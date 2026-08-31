/**
 * Web Worker for the cross-device stress test.
 *
 * Each ladder point runs in its own worker so that the driver can tell three
 * failure modes apart without losing the results collected so far:
 *   - the job throws (out of memory, allocation failure)  -> 'error'
 *   - the job outlives its time budget                    -> driver terminates -> 'timeout'
 *   - the worker itself dies (engine-level OOM kill)      -> 'worker-died'
 * Only a whole-tab crash escapes this, and the driver checkpoints to
 * localStorage so that case is recoverable too.
 */

import { parseXYZ } from './parseXYZ';
import { BasisSet } from './basisSet';
import { Molecular } from './molecular';
import { buildHF, type DFTConfig } from './builder';
import { initWasm, isWasmAvailable, getActiveBackend } from './eriWasm';
import type { EriBackend } from './hf';

export interface MemoryProbeRequest {
  type: 'memory-probe';
  /** Basis-function counts to try, ascending. */
  nbasisLadder: number[];
}

export type ProbeVerdict = 'resident' | 'swapping' | 'failed' | 'skipped';

export interface MemoryProbeResult {
  type: 'probe-result';
  nbasis: number;
  bytes: number;
  /** The allocation itself succeeded — see `verdict` for whether it was usable. */
  ok: boolean;
  ms: number;
  /** Commit throughput; collapses once the allocation is served by the page file. */
  gbPerSec: number;
  verdict: ProbeVerdict;
  detail?: string;
}

export interface MemoryProbeDone { type: 'probe-done'; }

export interface StressRequest {
  type: 'stress-run';
  xyzText: string;
  basisGBS: string;
  charge: number;
  eriBackend: EriBackend;
  baseUrl: string;
  dftConfig?: DFTConfig;
}

export interface StressProgress {
  type: 'progress';
  message: string;
  /** ms since the worker started the job — lets the driver show a live clock. */
  elapsedMs: number;
}

export interface StressDone {
  type: 'done';
  nbasis: number;
  natoms: number;
  energy: number;
  iterations: number;
  converged: boolean;
  setupMs: number;
  scfMs: number;
  totalMs: number;
  backend: string;
  /** Bytes of the unique-ERI array: (npair(npair+1)/2) doubles. */
  eriBytes: number;
  heapUsedBytes: number | null;
  heapLimitBytes: number | null;
}

export interface StressError {
  type: 'error';
  name: string;
  message: string;
  /** Set when the failure happened before the SCF loop began. */
  phase: 'setup' | 'eri' | 'scf' | 'unknown';
  nbasis: number | null;
  elapsedMs: number;
}

export type StressResponse = StressProgress | StressDone | StressError | MemoryProbeResult | MemoryProbeDone;
export type StressWorkerRequest = StressRequest | MemoryProbeRequest;

function heap(): { used: number | null; limit: number | null } {
  const p = performance as Performance & {
    memory?: { usedJSHeapSize: number; jsHeapSizeLimit: number };
  };
  if (!p.memory) return { used: null, limit: null };
  return { used: p.memory.usedJSHeapSize, limit: p.memory.jsHeapSizeLimit };
}

function eriBytesFor(n: number): number {
  const npair = n * (n + 1) / 2;
  return (npair * (npair + 1) / 2) * 8;
}

/**
 * Largest ERI array the device can hold *in real memory*.
 *
 * The SCF's dominant allocation is the unique-ERI array, so the memory ceiling
 * can be found without running any chemistry — which matters because an SCF at
 * that size would take tens of minutes.
 *
 * Plain allocation success is not the answer: on a desktop OS a multi-gigabyte
 * typed array is happily backed by the page file, so `new Float64Array` keeps
 * succeeding long past the point where the machine is thrashing (measured: a
 * 121 GB array "allocated" on a 32 GB desktop, at 1.4 GB/s and 84 s of disk
 * grinding). So the probe writes one value per 4 KiB page and watches the
 * commit throughput: once it collapses relative to the small-array baseline we
 * are paging, not computing, and the probe stops there.
 */
function probeMemory(req: MemoryProbeRequest): void {
  const post = (m: StressResponse) => (self as unknown as Worker).postMessage(m);

  const nav = navigator as Navigator & { deviceMemory?: number };
  // Stay well inside physical RAM: anything above it is swap by definition, and a
  // browser reporting 32 GB would otherwise let the probe commit tens of GB and
  // drag the whole machine into a paging storm.
  const hintGB = nav.deviceMemory ?? 4;
  const hardCapBytes = Math.min(Math.max(1, hintGB * 0.5), 8) * 1024 ** 3;

  let baselineGBps = 0;

  for (const nb of req.nbasisLadder) {
    const bytes = eriBytesFor(nb);
    const len = bytes / 8;

    if (bytes > hardCapBytes) {
      post({
        type: 'probe-result', nbasis: nb, bytes, ok: false, ms: 0, gbPerSec: 0,
        verdict: 'skipped',
        detail: `beyond the ${(hardCapBytes / 1024 ** 3).toFixed(0)} GB probe cap — far past any usable working set`,
      });
      break;
    }

    const t0 = performance.now();
    try {
      if (len > Number.MAX_SAFE_INTEGER) throw new RangeError('array length exceeds addressable range');
      let buf: Float64Array | null = new Float64Array(len);
      // One write per 4 KiB page (512 doubles) forces the pages to be committed.
      for (let i = 0; i < len; i += 512) buf[i] = 1;
      if (!Number.isFinite(buf[0])) throw new Error('unexpected buffer contents');
      buf = null;

      const ms = performance.now() - t0;
      const gbPerSec = (bytes / 1024 ** 3) / (ms / 1000);
      // Calibrate on the first array big enough to be measured reliably.
      if (baselineGBps === 0 && bytes >= 64 * 1024 ** 2 && ms > 5) baselineGBps = gbPerSec;

      const swapping = baselineGBps > 0 && gbPerSec < 0.4 * baselineGBps;
      post({
        type: 'probe-result', nbasis: nb, bytes, ok: true, ms, gbPerSec,
        verdict: swapping ? 'swapping' : 'resident',
        ...(swapping
          ? { detail: `commit throughput fell to ${gbPerSec.toFixed(2)} GB/s from a ${baselineGBps.toFixed(2)} GB/s baseline — this array no longer fits in real memory` }
          : {}),
      });
      if (swapping) break;
    } catch (err) {
      post({
        type: 'probe-result', nbasis: nb, bytes, ok: false,
        ms: performance.now() - t0, gbPerSec: 0, verdict: 'failed',
        detail: err instanceof Error ? `${err.name}: ${err.message}` : String(err),
      });
      break;
    }
  }
  post({ type: 'probe-done' });
}

self.onmessage = async (e: MessageEvent<StressWorkerRequest>) => {
  const req = e.data;
  if (req.type === 'memory-probe') { probeMemory(req); return; }
  if (req.type !== 'stress-run') return;

  const t0 = performance.now();
  let phase: StressError['phase'] = 'setup';
  let nbasis: number | null = null;
  const post = (m: StressResponse) => (self as unknown as Worker).postMessage(m);
  const progress = (message: string) =>
    post({ type: 'progress', message, elapsedMs: performance.now() - t0 });

  try {
    await initWasm(req.baseUrl);
    const backend = isWasmAvailable() ? (getActiveBackend() ?? 'wasm') : 'js';

    progress('Building molecule…');
    const basis = BasisSet.fromGBS(req.basisGBS);
    const atoms = parseXYZ(req.xyzText);
    const mol = new Molecular(atoms, basis, req.charge);
    nbasis = mol.numBasis;
    progress(`${atoms.length} atoms, ${nbasis} basis functions`);

    const setupMs = performance.now() - t0;
    const hf = buildHF(mol, 'RHF', req.dftConfig);

    let iterations = 0;
    let converged = false;
    phase = 'eri';
    const tScf = performance.now();
    const energy = await hf.solve({
      eriBackend: req.eriBackend,
      onProgress: (m) => {
        if (m.startsWith('Computing ERIs')) phase = 'eri';
        else if (m.startsWith('Converged')) converged = true;
        progress(m);
      },
      onIteration: (iter) => {
        phase = 'scf';
        iterations = iter;
        progress(`SCF iteration ${iter}`);
      },
    });

    const h = heap();
    post({
      type: 'done',
      nbasis,
      natoms: atoms.length,
      energy,
      iterations,
      converged,
      setupMs,
      scfMs: performance.now() - tScf,
      totalMs: performance.now() - t0,
      backend,
      eriBytes: eriBytesFor(nbasis),
      heapUsedBytes: h.used,
      heapLimitBytes: h.limit,
    });
  } catch (err) {
    post({
      type: 'error',
      name: err instanceof Error ? err.name : 'Error',
      message: err instanceof Error ? err.message : String(err),
      phase,
      nbasis,
      elapsedMs: performance.now() - t0,
    });
  }
};
