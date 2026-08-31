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
  baseUrl: string;
  /**
   * How long to keep the pair resident before releasing it. iOS terminates tabs
   * on sustained memory pressure rather than on a single allocation, so a probe
   * that allocates and frees within milliseconds passes where the real SCF —
   * which holds the same memory for minutes — is killed.
   */
  holdMs: number;
}

export type ProbeVerdict = 'ok' | 'failed' | 'unavailable';
/** Which half of the real peak footprint the probe was holding when it failed. */
export type ProbeStage = 'wasm' | 'js-copy';

export interface MemoryProbeResult {
  type: 'probe-result';
  nbasis: number;
  /** Size of one ERI array. The SCF holds two of these at once — see peakBytes. */
  bytes: number;
  /** What the calculation actually needs resident: the wasm array plus its JS copy. */
  peakBytes: number;
  ok: boolean;
  ms: number;
  verdict: ProbeVerdict;
  failedStage?: ProbeStage;
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
 * Minimal direct instantiation of the ERI module, exposing just the allocator.
 * `initWasm` wraps the exports and never hands out raw malloc, but the probe
 * needs to allocate inside the module's own linear memory to be meaningful.
 */
async function loadWasmAllocator(baseUrl: string): Promise<{
  malloc: (n: number, align: number) => number;
  free: (ptr: number, n: number, align: number) => void;
  memory: WebAssembly.Memory;
} | null> {
  for (const name of ['wasm_eri_simd_bg.wasm', 'wasm_eri_bg.wasm']) {
    try {
      const r = await fetch(new URL(`${baseUrl}wasm/${name}`, self.location.origin));
      if (!r.ok) continue;
      const mod = await WebAssembly.compile(await r.arrayBuffer());
      let instance: WebAssembly.Instance;
      instance = await WebAssembly.instantiate(mod, {
        './wasm_eri_bg.js': {
          __wbindgen_init_externref_table: () => {
            const t = (instance.exports as unknown as { __wbindgen_externrefs: WebAssembly.Table }).__wbindgen_externrefs;
            const o = t.grow(4);
            t.set(0, undefined); t.set(o + 0, undefined); t.set(o + 1, null); t.set(o + 2, true); t.set(o + 3, false);
          },
        },
      });
      const ex = instance.exports as unknown as Record<string, Function> & { memory: WebAssembly.Memory };
      if (ex.__wbindgen_start) ex.__wbindgen_start();
      return {
        malloc: ex.__wbindgen_malloc as (n: number, a: number) => number,
        free: ex.__wbindgen_free as (p: number, n: number, a: number) => void,
        memory: ex.memory,
      };
    } catch { /* try the next binary */ }
  }
  return null;
}

/**
 * Largest system the WASM backend can actually hold.
 *
 * This reproduces the calculation's real peak footprint, which is *two* copies of
 * the ERI array, not one: the Rust side allocates it inside wasm linear memory,
 * then `readF64` slices it into a JS `Float64Array` — and wasm linear memory
 * never shrinks, so both stay resident for the rest of the SCF.
 *
 * Getting this wrong is not academic. An earlier version probed only the wasm
 * allocation and reported the same 200-basis-function ceiling on an iPhone as on
 * a 32 GB desktop, while the iPhone's tab was in fact being killed at 175 — the
 * missing JS copy was the difference. Two ceilings therefore exist and the probe
 * has to distinguish them:
 *   - wasm32's 4 GiB address space, which is what stops a desktop; and
 *   - the device's own memory cap, which is what stops a phone.
 */
async function probeMemory(req: MemoryProbeRequest): Promise<void> {
  const post = (m: StressResponse) => (self as unknown as Worker).postMessage(m);

  const w = await loadWasmAllocator(req.baseUrl);
  if (!w) {
    post({
      type: 'probe-result', nbasis: 0, bytes: 0, peakBytes: 0, ok: false, ms: 0, verdict: 'unavailable',
      detail: 'could not instantiate the WASM module — this device runs the JS backend, whose ceiling is device RAM rather than the wasm32 address space',
    });
    post({ type: 'probe-done' });
    return;
  }

  for (const nb of req.nbasisLadder) {
    const bytes = eriBytesFor(nb);
    const t0 = performance.now();
    let stage: ProbeStage = 'wasm';
    let ptr = 0;
    let jsCopy: Float64Array | null = null;
    try {
      ptr = w.malloc(bytes, 8) >>> 0;
      if (!ptr) throw new Error('allocator returned null');
      // memory.buffer is detached by any growth, so take the view after malloc.
      const view = new Float64Array(w.memory.buffer, ptr, bytes / 8);
      for (let i = 0; i < view.length; i += 512) view[i] = 1;

      // The JS-side copy the real code makes, held at the same time.
      stage = 'js-copy';
      jsCopy = new Float64Array(bytes / 8);
      for (let i = 0; i < jsCopy.length; i += 512) jsCopy[i] = view[i];
      if (!Number.isFinite(jsCopy[0])) throw new Error('unexpected buffer contents');

      // Keep both resident, touching them so the pages cannot be reclaimed, for
      // long enough that a pressure-based killer has the chance it would get
      // during the real SCF.
      const holdUntil = performance.now() + req.holdMs;
      while (performance.now() < holdUntil) {
        for (let i = 0; i < jsCopy.length; i += 65536) jsCopy[i] += 1;
        for (let i = 0; i < view.length; i += 65536) view[i] += 1;
      }

      post({
        type: 'probe-result', nbasis: nb, bytes, peakBytes: bytes * 2,
        ok: true, ms: performance.now() - t0, verdict: 'ok',
      });
    } catch (err) {
      post({
        type: 'probe-result', nbasis: nb, bytes, peakBytes: bytes * 2,
        ok: false, ms: performance.now() - t0, verdict: 'failed', failedStage: stage,
        detail: err instanceof Error ? `${err.name}: ${err.message}` : String(err),
      });
      // A Rust abort leaves the module poisoned; nothing past this point is meaningful.
      break;
    } finally {
      jsCopy = null;
      if (ptr) { try { w.free(ptr, bytes, 8); } catch { /* poisoned module */ } }
    }
  }
  post({ type: 'probe-done' });
}

self.onmessage = async (e: MessageEvent<StressWorkerRequest>) => {
  const req = e.data;
  if (req.type === 'memory-probe') { await probeMemory(req); return; }
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
