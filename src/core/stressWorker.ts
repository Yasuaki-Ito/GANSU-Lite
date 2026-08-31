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

export type StressResponse = StressProgress | StressDone | StressError;

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

self.onmessage = async (e: MessageEvent<StressRequest>) => {
  const req = e.data;
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
