/**
 * Web Worker for Accuracy page — runs all 8 methods off the main thread.
 */

import { parseXYZ } from './parseXYZ';
import { BasisSet } from './basisSet';
import { Molecular } from './molecular';
import { buildHF } from './builder';
import type { DFTConfig } from './builder';
import { RHF } from './rhf';
import { computeMP2Energy } from './mp2';
import { computeMP3Energy } from './mp3';
import { computeCCSDEnergy } from './ccsd';
import {
  initWasm, isWasmAvailable, getActiveBackend,
  computeMP2EnergyWasm, computeMP3EnergyWasm, computeCCSDEnergyWasm,
} from './eriWasm';
import { getWasmPostHFThreshold } from './calibration';
import type { FunctionalName } from './xcFunctional';

// ── Message protocol ────────────────────────────────────────────────

export interface AccWorkerRequest {
  type: 'run';
  xyzText: string;
  basisGBS: string;
  charge: number;
  multiplicity: number;
  selectedMethods: string[];
  baseUrl: string;
}

/** Progress: method state change */
export interface AccWorkerMethodState {
  type: 'method-state';
  index: number;
  state: 'pending' | 'running' | 'done' | 'skipped';
  info?: string;
}

/** Progress: method info text update */
export interface AccWorkerMethodInfo {
  type: 'method-info';
  index: number;
  msg: string;
}

/** Progress: method completed with timing */
export interface AccWorkerMethodTime {
  type: 'method-time';
  index: number;
  ms: number;
}

/** Backend detection */
export interface AccWorkerBackend {
  type: 'backend';
  backend: 'wasm-simd' | 'wasm' | 'js' | 'none';
}

/** One method result */
export interface AccWorkerResult {
  type: 'method-result';
  methodIndex: number;
  totalEnergy: number;
  diffFromHF: number;
  timeMs: number;
}

/** All done */
export interface AccWorkerDone {
  type: 'done';
}

/** Error */
export interface AccWorkerError {
  type: 'error';
  message: string;
}

export type AccWorkerResponse =
  | AccWorkerMethodState
  | AccWorkerMethodInfo
  | AccWorkerMethodTime
  | AccWorkerBackend
  | AccWorkerResult
  | AccWorkerDone
  | AccWorkerError;

// ── Method definitions (must match accuracy.ts) ─────────────────────

interface MethodDef {
  id: string;
  type: 'hf' | 'postHF' | 'dft';
  functional?: FunctionalName;
}

const ALL_METHODS: MethodDef[] = [
  { id: 'HF',    type: 'hf' },
  { id: 'MP2',   type: 'postHF' },
  { id: 'MP3',   type: 'postHF' },
  { id: 'CCSD',  type: 'postHF' },
  { id: 'SVWN',  type: 'dft', functional: 'SVWN' },
  { id: 'BLYP',  type: 'dft', functional: 'BLYP' },
  { id: 'PBE',   type: 'dft', functional: 'PBE' },
  { id: 'B3LYP', type: 'dft', functional: 'B3LYP' },
];

// ── Worker handler ──────────────────────────────────────────────────

self.onmessage = async (e: MessageEvent<AccWorkerRequest>) => {
  const { xyzText, basisGBS, charge, multiplicity, selectedMethods, baseUrl } = e.data;
  const selected = new Set(selectedMethods);
  const post = (msg: AccWorkerResponse) => self.postMessage(msg);

  try {
    // Init WASM
    await initWasm(baseUrl);
    const backend = getActiveBackend();
    post({ type: 'backend', backend });

    // Parse molecule + basis
    const atoms = parseXYZ(xyzText);
    const basisSet = BasisSet.fromGBS(basisGBS);
    const betaToAlpha = Math.floor((multiplicity - 1) / 2);
    const molecular = new Molecular(atoms, basisSet, charge, betaToAlpha);

    // ── HF (always required) ──
    post({ type: 'method-state', index: 0, state: 'running', info: 'starting...' });

    const t0 = performance.now();
    const rhf = buildHF(molecular, 'RHF') as RHF;
    const eHF = await rhf.solve({
      eriBackend: 'auto',
      onProgress: (msg) => {
        // Detect backend from progress tags
        if (msg.includes('[WASM]') || msg.includes('[JS]')) {
          post({ type: 'backend', backend: getActiveBackend() });
        }
        post({ type: 'method-info', index: 0, msg });
      },
    });
    const hfTime = performance.now() - t0;

    post({ type: 'method-time', index: 0, ms: hfTime });
    post({ type: 'method-state', index: 0, state: 'done' });
    post({ type: 'method-result', methodIndex: 0, totalEnergy: eHF, diffFromHF: 0, timeMs: hfTime });

    // Extract RHF data for post-HF
    const C = rhf.coefficients;
    const eps = rhf.orbitalEnergies;
    const eri = rhf.eriStore;
    const nocc = molecular.numAlphaSpins;
    const nbasis = molecular.numBasis;
    const useWasmPH = isWasmAvailable() && nbasis >= getWasmPostHFThreshold();

    // ── MP2 ──
    if (selected.has('MP2')) {
      post({ type: 'method-state', index: 1, state: 'running', info: 'MO transform...' });

      const t1 = performance.now();
      let eCorrMP2: number;
      if (useWasmPH) {
        eCorrMP2 = computeMP2EnergyWasm(eri.data, C.data, eps, nocc, nbasis)
          ?? computeMP2Energy(C, eps, eri, nocc, nbasis);
      } else {
        eCorrMP2 = computeMP2Energy(C, eps, eri, nocc, nbasis);
      }
      const mp2Time = performance.now() - t1;

      post({ type: 'method-time', index: 1, ms: mp2Time });
      post({ type: 'method-state', index: 1, state: 'done' });
      post({ type: 'method-result', methodIndex: 1, totalEnergy: eHF + eCorrMP2, diffFromHF: eCorrMP2 * 1000, timeMs: mp2Time });
    }

    // ── MP3 ──
    if (selected.has('MP3')) {
      post({ type: 'method-state', index: 2, state: 'running', info: 'MO transform...' });

      const t2 = performance.now();
      let totalCorrMP3: number;
      if (useWasmPH) {
        const wasmResult = computeMP3EnergyWasm(eri.data, C.data, eps, nocc, nbasis);
        if (wasmResult) {
          totalCorrMP3 = wasmResult.mp2 + wasmResult.mp3;
        } else {
          const mp3Result = computeMP3Energy(C, eps, eri, nocc, nbasis);
          totalCorrMP3 = mp3Result.mp2 + mp3Result.mp3;
        }
      } else {
        const mp3Result = computeMP3Energy(C, eps, eri, nocc, nbasis);
        totalCorrMP3 = mp3Result.mp2 + mp3Result.mp3;
      }
      const mp3Time = performance.now() - t2;

      post({ type: 'method-time', index: 2, ms: mp3Time });
      post({ type: 'method-state', index: 2, state: 'done' });
      post({ type: 'method-result', methodIndex: 2, totalEnergy: eHF + totalCorrMP3, diffFromHF: totalCorrMP3 * 1000, timeMs: mp3Time });
    }

    // ── CCSD ──
    if (selected.has('CCSD')) {
      post({ type: 'method-state', index: 3, state: 'running', info: 'T/V amplitudes...' });

      const t3 = performance.now();
      let eCorrCCSD: number;
      if (useWasmPH) {
        eCorrCCSD = computeCCSDEnergyWasm(eri.data, C.data, eps, nocc, nbasis)
          ?? computeCCSDEnergy(C, eps, eri, nocc, nbasis);
      } else {
        eCorrCCSD = computeCCSDEnergy(C, eps, eri, nocc, nbasis);
      }
      const ccsdTime = performance.now() - t3;

      post({ type: 'method-time', index: 3, ms: ccsdTime });
      post({ type: 'method-state', index: 3, state: 'done' });
      post({ type: 'method-result', methodIndex: 3, totalEnergy: eHF + eCorrCCSD, diffFromHF: eCorrCCSD * 1000, timeMs: ccsdTime });
    }

    // ── DFT methods (4-7) ──
    for (let i = 4; i < ALL_METHODS.length; i++) {
      const meth = ALL_METHODS[i];
      if (!selected.has(meth.id)) continue;

      post({ type: 'method-state', index: i, state: 'running', info: 'starting...' });

      const td = performance.now();
      const dftConfig: DFTConfig = { functional: meth.functional!, gridLevel: 'medium' };
      const hf = buildHF(molecular, 'RHF', dftConfig);
      const energy = await hf.solve({
        eriBackend: 'auto',
        onProgress: (msg) => {
          if (msg.includes('[WASM]') || msg.includes('[JS]')) {
            post({ type: 'backend', backend: getActiveBackend() });
          }
          post({ type: 'method-info', index: i, msg });
        },
      });
      const dftTime = performance.now() - td;

      post({ type: 'method-time', index: i, ms: dftTime });
      post({ type: 'method-state', index: i, state: 'done' });
      post({ type: 'method-result', methodIndex: i, totalEnergy: energy, diffFromHF: (energy - eHF) * 1000, timeMs: dftTime });
    }

    post({ type: 'done' });
  } catch (err) {
    post({ type: 'error', message: err instanceof Error ? err.message : String(err) });
  }
};
