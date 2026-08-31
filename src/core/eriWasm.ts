/** WASM-accelerated ERI computation with JS fallback */

import type { PrimitiveShell } from './types';
import type { FloatArray } from '../linalg/matrix';

/** Pack PrimitiveShell[] into a flat Float64Array for WASM transfer.
 *  Layout: 7 fields per shell [exponent, coefficient, x, y, z, shellType, basisIndex] */
export function packShells(shells: PrimitiveShell[]): Float64Array {
  const flat = new Float64Array(shells.length * 7);
  for (let i = 0; i < shells.length; i++) {
    const off = i * 7;
    const s = shells[i];
    flat[off] = s.exponent;
    flat[off + 1] = s.coefficient;
    flat[off + 2] = s.coordinate.x;
    flat[off + 3] = s.coordinate.y;
    flat[off + 4] = s.coordinate.z;
    flat[off + 5] = s.shellType;
    flat[off + 6] = s.basisIndex;
  }
  return flat;
}

/** Pack PrimitiveShell[] into extended flat format (8 fields, includes atomIndex). */
export function packShellsExt(shells: PrimitiveShell[]): Float64Array {
  const flat = new Float64Array(shells.length * 8);
  for (let i = 0; i < shells.length; i++) {
    const off = i * 8;
    const s = shells[i];
    flat[off] = s.exponent;
    flat[off + 1] = s.coefficient;
    flat[off + 2] = s.coordinate.x;
    flat[off + 3] = s.coordinate.y;
    flat[off + 4] = s.coordinate.z;
    flat[off + 5] = s.shellType;
    flat[off + 6] = s.basisIndex;
    flat[off + 7] = s.atomIndex;
  }
  return flat;
}

type ComputeErisWasmFn = (
  shellsFlat: Float64Array,
  normFactors: Float64Array,
  numBasis: number,
  schwarzThreshold: number,
) => Float64Array;

type ComputeFockRhfWasmFn = (
  eriData: FloatArray, density: FloatArray,
  coreH: FloatArray, numBasis: number,
) => Float64Array;

type ComputeFockUhfWasmFn = (
  eriData: FloatArray, densityAlpha: FloatArray,
  densityBeta: FloatArray, densityTotal: FloatArray,
  coreH: FloatArray, numBasis: number,
) => Float64Array;

type ComputeMP2WasmFn = (
  eriData: FloatArray, coeff: FloatArray,
  epsilon: FloatArray, nocc: number, nbasis: number,
) => number;

type ComputeMP3WasmFn = (
  eriData: FloatArray, coeff: FloatArray,
  epsilon: FloatArray, nocc: number, nbasis: number,
) => { mp2: number; mp3: number };

type ComputeCCSDWasmFn = (
  eriData: FloatArray, coeff: FloatArray,
  epsilon: FloatArray, nocc: number, nbasis: number,
) => number;

type ComputeUMP2WasmFn = (
  eriData: FloatArray, ca: FloatArray, cb: FloatArray,
  epsilonA: FloatArray, epsilonB: FloatArray,
  noccA: number, noccB: number, nbasis: number,
) => number;

type ComputeUMP3WasmFn = (
  eriData: FloatArray, ca: FloatArray, cb: FloatArray,
  epsilonA: FloatArray, epsilonB: FloatArray,
  noccA: number, noccB: number, nbasis: number,
) => { mp2: number; mp3: number };

type ComputeUCCSDWasmFn = (
  eriData: FloatArray, ca: FloatArray, cb: FloatArray,
  epsilonA: FloatArray, epsilonB: FloatArray,
  noccA: number, noccB: number, nbasis: number,
) => number;

// RI WASM function types
type ComputeRISetupWasmFn = (
  primaryShellsFlat: FloatArray, primaryNorm: FloatArray,
  nbasis: number,
  auxShellsFlat: FloatArray, auxNorm: FloatArray,
  naux: number,
) => Float64Array;

type ComputeRIFockRhfWasmFn = (
  bMatrix: FloatArray, density: FloatArray,
  coreH: FloatArray, coefficients: FloatArray,
  naux: number, nbasis: number, nocc: number,
) => Float64Array;

type ComputeRIFockUhfWasmFn = (
  bMatrix: FloatArray, densityTotal: FloatArray,
  coreH: FloatArray,
  coeffAlpha: FloatArray, noccAlpha: number,
  coeffBeta: FloatArray, noccBeta: number,
  naux: number, nbasis: number,
) => Float64Array;

type ComputeRIMP2WasmFn = (
  bMatrix: FloatArray, coefficients: FloatArray,
  epsilon: FloatArray,
  naux: number, nocc: number, nbasis: number,
) => number;

type ComputeRIUMP2WasmFn = (
  bMatrix: FloatArray, ca: FloatArray, cb: FloatArray,
  epsilonA: FloatArray, epsilonB: FloatArray,
  naux: number, noccA: number, noccB: number, nbasis: number,
) => number;

let _wasmComputeEris: ComputeErisWasmFn | null = null;
let _wasmComputeFockRhf: ComputeFockRhfWasmFn | null = null;
let _wasmComputeFockUhf: ComputeFockUhfWasmFn | null = null;
let _wasmComputeMP2: ComputeMP2WasmFn | null = null;
let _wasmComputeMP3: ComputeMP3WasmFn | null = null;
let _wasmComputeCCSD: ComputeCCSDWasmFn | null = null;
let _wasmComputeUMP2: ComputeUMP2WasmFn | null = null;
let _wasmComputeUMP3: ComputeUMP3WasmFn | null = null;
let _wasmComputeUCCSD: ComputeUCCSDWasmFn | null = null;
let _wasmComputeRISetup: ComputeRISetupWasmFn | null = null;
let _wasmComputeRIFockRhf: ComputeRIFockRhfWasmFn | null = null;
let _wasmComputeRIFockUhf: ComputeRIFockUhfWasmFn | null = null;
let _wasmComputeRIMP2: ComputeRIMP2WasmFn | null = null;
let _wasmComputeRIUMP2: ComputeRIUMP2WasmFn | null = null;

type ComputeXCWasmFn = (
  shellsFlat: Float64Array, normFactors: Float64Array,
  densityA: Float64Array, densityB: Float64Array,
  gridFlat: Float64Array, nbasis: number,
  funcId: number, needGrad: boolean,
) => Float64Array;
let _wasmComputeXC: ComputeXCWasmFn | null = null;

// Gradient / Hessian / CPHF WASM function types
type Compute2eGradientWasmFn = (
  shellsFlat: Float64Array, density: Float64Array,
  norms: Float64Array, numBasis: number, numAtoms: number,
) => Float64Array;

type Compute2eHessianWasmFn = (
  shellsFlat: Float64Array, density: Float64Array,
  norms: Float64Array, numBasis: number, numAtoms: number,
) => Float64Array;

type TransformERItoMOWasmFn = (
  cFlat: Float64Array, eriSym: Float64Array, numBasis: number,
) => Float64Array;

type SolveCPHFWasmFn = (
  moERI: Float64Array, eps: Float64Array, rhs: Float64Array,
  nocc: number, nvir: number, nmo: number, nPert: number,
  tol: number, maxIter: number,
) => Float64Array;

let _wasmCompute2eGradient: Compute2eGradientWasmFn | null = null;
let _wasmCompute2eHessian: Compute2eHessianWasmFn | null = null;
let _wasmTransformERItoMO: TransformERItoMOWasmFn | null = null;
let _wasmSolveCPHF: SolveCPHFWasmFn | null = null;

let _wasmSimdAvailable = false;
let _wasmInitPromise: Promise<boolean> | null = null;

/** Detect WASM SIMD support by trying to compile a v128.const instruction */
async function detectWasmSimd(): Promise<boolean> {
  try {
    // Use WebAssembly.validate for fast feature detection
    // Module: (func (result v128) v128.const i64x2 0 0)
    const simdTest = new Uint8Array([
      0x00, 0x61, 0x73, 0x6d,  // \0asm magic
      0x01, 0x00, 0x00, 0x00,  // version 1
      // Type section (id=1)
      0x01, 0x05, 0x01,        // section, 5 bytes, 1 type
      0x60, 0x00, 0x01, 0x7b,  // func() -> v128
      // Function section (id=3)
      0x03, 0x02, 0x01, 0x00,  // section, 2 bytes, 1 func, type 0
      // Code section (id=10)
      0x0a, 0x16, 0x01,        // section, 22 bytes, 1 body
      0x14, 0x00,              // body: 20 bytes, 0 locals
      0xfd, 0x0c,              // v128.const
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x0b,                    // end
    ]);
    return WebAssembly.validate(simdTest);
  } catch {
    return false;
  }
}

/** Check if WASM SIMD backend is active */
export function isWasmSimdAvailable(): boolean {
  return _wasmSimdAvailable;
}

/** Get the active backend tier: 'wasm-simd', 'wasm', or 'js' */
export function getActiveBackend(): 'wasm-simd' | 'wasm' | 'js' {
  if (_wasmSimdAvailable) return 'wasm-simd';
  if (_wasmComputeEris !== null) return 'wasm';
  return 'js';
}

/** Initialize WASM module. Returns true if successful, false if fallback needed. */
export function initWasm(baseUrl?: string): Promise<boolean> {
  if (_wasmInitPromise) return _wasmInitPromise;

  _wasmInitPromise = (async () => {
    try {
      const base = baseUrl || '/';
      const simdSupported = await detectWasmSimd();

      // Try SIMD binary first, then fall back to standard WASM
      let wasmUrl: URL;
      if (simdSupported) {
        wasmUrl = new URL(`${base}wasm/wasm_eri_simd_bg.wasm`, self.location.origin);
        // If SIMD binary doesn't exist yet, fall back
        try {
          const probe = await fetch(wasmUrl, { method: 'HEAD' });
          if (!probe.ok) throw new Error('SIMD binary not found');
        } catch {
          wasmUrl = new URL(`${base}wasm/wasm_eri_bg.wasm`, self.location.origin);
        }
      } else {
        wasmUrl = new URL(`${base}wasm/wasm_eri_bg.wasm`, self.location.origin);
      }

      // Fetch and compile the WASM module
      const response = await fetch(wasmUrl);
      const bytes = await response.arrayBuffer();
      const module = await WebAssembly.compile(bytes);

      // Manually set up wasm-bindgen glue
      const imports = {
        './wasm_eri_bg.js': {
          __wbindgen_init_externref_table: () => {
            const table = (instance.exports as any).__wbindgen_externrefs;
            const offset = table.grow(4);
            table.set(0, undefined);
            table.set(offset + 0, undefined);
            table.set(offset + 1, null);
            table.set(offset + 2, true);
            table.set(offset + 3, false);
          },
        },
      };

      const instance = await WebAssembly.instantiate(module, imports);
      const exports = instance.exports as any;

      // Call __wbindgen_start to initialize
      if (exports.__wbindgen_start) {
        exports.__wbindgen_start();
      }

      // Wrap the raw WASM function with proper memory management
      const memory = exports.memory as WebAssembly.Memory;
      const malloc = exports.__wbindgen_malloc as (size: number, align: number) => number;
      const free = exports.__wbindgen_free as (ptr: number, size: number, align: number) => void;
      const wasmComputeEris = exports.compute_eris_wasm as (
        shellsPtr: number, shellsLen: number,
        normsPtr: number, normsLen: number,
        numBasis: number, threshold: number,
      ) => [number, number];

      // Helper: write f64 array to WASM memory, return [ptr, len]
      const writeF64 = (data: FloatArray): [number, number] => {
        const ptr = malloc(data.length * 8, 8) >>> 0;
        new Float64Array(memory.buffer).set(data, ptr / 8);
        return [ptr, data.length];
      };

      // Helper: read f64 result from WASM memory and free
      const readF64Result = (result: [number, number]): Float64Array => {
        const [ptr, len] = result;
        const output = new Float64Array(memory.buffer).slice(ptr / 8, ptr / 8 + len);
        free(ptr, len * 8, 8);
        return output;
      };

      _wasmComputeEris = (shellsFlat, normFactors, numBasis, schwarzThreshold) => {
        const [shellsPtr, shellsLen] = writeF64(shellsFlat);
        const [normsPtr, normsLen] = writeF64(normFactors);
        const result = wasmComputeEris(shellsPtr, shellsLen, normsPtr, normsLen, numBasis, schwarzThreshold);
        return readF64Result(result);
      };

      // Fock matrix WASM functions
      const wasmFockRhf = exports.compute_fock_rhf as (
        eriPtr: number, eriLen: number,
        densPtr: number, densLen: number,
        chPtr: number, chLen: number,
        numBasis: number,
      ) => [number, number];

      const wasmFockUhf = exports.compute_fock_uhf as (
        eriPtr: number, eriLen: number,
        daPtr: number, daLen: number,
        dbPtr: number, dbLen: number,
        dtPtr: number, dtLen: number,
        chPtr: number, chLen: number,
        numBasis: number,
      ) => [number, number];

      _wasmComputeFockRhf = (eriData, density, coreH, numBasis) => {
        const [eriPtr, eriLen] = writeF64(eriData);
        const [densPtr, densLen] = writeF64(density);
        const [chPtr, chLen] = writeF64(coreH);
        const result = wasmFockRhf(eriPtr, eriLen, densPtr, densLen, chPtr, chLen, numBasis);
        return readF64Result(result);
      };

      _wasmComputeFockUhf = (eriData, densityAlpha, densityBeta, densityTotal, coreH, numBasis) => {
        const [eriPtr, eriLen] = writeF64(eriData);
        const [daPtr, daLen] = writeF64(densityAlpha);
        const [dbPtr, dbLen] = writeF64(densityBeta);
        const [dtPtr, dtLen] = writeF64(densityTotal);
        const [chPtr, chLen] = writeF64(coreH);
        const result = wasmFockUhf(eriPtr, eriLen, daPtr, daLen, dbPtr, dbLen, dtPtr, dtLen, chPtr, chLen, numBasis);
        return readF64Result(result);
      };

      // MP2 WASM function
      const wasmMP2 = exports.compute_mp2_wasm as (
        eriPtr: number, eriLen: number,
        coeffPtr: number, coeffLen: number,
        epsPtr: number, epsLen: number,
        nocc: number, nbasis: number,
      ) => number;

      if (wasmMP2) {
        _wasmComputeMP2 = (eriData, coeff, epsilon, nocc, nbasis) => {
          const [eriPtr, eriLen] = writeF64(eriData);
          const [coeffPtr, coeffLen] = writeF64(coeff);
          const [epsPtr, epsLen] = writeF64(epsilon);
          return wasmMP2(eriPtr, eriLen, coeffPtr, coeffLen, epsPtr, epsLen, nocc, nbasis);
        };
      }

      // MP3 WASM function (returns Vec<f64> with 2 elements: [mp2, mp3])
      const wasmMP3 = exports.compute_mp3_wasm as (
        eriPtr: number, eriLen: number,
        coeffPtr: number, coeffLen: number,
        epsPtr: number, epsLen: number,
        nocc: number, nbasis: number,
      ) => [number, number];

      if (wasmMP3) {
        _wasmComputeMP3 = (eriData, coeff, epsilon, nocc, nbasis) => {
          const [eriPtr, eriLen] = writeF64(eriData);
          const [coeffPtr, coeffLen] = writeF64(coeff);
          const [epsPtr, epsLen] = writeF64(epsilon);
          const result = wasmMP3(eriPtr, eriLen, coeffPtr, coeffLen, epsPtr, epsLen, nocc, nbasis);
          const arr = readF64Result(result);
          return { mp2: arr[0], mp3: arr[1] };
        };
      }

      // CCSD WASM function (returns f64 scalar)
      const wasmCCSD = exports.compute_ccsd_wasm as (
        eriPtr: number, eriLen: number,
        coeffPtr: number, coeffLen: number,
        epsPtr: number, epsLen: number,
        nocc: number, nbasis: number,
      ) => number;

      if (wasmCCSD) {
        _wasmComputeCCSD = (eriData, coeff, epsilon, nocc, nbasis) => {
          const [eriPtr, eriLen] = writeF64(eriData);
          const [coeffPtr, coeffLen] = writeF64(coeff);
          const [epsPtr, epsLen] = writeF64(epsilon);
          return wasmCCSD(eriPtr, eriLen, coeffPtr, coeffLen, epsPtr, epsLen, nocc, nbasis);
        };
      }

      // UHF post-HF WASM functions
      const wasmUMP2 = exports.compute_ump2_wasm as (
        eriPtr: number, eriLen: number,
        caPtr: number, caLen: number,
        cbPtr: number, cbLen: number,
        epsAPtr: number, epsALen: number,
        epsBPtr: number, epsBLen: number,
        noccA: number, noccB: number, nbasis: number,
      ) => number;

      if (wasmUMP2) {
        _wasmComputeUMP2 = (eriData, ca, cb, epsilonA, epsilonB, noccA, noccB, nbasis) => {
          const [eriPtr, eriLen] = writeF64(eriData);
          const [caPtr, caLen] = writeF64(ca);
          const [cbPtr, cbLen] = writeF64(cb);
          const [epsAPtr, epsALen] = writeF64(epsilonA);
          const [epsBPtr, epsBLen] = writeF64(epsilonB);
          return wasmUMP2(eriPtr, eriLen, caPtr, caLen, cbPtr, cbLen, epsAPtr, epsALen, epsBPtr, epsBLen, noccA, noccB, nbasis);
        };
      }

      const wasmUMP3 = exports.compute_ump3_wasm as (
        eriPtr: number, eriLen: number,
        caPtr: number, caLen: number,
        cbPtr: number, cbLen: number,
        epsAPtr: number, epsALen: number,
        epsBPtr: number, epsBLen: number,
        noccA: number, noccB: number, nbasis: number,
      ) => [number, number];

      if (wasmUMP3) {
        _wasmComputeUMP3 = (eriData, ca, cb, epsilonA, epsilonB, noccA, noccB, nbasis) => {
          const [eriPtr, eriLen] = writeF64(eriData);
          const [caPtr, caLen] = writeF64(ca);
          const [cbPtr, cbLen] = writeF64(cb);
          const [epsAPtr, epsALen] = writeF64(epsilonA);
          const [epsBPtr, epsBLen] = writeF64(epsilonB);
          const result = wasmUMP3(eriPtr, eriLen, caPtr, caLen, cbPtr, cbLen, epsAPtr, epsALen, epsBPtr, epsBLen, noccA, noccB, nbasis);
          const arr = readF64Result(result);
          return { mp2: arr[0], mp3: arr[1] };
        };
      }

      const wasmUCCSD = exports.compute_uccsd_wasm as (
        eriPtr: number, eriLen: number,
        caPtr: number, caLen: number,
        cbPtr: number, cbLen: number,
        epsAPtr: number, epsALen: number,
        epsBPtr: number, epsBLen: number,
        noccA: number, noccB: number, nbasis: number,
      ) => number;

      if (wasmUCCSD) {
        _wasmComputeUCCSD = (eriData, ca, cb, epsilonA, epsilonB, noccA, noccB, nbasis) => {
          const [eriPtr, eriLen] = writeF64(eriData);
          const [caPtr, caLen] = writeF64(ca);
          const [cbPtr, cbLen] = writeF64(cb);
          const [epsAPtr, epsALen] = writeF64(epsilonA);
          const [epsBPtr, epsBLen] = writeF64(epsilonB);
          return wasmUCCSD(eriPtr, eriLen, caPtr, caLen, cbPtr, cbLen, epsAPtr, epsALen, epsBPtr, epsBLen, noccA, noccB, nbasis);
        };
      }

      // RI Setup WASM function
      const wasmRISetup = exports.compute_ri_setup_wasm as (
        pShellsPtr: number, pShellsLen: number,
        pNormPtr: number, pNormLen: number,
        nbasis: number,
        aShellsPtr: number, aShellsLen: number,
        aNormPtr: number, aNormLen: number,
        naux: number,
      ) => [number, number];

      if (wasmRISetup) {
        _wasmComputeRISetup = (primaryShellsFlat, primaryNorm, nbasis, auxShellsFlat, auxNorm, naux) => {
          const [pShellsPtr, pShellsLen] = writeF64(primaryShellsFlat);
          const [pNormPtr, pNormLen] = writeF64(primaryNorm);
          const [aShellsPtr, aShellsLen] = writeF64(auxShellsFlat);
          const [aNormPtr, aNormLen] = writeF64(auxNorm);
          const result = wasmRISetup(pShellsPtr, pShellsLen, pNormPtr, pNormLen, nbasis, aShellsPtr, aShellsLen, aNormPtr, aNormLen, naux);
          return readF64Result(result);
        };
      }

      // RI WASM functions
      const wasmRIFockRhf = exports.compute_ri_fock_rhf_wasm as (
        bPtr: number, bLen: number,
        densPtr: number, densLen: number,
        chPtr: number, chLen: number,
        coeffPtr: number, coeffLen: number,
        naux: number, nbasis: number, nocc: number,
      ) => [number, number];

      if (wasmRIFockRhf) {
        _wasmComputeRIFockRhf = (bMatrix, density, coreH, coefficients, naux, nbasis, nocc) => {
          const [bPtr, bLen] = writeF64(bMatrix);
          const [densPtr, densLen] = writeF64(density);
          const [chPtr, chLen] = writeF64(coreH);
          const [coeffPtr, coeffLen] = writeF64(coefficients);
          const result = wasmRIFockRhf(bPtr, bLen, densPtr, densLen, chPtr, chLen, coeffPtr, coeffLen, naux, nbasis, nocc);
          return readF64Result(result);
        };
      }

      const wasmRIFockUhf = exports.compute_ri_fock_uhf_wasm as (
        bPtr: number, bLen: number,
        dtPtr: number, dtLen: number,
        chPtr: number, chLen: number,
        caPtr: number, caLen: number,
        noccA: number,
        cbPtr: number, cbLen: number,
        noccB: number,
        naux: number, nbasis: number,
      ) => [number, number];

      if (wasmRIFockUhf) {
        _wasmComputeRIFockUhf = (bMatrix, densityTotal, coreH, coeffAlpha, noccAlpha, coeffBeta, noccBeta, naux, nbasis) => {
          const [bPtr, bLen] = writeF64(bMatrix);
          const [dtPtr, dtLen] = writeF64(densityTotal);
          const [chPtr, chLen] = writeF64(coreH);
          const [caPtr, caLen] = writeF64(coeffAlpha);
          const [cbPtr, cbLen] = writeF64(coeffBeta);
          const result = wasmRIFockUhf(bPtr, bLen, dtPtr, dtLen, chPtr, chLen, caPtr, caLen, noccAlpha, cbPtr, cbLen, noccBeta, naux, nbasis);
          return readF64Result(result);
        };
      }

      const wasmRIMP2 = exports.compute_ri_mp2_wasm as (
        bPtr: number, bLen: number,
        coeffPtr: number, coeffLen: number,
        epsPtr: number, epsLen: number,
        naux: number, nocc: number, nbasis: number,
      ) => number;

      if (wasmRIMP2) {
        _wasmComputeRIMP2 = (bMatrix, coefficients, epsilon, naux, nocc, nbasis) => {
          const [bPtr, bLen] = writeF64(bMatrix);
          const [coeffPtr, coeffLen] = writeF64(coefficients);
          const [epsPtr, epsLen] = writeF64(epsilon);
          return wasmRIMP2(bPtr, bLen, coeffPtr, coeffLen, epsPtr, epsLen, naux, nocc, nbasis);
        };
      }

      const wasmRIUMP2 = exports.compute_ri_ump2_wasm as (
        bPtr: number, bLen: number,
        caPtr: number, caLen: number,
        cbPtr: number, cbLen: number,
        epsAPtr: number, epsALen: number,
        epsBPtr: number, epsBLen: number,
        naux: number, noccA: number, noccB: number, nbasis: number,
      ) => number;

      if (wasmRIUMP2) {
        _wasmComputeRIUMP2 = (bMatrix, ca, cb, epsilonA, epsilonB, naux, noccA, noccB, nbasis) => {
          const [bPtr, bLen] = writeF64(bMatrix);
          const [caPtr, caLen] = writeF64(ca);
          const [cbPtr, cbLen] = writeF64(cb);
          const [epsAPtr, epsALen] = writeF64(epsilonA);
          const [epsBPtr, epsBLen] = writeF64(epsilonB);
          return wasmRIUMP2(bPtr, bLen, caPtr, caLen, cbPtr, cbLen, epsAPtr, epsALen, epsBPtr, epsBLen, naux, noccA, noccB, nbasis);
        };
      }

      // XC integration WASM
      const wasmComputeXC = exports.compute_xc_wasm as (
        shellsPtr: number, shellsLen: number,
        normsPtr: number, normsLen: number,
        daPtr: number, daLen: number,
        dbPtr: number, dbLen: number,
        gridPtr: number, gridLen: number,
        nbasis: number, funcId: number, needGrad: number,
      ) => [number, number];

      if (wasmComputeXC) {
        _wasmComputeXC = (shellsFlat, normFactors, densityA, densityB, gridFlat, nbasis, funcId, needGrad) => {
          const [shellsPtr, shellsLen] = writeF64(shellsFlat);
          const [normsPtr, normsLen] = writeF64(normFactors);
          const [daPtr, daLen] = writeF64(densityA);
          const [dbPtr, dbLen] = writeF64(densityB);
          const [gridPtr, gridLen] = writeF64(gridFlat);
          const result = wasmComputeXC(shellsPtr, shellsLen, normsPtr, normsLen,
            daPtr, daLen, dbPtr, dbLen, gridPtr, gridLen,
            nbasis, funcId, needGrad ? 1 : 0);
          return readF64Result(result);
        };
      }

      // 2-electron gradient WASM
      const wasmGrad = exports.compute_2e_gradient_wasm as (
        shellsPtr: number, shellsLen: number,
        densPtr: number, densLen: number,
        normsPtr: number, normsLen: number,
        nbasis: number, natoms: number,
      ) => [number, number];

      if (wasmGrad) {
        _wasmCompute2eGradient = (shellsFlat, density, norms, numBasis, numAtoms) => {
          const [shellsPtr, shellsLen] = writeF64(shellsFlat);
          const [densPtr, densLen] = writeF64(density);
          const [normsPtr, normsLen] = writeF64(norms);
          const result = wasmGrad(shellsPtr, shellsLen, densPtr, densLen, normsPtr, normsLen, numBasis, numAtoms);
          return readF64Result(result);
        };
      }

      // 2-electron Hessian WASM
      const wasmHess = exports.compute_2e_hessian_wasm as (
        shellsPtr: number, shellsLen: number,
        densPtr: number, densLen: number,
        normsPtr: number, normsLen: number,
        nbasis: number, natoms: number,
      ) => [number, number];

      if (wasmHess) {
        _wasmCompute2eHessian = (shellsFlat, density, norms, numBasis, numAtoms) => {
          const [shellsPtr, shellsLen] = writeF64(shellsFlat);
          const [densPtr, densLen] = writeF64(density);
          const [normsPtr, normsLen] = writeF64(norms);
          const result = wasmHess(shellsPtr, shellsLen, densPtr, densLen, normsPtr, normsLen, numBasis, numAtoms);
          return readF64Result(result);
        };
      }

      // MO-ERI transform WASM
      const wasmMOERI = exports.transform_eri_to_mo_wasm as (
        cPtr: number, cLen: number,
        eriPtr: number, eriLen: number,
        nbasis: number,
      ) => [number, number];

      if (wasmMOERI) {
        _wasmTransformERItoMO = (cFlat, eriSym, numBasis) => {
          const [cPtr, cLen] = writeF64(cFlat);
          const [eriPtr, eriLen] = writeF64(eriSym);
          const result = wasmMOERI(cPtr, cLen, eriPtr, eriLen, numBasis);
          return readF64Result(result);
        };
      }

      // CPHF solver WASM
      const wasmCPHF = exports.solve_cphf_wasm as (
        moEriPtr: number, moEriLen: number,
        epsPtr: number, epsLen: number,
        rhsPtr: number, rhsLen: number,
        nocc: number, nvir: number, nmo: number, nPert: number,
        tol: number, maxIter: number,
      ) => [number, number];

      if (wasmCPHF) {
        _wasmSolveCPHF = (moERI, eps, rhs, nocc, nvir, nmo, nPert, tol, maxIter) => {
          const [moEriPtr, moEriLen] = writeF64(moERI);
          const [epsPtr, epsLen] = writeF64(eps);
          const [rhsPtr, rhsLen] = writeF64(rhs);
          const result = wasmCPHF(moEriPtr, moEriLen, epsPtr, epsLen, rhsPtr, rhsLen, nocc, nvir, nmo, nPert, tol, maxIter);
          return readF64Result(result);
        };
      }

      // Detect if the loaded binary is the SIMD variant
      _wasmSimdAvailable = simdSupported && wasmUrl.href.includes('simd');

      return true;
    } catch (e) {
      console.warn('WASM ERI init failed, falling back to JS:', e);
      _wasmComputeEris = null;
      _wasmSimdAvailable = false;
      return false;
    }
  })();

  return _wasmInitPromise;
}

/** Check if WASM ERI is available */
export function isWasmAvailable(): boolean {
  return _wasmComputeEris !== null;
}

/** Compute ERIs using WASM. Returns null if WASM not available. */
export function computeERIsWasm(
  shells: PrimitiveShell[],
  normFactors: number[],
  numBasis: number,
  schwarzThreshold: number,
): FloatArray | null {
  if (!_wasmComputeEris) return null;

  const shellsFlat = packShells(shells);
  const normsF64 = new Float64Array(normFactors);
  return _wasmComputeEris(shellsFlat, normsF64, numBasis, schwarzThreshold);
}

/** Compute RHF Fock matrix using WASM. Returns null if WASM not available. */
export function computeFockRhfWasm(
  eriData: FloatArray, density: FloatArray,
  coreH: FloatArray, numBasis: number,
): FloatArray | null {
  if (!_wasmComputeFockRhf) return null;
  return _wasmComputeFockRhf(eriData, density, coreH, numBasis);
}

/** Compute UHF Fock matrices using WASM. Returns null if WASM not available. */
export function computeFockUhfWasm(
  eriData: FloatArray, densityAlpha: FloatArray,
  densityBeta: FloatArray, densityTotal: FloatArray,
  coreH: FloatArray, numBasis: number,
): { fockAlpha: FloatArray; fockBeta: FloatArray } | null {
  if (!_wasmComputeFockUhf) return null;
  const combined = _wasmComputeFockUhf(eriData, densityAlpha, densityBeta, densityTotal, coreH, numBasis);
  const nn = numBasis * numBasis;
  return {
    fockAlpha: combined.slice(0, nn),
    fockBeta: combined.slice(nn, 2 * nn),
  };
}

/** Compute MP2 correlation energy using WASM. Returns null if WASM not available. */
export function computeMP2EnergyWasm(
  eriData: FloatArray, coeff: FloatArray,
  epsilon: FloatArray, nocc: number, nbasis: number,
): number | null {
  if (!_wasmComputeMP2) return null;
  return _wasmComputeMP2(eriData, coeff, epsilon, nocc, nbasis);
}

/** Compute MP3 correlation energy using WASM. Returns null if WASM not available. */
export function computeMP3EnergyWasm(
  eriData: FloatArray, coeff: FloatArray,
  epsilon: FloatArray, nocc: number, nbasis: number,
): { mp2: number; mp3: number } | null {
  if (!_wasmComputeMP3) return null;
  return _wasmComputeMP3(eriData, coeff, epsilon, nocc, nbasis);
}

/** Compute CCSD correlation energy using WASM. Returns null if WASM not available. */
export function computeCCSDEnergyWasm(
  eriData: FloatArray, coeff: FloatArray,
  epsilon: FloatArray, nocc: number, nbasis: number,
): number | null {
  if (!_wasmComputeCCSD) return null;
  return _wasmComputeCCSD(eriData, coeff, epsilon, nocc, nbasis);
}

// ---------------------------------------------------------------------------
// UHF WASM wrappers
// ---------------------------------------------------------------------------

/** Compute UMP2 correlation energy using WASM. Returns null if WASM not available. */
export function computeUMP2EnergyWasm(
  eriData: FloatArray, ca: FloatArray, cb: FloatArray,
  epsilonA: FloatArray, epsilonB: FloatArray,
  noccA: number, noccB: number, nbasis: number,
): number | null {
  if (!_wasmComputeUMP2) return null;
  return _wasmComputeUMP2(eriData, ca, cb, epsilonA, epsilonB, noccA, noccB, nbasis);
}

/** Compute UMP3 correlation energy using WASM. Returns null if WASM not available. */
export function computeUMP3EnergyWasm(
  eriData: FloatArray, ca: FloatArray, cb: FloatArray,
  epsilonA: FloatArray, epsilonB: FloatArray,
  noccA: number, noccB: number, nbasis: number,
): { mp2: number; mp3: number } | null {
  if (!_wasmComputeUMP3) return null;
  return _wasmComputeUMP3(eriData, ca, cb, epsilonA, epsilonB, noccA, noccB, nbasis);
}

/** Compute UCCSD correlation energy using WASM. Returns null if WASM not available. */
export function computeUCCSDEnergyWasm(
  eriData: FloatArray, ca: FloatArray, cb: FloatArray,
  epsilonA: FloatArray, epsilonB: FloatArray,
  noccA: number, noccB: number, nbasis: number,
): number | null {
  if (!_wasmComputeUCCSD) return null;
  return _wasmComputeUCCSD(eriData, ca, cb, epsilonA, epsilonB, noccA, noccB, nbasis);
}

// ---------------------------------------------------------------------------
// ROHF WASM wrappers (semicanonical → UHF WASM)
// ---------------------------------------------------------------------------

import { Matrix } from '../linalg/matrix';
import { jacobiEigen } from '../linalg/eigendecomposition';
import { matmul, matmulAtB } from '../linalg/matmul';

/** Build semicanonical orbitals from ROHF MO coefficients and Fock matrix.
 *  Block-diagonalizes F_MO within occ [0,nocc) and vir [nocc,N) blocks. */
function semicanonicalize(
  C: Matrix, FMO: Matrix, nocc: number, nbasis: number,
): { Csc: Matrix; epsilon: Float64Array } {
  const N = nbasis;
  const nvir = N - nocc;
  const U = Matrix.identity(N);
  const eps = new Float64Array(N);

  // Occupied block
  if (nocc > 0) {
    const occ = new Matrix(nocc, nocc);
    for (let i = 0; i < nocc; i++)
      for (let j = 0; j < nocc; j++)
        occ.set(i, j, FMO.get(i, j));
    const eig = jacobiEigen(occ);
    for (let i = 0; i < nocc; i++) {
      eps[i] = eig.eigenvalues[i];
      for (let j = 0; j < nocc; j++)
        U.set(i, j, eig.eigenvectors.get(i, j));
    }
  }

  // Virtual block
  if (nvir > 0) {
    const vir = new Matrix(nvir, nvir);
    for (let i = 0; i < nvir; i++)
      for (let j = 0; j < nvir; j++)
        vir.set(i, j, FMO.get(nocc + i, nocc + j));
    const eig = jacobiEigen(vir);
    for (let i = 0; i < nvir; i++) {
      eps[nocc + i] = eig.eigenvalues[i];
      for (let j = 0; j < nvir; j++)
        U.set(nocc + i, nocc + j, eig.eigenvectors.get(i, j));
    }
  }

  return { Csc: matmul(C, U), epsilon: eps };
}

/** Compute ROMP2 energy using WASM (semicanonical → UMP2 WASM).
 *  Returns null if WASM not available. */
export function computeROMP2EnergyWasm(
  eriData: FloatArray, coeff: Matrix,
  Fa: Matrix, Fb: Matrix,
  noccA: number, noccB: number, nbasis: number,
): number | null {
  if (!_wasmComputeUMP2) return null;
  const FaMO = matmul(matmulAtB(coeff, Fa), coeff);
  const FbMO = matmul(matmulAtB(coeff, Fb), coeff);
  const { Csc: CaSc, epsilon: epsA } = semicanonicalize(coeff, FaMO, noccA, nbasis);
  const { Csc: CbSc, epsilon: epsB } = semicanonicalize(coeff, FbMO, noccB, nbasis);
  return _wasmComputeUMP2(eriData, CaSc.data, CbSc.data, epsA, epsB, noccA, noccB, nbasis);
}

/** Compute ROMP3 energy using WASM (semicanonical → UMP3 WASM).
 *  Returns null if WASM not available. */
export function computeROMP3EnergyWasm(
  eriData: FloatArray, coeff: Matrix,
  Fa: Matrix, Fb: Matrix,
  noccA: number, noccB: number, nbasis: number,
): { mp2: number; mp3: number } | null {
  if (!_wasmComputeUMP3) return null;
  const FaMO = matmul(matmulAtB(coeff, Fa), coeff);
  const FbMO = matmul(matmulAtB(coeff, Fb), coeff);
  const { Csc: CaSc, epsilon: epsA } = semicanonicalize(coeff, FaMO, noccA, nbasis);
  const { Csc: CbSc, epsilon: epsB } = semicanonicalize(coeff, FbMO, noccB, nbasis);
  return _wasmComputeUMP3(eriData, CaSc.data, CbSc.data, epsA, epsB, noccA, noccB, nbasis);
}

/** Compute ROCCSD energy using WASM (semicanonical → UCCSD WASM).
 *  Returns null if WASM not available. */
export function computeROCCSDEnergyWasm(
  eriData: FloatArray, coeff: Matrix,
  Fa: Matrix, Fb: Matrix,
  noccA: number, noccB: number, nbasis: number,
): number | null {
  if (!_wasmComputeUCCSD) return null;
  const FaMO = matmul(matmulAtB(coeff, Fa), coeff);
  const FbMO = matmul(matmulAtB(coeff, Fb), coeff);
  const { Csc: CaSc, epsilon: epsA } = semicanonicalize(coeff, FaMO, noccA, nbasis);
  const { Csc: CbSc, epsilon: epsB } = semicanonicalize(coeff, FbMO, noccB, nbasis);
  return _wasmComputeUCCSD(eriData, CaSc.data, CbSc.data, epsA, epsB, noccA, noccB, nbasis);
}

// ---------------------------------------------------------------------------
// RI WASM wrappers
// ---------------------------------------------------------------------------

/** Check if RI WASM functions are available */
export function isRIWasmAvailable(): boolean {
  return _wasmComputeRIFockRhf !== null;
}

/** Compute RI Setup (2c/3c integrals, Cholesky, B matrix) using WASM.
 *  Returns { naux, bMatrix } or null if WASM not available. */
export function computeRISetupWasm(
  primaryShellsFlat: FloatArray, primaryNorm: FloatArray, nbasis: number,
  auxShellsFlat: FloatArray, auxNorm: FloatArray, naux: number,
): { naux: number; bMatrix: Float64Array } | null {
  if (!_wasmComputeRISetup) return null;
  const result = _wasmComputeRISetup(primaryShellsFlat, primaryNorm, nbasis, auxShellsFlat, auxNorm, naux);
  const nauxEff = result[0];
  const bMatrix = result.slice(1);
  return { naux: nauxEff, bMatrix };
}

/** Compute RI-RHF Fock matrix using WASM. Returns null if WASM not available. */
export function computeRIFockRhfWasm(
  bMatrix: FloatArray, density: FloatArray,
  coreH: FloatArray, coefficients: FloatArray,
  naux: number, nbasis: number, nocc: number,
): Float64Array | null {
  if (!_wasmComputeRIFockRhf) return null;
  return _wasmComputeRIFockRhf(bMatrix, density, coreH, coefficients, naux, nbasis, nocc);
}

/** Compute RI-UHF Fock matrices using WASM. Returns null if WASM not available. */
export function computeRIFockUhfWasm(
  bMatrix: FloatArray, densityTotal: FloatArray,
  coreH: FloatArray,
  coeffAlpha: FloatArray, noccAlpha: number,
  coeffBeta: FloatArray, noccBeta: number,
  naux: number, nbasis: number,
): { fockAlpha: Float64Array; fockBeta: Float64Array } | null {
  if (!_wasmComputeRIFockUhf) return null;
  const combined = _wasmComputeRIFockUhf(bMatrix, densityTotal, coreH, coeffAlpha, noccAlpha, coeffBeta, noccBeta, naux, nbasis);
  const nn = nbasis * nbasis;
  return {
    fockAlpha: combined.slice(0, nn),
    fockBeta: combined.slice(nn, 2 * nn),
  };
}

/** Compute RI-MP2 energy using WASM. Returns null if WASM not available. */
export function computeRIMP2EnergyWasm(
  bMatrix: FloatArray, coefficients: FloatArray,
  epsilon: FloatArray,
  naux: number, nocc: number, nbasis: number,
): number | null {
  if (!_wasmComputeRIMP2) return null;
  return _wasmComputeRIMP2(bMatrix, coefficients, epsilon, naux, nocc, nbasis);
}

/** Compute RI-UMP2 energy using WASM. Returns null if WASM not available. */
export function computeRIUMP2EnergyWasm(
  bMatrix: FloatArray, ca: FloatArray, cb: FloatArray,
  epsilonA: FloatArray, epsilonB: FloatArray,
  naux: number, noccA: number, noccB: number, nbasis: number,
): number | null {
  if (!_wasmComputeRIUMP2) return null;
  return _wasmComputeRIUMP2(bMatrix, ca, cb, epsilonA, epsilonB, naux, noccA, noccB, nbasis);
}

// ---------------------------------------------------------------------------
// XC integration WASM wrapper
// ---------------------------------------------------------------------------

export interface XCWasmResult {
  exc: number;
  vxcA: Float64Array;
  vxcB?: Float64Array;
  numElectrons: number;
}

/** Check if XC WASM is available */
export function isXCWasmAvailable(): boolean {
  return _wasmComputeXC !== null;
}

/** Compute XC energy and V_xc matrix using WASM. Returns null if unavailable. */
export function computeXCWasm(
  shellsFlat: Float64Array, normFactors: Float64Array,
  densityA: Float64Array, densityB: Float64Array,
  gridFlat: Float64Array, nbasis: number,
  funcId: number, needGrad: boolean,
): XCWasmResult | null {
  if (!_wasmComputeXC) return null;
  const result = _wasmComputeXC(shellsFlat, normFactors, densityA, densityB, gridFlat, nbasis, funcId, needGrad);
  const nn = nbasis * nbasis;
  const unrestricted = densityB.length > 0;
  return {
    exc: result[0],
    numElectrons: result[1],
    vxcA: result.slice(2, 2 + nn),
    vxcB: unrestricted ? result.slice(2 + nn, 2 + 2 * nn) : undefined,
  };
}

// ---------------------------------------------------------------------------
// Gradient / Hessian / CPHF WASM wrappers
// ---------------------------------------------------------------------------

/** Check if 2-electron gradient WASM is available */
export function is2eGradientWasmAvailable(): boolean {
  return _wasmCompute2eGradient !== null;
}

/** Compute 2-electron gradient using WASM. Returns null if unavailable. */
export function compute2eGradientWasm(
  shellsFlat: Float64Array, density: Float64Array,
  norms: Float64Array, numBasis: number, numAtoms: number,
): Float64Array | null {
  if (!_wasmCompute2eGradient) return null;
  return _wasmCompute2eGradient(shellsFlat, density, norms, numBasis, numAtoms);
}

/** Check if 2-electron Hessian WASM is available */
export function is2eHessianWasmAvailable(): boolean {
  return _wasmCompute2eHessian !== null;
}

/** Compute 2-electron skeleton Hessian using WASM. Returns null if unavailable. */
export function compute2eHessianWasm(
  shellsFlat: Float64Array, density: Float64Array,
  norms: Float64Array, numBasis: number, numAtoms: number,
): Float64Array | null {
  if (!_wasmCompute2eHessian) return null;
  return _wasmCompute2eHessian(shellsFlat, density, norms, numBasis, numAtoms);
}

/** Check if MO-ERI transform WASM is available */
export function isTransformERItoMOWasmAvailable(): boolean {
  return _wasmTransformERItoMO !== null;
}

/** Transform AO ERIs to MO basis using WASM. Returns null if unavailable. */
export function transformERItoMOWasm(
  cFlat: Float64Array, eriSym: Float64Array, numBasis: number,
): Float64Array | null {
  if (!_wasmTransformERItoMO) return null;
  return _wasmTransformERItoMO(cFlat, eriSym, numBasis);
}

/** Check if CPHF solver WASM is available */
export function isCPHFWasmAvailable(): boolean {
  return _wasmSolveCPHF !== null;
}

/** Solve CPHF equations using WASM. Returns null if unavailable. */
export function solveCPHFWasm(
  moERI: Float64Array, eps: Float64Array, rhs: Float64Array,
  nocc: number, nvir: number, nmo: number, nPert: number,
  tol: number, maxIter: number,
): Float64Array | null {
  if (!_wasmSolveCPHF) return null;
  return _wasmSolveCPHF(moERI, eps, rhs, nocc, nvir, nmo, nPert, tol, maxIter);
}
