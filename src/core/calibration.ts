/** Benchmark-based backend calibration.
 *
 *  bench.html runs ERI/Fock benchmarks and saves crossover thresholds
 *  to localStorage.  The auto backend selector reads them here instead
 *  of using hardcoded constants. */

const STORAGE_KEY = 'gansu-calibration';

/** Default thresholds (used when no calibration data exists) */
const DEFAULT_WASM_POSTHF_THRESHOLD = 5;

export interface CalibrationData {
  /** ISO timestamp of calibration run */
  timestamp: string;
  /** Whether WASM was available during calibration */
  wasmAvailable: boolean;
  /** numBasis above which WASM post-HF (MP2/MP3/CCSD) is faster than JS */
  wasmPostHFThreshold?: number;
}

let _cached: CalibrationData | null | undefined;

/** Inject calibration data directly (for Web Worker context where
 *  localStorage is unavailable).  Call before any SCF computation. */
export function setCalibration(data: CalibrationData): void {
  _cached = data;
}

/** Read calibration from localStorage (cached after first call).
 *  Returns null in Worker context or when no data saved. */
export function getCalibration(): CalibrationData | null {
  if (_cached !== undefined) return _cached;
  try {
    if (typeof localStorage === 'undefined') { _cached = null; return null; }
    const raw = localStorage.getItem(STORAGE_KEY);
    if (!raw) { _cached = null; return null; }
    const data = JSON.parse(raw) as CalibrationData;
    if (typeof data.wasmAvailable !== 'boolean') {
      _cached = null;
      return null;
    }
    _cached = data;
    return data;
  } catch {
    _cached = null;
    return null;
  }
}

/** Save calibration data to localStorage */
export function saveCalibration(data: CalibrationData): void {
  localStorage.setItem(STORAGE_KEY, JSON.stringify(data));
  _cached = data;
}

/** Clear cached calibration (for testing) */
export function clearCalibrationCache(): void {
  _cached = undefined;
}

/** WASM threshold for post-HF computation (calibrated or default) */
export function getWasmPostHFThreshold(): number {
  const cal = getCalibration();
  return cal?.wasmPostHFThreshold ?? DEFAULT_WASM_POSTHF_THRESHOLD;
}

/** Determine WASM threshold from benchmark results. */
export function computeThresholds(
  results: Array<{
    numBasis: number;
    mp2JS?: number; mp2WASM?: number;
  }>,
): { wasmPostHFThreshold: number } {
  const sorted = [...results].sort((a, b) => a.numBasis - b.numBasis);

  function findCrossover<T extends { numBasis: number }>(
    entries: T[],
    getSlow: (r: T) => number,
    getFast: (r: T) => number,
  ): number {
    for (let i = 0; i < entries.length; i++) {
      const s = getSlow(entries[i]);
      const f = getFast(entries[i]);
      if (s <= 0 || f <= 0) continue;
      if (f < s) {
        let consistent = true;
        for (let j = i + 1; j < entries.length; j++) {
          const sj = getSlow(entries[j]);
          const fj = getFast(entries[j]);
          if (sj > 0 && fj > 0 && fj > sj * 1.1) {
            consistent = false;
            break;
          }
        }
        if (consistent) return entries[i].numBasis;
      }
    }
    return Infinity;
  }

  return {
    wasmPostHFThreshold: findCrossover(
      sorted.filter(r => (r.mp2JS ?? 0) > 0 && (r.mp2WASM ?? 0) > 0),
      r => r.mp2JS ?? 0, r => r.mp2WASM ?? 0,
    ),
  };
}
