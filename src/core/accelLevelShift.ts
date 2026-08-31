/** Level Shifting accelerator: F' = F + b*(S - S*P*S) */

import type { SCFAccelerator, SCFAccelInput, SCFAccelOutput } from './scfAccelerator';
import { copyToFloatArray, type FloatArray } from '../linalg/matrix';

export class AccelLevelShift implements SCFAccelerator {
  readonly name = 'Level Shift';
  private readonly shift: number;

  constructor(shift?: number) {
    this.shift = shift ?? 0.5;
  }

  accelerate(input: SCFAccelInput): SCFAccelOutput {
    const n = input.numBasis;
    const S = input.overlap;
    const b = this.shift;
    const nCh = input.fock.length;
    const fockOut: FloatArray[] = [];

    for (let ch = 0; ch < nCh; ch++) {
      const F = input.fock[ch];
      const P = input.density[ch]; // idempotent density (P/2 for RHF, P for UHF)
      // Compute S*P*S
      const SP = new Float64Array(n * n);
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          let sum = 0;
          for (let k = 0; k < n; k++) sum += S[i * n + k] * P[k * n + j];
          SP[i * n + j] = sum;
        }
      }
      const SPS = new Float64Array(n * n);
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          let sum = 0;
          for (let k = 0; k < n; k++) sum += SP[i * n + k] * S[k * n + j];
          SPS[i * n + j] = sum;
        }
      }
      // F' = F + b*(S - S*P*S)
      const shifted = copyToFloatArray(F);
      for (let k = 0; k < n * n; k++) {
        shifted[k] += b * (S[k] - SPS[k]);
      }
      fockOut.push(shifted);
    }
    return { fock: fockOut };
  }

  reset() { /* stateless */ }
}
