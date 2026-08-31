/** Damping accelerator: F = α*F_new + (1-α)*F_old */

import type { SCFAccelerator, SCFAccelInput, SCFAccelOutput } from './scfAccelerator';
import { copyToFloatArray, type FloatArray } from '../linalg/matrix';

export class AccelDamping implements SCFAccelerator {
  readonly name = 'Damping';
  private readonly alpha: number;
  private prevFock: FloatArray[] | null = null;

  constructor(alpha?: number) {
    this.alpha = alpha ?? 0.5;
  }

  accelerate(input: SCFAccelInput): SCFAccelOutput {
    const nCh = input.fock.length;
    const fockOut: FloatArray[] = [];

    for (let ch = 0; ch < nCh; ch++) {
      const fNew = input.fock[ch];
      if (this.prevFock && this.prevFock[ch]) {
        const fOld = this.prevFock[ch];
        const mixed = copyToFloatArray(fNew);
        for (let k = 0; k < mixed.length; k++) {
          mixed[k] = this.alpha * fNew[k] + (1 - this.alpha) * fOld[k];
        }
        fockOut.push(mixed);
      } else {
        fockOut.push(copyToFloatArray(fNew));
      }
    }
    this.prevFock = fockOut.map(f => copyToFloatArray(f));
    return { fock: fockOut };
  }

  reset() {
    this.prevFock = null;
  }
}
