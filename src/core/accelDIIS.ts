/** DIIS accelerator — wraps existing DIIS class, one instance per channel */

import { DIIS } from './diis';
import type { SCFAccelerator, SCFAccelInput, SCFAccelOutput } from './scfAccelerator';
import { copyToFloatArray } from '../linalg/matrix';

export class AccelDIIS implements SCFAccelerator {
  readonly name = 'DIIS';
  private readonly maxHist: number;
  private instances: DIIS[] = [];

  constructor(maxHist?: number) {
    this.maxHist = maxHist ?? 8;
  }

  accelerate(input: SCFAccelInput): SCFAccelOutput {
    const nCh = input.fock.length;
    // Lazily create one DIIS instance per channel
    while (this.instances.length < nCh) {
      this.instances.push(new DIIS(this.maxHist, 2));
    }

    const fockOut: typeof input.fock = [];
    for (let ch = 0; ch < nCh; ch++) {
      const diis = this.instances[ch];
      diis.push(copyToFloatArray(input.fock[ch]), copyToFloatArray(input.error[ch]));
      if (diis.canExtrapolate()) {
        fockOut.push(diis.extrapolate());
      } else {
        fockOut.push(copyToFloatArray(input.fock[ch]));
      }
    }
    return { fock: fockOut };
  }

  reset() {
    for (const d of this.instances) d.reset();
    this.instances.length = 0;
  }
}
