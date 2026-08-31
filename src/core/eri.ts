/** ERI storage with 8-fold symmetry index */

import { eri1DIndex } from './integrals2e';
import type { FloatArray } from '../linalg/matrix';

export class ERIStored {
  readonly data: FloatArray;
  private readonly numBasis: number;

  constructor(data: FloatArray, numBasis: number) {
    this.data = data;
    this.numBasis = numBasis;
  }

  /** Get (ij|kl) integral value */
  get(i: number, j: number, k: number, l: number): number {
    return this.data[eri1DIndex(i, j, k, l, this.numBasis)];
  }
}
