/** DIIS (Direct Inversion in the Iterative Subspace) — port of GANSU diis.hpp */

import { solveLinearSystem } from '../linalg/solve';
import { copyToFloatArray, type FloatArray } from '../linalg/matrix';

export class DIIS {
  private readonly maxSize: number;
  private readonly minSize: number;
  private readonly fockHistory: FloatArray[] = [];
  private readonly errorHistory: FloatArray[] = [];

  constructor(maxSize = 8, minSize = 2) {
    this.maxSize = maxSize;
    this.minSize = minSize;
  }

  push(fock: FloatArray, error: FloatArray) {
    this.fockHistory.push(copyToFloatArray(fock));
    this.errorHistory.push(copyToFloatArray(error));

    if (this.fockHistory.length > this.maxSize) {
      this.fockHistory.shift();
      this.errorHistory.shift();
    }
  }

  canExtrapolate(): boolean {
    return this.fockHistory.length >= this.minSize;
  }

  /** Extrapolate a new Fock matrix from the history */
  extrapolate(): FloatArray {
    const m = this.fockHistory.length;
    const n = m + 1; // B matrix is (m+1) x (m+1) with Lagrange constraint
    const vecSize = this.fockHistory[0].length;

    // Build B matrix: B[i][j] = <e_i | e_j>, plus Lagrange multiplier row/col
    // B-matrix solve stays Float64Array for numerical stability (small matrix)
    const B = new Float64Array(n * n);
    const rhs = new Float64Array(n);

    for (let i = 0; i < m; i++) {
      for (let j = i; j < m; j++) {
        let dot = 0;
        for (let k = 0; k < this.errorHistory[i].length; k++) {
          dot += this.errorHistory[i][k] * this.errorHistory[j][k];
        }
        B[i * n + j] = dot;
        B[j * n + i] = dot;
      }
      // Lagrange constraint
      B[i * n + m] = -1;
      B[m * n + i] = -1;
    }
    B[m * n + m] = 0;
    rhs[m] = -1;

    // Solve B * c = rhs
    const c = solveLinearSystem(B, rhs, n);

    // New Fock = sum c_i * F_i
    const newFock = copyToFloatArray(new Float64Array(vecSize));
    for (let i = 0; i < m; i++) {
      for (let k = 0; k < vecSize; k++) {
        newFock[k] += c[i] * this.fockHistory[i][k];
      }
    }

    return newFock;
  }

  reset() {
    this.fockHistory.length = 0;
    this.errorHistory.length = 0;
  }
}
