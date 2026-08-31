/** EDIIS (Energy-DIIS) accelerator — minimizes energy on simplex via projected gradient */

import type { SCFAccelerator, SCFAccelInput, SCFAccelOutput } from './scfAccelerator';
import { copyToFloatArray, type FloatArray } from '../linalg/matrix';

export class AccelEDIIS implements SCFAccelerator {
  readonly name = 'EDIIS';
  private readonly maxHist: number;
  private fockHist: FloatArray[][] = [];

  constructor(maxHist?: number) {
    this.maxHist = maxHist ?? 8;
  }    // [iteration][channel]
  private densHist: FloatArray[][] = [];    // [iteration][channel]
  private energyHist: number[] = [];

  accelerate(input: SCFAccelInput): SCFAccelOutput {
    const nCh = input.fock.length;
    const n2 = input.numBasis * input.numBasis;

    // Store current state
    this.fockHist.push(input.fock.map(f => copyToFloatArray(f)));
    this.densHist.push(input.density.map(d => copyToFloatArray(d)));
    this.energyHist.push(input.energy);

    // Trim history
    if (this.fockHist.length > this.maxHist) {
      this.fockHist.shift();
      this.densHist.shift();
      this.energyHist.shift();
    }

    const m = this.fockHist.length;
    if (m < 2) return { fock: input.fock.map(f => copyToFloatArray(f)) };

    // Compute EDIIS objective gradient terms:
    // E(c) = Σ_i c_i * E_i + Σ_{ij} c_i * c_j * Tr((D_i - D_j)(F_i - F_j))
    // We solve for c on the simplex Σc_i=1, c_i>=0 via projected gradient.

    // Build quadratic matrix B[i][j] = Σ_ch Tr((D_i-D_j)(F_i-F_j))
    const B = new Float64Array(m * m);
    for (let i = 0; i < m; i++) {
      for (let j = i; j < m; j++) {
        let val = 0;
        for (let ch = 0; ch < nCh; ch++) {
          const Di = this.densHist[i][ch];
          const Dj = this.densHist[j][ch];
          const Fi = this.fockHist[i][ch];
          const Fj = this.fockHist[j][ch];
          for (let k = 0; k < n2; k++) {
            val += (Di[k] - Dj[k]) * (Fi[k] - Fj[k]);
          }
        }
        B[i * m + j] = val;
        B[j * m + i] = val;
      }
    }

    // Projected gradient on simplex
    const c = new Float64Array(m);
    c.fill(1.0 / m);

    for (let iter = 0; iter < 200; iter++) {
      // Compute gradient: g_i = E_i + 2*Σ_j B[i][j]*c_j
      const g = new Float64Array(m);
      for (let i = 0; i < m; i++) {
        g[i] = this.energyHist[i];
        for (let j = 0; j < m; j++) {
          g[i] += 2 * B[i * m + j] * c[j];
        }
      }

      // Project gradient onto simplex tangent plane: g - (Σg_i/m)
      let gMean = 0;
      for (let i = 0; i < m; i++) gMean += g[i];
      gMean /= m;
      for (let i = 0; i < m; i++) g[i] -= gMean;

      // Step
      const step = 0.01;
      for (let i = 0; i < m; i++) c[i] -= step * g[i];

      // Project back to simplex (clip + normalize)
      for (let i = 0; i < m; i++) c[i] = Math.max(c[i], 0);
      let cSum = 0;
      for (let i = 0; i < m; i++) cSum += c[i];
      if (cSum < 1e-15) { c.fill(1.0 / m); break; }
      for (let i = 0; i < m; i++) c[i] /= cSum;
    }

    // Combine: F = Σ c_i * F_i
    const fockOut: FloatArray[] = [];
    for (let ch = 0; ch < nCh; ch++) {
      const combined = copyToFloatArray(new Float64Array(n2));
      for (let i = 0; i < m; i++) {
        const w = c[i];
        if (Math.abs(w) < 1e-15) continue;
        const fi = this.fockHist[i][ch];
        for (let k = 0; k < n2; k++) combined[k] += w * fi[k];
      }
      fockOut.push(combined);
    }
    return { fock: fockOut };
  }

  reset() {
    this.fockHist.length = 0;
    this.densHist.length = 0;
    this.energyHist.length = 0;
  }
}
