/** ADIIS+DIIS hybrid — ADIIS while error is large, switches to DIIS when converging */

import type { SCFAccelerator, SCFAccelInput, SCFAccelOutput } from './scfAccelerator';
import { copyToFloatArray, type FloatArray } from '../linalg/matrix';
import { DIIS } from './diis';

export class AccelADIIS implements SCFAccelerator {
  readonly name = 'ADIIS+DIIS';
  private readonly maxHist: number;
  private readonly switchThreshold: number;
  private fockHist: FloatArray[][] = [];
  private densHist: FloatArray[][] = [];
  private energyHist: number[] = [];
  private diisInstances: DIIS[] = [];

  constructor(maxHist?: number, switchThreshold?: number) {
    this.maxHist = maxHist ?? 8;
    this.switchThreshold = switchThreshold ?? 0.1;
  }

  accelerate(input: SCFAccelInput): SCFAccelOutput {
    const nCh = input.fock.length;
    const n2 = input.numBasis * input.numBasis;

    // Ensure DIIS instances exist
    while (this.diisInstances.length < nCh) this.diisInstances.push(new DIIS(this.maxHist, 2));

    // Feed DIIS
    for (let ch = 0; ch < nCh; ch++) {
      this.diisInstances[ch].push(copyToFloatArray(input.fock[ch]), copyToFloatArray(input.error[ch]));
    }

    // Store ADIIS history
    this.fockHist.push(input.fock.map(f => copyToFloatArray(f)));
    this.densHist.push(input.density.map(d => copyToFloatArray(d)));
    this.energyHist.push(input.energy);
    if (this.fockHist.length > this.maxHist) {
      this.fockHist.shift();
      this.densHist.shift();
      this.energyHist.shift();
    }

    // Determine max error norm
    let maxErr = 0;
    for (let ch = 0; ch < nCh; ch++) {
      let norm = 0;
      for (let k = 0; k < input.error[ch].length; k++) {
        norm += input.error[ch][k] * input.error[ch][k];
      }
      maxErr = Math.max(maxErr, Math.sqrt(norm));
    }

    // Use DIIS if error is small enough
    if (maxErr < this.switchThreshold && this.diisInstances[0].canExtrapolate()) {
      const fockOut: FloatArray[] = [];
      for (let ch = 0; ch < nCh; ch++) {
        fockOut.push(this.diisInstances[ch].extrapolate());
      }
      return { fock: fockOut };
    }

    // ADIIS: minimize Tr(F*ΔD) based energy estimate on simplex
    const m = this.fockHist.length;
    if (m < 2) return { fock: input.fock.map(f => copyToFloatArray(f)) };

    // Reference = last entry
    const ref = m - 1;
    // ΔD_i = D_i - D_ref, ΔF_i = F_i - F_ref
    // E(c) ≈ E_ref + Σ c_i Tr(F_ref * ΔD_i) + Σ c_i c_j Tr(ΔF_i * ΔD_j)

    // Linear term
    const lin = new Float64Array(m);
    for (let i = 0; i < m; i++) {
      for (let ch = 0; ch < nCh; ch++) {
        const Fref = this.fockHist[ref][ch];
        const Di = this.densHist[i][ch];
        const Dref = this.densHist[ref][ch];
        for (let k = 0; k < n2; k++) {
          lin[i] += Fref[k] * (Di[k] - Dref[k]);
        }
      }
    }

    // Quadratic term
    const quad = new Float64Array(m * m);
    for (let i = 0; i < m; i++) {
      for (let j = i; j < m; j++) {
        let val = 0;
        for (let ch = 0; ch < nCh; ch++) {
          const Fi = this.fockHist[i][ch];
          const Fref = this.fockHist[ref][ch];
          const Dj = this.densHist[j][ch];
          const Dref = this.densHist[ref][ch];
          for (let k = 0; k < n2; k++) {
            val += (Fi[k] - Fref[k]) * (Dj[k] - Dref[k]);
          }
        }
        quad[i * m + j] = val;
        quad[j * m + i] = val;
      }
    }

    // Projected gradient on simplex
    const c = new Float64Array(m);
    c.fill(1.0 / m);

    for (let iter = 0; iter < 200; iter++) {
      const g = new Float64Array(m);
      for (let i = 0; i < m; i++) {
        g[i] = lin[i];
        for (let j = 0; j < m; j++) g[i] += 2 * quad[i * m + j] * c[j];
      }
      let gMean = 0;
      for (let i = 0; i < m; i++) gMean += g[i];
      gMean /= m;
      for (let i = 0; i < m; i++) g[i] -= gMean;

      const step = 0.01;
      for (let i = 0; i < m; i++) c[i] -= step * g[i];

      for (let i = 0; i < m; i++) c[i] = Math.max(c[i], 0);
      let cSum = 0;
      for (let i = 0; i < m; i++) cSum += c[i];
      if (cSum < 1e-15) { c.fill(1.0 / m); break; }
      for (let i = 0; i < m; i++) c[i] /= cSum;
    }

    // Combine Fock
    const fockOut: FloatArray[] = [];
    for (let ch = 0; ch < nCh; ch++) {
      const combined = copyToFloatArray(new Float64Array(n2));
      for (let i = 0; i < m; i++) {
        if (Math.abs(c[i]) < 1e-15) continue;
        const fi = this.fockHist[i][ch];
        for (let k = 0; k < n2; k++) combined[k] += c[i] * fi[k];
      }
      fockOut.push(combined);
    }
    return { fock: fockOut };
  }

  reset() {
    this.fockHist.length = 0;
    this.densHist.length = 0;
    this.energyHist.length = 0;
    for (const d of this.diisInstances) d.reset();
    this.diisInstances.length = 0;
  }
}
