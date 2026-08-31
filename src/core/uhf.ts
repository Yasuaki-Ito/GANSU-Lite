/** UHF/UKS (Unrestricted Hartree-Fock / Unrestricted Kohn-Sham) */

import { HF } from './hf';
import type { Atom, PrimitiveShell, ShellTypeInfo, BasisRange } from './types';
import { Matrix, type FloatArray, copyToFloatArray } from '../linalg/matrix';
import { jacobiEigen } from '../linalg/eigendecomposition';
import { matmul, matmulAtB, traceProduct } from '../linalg/matmul';
import { computeERIs } from './integrals2e';
import type { EriBackend } from './hf';
import { ERIStored } from './eri';
import { isWasmAvailable, computeFockUhfWasm } from './eriWasm';
import type { RIData } from './ri';
import { computeXC } from './xcIntegration';

export class UHF extends HF {
  private fockAlpha!: Matrix;
  private fockBeta!: Matrix;
  private densityAlpha!: Matrix;
  private densityBeta!: Matrix;
  private densityTotal!: Matrix;
  private coeffAlpha!: Matrix;
  private coeffBeta!: Matrix;
  private _orbitalEnergiesAlpha!: FloatArray;
  private _orbitalEnergiesBeta!: FloatArray;
  private eri!: ERIStored;
  private _riData?: RIData;
  private energy = 0;

  async precomputeERI(onProgress?: (percent: number, msg?: string) => void, backend?: EriBackend) {
    if (this._riData) {
      if (onProgress) onProgress(100, 'RI');
      return;
    }

    const eriData = await computeERIs(
      this.primitiveShells, this.normFactors, this.numBasis, 1e-10, onProgress, backend,
    );
    this.eri = new ERIStored(eriData, this.numBasis);
  }

  setRIData(ri: RIData) {
    this._riData = ri;
  }

  get riData(): RIData | undefined { return this._riData; }

  private initialAlphaDensity?: Matrix;
  private initialBetaDensity?: Matrix;

  /** Set separate alpha/beta initial density (broken-symmetry guess) */
  setInitialDensityGuessAlphaBeta(Pa: Matrix, Pb: Matrix) {
    this.initialAlphaDensity = Pa;
    this.initialBetaDensity = Pb;
    // Set total density so parent solve() triggers the density-guess branch
    const Pt = new Matrix(Pa.rows, Pa.cols);
    for (let k = 0; k < Pa.data.length; k++) Pt.data[k] = Pa.data[k] + Pb.data[k];
    this.setInitialDensityGuess(Pt);
  }

  protected setInitialDensity(P: Matrix) {
    if (this.initialAlphaDensity && this.initialBetaDensity) {
      this.densityAlpha = this.initialAlphaDensity;
      this.densityBeta = this.initialBetaDensity;
      const n = this.numBasis;
      const Pt = new Matrix(n, n);
      for (let k = 0; k < n * n; k++) Pt.data[k] = this.initialAlphaDensity.data[k] + this.initialBetaDensity.data[k];
      this.densityTotal = Pt;
      this.initialAlphaDensity = undefined;
      this.initialBetaDensity = undefined;
      return;
    }
    const n = this.numBasis;
    const Pa = new Matrix(n, n);
    const Pb = new Matrix(n, n);
    for (let k = 0; k < n * n; k++) {
      Pa.data[k] = P.data[k] * 0.5;
      Pb.data[k] = P.data[k] * 0.5;
    }
    this.densityAlpha = Pa;
    this.densityBeta = Pb;
    this.densityTotal = P.clone();
  }

  /** Initial guess: F_alpha = F_beta = H_core */
  guessInitialFockMatrix() {
    this.fockAlpha = this.coreHamiltonianMatrix.clone();
    this.fockBeta = this.coreHamiltonianMatrix.clone();
  }

  guessInitialFockMatrixGWH() {
    const gwh = this.buildGWHFock();
    this.fockAlpha = gwh;
    this.fockBeta = gwh.clone();
  }

  /** Diagonalize Fock matrices in orthogonal basis for alpha and beta separately */
  computeCoefficientMatrix() {
    // Alpha
    const FpA = matmul(matmulAtB(this.transformMatrix, this.fockAlpha), this.transformMatrix);
    const eigA = jacobiEigen(FpA);
    this.coeffAlpha = matmul(this.transformMatrix, eigA.eigenvectors);
    this._orbitalEnergiesAlpha = eigA.eigenvalues;

    // Beta
    const FpB = matmul(matmulAtB(this.transformMatrix, this.fockBeta), this.transformMatrix);
    const eigB = jacobiEigen(FpB);
    this.coeffBeta = matmul(this.transformMatrix, eigB.eigenvectors);
    this._orbitalEnergiesBeta = eigB.eigenvalues;
  }

  /** P_alpha = C_alpha_occ * C_alpha_occ^T (no factor of 2)
   *  P_beta  = C_beta_occ  * C_beta_occ^T
   *  P_total = P_alpha + P_beta */
  computeDensityMatrix() {
    const n = this.numBasis;

    // Alpha density
    const Pa = new Matrix(n, n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        let sum = 0;
        for (let a = 0; a < this.numAlphaSpins; a++) {
          sum += this.coeffAlpha.get(i, a) * this.coeffAlpha.get(j, a);
        }
        Pa.set(i, j, sum);
      }
    }
    this.densityAlpha = Pa;

    // Beta density
    const Pb = new Matrix(n, n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        let sum = 0;
        for (let a = 0; a < this.numBetaSpins; a++) {
          sum += this.coeffBeta.get(i, a) * this.coeffBeta.get(j, a);
        }
        Pb.set(i, j, sum);
      }
    }
    this.densityBeta = Pb;

    // Total density
    const Pt = new Matrix(n, n);
    for (let k = 0; k < n * n; k++) {
      Pt.data[k] = Pa.data[k] + Pb.data[k];
    }
    this.densityTotal = Pt;
  }

  private _lastVxcA?: Matrix;
  private _lastVxcB?: Matrix;

  /** Fock matrix: HF or KS */
  async computeFockMatrix() {
    const n = this.numBasis;
    const xc = this._xcFunctional;
    const hfExFrac = xc ? xc.exactExchangeFraction : 1.0;

    // DFT: compute V_xc
    if (xc && this._grid && this.densityAlpha) {
      const xcResult = computeXC(
        xc, this._grid, this.primitiveShells, this.normFactors,
        n, this.densityAlpha.data as Float64Array,
        this.densityBeta.data as Float64Array,
      );
      this._xcEnergy = xcResult.exc;
      this._lastVxcA = new Matrix(n, n, xcResult.vxcA);
      this._lastVxcB = xcResult.vxcB ? new Matrix(n, n, xcResult.vxcB) : this._lastVxcA;
    }

    // RI-JK path
    if (this._riData && this.coeffAlpha && this.coeffBeta) {
      if (xc && hfExFrac === 0) {
        const J = this._riData.buildCoulomb(this.densityTotal);
        const Fa = new Matrix(n, n);
        const Fb = new Matrix(n, n);
        for (let k = 0; k < n * n; k++) {
          const hj = this.coreHamiltonianMatrix.data[k] + J.data[k];
          Fa.data[k] = hj + (this._lastVxcA ? this._lastVxcA.data[k] : 0);
          Fb.data[k] = hj + (this._lastVxcB ? this._lastVxcB.data[k] : 0);
        }
        this.fockAlpha = Fa;
        this.fockBeta = Fb;
      } else {
        const result = this._riData.buildFockUHF(
          this.densityTotal, this.coreHamiltonianMatrix,
          this.coeffAlpha, this.numAlphaSpins,
          this.coeffBeta, this.numBetaSpins,
        );
        this.fockAlpha = result.fockAlpha;
        this.fockBeta = result.fockBeta;
        if (xc) this.addVxcToFock();
      }
      return;
    }
    // RI first iteration
    if (this._riData && this.densityTotal) {
      const J = this._riData.buildCoulomb(this.densityTotal);
      const Fa = new Matrix(n, n);
      const Fb = new Matrix(n, n);
      for (let k = 0; k < n * n; k++) {
        const hj = this.coreHamiltonianMatrix.data[k] + J.data[k];
        Fa.data[k] = hj + (this._lastVxcA ? this._lastVxcA.data[k] : 0);
        Fb.data[k] = hj + (this._lastVxcB ? this._lastVxcB.data[k] : 0);
      }
      this.fockAlpha = Fa;
      this.fockBeta = Fb;
      return;
    }

    // Pure DFT with stored ERIs: J only
    if (xc && hfExFrac === 0) {
      const Fa = new Matrix(n, n);
      const Fb = new Matrix(n, n);
      for (let i = 0; i < n; i++) {
        for (let j = i; j < n; j++) {
          let J = 0;
          for (let k = 0; k < n; k++) {
            for (let l = 0; l < n; l++) {
              const Ptot = this.densityTotal.get(k, l);
              if (Math.abs(Ptot) < 1e-10) continue;
              J += Ptot * this.eri.get(i, j, k, l);
            }
          }
          const hj = this.coreHamiltonianMatrix.get(i, j) + J;
          Fa.set(i, j, hj); Fa.set(j, i, hj);
          Fb.set(i, j, hj); Fb.set(j, i, hj);
        }
      }
      if (this._lastVxcA) for (let k = 0; k < n * n; k++) Fa.data[k] += this._lastVxcA.data[k];
      if (this._lastVxcB) for (let k = 0; k < n * n; k++) Fb.data[k] += this._lastVxcB.data[k];
      this.fockAlpha = Fa;
      this.fockBeta = Fb;
      return;
    }

    // Hybrid DFT with stored ERIs: J + scaled K
    if (xc && hfExFrac > 0) {
      const Fa = new Matrix(n, n);
      const Fb = new Matrix(n, n);
      for (let i = 0; i < n; i++) {
        for (let j = i; j < n; j++) {
          let J = 0, Ka = 0, Kb = 0;
          for (let k = 0; k < n; k++) {
            for (let l = 0; l < n; l++) {
              const Ptot = this.densityTotal.get(k, l);
              const Pa = this.densityAlpha.get(k, l);
              const Pb = this.densityBeta.get(k, l);
              const Jij = this.eri.get(i, j, k, l);
              const Kij = this.eri.get(i, k, j, l);
              J += Ptot * Jij;
              Ka += Pa * Kij;
              Kb += Pb * Kij;
            }
          }
          const ha = this.coreHamiltonianMatrix.get(i, j) + J - hfExFrac * Ka;
          const hb = this.coreHamiltonianMatrix.get(i, j) + J - hfExFrac * Kb;
          Fa.set(i, j, ha); Fa.set(j, i, ha);
          Fb.set(i, j, hb); Fb.set(j, i, hb);
        }
      }
      if (this._lastVxcA) for (let k = 0; k < n * n; k++) Fa.data[k] += this._lastVxcA.data[k];
      if (this._lastVxcB) for (let k = 0; k < n * n; k++) Fb.data[k] += this._lastVxcB.data[k];
      this.fockAlpha = Fa;
      this.fockBeta = Fb;
      return;
    }

    // Pure HF: try WASM
    if (this.eriBackend !== 'js' && isWasmAvailable()) {
      const result = computeFockUhfWasm(
        this.eri.data, this.densityAlpha.data,
        this.densityBeta.data, this.densityTotal.data,
        this.coreHamiltonianMatrix.data, n,
      );
      if (result) {
        this.fockAlpha = new Matrix(n, n, result.fockAlpha);
        this.fockBeta = new Matrix(n, n, result.fockBeta);
        return;
      }
    }

    // JS fallback (pure HF)
    const Fa = new Matrix(n, n);
    const Fb = new Matrix(n, n);

    for (let i = 0; i < n; i++) {
      for (let j = i; j < n; j++) {
        let Ga = 0, Gb = 0;
        for (let k = 0; k < n; k++) {
          for (let l = 0; l < n; l++) {
            const Ptot = this.densityTotal.get(k, l);
            const Pa = this.densityAlpha.get(k, l);
            const Pb = this.densityBeta.get(k, l);
            if (Math.abs(Ptot) < 1e-10 && Math.abs(Pa) < 1e-10) continue;
            const Jij = this.eri.get(i, j, k, l);
            const Kij = this.eri.get(i, k, j, l);
            Ga += Ptot * Jij - Pa * Kij;
            Gb += Ptot * Jij - Pb * Kij;
          }
        }
        const Ha = this.coreHamiltonianMatrix.get(i, j) + Ga;
        const Hb = this.coreHamiltonianMatrix.get(i, j) + Gb;
        Fa.set(i, j, Ha); Fa.set(j, i, Ha);
        Fb.set(i, j, Hb); Fb.set(j, i, Hb);
      }
    }

    this.fockAlpha = Fa;
    this.fockBeta = Fb;
  }

  private addVxcToFock() {
    const n2 = this.numBasis * this.numBasis;
    if (this._lastVxcA) for (let k = 0; k < n2; k++) this.fockAlpha.data[k] += this._lastVxcA.data[k];
    if (this._lastVxcB) for (let k = 0; k < n2; k++) this.fockBeta.data[k] += this._lastVxcB.data[k];
  }

  /** E = 0.5 * (Tr(Pα*(H+Fα)) + Tr(Pβ*(H+Fβ))), with DFT correction */
  computeEnergy(): number {
    const HpFa = this.coreHamiltonianMatrix.addMatrix(this.fockAlpha);
    const HpFb = this.coreHamiltonianMatrix.addMatrix(this.fockBeta);
    this.energy = 0.5 * (traceProduct(this.densityAlpha, HpFa) + traceProduct(this.densityBeta, HpFb));

    if (this._xcFunctional && this._lastVxcA && this._lastVxcB) {
      const trPVxc = traceProduct(this.densityAlpha, this._lastVxcA) +
                     traceProduct(this.densityBeta, this._lastVxcB);
      this.energy += this._xcEnergy - 0.5 * trPVxc;
    }

    return this.energy;
  }

  /** Update Fock matrices using the selected SCF accelerator */
  updateFockMatrix() {
    const n = this.numBasis;

    // Commutator errors for alpha and beta
    const FPSa = matmul(matmul(this.fockAlpha, this.densityAlpha), this.overlapMatrix);
    const SPFa = matmul(matmul(this.overlapMatrix, this.densityAlpha), this.fockAlpha);
    const errorA = FPSa.subMatrix(SPFa);

    const FPSb = matmul(matmul(this.fockBeta, this.densityBeta), this.overlapMatrix);
    const SPFb = matmul(matmul(this.overlapMatrix, this.densityBeta), this.fockBeta);
    const errorB = FPSb.subMatrix(SPFb);

    if (this.scfAccelerator) {
      const result = this.scfAccelerator.accelerate({
        fock: [copyToFloatArray(this.fockAlpha.data), copyToFloatArray(this.fockBeta.data)],
        density: [copyToFloatArray(this.densityAlpha.data), copyToFloatArray(this.densityBeta.data)],
        error: [copyToFloatArray(errorA.data), copyToFloatArray(errorB.data)],
        energy: this.energy,
        overlap: this.overlapMatrix.data,
        numBasis: n,
        iteration: 0,
      });
      if (result.fock) {
        this.fockAlpha = new Matrix(n, n, result.fock[0]);
        this.fockBeta = new Matrix(n, n, result.fock[1]);
      }
    }
  }

  getTotalEnergy(): number {
    return this.energy + this.nuclearRepulsionEnergy;
  }

  // Public getters
  get orbitalEnergiesAlpha(): FloatArray { return this._orbitalEnergiesAlpha; }
  get orbitalEnergiesBeta(): FloatArray { return this._orbitalEnergiesBeta; }
  get coefficientsAlpha(): Matrix { return this.coeffAlpha; }
  get coefficientsBeta(): Matrix { return this.coeffBeta; }
  get densityAlphaMatrix(): Matrix { return this.densityAlpha; }
  get densityBetaMatrix(): Matrix { return this.densityBeta; }
  get density(): Matrix { return this.densityTotal; }
  get eriStore(): ERIStored { return this.eri; }
  get fockAlphaMatrix(): Matrix { return this.fockAlpha; }
  get fockBetaMatrix(): Matrix { return this.fockBeta; }
}
