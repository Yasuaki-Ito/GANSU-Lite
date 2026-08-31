/** ROHF (Restricted Open-shell Hartree-Fock) — Guest-Saunders effective Fock matrix */

import { HF } from './hf';
import { Matrix, type FloatArray, copyToFloatArray } from '../linalg/matrix';
import { jacobiEigen } from '../linalg/eigendecomposition';
import { matmul, matmulAtB, matmulABt, traceProduct } from '../linalg/matmul';
import { computeERIs } from './integrals2e';
import type { EriBackend } from './hf';
import { ERIStored } from './eri';
import { isWasmAvailable, computeFockUhfWasm } from './eriWasm';
import type { RIData } from './ri';
import { computeXC } from './xcIntegration';

export class ROHF extends HF {
  private fockAlpha!: Matrix;
  private fockBeta!: Matrix;
  private densAlpha!: Matrix;
  private densBeta!: Matrix;
  private densTotal!: Matrix;
  private coeff!: Matrix;
  private coeffOrth!: Matrix; // MO coefficients in orthogonal basis
  private _orbitalEnergies!: FloatArray;
  private eri!: ERIStored;
  private _riData?: RIData;
  private energy = 0;
  private scfIter = 0;
  private _lastVxcA?: Matrix;
  private _lastVxcB?: Matrix;

  /** Number of doubly occupied (closed-shell) orbitals */
  private get numClosed(): number { return this.numBetaSpins; }
  /** Number of singly occupied (open-shell) orbitals */
  private get numOpen(): number { return this.numAlphaSpins - this.numBetaSpins; }

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

  protected setInitialDensity(P: Matrix) {
    const n = this.numBasis;
    const nel = this.numAlphaSpins + this.numBetaSpins;
    const fracA = nel > 0 ? this.numAlphaSpins / nel : 0.5;
    const fracB = nel > 0 ? this.numBetaSpins / nel : 0.5;
    const Pa = new Matrix(n, n);
    const Pb = new Matrix(n, n);
    for (let k = 0; k < n * n; k++) {
      Pa.data[k] = P.data[k] * fracA;
      Pb.data[k] = P.data[k] * fracB;
    }
    this.densAlpha = Pa;
    this.densBeta = Pb;
    this.densTotal = P.clone();
  }

  guessInitialFockMatrix() {
    this.fockAlpha = this.coreHamiltonianMatrix.clone();
    this.fockBeta = this.coreHamiltonianMatrix.clone();
  }

  guessInitialFockMatrixGWH() {
    const gwh = this.buildGWHFock();
    this.fockAlpha = gwh;
    this.fockBeta = gwh.clone();
  }

  /** Build effective Fock matrix and diagonalize.
   *  First iteration: diagonalize average Fock.
   *  Subsequent iterations: Guest-Saunders coupling in MO basis. */
  computeCoefficientMatrix() {
    const n = this.numBasis;
    const X = this.transformMatrix;

    // Transform Fα, Fβ to orthogonal basis
    const FaOrth = matmul(matmulAtB(X, this.fockAlpha), X);
    const FbOrth = matmul(matmulAtB(X, this.fockBeta), X);

    let FeffOrth: Matrix;

    if (!this.coeffOrth) {
      // First iteration: use average Fock matrix
      FeffOrth = new Matrix(n, n);
      for (let k = 0; k < n * n; k++) {
        FeffOrth.data[k] = 0.5 * (FaOrth.data[k] + FbOrth.data[k]);
      }
    } else {
      // Transform to MO basis using previous orthogonal coefficients
      const Cp = this.coeffOrth;
      const FaMO = matmul(matmulAtB(Cp, FaOrth), Cp); // C'^T Fα' C'
      const FbMO = matmul(matmulAtB(Cp, FbOrth), Cp); // C'^T Fβ' C'

      // Build effective Fock in MO basis (Guest-Saunders coupling)
      const nc = this.numClosed;
      const no = this.numOpen;
      const FeffMO = new Matrix(n, n);

      for (let i = 0; i < n; i++) {
        for (let j = i; j < n; j++) {
          // Determine subspace: 0=closed, 1=open, 2=virtual
          const si = i < nc ? 0 : (i < nc + no ? 1 : 2);
          const sj = j < nc ? 0 : (j < nc + no ? 1 : 2);
          const smin = Math.min(si, sj);
          const smax = Math.max(si, sj);

          const fa = FaMO.get(i, j);
          const fb = FbMO.get(i, j);

          let val: number;
          if (smin === 0 && smax === 0)      val = 0.5 * (fa + fb); // closed-closed
          else if (smin === 0 && smax === 1) val = fb;               // closed-open
          else if (smin === 0 && smax === 2) val = 0.5 * (fa + fb); // closed-virtual
          else if (smin === 1 && smax === 1) val = fa;               // open-open
          else if (smin === 1 && smax === 2) val = fa;               // open-virtual
          else                               val = 0.5 * (fa + fb); // virtual-virtual

          FeffMO.set(i, j, val);
          FeffMO.set(j, i, val);
        }
      }

      // Transform back to orthogonal basis: Feff' = C' FeffMO C'^T
      FeffOrth = matmulABt(matmul(Cp, FeffMO), Cp);

      // Accelerator on effective Fock using commutator [Feff', P'_total] as error
      // P'_total in orthogonal basis from previous MO coefficients
      const PtOrth = new Matrix(n, n);
      for (let i = 0; i < n; i++) {
        for (let j = i; j < n; j++) {
          let sum = 0;
          for (let a = 0; a < this.numAlphaSpins; a++) {
            sum += Cp.get(i, a) * Cp.get(j, a);
          }
          for (let a = 0; a < this.numBetaSpins; a++) {
            sum += Cp.get(i, a) * Cp.get(j, a);
          }
          PtOrth.set(i, j, sum);
          PtOrth.set(j, i, sum);
        }
      }
      // error = [Feff', P'] = Feff'*P' - P'*Feff' (zero at convergence)
      const FP = matmul(FeffOrth, PtOrth);
      const PF = matmul(PtOrth, FeffOrth);
      const error = FP.subMatrix(PF);

      if (this.scfAccelerator) {
        const result = this.scfAccelerator.accelerate({
          fock: [copyToFloatArray(FeffOrth.data)],
          density: [copyToFloatArray(PtOrth.data)],
          error: [copyToFloatArray(error.data)],
          energy: this.energy,
          overlap: this.overlapMatrix.data,
          numBasis: n,
          iteration: this.scfIter++,
        });
        if (result.fock) {
          FeffOrth = new Matrix(n, n, result.fock[0]);
        }
      }
    }

    // Diagonalize effective Fock in orthogonal basis
    const eig = jacobiEigen(FeffOrth);
    this.coeffOrth = eig.eigenvectors;
    this.coeff = matmul(X, eig.eigenvectors);
    this._orbitalEnergies = eig.eigenvalues;
  }

  /** Build Pα, Pβ, Ptotal from single coefficient matrix C.
   *  Pα uses numAlphaSpins occupied orbitals, Pβ uses numBetaSpins. */
  computeDensityMatrix() {
    const n = this.numBasis;
    const C = this.coeff;

    // Alpha density
    const Pa = new Matrix(n, n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        let sum = 0;
        for (let a = 0; a < this.numAlphaSpins; a++) {
          sum += C.get(i, a) * C.get(j, a);
        }
        Pa.set(i, j, sum);
      }
    }
    this.densAlpha = Pa;

    // Beta density
    const Pb = new Matrix(n, n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        let sum = 0;
        for (let a = 0; a < this.numBetaSpins; a++) {
          sum += C.get(i, a) * C.get(j, a);
        }
        Pb.set(i, j, sum);
      }
    }
    this.densBeta = Pb;

    // Total density
    const Pt = new Matrix(n, n);
    for (let k = 0; k < n * n; k++) {
      Pt.data[k] = Pa.data[k] + Pb.data[k];
    }
    this.densTotal = Pt;
  }

  /** Fock matrix construction — same structure as UHF, with DFT/KS support */
  async computeFockMatrix() {
    const n = this.numBasis;
    const xc = this._xcFunctional;
    const hfExFrac = xc ? xc.exactExchangeFraction : 1.0;

    // DFT: compute V_xc from α/β densities
    if (xc && this._grid && this.densAlpha) {
      const xcResult = computeXC(
        xc, this._grid, this.primitiveShells, this.normFactors,
        n, this.densAlpha.data as Float64Array,
        this.densBeta.data as Float64Array,
      );
      this._xcEnergy = xcResult.exc;
      this._lastVxcA = new Matrix(n, n, xcResult.vxcA);
      this._lastVxcB = xcResult.vxcB ? new Matrix(n, n, xcResult.vxcB) : this._lastVxcA;
    }

    // RI-JK path (requires coefficients for exchange)
    if (this._riData && this.coeff) {
      if (xc && hfExFrac === 0) {
        // Pure DFT: J only
        const J = this._riData.buildCoulomb(this.densTotal);
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
          this.densTotal, this.coreHamiltonianMatrix,
          this.coeff, this.numAlphaSpins,
          this.coeff, this.numBetaSpins,
        );
        this.fockAlpha = result.fockAlpha;
        this.fockBeta = result.fockBeta;
        if (xc) this.addVxcToFock();
      }
      return;
    }
    // RI first iteration with SAD density but no coefficients: J-only approximation
    if (this._riData && this.densTotal) {
      const J = this._riData.buildCoulomb(this.densTotal);
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
              const Ptot = this.densTotal.get(k, l);
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
              const Ptot = this.densTotal.get(k, l);
              const Pa = this.densAlpha.get(k, l);
              const Pb = this.densBeta.get(k, l);
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
        this.eri.data, this.densAlpha.data,
        this.densBeta.data, this.densTotal.data,
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
            const Ptot = this.densTotal.get(k, l);
            const Pa = this.densAlpha.get(k, l);
            const Pb = this.densBeta.get(k, l);
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
    this.energy = 0.5 * (traceProduct(this.densAlpha, HpFa) + traceProduct(this.densBeta, HpFb));

    if (this._xcFunctional && this._lastVxcA && this._lastVxcB) {
      const trPVxc = traceProduct(this.densAlpha, this._lastVxcA) +
                     traceProduct(this.densBeta, this._lastVxcB);
      this.energy += this._xcEnergy - 0.5 * trPVxc;
    }

    return this.energy;
  }

  /** Accelerator is applied to effective Fock in computeCoefficientMatrix() */
  updateFockMatrix() {
    // no-op: ROHF applies accelerator to the effective Fock, not Fα/Fβ individually
  }

  getTotalEnergy(): number {
    return this.energy + this.nuclearRepulsionEnergy;
  }

  // Public getters
  get orbitalEnergies(): FloatArray { return this._orbitalEnergies; }
  get coefficients(): Matrix { return this.coeff; }
  get density(): Matrix { return this.densTotal; }
  get densityAlphaMatrix(): Matrix { return this.densAlpha; }
  get densityBetaMatrix(): Matrix { return this.densBeta; }
  get eriStore(): ERIStored { return this.eri; }
  get fockAlphaMatrix(): Matrix { return this.fockAlpha; }
  get fockBetaMatrix(): Matrix { return this.fockBeta; }
}
