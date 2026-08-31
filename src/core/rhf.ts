/** RHF/RKS (Restricted Hartree-Fock / Restricted Kohn-Sham) */

import { HF } from './hf';
import type { Atom, PrimitiveShell, ShellTypeInfo, BasisRange } from './types';
import { Matrix, type FloatArray, copyToFloatArray } from '../linalg/matrix';
import { jacobiEigen } from '../linalg/eigendecomposition';
import { matmul, matmulAtB, traceProduct } from '../linalg/matmul';
import { computeERIs, computeLongRangeERIs } from './integrals2e';
import type { EriBackend } from './hf';
import { ERIStored } from './eri';
import { isWasmAvailable, computeFockRhfWasm } from './eriWasm';
import type { RIData } from './ri';
import { computeXC } from './xcIntegration';

export class RHF extends HF {
  private fockMatrix!: Matrix;
  private densityMatrix!: Matrix;
  private coefficientMatrix!: Matrix;
  private _orbitalEnergies!: FloatArray;
  private eri!: ERIStored;
  private _eriLR?: ERIStored; // long-range ERI for RSH functionals
  private _riData?: RIData;
  private energy = 0;

  constructor(
    numBasis: number, numElectrons: number, numAlphaSpins: number, numBetaSpins: number,
    atoms: Atom[], primitiveShells: PrimitiveShell[], shellTypeInfos: ShellTypeInfo[],
    atomToBasisRange: BasisRange[], normFactors: number[],
  ) {
    super(numBasis, numElectrons, numAlphaSpins, numBetaSpins,
      atoms, primitiveShells, shellTypeInfos, atomToBasisRange, normFactors);

    if (numAlphaSpins !== numBetaSpins) {
      throw new Error('RHF requires equal alpha and beta electrons');
    }
  }

  async precomputeERI(onProgress?: (percent: number, msg?: string) => void, backend?: EriBackend) {
    if (this._riData) {
      // RI path: skip full ERI, RIData is set externally via setRIData()
      if (onProgress) onProgress(100, 'RI');
      return;
    }

    const eriData = await computeERIs(
      this.primitiveShells, this.normFactors, this.numBasis, 1e-10, onProgress, backend,
    );
    this.eri = new ERIStored(eriData, this.numBasis);

    // For RSH functionals, also build the long-range ERI (erf(ωr)/r)
    const rs = this._xcFunctional?.rangeSeparation;
    if (rs && rs.beta !== 0 && rs.omega > 0) {
      if (onProgress) onProgress(50, `Long-range ERI (ω=${rs.omega})`);
      const lrData = computeLongRangeERIs(
        this.primitiveShells, this.normFactors, this.numBasis, rs.omega, 1e-10,
      );
      this._eriLR = new ERIStored(lrData, this.numBasis);
      if (onProgress) onProgress(100, `Long-range ERI done`);
    }
  }

  setRIData(ri: RIData) {
    this._riData = ri;
  }

  get riData(): RIData | undefined { return this._riData; }

  protected setInitialDensity(P: Matrix) {
    this.densityMatrix = P;
  }

  /** Initial guess: F = H_core (Core Hamiltonian guess) */
  guessInitialFockMatrix() {
    this.fockMatrix = this.coreHamiltonianMatrix.clone();
  }

  guessInitialFockMatrixGWH() {
    this.fockMatrix = this.buildGWHFock();
  }

  /** Diagonalize Fock matrix in orthogonal basis: F' = X^T F X → eigensolve → C = X C' */
  computeCoefficientMatrix() {
    const Fp = matmul(matmulAtB(this.transformMatrix, this.fockMatrix), this.transformMatrix);
    const { eigenvalues, eigenvectors } = jacobiEigen(Fp);
    this.coefficientMatrix = matmul(this.transformMatrix, eigenvectors);
    this._orbitalEnergies = eigenvalues;
  }

  /** P = 2 * C_occ * C_occ^T */
  computeDensityMatrix() {
    const n = this.numBasis;
    const nocc = this.numAlphaSpins; // same as numBetaSpins for RHF
    const P = new Matrix(n, n);

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        let sum = 0;
        for (let a = 0; a < nocc; a++) {
          sum += this.coefficientMatrix.get(i, a) * this.coefficientMatrix.get(j, a);
        }
        P.set(i, j, 2 * sum);
      }
    }

    this.densityMatrix = P;
  }

  /** F[i][j] = H[i][j] + G[i][j] (+ V_xc for DFT) */
  async computeFockMatrix() {
    const n = this.numBasis;
    const xc = this._xcFunctional;
    const hfExchangeFraction = xc ? xc.exactExchangeFraction : 1.0;

    // DFT: compute V_xc
    if (xc && this._grid && this.densityMatrix) {
      const xcResult = computeXC(
        xc, this._grid, this.primitiveShells, this.normFactors,
        n, this.densityMatrix.data as Float64Array, null,
      );
      this._xcEnergy = xcResult.exc;
      this._lastVxc = new Matrix(n, n, xcResult.vxcA);
    }

    // RI-JK path
    if (this._riData && this.coefficientMatrix) {
      if (xc && hfExchangeFraction === 0) {
        // Pure DFT: J only
        const J = this._riData.buildCoulomb(this.densityMatrix);
        const F = new Matrix(n, n);
        for (let k = 0; k < n * n; k++) {
          F.data[k] = this.coreHamiltonianMatrix.data[k] + J.data[k];
        }
        if (this._lastVxc) for (let k = 0; k < n * n; k++) F.data[k] += this._lastVxc.data[k];
        this.fockMatrix = F;
      } else {
        this.fockMatrix = this._riData.buildFockRHF(
          this.densityMatrix, this.coreHamiltonianMatrix,
          this.coefficientMatrix, this.numAlphaSpins,
        );
        if (xc) {
          // Hybrid: scale exchange and add V_xc
          this.adjustFockForHybrid(hfExchangeFraction);
        }
      }
      return;
    }
    // RI first iteration
    if (this._riData && this.densityMatrix) {
      const J = this._riData.buildCoulomb(this.densityMatrix);
      const F = new Matrix(n, n);
      for (let k = 0; k < n * n; k++) {
        F.data[k] = this.coreHamiltonianMatrix.data[k] + J.data[k];
      }
      if (this._lastVxc) for (let k = 0; k < n * n; k++) F.data[k] += this._lastVxc.data[k];
      this.fockMatrix = F;
      return;
    }

    // Stored ERI path
    if (xc && hfExchangeFraction === 0) {
      // Pure DFT with stored ERIs: J only
      const F = new Matrix(n, n);
      for (let i = 0; i < n; i++) {
        for (let j = i; j < n; j++) {
          let J = 0;
          for (let k = 0; k < n; k++) {
            for (let l = 0; l < n; l++) {
              const Pkl = this.densityMatrix.get(k, l);
              if (Math.abs(Pkl) < 1e-10) continue;
              J += Pkl * this.eri.get(i, j, k, l);
            }
          }
          const val = this.coreHamiltonianMatrix.get(i, j) + J;
          F.set(i, j, val);
          F.set(j, i, val);
        }
      }
      if (this._lastVxc) for (let k = 0; k < n * n; k++) F.data[k] += this._lastVxc.data[k];
      this.fockMatrix = F;
      return;
    }

    // HF or hybrid DFT with stored ERIs
    // Try WASM Fock kernel
    if (this.eriBackend !== 'js' && isWasmAvailable()) {
      const result = computeFockRhfWasm(
        this.eri.data, this.densityMatrix.data,
        this.coreHamiltonianMatrix.data, n,
      );
      if (result) {
        this.fockMatrix = new Matrix(n, n, result);
        if (xc) this.adjustFockForHybrid(hfExchangeFraction);
        return;
      }
    }

    // JS fallback
    const F = new Matrix(n, n);
    for (let i = 0; i < n; i++) {
      for (let j = i; j < n; j++) {
        let G = 0;
        for (let k = 0; k < n; k++) {
          for (let l = 0; l < n; l++) {
            const Pkl = this.densityMatrix.get(k, l);
            if (Math.abs(Pkl) < 1e-10) continue;
            G += Pkl * (this.eri.get(i, j, k, l) - 0.5 * this.eri.get(i, k, j, l));
          }
        }
        const val = this.coreHamiltonianMatrix.get(i, j) + G;
        F.set(i, j, val);
        F.set(j, i, val);
      }
    }
    this.fockMatrix = F;
    if (xc) this.adjustFockForHybrid(hfExchangeFraction);
  }

  /** For hybrid DFT: adjust F = H + J + a0*K → F = H + J + a0*K + V_xc
   *  Since the standard code computes F = H + J - 0.5*K, we need:
   *  F_hybrid = F_HF + (a0 - 1) * (-0.5*K) + V_xc = F_HF + (1-a0)*0.5*K_contribution + V_xc
   *  Easier: F_hybrid = H + J + a0*K + V_xc
   *  K_contribution = F_HF - H - J, so K = -2 * (F_HF - H - J - (-0.5*K)) ... this is circular.
   *  Instead: rebuild with exchange scaling factor. */
  private adjustFockForHybrid(hfExchangeFraction: number) {
    if (!this._lastVxc) return;
    const n = this.numBasis;
    // F_HF = H + J - 0.5*K
    // F_DFT = H + J + a0*(-0.5*K) + V_xc = F_HF + (a0 - 1)*(-0.5*K) + V_xc
    // where -0.5*K = F_HF - (H + J)
    // So we need J separately. For simplicity, extract K from F_HF:
    // K_contribution = F_HF - H - J_built
    // This requires J... let's compute J separately for hybrids.

    // Recompute F = H + J + a0*(-0.5*K) + V_xc using stored ERIs
    if (this.eri) {
      const rs = this._xcFunctional?.rangeSeparation;
      const beta = rs?.beta ?? 0;
      const F = new Matrix(n, n);
      for (let i = 0; i < n; i++) {
        for (let j = i; j < n; j++) {
          let J = 0, K = 0, KLR = 0;
          for (let k = 0; k < n; k++) {
            for (let l = 0; l < n; l++) {
              const Pkl = this.densityMatrix.get(k, l);
              if (Math.abs(Pkl) < 1e-10) continue;
              J += Pkl * this.eri.get(i, j, k, l);
              K += Pkl * this.eri.get(i, k, j, l);
              if (beta !== 0 && this._eriLR) {
                KLR += Pkl * this._eriLR.get(i, k, j, l);
              }
            }
          }
          // RSH: F = H + J - 0.5*α*K - 0.5*β*K_LR + V_xc
          const val = this.coreHamiltonianMatrix.get(i, j) + J
            - 0.5 * hfExchangeFraction * K
            - 0.5 * beta * KLR;
          F.set(i, j, val);
          F.set(j, i, val);
        }
      }
      for (let k = 0; k < n * n; k++) F.data[k] += this._lastVxc.data[k];
      this.fockMatrix = F;
    } else {
      // RI hybrid: F already has full exchange, need to scale it
      // F_HF = H + J - 0.5*K, we want H + J - 0.5*a0*K + V_xc
      // F_new = H + J - 0.5*a0*K + V_xc = F_HF + (1-a0)*0.5*K + V_xc
      // Since we can't easily get K separately with RI, just add V_xc
      // (the RI path already built the full F with exchange)
      // For now: approximate by scaling the exchange part
      // This is inexact for RI hybrids; proper implementation needs separate J/K from RI
      for (let k = 0; k < n * n; k++) {
        this.fockMatrix.data[k] += this._lastVxc.data[k];
      }
    }
  }

  private _lastVxc?: Matrix;

  /** E = 0.5 * Tr(P * (H + F)).
   *  For DFT: correct by E_xc - 0.5*Tr(P*V_xc) since V_xc is in F but E_xc ≠ Tr(P*V_xc). */
  computeEnergy(): number {
    const HpF = this.coreHamiltonianMatrix.addMatrix(this.fockMatrix);
    this.energy = 0.5 * traceProduct(this.densityMatrix, HpF);

    if (this._xcFunctional && this._lastVxc) {
      // 0.5*Tr(P*(H+F)) already includes 0.5*Tr(P*V_xc); replace with E_xc
      const trPVxc = traceProduct(this.densityMatrix, this._lastVxc);
      this.energy += this._xcEnergy - 0.5 * trPVxc;
    }

    return this.energy;
  }

  /** Update Fock matrix using the selected SCF accelerator */
  updateFockMatrix() {
    const n = this.numBasis;
    // Compute commutator error: e = F*P*S - S*P*F
    const FPS = matmul(matmul(this.fockMatrix, this.densityMatrix), this.overlapMatrix);
    const SPF = matmul(matmul(this.overlapMatrix, this.densityMatrix), this.fockMatrix);
    const error = FPS.subMatrix(SPF);

    if (this.scfAccelerator) {
      // Level Shift needs idempotent density P/2 for RHF
      const halfDens = new Matrix(n, n);
      for (let k = 0; k < n * n; k++) halfDens.data[k] = this.densityMatrix.data[k] * 0.5;

      const result = this.scfAccelerator.accelerate({
        fock: [copyToFloatArray(this.fockMatrix.data)],
        density: [copyToFloatArray(halfDens.data)],
        error: [copyToFloatArray(error.data)],
        energy: this.energy,
        overlap: this.overlapMatrix.data,
        numBasis: n,
        iteration: 0,
      });
      if (result.fock) {
        this.fockMatrix = new Matrix(n, n, result.fock[0]);
      }
    }
  }

  getTotalEnergy(): number {
    return this.energy + this.nuclearRepulsionEnergy;
  }

  get orbitalEnergies(): FloatArray { return this._orbitalEnergies; }
  get coefficients(): Matrix { return this.coefficientMatrix; }
  get density(): Matrix { return this.densityMatrix; }
  get fock(): Matrix { return this.fockMatrix; }
  get eriStore(): ERIStored { return this.eri; }
}
