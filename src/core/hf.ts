/** Abstract HF/KS base class with SCF loop — port of GANSU hf.cu */

import type { Atom, PrimitiveShell, ShellTypeInfo, BasisRange } from './types';
import { Matrix } from '../linalg/matrix';
import { jacobiEigen } from '../linalg/eigendecomposition';
import { matmul, matmulAtB } from '../linalg/matmul';
import { computeOneElectronIntegrals } from './integrals1e';
import type { SCFAccelMethod, SCFAccelParams, SCFAccelerator } from './scfAccelerator';
import { createAccelerator } from './accelFactory';
import type { XCFunctional } from './xcFunctional';
import type { GridPoint, GridLevel } from './grid';

export type EriBackend = 'auto' | 'js' | 'wasm';
export type EriAlgorithm = 'auto' | 'md' | 'os' | 'rys' | 'hgp';

export interface SCFCallbacks {
  onIteration?: (iter: number, energy: number, deltaE: number) => void;
  onProgress?: (message: string) => void;
  onPhase?: (phase: string, status: 'start' | 'done') => void;
  eriBackend?: EriBackend;
  scfAccelMethod?: SCFAccelMethod;
  scfAccelParams?: SCFAccelParams;
}

export abstract class HF {
  readonly numBasis: number;
  readonly numElectrons: number;
  readonly numAlphaSpins: number;
  readonly numBetaSpins: number;
  readonly atoms: Atom[];
  readonly primitiveShells: PrimitiveShell[];
  readonly shellTypeInfos: ShellTypeInfo[];
  readonly atomToBasisRange: BasisRange[];
  readonly normFactors: number[];

  protected overlapMatrix!: Matrix;
  protected coreHamiltonianMatrix!: Matrix;
  protected kineticMatrix!: Matrix;
  protected transformMatrix!: Matrix;
  protected nuclearRepulsionEnergy = 0;

  get overlap(): Matrix { return this.overlapMatrix; }
  get coreHamiltonian(): Matrix { return this.coreHamiltonianMatrix; }
  get kinetic(): Matrix { return this.kineticMatrix; }
  get nuclearRepulsion(): number { return this.nuclearRepulsionEnergy; }

  protected maxIter = 200;
  protected eriBackend: EriBackend = 'auto';
  protected scfAccelerator?: SCFAccelerator;

  protected readonly convergenceThreshold = 1e-8;

  // DFT configuration (null = pure HF)
  protected _xcFunctional?: XCFunctional;
  protected _grid?: GridPoint[];
  protected _coarseGrid?: GridPoint[];
  protected _targetGrid?: GridPoint[];
  protected _gridLevel: GridLevel = 'medium';
  protected _xcEnergy = 0;
  protected _gridElectrons?: number;
  protected _gridSwitched = false;

  get xcFunctional(): XCFunctional | undefined { return this._xcFunctional; }
  get grid(): GridPoint[] | undefined { return this._grid; }
  get xcEnergy(): number { return this._xcEnergy; }
  /** ∫ρ dr evaluated on the XC grid (DFT only). Deviation from the electron
   *  count is the standard sanity check on grid quality. */
  get gridElectrons(): number | undefined { return this._gridElectrons; }

  setDFT(functional: XCFunctional, grid: GridPoint[], coarseGrid?: GridPoint[]) {
    this._xcFunctional = functional;
    this._targetGrid = grid;
    this._coarseGrid = coarseGrid;
    // Start with coarse grid if available, otherwise use target
    this._grid = coarseGrid ?? grid;
    this._gridSwitched = !coarseGrid;
  }

  /** Switch to target (fine) grid when SCF is nearly converged */
  protected switchToTargetGrid(): void {
    if (!this._gridSwitched && this._targetGrid) {
      this._grid = this._targetGrid;
      this._gridSwitched = true;
      // The DIIS/ADIIS history holds Fock matrices built on the coarse grid. They are
      // inconsistent with the new objective, and when the density can no longer change
      // (H2/STO-3G: fixed by symmetry) the error vectors vanish and the extrapolation
      // returns an arbitrary mixture of coarse- and target-grid Fock matrices — the
      // energy still converges but orbital energies come out 0.5 eV off. Start afresh.
      this.scfAccelerator?.reset();
    }
  }

  get isDFT(): boolean { return this._xcFunctional !== undefined; }

  constructor(
    numBasis: number, numElectrons: number, numAlphaSpins: number, numBetaSpins: number,
    atoms: Atom[], primitiveShells: PrimitiveShell[], shellTypeInfos: ShellTypeInfo[],
    atomToBasisRange: BasisRange[], normFactors: number[],
  ) {
    this.numBasis = numBasis;
    this.numElectrons = numElectrons;
    this.numAlphaSpins = numAlphaSpins;
    this.numBetaSpins = numBetaSpins;
    this.atoms = atoms;
    this.primitiveShells = primitiveShells;
    this.shellTypeInfos = shellTypeInfos;
    this.atomToBasisRange = atomToBasisRange;
    this.normFactors = normFactors;
  }

  /** Compute nuclear repulsion energy E_nuc = Σ_i Σ_{j>i} Z_i*Z_j / r_ij */
  computeNuclearRepulsionEnergy(): number {
    let E = 0;
    for (let i = 0; i < this.atoms.length; i++) {
      for (let j = i + 1; j < this.atoms.length; j++) {
        const ai = this.atoms[i], aj = this.atoms[j];
        const dx = ai.coordinate.x - aj.coordinate.x;
        const dy = ai.coordinate.y - aj.coordinate.y;
        const dz = ai.coordinate.z - aj.coordinate.z;
        const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
        E += ai.atomicNumber * aj.atomicNumber / r;
      }
    }
    this.nuclearRepulsionEnergy = E;
    return E;
  }

  /** Compute overlap and core Hamiltonian matrices */
  computeCoreHamiltonianMatrix() {
    const { overlap, coreHamiltonian, kinetic } = computeOneElectronIntegrals(
      this.primitiveShells, this.atoms, this.normFactors, this.numBasis,
    );
    this.overlapMatrix = overlap;
    this.coreHamiltonianMatrix = coreHamiltonian;
    this.kineticMatrix = kinetic;
  }

  /**
   * Cartesian → spherical map T (numBasis × nSph), set when the basis is used with
   * pure d/f functions. Every integral stays Cartesian; T only restricts the
   * variational space, inside computeTransformMatrix.
   */
  protected sphericalTransform: Matrix | null = null;

  setSphericalTransform(T: Matrix) {
    if (T.rows !== this.numBasis) throw new Error(`Spherical map has ${T.rows} rows, basis has ${this.numBasis}`);
    this.sphericalTransform = T;
  }

  /**
   * Columns of X that span nothing: linearly dependent directions, and in a
   * spherical basis the Cartesian-only combinations (the s hidden in a d shell, the
   * p in an f shell). Diagonalising X^T F X leaves one fake orbital per such column.
   */
  protected nullDirections: number[] = [];

  /**
   * Energy given to the fake orbitals so they sort above every real one. They have
   * identically zero coefficients, so any correlation method sees them as virtuals
   * whose integrals all vanish and that therefore contribute nothing.
   */
  static readonly NULL_ORBITAL_ENERGY = 1e4;

  /** Real molecular orbitals: the first numMO columns of C (and entries of ε). */
  get numMO(): number { return this.numBasis - this.nullDirections.length; }

  /**
   * Compute the orthogonalising transformation X (numBasis × numBasis).
   *
   * Cartesian basis: X = U s^{-1/2}, with near-singular directions zeroed.
   * Spherical basis: X = T U' s'^{-1/2}, where U', s' diagonalise T^T S T. Only the
   * first nSph columns are non-zero, so every matrix keeps its Cartesian size —
   * nothing downstream has to know about two dimensions — while C = X C' can only
   * ever contain pure functions.
   */
  computeTransformMatrix() {
    const n = this.numBasis;
    const T = this.sphericalTransform;
    const S = T ? matmul(matmulAtB(T, this.overlapMatrix), T) : this.overlapMatrix;
    const m = S.rows;
    const { eigenvalues, eigenvectors } = jacobiEigen(S);

    // U s^{-1/2} (filter out near-zero eigenvalues for linear dependence)
    const sInvSqrt = new Matrix(m, m);
    for (let i = 0; i < m; i++) {
      if (eigenvalues[i] > 1e-6) {
        sInvSqrt.set(i, i, 1.0 / Math.sqrt(eigenvalues[i]));
      }
    }
    const Xsub = matmul(eigenvectors, sInvSqrt);

    if (!T) {
      this.transformMatrix = Xsub;
    } else {
      const Xsph = matmul(T, Xsub); // n × m
      const X = new Matrix(n, n);
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < m; j++) X.data[i * n + j] = Xsph.data[i * m + j];
      }
      this.transformMatrix = X;
    }

    this.nullDirections = [];
    const X = this.transformMatrix;
    for (let j = 0; j < n; j++) {
      let zero = true;
      for (let i = 0; i < n && zero; i++) if (X.data[i * n + j] !== 0) zero = false;
      if (zero) this.nullDirections.push(j);
    }
  }

  /**
   * F' = X^T F X, with the null directions pushed to NULL_ORBITAL_ENERGY. Without
   * the shift they diagonalise to ε = 0, which sits below a positive LUMO and would
   * be reported as one.
   */
  protected toOrthogonalBasis(F: Matrix): Matrix {
    const X = this.transformMatrix;
    const Fp = matmul(matmulAtB(X, F), X);
    for (const k of this.nullDirections) Fp.set(k, k, HF.NULL_ORBITAL_ENERGY);
    return Fp;
  }

  /**
   * DIIS error for a restricted variational space. FPS - SPF in the AO basis does
   * not vanish at convergence when F couples the occupied space to functions
   * outside it, so the commutator is measured in the orthogonal basis instead.
   * A full Cartesian basis is returned unchanged, keeping that path bit-identical.
   */
  protected projectError(e: Matrix): Matrix {
    if (this.nullDirections.length === 0) return e;
    const X = this.transformMatrix;
    return matmul(matmulAtB(X, e), X);
  }

  /** Abstract methods implemented by derived classes */
  abstract precomputeERI(onProgress?: (percent: number, msg?: string) => void, backend?: EriBackend): Promise<void>;
  abstract guessInitialFockMatrix(): void;
  abstract computeCoefficientMatrix(): void;
  abstract computeDensityMatrix(): void;
  abstract computeFockMatrix(): void | Promise<void>;
  abstract computeEnergy(): number;
  abstract updateFockMatrix(): void;
  abstract getTotalEnergy(): number;

  /** Optional SAD initial density (set before solve) */
  protected initialDensityGuess?: Matrix;
  /** Initial guess type: 'core' (default) | 'gwh' */
  protected initialGuessType: 'core' | 'gwh' = 'core';

  setInitialDensityGuess(density: Matrix) {
    this.initialDensityGuess = density;
  }

  setInitialGuessType(type: 'core' | 'gwh') {
    this.initialGuessType = type;
  }

  /** Override in subclasses to set the density matrix from SAD guess */
  protected setInitialDensity(_density: Matrix): void {
    // default no-op; subclasses override
  }

  /** Build GWH (Generalized Wolfsberg-Helmholz) Fock matrix.
   *  Diagonal: F_ii = H_ii, Off-diagonal: F_ij = K * S_ij * (H_ii + H_jj) / 2 */
  protected buildGWHFock(): Matrix {
    const K = 1.75;
    const n = this.numBasis;
    const F = new Matrix(n, n);
    for (let i = 0; i < n; i++) {
      F.set(i, i, this.coreHamiltonianMatrix.get(i, i));
      for (let j = i + 1; j < n; j++) {
        const val = K * this.overlapMatrix.get(i, j) *
          (this.coreHamiltonianMatrix.get(i, i) + this.coreHamiltonianMatrix.get(j, j)) / 2;
        F.set(i, j, val);
        F.set(j, i, val);
      }
    }
    return F;
  }

  /** GWH initial Fock guess — subclasses override to set their Fock field(s) */
  guessInitialFockMatrixGWH(): void {
    // fallback: use core Hamiltonian (subclasses should override)
    this.guessInitialFockMatrix();
  }

  /** Run the SCF procedure */
  async solve(callbacks?: SCFCallbacks): Promise<number> {
    if (callbacks?.eriBackend) this.eriBackend = callbacks.eriBackend;
    const accelMethod = callbacks?.scfAccelMethod ?? 'diis';
    this.scfAccelerator = createAccelerator(accelMethod, callbacks?.scfAccelParams);

    callbacks?.onPhase?.('integrals', 'start');

    callbacks?.onProgress?.('Computing nuclear repulsion energy...');
    this.computeNuclearRepulsionEnergy();

    callbacks?.onProgress?.('Computing one-electron integrals...');
    this.computeCoreHamiltonianMatrix();

    callbacks?.onProgress?.('Computing two-electron integrals...');
    await this.precomputeERI((p, backendInfo) => {
      if (backendInfo) {
        callbacks?.onProgress?.(`Computing ERIs... ${p.toFixed(0)}% [${backendInfo}]`);
      } else {
        callbacks?.onProgress?.(`Computing ERIs... ${p.toFixed(0)}%`);
      }
    }, callbacks?.eriBackend);

    callbacks?.onProgress?.('Computing transformation matrix...');
    this.computeTransformMatrix();

    if (this.initialDensityGuess) {
      callbacks?.onProgress?.('Setting SAD initial density...');
      this.setInitialDensity(this.initialDensityGuess);
      await this.computeFockMatrix();
    } else if (this.initialGuessType === 'gwh') {
      callbacks?.onProgress?.('Computing GWH initial Fock matrix...');
      this.guessInitialFockMatrixGWH();
    } else {
      callbacks?.onProgress?.('Guessing initial Fock matrix...');
      this.guessInitialFockMatrix();
    }

    callbacks?.onPhase?.('integrals', 'done');
    callbacks?.onPhase?.('scf', 'start');

    let prevEnergy = 0;

    for (let iter = 0; iter < this.maxIter; iter++) {
      this.computeCoefficientMatrix();
      this.computeDensityMatrix();
      await this.computeFockMatrix();
      const energy = this.computeEnergy();

      const deltaE = energy - prevEnergy;
      const totalEnergy = energy + this.nuclearRepulsionEnergy;

      callbacks?.onIteration?.(iter, totalEnergy, deltaE);

      // Switch to target grid when close to convergence (DFT only)
      let switchedThisIter = false;
      if (this._xcFunctional && !this._gridSwitched && iter > 0 && Math.abs(deltaE) < 1e-4) {
        this.switchToTargetGrid();
        switchedThisIter = true;
        callbacks?.onProgress?.(`Switching to target grid (${this._grid!.length} points)`);
      }

      // The energy of the switching iteration was evaluated on the coarse grid, so it
      // must never satisfy the convergence test: small systems (e.g. H2/STO-3G) converge
      // on the coarse grid within the same iteration and would otherwise return without
      // a single target-grid evaluation (1e-4 Eh error, ∫ρ off by 2e-4).
      if (!switchedThisIter && iter > 0 && Math.abs(deltaE) < this.convergenceThreshold) {
        callbacks?.onProgress?.(`Converged after ${iter + 1} iterations`);
        callbacks?.onPhase?.('scf', 'done');
        return totalEnergy;
      }

      this.updateFockMatrix();
      prevEnergy = energy;

      // Yield to browser event loop
      await new Promise<void>(r => setTimeout(r, 0));
    }

    callbacks?.onProgress?.(`Warning: SCF did not converge in ${this.maxIter} iterations`);
    callbacks?.onPhase?.('scf', 'done');
    return this.getTotalEnergy();
  }
}
