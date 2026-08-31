/**
 * TDDFT (Time-Dependent DFT) excited states.
 *   - Tamm-Dancoff Approximation (TDA): solves A·X = ω·X
 *   - Full Casida (RPA): solves the [A,B;-B,-A] super-eigenproblem
 *
 * Casida matrices in the singles space (occ → vir, Mulliken notation):
 *
 *   A_{ia,jb} = δ_{ij}δ_{ab}(ε_a - ε_i)
 *             + 2(ia|jb) - α(ij|ab) - β(ij|ab)_LR + 2(ia|f^S_xc|jb)   [singlet]
 *             - α(ij|ab) - β(ij|ab)_LR + 1(ia|f^T_xc|jb)              [triplet]
 *
 *   B_{ia,jb} = 2(ia|bj) - α(ib|aj) - β(ib|aj)_LR + 2(ia|f^S_xc|jb)   [singlet]
 *             - α(ib|aj) - β(ib|aj)_LR + 1(ia|f^T_xc|jb)              [triplet]
 *
 * The full Casida problem reduces to the Hermitian form (when A-B > 0):
 *   (A-B)^{1/2} (A+B) (A-B)^{1/2} · Q = ω² · Q
 * with X+Y = (A-B)^{-1/2}·Q used for oscillator strengths.
 *
 * For pure DFT (α = 0): the HF-exchange terms vanish; kernel + Coulomb provide response.
 * RSH long-range exchange uses a separate LR ERI store passed in.
 *
 * Implementation notes:
 *   - LDA-like (ALDA) kernel: numerical 2nd derivative of ε_xc·ρ w.r.t. spin densities.
 *   - GGA/meta-GGA gradient/τ kernel pieces are omitted; ALDA captures the dominant part.
 *   - When A-B has negative eigenvalues (triplet instability), Casida gives imaginary ω;
 *     we report |ω| with a sign flag in the message.
 */

import type { ERIStored } from './eri';
import type { PrimitiveShell } from './types';
import type { XCFunctional, XCInput } from './xcFunctional';
import type { GridPoint } from './grid';
import { Matrix, type FloatArray } from '../linalg/matrix';
import { jacobiEigen } from '../linalg/eigendecomposition';
import { computeDipoleIntegrals } from './integralsDipole';
import { buildShellNormCache, evaluateBasisFunctionsScreened } from './xcIntegration';
import { transformERItoMO } from './cphf';
import { transformERItoMOWasm, isTransformERItoMOWasmAvailable } from './eriWasm';
import type { CISExcitedState } from './cis';

const HA_TO_EV = 27.211386245988;

export interface TDDFTResult {
  states: CISExcitedState[];
  isTriplet: boolean;
  method?: 'TDA' | 'Casida';
}

/** Compute LDA-like XC kernel on a grid point.
 *  Returns the singlet kernel: K^S = ∂²(ε·ρ)/∂ρ² evaluated at total density,
 *  and the triplet kernel: K^T = (1/ρ_α) ∂(ε_x · ρ)/∂ρ_α — spin-flip part.
 *  Both via central FD on the functional. */
function xcKernelAtPoint(
  fn: XCFunctional, rho: number, isTriplet: boolean,
): number {
  if (rho < 1e-12) return 0;
  // Closed-shell partition: ρ_α = ρ_β = ρ/2
  const rhoSpin = rho / 2;
  const eps = Math.max(rho * 1e-5, 1e-9);

  if (!isTriplet) {
    // Singlet: f^S = (1/2) (f_αα + f_αβ + f_βα + f_ββ) ≈ d²(ε·ρ)/dρ² at uniform spin
    // Use total density variation: ρ → ρ ± h, ρ_α=ρ_β=ρ/2 always
    const inp = (rho_total: number): XCInput => ({ rhoA: rho_total / 2, rhoB: rho_total / 2 });
    const f = (r: number): number => fn.evaluate(inp(r)).exc * r;
    const f0 = f(rho);
    const fp = f(rho + eps);
    const fm = f(Math.max(rho - eps, 1e-15));
    return (fp - 2 * f0 + fm) / (eps * eps);
  } else {
    // Triplet kernel: f^T = ∂²(ε·ρ)/∂s², s = ρ_α - ρ_β at fixed ρ
    // For closed-shell reference at ρ_α=ρ_β=ρ/2, vary s while keeping ρ fixed
    const f = (s: number): number => {
      const rhoA = rhoSpin + 0.5 * s;
      const rhoB = rhoSpin - 0.5 * s;
      if (rhoA < 0 || rhoB < 0) return 0;
      return fn.evaluate({ rhoA, rhoB }).exc * (rhoA + rhoB);
    };
    const f0 = f(0);
    const fp = f(eps);
    const fm = f(-eps);
    return (fp - 2 * f0 + fm) / (eps * eps);
  }
}

/** Compute TDDFT excited states (TDA or Full Casida).
 *  C, epsilon: KS (or HF) orbitals
 *  eri: full AO ERIs
 *  eriLR: optional long-range AO ERIs (for RSH); null otherwise
 *  xcFunctional: the XC functional (used only for the kernel)
 *  grid: integration grid points
 *  hfFraction: short-range HF coefficient α (= exactExchangeFraction)
 *  rsBeta: additional long-range HF coefficient β (0 for non-RSH)
 *  useFullCasida: if true, build B and solve full Casida equations (default: TDA only) */
export function computeTDDFT(
  C: Matrix,
  epsilon: FloatArray,
  eri: ERIStored,
  eriLR: ERIStored | null,
  xcFunctional: XCFunctional | null,
  grid: GridPoint[] | null,
  shells: PrimitiveShell[],
  normFactors: number[],
  nocc: number,
  nbasis: number,
  nStates: number,
  isTriplet: boolean,
  hfFraction: number,
  rsBeta: number,
  onProgress?: (msg: string) => void,
  useFullCasida: boolean = false,
): TDDFTResult {
  const nvir = nbasis - nocc;
  const dim = nocc * nvir;
  const N = nbasis;
  const tag = useFullCasida ? 'TDDFT-Casida' : 'TDDFT-TDA';

  onProgress?.(`${tag} ${isTriplet ? 'triplet' : 'singlet'}: building A (${dim}×${dim})...`);

  // Pre-compute MO orbital values on the grid (only occupied + virtual we care about)
  let phiMOgrid: Float64Array | null = null;
  let kxcGrid: Float64Array | null = null;

  if (xcFunctional && grid) {
    onProgress?.('TDDFT-TDA: evaluating MOs and XC kernel on grid...');
    phiMOgrid = new Float64Array(grid.length * nbasis); // φ_p^MO at each grid point
    kxcGrid = new Float64Array(grid.length);
    const shellNorms = buildShellNormCache(shells);
    for (let g = 0; g < grid.length; g++) {
      const gp = grid[g];
      const { phi, sigIdx, nsig } = evaluateBasisFunctionsScreened(
        shells, shellNorms, normFactors, nbasis, gp.x, gp.y, gp.z, false,
      );
      // Density at this grid point (closed-shell)
      let rho = 0;
      for (let ii = 0; ii < nsig; ii++) {
        const mu = sigIdx[ii];
        const phiMu = phi[mu];
        const muOff = mu * nbasis;
        // Only occupied orbitals contribute to ρ (closed-shell: ρ = 2 Σ_i |φ_i^MO|²)
        for (let i = 0; i < nocc; i++) {
          let phi_i_mo = 0;
          for (let jj = 0; jj < nsig; jj++) {
            const nu = sigIdx[jj];
            phi_i_mo += C.get(nu, i) * phi[nu];
          }
          // Only need to compute once per grid point per orbital — restructure
          break;
        }
        // We compute MO values per orbital outside this inner loop instead
      }
      // Cleaner approach: compute MO values per orbital
      for (let p = 0; p < nbasis; p++) {
        let mo = 0;
        for (let ii = 0; ii < nsig; ii++) {
          const mu = sigIdx[ii];
          mo += C.get(mu, p) * phi[mu];
        }
        phiMOgrid[g * nbasis + p] = mo;
      }
      // Density from occupied MOs
      rho = 0;
      for (let i = 0; i < nocc; i++) {
        const mo = phiMOgrid[g * nbasis + i];
        rho += 2 * mo * mo;
      }
      kxcGrid[g] = rho > 1e-12 ? xcKernelAtPoint(xcFunctional, rho, isTriplet) : 0;
    }
  }

  // ── Pre-compute MO ERIs (full and long-range) using WASM where available ──
  onProgress?.('TDDFT-TDA: AO → MO ERI transform...');
  const Cflat = new Float64Array(C.data);
  let moERI: Float64Array;
  if (isTransformERItoMOWasmAvailable()) {
    const r = transformERItoMOWasm(Cflat, new Float64Array(eri.data), N);
    moERI = r ?? transformERItoMO(C, eri, N);
  } else {
    moERI = transformERItoMO(C, eri, N);
  }
  let moERILR: Float64Array | null = null;
  if (eriLR && rsBeta !== 0) {
    onProgress?.('TDDFT-TDA: AO → MO long-range ERI transform...');
    if (isTransformERItoMOWasmAvailable()) {
      const r = transformERItoMOWasm(Cflat, new Float64Array(eriLR.data), N);
      moERILR = r ?? transformERItoMO(C, eriLR, N);
    } else {
      moERILR = transformERItoMO(C, eriLR, N);
    }
  }

  // MO-ERI accessor: moERI is laid out as ((p*N + q)*N + r)*N + s = (pq|rs) (Mulliken)
  const N2 = N * N, N3 = N2 * N;
  const moGet = (p: number, q: number, r: number, s: number, table: Float64Array): number =>
    table[p * N3 + q * N2 + r * N + s];

  // Helper: build the A or B matrix in the singles basis.
  // For A: Coulomb (ia|jb), exchange (ij|ab), kernel (ia|f|jb)
  // For B: Coulomb (ia|bj), exchange (ib|aj), kernel (ia|f|jb)  (kernel symmetric in a↔b swap)
  function buildABMatrix(isB: boolean): Matrix {
    const M = new Matrix(dim, dim);
    for (let i = 0; i < nocc; i++) {
      for (let a = 0; a < nvir; a++) {
        const ia = i * nvir + a;
        const aMO = nocc + a;
        for (let j = 0; j < nocc; j++) {
          for (let b = 0; b < nvir; b++) {
            const jb = j * nvir + b;
            const bMO = nocc + b;
            let val = 0;

            // Diagonal (only A)
            if (!isB && i === j && a === b) val += epsilon[aMO] - epsilon[i];

            // Coulomb (singlet only)
            //   A: 2 (ia|jb)_MO,  B: 2 (ia|bj)_MO  (in Mulliken, swap j↔b in 2nd pair)
            if (!isTriplet) {
              val += 2 * (isB ? moGet(i, aMO, bMO, j, moERI) : moGet(i, aMO, j, bMO, moERI));
            }

            // HF exchange
            //   A: -α (ij|ab)_MO - β (ij|ab)_MO_LR
            //   B: -α (ib|aj)_MO - β (ib|aj)_MO_LR
            if (hfFraction !== 0) {
              val -= hfFraction * (isB ? moGet(i, bMO, aMO, j, moERI) : moGet(i, j, aMO, bMO, moERI));
            }
            if (moERILR && rsBeta !== 0) {
              val -= rsBeta * (isB ? moGet(i, bMO, aMO, j, moERILR) : moGet(i, j, aMO, bMO, moERILR));
            }

            // XC kernel (ALDA): same for A and B (kernel is symmetric in (i,a)↔(j,b))
            if (xcFunctional && grid && phiMOgrid && kxcGrid) {
              let kernel = 0;
              for (let g = 0; g < grid.length; g++) {
                const phi_i = phiMOgrid[g * nbasis + i];
                const phi_a = phiMOgrid[g * nbasis + aMO];
                const phi_j = phiMOgrid[g * nbasis + j];
                const phi_b = phiMOgrid[g * nbasis + bMO];
                kernel += grid[g].weight * phi_i * phi_a * kxcGrid[g] * phi_j * phi_b;
              }
              val += (isTriplet ? 1.0 : 2.0) * kernel;
            }

            M.set(ia, jb, val);
          }
        }
      }
    }
    return M;
  }

  const A = buildABMatrix(false);

  // ── Diagonalise: TDA solves A·X=ωX directly; Casida builds B and solves the
  //    Hermitian reduced form M·Q = ω²·Q with M = (A-B)^{1/2}(A+B)(A-B)^{1/2}.
  let eigenvalues: Float64Array;
  let eigenvectors: Matrix;          // for TDA: cols are X; for Casida: cols are X+Y
  let casidaImaginary: number = 0;   // count of states with ω²<0 (instability)

  if (!useFullCasida) {
    onProgress?.(`${tag}: diagonalising (${dim}×${dim})...`);
    const eig = jacobiEigen(A);
    eigenvalues = eig.eigenvalues;
    eigenvectors = eig.eigenvectors;
  } else {
    onProgress?.(`${tag}: building B (${dim}×${dim})...`);
    const B = buildABMatrix(true);
    // (A±B) — both are symmetric
    const APlusB = new Matrix(dim, dim);
    const AMinusB = new Matrix(dim, dim);
    for (let i = 0; i < dim * dim; i++) {
      APlusB.data[i] = A.data[i] + B.data[i];
      AMinusB.data[i] = A.data[i] - B.data[i];
    }
    onProgress?.(`${tag}: eigendecomposing A-B...`);
    const { eigenvalues: dAMB, eigenvectors: U_AMB } = jacobiEigen(AMinusB);
    // Capture instability (ω² < 0 will follow from negative eigenvalues of A-B at ground state)
    let nNeg = 0;
    const sqD = new Float64Array(dim);
    for (let k = 0; k < dim; k++) {
      if (dAMB[k] < 0) { nNeg++; sqD[k] = -Math.sqrt(-dAMB[k]); }
      else sqD[k] = Math.sqrt(dAMB[k]);
    }
    casidaImaginary = nNeg;
    onProgress?.(`${tag}: forming reduced M' = sqD·U^T·(A+B)·U·sqD ...`);
    // tmp = (A+B) · U
    const tmp = new Matrix(dim, dim);
    for (let i = 0; i < dim; i++) for (let j = 0; j < dim; j++) {
      let s = 0; for (let k = 0; k < dim; k++) s += APlusB.data[i*dim+k] * U_AMB.data[k*dim+j];
      tmp.data[i*dim+j] = s;
    }
    // M'_ij = sqD_i · (Σ_k U_ki · tmp_kj) · sqD_j
    const Mp = new Matrix(dim, dim);
    for (let i = 0; i < dim; i++) {
      const si = sqD[i];
      for (let j = 0; j < dim; j++) {
        const sj = sqD[j];
        let s = 0; for (let k = 0; k < dim; k++) s += U_AMB.data[k*dim+i] * tmp.data[k*dim+j];
        Mp.data[i*dim+j] = si * s * sj;
      }
    }
    // Symmetrize Mp (small numerical asymmetry from finite precision)
    for (let i = 0; i < dim; i++) for (let j = i+1; j < dim; j++) {
      const v = 0.5 * (Mp.data[i*dim+j] + Mp.data[j*dim+i]);
      Mp.data[i*dim+j] = v; Mp.data[j*dim+i] = v;
    }

    onProgress?.(`${tag}: diagonalising M' (${dim}×${dim})...`);
    const { eigenvalues: omega2, eigenvectors: Qp } = jacobiEigen(Mp);

    // Excitation energies: ω = sign(ω²) · sqrt(|ω²|)  (negative ω² → instability)
    eigenvalues = new Float64Array(dim);
    for (let n = 0; n < dim; n++) {
      eigenvalues[n] = omega2[n] >= 0 ? Math.sqrt(omega2[n]) : -Math.sqrt(-omega2[n]);
    }

    // Recover (X+Y) eigenvectors in the original singles basis.
    //   M·q = ω²·q with q = (A-B)^{1/2}·(X-Y) and Z·Z = 1 (Jacobi normalization).
    //   Then  (X+Y) = (1/ω)·(A-B)^{1/2}·q.
    //   For Casida normalization <X+Y|X-Y>=1, the proper (X+Y) is rescaled by sqrt(ω).
    //   Storing (X+Y)·sqrt(ω) lets the oscillator strength code share the TDA formula
    //   f = (2/3)·ω·|<vec, μ>|² (with √2 spin factor) without modification.
    eigenvectors = new Matrix(dim, dim);
    for (let n = 0; n < dim; n++) {
      const w = eigenvalues[n];
      const aw = Math.abs(w);
      const sqrtwInv = aw > 1e-10 ? 1.0 / Math.sqrt(aw) : 0;
      for (let i = 0; i < dim; i++) {
        let s = 0;
        for (let k = 0; k < dim; k++) s += U_AMB.data[i*dim+k] * sqD[k] * Qp.data[k*dim+n];
        eigenvectors.data[i*dim+n] = s * sqrtwInv;
      }
    }

    if (casidaImaginary > 0) {
      onProgress?.(`${tag}: ${casidaImaginary} state(s) with imaginary ω (instability)`);
    }
  }

  // Oscillator strengths (singlet only)
  let dipX: Float64Array | null = null, dipY: Float64Array | null = null, dipZ: Float64Array | null = null;
  if (!isTriplet) {
    const { Dx, Dy, Dz } = computeDipoleIntegrals(shells, normFactors, nbasis);
    dipX = new Float64Array(dim); dipY = new Float64Array(dim); dipZ = new Float64Array(dim);
    for (let i = 0; i < nocc; i++) {
      for (let a = 0; a < nvir; a++) {
        const ia = i * nvir + a;
        let dx = 0, dy = 0, dz = 0;
        for (let mu = 0; mu < N; mu++) {
          const Cmi = C.get(mu, i);
          if (Math.abs(Cmi) < 1e-15) continue;
          for (let nu = 0; nu < N; nu++) {
            const w = Cmi * C.get(nu, nocc + a);
            dx += w * Dx.get(mu, nu);
            dy += w * Dy.get(mu, nu);
            dz += w * Dz.get(mu, nu);
          }
        }
        dipX[ia] = dx; dipY[ia] = dy; dipZ[ia] = dz;
      }
    }
  }

  const states: CISExcitedState[] = [];
  const actualStates = Math.min(nStates, dim);
  for (let n = 0; n < actualStates; n++) {
    let f = 0;
    if (!isTriplet && dipX && dipY && dipZ) {
      let tdx = 0, tdy = 0, tdz = 0;
      for (let ia = 0; ia < dim; ia++) {
        const c = eigenvectors.get(ia, n);
        tdx += c * dipX[ia]; tdy += c * dipY[ia]; tdz += c * dipZ[ia];
      }
      tdx *= Math.SQRT2; tdy *= Math.SQRT2; tdz *= Math.SQRT2;
      f = (2 / 3) * eigenvalues[n] * (tdx * tdx + tdy * tdy + tdz * tdz);
    }
    const transitions: Array<{ i: number; a: number; coeff: number }> = [];
    for (let i = 0; i < nocc; i++) {
      for (let a = 0; a < nvir; a++) {
        const c = eigenvectors.get(i * nvir + a, n);
        if (Math.abs(c) > 0.1) transitions.push({ i, a: a + nocc, coeff: c });
      }
    }
    transitions.sort((x, y) => Math.abs(y.coeff) - Math.abs(x.coeff));
    states.push({
      energy: eigenvalues[n],
      energyEV: eigenvalues[n] * HA_TO_EV,
      oscillatorStrength: Math.max(0, f),
      dominantTransitions: transitions.slice(0, 3),
    });
  }

  onProgress?.(`${tag}: ${actualStates} states computed`);
  return { states, isTriplet, method: useFullCasida ? 'Casida' : 'TDA' };
}
