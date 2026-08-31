/** CIS (Configuration Interaction Singles) excited state calculation.
 *
 *  CIS matrix: A_{ia,jb} = δ_{ij}δ_{ab}(ε_a - ε_i) + 2(ia|jb) - (ij|ab)   [singlet]
 *              A_{ia,jb} = δ_{ij}δ_{ab}(ε_a - ε_i) - (ij|ab)               [triplet]
 *
 *  Diagonalize A to get excitation energies and transition amplitudes.
 *  Oscillator strengths: f_n = (2/3) ω_n |<0|μ|n>|²
 */

import type { ERIStored } from './eri';
import { Matrix, type FloatArray } from '../linalg/matrix';
import { jacobiEigen } from '../linalg/eigendecomposition';
import { computeDipoleIntegrals } from './integralsDipole';
import type { PrimitiveShell } from './types';

export interface CISExcitedState {
  energy: number;           // excitation energy in Hartree
  energyEV: number;         // excitation energy in eV
  oscillatorStrength: number;
  dominantTransitions: Array<{ i: number; a: number; coeff: number }>;
}

export interface CISResult {
  states: CISExcitedState[];
  isTriplet: boolean;
}

const HA_TO_EV = 27.211386245988;

/** Build half-transformed integrals: half[p][q,λ,σ] = (pq|λσ) in mixed basis.
 *  p runs over orbitals pStart..pStart+np-1, q over qStart..qStart+nq-1. */
function buildHalfTransform(
  C: Matrix, eri: ERIStored, N: number,
  pStart: number, np: number, qStart: number, nq: number,
): Float64Array[] {
  const result: Float64Array[] = new Array(np);
  const buf = new Float64Array(N * nq);

  for (let p = 0; p < np; p++) {
    const tmp = new Float64Array(nq * N * N);
    result[p] = tmp;
    const pIdx = pStart + p;

    for (let lam = 0; lam < N; lam++) {
      for (let sig = 0; sig < N; sig++) {
        buf.fill(0);
        for (let mu = 0; mu < N; mu++) {
          for (let nu = 0; nu < N; nu++) {
            const eriVal = eri.get(mu, nu, lam, sig);
            if (Math.abs(eriVal) < 1e-15) continue;
            for (let q = 0; q < nq; q++) {
              buf[mu * nq + q] += C.get(nu, qStart + q) * eriVal;
            }
          }
        }
        for (let mu = 0; mu < N; mu++) {
          const Cmp = C.get(mu, pIdx);
          if (Math.abs(Cmp) < 1e-15) continue;
          for (let q = 0; q < nq; q++) {
            tmp[q * N * N + lam * N + sig] += Cmp * buf[mu * nq + q];
          }
        }
      }
    }
  }
  return result;
}

/** Compute CIS excited states from converged RHF.
 *  @param C           MO coefficient matrix (AO × MO)
 *  @param epsilon     Orbital energies (sorted)
 *  @param eri         AO ERI storage
 *  @param nocc        Number of occupied orbitals
 *  @param nbasis      Number of basis functions
 *  @param nStates     Number of excited states to return
 *  @param isTriplet   If true, compute triplet states
 *  @param primitiveShells  Primitive shells for dipole integrals
 *  @param normFactors      CGTO normalization factors
 *  @param onProgress  Optional progress callback
 */
export function computeCIS(
  C: Matrix,
  epsilon: FloatArray,
  eri: ERIStored,
  nocc: number,
  nbasis: number,
  nStates: number,
  isTriplet: boolean,
  primitiveShells: PrimitiveShell[],
  normFactors: number[],
  onProgress?: (msg: string) => void,
): CISResult {
  const nvir = nbasis - nocc;
  const dim = nocc * nvir;
  const N = nbasis;

  onProgress?.(`Building CIS matrix (${dim}×${dim})...`);

  // Half-transform: halfOV[i][a,λ,σ] = (ia|λσ), i∈occ, a∈vir
  const halfOV = buildHalfTransform(C, eri, N, 0, nocc, nocc, nvir);

  // Half-transform: halfOO[i][j,λ,σ] = (ij|λσ), i,j∈occ
  const halfOO = buildHalfTransform(C, eri, N, 0, nocc, 0, nocc);

  onProgress?.('Assembling CIS Hamiltonian...');

  // Build CIS matrix
  const A = new Matrix(dim, dim);

  for (let i = 0; i < nocc; i++) {
    for (let a = 0; a < nvir; a++) {
      const ia = i * nvir + a;

      for (let j = 0; j < nocc; j++) {
        for (let b = 0; b < nvir; b++) {
          const jb = j * nvir + b;

          let val = 0.0;

          // Diagonal: (ε_a - ε_i) δ_{ij} δ_{ab}
          if (i === j && a === b) {
            val += epsilon[nocc + a] - epsilon[i];
          }

          // (ia|jb) from halfOV[i]: complete second half-transform
          // (ia|jb) = Σ_{λσ} C[λ,j] * C[σ,b+nocc] * halfOV[i][a,λ,σ]
          if (!isTriplet) {
            let eri_iajb = 0.0;
            for (let lam = 0; lam < N; lam++) {
              const Clj = C.get(lam, j);
              if (Math.abs(Clj) < 1e-15) continue;
              for (let sig = 0; sig < N; sig++) {
                eri_iajb += Clj * C.get(sig, nocc + b) * halfOV[i][a * N * N + lam * N + sig];
              }
            }
            val += 2.0 * eri_iajb;
          }

          // (ij|ab) from halfOO[i]: complete second half-transform
          // (ij|ab) = Σ_{λσ} C[λ,a+nocc] * C[σ,b+nocc] * halfOO[i][j,λ,σ]
          let eri_ijab = 0.0;
          for (let lam = 0; lam < N; lam++) {
            const Cla = C.get(lam, nocc + a);
            if (Math.abs(Cla) < 1e-15) continue;
            for (let sig = 0; sig < N; sig++) {
              eri_ijab += Cla * C.get(sig, nocc + b) * halfOO[i][j * N * N + lam * N + sig];
            }
          }
          val -= eri_ijab;

          A.set(ia, jb, val);
        }
      }
    }
  }

  onProgress?.('Diagonalizing CIS matrix...');

  // Diagonalize
  const { eigenvalues, eigenvectors } = jacobiEigen(A);

  // Compute oscillator strengths (singlet only; triplet transitions are spin-forbidden)
  let oscStrengths: number[];
  if (isTriplet) {
    oscStrengths = new Array(Math.min(nStates, dim)).fill(0.0);
  } else {
    onProgress?.('Computing oscillator strengths...');
    const { Dx, Dy, Dz } = computeDipoleIntegrals(primitiveShells, normFactors, nbasis);

    // MO dipole integrals: μ_{ia} = Σ_{μν} C[μ,i] * C[ν,a+nocc] * D[μ,ν]
    const dipX = new Float64Array(dim);
    const dipY = new Float64Array(dim);
    const dipZ = new Float64Array(dim);

    for (let i = 0; i < nocc; i++) {
      for (let a = 0; a < nvir; a++) {
        const ia = i * nvir + a;
        let dx = 0.0, dy = 0.0, dz = 0.0;
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
        dipX[ia] = dx;
        dipY[ia] = dy;
        dipZ[ia] = dz;
      }
    }

    oscStrengths = [];
    for (let n = 0; n < Math.min(nStates, dim); n++) {
      let tdx = 0.0, tdy = 0.0, tdz = 0.0;
      for (let ia = 0; ia < dim; ia++) {
        const c = eigenvectors.get(ia, n);
        tdx += c * dipX[ia];
        tdy += c * dipY[ia];
        tdz += c * dipZ[ia];
      }
      // √2 factor for RHF (α+β spin contributions)
      tdx *= Math.SQRT2;
      tdy *= Math.SQRT2;
      tdz *= Math.SQRT2;
      const tdm2 = tdx * tdx + tdy * tdy + tdz * tdz;
      const omega = eigenvalues[n];
      const f = (2.0 / 3.0) * omega * tdm2;
      oscStrengths.push(Math.max(0, f));
    }
  }

  // Extract results
  const actualStates = Math.min(nStates, dim);
  const states: CISExcitedState[] = [];

  for (let n = 0; n < actualStates; n++) {
    const transitions: Array<{ i: number; a: number; coeff: number }> = [];
    for (let i = 0; i < nocc; i++) {
      for (let a = 0; a < nvir; a++) {
        const ia = i * nvir + a;
        const c = eigenvectors.get(ia, n);
        if (Math.abs(c) > 0.1) {
          transitions.push({ i, a: a + nocc, coeff: c });
        }
      }
    }
    transitions.sort((x, y) => Math.abs(y.coeff) - Math.abs(x.coeff));

    states.push({
      energy: eigenvalues[n],
      energyEV: eigenvalues[n] * HA_TO_EV,
      oscillatorStrength: oscStrengths[n],
      dominantTransitions: transitions.slice(0, 3),
    });
  }

  onProgress?.(`CIS: ${actualStates} excited states computed`);
  return { states, isTriplet };
}
