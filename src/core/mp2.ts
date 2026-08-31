/** MP2 (second-order Moller-Plesset perturbation theory) correlation energy.
 *
 *  E_MP2 = Σ_{i<j,a<b} (ia|jb)[2(ia|jb) - (ib|ja)] / (ε_i + ε_j - ε_a - ε_b)
 *
 *  MO ERI via two half-transformations: O(N^5) total.
 *
 *  Step 1 (first half): for each i, build tmp_i[a,λ,σ] = Σ_μν C[μ,i]*C[ν,a]*(μν|λσ)
 *    Decomposed as:
 *      B[μ,a,λ,σ] = Σ_ν C[ν,a+nocc] * (μν|λσ)     ... O(N * nvir * N^2) per μ
 *      tmp_i[a,λ,σ] = Σ_μ C[μ,i] * B[μ,a,λ,σ]      ... O(N * nvir * N^2)
 *    Total step 1: O(nocc * N^2 * N * nvir) = O(N^5)
 *
 *  Step 2 (second half): (ia|jb) = Σ_λσ C[λ,j] * C[σ,b+nocc] * tmp_i[a,λ,σ]
 *    For each (i,j) pair, contract tmp_i and tmp_j to get (ia|jb) for all a,b
 *    Total step 2: O(nocc^2 * nvir^2 * N^2) ≈ O(N^4) to O(N^5)
 */

import type { ERIStored } from './eri';
import { Matrix, type FloatArray } from '../linalg/matrix';
import { matmul, matmulAtB } from '../linalg/matmul';
import { jacobiEigen } from '../linalg/eigendecomposition';

/** Compute MP2 correlation energy from converged RHF results.
 *  @param C       MO coefficient matrix (AO x MO), columns sorted by orbital energy
 *  @param epsilon Orbital energies (sorted)
 *  @param eri     AO basis ERI storage
 *  @param nocc    Number of occupied orbitals
 *  @param nbasis  Number of basis functions
 */
export function computeMP2Energy(
  C: Matrix,
  epsilon: FloatArray,
  eri: ERIStored,
  nocc: number,
  nbasis: number,
): number {
  const nvir = nbasis - nocc;
  const N = nbasis;

  // Step 1: Build half-transformed integrals for each occupied orbital i.
  // half[i] is a Float64Array of size [nvir × N × N], storing tmp_i[a, λ, σ].
  //
  // tmp_i[a,λ,σ] = Σ_μ C[μ,i] * Σ_ν C[ν,a+nocc] * (μν|λσ)
  //
  // We decompose into two steps to avoid O(N^4) per (i,a,λ,σ):
  //   For each (λ,σ): buf[μ,a] = Σ_ν C[ν,a+nocc] * (μν|λσ)  ... O(N^2 * nvir * N) total over all λ,σ
  //                    tmp_i[a,λ,σ] = Σ_μ C[μ,i] * buf[μ,a]   ... O(nocc * nvir * N) per (λ,σ)

  const halfTransformed: Float64Array[] = new Array(nocc);
  const buf = new Float64Array(N * nvir); // buf[mu * nvir + a]

  for (let i = 0; i < nocc; i++) {
    const tmp = new Float64Array(nvir * N * N); // tmp[a * N*N + lam * N + sig]
    halfTransformed[i] = tmp;

    for (let lam = 0; lam < N; lam++) {
      for (let sig = 0; sig < N; sig++) {
        // First quarter-transform: buf[mu, a] = Σ_ν C[ν, a+nocc] * (μν|λσ)
        buf.fill(0);
        for (let mu = 0; mu < N; mu++) {
          for (let nu = 0; nu < N; nu++) {
            const eriVal = eri.get(mu, nu, lam, sig);
            if (Math.abs(eriVal) < 1e-15) continue;
            for (let a = 0; a < nvir; a++) {
              buf[mu * nvir + a] += C.get(nu, nocc + a) * eriVal;
            }
          }
        }

        // Second quarter-transform: tmp[a, lam, sig] = Σ_μ C[μ,i] * buf[μ, a]
        for (let mu = 0; mu < N; mu++) {
          const Cmi = C.get(mu, i);
          if (Math.abs(Cmi) < 1e-15) continue;
          for (let a = 0; a < nvir; a++) {
            tmp[a * N * N + lam * N + sig] += Cmi * buf[mu * nvir + a];
          }
        }
      }
    }
  }

  // Step 2: Contract half-transformed integrals to get (ia|jb) and accumulate MP2 energy.
  // (ia|jb) = Σ_λσ C[λ,j] * C[σ,b+nocc] * half[i][a, λ, σ]
  //
  // Optimization: precompute second half for j as well:
  //   half[j][b, λ, σ] is already computed, so:
  //   (ia|jb) = Σ_λσ half_i[a,λ,σ] * ... but half_i already has the bra indices transformed.
  //   We need: (ia|jb) = Σ_λσ tmp_i[a,λ,σ] * C[λ,j] * C[σ,b+nocc]
  //
  // Actually, tmp_i[a,λ,σ] has bra fully transformed (i,a) and ket in AO (λ,σ).
  // So (ia|jb) = Σ_λ C[λ,j] * Σ_σ C[σ,b+nocc] * tmp_i[a,λ,σ]

  let emp2 = 0;

  for (let i = 0; i < nocc; i++) {
    const tmp_i = halfTransformed[i];
    for (let j = i; j < nocc; j++) {
      // For this (i,j) pair, compute all (ia|jb) and (ib|ja)
      for (let a = 0; a < nvir; a++) {
        for (let b = 0; b < nvir; b++) {
          // (ia|jb) = Σ_λσ C[λ,j] * C[σ,b+nocc] * tmp_i[a,λ,σ]
          let iajb = 0;
          for (let lam = 0; lam < N; lam++) {
            const Clj = C.get(lam, j);
            if (Math.abs(Clj) < 1e-15) continue;
            for (let sig = 0; sig < N; sig++) {
              iajb += Clj * C.get(sig, nocc + b) * tmp_i[a * N * N + lam * N + sig];
            }
          }

          // (ib|ja) = Σ_λσ C[λ,j] * C[σ,a+nocc] * tmp_i[b,λ,σ]
          // Note: tmp_i[b,...] already has the i index baked in
          let ibja = 0;
          for (let lam = 0; lam < N; lam++) {
            const Clj = C.get(lam, j);
            if (Math.abs(Clj) < 1e-15) continue;
            for (let sig = 0; sig < N; sig++) {
              ibja += Clj * C.get(sig, nocc + a) * tmp_i[b * N * N + lam * N + sig];
            }
          }

          const aa = nocc + a;
          const bb = nocc + b;
          const denom = epsilon[i] + epsilon[j] - epsilon[aa] - epsilon[bb];
          const factor = (i === j) ? 1.0 : 2.0; // account for i<j symmetry
          emp2 += factor * iajb * (2 * iajb - ibja) / denom;
        }
      }
    }
  }

  return emp2;
}

// ---------------------------------------------------------------------------
// UMP2 / ROMP2
// ---------------------------------------------------------------------------

/** Half-transform AO ERIs with coefficients Cocc (occupied) and Cvir (virtual).
 *  Returns half[i][a, λ, σ] for each occupied orbital i.
 *  bra is transformed: (i_occ, a_vir | λ_AO, σ_AO)
 */
function halfTransform(
  Cocc: Matrix, noccOcc: number,
  Cvir: Matrix, noccVir: number, nvirVir: number,
  eri: ERIStored, nbasis: number,
): Float64Array[] {
  const N = nbasis;
  const half: Float64Array[] = new Array(noccOcc);
  const buf = new Float64Array(N * nvirVir);

  for (let i = 0; i < noccOcc; i++) {
    const tmp = new Float64Array(nvirVir * N * N);
    half[i] = tmp;

    for (let lam = 0; lam < N; lam++) {
      for (let sig = 0; sig < N; sig++) {
        buf.fill(0);
        for (let mu = 0; mu < N; mu++) {
          for (let nu = 0; nu < N; nu++) {
            const eriVal = eri.get(mu, nu, lam, sig);
            if (Math.abs(eriVal) < 1e-15) continue;
            for (let a = 0; a < nvirVir; a++) {
              buf[mu * nvirVir + a] += Cvir.get(nu, noccVir + a) * eriVal;
            }
          }
        }
        for (let mu = 0; mu < N; mu++) {
          const Cmi = Cocc.get(mu, i);
          if (Math.abs(Cmi) < 1e-15) continue;
          for (let a = 0; a < nvirVir; a++) {
            tmp[a * N * N + lam * N + sig] += Cmi * buf[mu * nvirVir + a];
          }
        }
      }
    }
  }
  return half;
}

/** Contract half-transformed integral (ia|λσ) with Cj and Cb to get (ia|jb).
 *  Returns the integral value.
 */
function contractMOInt(
  halfI: Float64Array, a: number,
  Cj: Matrix, j: number,
  Cb: Matrix, bOffset: number, b: number,
  N: number,
): number {
  let val = 0;
  for (let lam = 0; lam < N; lam++) {
    const clj = Cj.get(lam, j);
    if (Math.abs(clj) < 1e-15) continue;
    for (let sig = 0; sig < N; sig++) {
      val += clj * Cb.get(sig, bOffset + b) * halfI[a * N * N + lam * N + sig];
    }
  }
  return val;
}

/** Compute UMP2 correlation energy from converged UHF results.
 *  E = E_αα + E_ββ + E_αβ
 */
export function computeUMP2Energy(
  Ca: Matrix, Cb: Matrix,
  epsilonA: FloatArray, epsilonB: FloatArray,
  eri: ERIStored,
  noccA: number, noccB: number,
  nbasis: number,
): number {
  const N = nbasis;
  const nvirA = N - noccA;
  const nvirB = N - noccB;

  // Half-transform: bra=(i_α, a_α), ket=AO
  const halfA = halfTransform(Ca, noccA, Ca, noccA, nvirA, eri, N);
  // Half-transform: bra=(i_β, a_β), ket=AO
  const halfB = halfTransform(Cb, noccB, Cb, noccB, nvirB, eri, N);

  let emp2 = 0;

  // αα contribution: E = Σ_{i<j} Σ_{a<b} |(ia|jb)-(ib|ja)|² / denom
  for (let i = 0; i < noccA; i++) {
    for (let j = i + 1; j < noccA; j++) {
      for (let a = 0; a < nvirA; a++) {
        for (let b = a + 1; b < nvirA; b++) {
          const iajb = contractMOInt(halfA[i], a, Ca, j, Ca, noccA, b, N);
          const ibja = contractMOInt(halfA[i], b, Ca, j, Ca, noccA, a, N);
          const antisym = iajb - ibja;
          const denom = epsilonA[i] + epsilonA[j] - epsilonA[noccA + a] - epsilonA[noccA + b];
          emp2 += antisym * antisym / denom;
        }
      }
    }
  }

  // ββ contribution
  for (let i = 0; i < noccB; i++) {
    for (let j = i + 1; j < noccB; j++) {
      for (let a = 0; a < nvirB; a++) {
        for (let b = a + 1; b < nvirB; b++) {
          const iajb = contractMOInt(halfB[i], a, Cb, j, Cb, noccB, b, N);
          const ibja = contractMOInt(halfB[i], b, Cb, j, Cb, noccB, a, N);
          const antisym = iajb - ibja;
          const denom = epsilonB[i] + epsilonB[j] - epsilonB[noccB + a] - epsilonB[noccB + b];
          emp2 += antisym * antisym / denom;
        }
      }
    }
  }

  // αβ contribution: E = Σ_{iα,jβ,aα,bβ} |(iα aα|jβ bβ)|² / denom
  // half_iα has bra=(iα, aα) with Cα, ket still in AO.
  // Need to contract ket with Cβ for j and b.
  for (let i = 0; i < noccA; i++) {
    for (let j = 0; j < noccB; j++) {
      for (let a = 0; a < nvirA; a++) {
        for (let b = 0; b < nvirB; b++) {
          const iajb = contractMOInt(halfA[i], a, Cb, j, Cb, noccB, b, N);
          const denom = epsilonA[i] + epsilonB[j] - epsilonA[noccA + a] - epsilonB[noccB + b];
          emp2 += iajb * iajb / denom;
        }
      }
    }
  }

  return emp2;
}

/** Compute ROMP2 (semicanonical) correlation energy from converged ROHF results.
 *  1. Build semicanonical orbitals by block-diagonalizing Fα and Fβ in MO basis
 *  2. Call UMP2 with the resulting coefficients and orbital energies
 */
export function computeROMP2Energy(
  C: Matrix,
  Fa: Matrix, Fb: Matrix,
  eri: ERIStored,
  noccA: number, noccB: number,
  nbasis: number,
): number {
  const N = nbasis;

  // Transform Fock matrices to MO basis: F_MO = C^T F C
  const FaMO = matmul(matmulAtB(C, Fa), C);
  const FbMO = matmul(matmulAtB(C, Fb), C);

  // 2-block semicanonical: α occ=[0,noccA) vir=[noccA,N), β occ=[0,noccB) vir=[noccB,N)
  function semicanonicalize(FMO: Matrix, nocc: number): { U: Matrix; epsilon: Float64Array } {
    const nvir = N - nocc;
    const U = Matrix.identity(N);
    const eps = new Float64Array(N);

    // Occupied block
    if (nocc > 0) {
      const occ = new Matrix(nocc, nocc);
      for (let i = 0; i < nocc; i++)
        for (let j = 0; j < nocc; j++)
          occ.set(i, j, FMO.get(i, j));
      const eig = jacobiEigen(occ);
      for (let i = 0; i < nocc; i++) {
        eps[i] = eig.eigenvalues[i];
        for (let j = 0; j < nocc; j++)
          U.set(i, j, eig.eigenvectors.get(i, j));
      }
    }

    // Virtual block
    if (nvir > 0) {
      const vir = new Matrix(nvir, nvir);
      for (let i = 0; i < nvir; i++)
        for (let j = 0; j < nvir; j++)
          vir.set(i, j, FMO.get(nocc + i, nocc + j));
      const eig = jacobiEigen(vir);
      for (let i = 0; i < nvir; i++) {
        eps[nocc + i] = eig.eigenvalues[i];
        for (let j = 0; j < nvir; j++)
          U.set(nocc + i, nocc + j, eig.eigenvectors.get(i, j));
      }
    }

    return { U, epsilon: eps };
  }

  const { U: Ua, epsilon: epsA } = semicanonicalize(FaMO, noccA);
  const { U: Ub, epsilon: epsB } = semicanonicalize(FbMO, noccB);

  // Semicanonical coefficients: C_sc = C * U
  const CaScf = matmul(C, Ua);
  const CbScf = matmul(C, Ub);

  return computeUMP2Energy(CaScf, CbScf, epsA, epsB, eri, noccA, noccB, N);
}
