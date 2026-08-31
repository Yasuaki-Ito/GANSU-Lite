/** RI-MP2: Resolution of Identity approximation for MP2.
 *
 *  (ia|jb) ≈ Σ_P B_ia^P * B_jb^P
 *
 *  E_MP2 = Σ_{i<j,a<b} (ia|jb)[2(ia|jb) - (ib|ja)] / (ε_i + ε_j - ε_a - ε_b)
 *
 *  Cost: O(N^2 * N_aux * N_occ) for MO transform, O(N_occ^2 * N_vir^2 * N_aux) for energy.
 */

import type { FloatArray } from '../linalg/matrix';
import type { Matrix } from '../linalg/matrix';
import type { RIData } from './ri';
import { computeRIMP2EnergyWasm, computeRIUMP2EnergyWasm } from './eriWasm';

/** Compute RI-MP2 correlation energy from converged RHF results. */
export function computeRIMP2Energy(
  C: Matrix,
  epsilon: FloatArray,
  ri: RIData,
  nocc: number,
  nbasis: number,
): number {
  // Try WASM path
  const wasmResult = computeRIMP2EnergyWasm(
    ri.bMatrix, C.data, epsilon, ri.naux, nocc, nbasis,
  );
  if (wasmResult !== null) return wasmResult;

  // JS fallback
  const nvir = nbasis - nocc;
  const naux = ri.naux;

  // B_ia^P = Σ_μν C[μ,i] C[ν,a] B_P^{μν}
  const bmo = ri.transformToMO(C, nocc);

  // E_MP2 = Σ_{ij} Σ_{ab} (ia|jb) [2(ia|jb) - (ib|ja)] / (ε_i + ε_j - ε_a - ε_b)
  let emp2 = 0;

  for (let i = 0; i < nocc; i++) {
    for (let j = i; j < nocc; j++) {
      for (let a = 0; a < nvir; a++) {
        for (let b = 0; b < nvir; b++) {
          // (ia|jb) = Σ_P B_ia^P * B_jb^P
          let iajb = 0;
          const iaOff = i * nvir * naux + a * naux;
          const jbOff = j * nvir * naux + b * naux;
          for (let P = 0; P < naux; P++) {
            iajb += bmo[iaOff + P] * bmo[jbOff + P];
          }

          // (ib|ja) = Σ_P B_ib^P * B_ja^P
          let ibja = 0;
          const ibOff = i * nvir * naux + b * naux;
          const jaOff = j * nvir * naux + a * naux;
          for (let P = 0; P < naux; P++) {
            ibja += bmo[ibOff + P] * bmo[jaOff + P];
          }

          const denom = epsilon[i] + epsilon[j] - epsilon[nocc + a] - epsilon[nocc + b];
          const factor = (i === j) ? 1.0 : 2.0;
          emp2 += factor * iajb * (2 * iajb - ibja) / denom;
        }
      }
    }
  }

  return emp2;
}

/** Compute RI-UMP2 correlation energy from converged UHF results. */
export function computeRIUMP2Energy(
  Ca: Matrix, Cb: Matrix,
  epsilonA: FloatArray, epsilonB: FloatArray,
  ri: RIData,
  noccA: number, noccB: number,
  nbasis: number,
): number {
  // Try WASM path
  const wasmResult = computeRIUMP2EnergyWasm(
    ri.bMatrix, Ca.data, Cb.data, epsilonA, epsilonB,
    ri.naux, noccA, noccB, nbasis,
  );
  if (wasmResult !== null) return wasmResult;

  // JS fallback
  const nvirA = nbasis - noccA;
  const nvirB = nbasis - noccB;
  const naux = ri.naux;

  const bmoA = ri.transformToMO(Ca, noccA);
  const bmoB = ri.transformToMO(Cb, noccB);

  let emp2 = 0;

  // αα contribution
  for (let i = 0; i < noccA; i++) {
    for (let j = i + 1; j < noccA; j++) {
      for (let a = 0; a < nvirA; a++) {
        for (let b = a + 1; b < nvirA; b++) {
          let iajb = 0, ibja = 0;
          const iaOff = i * nvirA * naux + a * naux;
          const jbOff = j * nvirA * naux + b * naux;
          const ibOff = i * nvirA * naux + b * naux;
          const jaOff = j * nvirA * naux + a * naux;
          for (let P = 0; P < naux; P++) {
            iajb += bmoA[iaOff + P] * bmoA[jbOff + P];
            ibja += bmoA[ibOff + P] * bmoA[jaOff + P];
          }
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
          let iajb = 0, ibja = 0;
          const iaOff = i * nvirB * naux + a * naux;
          const jbOff = j * nvirB * naux + b * naux;
          const ibOff = i * nvirB * naux + b * naux;
          const jaOff = j * nvirB * naux + a * naux;
          for (let P = 0; P < naux; P++) {
            iajb += bmoB[iaOff + P] * bmoB[jbOff + P];
            ibja += bmoB[ibOff + P] * bmoB[jaOff + P];
          }
          const antisym = iajb - ibja;
          const denom = epsilonB[i] + epsilonB[j] - epsilonB[noccB + a] - epsilonB[noccB + b];
          emp2 += antisym * antisym / denom;
        }
      }
    }
  }

  // αβ contribution
  for (let i = 0; i < noccA; i++) {
    for (let j = 0; j < noccB; j++) {
      for (let a = 0; a < nvirA; a++) {
        for (let b = 0; b < nvirB; b++) {
          let iajb = 0;
          const iaOff = i * nvirA * naux + a * naux;
          const jbOff = j * nvirB * naux + b * naux;
          for (let P = 0; P < naux; P++) {
            iajb += bmoA[iaOff + P] * bmoB[jbOff + P];
          }
          const denom = epsilonA[i] + epsilonB[j] - epsilonA[noccA + a] - epsilonB[noccB + b];
          emp2 += iajb * iajb / denom;
        }
      }
    }
  }

  return emp2;
}
