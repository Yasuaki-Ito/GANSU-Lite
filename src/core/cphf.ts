/** Coupled-Perturbed Hartree-Fock (CPHF) solver — ported from GANSU.
 *
 *  Solves A U^x = B^x for each perturbation x using preconditioned CG.
 *
 *  A_{ai,bj} = (ε_a - ε_i) δ_{ab}δ_{ij} + 4(ai|bj) - (ab|ij) - (aj|ib)
 *  B^x_{ai} = F^x_{ai} - ε_i S^x_{ai}   (skeleton Fock/overlap derivatives)
 *
 *  Returns U^x for each perturbation, from which dD/dR can be constructed.
 */

import type { ERIStored } from './eri';
import { Matrix, type FloatArray } from '../linalg/matrix';

/** CPHF operator: A * U where U is a vector of size nocc*nvir.
 *  Uses MO-basis ERIs (full 4-index, precomputed). */
class CPHFOperator {
  private readonly moERI: Float64Array; // (pq|rs) in MO basis, p*nmo³+q*nmo²+r*nmo+s
  private readonly eps: FloatArray;
  private readonly nocc: number;
  private readonly nvir: number;
  private readonly nmo: number;
  private readonly diagonal: Float64Array;

  constructor(moERI: Float64Array, eps: FloatArray, nocc: number, nvir: number, nmo: number) {
    this.moERI = moERI;
    this.eps = eps;
    this.nocc = nocc;
    this.nvir = nvir;
    this.nmo = nmo;

    // Precompute diagonal: ε_a - ε_i
    this.diagonal = new Float64Array(nocc * nvir);
    for (let i = 0; i < nocc; i++)
      for (let a = 0; a < nvir; a++)
        this.diagonal[i * nvir + a] = eps[nocc + a] - eps[i];
  }

  get dimension(): number { return this.nocc * this.nvir; }

  /** Apply A to input vector, store in output. */
  apply(input: Float64Array, output: Float64Array): void {
    const { nocc, nvir, nmo, eps, moERI, diagonal } = this;
    const n = nocc * nvir;

    // Diagonal part: (ε_a - ε_i) U_{ia}
    for (let idx = 0; idx < n; idx++) {
      output[idx] = diagonal[idx] * input[idx];
    }

    // 2-electron part: Σ_{bj} [4(ai|bj) - (ab|ij) - (aj|ib)] U_{bj}
    for (let i = 0; i < nocc; i++) {
      for (let a = 0; a < nvir; a++) {
        const ia = i * nvir + a;
        const aMO = nocc + a;
        let sum = 0;
        for (let j = 0; j < nocc; j++) {
          for (let b = 0; b < nvir; b++) {
            const bMO = nocc + b;
            const Ubj = input[j * nvir + b];
            if (Math.abs(Ubj) < 1e-15) continue;

            // Mulliken notation: (pq|rs) = moERI[p*nmo³+q*nmo²+r*nmo+s]
            const aibj = moERI[((aMO * nmo + i) * nmo + bMO) * nmo + j];
            const abij = moERI[((aMO * nmo + bMO) * nmo + i) * nmo + j];
            const ajib = moERI[((aMO * nmo + j) * nmo + i) * nmo + bMO];

            sum += (4 * aibj - abij - ajib) * Ubj;
          }
        }
        output[ia] += sum;
      }
    }
  }

  /** Preconditioner: M⁻¹ r = r / (ε_a - ε_i) */
  applyPreconditioner(input: Float64Array, output: Float64Array): void {
    for (let idx = 0; idx < this.dimension; idx++) {
      const d = this.diagonal[idx];
      output[idx] = Math.abs(d) > 1e-12 ? input[idx] / d : input[idx];
    }
  }
}

/** Transform AO ERIs to full MO basis (pq|rs).
 *  O(N⁵) — feasible for small basis sets. */
export function transformERItoMO(C: Matrix, eri: ERIStored, n: number): Float64Array {
  const n4 = n * n * n * n;
  const moERI = new Float64Array(n4);

  // Quarter transform: (pν|λσ) = Σ_μ C_{μp} (μν|λσ)
  const half1 = new Float64Array(n4);
  for (let p = 0; p < n; p++)
    for (let nu = 0; nu < n; nu++)
      for (let lam = 0; lam < n; lam++)
        for (let sig = 0; sig < n; sig++) {
          let val = 0;
          for (let mu = 0; mu < n; mu++)
            val += C.get(mu, p) * eri.get(mu, nu, lam, sig);
          half1[((p * n + nu) * n + lam) * n + sig] = val;
        }

  // Second quarter: (pq|λσ) = Σ_ν C_{νq} (pν|λσ)
  const half2 = new Float64Array(n4);
  for (let p = 0; p < n; p++)
    for (let q = 0; q < n; q++)
      for (let lam = 0; lam < n; lam++)
        for (let sig = 0; sig < n; sig++) {
          let val = 0;
          for (let nu = 0; nu < n; nu++)
            val += C.get(nu, q) * half1[((p * n + nu) * n + lam) * n + sig];
          half2[((p * n + q) * n + lam) * n + sig] = val;
        }

  // Third quarter: (pq|rσ) = Σ_λ C_{λr} (pq|λσ)
  const half3 = new Float64Array(n4);
  for (let p = 0; p < n; p++)
    for (let q = 0; q < n; q++)
      for (let r = 0; r < n; r++)
        for (let sig = 0; sig < n; sig++) {
          let val = 0;
          for (let lam = 0; lam < n; lam++)
            val += C.get(lam, r) * half2[((p * n + q) * n + lam) * n + sig];
          half3[((p * n + q) * n + r) * n + sig] = val;
        }

  // Fourth quarter: (pq|rs) = Σ_σ C_{σs} (pq|rσ)
  for (let p = 0; p < n; p++)
    for (let q = 0; q < n; q++)
      for (let r = 0; r < n; r++)
        for (let s = 0; s < n; s++) {
          let val = 0;
          for (let sig = 0; sig < n; sig++)
            val += C.get(sig, s) * half3[((p * n + q) * n + r) * n + sig];
          moERI[((p * n + q) * n + r) * n + s] = val;
        }

  return moERI;
}

/** Solve CPHF equations: A U^x = B^x for all perturbations.
 *  @param moERI Full MO-basis ERIs
 *  @param eps Orbital energies
 *  @param rhs RHS vectors B^x, flat array [nPert × dim]
 *  @param nocc Number of occupied orbitals
 *  @param nvir Number of virtual orbitals
 *  @param nmo Total number of MOs
 *  @param nPert Number of perturbation directions (3N)
 *  @param tol Convergence tolerance
 *  @param maxIter Maximum CG iterations
 *  @returns Solution vectors U^x, flat array [nPert × dim]
 */
export function solveCPHF(
  moERI: Float64Array, eps: FloatArray,
  rhs: Float64Array,
  nocc: number, nvir: number, nmo: number,
  nPert: number,
  tol: number = 1e-8,
  maxIter: number = 200,
  onProgress?: (msg: string) => void,
): Float64Array {
  const op = new CPHFOperator(moERI, eps, nocc, nvir, nmo);
  const dim = op.dimension;
  const U = new Float64Array(nPert * dim);

  const r = new Float64Array(dim);
  const z = new Float64Array(dim);
  const p = new Float64Array(dim);
  const Ap = new Float64Array(dim);

  for (let pert = 0; pert < nPert; pert++) {
    const b = rhs.subarray(pert * dim, (pert + 1) * dim);
    const x = U.subarray(pert * dim, (pert + 1) * dim);

    // x = 0, r = b
    x.fill(0);
    r.set(b);

    // z = M⁻¹ r
    op.applyPreconditioner(r, z);

    // p = z
    p.set(z);

    // rz = r·z
    let rz = dot(r, z, dim);
    const bNorm = norm(r, dim);
    if (bNorm < 1e-15) continue;

    for (let iter = 0; iter < maxIter; iter++) {
      // Ap = A p
      op.apply(p, Ap);

      // α = rz / (p·Ap)
      const pAp = dot(p, Ap, dim);
      const alpha = rz / (Math.abs(pAp) > 1e-30 ? pAp : 1e-30);

      // x += α p, r -= α Ap
      for (let i = 0; i < dim; i++) { x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; }

      // Convergence check
      const rNorm = norm(r, dim);
      if (rNorm / bNorm < tol) {
        if (pert === 0) onProgress?.(`CPHF: pert 0 converged in ${iter + 1} iters`);
        break;
      }

      // z = M⁻¹ r
      op.applyPreconditioner(r, z);

      // β = (r·z_new) / rz_old
      const rzNew = dot(r, z, dim);
      const beta = rzNew / (Math.abs(rz) > 1e-30 ? rz : 1e-30);

      // p = z + β p
      for (let i = 0; i < dim; i++) p[i] = z[i] + beta * p[i];

      rz = rzNew;
    }
  }

  onProgress?.(`CPHF: ${nPert} perturbations solved`);
  return U;
}

function dot(a: Float64Array, b: Float64Array, n: number): number {
  let s = 0;
  for (let i = 0; i < n; i++) s += a[i] * b[i];
  return s;
}

function norm(a: Float64Array, n: number): number {
  return Math.sqrt(dot(a, a, n));
}
