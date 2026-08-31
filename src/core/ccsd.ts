/** CCSD (Coupled Cluster Singles and Doubles) correlation energy.
 *
 *  RHF: Spatial-orbital formulation (Crawdad notation, ~16x faster than spin-orbital).
 *  UHF: Spin-orbital with blocked ordering [occα, occβ, virα, virβ].
 *  ROHF: Semicanonical orbitals → UHF solver.
 *
 *  Reference: Stanton, Gauss, Watts, Bartlett (JCP 94, 4334, 1991).
 */

import type { ERIStored } from './eri';
import { Matrix, type FloatArray } from '../linalg/matrix';
import { DIIS } from './diis';
import { jacobiEigen } from '../linalg/eigendecomposition';
import { matmul, matmulAtB } from '../linalg/matmul';

// ---------------------------------------------------------------------------
// Full 4-index AO → MO integral transform
// ---------------------------------------------------------------------------

/** 4-step quarter-transform AO ERIs → spatial MO integrals (pq|rs).
 *  Returns flat Float64Array indexed as [p*N³ + q*N² + r*N + s]. */
export function fullMOTransform(C1: Matrix, C2: Matrix, C3: Matrix, C4: Matrix,
                         eri: ERIStored, N: number): Float64Array {
  const N2 = N * N, N3 = N2 * N;
  // Step 1: (pν|λσ) = Σ_μ C1[μ,p] * (μν|λσ)
  const step1 = new Float64Array(N * N3);
  for (let p = 0; p < N; p++) {
    for (let nu = 0; nu < N; nu++) {
      for (let lam = 0; lam < N; lam++) {
        for (let sig = 0; sig < N; sig++) {
          let val = 0;
          for (let mu = 0; mu < N; mu++) {
            const c = C1.get(mu, p);
            if (Math.abs(c) < 1e-15) continue;
            val += c * eri.get(mu, nu, lam, sig);
          }
          step1[p * N3 + nu * N2 + lam * N + sig] = val;
        }
      }
    }
  }
  // Step 2: (pq|λσ) = Σ_ν C2[ν,q] * (pν|λσ)
  const step2 = new Float64Array(N * N3);
  for (let p = 0; p < N; p++) {
    for (let q = 0; q < N; q++) {
      for (let lam = 0; lam < N; lam++) {
        for (let sig = 0; sig < N; sig++) {
          let val = 0;
          for (let nu = 0; nu < N; nu++) {
            const c = C2.get(nu, q);
            if (Math.abs(c) < 1e-15) continue;
            val += c * step1[p * N3 + nu * N2 + lam * N + sig];
          }
          step2[p * N3 + q * N2 + lam * N + sig] = val;
        }
      }
    }
  }
  // Step 3: (pq|rσ) = Σ_λ C3[λ,r] * (pq|λσ)
  const step3 = new Float64Array(N * N3);
  for (let p = 0; p < N; p++) {
    for (let q = 0; q < N; q++) {
      for (let r = 0; r < N; r++) {
        for (let sig = 0; sig < N; sig++) {
          let val = 0;
          for (let lam = 0; lam < N; lam++) {
            const c = C3.get(lam, r);
            if (Math.abs(c) < 1e-15) continue;
            val += c * step2[p * N3 + q * N2 + lam * N + sig];
          }
          step3[p * N3 + q * N2 + r * N + sig] = val;
        }
      }
    }
  }
  // Step 4: (pq|rs) = Σ_σ C4[σ,s] * (pq|rσ)
  const moInts = new Float64Array(N * N3);
  for (let p = 0; p < N; p++) {
    for (let q = 0; q < N; q++) {
      for (let r = 0; r < N; r++) {
        for (let s = 0; s < N; s++) {
          let val = 0;
          for (let sig = 0; sig < N; sig++) {
            const c = C4.get(sig, s);
            if (Math.abs(c) < 1e-15) continue;
            val += c * step3[p * N3 + q * N2 + r * N + sig];
          }
          moInts[p * N3 + q * N2 + r * N + s] = val;
        }
      }
    }
  }
  return moInts;
}

// ---------------------------------------------------------------------------
// T2 indexing helper
// ---------------------------------------------------------------------------
export function t2idx(i: number, j: number, a: number, b: number, nocc: number, nvir: number): number {
  return ((i * nocc + j) * nvir + a) * nvir + b;
}

// ---------------------------------------------------------------------------
// RHF spatial-orbital CCSD (Crawdad notation)
// ---------------------------------------------------------------------------
// v(p,q,r,s) = (pr|qs)   — Coulomb integral
// w(p,q,r,s) = 2*v(p,q,r,s) - v(p,q,s,r)  — antisymmetrized
// P(ia,jb) f = f(i,a,j,b) + f(j,b,i,a)

function ccsdSpatialLoop(
  moInts: Float64Array,
  eps: FloatArray,
  nocc: number, N: number,
  onProgress?: (msg: string) => void,
): number {
  const nvir = N - nocc;
  const N2 = N * N, N3 = N2 * N;

  // v(p,q,r,s) = (pr|qs) = moInts[p*N³ + r*N² + q*N + s]
  const v = (p: number, q: number, r: number, s: number): number =>
    moInts[p * N3 + r * N2 + q * N + s];
  const w = (p: number, q: number, r: number, s: number): number =>
    2.0 * v(p, q, r, s) - v(p, q, s, r);

  const T2 = (i: number, j: number, a: number, b: number): number =>
    ((i * nocc + j) * nvir + a) * nvir + b;

  // Denominators
  const t1Size = nocc * nvir;
  const t2Size = nocc * nocc * nvir * nvir;

  const Dia = new Float64Array(t1Size);
  for (let i = 0; i < nocc; i++)
    for (let a = 0; a < nvir; a++)
      Dia[i * nvir + a] = eps[i] - eps[nocc + a];

  const Dijab = new Float64Array(t2Size);
  for (let i = 0; i < nocc; i++)
    for (let j = 0; j < nocc; j++)
      for (let a = 0; a < nvir; a++)
        for (let b = 0; b < nvir; b++)
          Dijab[T2(i, j, a, b)] = eps[i] + eps[j] - eps[nocc + a] - eps[nocc + b];

  // MP2 initial guess: T1=0, t2(i,j,a,b) = v(i,j,A,B) / D
  let t1 = new Float64Array(t1Size);
  let t2 = new Float64Array(t2Size);
  for (let i = 0; i < nocc; i++)
    for (let j = 0; j < nocc; j++)
      for (let a = 0; a < nvir; a++)
        for (let b = 0; b < nvir; b++) {
          const idx = T2(i, j, a, b);
          t2[idx] = v(i, j, nocc + a, nocc + b) / Dijab[idx];
        }

  // Energy: E = Σ_{ijab} w(i,j,A,B) * (t2(i,j,a,b) + t1(i,a)*t1(j,b))
  function energy(): number {
    let E = 0;
    for (let i = 0; i < nocc; i++)
      for (let j = 0; j < nocc; j++)
        for (let a = 0; a < nvir; a++)
          for (let b = 0; b < nvir; b++) {
            const A = nocc + a, B = nocc + b;
            E += w(i, j, A, B) * (t2[T2(i, j, a, b)] + t1[i * nvir + a] * t1[j * nvir + b]);
          }
    return E;
  }

  let Ecc = energy();
  onProgress?.(`CCSD iter  0: E = ${Ecc.toFixed(12)} (MP2 initial guess)`);

  const diis = new DIIS(8, 4);
  const MAX_ITER = 100;
  const CONV = 1e-10;

  for (let iter = 1; iter <= MAX_ITER; iter++) {

    // ---- F intermediates ----
    // Fkc[k,c] = Σ_{ld} w(k,l,C,D) * t1(l,d)
    const Fkc = new Float64Array(nocc * nvir);
    for (let k = 0; k < nocc; k++)
      for (let c = 0; c < nvir; c++) {
        let val = 0;
        for (let l = 0; l < nocc; l++)
          for (let d = 0; d < nvir; d++)
            val += w(k, l, nocc + c, nocc + d) * t1[l * nvir + d];
        Fkc[k * nvir + c] = val;
      }

    // Fki[k,i] = Σ_{lcd} w(k,l,C,D) * (t2(i,l,c,d) + t1(i,c)*t1(l,d))
    const Fki = new Float64Array(nocc * nocc);
    for (let k = 0; k < nocc; k++)
      for (let i = 0; i < nocc; i++) {
        let val = 0;
        for (let l = 0; l < nocc; l++)
          for (let c = 0; c < nvir; c++)
            for (let d = 0; d < nvir; d++) {
              const ww = w(k, l, nocc + c, nocc + d);
              val += ww * (t2[T2(i, l, c, d)] + t1[i * nvir + c] * t1[l * nvir + d]);
            }
        Fki[k * nocc + i] = val;
      }

    // Fac[a,c] = -Σ_{kld} w(k,l,C,D) * (t2(k,l,a,d) + t1(k,a)*t1(l,d))
    const Fac = new Float64Array(nvir * nvir);
    for (let a = 0; a < nvir; a++)
      for (let c = 0; c < nvir; c++) {
        let val = 0;
        for (let k = 0; k < nocc; k++)
          for (let l = 0; l < nocc; l++)
            for (let d = 0; d < nvir; d++) {
              const ww = w(k, l, nocc + c, nocc + d);
              val -= ww * (t2[T2(k, l, a, d)] + t1[k * nvir + a] * t1[l * nvir + d]);
            }
        Fac[a * nvir + c] = val;
      }

    // ---- L intermediates ----
    // Lki[k,i] = Fki[k,i] + Σ_{lc} w(l,k,C,i) * t1(l,c)
    const Lki = new Float64Array(nocc * nocc);
    for (let k = 0; k < nocc; k++)
      for (let i = 0; i < nocc; i++) {
        let val = Fki[k * nocc + i];
        for (let l = 0; l < nocc; l++)
          for (let c = 0; c < nvir; c++)
            val += w(l, k, nocc + c, i) * t1[l * nvir + c];
        Lki[k * nocc + i] = val;
      }

    // Lac[a,c] = Fac[a,c] + Σ_{kd} w(k,A,D,C) * t1(k,d)
    const Lac = new Float64Array(nvir * nvir);
    for (let a = 0; a < nvir; a++)
      for (let c = 0; c < nvir; c++) {
        let val = Fac[a * nvir + c];
        for (let k = 0; k < nocc; k++)
          for (let d = 0; d < nvir; d++)
            val += w(k, nocc + a, nocc + d, nocc + c) * t1[k * nvir + d];
        Lac[a * nvir + c] = val;
      }

    // ---- W^{kl}_{ij} ----
    const Wklij = new Float64Array(nocc * nocc * nocc * nocc);
    for (let k = 0; k < nocc; k++)
      for (let l = 0; l < nocc; l++)
        for (let i = 0; i < nocc; i++)
          for (let j = 0; j < nocc; j++) {
            let val = v(k, l, i, j);
            for (let c = 0; c < nvir; c++) {
              const C = nocc + c;
              val += v(l, k, C, i) * t1[j * nvir + c];
              val += v(k, l, C, j) * t1[i * nvir + c];
            }
            for (let c = 0; c < nvir; c++)
              for (let d = 0; d < nvir; d++) {
                const C = nocc + c, D = nocc + d;
                val += v(k, l, C, D) * (t2[T2(i, j, c, d)] + t1[i * nvir + c] * t1[j * nvir + d]);
              }
            Wklij[((k * nocc + l) * nocc + i) * nocc + j] = val;
          }

    // ---- W^{ab}_{cd} ----
    const vv = nvir * nvir;
    const Wabcd = new Float64Array(vv * vv);
    for (let a = 0; a < nvir; a++)
      for (let b = 0; b < nvir; b++) {
        const A = nocc + a, B = nocc + b;
        for (let c = 0; c < nvir; c++)
          for (let d = 0; d < nvir; d++) {
            const C = nocc + c, D = nocc + d;
            let val = v(A, B, C, D);
            for (let k = 0; k < nocc; k++) {
              val -= v(k, A, D, C) * t1[k * nvir + b];
              val -= v(k, B, C, D) * t1[k * nvir + a];
            }
            Wabcd[(a * nvir + b) * vv + c * nvir + d] = val;
          }
      }

    // ---- W^{ak}_{ic} and W^{ak}_{ci} (exchange intermediates) ----
    const Wakic = new Float64Array(nvir * nocc * nocc * nvir);
    for (let a = 0; a < nvir; a++)
      for (let k = 0; k < nocc; k++)
        for (let i = 0; i < nocc; i++)
          for (let c = 0; c < nvir; c++) {
            const A = nocc + a, C = nocc + c;
            let val = v(A, k, i, C);
            for (let l = 0; l < nocc; l++)
              val -= v(k, l, C, i) * t1[l * nvir + a];
            for (let d = 0; d < nvir; d++) {
              const D = nocc + d;
              val += v(k, A, C, D) * t1[i * nvir + d];
            }
            for (let l = 0; l < nocc; l++)
              for (let d = 0; d < nvir; d++) {
                const D = nocc + d;
                const vlk = v(l, k, D, C);
                val -= 0.5 * vlk * t2[T2(i, l, d, a)];
                val -= vlk * t1[i * nvir + d] * t1[l * nvir + a];
                val += 0.5 * w(l, k, D, C) * t2[T2(i, l, a, d)];
              }
            Wakic[((a * nocc + k) * nocc + i) * nvir + c] = val;
          }

    const Wakci = new Float64Array(nvir * nocc * nvir * nocc);
    for (let a = 0; a < nvir; a++)
      for (let k = 0; k < nocc; k++)
        for (let c = 0; c < nvir; c++)
          for (let i = 0; i < nocc; i++) {
            const A = nocc + a, C = nocc + c;
            let val = v(A, k, C, i);
            for (let l = 0; l < nocc; l++)
              val -= v(l, k, C, i) * t1[l * nvir + a];
            for (let d = 0; d < nvir; d++) {
              const D = nocc + d;
              val += v(k, A, D, C) * t1[i * nvir + d];
            }
            for (let l = 0; l < nocc; l++)
              for (let d = 0; d < nvir; d++) {
                const D = nocc + d;
                const vlk = v(l, k, C, D);
                val -= 0.5 * vlk * t2[T2(i, l, d, a)];
                val -= vlk * t1[i * nvir + d] * t1[l * nvir + a];
              }
            Wakci[((a * nocc + k) * nvir + c) * nocc + i] = val;
          }

    // ---- T1 update ----
    const newT1 = new Float64Array(t1Size);
    for (let i = 0; i < nocc; i++)
      for (let a = 0; a < nvir; a++) {
        const A = nocc + a;
        let val = 0;
        // Fac * t1
        for (let c = 0; c < nvir; c++)
          val += Fac[a * nvir + c] * t1[i * nvir + c];
        // -Fki * t1
        for (let k = 0; k < nocc; k++)
          val -= Fki[k * nocc + i] * t1[k * nvir + a];
        // 2*Fkc*t2(ki,ca) - Fkc*t2(ik,ca) + Fkc*t1(i,c)*t1(k,a)
        for (let k = 0; k < nocc; k++)
          for (let c = 0; c < nvir; c++) {
            const fc = Fkc[k * nvir + c];
            val += fc * (2.0 * t2[T2(k, i, c, a)] - t2[T2(i, k, c, a)] + t1[i * nvir + c] * t1[k * nvir + a]);
          }
        // w(A,k,i,C) * t1(k,c)
        for (let k = 0; k < nocc; k++)
          for (let c = 0; c < nvir; c++)
            val += w(A, k, i, nocc + c) * t1[k * nvir + c];
        // w(A,k,C,D) * (t2(ik,cd) + t1(i,c)*t1(k,d))
        for (let k = 0; k < nocc; k++)
          for (let c = 0; c < nvir; c++)
            for (let d = 0; d < nvir; d++)
              val += w(A, k, nocc + c, nocc + d) * (t2[T2(i, k, c, d)] + t1[i * nvir + c] * t1[k * nvir + d]);
        // -w(k,l,i,C) * (t2(kl,ac) + t1(k,a)*t1(l,c))
        for (let k = 0; k < nocc; k++)
          for (let l = 0; l < nocc; l++)
            for (let c = 0; c < nvir; c++)
              val -= w(k, l, i, nocc + c) * (t2[T2(k, l, a, c)] + t1[k * nvir + a] * t1[l * nvir + c]);
        newT1[i * nvir + a] = val / Dia[i * nvir + a];
      }

    // ---- T2 update ----
    // tau[ij,ab] = t2[i,j,a,b] + t1[i,a]*t1[j,b]
    const tau = new Float64Array(t2Size);
    for (let i = 0; i < nocc; i++)
      for (let j = 0; j < nocc; j++)
        for (let a = 0; a < nvir; a++)
          for (let b = 0; b < nvir; b++) {
            const idx = T2(i, j, a, b);
            tau[idx] = t2[idx] + t1[i * nvir + a] * t1[j * nvir + b];
          }

    // raw(i,a,j,b): compute per (i,a,j,b), then symmetrize
    const raw = new Float64Array(t2Size);

    // Batch: 0.5 * tau[ij,cd] * Wabcd[ab,cd]^T  and  0.5 * Wklij^T[ij,kl] * tau[kl,ab]
    // Using layout: raw[T2(i,j,a,b)], tau[T2(i,j,c,d)], Wklij[(k*nocc+l)*nocc^2 + i*nocc+j]
    const oo = nocc * nocc;
    // 0.5 * Σ_{cd} tau[ij,cd] * Wabcd[ab,cd]
    for (let ij = 0; ij < oo; ij++) {
      const ijOff = ij * vv;
      for (let ab = 0; ab < vv; ab++) {
        let sum = 0;
        const abOff = ab * vv;
        for (let cd = 0; cd < vv; cd++)
          sum += tau[ijOff + cd] * Wabcd[abOff + cd];
        raw[ijOff + ab] += 0.5 * sum;
      }
    }
    // 0.5 * Σ_{kl} Wklij[kl,ij] * tau[kl,ab]
    for (let ij = 0; ij < oo; ij++) {
      for (let ab = 0; ab < vv; ab++) {
        let sum = 0;
        for (let kl = 0; kl < oo; kl++)
          sum += Wklij[kl * oo + ij] * tau[kl * vv + ab];
        raw[ij * vv + ab] += 0.5 * sum;
      }
    }

    // Add remaining terms to raw
    for (let i = 0; i < nocc; i++)
      for (let a = 0; a < nvir; a++) {
        const A = nocc + a;
        for (let j = 0; j < nocc; j++)
          for (let b = 0; b < nvir; b++) {
            const B = nocc + b;
            let val = 0.5 * v(i, j, A, B);

            // Lac * t2(i,j,c,b)
            for (let c = 0; c < nvir; c++)
              val += Lac[a * nvir + c] * t2[T2(i, j, c, b)];

            // -Lki * t2(k,j,a,b)
            for (let k = 0; k < nocc; k++)
              val -= Lki[k * nocc + i] * t2[T2(k, j, a, b)];

            // v(A,B,i,C) * t1(j,c)
            for (let c = 0; c < nvir; c++)
              val += v(A, B, i, nocc + c) * t1[j * nvir + c];

            // -v(k,B,i,C) * t1(k,a) * t1(j,c)
            for (let k = 0; k < nocc; k++)
              for (let c = 0; c < nvir; c++)
                val -= v(k, B, i, nocc + c) * t1[k * nvir + a] * t1[j * nvir + c];

            // -v(A,k,i,j) * t1(k,b)
            for (let k = 0; k < nocc; k++)
              val -= v(A, k, i, j) * t1[k * nvir + b];

            // -v(A,k,i,C) * t1(j,c) * t1(k,b)
            for (let k = 0; k < nocc; k++)
              for (let c = 0; c < nvir; c++)
                val -= v(A, k, i, nocc + c) * t1[j * nvir + c] * t1[k * nvir + b];

            // W exchange terms:
            // 2*Wakic*t2(kj,cb) - Wakci*t2(kj,cb) - Wakic*t2(kj,bc) - Wbkci*t2(kj,ac)
            for (let k = 0; k < nocc; k++)
              for (let c = 0; c < nvir; c++) {
                const w1 = Wakic[((a * nocc + k) * nocc + i) * nvir + c];
                const w2 = Wakci[((a * nocc + k) * nvir + c) * nocc + i];
                const w3 = Wakci[((b * nocc + k) * nvir + c) * nocc + i];
                val += 2.0 * w1 * t2[T2(k, j, c, b)];
                val -= w2 * t2[T2(k, j, c, b)];
                val -= w1 * t2[T2(k, j, b, c)];
                val -= w3 * t2[T2(k, j, a, c)];
              }

            raw[T2(i, j, a, b)] += val;
          }
      }

    // Symmetrize: t2_new(i,j,a,b) = [raw(i,j,a,b) + raw(j,i,b,a)] / D
    const newT2 = new Float64Array(t2Size);
    for (let i = 0; i < nocc; i++)
      for (let j = 0; j < nocc; j++)
        for (let a = 0; a < nvir; a++)
          for (let b = 0; b < nvir; b++) {
            const idx = T2(i, j, a, b);
            newT2[idx] = (raw[idx] + raw[T2(j, i, b, a)]) / Dijab[idx];
          }

    // DIIS extrapolation
    const ampVec = new Float64Array(t1Size + t2Size);
    const errVec = new Float64Array(t1Size + t2Size);
    for (let k = 0; k < t1Size; k++) {
      ampVec[k] = newT1[k];
      errVec[k] = newT1[k] - t1[k];
    }
    for (let k = 0; k < t2Size; k++) {
      ampVec[t1Size + k] = newT2[k];
      errVec[t1Size + k] = newT2[k] - t2[k];
    }
    diis.push(ampVec, errVec);
    if (diis.canExtrapolate()) {
      const extrap = diis.extrapolate();
      for (let k = 0; k < t1Size; k++) newT1[k] = extrap[k];
      for (let k = 0; k < t2Size; k++) newT2[k] = extrap[t1Size + k] as number;
    }

    t1 = newT1;
    t2 = newT2;

    const newEcc = energy();
    const deltaE = newEcc - Ecc;
    Ecc = newEcc;

    onProgress?.(`CCSD iter ${iter.toString().padStart(2)}: E = ${Ecc.toFixed(12)}, ΔE = ${deltaE.toExponential(4)}`);

    if (Math.abs(deltaE) < CONV) {
      onProgress?.(`CCSD converged after ${iter} iterations`);
      return Ecc;
    }
  }

  onProgress?.(`CCSD: WARNING — not converged in ${MAX_ITER} iterations`);
  return Ecc;
}

// ---------------------------------------------------------------------------
// Spin-orbital CCSD loop (for UHF/ROHF)
// ---------------------------------------------------------------------------

function ccsdLoop(
  getInt: (P: number, Q: number, R: number, S: number) => number,
  fock: Float64Array,
  nocc: number, nvir: number, nso: number,
  onProgress?: (msg: string) => void,
): number {
  // Denominators
  const Dia = new Float64Array(nocc * nvir);
  for (let i = 0; i < nocc; i++) {
    for (let a = 0; a < nvir; a++) {
      Dia[i * nvir + a] = fock[i * nso + i] - fock[(nocc + a) * nso + (nocc + a)];
    }
  }
  const t2Size = nocc * nocc * nvir * nvir;
  const Dijab = new Float64Array(t2Size);
  for (let i = 0; i < nocc; i++) {
    for (let j = 0; j < nocc; j++) {
      for (let a = 0; a < nvir; a++) {
        for (let b = 0; b < nvir; b++) {
          Dijab[t2idx(i, j, a, b, nocc, nvir)] =
            fock[i * nso + i] + fock[j * nso + j]
            - fock[(nocc + a) * nso + (nocc + a)]
            - fock[(nocc + b) * nso + (nocc + b)];
        }
      }
    }
  }

  // MP2 initial guess: T1=0, T2 = <ij||ab>/D
  const t1Size = nocc * nvir;
  let t1 = new Float64Array(t1Size);
  let t2 = new Float64Array(t2Size);
  for (let i = 0; i < nocc; i++) {
    for (let j = 0; j < nocc; j++) {
      for (let a = 0; a < nvir; a++) {
        for (let b = 0; b < nvir; b++) {
          const idx = t2idx(i, j, a, b, nocc, nvir);
          t2[idx] = getInt(i, j, nocc + a, nocc + b) / Dijab[idx];
        }
      }
    }
  }

  // Energy function
  function energy(): number {
    let E = 0;
    for (let i = 0; i < nocc; i++) {
      for (let a = 0; a < nvir; a++) {
        const A = nocc + a;
        E += fock[i * nso + A] * t1[i * nvir + a];
        for (let j = 0; j < nocc; j++) {
          for (let b = 0; b < nvir; b++) {
            const ijab = getInt(i, j, A, nocc + b);
            E += 0.25 * ijab * t2[t2idx(i, j, a, b, nocc, nvir)];
            E += 0.5 * ijab * t1[i * nvir + a] * t1[j * nvir + b];
          }
        }
      }
    }
    return E;
  }

  let Ecc = energy();
  onProgress?.(`CCSD iter  0: E = ${Ecc.toFixed(12)} (MP2 initial guess)`);

  const diis = new DIIS(8, 4);
  const MAX_ITER = 100;
  const CONV = 1e-10;

  for (let iter = 1; iter <= MAX_ITER; iter++) {
    // τ̃ and τ
    const tauTilde = new Float64Array(t2Size);
    const tau = new Float64Array(t2Size);
    for (let i = 0; i < nocc; i++) {
      for (let j = 0; j < nocc; j++) {
        for (let a = 0; a < nvir; a++) {
          for (let b = 0; b < nvir; b++) {
            const idx = t2idx(i, j, a, b, nocc, nvir);
            const ia_jb = t1[i * nvir + a] * t1[j * nvir + b];
            const ib_ja = t1[i * nvir + b] * t1[j * nvir + a];
            tauTilde[idx] = t2[idx] + 0.5 * (ia_jb - ib_ja);
            tau[idx] = t2[idx] + ia_jb - ib_ja;
          }
        }
      }
    }

    // --- Stanton Eq.3: F_ae ---
    const Fae = new Float64Array(nvir * nvir);
    for (let a = 0; a < nvir; a++) {
      for (let e = 0; e < nvir; e++) {
        const A = nocc + a, E = nocc + e;
        let val = (a !== e) ? fock[A * nso + E] : 0;
        for (let m = 0; m < nocc; m++) val -= 0.5 * fock[m * nso + E] * t1[m * nvir + a];
        for (let m = 0; m < nocc; m++)
          for (let f = 0; f < nvir; f++)
            val += t1[m * nvir + f] * getInt(m, A, nocc + f, E);
        for (let m = 0; m < nocc; m++)
          for (let n = 0; n < nocc; n++)
            for (let f = 0; f < nvir; f++)
              val -= 0.5 * tauTilde[t2idx(m, n, a, f, nocc, nvir)] * getInt(m, n, E, nocc + f);
        Fae[a * nvir + e] = val;
      }
    }

    // --- Stanton Eq.4: F_mi ---
    const Fmi = new Float64Array(nocc * nocc);
    for (let m = 0; m < nocc; m++) {
      for (let i = 0; i < nocc; i++) {
        let val = (m !== i) ? fock[m * nso + i] : 0;
        for (let e = 0; e < nvir; e++) val += 0.5 * fock[m * nso + (nocc + e)] * t1[i * nvir + e];
        for (let n = 0; n < nocc; n++)
          for (let e = 0; e < nvir; e++)
            val += t1[n * nvir + e] * getInt(m, n, i, nocc + e);
        for (let n = 0; n < nocc; n++)
          for (let e = 0; e < nvir; e++)
            for (let f = 0; f < nvir; f++)
              val += 0.5 * tauTilde[t2idx(i, n, e, f, nocc, nvir)] * getInt(m, n, nocc + e, nocc + f);
        Fmi[m * nocc + i] = val;
      }
    }

    // --- Stanton Eq.5: F_me ---
    const Fme = new Float64Array(nocc * nvir);
    for (let m = 0; m < nocc; m++) {
      for (let e = 0; e < nvir; e++) {
        let val = fock[m * nso + (nocc + e)];
        for (let n = 0; n < nocc; n++)
          for (let f = 0; f < nvir; f++)
            val += t1[n * nvir + f] * getInt(m, n, nocc + e, nocc + f);
        Fme[m * nvir + e] = val;
      }
    }

    // --- Stanton Eq.6: W_mnij ---
    const Wmnij = new Float64Array(nocc * nocc * nocc * nocc);
    for (let m = 0; m < nocc; m++) {
      for (let n = 0; n < nocc; n++) {
        for (let i = 0; i < nocc; i++) {
          for (let j = 0; j < nocc; j++) {
            let val = getInt(m, n, i, j);
            for (let e = 0; e < nvir; e++) {
              val += t1[j * nvir + e] * getInt(m, n, i, nocc + e);
              val -= t1[i * nvir + e] * getInt(m, n, j, nocc + e);
            }
            for (let e = 0; e < nvir; e++)
              for (let f = 0; f < nvir; f++)
                val += 0.25 * tau[t2idx(i, j, e, f, nocc, nvir)] * getInt(m, n, nocc + e, nocc + f);
            Wmnij[((m * nocc + n) * nocc + i) * nocc + j] = val;
          }
        }
      }
    }

    // --- Stanton Eq.7: W_abef ---
    const Wabef = new Float64Array(nvir * nvir * nvir * nvir);
    for (let a = 0; a < nvir; a++) {
      for (let b = 0; b < nvir; b++) {
        for (let e = 0; e < nvir; e++) {
          for (let f = 0; f < nvir; f++) {
            const A = nocc + a, B = nocc + b, E = nocc + e, F = nocc + f;
            let val = getInt(A, B, E, F);
            for (let m = 0; m < nocc; m++) {
              val -= t1[m * nvir + b] * getInt(A, m, E, F);
              val += t1[m * nvir + a] * getInt(B, m, E, F);
            }
            for (let m = 0; m < nocc; m++)
              for (let n = 0; n < nocc; n++)
                val += 0.25 * tau[t2idx(m, n, a, b, nocc, nvir)] * getInt(m, n, E, F);
            Wabef[((a * nvir + b) * nvir + e) * nvir + f] = val;
          }
        }
      }
    }

    // --- Stanton Eq.8: W_mbej ---
    const Wmbej = new Float64Array(nocc * nvir * nvir * nocc);
    for (let m = 0; m < nocc; m++) {
      for (let b = 0; b < nvir; b++) {
        for (let e = 0; e < nvir; e++) {
          for (let j = 0; j < nocc; j++) {
            const B = nocc + b, E = nocc + e;
            let val = getInt(m, B, E, j);
            for (let f = 0; f < nvir; f++)
              val += t1[j * nvir + f] * getInt(m, B, E, nocc + f);
            for (let n = 0; n < nocc; n++)
              val -= t1[n * nvir + b] * getInt(m, n, E, j);
            for (let n = 0; n < nocc; n++)
              for (let f = 0; f < nvir; f++)
                val -= (0.5 * t2[t2idx(j, n, f, b, nocc, nvir)] + t1[j * nvir + f] * t1[n * nvir + b])
                     * getInt(m, n, E, nocc + f);
            Wmbej[((m * nvir + b) * nvir + e) * nocc + j] = val;
          }
        }
      }
    }

    // --- Stanton Eq.1: T1 update ---
    const newT1 = new Float64Array(t1Size);
    for (let i = 0; i < nocc; i++) {
      for (let a = 0; a < nvir; a++) {
        const A = nocc + a;
        let val = fock[i * nso + A];
        for (let e = 0; e < nvir; e++) val += t1[i * nvir + e] * Fae[a * nvir + e];
        for (let m = 0; m < nocc; m++) val -= t1[m * nvir + a] * Fmi[m * nocc + i];
        for (let m = 0; m < nocc; m++)
          for (let e = 0; e < nvir; e++)
            val += t2[t2idx(i, m, a, e, nocc, nvir)] * Fme[m * nvir + e];
        for (let n = 0; n < nocc; n++)
          for (let f = 0; f < nvir; f++)
            val -= t1[n * nvir + f] * getInt(n, A, i, nocc + f);
        for (let m = 0; m < nocc; m++)
          for (let e = 0; e < nvir; e++)
            for (let f = 0; f < nvir; f++)
              val -= 0.5 * t2[t2idx(i, m, e, f, nocc, nvir)] * getInt(m, A, nocc + e, nocc + f);
        for (let m = 0; m < nocc; m++)
          for (let n = 0; n < nocc; n++)
            for (let e = 0; e < nvir; e++)
              val -= 0.5 * t2[t2idx(m, n, a, e, nocc, nvir)] * getInt(n, m, nocc + e, i);
        newT1[i * nvir + a] = val / Dia[i * nvir + a];
      }
    }

    // --- Stanton Eq.2: T2 update ---
    const newT2 = new Float64Array(t2Size);
    for (let i = 0; i < nocc; i++) {
      for (let j = 0; j < nocc; j++) {
        for (let a = 0; a < nvir; a++) {
          for (let b = 0; b < nvir; b++) {
            const A = nocc + a, B = nocc + b;
            const idx = t2idx(i, j, a, b, nocc, nvir);
            let val = getInt(i, j, A, B);

            for (let e = 0; e < nvir; e++) {
              let Fb = Fae[b * nvir + e];
              for (let m = 0; m < nocc; m++) Fb -= 0.5 * t1[m * nvir + b] * Fme[m * nvir + e];
              val += t2[t2idx(i, j, a, e, nocc, nvir)] * Fb;
              let Fa = Fae[a * nvir + e];
              for (let m = 0; m < nocc; m++) Fa -= 0.5 * t1[m * nvir + a] * Fme[m * nvir + e];
              val -= t2[t2idx(i, j, b, e, nocc, nvir)] * Fa;
            }

            for (let m = 0; m < nocc; m++) {
              let Fj = Fmi[m * nocc + j];
              for (let e = 0; e < nvir; e++) Fj += 0.5 * t1[j * nvir + e] * Fme[m * nvir + e];
              val -= t2[t2idx(i, m, a, b, nocc, nvir)] * Fj;
              let Fi = Fmi[m * nocc + i];
              for (let e = 0; e < nvir; e++) Fi += 0.5 * t1[i * nvir + e] * Fme[m * nvir + e];
              val += t2[t2idx(j, m, a, b, nocc, nvir)] * Fi;
            }

            for (let m = 0; m < nocc; m++)
              for (let n = 0; n < nocc; n++)
                val += 0.5 * tau[t2idx(m, n, a, b, nocc, nvir)]
                     * Wmnij[((m * nocc + n) * nocc + i) * nocc + j];

            for (let e = 0; e < nvir; e++)
              for (let f = 0; f < nvir; f++)
                val += 0.5 * tau[t2idx(i, j, e, f, nocc, nvir)]
                     * Wabef[((a * nvir + b) * nvir + e) * nvir + f];

            for (let m = 0; m < nocc; m++) {
              for (let e = 0; e < nvir; e++) {
                const E = nocc + e;
                val += t2[t2idx(i, m, a, e, nocc, nvir)] * Wmbej[((m * nvir + b) * nvir + e) * nocc + j];
                val -= t1[i * nvir + e] * t1[m * nvir + a] * getInt(m, B, E, j);
                val -= t2[t2idx(j, m, a, e, nocc, nvir)] * Wmbej[((m * nvir + b) * nvir + e) * nocc + i];
                val += t1[j * nvir + e] * t1[m * nvir + a] * getInt(m, B, E, i);
                val -= t2[t2idx(i, m, b, e, nocc, nvir)] * Wmbej[((m * nvir + a) * nvir + e) * nocc + j];
                val += t1[i * nvir + e] * t1[m * nvir + b] * getInt(m, A, E, j);
                val += t2[t2idx(j, m, b, e, nocc, nvir)] * Wmbej[((m * nvir + a) * nvir + e) * nocc + i];
                val -= t1[j * nvir + e] * t1[m * nvir + b] * getInt(m, A, E, i);
              }
            }

            for (let e = 0; e < nvir; e++) {
              val += t1[i * nvir + e] * getInt(A, B, nocc + e, j);
              val -= t1[j * nvir + e] * getInt(A, B, nocc + e, i);
            }

            for (let m = 0; m < nocc; m++) {
              val -= t1[m * nvir + a] * getInt(m, B, i, j);
              val += t1[m * nvir + b] * getInt(m, A, i, j);
            }

            newT2[idx] = val / Dijab[idx];
          }
        }
      }
    }

    // DIIS extrapolation
    const ampVec = new Float64Array(t1Size + t2Size);
    const errVec = new Float64Array(t1Size + t2Size);
    for (let k = 0; k < t1Size; k++) {
      ampVec[k] = newT1[k];
      errVec[k] = newT1[k] - t1[k];
    }
    for (let k = 0; k < t2Size; k++) {
      ampVec[t1Size + k] = newT2[k];
      errVec[t1Size + k] = newT2[k] - t2[k];
    }
    diis.push(ampVec, errVec);
    if (diis.canExtrapolate()) {
      const extrap = diis.extrapolate();
      for (let k = 0; k < t1Size; k++) newT1[k] = extrap[k];
      for (let k = 0; k < t2Size; k++) newT2[k] = extrap[t1Size + k] as number;
    }

    t1 = newT1;
    t2 = newT2;

    const newEcc = energy();
    const deltaE = newEcc - Ecc;
    Ecc = newEcc;

    onProgress?.(`CCSD iter ${iter.toString().padStart(2)}: E = ${Ecc.toFixed(12)}, ΔE = ${deltaE.toExponential(4)}`);

    if (Math.abs(deltaE) < CONV) {
      onProgress?.(`CCSD converged after ${iter} iterations`);
      return Ecc;
    }
  }

  onProgress?.(`CCSD: WARNING — not converged in ${MAX_ITER} iterations`);
  return Ecc;
}

// ---------------------------------------------------------------------------
// RHF-CCSD entry point (spatial orbital)
// ---------------------------------------------------------------------------

/** Compute CCSD correlation energy from converged RHF results.
 *  Uses spatial-orbital formulation (~16x faster than spin-orbital). */
export function computeCCSDEnergy(
  C: Matrix, epsilon: FloatArray, eri: ERIStored,
  nocc: number, nbasis: number,
  onProgress?: (msg: string) => void,
): number {
  onProgress?.('CCSD: Full MO integral transformation...');
  const moInts = fullMOTransform(C, C, C, C, eri, nbasis);

  return ccsdSpatialLoop(moInts, epsilon, nocc, nbasis, onProgress);
}

// ---------------------------------------------------------------------------
// UHF-CCSD entry point (spin-orbital, blocked ordering)
// ---------------------------------------------------------------------------

/** Compute CCSD correlation energy from converged UHF results.
 *  Uses blocked spin-orbital ordering:
 *    [0, noccA) = occupied α,  [noccA, noccA+noccB) = occupied β,
 *    [nocc_s, nocc_s+nvirA) = virtual α,  [nocc_s+nvirA, nso) = virtual β
 */
export function computeUCCSDEnergy(
  Ca: Matrix, Cb: Matrix,
  epsilonA: FloatArray, epsilonB: FloatArray,
  eri: ERIStored,
  noccA: number, noccB: number, nbasis: number,
  onProgress?: (msg: string) => void,

): number {
  const N = nbasis;
  const nvirA = N - noccA, nvirB = N - noccB;
  const nocc_s = noccA + noccB;
  const nvir_s = nvirA + nvirB;
  const nso = nocc_s + nvir_s;

  onProgress?.('UCCSD: MO integral transformation (3 sets)...');
  const moAA = fullMOTransform(Ca, Ca, Ca, Ca, eri, N);
  const moBB = fullMOTransform(Cb, Cb, Cb, Cb, eri, N);
  const moAB = fullMOTransform(Ca, Ca, Cb, Cb, eri, N);

  const N3 = N * N * N, N2 = N * N;

  // Spin type and spatial index from blocked spin-orbital index
  function spinOf(P: number): number {
    if (P < noccA) return 0;
    if (P < nocc_s) return 1;
    if (P < nocc_s + nvirA) return 0;
    return 1;
  }
  function spatOf(P: number): number {
    if (P < noccA) return P;
    if (P < nocc_s) return P - noccA;
    if (P < nocc_s + nvirA) return noccA + (P - nocc_s);
    return noccB + (P - nocc_s - nvirA);
  }

  function getInt(P: number, Q: number, R: number, S: number): number {
    const sP = spinOf(P), sQ = spinOf(Q), sR = spinOf(R), sS = spinOf(S);
    const p = spatOf(P), q = spatOf(Q), r = spatOf(R), s = spatOf(S);
    let val = 0;
    // Direct: (pr|qs) δ(σP,σR) δ(σQ,σS)
    if (sP === sR && sQ === sS) {
      if (sP === 0 && sQ === 0)      val += moAA[p * N3 + r * N2 + q * N + s];
      else if (sP === 1 && sQ === 1) val += moBB[p * N3 + r * N2 + q * N + s];
      else if (sP === 0 && sQ === 1) val += moAB[p * N3 + r * N2 + q * N + s];
      else                           val += moAB[q * N3 + s * N2 + p * N + r]; // (βα|βα)→transpose
    }
    // Exchange: -(ps|qr) δ(σP,σS) δ(σQ,σR)
    if (sP === sS && sQ === sR) {
      if (sP === 0 && sQ === 0)      val -= moAA[p * N3 + s * N2 + q * N + r];
      else if (sP === 1 && sQ === 1) val -= moBB[p * N3 + s * N2 + q * N + r];
      else if (sP === 0 && sQ === 1) val -= moAB[p * N3 + s * N2 + q * N + r];
      else                           val -= moAB[q * N3 + r * N2 + p * N + s];
    }
    return val;
  }

  // Blocked spin-orbital Fock
  const fock = new Float64Array(nso * nso);
  for (let i = 0; i < noccA; i++) fock[i * nso + i] = epsilonA[i];
  for (let i = 0; i < noccB; i++) fock[(noccA + i) * nso + (noccA + i)] = epsilonB[i];
  for (let a = 0; a < nvirA; a++) fock[(nocc_s + a) * nso + (nocc_s + a)] = epsilonA[noccA + a];
  for (let a = 0; a < nvirB; a++) fock[(nocc_s + nvirA + a) * nso + (nocc_s + nvirA + a)] = epsilonB[noccB + a];

  return ccsdLoop(getInt, fock, nocc_s, nvir_s, nso, onProgress);
}

// ---------------------------------------------------------------------------
// ROHF-CCSD entry point (semicanonical orbitals)
// ---------------------------------------------------------------------------

/** Compute CCSD correlation energy from converged ROHF results.
 *  Builds semicanonical orbitals by diagonalizing Fα and Fβ within
 *  occ-occ and vir-vir blocks of the ROHF MO basis, then delegates
 *  to the UCCSD solver. */
export function computeROCCSDEnergy(
  C: Matrix, Fa: Matrix, Fb: Matrix,
  eri: ERIStored,
  noccA: number, noccB: number, nbasis: number,
  onProgress?: (msg: string) => void,

): number {
  const N = nbasis;

  onProgress?.('ROCCSD: Building semicanonical orbitals...');

  // Transform Fα, Fβ to ROHF MO basis
  const FaMO = matmul(matmulAtB(C, Fa), C);
  const FbMO = matmul(matmulAtB(C, Fb), C);

  // --- Semicanonical α ---
  const occA = new Matrix(noccA, noccA);
  for (let i = 0; i < noccA; i++)
    for (let j = 0; j < noccA; j++)
      occA.set(i, j, FaMO.get(i, j));
  const eigOccA = jacobiEigen(occA);

  const nvirA = N - noccA;
  const virA = new Matrix(nvirA, nvirA);
  for (let i = 0; i < nvirA; i++)
    for (let j = 0; j < nvirA; j++)
      virA.set(i, j, FaMO.get(noccA + i, noccA + j));
  const eigVirA = jacobiEigen(virA);

  // Block-diagonal Uα
  const Ua = new Matrix(N, N);
  for (let i = 0; i < noccA; i++)
    for (let j = 0; j < noccA; j++)
      Ua.set(i, j, eigOccA.eigenvectors.get(i, j));
  for (let i = 0; i < nvirA; i++)
    for (let j = 0; j < nvirA; j++)
      Ua.set(noccA + i, noccA + j, eigVirA.eigenvectors.get(i, j));

  const CaSc = matmul(C, Ua);
  const epsilonA = new Float64Array(N);
  for (let i = 0; i < noccA; i++) epsilonA[i] = eigOccA.eigenvalues[i];
  for (let i = 0; i < nvirA; i++) epsilonA[noccA + i] = eigVirA.eigenvalues[i];

  // --- Semicanonical β ---
  const occB = new Matrix(noccB, noccB);
  for (let i = 0; i < noccB; i++)
    for (let j = 0; j < noccB; j++)
      occB.set(i, j, FbMO.get(i, j));
  const eigOccB = jacobiEigen(occB);

  const nvirB = N - noccB;
  const virB = new Matrix(nvirB, nvirB);
  for (let i = 0; i < nvirB; i++)
    for (let j = 0; j < nvirB; j++)
      virB.set(i, j, FbMO.get(noccB + i, noccB + j));
  const eigVirB = jacobiEigen(virB);

  // Block-diagonal Uβ
  const Ub = new Matrix(N, N);
  for (let i = 0; i < noccB; i++)
    for (let j = 0; j < noccB; j++)
      Ub.set(i, j, eigOccB.eigenvectors.get(i, j));
  for (let i = 0; i < nvirB; i++)
    for (let j = 0; j < nvirB; j++)
      Ub.set(noccB + i, noccB + j, eigVirB.eigenvectors.get(i, j));

  const CbSc = matmul(C, Ub);
  const epsilonB = new Float64Array(N);
  for (let i = 0; i < noccB; i++) epsilonB[i] = eigOccB.eigenvalues[i];
  for (let i = 0; i < nvirB; i++) epsilonB[noccB + i] = eigVirB.eigenvalues[i];

  // Delegate to UCCSD with semicanonical orbitals
  return computeUCCSDEnergy(CaSc, CbSc, epsilonA, epsilonB, eri, noccA, noccB, nbasis, onProgress);
}
