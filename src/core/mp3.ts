/** MP3 (Third-order Møller-Plesset perturbation theory) correlation energy.
 *
 *  RHF uses spatial-orbital (chemist notation) formulas directly:
 *    PP (2h4p): Σ (ia|jb)(ac|bd)(2(ic|jd)-(id|jc)) / (D_ij^ab × D_ij^cd)
 *    HH (4h2p): Σ (ia|jb)(ik|jl)(2(ka|lb)-(kb|la)) / (D_ij^ab × D_kl^ab)
 *    PH (3h3p): Σ ((2(ia|jb)-(ij|ab))(2(kc|ia)-(ka|ic))(2(kc|jb)-(kb|jc))
 *                   - 3(ij|ab)(ka|ic)(kb|jc)) / (D_ik^ac × D_kj^bc)
 *
 *  UHF/ROHF uses spatial-orbital formulas with explicit spin channels (αα, ββ, αβ).
 *  Reuses full 4-index MO transform from ccsd.ts.
 */

import type { ERIStored } from './eri';
import { Matrix, type FloatArray } from '../linalg/matrix';
import { fullMOTransform } from './ccsd';
import { jacobiEigen } from '../linalg/eigendecomposition';
import { matmul, matmulAtB } from '../linalg/matmul';

// ---------------------------------------------------------------------------
// RHF-MP3 — spatial orbital (chemist notation) formulas
// ---------------------------------------------------------------------------

/** Compute MP3 correlation energy from converged RHF results.
 *  Uses spatial-orbital formulas from GANSU C++ (proven correct).
 *  Returns both MP2 and MP3 correlation energies. */
export function computeMP3Energy(
  C: Matrix, epsilon: FloatArray, eri: ERIStored,
  nocc: number, nbasis: number,
  onProgress?: (msg: string) => void,
): { mp2: number; mp3: number } {
  const N = nbasis;

  onProgress?.('MP3: Full MO integral transformation...');
  const moInts = fullMOTransform(C, C, C, C, eri, N);

  // Chemist-notation MO integral accessor: (pq|rs)
  const N3 = N * N * N, N2 = N * N;
  function mo(p: number, q: number, r: number, s: number): number {
    return moInts[p * N3 + q * N2 + r * N + s];
  }

  // --- MP2 energy (spatial orbital, spin-traced) ---
  let Emp2 = 0;
  for (let i = 0; i < nocc; i++) {
    for (let j = 0; j < nocc; j++) {
      for (let a = nocc; a < N; a++) {
        for (let b = nocc; b < N; b++) {
          const iajb = mo(i, a, j, b);
          const ibja = mo(i, b, j, a);
          const D = epsilon[i] + epsilon[j] - epsilon[a] - epsilon[b];
          Emp2 += iajb * (2 * iajb - ibja) / D;
        }
      }
    }
  }
  onProgress?.(`MP3: MP2 correlation energy = ${Emp2.toFixed(10)}`);
  onProgress?.('MP3: Computing third-order correction...');

  // --- PP (2h4p): particle-particle ladder ---
  let Epp = 0;
  for (let i = 0; i < nocc; i++) {
    for (let j = 0; j < nocc; j++) {
      for (let a = nocc; a < N; a++) {
        for (let b = nocc; b < N; b++) {
          const D1 = epsilon[i] + epsilon[j] - epsilon[a] - epsilon[b];
          const iajb = mo(i, a, j, b);
          for (let c = nocc; c < N; c++) {
            for (let d = nocc; d < N; d++) {
              const D2 = epsilon[i] + epsilon[j] - epsilon[c] - epsilon[d];
              Epp += iajb * mo(a, c, b, d) * (2 * mo(i, c, j, d) - mo(i, d, j, c)) / (D1 * D2);
            }
          }
        }
      }
    }
  }

  // --- HH (4h2p): hole-hole ladder ---
  let Ehh = 0;
  for (let i = 0; i < nocc; i++) {
    for (let j = 0; j < nocc; j++) {
      for (let k = 0; k < nocc; k++) {
        for (let l = 0; l < nocc; l++) {
          const ikjl = mo(i, k, j, l);
          for (let a = nocc; a < N; a++) {
            for (let b = nocc; b < N; b++) {
              const D1 = epsilon[i] + epsilon[j] - epsilon[a] - epsilon[b];
              const D2 = epsilon[k] + epsilon[l] - epsilon[a] - epsilon[b];
              Ehh += mo(i, a, j, b) * ikjl * (2 * mo(k, a, l, b) - mo(k, b, l, a)) / (D1 * D2);
            }
          }
        }
      }
    }
  }

  // --- PH (3h3p): particle-hole ring ---
  let Eph = 0;
  for (let i = 0; i < nocc; i++) {
    for (let j = 0; j < nocc; j++) {
      for (let k = 0; k < nocc; k++) {
        for (let a = nocc; a < N; a++) {
          for (let b = nocc; b < N; b++) {
            for (let c = nocc; c < N; c++) {
              const D1 = epsilon[i] + epsilon[k] - epsilon[a] - epsilon[c];
              const D2 = epsilon[k] + epsilon[j] - epsilon[b] - epsilon[c];
              const denom = D1 * D2;
              if (Math.abs(denom) < 1e-14) continue;
              const iajb = mo(i, a, j, b);
              const ijab = mo(i, j, a, b);
              const kcia = mo(k, c, i, a);
              const kaic = mo(k, a, i, c);
              const kcjb = mo(k, c, j, b);
              const kbjc = mo(k, b, j, c);
              Eph += ((2 * iajb - ijab) * (2 * kcia - kaic) * (2 * kcjb - kbjc)
                      - 3 * ijab * kaic * kbjc) / denom;
            }
          }
        }
      }
    }
  }

  const Emp3 = Epp + Ehh + Eph;
  onProgress?.(`MP3: pp=${Epp.toFixed(10)}, hh=${Ehh.toFixed(10)}, ph=${Eph.toFixed(10)}`);
  onProgress?.(`MP3 correction energy: ${Emp3.toFixed(10)} Hartree`);

  return { mp2: Emp2, mp3: Emp3 };
}

// ---------------------------------------------------------------------------
// Same-spin helpers for UHF-MP3 spatial orbital
// ---------------------------------------------------------------------------

/** Same-spin PP ladder: 0.5 Σ (acbd)(cidj)[(iajb)-(ibja)] / (D_ij^ab D_ij^cd) */
function ppSameSpin(mo: Float64Array, eps: FloatArray, nocc: number, N: number): number {
  const N3 = N * N * N, N2 = N * N;
  let E = 0;
  for (let i = 0; i < nocc; i++)
    for (let j = 0; j < nocc; j++)
      for (let a = nocc; a < N; a++)
        for (let b = nocc; b < N; b++) {
          const D1 = eps[i] + eps[j] - eps[a] - eps[b];
          const X = mo[i * N3 + a * N2 + j * N + b] - mo[i * N3 + b * N2 + j * N + a];
          if (Math.abs(X) < 1e-15) continue;
          for (let c = nocc; c < N; c++)
            for (let d = nocc; d < N; d++) {
              const D2 = eps[i] + eps[j] - eps[c] - eps[d];
              E += X * mo[a * N3 + c * N2 + b * N + d]
                     * mo[c * N3 + i * N2 + d * N + j] / (D1 * D2);
            }
        }
  return 0.5 * E;
}

/** Same-spin HH ladder: 0.5 Σ (kilj)(akbl)[(iajb)-(ibja)] / (D_ij^ab D_kl^ab) */
function hhSameSpin(mo: Float64Array, eps: FloatArray, nocc: number, N: number): number {
  const N3 = N * N * N, N2 = N * N;
  let E = 0;
  for (let i = 0; i < nocc; i++)
    for (let j = 0; j < nocc; j++)
      for (let a = nocc; a < N; a++)
        for (let b = nocc; b < N; b++) {
          const D1 = eps[i] + eps[j] - eps[a] - eps[b];
          const X = mo[i * N3 + a * N2 + j * N + b] - mo[i * N3 + b * N2 + j * N + a];
          if (Math.abs(X) < 1e-15) continue;
          for (let k = 0; k < nocc; k++)
            for (let l = 0; l < nocc; l++) {
              const D2 = eps[k] + eps[l] - eps[a] - eps[b];
              E += X * mo[k * N3 + i * N2 + l * N + j]
                     * mo[a * N3 + k * N2 + b * N + l] / (D1 * D2);
            }
        }
  return 0.5 * E;
}

/** Same-spin PH ring: Σ [(iajb)-(ibja)][(aikc)-(acki)][(cjbk)-(ckbj)] / (D_ij^ab D_kj^cb) */
function phSameSpin(mo: Float64Array, eps: FloatArray, nocc: number, N: number): number {
  const N3 = N * N * N, N2 = N * N;
  let E = 0;
  for (let i = 0; i < nocc; i++)
    for (let j = 0; j < nocc; j++)
      for (let a = nocc; a < N; a++)
        for (let b = nocc; b < N; b++) {
          const D1 = eps[i] + eps[j] - eps[a] - eps[b];
          const X1 = mo[i * N3 + a * N2 + j * N + b] - mo[i * N3 + b * N2 + j * N + a];
          if (Math.abs(X1) < 1e-15) continue;
          for (let k = 0; k < nocc; k++)
            for (let c = nocc; c < N; c++) {
              const D2 = eps[k] + eps[j] - eps[c] - eps[b];
              const X2 = mo[a * N3 + i * N2 + k * N + c] - mo[a * N3 + c * N2 + k * N + i];
              const X3 = mo[c * N3 + j * N2 + b * N + k] - mo[c * N3 + k * N2 + b * N + j];
              E += X1 * X2 * X3 / (D1 * D2);
            }
        }
  return E;
}

// ---------------------------------------------------------------------------
// UHF-MP3 — spatial orbital (chemist notation) formulas
// ---------------------------------------------------------------------------

/** Compute MP3 correlation energy from converged UHF results.
 *  Uses spatial-orbital formulas with explicit spin channels (αα, ββ, αβ). */
export function computeUMP3Energy(
  Ca: Matrix, Cb: Matrix,
  epsilonA: FloatArray, epsilonB: FloatArray,
  eri: ERIStored,
  noccA: number, noccB: number, nbasis: number,
  onProgress?: (msg: string) => void,
): { mp2: number; mp3: number } {
  const N = nbasis;
  const N3 = N * N * N, N2 = N * N;
  const oA = noccA, oB = noccB;
  const eA = epsilonA, eB = epsilonB;

  onProgress?.('UMP3: MO integral transformation (3 sets)...');
  const moAA = fullMOTransform(Ca, Ca, Ca, Ca, eri, N);
  const moBB = fullMOTransform(Cb, Cb, Cb, Cb, eri, N);
  const moAB = fullMOTransform(Ca, Ca, Cb, Cb, eri, N);

  // Chemist-notation accessors:
  //   aa(p,q,r,s) = (pα qα|rα sα) = moAA[p·N³+q·N²+r·N+s]
  //   bb(p,q,r,s) = (pβ qβ|rβ sβ) = moBB[p·N³+q·N²+r·N+s]
  //   ab(p,q,r,s) = (pα qα|rβ sβ) = moAB[p·N³+q·N²+r·N+s]
  //   ba(p,q,r,s) = (pβ qβ|rα sα) = moAB[r·N³+s·N²+p·N+q]  (by (pq|rs)=(rs|pq))
  const _aa = (p: number, q: number, r: number, s: number) => moAA[p * N3 + q * N2 + r * N + s];
  const _bb = (p: number, q: number, r: number, s: number) => moBB[p * N3 + q * N2 + r * N + s];
  const _ab = (p: number, q: number, r: number, s: number) => moAB[p * N3 + q * N2 + r * N + s];
  const _ba = (p: number, q: number, r: number, s: number) => moAB[r * N3 + s * N2 + p * N + q];

  // ---- MP2 energy (spatial orbital, 3 spin channels) ----
  let Emp2 = 0;
  // αα: 0.5 Σ (aibj)[(iajb)-(ibja)] / D
  for (let i = 0; i < oA; i++)
    for (let j = 0; j < oA; j++)
      for (let a = oA; a < N; a++)
        for (let b = oA; b < N; b++) {
          const D = eA[i] + eA[j] - eA[a] - eA[b];
          Emp2 += 0.5 * _aa(a, i, b, j) * (_aa(i, a, j, b) - _aa(i, b, j, a)) / D;
        }
  // ββ: 0.5 Σ (aibj)[(iajb)-(ibja)] / D
  for (let i = 0; i < oB; i++)
    for (let j = 0; j < oB; j++)
      for (let a = oB; a < N; a++)
        for (let b = oB; b < N; b++) {
          const D = eB[i] + eB[j] - eB[a] - eB[b];
          Emp2 += 0.5 * _bb(a, i, b, j) * (_bb(i, a, j, b) - _bb(i, b, j, a)) / D;
        }
  // αβ: Σ |(iβ aβ|jα bα)|² / D
  for (let i = 0; i < oB; i++)
    for (let j = 0; j < oA; j++)
      for (let a = oB; a < N; a++)
        for (let b = oA; b < N; b++) {
          const D = eB[i] + eA[j] - eB[a] - eA[b];
          const v = _ba(i, a, j, b);
          Emp2 += v * v / D;
        }

  onProgress?.(`UMP3: MP2 = ${Emp2.toFixed(10)}`);
  onProgress?.('UMP3: Computing third-order correction...');

  // ---- PP channel ----
  // αα + ββ (same-spin antisymmetric)
  let Epp = ppSameSpin(moAA, eA, oA, N) + ppSameSpin(moBB, eB, oB, N);
  // αβ: Σ (iβaβ|jαbα)(aβcβ|bαdα)(cβiβ|dαjα) / (D_βα D_βα)
  for (let i = 0; i < oB; i++)
    for (let j = 0; j < oA; j++)
      for (let a = oB; a < N; a++)
        for (let b = oA; b < N; b++) {
          const D1 = eB[i] + eA[j] - eB[a] - eA[b];
          const v1 = _ba(i, a, j, b);
          if (Math.abs(v1) < 1e-15) continue;
          for (let c = oB; c < N; c++)
            for (let d = oA; d < N; d++) {
              const D2 = eB[i] + eA[j] - eB[c] - eA[d];
              Epp += v1 * _ba(a, c, b, d) * _ba(c, i, d, j) / (D1 * D2);
            }
        }

  // ---- HH channel ----
  // αα + ββ (same-spin antisymmetric)
  let Ehh = hhSameSpin(moAA, eA, oA, N) + hhSameSpin(moBB, eB, oB, N);
  // αβ: Σ (iβaβ|jαbα)(kβiβ|lαjα)(aβkβ|bαlα) / (D_βα D_βα)
  for (let i = 0; i < oB; i++)
    for (let j = 0; j < oA; j++)
      for (let a = oB; a < N; a++)
        for (let b = oA; b < N; b++) {
          const D1 = eB[i] + eA[j] - eB[a] - eA[b];
          const v1 = _ba(i, a, j, b);
          if (Math.abs(v1) < 1e-15) continue;
          for (let k = 0; k < oB; k++)
            for (let l = 0; l < oA; l++) {
              const D2 = eB[k] + eA[l] - eB[a] - eA[b];
              Ehh += v1 * _ba(k, i, l, j) * _ba(a, k, b, l) / (D1 * D2);
            }
        }

  // ---- PH channel ----
  // Same-spin αα + ββ: Σ [(iajb)-(ibja)][(aikc)-(acki)][(cjbk)-(ckbj)] / (D1 D2)
  let Eph = phSameSpin(moAA, eA, oA, N) + phSameSpin(moBB, eB, oB, N);

  // PH group 3: i,j∈occA; k∈occB; a,b∈virA; c∈virB
  //   -[(iajb)_aa - (ibja)_aa] × (aikc)_ab × (ckbj)_ba / (D_αα D_βα)
  for (let i = 0; i < oA; i++)
    for (let j = 0; j < oA; j++) {
      for (let a = oA; a < N; a++)
        for (let b = oA; b < N; b++) {
          const D1 = eA[i] + eA[j] - eA[a] - eA[b];
          const X = _aa(i, a, j, b) - _aa(i, b, j, a);
          if (Math.abs(X) < 1e-15) continue;
          for (let k = 0; k < oB; k++)
            for (let c = oB; c < N; c++) {
              const D2 = eB[k] + eA[j] - eB[c] - eA[b];
              Eph -= X * _ab(a, i, k, c) * _ba(c, k, b, j) / (D1 * D2);
            }
        }
    }

  // PH group 4: i∈occB; j,k∈occA; a∈virB; b,c∈virA
  //   (iajb)_ba × (aikc)_ba × [(cjbk)_aa - (ckbj)_aa] / (D_βα D_αα)
  for (let i = 0; i < oB; i++)
    for (let j = 0; j < oA; j++)
      for (let a = oB; a < N; a++)
        for (let b = oA; b < N; b++) {
          const D1 = eB[i] + eA[j] - eB[a] - eA[b];
          const v1 = _ba(i, a, j, b);
          if (Math.abs(v1) < 1e-15) continue;
          for (let k = 0; k < oA; k++)
            for (let c = oA; c < N; c++) {
              const D2 = eA[k] + eA[j] - eA[c] - eA[b];
              const X3 = _aa(c, j, b, k) - _aa(c, k, b, j);
              Eph += v1 * _ba(a, i, k, c) * X3 / (D1 * D2);
            }
        }

  // PH group 5: i,k∈occA; j∈occB; a,c∈virA; b∈virB
  //   (iajb)_ab × [(acki)_aa - (aikc)_aa] × (ckbj)_ab / (D_αβ D_αβ)
  for (let i = 0; i < oA; i++)
    for (let j = 0; j < oB; j++)
      for (let a = oA; a < N; a++)
        for (let b = oB; b < N; b++) {
          const D1 = eA[i] + eB[j] - eA[a] - eB[b];
          const v1 = _ab(i, a, j, b);
          if (Math.abs(v1) < 1e-15) continue;
          for (let k = 0; k < oA; k++)
            for (let c = oA; c < N; c++) {
              const D2 = eA[k] + eB[j] - eA[c] - eB[b];
              const X2 = _aa(a, c, k, i) - _aa(a, i, k, c);
              Eph += v1 * X2 * _ab(c, k, b, j) / (D1 * D2);
            }
        }

  // PH group 6: i∈occA; j,k∈occB; a∈virA; b,c∈virB
  //   (iajb)_ab × (aikc)_ab × [(cjbk)_bb - (ckbj)_bb] / (D_αβ D_ββ)
  for (let i = 0; i < oA; i++)
    for (let j = 0; j < oB; j++)
      for (let a = oA; a < N; a++)
        for (let b = oB; b < N; b++) {
          const D1 = eA[i] + eB[j] - eA[a] - eB[b];
          const v1 = _ab(i, a, j, b);
          if (Math.abs(v1) < 1e-15) continue;
          for (let k = 0; k < oB; k++)
            for (let c = oB; c < N; c++) {
              const D2 = eB[k] + eB[j] - eB[c] - eB[b];
              const X3 = _bb(c, j, b, k) - _bb(c, k, b, j);
              Eph += v1 * _ab(a, i, k, c) * X3 / (D1 * D2);
            }
        }

  // PH group 7: i,k∈occB; j∈occA; a,c∈virB; b∈virA
  //   (iajb)_ba × [(acki)_bb - (aikc)_bb] × (ckbj)_ba / (D_βα D_βα)
  for (let i = 0; i < oB; i++)
    for (let j = 0; j < oA; j++)
      for (let a = oB; a < N; a++)
        for (let b = oA; b < N; b++) {
          const D1 = eB[i] + eA[j] - eB[a] - eA[b];
          const v1 = _ba(i, a, j, b);
          if (Math.abs(v1) < 1e-15) continue;
          for (let k = 0; k < oB; k++)
            for (let c = oB; c < N; c++) {
              const D2 = eB[k] + eA[j] - eB[c] - eA[b];
              const X2 = _bb(a, c, k, i) - _bb(a, i, k, c);
              Eph += v1 * X2 * _ba(c, k, b, j) / (D1 * D2);
            }
        }

  // PH group 8: i,j∈occB; k∈occA; a,b∈virB; c∈virA
  //   -[(iajb)_bb - (ibja)_bb] × (aikc)_ba × (ckbj)_ab / (D_ββ D_αβ)
  for (let i = 0; i < oB; i++)
    for (let j = 0; j < oB; j++)
      for (let a = oB; a < N; a++)
        for (let b = oB; b < N; b++) {
          const D1 = eB[i] + eB[j] - eB[a] - eB[b];
          const X = _bb(i, a, j, b) - _bb(i, b, j, a);
          if (Math.abs(X) < 1e-15) continue;
          for (let k = 0; k < oA; k++)
            for (let c = oA; c < N; c++) {
              const D2 = eA[k] + eB[j] - eA[c] - eB[b];
              Eph -= X * _ba(a, i, k, c) * _ab(c, k, b, j) / (D1 * D2);
            }
        }

  // PH group 9: i,k∈occA; j∈occB; a,c∈virB; b∈virA  (single term)
  //   +(ibja)_ab × (acki)_ba × (cjbk)_ba / (D_{iα jβ}^{aβ bα} × D_{kα jβ}^{cβ bα})
  for (let i = 0; i < oA; i++)
    for (let j = 0; j < oB; j++)
      for (let a = oB; a < N; a++)
        for (let b = oA; b < N; b++) {
          const D1 = eA[i] + eB[j] - eB[a] - eA[b];
          const v1 = _ab(i, b, j, a);  // (iα bα|jβ aβ)
          if (Math.abs(v1) < 1e-15) continue;
          for (let k = 0; k < oA; k++)
            for (let c = oB; c < N; c++) {
              const D2 = eA[k] + eB[j] - eB[c] - eA[b];
              // (aβ cβ|kα iα) = ba(a,c,k,i) = ab(k,i,a,c)
              // (cβ jβ|bα kα) = ba(c,j,b,k) = ab(b,k,c,j)
              Eph += v1 * _ab(k, i, a, c) * _ab(b, k, c, j) / (D1 * D2);
            }
        }

  // PH group 10: i,k∈occB; j∈occA; a,c∈virA; b∈virB  (single term)
  //   +(ibja)_ba × (acki)_ab × (cjbk)_ab / (D_{iβ jα}^{aα bβ} × D_{kβ jα}^{cα bβ})
  for (let i = 0; i < oB; i++)
    for (let j = 0; j < oA; j++)
      for (let a = oA; a < N; a++)
        for (let b = oB; b < N; b++) {
          const D1 = eB[i] + eA[j] - eA[a] - eB[b];
          const v1 = _ba(i, b, j, a);  // (iβ bβ|jα aα)
          if (Math.abs(v1) < 1e-15) continue;
          for (let k = 0; k < oB; k++)
            for (let c = oA; c < N; c++) {
              const D2 = eB[k] + eA[j] - eA[c] - eB[b];
              // (aα cα|kβ iβ) = ab(a,c,k,i)
              // (cα jα|bβ kβ) = ab(c,j,b,k)
              Eph += v1 * _ab(a, c, k, i) * _ab(c, j, b, k) / (D1 * D2);
            }
        }

  const Emp3 = Epp + Ehh + Eph;
  onProgress?.(`UMP3: pp=${Epp.toFixed(10)}, hh=${Ehh.toFixed(10)}, ph=${Eph.toFixed(10)}`);
  onProgress?.(`UMP3 correction energy: ${Emp3.toFixed(10)} Hartree`);

  return { mp2: Emp2, mp3: Emp3 };
}

// ---------------------------------------------------------------------------
// ROHF-MP3 (semicanonical orbitals)
// ---------------------------------------------------------------------------

/** Compute MP3 correlation energy from converged ROHF results.
 *  Builds semicanonical orbitals (same pattern as ROCCSD/ROMP2),
 *  then delegates to UMP3. */
export function computeROMP3Energy(
  C: Matrix, Fa: Matrix, Fb: Matrix,
  eri: ERIStored,
  noccA: number, noccB: number, nbasis: number,
  onProgress?: (msg: string) => void,
): { mp2: number; mp3: number } {
  const N = nbasis;

  onProgress?.('ROMP3: Building semicanonical orbitals...');

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

  return computeUMP3Energy(CaSc, CbSc, epsilonA, epsilonB, eri, noccA, noccB, nbasis, onProgress);
}
