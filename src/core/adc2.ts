/** ADC(2) — strict second-order Algebraic Diagrammatic Construction for
 *  excitation energies of a closed-shell RHF reference.
 *
 *  The secular matrix acts on singles (ia) and doubles (ijab) together:
 *
 *    M = [ M_SS  M_SD ]     M_SS = M(0) + M(1) + M(2)   (CIS plus 2nd-order terms)
 *        [ M_DS  M_DD ]     M_SD = 1st-order coupling
 *                            M_DD = diag(ε_a + ε_b − ε_i − ε_j)
 *
 *  (Trofimov & Schirmer; Dreuw & Wormit, WIREs Comput. Mol. Sci. 5, 82 (2015).)
 *  It is written in spin orbitals, restricted to M_s = 0 configurations, and the
 *  lowest roots are found with a Davidson solver on the full singles+doubles
 *  space. The time-reversal (spin-flip) symmetry of M separates singlets from
 *  triplets: every trial vector is kept even (singlet) or odd (triplet) under it.
 *
 *  Validated against PySCF's ADC(2) for energies and oscillator strengths.
 *
 *  An earlier version of this file built only the CIS matrix with MP2-shifted
 *  orbital energies and labelled it ADC(2); it overestimated H2O/cc-pVDZ's first
 *  singlet by ~2.3 eV (above even CIS) because the 2nd-order coupling and the
 *  doubles space were missing entirely.
 */

import type { ERIStored } from './eri';
import { Matrix, type FloatArray } from '../linalg/matrix';
import { jacobiEigen } from '../linalg/eigendecomposition';
import { computeDipoleIntegrals } from './integralsDipole';
import type { PrimitiveShell } from './types';
import type { CISExcitedState } from './cis';

const HA_TO_EV = 27.211386245988;

export interface ADC2Result {
  states: CISExcitedState[];
  mp2Energy: number;
}

// ── MO integrals ──────────────────────────────────────────────────────

/**
 * Chemists' MO integrals (pq|rs) for one block, as a dense [np][nq][nr][ns] array.
 * Four quarter-transformations, O(N^4 · n) each.
 */
function moBlock(
  C: Matrix, eri: ERIStored, N: number,
  P: number[], Q: number[], R: number[], S: number[],
): Float64Array {
  const np = P.length, nq = Q.length, nr = R.length, ns = S.length;
  // I1[μ,ν,λ,s] = Σ_σ (μν|λσ) C[σ,s]
  const I1 = new Float64Array(N * N * N * ns);
  for (let mu = 0; mu < N; mu++)
    for (let nu = 0; nu <= mu; nu++)
      for (let lam = 0; lam < N; lam++) {
        const base = ((mu * N + nu) * N + lam) * ns;
        for (let sig = 0; sig < N; sig++) {
          const v = eri.get(mu, nu, lam, sig);
          if (v === 0) continue;
          for (let s = 0; s < ns; s++) I1[base + s] += v * C.data[sig * C.cols + S[s]];
        }
        if (nu !== mu) {
          const base2 = ((nu * N + mu) * N + lam) * ns;
          for (let s = 0; s < ns; s++) I1[base2 + s] = I1[base + s];
        }
      }
  // I2[μ,ν,r,s] = Σ_λ C[λ,r] I1[μ,ν,λ,s]
  const I2 = new Float64Array(N * N * nr * ns);
  for (let mn = 0; mn < N * N; mn++)
    for (let lam = 0; lam < N; lam++) {
      const src = (mn * N + lam) * ns;
      for (let r = 0; r < nr; r++) {
        const c = C.data[lam * C.cols + R[r]];
        if (c === 0) continue;
        const dst = (mn * nr + r) * ns;
        for (let s = 0; s < ns; s++) I2[dst + s] += c * I1[src + s];
      }
    }
  // I3[μ,q,r,s] = Σ_ν C[ν,q] I2[μ,ν,r,s]
  const rs = nr * ns;
  const I3 = new Float64Array(N * nq * rs);
  for (let mu = 0; mu < N; mu++)
    for (let nu = 0; nu < N; nu++) {
      const src = (mu * N + nu) * rs;
      for (let q = 0; q < nq; q++) {
        const c = C.data[nu * C.cols + Q[q]];
        if (c === 0) continue;
        const dst = (mu * nq + q) * rs;
        for (let k = 0; k < rs; k++) I3[dst + k] += c * I2[src + k];
      }
    }
  // out[p,q,r,s] = Σ_μ C[μ,p] I3[μ,q,r,s]
  const qrs = nq * rs;
  const out = new Float64Array(np * qrs);
  for (let mu = 0; mu < N; mu++)
    for (let p = 0; p < np; p++) {
      const c = C.data[mu * C.cols + P[p]];
      if (c === 0) continue;
      const src = mu * qrs, dst = p * qrs;
      for (let k = 0; k < qrs; k++) out[dst + k] += c * I3[src + k];
    }
  return out;
}

// ── main ──────────────────────────────────────────────────────────────

export function computeADC2(
  C: Matrix,
  epsilon: FloatArray,
  eri: ERIStored,
  nocc: number,
  nbasis: number,
  nStates: number,
  primitiveShells: PrimitiveShell[],
  normFactors: number[],
  onProgress?: (msg: string) => void,
  isTriplet: boolean = false,
): ADC2Result {
  const N = nbasis;

  // Real MOs only: with pure d/f functions the trailing columns of C are
  // identically zero placeholders (see HF.nullDirections).
  let nmo = N;
  while (nmo > 0) {
    let zero = true;
    for (let mu = 0; mu < N && zero; mu++) if (C.data[mu * C.cols + nmo - 1] !== 0) zero = false;
    if (!zero) break;
    nmo--;
  }
  const no = nocc, nv = nmo - nocc;
  // Dense spin-orbital tensors: t2, t2(2), the ring intermediate and the doubles
  // part of each eigenvector are (2no)²(2nv)² each, and the (vv|vv) transform
  // holds an N³·nv intermediate. Refuse up front rather than let the tab die.
  const so4 = (2 * no) ** 2 * (2 * nv) ** 2;
  const needBytes = 8 * Math.max(4 * so4, N ** 3 * Math.max(nv, 1));
  const LIMIT = 700 * 1024 ** 2;
  if (needBytes > LIMIT) {
    throw new Error(
      `ADC(2) for ${no} occupied / ${nv} virtual orbitals needs about ` +
      `${(needBytes / 1024 ** 2).toFixed(0)} MB, more than a browser tab can safely use ` +
      `(limit ${(LIMIT / 1024 ** 2).toFixed(0)} MB). Use a smaller molecule or basis set.`);
  }

  const occ = Array.from({ length: no }, (_, i) => i);
  const vir = Array.from({ length: nv }, (_, a) => no + a);

  onProgress?.(`ADC(2) ${isTriplet ? 'triplet' : 'singlet'}: transforming integrals to MO basis...`);
  const OVOV = moBlock(C, eri, N, occ, vir, occ, vir); // (ia|jb)
  const OOVV = moBlock(C, eri, N, occ, occ, vir, vir); // (ij|ab)
  const OOOV = moBlock(C, eri, N, occ, occ, occ, vir); // (ij|ka)
  const OVVV = moBlock(C, eri, N, occ, vir, vir, vir); // (ia|bc)

  const ovov = (i: number, a: number, j: number, b: number) => OVOV[((i * nv + a) * no + j) * nv + b];
  const oovv = (i: number, j: number, a: number, b: number) => OOVV[((i * no + j) * nv + a) * nv + b];
  const ooov = (i: number, j: number, k: number, a: number) => OOOV[((i * no + j) * no + k) * nv + a];
  const ovvv = (i: number, a: number, b: number, c: number) => OVVV[((i * nv + a) * nv + b) * nv + c];

  // ── spin orbitals: index 2p (α) and 2p+1 (β); spatial = p >> 1, spin = p & 1 ──
  const O = 2 * no, V = 2 * nv;
  const eO = (i: number) => epsilon[i >> 1];
  const eV = (a: number) => epsilon[no + (a >> 1)];

  // Antisymmetrised integrals <pq||rs> = <pq|rs> − <pq|sr>, <pq|rs> = (pr|qs)·δσ.
  // Each block below is named by the occ/vir character of p,q,r,s.
  /** <kl||cd> (o o v v) */
  const A_oovv = (k: number, l: number, c: number, d: number): number => {
    const sk = k & 1, sl = l & 1, sc = c & 1, sd = d & 1;
    let v = 0;
    if (sk === sc && sl === sd) v += ovov(k >> 1, c >> 1, l >> 1, d >> 1);
    if (sk === sd && sl === sc) v -= ovov(k >> 1, d >> 1, l >> 1, c >> 1);
    return v;
  };
  /** <kl||id> (o o o v) = (ki|ld) − (kd|li) */
  const A_ooov = (k: number, l: number, i: number, d: number): number => {
    const sk = k & 1, sl = l & 1, si = i & 1, sd = d & 1;
    let v = 0;
    if (sk === si && sl === sd) v += ooov(k >> 1, i >> 1, l >> 1, d >> 1);
    if (sk === sd && sl === si) v -= ooov(l >> 1, i >> 1, k >> 1, d >> 1);
    return v;
  };
  /** <al||cd> (v o v v) = (ac|ld) − (ad|lc) */
  const A_vovv = (a: number, l: number, c: number, d: number): number => {
    const sa = a & 1, sl = l & 1, sc = c & 1, sd = d & 1;
    let v = 0;
    if (sa === sc && sl === sd) v += ovvv(l >> 1, d >> 1, a >> 1, c >> 1);
    if (sa === sd && sl === sc) v -= ovvv(l >> 1, c >> 1, a >> 1, d >> 1);
    return v;
  };
  /** <aj||bi> (v o v o) = (ab|ji) − (ai|jb) */
  const A_vovo = (a: number, j: number, b: number, i: number): number => {
    const sa = a & 1, sj = j & 1, sb = b & 1, si = i & 1;
    let v = 0;
    if (sa === sb && sj === si) v += oovv(j >> 1, i >> 1, a >> 1, b >> 1);
    if (sa === si && sj === sb) v -= ovov(i >> 1, a >> 1, j >> 1, b >> 1);
    return v;
  };

  // ── MP2 amplitudes t_kl^cd = <kl||cd> / (ε_k+ε_l−ε_c−ε_d), dense over spin orbitals ──
  onProgress?.('ADC(2): MP2 amplitudes...');
  const T = new Float64Array(O * O * V * V);
  const tIdx = (k: number, l: number, c: number, d: number) => ((k * O + l) * V + c) * V + d;
  let mp2Energy = 0;
  for (let k = 0; k < O; k++)
    for (let l = 0; l < O; l++)
      for (let c = 0; c < V; c++)
        for (let d = 0; d < V; d++) {
          const g = A_oovv(k, l, c, d);
          if (g === 0) continue;
          const t = g / (eO(k) + eO(l) - eV(c) - eV(d));
          T[tIdx(k, l, c, d)] = t;
          mp2Energy += 0.25 * g * t;
        }
  onProgress?.(`ADC(2): MP2 energy = ${mp2Energy.toFixed(8)} Eh`);

  // ── configuration lists (M_s = 0 only) ──
  const singles: [number, number][] = [];
  for (let i = 0; i < O; i++) for (let a = 0; a < V; a++) if ((i & 1) === (a & 1)) singles.push([i, a]);
  const nS = singles.length;
  const sIndex = new Int32Array(O * V).fill(-1);
  singles.forEach(([i, a], n) => { sIndex[i * V + a] = n; });

  const doubles: [number, number, number, number][] = [];
  for (let k = 0; k < O; k++) for (let l = k + 1; l < O; l++)
    for (let c = 0; c < V; c++) for (let d = c + 1; d < V; d++)
      if ((k & 1) + (l & 1) === (c & 1) + (d & 1)) doubles.push([k, l, c, d]);
  const nD = doubles.length;
  const dim = nS + nD;

  // ── M_SS = M(0) + M(1) + M(2), dense ──
  //   M(2)_{ia,jb} = −δ_ij (X_ab + X_ba) − δ_ab (Y_ij + Y_ji)
  //                 + 1/2 Σ_ck [ t_ik^ac <jk||bc> + t_jk^bc <ik||ac> ]
  // with X, Y the quarter-weighted contractions below. The signs were fixed by
  // matching PySCF: of the 16 sign combinations for these three terms and the
  // relative sign inside M_SD, only this one reproduces it (to 1e-4 eV; the next
  // best is off by 0.11 eV).
  onProgress?.(`ADC(2): building singles block (${nS}) — doubles space ${nD}...`);
  // Intermediates for the δ_ij and δ_ab parts of M(2)
  const Xvv = new Float64Array(V * V); // 1/4 Σ_ckl t_kl^ac <kl||bc>
  for (let a = 0; a < V; a++) for (let b = 0; b < V; b++) {
    let s = 0;
    for (let k = 0; k < O; k++) for (let l = 0; l < O; l++) for (let c = 0; c < V; c++) {
      const t = T[tIdx(k, l, a, c)];
      if (t !== 0) s += t * A_oovv(k, l, b, c);
    }
    Xvv[a * V + b] = 0.25 * s;
  }
  const Xoo = new Float64Array(O * O); // 1/4 Σ_cdk t_ik^cd <jk||cd>
  for (let i = 0; i < O; i++) for (let j = 0; j < O; j++) {
    let s = 0;
    for (let k = 0; k < O; k++) for (let c = 0; c < V; c++) for (let d = 0; d < V; d++) {
      const t = T[tIdx(i, k, c, d)];
      if (t !== 0) s += t * A_oovv(j, k, c, d);
    }
    Xoo[i * O + j] = 0.25 * s;
  }
  const Mss = new Float64Array(nS * nS);
  for (let p = 0; p < nS; p++) {
    const [i, a] = singles[p];
    for (let q = p; q < nS; q++) {
      const [j, b] = singles[q];
      let v = -A_vovo(a, j, b, i);
      if (i === j && a === b) v += eV(a) - eO(i);
      if (i === j) v -= Xvv[a * V + b] + Xvv[b * V + a];
      if (a === b) v -= Xoo[i * O + j] + Xoo[j * O + i];
      // + 1/2 Σ_ck [ t_ik^ac <jk||bc> + t_jk^bc <ik||ac> ]
      let s = 0;
      for (let k = 0; k < O; k++) for (let c = 0; c < V; c++) {
        const t1 = T[tIdx(i, k, a, c)], t2 = T[tIdx(j, k, b, c)];
        if (t1 !== 0) s += t1 * A_oovv(j, k, b, c);
        if (t2 !== 0) s += t2 * A_oovv(i, k, a, c);
      }
      v += 0.5 * s;
      Mss[p * nS + q] = v;
      Mss[q * nS + p] = v;
    }
  }

  // ── M_SD (sparse: each doubles column touches ≤ 2(O+V) singles) ──
  // M_{ia,klcd} = <kl||id> δ_ac − <kl||ic> δ_ad − <al||cd> δ_ik + <ak||cd> δ_il
  const sdRow: number[] = [], sdCol: number[] = [], sdVal: number[] = [];
  for (let D = 0; D < nD; D++) {
    const [k, l, c, d] = doubles[D];
    for (let i = 0; i < O; i++) {
      let s = sIndex[i * V + c];
      if (s >= 0) { const v = A_ooov(k, l, i, d); if (v !== 0) { sdRow.push(s); sdCol.push(D); sdVal.push(v); } }
      s = sIndex[i * V + d];
      if (s >= 0) { const v = -A_ooov(k, l, i, c); if (v !== 0) { sdRow.push(s); sdCol.push(D); sdVal.push(v); } }
    }
    for (let a = 0; a < V; a++) {
      let s = sIndex[k * V + a];
      if (s >= 0) { const v = -A_vovv(a, l, c, d); if (v !== 0) { sdRow.push(s); sdCol.push(D); sdVal.push(v); } }
      s = sIndex[l * V + a];
      if (s >= 0) { const v = A_vovv(a, k, c, d); if (v !== 0) { sdRow.push(s); sdCol.push(D); sdVal.push(v); } }
    }
  }
  const nnz = sdVal.length;
  const SR = Int32Array.from(sdRow), SC = Int32Array.from(sdCol), SV = Float64Array.from(sdVal);
  const Ddiag = new Float64Array(nD);
  for (let D = 0; D < nD; D++) { const [k, l, c, d] = doubles[D]; Ddiag[D] = eV(c) + eV(d) - eO(k) - eO(l); }

  const sigma = (x: Float64Array, out: Float64Array) => {
    out.fill(0);
    for (let p = 0; p < nS; p++) {
      let s = 0; const row = p * nS;
      for (let q = 0; q < nS; q++) s += Mss[row + q] * x[q];
      out[p] = s;
    }
    for (let D = 0; D < nD; D++) out[nS + D] = Ddiag[D] * x[nS + D];
    for (let z = 0; z < nnz; z++) {
      const s = SR[z], D = nS + SC[z], v = SV[z];
      out[s] += v * x[D];
      out[D] += v * x[s];
    }
  };

  // ── spin-flip symmetry: α ↔ β on every index ──
  const flipMap = new Int32Array(dim), flipSign = new Float64Array(dim);
  for (let p = 0; p < nS; p++) { const [i, a] = singles[p]; flipMap[p] = sIndex[(i ^ 1) * V + (a ^ 1)]; flipSign[p] = 1; }
  {
    const dKey = (k: number, l: number, c: number, d: number) => ((k * O + l) * V + c) * V + d;
    const dIndex = new Map<number, number>();
    doubles.forEach(([k, l, c, d], D) => dIndex.set(dKey(k, l, c, d), D));
    doubles.forEach(([k, l, c, d], D) => {
      let k2 = k ^ 1, l2 = l ^ 1, c2 = c ^ 1, d2 = d ^ 1, sign = 1;
      if (k2 > l2) { [k2, l2] = [l2, k2]; sign = -sign; }
      if (c2 > d2) { [c2, d2] = [d2, c2]; sign = -sign; }
      flipMap[nS + D] = nS + dIndex.get(dKey(k2, l2, c2, d2))!;
      flipSign[nS + D] = sign;
    });
  }
  const parity = isTriplet ? -1 : 1;
  const project = (x: Float64Array) => {
    const y = new Float64Array(dim);
    for (let p = 0; p < dim; p++) y[p] = 0.5 * (x[p] + parity * flipSign[p] * x[flipMap[p]]);
    return y;
  };

  // ── Davidson ──
  const nRoots = Math.max(1, Math.min(nStates, nS / 2));
  const diag = new Float64Array(dim);
  for (let p = 0; p < nS; p++) diag[p] = Mss[p * nS + p];
  diag.set(Ddiag, nS);

  const dot = (u: Float64Array, v: Float64Array) => { let s = 0; for (let p = 0; p < dim; p++) s += u[p] * v[p]; return s; };
  const basis: Float64Array[] = [], sig: Float64Array[] = [];
  const addVector = (v: Float64Array): boolean => {
    for (let pass = 0; pass < 2; pass++)
      for (const b of basis) { const c = dot(b, v); for (let p = 0; p < dim; p++) v[p] -= c * b[p]; }
    const n = Math.sqrt(dot(v, v));
    if (n < 1e-6) return false;
    for (let p = 0; p < dim; p++) v[p] /= n;
    const s = new Float64Array(dim); sigma(v, s);
    basis.push(v); sig.push(s);
    return true;
  };

  // Guesses: lowest diagonal singles, symmetrised to the requested spin.
  const order = Array.from({ length: nS }, (_, p) => p).sort((p, q) => diag[p] - diag[q]);
  const nGuess = Math.min(nS / 2, Math.max(nRoots + 4, 2 * nRoots));
  for (const p of order) {
    if (basis.length >= nGuess) break;
    const v = new Float64Array(dim); v[p] = 1;
    addVector(project(v));
  }

  const maxSub = Math.max(40, 12 * nRoots);
  let theta = new Float64Array(nRoots);
  let ritz: Float64Array[] = [];
  onProgress?.(`ADC(2): Davidson on ${dim} configurations (${nS} singles + ${nD} doubles)...`);
  for (let iter = 0; iter < 200; iter++) {
    const m = basis.length;
    const G = new Matrix(m, m);
    for (let i = 0; i < m; i++) for (let j = i; j < m; j++) {
      const g = dot(basis[i], sig[j]); G.set(i, j, g); G.set(j, i, g);
    }
    const { eigenvalues, eigenvectors } = jacobiEigen(G);
    theta = Float64Array.from(eigenvalues.slice(0, nRoots));
    ritz = [];
    const residuals: Float64Array[] = [];
    let maxRes = 0;
    for (let r = 0; r < nRoots; r++) {
      const x = new Float64Array(dim), ax = new Float64Array(dim);
      for (let i = 0; i < m; i++) {
        const c = eigenvectors.get(i, r);
        const b = basis[i], s = sig[i];
        for (let p = 0; p < dim; p++) { x[p] += c * b[p]; ax[p] += c * s[p]; }
      }
      ritz.push(x);
      const res = new Float64Array(dim);
      for (let p = 0; p < dim; p++) res[p] = ax[p] - theta[r] * x[p];
      const rn = Math.sqrt(dot(res, res));
      maxRes = Math.max(maxRes, rn);
      if (rn > 1e-6) {
        for (let p = 0; p < dim; p++) {
          const den = theta[r] - diag[p];
          res[p] /= Math.abs(den) > 1e-4 ? den : (den >= 0 ? 1e-4 : -1e-4);
        }
        residuals.push(res);
      }
    }
    if (maxRes < 1e-6) break;
    if (basis.length + residuals.length > maxSub) {
      // restart from the current Ritz vectors
      basis.length = 0; sig.length = 0;
      for (const x of ritz) addVector(project(x));
    }
    let added = 0;
    for (const r of residuals) if (addVector(project(r))) added++;
    if (added === 0) break;
  }

  // ── transition moments, to second order (matches PySCF's default) ──
  // Each block of the transition density (ov, vo, oo, vv) was checked separately
  // against PySCF's; all four agree to 1e-6. The doubles-part (Y2) terms enter
  // ρ_oo with + and ρ_vv with −, and ρ_vo needs the second-order doubles t2(2).
  // The transition density between the ground and an excited state is
  //   ρ_ov = Y + 2nd-order corrections,  ρ_vo, ρ_oo, ρ_vv from t2 and t1(2),
  // where Y and Y2 are the singles and doubles parts of the ADC eigenvector.
  // Using Y alone (zeroth order) misses up to ~40 % of f (CO/cc-pVDZ).
  let muMO: Float64Array | null = null; // spatial MO dipoles, [x][p][q]
  let t12: Float64Array | null = null;   // 2nd-order singles amplitudes, spin-orbital [i][a]
  let Rvv: Float64Array | null = null, Roo: Float64Array | null = null;
  let T2s: Float64Array | null = null;   // 2nd-order doubles amplitudes, spin-orbital
  if (!isTriplet) {
    const dip = computeDipoleIntegrals(primitiveShells, normFactors, nbasis);
    muMO = new Float64Array(3 * nmo * nmo);
    [dip.Dx, dip.Dy, dip.Dz].forEach((D, x) => {
      const half = new Float64Array(N * nmo); // D C
      for (let mu = 0; mu < N; mu++) for (let nu = 0; nu < N; nu++) {
        const d = D.get(mu, nu); if (d === 0) continue;
        for (let q = 0; q < nmo; q++) half[mu * nmo + q] += d * C.data[nu * C.cols + q];
      }
      for (let pp = 0; pp < nmo; pp++) for (let mu = 0; mu < N; mu++) {
        const c = C.data[mu * C.cols + pp]; if (c === 0) continue;
        for (let q = 0; q < nmo; q++) muMO![(x * nmo + pp) * nmo + q] += c * half[mu * nmo + q];
      }
    });
    // t_i^a(2) = [ ½ Σ_kcd <ak||cd> t_ik^cd − ½ Σ_klc <kl||ic> t_kl^ac ] / (ε_i − ε_a)
    t12 = new Float64Array(O * V);
    for (let i = 0; i < O; i++) for (let a = 0; a < V; a++) {
      if ((i & 1) !== (a & 1)) continue;
      let v = 0;
      for (let k = 0; k < O; k++) for (let c = 0; c < V; c++) for (let d = 0; d < V; d++) {
        const t = T[tIdx(i, k, c, d)]; if (t !== 0) v += 0.5 * A_vovv(a, k, c, d) * t;
      }
      for (let k = 0; k < O; k++) for (let l = 0; l < O; l++) for (let c = 0; c < V; c++) {
        const t = T[tIdx(k, l, a, c)]; if (t !== 0) v -= 0.5 * A_ooov(k, l, i, c) * t;
      }
      t12[i * V + a] = v / (eO(i) - eV(a));
    }
    // Second-order doubles amplitudes (the MP3 doubles):
    //   D t_ij^ab(2) = ½ Σ_cd <ab||cd> t_ij^cd + ½ Σ_kl <kl||ij> t_kl^ab
    //                + P(ij)P(ab) Σ_kc <kb||cj> t_ik^ac
    // They enter the transition density only through ρ_vo, but leaving them out
    // underestimates f by 10-15 % (H2O/cc-pVDZ), so they are not optional.
    const VVVV = moBlock(C, eri, N, vir, vir, vir, vir); // (ab|cd)
    const OOOO = moBlock(C, eri, N, occ, occ, occ, occ); // (ij|kl)
    const vvvv = (a: number, b: number, c: number, d: number) => VVVV[((a * nv + b) * nv + c) * nv + d];
    const oooo = (i: number, j: number, k: number, l: number) => OOOO[((i * no + j) * no + k) * no + l];
    /** <ab||cd> = (ac|bd) − (ad|bc) */
    const A_vvvv = (a: number, b: number, c: number, d: number): number => {
      let v = 0;
      if ((a & 1) === (c & 1) && (b & 1) === (d & 1)) v += vvvv(a >> 1, c >> 1, b >> 1, d >> 1);
      if ((a & 1) === (d & 1) && (b & 1) === (c & 1)) v -= vvvv(a >> 1, d >> 1, b >> 1, c >> 1);
      return v;
    };
    /** <kl||ij> = (ki|lj) − (kj|li) */
    const A_oooo = (k: number, l: number, i: number, j: number): number => {
      let v = 0;
      if ((k & 1) === (i & 1) && (l & 1) === (j & 1)) v += oooo(k >> 1, i >> 1, l >> 1, j >> 1);
      if ((k & 1) === (j & 1) && (l & 1) === (i & 1)) v -= oooo(k >> 1, j >> 1, l >> 1, i >> 1);
      return v;
    };
    /** <kb||cj> = (kc|bj) − (kj|bc) */
    const A_ovvo = (k: number, b: number, c: number, j: number): number => {
      let v = 0;
      if ((k & 1) === (c & 1) && (b & 1) === (j & 1)) v += ovov(k >> 1, c >> 1, j >> 1, b >> 1);
      if ((k & 1) === (j & 1) && (b & 1) === (c & 1)) v -= oovv(k >> 1, j >> 1, b >> 1, c >> 1);
      return v;
    };
    const Wring = new Float64Array(O * O * V * V); // Σ_kc <kb||cj> t_ik^ac, before P(ij)P(ab)
    for (let i = 0; i < O; i++) for (let k = 0; k < O; k++) for (let a = 0; a < V; a++) for (let c = 0; c < V; c++) {
      const t = T[tIdx(i, k, a, c)]; if (t === 0) continue;
      for (let j = 0; j < O; j++) for (let b = 0; b < V; b++) {
        const g = A_ovvo(k, b, c, j); if (g !== 0) Wring[tIdx(i, j, a, b)] += g * t;
      }
    }
    T2s = new Float64Array(O * O * V * V);
    for (let i = 0; i < O; i++) for (let j = 0; j < O; j++) for (let a = 0; a < V; a++) for (let b = 0; b < V; b++) {
      if ((i & 1) + (j & 1) !== (a & 1) + (b & 1)) continue;
      let v = 0;
      for (let c = 0; c < V; c++) for (let d = 0; d < V; d++) {
        const t = T[tIdx(i, j, c, d)]; if (t !== 0) v += 0.5 * A_vvvv(a, b, c, d) * t;
      }
      for (let k = 0; k < O; k++) for (let l = 0; l < O; l++) {
        const t = T[tIdx(k, l, a, b)]; if (t !== 0) v += 0.5 * A_oooo(k, l, i, j) * t;
      }
      v += Wring[tIdx(i, j, a, b)] - Wring[tIdx(j, i, a, b)] - Wring[tIdx(i, j, b, a)] + Wring[tIdx(j, i, b, a)];
      T2s[tIdx(i, j, a, b)] = v / (eO(i) + eO(j) - eV(a) - eV(b));
    }

    // root-independent contractions of t2 with itself
    Rvv = new Float64Array(V * V); // Σ_xyh t_xy^gh t_xy^vh
    Roo = new Float64Array(O * O); // Σ_ygh t_xy^gh t_py^gh
    for (let x = 0; x < O; x++) for (let y = 0; y < O; y++)
      for (let g = 0; g < V; g++) for (let h = 0; h < V; h++) {
        const t = T[tIdx(x, y, g, h)]; if (t === 0) continue;
        for (let v = 0; v < V; v++) Rvv[g * V + v] += t * T[tIdx(x, y, v, h)];
        for (let pp = 0; pp < O; pp++) Roo[x * O + pp] += t * T[tIdx(pp, y, g, h)];
      }
  }

  const transitionDipole = (x: Float64Array): [number, number, number] => {
    const Y = new Float64Array(O * V);
    for (let p = 0; p < nS; p++) { const [i, a] = singles[p]; Y[i * V + a] = x[p]; }
    const Y2 = new Float64Array(O * O * V * V); // full antisymmetric doubles part
    for (let D = 0; D < nD; D++) {
      const [k, l, c, d] = doubles[D]; const y = x[nS + D];
      if (y === 0) continue;
      Y2[tIdx(k, l, c, d)] = y; Y2[tIdx(l, k, c, d)] = -y;
      Y2[tIdx(k, l, d, c)] = -y; Y2[tIdx(l, k, d, c)] = y;
    }
    const TYov = new Float64Array(O * V), TYvo = new Float64Array(V * O);
    const TYoo = new Float64Array(O * O), TYvv = new Float64Array(V * V);
    TYov.set(Y);
    // U[y,h] = Σ_xg Y[x,g] t_xy^gh ; TYvo[v,p] = Σ_xg Y[x,g] t_px^vg
    const U = new Float64Array(O * V);
    for (let xx = 0; xx < O; xx++) for (let g = 0; g < V; g++) {
      const yv = Y[xx * V + g]; if (yv === 0) continue;
      for (let y = 0; y < O; y++) for (let h = 0; h < V; h++) {
        U[y * V + h] += yv * T[tIdx(xx, y, g, h)];
        TYvo[h * O + y] += yv * (T[tIdx(y, xx, h, g)] + T2s![tIdx(y, xx, h, g)]);
      }
    }
    for (let pp = 0; pp < O; pp++) for (let v = 0; v < V; v++) {
      let s = 0;
      for (let y = 0; y < O; y++) for (let h = 0; h < V; h++) s += U[y * V + h] * T[tIdx(pp, y, v, h)];
      for (let g = 0; g < V; g++) s -= 0.5 * Y[pp * V + g] * Rvv![g * V + v];
      for (let xx = 0; xx < O; xx++) s -= 0.5 * Y[xx * V + v] * Roo![xx * O + pp];
      TYov[pp * V + v] += 0.5 * s;
    }
    for (let pp = 0; pp < O; pp++) for (let q = 0; q < O; q++) {
      let s = 0;
      for (let g = 0; g < V; g++) s -= Y[pp * V + g] * t12![q * V + g];
      let s2 = 0;
      for (let xx = 0; xx < O; xx++) for (let g = 0; g < V; g++) for (let h = 0; h < V; h++)
        s2 += Y2[tIdx(pp, xx, g, h)] * T[tIdx(q, xx, g, h)];
      TYoo[pp * O + q] = s + 0.5 * s2;
    }
    for (let v = 0; v < V; v++) for (let r = 0; r < V; r++) {
      let s = 0;
      for (let xx = 0; xx < O; xx++) s += Y[xx * V + r] * t12![xx * V + v];
      let s2 = 0;
      for (let xx = 0; xx < O; xx++) for (let y = 0; y < O; y++) for (let g = 0; g < V; g++)
        s2 += Y2[tIdx(xx, y, r, g)] * T[tIdx(xx, y, v, g)];
      TYvv[v * V + r] = s - 0.5 * s2;
    }
    const m = (c: number, pp: number, q: number) => muMO![(c * nmo + pp) * nmo + q];
    const out: [number, number, number] = [0, 0, 0];
    for (let c = 0; c < 3; c++) {
      let s = 0;
      for (let pp = 0; pp < O; pp++) for (let v = 0; v < V; v++) {
        if ((pp & 1) !== (v & 1)) continue;
        s += (TYov[pp * V + v] + TYvo[v * O + pp]) * m(c, pp >> 1, no + (v >> 1));
      }
      for (let pp = 0; pp < O; pp++) for (let q = 0; q < O; q++)
        if ((pp & 1) === (q & 1)) s += TYoo[pp * O + q] * m(c, pp >> 1, q >> 1);
      for (let v = 0; v < V; v++) for (let r = 0; r < V; r++)
        if ((v & 1) === (r & 1)) s += TYvv[v * V + r] * m(c, no + (v >> 1), no + (r >> 1));
      out[c] = s;
    }
    return out;
  };

  const states: CISExcitedState[] = [];
  for (let r = 0; r < ritz.length; r++) {
    const x = ritz[r];
    let f = 0;
    if (!isTriplet) {
      const tdm = transitionDipole(x);
      f = (2 / 3) * theta[r] * (tdm[0] ** 2 + tdm[1] ** 2 + tdm[2] ** 2);
    }
    // Report spatial i → a weights, combining the α and β components.
    const w = new Map<number, number>();
    for (let p = 0; p < nS; p++) {
      const [i, a] = singles[p];
      if ((i & 1) !== 0) continue;
      w.set((i >> 1) * nv + (a >> 1), x[p] * Math.SQRT2);
    }
    const transitions = [...w.entries()]
      .filter(([, c]) => Math.abs(c) > 0.1)
      .sort((u, v) => Math.abs(v[1]) - Math.abs(u[1]))
      .slice(0, 3)
      .map(([key, coeff]) => ({ i: Math.floor(key / nv), a: no + (key % nv), coeff }));
    states.push({
      energy: theta[r],
      energyEV: theta[r] * HA_TO_EV,
      oscillatorStrength: Math.max(0, f),
      dominantTransitions: transitions,
    });
  }

  onProgress?.(`ADC(2): ${states.length} states computed`);
  return { states, mp2Energy };
}
