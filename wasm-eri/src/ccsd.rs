/// WASM-accelerated CCSD correlation energy.
///
/// RHF: Spatial-orbital formulation (Crawdad notation, ~16x faster than spin-orbital).
/// UHF: Spin-orbital with blocked ordering [occα, occβ, virα, virβ].
///
/// Reference: Stanton, Gauss, Watts, Bartlett (JCP 94, 4334, 1991).

use crate::mp3::{full_mo_transform, full_mo_transform_4c};
use crate::simd_utils;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

#[inline(always)]
fn t2idx(i: usize, j: usize, a: usize, b: usize, nocc: usize, nvir: usize) -> usize {
    ((i * nocc + j) * nvir + a) * nvir + b
}

// ---------------------------------------------------------------------------
// DIIS
// ---------------------------------------------------------------------------

struct Diis {
    max_size: usize,
    min_size: usize,
    amp_history: Vec<Vec<f64>>,
    err_history: Vec<Vec<f64>>,
}

impl Diis {
    fn new(max_size: usize, min_size: usize) -> Self {
        Diis {
            max_size,
            min_size,
            amp_history: Vec::new(),
            err_history: Vec::new(),
        }
    }

    fn push(&mut self, amp: &[f64], err: &[f64]) {
        self.amp_history.push(amp.to_vec());
        self.err_history.push(err.to_vec());
        if self.amp_history.len() > self.max_size {
            self.amp_history.remove(0);
            self.err_history.remove(0);
        }
    }

    fn can_extrapolate(&self) -> bool {
        self.amp_history.len() >= self.min_size
    }

    fn extrapolate(&self) -> Vec<f64> {
        let m = self.amp_history.len();
        let n = m + 1;
        let vec_size = self.amp_history[0].len();

        let mut b_mat = vec![0.0_f64; n * n];
        let mut rhs = vec![0.0_f64; n];

        for i in 0..m {
            for j in i..m {
                let dot = simd_utils::dot(&self.err_history[i], &self.err_history[j]);
                b_mat[i * n + j] = dot;
                b_mat[j * n + i] = dot;
            }
            b_mat[i * n + m] = -1.0;
            b_mat[m * n + i] = -1.0;
        }
        b_mat[m * n + m] = 0.0;
        rhs[m] = -1.0;

        let c = solve_linear_system(&b_mat, &rhs, n);

        let mut result = vec![0.0_f64; vec_size];
        for i in 0..m {
            simd_utils::daxpy(c[i], &self.amp_history[i], &mut result);
        }
        result
    }
}

fn solve_linear_system(a_in: &[f64], b_in: &[f64], n: usize) -> Vec<f64> {
    let mut a = a_in.to_vec();
    let mut x = b_in.to_vec();

    for col in 0..n {
        let mut max_val = a[col * n + col].abs();
        let mut max_row = col;
        for row in (col + 1)..n {
            let v = a[row * n + col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }
        if max_row != col {
            for j in 0..n {
                a.swap(col * n + j, max_row * n + j);
            }
            x.swap(col, max_row);
        }
        let pivot = a[col * n + col];
        if pivot.abs() < 1e-30 {
            continue;
        }
        for row in (col + 1)..n {
            let factor = a[row * n + col] / pivot;
            for j in col..n {
                a[row * n + j] -= factor * a[col * n + j];
            }
            x[row] -= factor * x[col];
        }
    }

    for row in (0..n).rev() {
        for j in (row + 1)..n {
            x[row] -= a[row * n + j] * x[j];
        }
        let diag = a[row * n + row];
        if diag.abs() > 1e-30 {
            x[row] /= diag;
        } else {
            x[row] = 0.0;
        }
    }
    x
}

// ---------------------------------------------------------------------------
// RHF spatial-orbital CCSD (Crawdad notation)
// ---------------------------------------------------------------------------
// v(p,q,r,s) = (pr|qs)   — Coulomb integral
// w(p,q,r,s) = 2*v(p,q,r,s) - v(p,q,s,r)  — antisymmetrized
// P(ia,jb) f = f(i,a,j,b) + f(j,b,i,a)

fn ccsd_spatial_loop(
    mo_ints: &[f64],
    eps: &[f64],
    nocc: usize,
    n: usize,
) -> f64 {
    let nvir = n - nocc;
    let n2 = n * n;
    let n3 = n2 * n;

    // v(p,q,r,s) = (pr|qs) = mo_ints[p*N³ + r*N² + q*N + s]
    let v = |p: usize, q: usize, r: usize, s: usize| -> f64 {
        mo_ints[p * n3 + r * n2 + q * n + s]
    };
    let w = |p: usize, q: usize, r: usize, s: usize| -> f64 {
        2.0 * v(p, q, r, s) - v(p, q, s, r)
    };

    let t2i = |i: usize, j: usize, a: usize, b: usize| -> usize {
        ((i * nocc + j) * nvir + a) * nvir + b
    };

    // Denominators
    let t1_size = nocc * nvir;
    let t2_size = nocc * nocc * nvir * nvir;

    let mut dia = vec![0.0_f64; t1_size];
    for i in 0..nocc {
        for a in 0..nvir {
            dia[i * nvir + a] = eps[i] - eps[nocc + a];
        }
    }

    let mut dijab = vec![0.0_f64; t2_size];
    for i in 0..nocc {
        for j in 0..nocc {
            for a in 0..nvir {
                for b in 0..nvir {
                    dijab[t2i(i, j, a, b)] = eps[i] + eps[j] - eps[nocc + a] - eps[nocc + b];
                }
            }
        }
    }

    // MP2 initial guess: T1=0, t2(i,j,a,b) = v(i,j,A,B) / D
    let mut t1 = vec![0.0_f64; t1_size];
    let mut t2 = vec![0.0_f64; t2_size];
    for i in 0..nocc {
        for j in 0..nocc {
            for a in 0..nvir {
                for b in 0..nvir {
                    let idx = t2i(i, j, a, b);
                    t2[idx] = v(i, j, nocc + a, nocc + b) / dijab[idx];
                }
            }
        }
    }

    // Energy: E = Σ_{ijab} w(i,j,A,B) * (t2(i,j,a,b) + t1(i,a)*t1(j,b))
    let energy = |t1: &[f64], t2: &[f64]| -> f64 {
        let mut e = 0.0;
        for i in 0..nocc {
            for j in 0..nocc {
                for a in 0..nvir {
                    for b in 0..nvir {
                        let aa = nocc + a;
                        let bb = nocc + b;
                        e += w(i, j, aa, bb) * (t2[t2i(i, j, a, b)] + t1[i * nvir + a] * t1[j * nvir + b]);
                    }
                }
            }
        }
        e
    };

    let mut ecc = energy(&t1, &t2);

    let mut diis = Diis::new(8, 4);
    let max_iter = 100;
    let conv = 1e-10;

    for _iter in 1..=max_iter {

        // ---- F intermediates ----
        let mut fkc = vec![0.0_f64; nocc * nvir];
        for k in 0..nocc {
            for c in 0..nvir {
                let mut val = 0.0;
                for l in 0..nocc {
                    for d in 0..nvir {
                        val += w(k, l, nocc + c, nocc + d) * t1[l * nvir + d];
                    }
                }
                fkc[k * nvir + c] = val;
            }
        }

        let mut fki = vec![0.0_f64; nocc * nocc];
        for k in 0..nocc {
            for i in 0..nocc {
                let mut val = 0.0;
                for l in 0..nocc {
                    for c in 0..nvir {
                        for d in 0..nvir {
                            let ww = w(k, l, nocc + c, nocc + d);
                            val += ww * (t2[t2i(i, l, c, d)] + t1[i * nvir + c] * t1[l * nvir + d]);
                        }
                    }
                }
                fki[k * nocc + i] = val;
            }
        }

        let mut fac = vec![0.0_f64; nvir * nvir];
        for a in 0..nvir {
            for c in 0..nvir {
                let mut val = 0.0;
                for k in 0..nocc {
                    for l in 0..nocc {
                        for d in 0..nvir {
                            let ww = w(k, l, nocc + c, nocc + d);
                            val -= ww * (t2[t2i(k, l, a, d)] + t1[k * nvir + a] * t1[l * nvir + d]);
                        }
                    }
                }
                fac[a * nvir + c] = val;
            }
        }

        // ---- L intermediates ----
        let mut lki = vec![0.0_f64; nocc * nocc];
        for k in 0..nocc {
            for i in 0..nocc {
                let mut val = fki[k * nocc + i];
                for l in 0..nocc {
                    for c in 0..nvir {
                        val += w(l, k, nocc + c, i) * t1[l * nvir + c];
                    }
                }
                lki[k * nocc + i] = val;
            }
        }

        let mut lac = vec![0.0_f64; nvir * nvir];
        for a in 0..nvir {
            for c in 0..nvir {
                let mut val = fac[a * nvir + c];
                for k in 0..nocc {
                    for d in 0..nvir {
                        val += w(k, nocc + a, nocc + d, nocc + c) * t1[k * nvir + d];
                    }
                }
                lac[a * nvir + c] = val;
            }
        }

        // ---- W^{kl}_{ij} ----
        let mut wklij = vec![0.0_f64; nocc * nocc * nocc * nocc];
        for k in 0..nocc {
            for l in 0..nocc {
                for i in 0..nocc {
                    for j in 0..nocc {
                        let mut val = v(k, l, i, j);
                        for c in 0..nvir {
                            let cc = nocc + c;
                            val += v(l, k, cc, i) * t1[j * nvir + c];
                            val += v(k, l, cc, j) * t1[i * nvir + c];
                        }
                        for c in 0..nvir {
                            for d in 0..nvir {
                                let cc = nocc + c;
                                let dd = nocc + d;
                                val += v(k, l, cc, dd) * (t2[t2i(i, j, c, d)] + t1[i * nvir + c] * t1[j * nvir + d]);
                            }
                        }
                        wklij[((k * nocc + l) * nocc + i) * nocc + j] = val;
                    }
                }
            }
        }

        // ---- W^{ab}_{cd} ----
        let vv = nvir * nvir;
        let mut wabcd = vec![0.0_f64; vv * vv];
        for a in 0..nvir {
            for b in 0..nvir {
                let aa = nocc + a;
                let bb = nocc + b;
                for c in 0..nvir {
                    for d in 0..nvir {
                        let cc = nocc + c;
                        let dd = nocc + d;
                        let mut val = v(aa, bb, cc, dd);
                        for k in 0..nocc {
                            val -= v(k, aa, dd, cc) * t1[k * nvir + b];
                            val -= v(k, bb, cc, dd) * t1[k * nvir + a];
                        }
                        wabcd[(a * nvir + b) * vv + c * nvir + d] = val;
                    }
                }
            }
        }

        // ---- W^{ak}_{ic} and W^{ak}_{ci} ----
        let mut wakic = vec![0.0_f64; nvir * nocc * nocc * nvir];
        for a in 0..nvir {
            for k in 0..nocc {
                for i in 0..nocc {
                    for c in 0..nvir {
                        let aa = nocc + a;
                        let cc = nocc + c;
                        let mut val = v(aa, k, i, cc);
                        for l in 0..nocc {
                            val -= v(k, l, cc, i) * t1[l * nvir + a];
                        }
                        for d in 0..nvir {
                            let dd = nocc + d;
                            val += v(k, aa, cc, dd) * t1[i * nvir + d];
                        }
                        for l in 0..nocc {
                            for d in 0..nvir {
                                let dd = nocc + d;
                                let vlk = v(l, k, dd, cc);
                                val -= 0.5 * vlk * t2[t2i(i, l, d, a)];
                                val -= vlk * t1[i * nvir + d] * t1[l * nvir + a];
                                val += 0.5 * w(l, k, dd, cc) * t2[t2i(i, l, a, d)];
                            }
                        }
                        wakic[((a * nocc + k) * nocc + i) * nvir + c] = val;
                    }
                }
            }
        }

        let mut wakci = vec![0.0_f64; nvir * nocc * nvir * nocc];
        for a in 0..nvir {
            for k in 0..nocc {
                for c in 0..nvir {
                    for i in 0..nocc {
                        let aa = nocc + a;
                        let cc = nocc + c;
                        let mut val = v(aa, k, cc, i);
                        for l in 0..nocc {
                            val -= v(l, k, cc, i) * t1[l * nvir + a];
                        }
                        for d in 0..nvir {
                            let dd = nocc + d;
                            val += v(k, aa, dd, cc) * t1[i * nvir + d];
                        }
                        for l in 0..nocc {
                            for d in 0..nvir {
                                let dd = nocc + d;
                                let vlk = v(l, k, cc, dd);
                                val -= 0.5 * vlk * t2[t2i(i, l, d, a)];
                                val -= vlk * t1[i * nvir + d] * t1[l * nvir + a];
                            }
                        }
                        wakci[((a * nocc + k) * nvir + c) * nocc + i] = val;
                    }
                }
            }
        }

        // ---- T1 update ----
        let mut new_t1 = vec![0.0_f64; t1_size];
        for i in 0..nocc {
            for a in 0..nvir {
                let aa = nocc + a;
                let mut val = 0.0;
                for c in 0..nvir {
                    val += fac[a * nvir + c] * t1[i * nvir + c];
                }
                for k in 0..nocc {
                    val -= fki[k * nocc + i] * t1[k * nvir + a];
                }
                for k in 0..nocc {
                    for c in 0..nvir {
                        let fc = fkc[k * nvir + c];
                        val += fc * (2.0 * t2[t2i(k, i, c, a)] - t2[t2i(i, k, c, a)] + t1[i * nvir + c] * t1[k * nvir + a]);
                    }
                }
                for k in 0..nocc {
                    for c in 0..nvir {
                        val += w(aa, k, i, nocc + c) * t1[k * nvir + c];
                    }
                }
                for k in 0..nocc {
                    for c in 0..nvir {
                        for d in 0..nvir {
                            val += w(aa, k, nocc + c, nocc + d) * (t2[t2i(i, k, c, d)] + t1[i * nvir + c] * t1[k * nvir + d]);
                        }
                    }
                }
                for k in 0..nocc {
                    for l in 0..nocc {
                        for c in 0..nvir {
                            val -= w(k, l, i, nocc + c) * (t2[t2i(k, l, a, c)] + t1[k * nvir + a] * t1[l * nvir + c]);
                        }
                    }
                }
                new_t1[i * nvir + a] = val / dia[i * nvir + a];
            }
        }

        // ---- T2 update ----
        let mut tau = vec![0.0_f64; t2_size];
        for i in 0..nocc {
            for j in 0..nocc {
                for a in 0..nvir {
                    for b in 0..nvir {
                        let idx = t2i(i, j, a, b);
                        tau[idx] = t2[idx] + t1[i * nvir + a] * t1[j * nvir + b];
                    }
                }
            }
        }

        let oo = nocc * nocc;
        let mut raw = vec![0.0_f64; t2_size];

        // 0.5 * Σ_{cd} tau[ij,cd] * Wabcd[ab,cd]
        for ij in 0..oo {
            let ij_off = ij * vv;
            for ab in 0..vv {
                let mut sum = 0.0;
                let ab_off = ab * vv;
                for cd in 0..vv {
                    sum += tau[ij_off + cd] * wabcd[ab_off + cd];
                }
                raw[ij_off + ab] += 0.5 * sum;
            }
        }

        // 0.5 * Σ_{kl} Wklij[kl,ij] * tau[kl,ab]
        for ij in 0..oo {
            for ab in 0..vv {
                let mut sum = 0.0;
                for kl in 0..oo {
                    sum += wklij[kl * oo + ij] * tau[kl * vv + ab];
                }
                raw[ij * vv + ab] += 0.5 * sum;
            }
        }

        // Add remaining terms
        for i in 0..nocc {
            for a in 0..nvir {
                let aa = nocc + a;
                for j in 0..nocc {
                    for b in 0..nvir {
                        let bb = nocc + b;
                        let mut val = 0.5 * v(i, j, aa, bb);

                        for c in 0..nvir {
                            val += lac[a * nvir + c] * t2[t2i(i, j, c, b)];
                        }
                        for k in 0..nocc {
                            val -= lki[k * nocc + i] * t2[t2i(k, j, a, b)];
                        }
                        for c in 0..nvir {
                            val += v(aa, bb, i, nocc + c) * t1[j * nvir + c];
                        }
                        for k in 0..nocc {
                            for c in 0..nvir {
                                val -= v(k, bb, i, nocc + c) * t1[k * nvir + a] * t1[j * nvir + c];
                            }
                        }
                        for k in 0..nocc {
                            val -= v(aa, k, i, j) * t1[k * nvir + b];
                        }
                        for k in 0..nocc {
                            for c in 0..nvir {
                                val -= v(aa, k, i, nocc + c) * t1[j * nvir + c] * t1[k * nvir + b];
                            }
                        }
                        for k in 0..nocc {
                            for c in 0..nvir {
                                let w1 = wakic[((a * nocc + k) * nocc + i) * nvir + c];
                                let w2 = wakci[((a * nocc + k) * nvir + c) * nocc + i];
                                let w3 = wakci[((b * nocc + k) * nvir + c) * nocc + i];
                                val += 2.0 * w1 * t2[t2i(k, j, c, b)];
                                val -= w2 * t2[t2i(k, j, c, b)];
                                val -= w1 * t2[t2i(k, j, b, c)];
                                val -= w3 * t2[t2i(k, j, a, c)];
                            }
                        }

                        raw[t2i(i, j, a, b)] += val;
                    }
                }
            }
        }

        // Symmetrize: t2_new(i,j,a,b) = [raw(i,j,a,b) + raw(j,i,b,a)] / D
        let mut new_t2 = vec![0.0_f64; t2_size];
        for i in 0..nocc {
            for j in 0..nocc {
                for a in 0..nvir {
                    for b in 0..nvir {
                        let idx = t2i(i, j, a, b);
                        new_t2[idx] = (raw[idx] + raw[t2i(j, i, b, a)]) / dijab[idx];
                    }
                }
            }
        }

        // DIIS
        let total_size = t1_size + t2_size;
        let mut amp_vec = vec![0.0_f64; total_size];
        let mut err_vec = vec![0.0_f64; total_size];
        for k in 0..t1_size {
            amp_vec[k] = new_t1[k];
            err_vec[k] = new_t1[k] - t1[k];
        }
        for k in 0..t2_size {
            amp_vec[t1_size + k] = new_t2[k];
            err_vec[t1_size + k] = new_t2[k] - t2[k];
        }
        diis.push(&amp_vec, &err_vec);
        if diis.can_extrapolate() {
            let extrap = diis.extrapolate();
            for k in 0..t1_size {
                new_t1[k] = extrap[k];
            }
            for k in 0..t2_size {
                new_t2[k] = extrap[t1_size + k];
            }
        }

        t1 = new_t1;
        t2 = new_t2;

        let new_ecc = energy(&t1, &t2);
        let delta_e = new_ecc - ecc;
        ecc = new_ecc;

        if delta_e.abs() < conv {
            return ecc;
        }
    }

    ecc
}

// ---------------------------------------------------------------------------
// Spin-orbital CCSD loop (for UHF)
// ---------------------------------------------------------------------------

fn ccsd_loop<F: Fn(usize, usize, usize, usize) -> f64>(
    get_int: &F,
    fock: &[f64],
    nocc: usize,
    nvir: usize,
    nso: usize,
) -> f64 {

    let mut dia = vec![0.0_f64; nocc * nvir];
    for i in 0..nocc {
        for a in 0..nvir {
            dia[i * nvir + a] = fock[i * nso + i] - fock[(nocc + a) * nso + (nocc + a)];
        }
    }

    let t2_size = nocc * nocc * nvir * nvir;
    let mut dijab = vec![0.0_f64; t2_size];
    for i in 0..nocc {
        for j in 0..nocc {
            for a in 0..nvir {
                for b in 0..nvir {
                    dijab[t2idx(i, j, a, b, nocc, nvir)] =
                        fock[i * nso + i] + fock[j * nso + j]
                        - fock[(nocc + a) * nso + (nocc + a)]
                        - fock[(nocc + b) * nso + (nocc + b)];
                }
            }
        }
    }

    let t1_size = nocc * nvir;
    let mut t1 = vec![0.0_f64; t1_size];
    let mut t2 = vec![0.0_f64; t2_size];
    for i in 0..nocc {
        for j in 0..nocc {
            for a in 0..nvir {
                for b in 0..nvir {
                    let idx = t2idx(i, j, a, b, nocc, nvir);
                    t2[idx] = get_int(i, j, nocc + a, nocc + b) / dijab[idx];
                }
            }
        }
    }

    let energy = |t1: &[f64], t2: &[f64]| -> f64 {
        let mut e = 0.0;
        for i in 0..nocc {
            for a in 0..nvir {
                let aa = nocc + a;
                e += fock[i * nso + aa] * t1[i * nvir + a];
                for j in 0..nocc {
                    for b in 0..nvir {
                        let ijab = get_int(i, j, aa, nocc + b);
                        e += 0.25 * ijab * t2[t2idx(i, j, a, b, nocc, nvir)];
                        e += 0.5 * ijab * t1[i * nvir + a] * t1[j * nvir + b];
                    }
                }
            }
        }
        e
    };

    let mut ecc = energy(&t1, &t2);

    let mut diis = Diis::new(8, 4);
    let max_iter = 100;
    let conv = 1e-10;

    for _iter in 1..=max_iter {
        let mut tau_tilde = vec![0.0_f64; t2_size];
        let mut tau = vec![0.0_f64; t2_size];
        for i in 0..nocc {
            for j in 0..nocc {
                for a in 0..nvir {
                    for b in 0..nvir {
                        let idx = t2idx(i, j, a, b, nocc, nvir);
                        let ia_jb = t1[i * nvir + a] * t1[j * nvir + b];
                        let ib_ja = t1[i * nvir + b] * t1[j * nvir + a];
                        tau_tilde[idx] = t2[idx] + 0.5 * (ia_jb - ib_ja);
                        tau[idx] = t2[idx] + ia_jb - ib_ja;
                    }
                }
            }
        }

        // Stanton Eq.3: F_ae
        let mut f_ae = vec![0.0_f64; nvir * nvir];
        for a in 0..nvir {
            for e in 0..nvir {
                let aa = nocc + a;
                let ee = nocc + e;
                let mut val = if a != e { fock[aa * nso + ee] } else { 0.0 };
                for m in 0..nocc { val -= 0.5 * fock[m * nso + ee] * t1[m * nvir + a]; }
                for m in 0..nocc {
                    for f in 0..nvir { val += t1[m * nvir + f] * get_int(m, aa, nocc + f, ee); }
                }
                for m in 0..nocc {
                    for nn in 0..nocc {
                        for f in 0..nvir {
                            val -= 0.5 * tau_tilde[t2idx(m, nn, a, f, nocc, nvir)] * get_int(m, nn, ee, nocc + f);
                        }
                    }
                }
                f_ae[a * nvir + e] = val;
            }
        }

        // Stanton Eq.4: F_mi
        let mut f_mi = vec![0.0_f64; nocc * nocc];
        for m in 0..nocc {
            for i in 0..nocc {
                let mut val = if m != i { fock[m * nso + i] } else { 0.0 };
                for e in 0..nvir { val += 0.5 * fock[m * nso + (nocc + e)] * t1[i * nvir + e]; }
                for nn in 0..nocc {
                    for e in 0..nvir { val += t1[nn * nvir + e] * get_int(m, nn, i, nocc + e); }
                }
                for nn in 0..nocc {
                    for e in 0..nvir {
                        for f in 0..nvir {
                            val += 0.5 * tau_tilde[t2idx(i, nn, e, f, nocc, nvir)] * get_int(m, nn, nocc + e, nocc + f);
                        }
                    }
                }
                f_mi[m * nocc + i] = val;
            }
        }

        // Stanton Eq.5: F_me
        let mut f_me = vec![0.0_f64; nocc * nvir];
        for m in 0..nocc {
            for e in 0..nvir {
                let mut val = fock[m * nso + (nocc + e)];
                for nn in 0..nocc {
                    for f in 0..nvir { val += t1[nn * nvir + f] * get_int(m, nn, nocc + e, nocc + f); }
                }
                f_me[m * nvir + e] = val;
            }
        }

        // Stanton Eq.6: W_mnij
        let mut w_mnij = vec![0.0_f64; nocc * nocc * nocc * nocc];
        for m in 0..nocc {
            for nn in 0..nocc {
                for i in 0..nocc {
                    for j in 0..nocc {
                        let mut val = get_int(m, nn, i, j);
                        for e in 0..nvir {
                            val += t1[j * nvir + e] * get_int(m, nn, i, nocc + e);
                            val -= t1[i * nvir + e] * get_int(m, nn, j, nocc + e);
                        }
                        for e in 0..nvir {
                            for f in 0..nvir {
                                val += 0.25 * tau[t2idx(i, j, e, f, nocc, nvir)] * get_int(m, nn, nocc + e, nocc + f);
                            }
                        }
                        w_mnij[((m * nocc + nn) * nocc + i) * nocc + j] = val;
                    }
                }
            }
        }

        // Stanton Eq.7: W_abef
        let mut w_abef = vec![0.0_f64; nvir * nvir * nvir * nvir];
        for a in 0..nvir {
            for b in 0..nvir {
                for e in 0..nvir {
                    for f in 0..nvir {
                        let aa = nocc + a; let bb = nocc + b; let ee = nocc + e; let ff = nocc + f;
                        let mut val = get_int(aa, bb, ee, ff);
                        for m in 0..nocc {
                            val -= t1[m * nvir + b] * get_int(aa, m, ee, ff);
                            val += t1[m * nvir + a] * get_int(bb, m, ee, ff);
                        }
                        for m in 0..nocc {
                            for nn in 0..nocc {
                                val += 0.25 * tau[t2idx(m, nn, a, b, nocc, nvir)] * get_int(m, nn, ee, ff);
                            }
                        }
                        w_abef[((a * nvir + b) * nvir + e) * nvir + f] = val;
                    }
                }
            }
        }

        // Stanton Eq.8: W_mbej
        let mut w_mbej = vec![0.0_f64; nocc * nvir * nvir * nocc];
        for m in 0..nocc {
            for b in 0..nvir {
                for e in 0..nvir {
                    for j in 0..nocc {
                        let bb = nocc + b; let ee = nocc + e;
                        let mut val = get_int(m, bb, ee, j);
                        for f in 0..nvir { val += t1[j * nvir + f] * get_int(m, bb, ee, nocc + f); }
                        for nn in 0..nocc { val -= t1[nn * nvir + b] * get_int(m, nn, ee, j); }
                        for nn in 0..nocc {
                            for f in 0..nvir {
                                val -= (0.5 * t2[t2idx(j, nn, f, b, nocc, nvir)] + t1[j * nvir + f] * t1[nn * nvir + b])
                                    * get_int(m, nn, ee, nocc + f);
                            }
                        }
                        w_mbej[((m * nvir + b) * nvir + e) * nocc + j] = val;
                    }
                }
            }
        }

        // Stanton Eq.1: T1 update
        let mut new_t1 = vec![0.0_f64; t1_size];
        for i in 0..nocc {
            for a in 0..nvir {
                let aa = nocc + a;
                let mut val = fock[i * nso + aa];
                for e in 0..nvir { val += t1[i * nvir + e] * f_ae[a * nvir + e]; }
                for m in 0..nocc { val -= t1[m * nvir + a] * f_mi[m * nocc + i]; }
                for m in 0..nocc {
                    for e in 0..nvir { val += t2[t2idx(i, m, a, e, nocc, nvir)] * f_me[m * nvir + e]; }
                }
                for nn in 0..nocc {
                    for f in 0..nvir { val -= t1[nn * nvir + f] * get_int(nn, aa, i, nocc + f); }
                }
                for m in 0..nocc {
                    for e in 0..nvir {
                        for f in 0..nvir {
                            val -= 0.5 * t2[t2idx(i, m, e, f, nocc, nvir)] * get_int(m, aa, nocc + e, nocc + f);
                        }
                    }
                }
                for m in 0..nocc {
                    for nn in 0..nocc {
                        for e in 0..nvir {
                            val -= 0.5 * t2[t2idx(m, nn, a, e, nocc, nvir)] * get_int(nn, m, nocc + e, i);
                        }
                    }
                }
                new_t1[i * nvir + a] = val / dia[i * nvir + a];
            }
        }

        // Stanton Eq.2: T2 update
        let mut new_t2 = vec![0.0_f64; t2_size];
        for i in 0..nocc {
            for j in 0..nocc {
                for a in 0..nvir {
                    for b in 0..nvir {
                        let aa = nocc + a; let bb = nocc + b;
                        let idx = t2idx(i, j, a, b, nocc, nvir);
                        let mut val = get_int(i, j, aa, bb);

                        for e in 0..nvir {
                            let mut fb = f_ae[b * nvir + e];
                            for m in 0..nocc { fb -= 0.5 * t1[m * nvir + b] * f_me[m * nvir + e]; }
                            val += t2[t2idx(i, j, a, e, nocc, nvir)] * fb;
                            let mut fa = f_ae[a * nvir + e];
                            for m in 0..nocc { fa -= 0.5 * t1[m * nvir + a] * f_me[m * nvir + e]; }
                            val -= t2[t2idx(i, j, b, e, nocc, nvir)] * fa;
                        }

                        for m in 0..nocc {
                            let mut fj = f_mi[m * nocc + j];
                            for e in 0..nvir { fj += 0.5 * t1[j * nvir + e] * f_me[m * nvir + e]; }
                            val -= t2[t2idx(i, m, a, b, nocc, nvir)] * fj;
                            let mut fi = f_mi[m * nocc + i];
                            for e in 0..nvir { fi += 0.5 * t1[i * nvir + e] * f_me[m * nvir + e]; }
                            val += t2[t2idx(j, m, a, b, nocc, nvir)] * fi;
                        }

                        for m in 0..nocc {
                            for nn in 0..nocc {
                                val += 0.5 * tau[t2idx(m, nn, a, b, nocc, nvir)]
                                    * w_mnij[((m * nocc + nn) * nocc + i) * nocc + j];
                            }
                        }

                        for e in 0..nvir {
                            for f in 0..nvir {
                                val += 0.5 * tau[t2idx(i, j, e, f, nocc, nvir)]
                                    * w_abef[((a * nvir + b) * nvir + e) * nvir + f];
                            }
                        }

                        for m in 0..nocc {
                            for e in 0..nvir {
                                let ee = nocc + e;
                                val += t2[t2idx(i, m, a, e, nocc, nvir)] * w_mbej[((m * nvir + b) * nvir + e) * nocc + j];
                                val -= t1[i * nvir + e] * t1[m * nvir + a] * get_int(m, bb, ee, j);
                                val -= t2[t2idx(j, m, a, e, nocc, nvir)] * w_mbej[((m * nvir + b) * nvir + e) * nocc + i];
                                val += t1[j * nvir + e] * t1[m * nvir + a] * get_int(m, bb, ee, i);
                                val -= t2[t2idx(i, m, b, e, nocc, nvir)] * w_mbej[((m * nvir + a) * nvir + e) * nocc + j];
                                val += t1[i * nvir + e] * t1[m * nvir + b] * get_int(m, aa, ee, j);
                                val += t2[t2idx(j, m, b, e, nocc, nvir)] * w_mbej[((m * nvir + a) * nvir + e) * nocc + i];
                                val -= t1[j * nvir + e] * t1[m * nvir + b] * get_int(m, aa, ee, i);
                            }
                        }

                        for e in 0..nvir {
                            val += t1[i * nvir + e] * get_int(aa, bb, nocc + e, j);
                            val -= t1[j * nvir + e] * get_int(aa, bb, nocc + e, i);
                        }

                        for m in 0..nocc {
                            val -= t1[m * nvir + a] * get_int(m, bb, i, j);
                            val += t1[m * nvir + b] * get_int(m, aa, i, j);
                        }

                        new_t2[idx] = val / dijab[idx];
                    }
                }
            }
        }

        // DIIS
        let total_size = t1_size + t2_size;
        let mut amp_vec = vec![0.0_f64; total_size];
        let mut err_vec = vec![0.0_f64; total_size];
        for k in 0..t1_size {
            amp_vec[k] = new_t1[k];
            err_vec[k] = new_t1[k] - t1[k];
        }
        for k in 0..t2_size {
            amp_vec[t1_size + k] = new_t2[k];
            err_vec[t1_size + k] = new_t2[k] - t2[k];
        }
        diis.push(&amp_vec, &err_vec);
        if diis.can_extrapolate() {
            let extrap = diis.extrapolate();
            for k in 0..t1_size { new_t1[k] = extrap[k]; }
            for k in 0..t2_size { new_t2[k] = extrap[t1_size + k]; }
        }

        t1 = new_t1;
        t2 = new_t2;

        let new_ecc = energy(&t1, &t2);
        let delta_e = new_ecc - ecc;
        ecc = new_ecc;

        if delta_e.abs() < conv {
            return ecc;
        }
    }

    ecc
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

/// Compute CCSD correlation energy from converged RHF results.
/// Uses spatial-orbital formulation (~16x faster than spin-orbital).
pub fn compute_ccsd(
    eri: &[f64],
    coeff: &[f64],
    epsilon: &[f64],
    nocc: usize,
    n: usize,
) -> f64 {
    let mo_ints = full_mo_transform(coeff, eri, n);
    ccsd_spatial_loop(&mo_ints, epsilon, nocc, n)
}

/// Compute CCSD correlation energy from converged UHF results.
/// Uses blocked spin-orbital ordering.
pub fn compute_uccsd(
    eri: &[f64],
    ca: &[f64],
    cb: &[f64],
    epsilon_a: &[f64],
    epsilon_b: &[f64],
    nocc_a: usize,
    nocc_b: usize,
    n: usize,
) -> f64 {
    let nvir_a = n - nocc_a;
    let nvir_b = n - nocc_b;
    let nocc_s = nocc_a + nocc_b;
    let nvir_s = nvir_a + nvir_b;
    let nso = nocc_s + nvir_s;

    let n2 = n * n;
    let n3 = n2 * n;

    let mo_aa = full_mo_transform_4c(ca, ca, ca, ca, eri, n);
    let mo_bb = full_mo_transform_4c(cb, cb, cb, cb, eri, n);
    let mo_ab = full_mo_transform_4c(ca, ca, cb, cb, eri, n);

    let spin_of = |p: usize| -> usize {
        if p < nocc_a { 0 }
        else if p < nocc_s { 1 }
        else if p < nocc_s + nvir_a { 0 }
        else { 1 }
    };
    let spat_of = |p: usize| -> usize {
        if p < nocc_a { p }
        else if p < nocc_s { p - nocc_a }
        else if p < nocc_s + nvir_a { nocc_a + (p - nocc_s) }
        else { nocc_b + (p - nocc_s - nvir_a) }
    };

    let get_int = |pp: usize, qq: usize, rr: usize, ss: usize| -> f64 {
        let sp = spin_of(pp); let sq = spin_of(qq);
        let sr = spin_of(rr); let sss = spin_of(ss);
        let p = spat_of(pp); let q = spat_of(qq);
        let r = spat_of(rr); let s = spat_of(ss);
        let mut val = 0.0;
        if sp == sr && sq == sss {
            if sp == 0 && sq == 0      { val += mo_aa[p * n3 + r * n2 + q * n + s]; }
            else if sp == 1 && sq == 1 { val += mo_bb[p * n3 + r * n2 + q * n + s]; }
            else if sp == 0 && sq == 1 { val += mo_ab[p * n3 + r * n2 + q * n + s]; }
            else                       { val += mo_ab[q * n3 + s * n2 + p * n + r]; }
        }
        if sp == sss && sq == sr {
            if sp == 0 && sq == 0      { val -= mo_aa[p * n3 + s * n2 + q * n + r]; }
            else if sp == 1 && sq == 1 { val -= mo_bb[p * n3 + s * n2 + q * n + r]; }
            else if sp == 0 && sq == 1 { val -= mo_ab[p * n3 + s * n2 + q * n + r]; }
            else                       { val -= mo_ab[q * n3 + r * n2 + p * n + s]; }
        }
        val
    };

    let mut fock = vec![0.0_f64; nso * nso];
    for i in 0..nocc_a { fock[i * nso + i] = epsilon_a[i]; }
    for i in 0..nocc_b { fock[(nocc_a + i) * nso + (nocc_a + i)] = epsilon_b[i]; }
    for a in 0..nvir_a { fock[(nocc_s + a) * nso + (nocc_s + a)] = epsilon_a[nocc_a + a]; }
    for a in 0..nvir_b { fock[(nocc_s + nvir_a + a) * nso + (nocc_s + nvir_a + a)] = epsilon_b[nocc_b + a]; }

    ccsd_loop(&get_int, &fock, nocc_s, nvir_s, nso)
}
