/// WASM-accelerated MP3 correlation energy — port of mp3.ts (RHF spatial-orbital)

use crate::integrals::eri_1d_index;

/// Full 4-index MO integral transformation with 4 independent coefficient matrices.
/// c1..c4: N*N row-major (AO x MO), eri: 8-fold symmetric.
/// Step k uses ck. Returns (pq|rs) chemist notation, N^4 elements.
pub(crate) fn full_mo_transform_4c(
    c1: &[f64], c2: &[f64], c3: &[f64], c4: &[f64],
    eri: &[f64], n: usize,
) -> Vec<f64> {
    let n2 = n * n;
    let n3 = n2 * n;
    let n4 = n3 * n;

    // Step 1: (pν|λσ) = Σ_μ C1[μ,p] * (μν|λσ)
    let mut step1 = vec![0.0_f64; n4];
    for p in 0..n {
        for nu in 0..n {
            for lam in 0..n {
                for sig in 0..n {
                    let mut val = 0.0;
                    for mu in 0..n {
                        let c = c1[mu * n + p];
                        if c.abs() < 1e-15 {
                            continue;
                        }
                        val += c * eri[eri_1d_index(mu, nu, lam, sig)];
                    }
                    step1[p * n3 + nu * n2 + lam * n + sig] = val;
                }
            }
        }
    }

    // Step 2: (pq|λσ) = Σ_ν C2[ν,q] * (pν|λσ)
    let mut step2 = vec![0.0_f64; n4];
    for p in 0..n {
        for q in 0..n {
            for lam in 0..n {
                for sig in 0..n {
                    let mut val = 0.0;
                    for nu in 0..n {
                        let c = c2[nu * n + q];
                        if c.abs() < 1e-15 {
                            continue;
                        }
                        val += c * step1[p * n3 + nu * n2 + lam * n + sig];
                    }
                    step2[p * n3 + q * n2 + lam * n + sig] = val;
                }
            }
        }
    }
    drop(step1);

    // Step 3: (pq|rσ) = Σ_λ C3[λ,r] * (pq|λσ)
    let mut step3 = vec![0.0_f64; n4];
    for p in 0..n {
        for q in 0..n {
            for r in 0..n {
                for sig in 0..n {
                    let mut val = 0.0;
                    for lam in 0..n {
                        let c = c3[lam * n + r];
                        if c.abs() < 1e-15 {
                            continue;
                        }
                        val += c * step2[p * n3 + q * n2 + lam * n + sig];
                    }
                    step3[p * n3 + q * n2 + r * n + sig] = val;
                }
            }
        }
    }
    drop(step2);

    // Step 4: (pq|rs) = Σ_σ C4[σ,s] * (pq|rσ)
    let mut mo_ints = vec![0.0_f64; n4];
    for p in 0..n {
        for q in 0..n {
            for r in 0..n {
                for s in 0..n {
                    let mut val = 0.0;
                    for sig in 0..n {
                        let c = c4[sig * n + s];
                        if c.abs() < 1e-15 {
                            continue;
                        }
                        val += c * step3[p * n3 + q * n2 + r * n + sig];
                    }
                    mo_ints[p * n3 + q * n2 + r * n + s] = val;
                }
            }
        }
    }

    mo_ints
}

/// Full 4-index MO integral transformation (RHF: same C for all 4 indices).
pub(crate) fn full_mo_transform(coeff: &[f64], eri: &[f64], n: usize) -> Vec<f64> {
    full_mo_transform_4c(coeff, coeff, coeff, coeff, eri, n)
}

/// Compute MP3 correlation energy from converged RHF results.
///
/// `eri`: 8-fold symmetric ERI data
/// `coeff`: MO coefficient matrix (N*N, row-major)
/// `epsilon`: orbital energies (N elements)
/// `nocc`: number of occupied spatial orbitals
/// `n`: number of basis functions
///
/// Returns [mp2, mp3] as a 2-element Vec.
pub fn compute_mp3(
    eri: &[f64],
    coeff: &[f64],
    epsilon: &[f64],
    nocc: usize,
    n: usize,
) -> Vec<f64> {
    let mo_ints = full_mo_transform(coeff, eri, n);

    let n3 = n * n * n;
    let n2 = n * n;

    // MO integral accessor (chemist notation)
    let mo = |p: usize, q: usize, r: usize, s: usize| -> f64 {
        mo_ints[p * n3 + q * n2 + r * n + s]
    };

    // --- MP2 energy (spatial orbital, spin-traced) ---
    let mut emp2 = 0.0_f64;
    for i in 0..nocc {
        for j in 0..nocc {
            for a in nocc..n {
                for b in nocc..n {
                    let iajb = mo(i, a, j, b);
                    let ibja = mo(i, b, j, a);
                    let d = epsilon[i] + epsilon[j] - epsilon[a] - epsilon[b];
                    emp2 += iajb * (2.0 * iajb - ibja) / d;
                }
            }
        }
    }

    // --- PP (2h4p): particle-particle ladder ---
    let mut epp = 0.0_f64;
    for i in 0..nocc {
        for j in 0..nocc {
            for a in nocc..n {
                for b in nocc..n {
                    let d1 = epsilon[i] + epsilon[j] - epsilon[a] - epsilon[b];
                    let iajb = mo(i, a, j, b);
                    for c in nocc..n {
                        for d in nocc..n {
                            let d2 = epsilon[i] + epsilon[j] - epsilon[c] - epsilon[d];
                            epp += iajb * mo(a, c, b, d)
                                * (2.0 * mo(i, c, j, d) - mo(i, d, j, c))
                                / (d1 * d2);
                        }
                    }
                }
            }
        }
    }

    // --- HH (4h2p): hole-hole ladder ---
    let mut ehh = 0.0_f64;
    for i in 0..nocc {
        for j in 0..nocc {
            for k in 0..nocc {
                for l in 0..nocc {
                    let ikjl = mo(i, k, j, l);
                    for a in nocc..n {
                        for b in nocc..n {
                            let d1 = epsilon[i] + epsilon[j] - epsilon[a] - epsilon[b];
                            let d2 = epsilon[k] + epsilon[l] - epsilon[a] - epsilon[b];
                            ehh += mo(i, a, j, b) * ikjl
                                * (2.0 * mo(k, a, l, b) - mo(k, b, l, a))
                                / (d1 * d2);
                        }
                    }
                }
            }
        }
    }

    // --- PH (3h3p): particle-hole ring ---
    let mut eph = 0.0_f64;
    for i in 0..nocc {
        for j in 0..nocc {
            for k in 0..nocc {
                for a in nocc..n {
                    for b in nocc..n {
                        for c in nocc..n {
                            let d1 = epsilon[i] + epsilon[k] - epsilon[a] - epsilon[c];
                            let d2 = epsilon[k] + epsilon[j] - epsilon[b] - epsilon[c];
                            let denom = d1 * d2;
                            if denom.abs() < 1e-14 {
                                continue;
                            }
                            let iajb = mo(i, a, j, b);
                            let ijab = mo(i, j, a, b);
                            let kcia = mo(k, c, i, a);
                            let kaic = mo(k, a, i, c);
                            let kcjb = mo(k, c, j, b);
                            let kbjc = mo(k, b, j, c);
                            eph += ((2.0 * iajb - ijab)
                                * (2.0 * kcia - kaic)
                                * (2.0 * kcjb - kbjc)
                                - 3.0 * ijab * kaic * kbjc)
                                / denom;
                        }
                    }
                }
            }
        }
    }

    let emp3 = epp + ehh + eph;
    vec![emp2, emp3]
}

// ---------------------------------------------------------------------------
// Same-spin helpers for UHF-MP3 spatial orbital
// ---------------------------------------------------------------------------

/// Same-spin PP ladder: 0.5 Σ (acbd)(cidj)[(iajb)-(ibja)] / (D_ij^ab D_ij^cd)
fn pp_same_spin(mo: &[f64], eps: &[f64], nocc: usize, n: usize) -> f64 {
    let n3 = n * n * n;
    let n2 = n * n;
    let mut e = 0.0_f64;
    for i in 0..nocc {
        for j in 0..nocc {
            for a in nocc..n {
                for b in nocc..n {
                    let d1 = eps[i] + eps[j] - eps[a] - eps[b];
                    let x = mo[i * n3 + a * n2 + j * n + b] - mo[i * n3 + b * n2 + j * n + a];
                    if x.abs() < 1e-15 { continue; }
                    for c in nocc..n {
                        for d in nocc..n {
                            let d2 = eps[i] + eps[j] - eps[c] - eps[d];
                            e += x * mo[a * n3 + c * n2 + b * n + d]
                                   * mo[c * n3 + i * n2 + d * n + j] / (d1 * d2);
                        }
                    }
                }
            }
        }
    }
    0.5 * e
}

/// Same-spin HH ladder: 0.5 Σ (kilj)(akbl)[(iajb)-(ibja)] / (D_ij^ab D_kl^ab)
fn hh_same_spin(mo: &[f64], eps: &[f64], nocc: usize, n: usize) -> f64 {
    let n3 = n * n * n;
    let n2 = n * n;
    let mut e = 0.0_f64;
    for i in 0..nocc {
        for j in 0..nocc {
            for a in nocc..n {
                for b in nocc..n {
                    let d1 = eps[i] + eps[j] - eps[a] - eps[b];
                    let x = mo[i * n3 + a * n2 + j * n + b] - mo[i * n3 + b * n2 + j * n + a];
                    if x.abs() < 1e-15 { continue; }
                    for k in 0..nocc {
                        for l in 0..nocc {
                            let d2 = eps[k] + eps[l] - eps[a] - eps[b];
                            e += x * mo[k * n3 + i * n2 + l * n + j]
                                   * mo[a * n3 + k * n2 + b * n + l] / (d1 * d2);
                        }
                    }
                }
            }
        }
    }
    0.5 * e
}

/// Same-spin PH ring: Σ [(iajb)-(ibja)][(aikc)-(acki)][(cjbk)-(ckbj)] / (D_ij^ab D_kj^cb)
fn ph_same_spin(mo: &[f64], eps: &[f64], nocc: usize, n: usize) -> f64 {
    let n3 = n * n * n;
    let n2 = n * n;
    let mut e = 0.0_f64;
    for i in 0..nocc {
        for j in 0..nocc {
            for a in nocc..n {
                for b in nocc..n {
                    let d1 = eps[i] + eps[j] - eps[a] - eps[b];
                    let x1 = mo[i * n3 + a * n2 + j * n + b] - mo[i * n3 + b * n2 + j * n + a];
                    if x1.abs() < 1e-15 { continue; }
                    for k in 0..nocc {
                        for c in nocc..n {
                            let d2 = eps[k] + eps[j] - eps[c] - eps[b];
                            let x2 = mo[a * n3 + i * n2 + k * n + c] - mo[a * n3 + c * n2 + k * n + i];
                            let x3 = mo[c * n3 + j * n2 + b * n + k] - mo[c * n3 + k * n2 + b * n + j];
                            e += x1 * x2 * x3 / (d1 * d2);
                        }
                    }
                }
            }
        }
    }
    e
}

// ---------------------------------------------------------------------------
// UHF MP3 — spatial orbital with explicit spin channels (αα, ββ, αβ)
// ---------------------------------------------------------------------------

/// Compute UMP3 correlation energy from converged UHF results.
/// Returns [mp2, mp3] as a 2-element Vec.
pub fn compute_ump3(
    eri: &[f64],
    ca: &[f64],
    cb: &[f64],
    epsilon_a: &[f64],
    epsilon_b: &[f64],
    nocc_a: usize,
    nocc_b: usize,
    n: usize,
) -> Vec<f64> {
    let n2 = n * n;
    let n3 = n2 * n;
    let o_a = nocc_a;
    let o_b = nocc_b;
    let ea = epsilon_a;
    let eb = epsilon_b;

    // 3 MO integral blocks
    let mo_aa = full_mo_transform_4c(ca, ca, ca, ca, eri, n);
    let mo_bb = full_mo_transform_4c(cb, cb, cb, cb, eri, n);
    let mo_ab = full_mo_transform_4c(ca, ca, cb, cb, eri, n);

    // Chemist-notation accessors:
    //   aa(p,q,r,s) = (pα qα|rα sα) = mo_aa[p·n³+q·n²+r·n+s]
    //   bb(p,q,r,s) = (pβ qβ|rβ sβ) = mo_bb[p·n³+q·n²+r·n+s]
    //   ab(p,q,r,s) = (pα qα|rβ sβ) = mo_ab[p·n³+q·n²+r·n+s]
    //   ba(p,q,r,s) = (pβ qβ|rα sα) = mo_ab[r·n³+s·n²+p·n+q]
    let aa = |p: usize, q: usize, r: usize, s: usize| mo_aa[p * n3 + q * n2 + r * n + s];
    let bb = |p: usize, q: usize, r: usize, s: usize| mo_bb[p * n3 + q * n2 + r * n + s];
    let ab = |p: usize, q: usize, r: usize, s: usize| mo_ab[p * n3 + q * n2 + r * n + s];
    let ba = |p: usize, q: usize, r: usize, s: usize| mo_ab[r * n3 + s * n2 + p * n + q];

    // ---- MP2 energy (spatial orbital, 3 spin channels) ----
    let mut emp2 = 0.0_f64;
    // αα
    for i in 0..o_a {
        for j in 0..o_a {
            for a in o_a..n {
                for b in o_a..n {
                    let d = ea[i] + ea[j] - ea[a] - ea[b];
                    emp2 += 0.5 * aa(a, i, b, j) * (aa(i, a, j, b) - aa(i, b, j, a)) / d;
                }
            }
        }
    }
    // ββ
    for i in 0..o_b {
        for j in 0..o_b {
            for a in o_b..n {
                for b in o_b..n {
                    let d = eb[i] + eb[j] - eb[a] - eb[b];
                    emp2 += 0.5 * bb(a, i, b, j) * (bb(i, a, j, b) - bb(i, b, j, a)) / d;
                }
            }
        }
    }
    // αβ
    for i in 0..o_b {
        for j in 0..o_a {
            for a in o_b..n {
                for b in o_a..n {
                    let d = eb[i] + ea[j] - eb[a] - ea[b];
                    let v = ba(i, a, j, b);
                    emp2 += v * v / d;
                }
            }
        }
    }

    // ---- PP channel ----
    let mut epp = pp_same_spin(&mo_aa, ea, o_a, n) + pp_same_spin(&mo_bb, eb, o_b, n);
    // αβ: Σ ba(i,a,j,b) * ba(a,c,b,d) * ba(c,i,d,j) / (D_βα D_βα)
    for i in 0..o_b {
        for j in 0..o_a {
            for a in o_b..n {
                for b in o_a..n {
                    let d1 = eb[i] + ea[j] - eb[a] - ea[b];
                    let v1 = ba(i, a, j, b);
                    if v1.abs() < 1e-15 { continue; }
                    for c in o_b..n {
                        for d in o_a..n {
                            let d2 = eb[i] + ea[j] - eb[c] - ea[d];
                            epp += v1 * ba(a, c, b, d) * ba(c, i, d, j) / (d1 * d2);
                        }
                    }
                }
            }
        }
    }

    // ---- HH channel ----
    let mut ehh = hh_same_spin(&mo_aa, ea, o_a, n) + hh_same_spin(&mo_bb, eb, o_b, n);
    // αβ: Σ ba(i,a,j,b) * ba(k,i,l,j) * ba(a,k,b,l) / (D_βα D_βα)
    for i in 0..o_b {
        for j in 0..o_a {
            for a in o_b..n {
                for b in o_a..n {
                    let d1 = eb[i] + ea[j] - eb[a] - ea[b];
                    let v1 = ba(i, a, j, b);
                    if v1.abs() < 1e-15 { continue; }
                    for k in 0..o_b {
                        for l in 0..o_a {
                            let d2 = eb[k] + ea[l] - eb[a] - ea[b];
                            ehh += v1 * ba(k, i, l, j) * ba(a, k, b, l) / (d1 * d2);
                        }
                    }
                }
            }
        }
    }

    // ---- PH channel ----
    let mut eph = ph_same_spin(&mo_aa, ea, o_a, n) + ph_same_spin(&mo_bb, eb, o_b, n);

    // PH group 3: i,j∈occA; k∈occB; a,b∈virA; c∈virB
    for i in 0..o_a {
        for j in 0..o_a {
            for a in o_a..n {
                for b in o_a..n {
                    let d1 = ea[i] + ea[j] - ea[a] - ea[b];
                    let x = aa(i, a, j, b) - aa(i, b, j, a);
                    if x.abs() < 1e-15 { continue; }
                    for k in 0..o_b {
                        for c in o_b..n {
                            let d2 = eb[k] + ea[j] - eb[c] - ea[b];
                            eph -= x * ab(a, i, k, c) * ba(c, k, b, j) / (d1 * d2);
                        }
                    }
                }
            }
        }
    }

    // PH group 4: i∈occB; j,k∈occA; a∈virB; b,c∈virA
    for i in 0..o_b {
        for j in 0..o_a {
            for a in o_b..n {
                for b in o_a..n {
                    let d1 = eb[i] + ea[j] - eb[a] - ea[b];
                    let v1 = ba(i, a, j, b);
                    if v1.abs() < 1e-15 { continue; }
                    for k in 0..o_a {
                        for c in o_a..n {
                            let d2 = ea[k] + ea[j] - ea[c] - ea[b];
                            let x3 = aa(c, j, b, k) - aa(c, k, b, j);
                            eph += v1 * ba(a, i, k, c) * x3 / (d1 * d2);
                        }
                    }
                }
            }
        }
    }

    // PH group 5: i,k∈occA; j∈occB; a,c∈virA; b∈virB
    for i in 0..o_a {
        for j in 0..o_b {
            for a in o_a..n {
                for b in o_b..n {
                    let d1 = ea[i] + eb[j] - ea[a] - eb[b];
                    let v1 = ab(i, a, j, b);
                    if v1.abs() < 1e-15 { continue; }
                    for k in 0..o_a {
                        for c in o_a..n {
                            let d2 = ea[k] + eb[j] - ea[c] - eb[b];
                            let x2 = aa(a, c, k, i) - aa(a, i, k, c);
                            eph += v1 * x2 * ab(c, k, b, j) / (d1 * d2);
                        }
                    }
                }
            }
        }
    }

    // PH group 6: i∈occA; j,k∈occB; a∈virA; b,c∈virB
    for i in 0..o_a {
        for j in 0..o_b {
            for a in o_a..n {
                for b in o_b..n {
                    let d1 = ea[i] + eb[j] - ea[a] - eb[b];
                    let v1 = ab(i, a, j, b);
                    if v1.abs() < 1e-15 { continue; }
                    for k in 0..o_b {
                        for c in o_b..n {
                            let d2 = eb[k] + eb[j] - eb[c] - eb[b];
                            let x3 = bb(c, j, b, k) - bb(c, k, b, j);
                            eph += v1 * ab(a, i, k, c) * x3 / (d1 * d2);
                        }
                    }
                }
            }
        }
    }

    // PH group 7: i,k∈occB; j∈occA; a,c∈virB; b∈virA
    for i in 0..o_b {
        for j in 0..o_a {
            for a in o_b..n {
                for b in o_a..n {
                    let d1 = eb[i] + ea[j] - eb[a] - ea[b];
                    let v1 = ba(i, a, j, b);
                    if v1.abs() < 1e-15 { continue; }
                    for k in 0..o_b {
                        for c in o_b..n {
                            let d2 = eb[k] + ea[j] - eb[c] - ea[b];
                            let x2 = bb(a, c, k, i) - bb(a, i, k, c);
                            eph += v1 * x2 * ba(c, k, b, j) / (d1 * d2);
                        }
                    }
                }
            }
        }
    }

    // PH group 8: i,j∈occB; k∈occA; a,b∈virB; c∈virA
    for i in 0..o_b {
        for j in 0..o_b {
            for a in o_b..n {
                for b in o_b..n {
                    let d1 = eb[i] + eb[j] - eb[a] - eb[b];
                    let x = bb(i, a, j, b) - bb(i, b, j, a);
                    if x.abs() < 1e-15 { continue; }
                    for k in 0..o_a {
                        for c in o_a..n {
                            let d2 = ea[k] + eb[j] - ea[c] - eb[b];
                            eph -= x * ba(a, i, k, c) * ab(c, k, b, j) / (d1 * d2);
                        }
                    }
                }
            }
        }
    }

    // PH group 9: i,k∈occA; j∈occB; a,c∈virB; b∈virA  (single term)
    for i in 0..o_a {
        for j in 0..o_b {
            for a in o_b..n {
                for b in o_a..n {
                    let d1 = ea[i] + eb[j] - eb[a] - ea[b];
                    let v1 = ab(i, b, j, a); // (iα bα|jβ aβ)
                    if v1.abs() < 1e-15 { continue; }
                    for k in 0..o_a {
                        for c in o_b..n {
                            let d2 = ea[k] + eb[j] - eb[c] - ea[b];
                            eph += v1 * ab(k, i, a, c) * ab(b, k, c, j) / (d1 * d2);
                        }
                    }
                }
            }
        }
    }

    // PH group 10: i,k∈occB; j∈occA; a,c∈virA; b∈virB  (single term)
    for i in 0..o_b {
        for j in 0..o_a {
            for a in o_a..n {
                for b in o_b..n {
                    let d1 = eb[i] + ea[j] - ea[a] - eb[b];
                    let v1 = ba(i, b, j, a); // (iβ bβ|jα aα)
                    if v1.abs() < 1e-15 { continue; }
                    for k in 0..o_b {
                        for c in o_a..n {
                            let d2 = eb[k] + ea[j] - ea[c] - eb[b];
                            eph += v1 * ab(a, c, k, i) * ab(c, j, b, k) / (d1 * d2);
                        }
                    }
                }
            }
        }
    }

    let emp3 = epp + ehh + eph;
    vec![emp2, emp3]
}
