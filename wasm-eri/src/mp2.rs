/// WASM-accelerated MP2 correlation energy — port of mp2.ts

use crate::integrals::eri_1d_index;
use crate::mp3::full_mo_transform_4c;
use crate::simd_utils;

/// Compute MP2 correlation energy from converged RHF results.
///
/// `eri`: 8-fold symmetric ERI data (flat f64)
/// `coeff`: MO coefficient matrix (N*N, row-major, AO x MO)
/// `epsilon`: orbital energies (N elements, sorted)
/// `nocc`: number of occupied spatial orbitals
/// `n`: number of basis functions
///
/// Algorithm: two half-transforms + contraction, O(N^5) total.
pub fn compute_mp2(
    eri: &[f64],
    coeff: &[f64],
    epsilon: &[f64],
    nocc: usize,
    n: usize,
) -> f64 {
    let nvir = n - nocc;

    // Pre-extract virtual orbital coefficient columns (contiguous for SIMD)
    // c_vir[a][sig] = coeff[sig * n + (nocc + a)]
    let mut c_vir: Vec<Vec<f64>> = Vec::with_capacity(nvir);
    for a in 0..nvir {
        let mut col = vec![0.0_f64; n];
        for sig in 0..n {
            col[sig] = coeff[sig * n + (nocc + a)];
        }
        c_vir.push(col);
    }

    // Pre-extract occupied orbital coefficient columns
    // c_occ[j][lam] = coeff[lam * n + j]
    let mut c_occ: Vec<Vec<f64>> = Vec::with_capacity(nocc);
    for j in 0..nocc {
        let mut col = vec![0.0_f64; n];
        for lam in 0..n {
            col[lam] = coeff[lam * n + j];
        }
        c_occ.push(col);
    }

    // Step 1: Build half-transformed integrals for each occupied orbital i.
    let mut half_transformed: Vec<Vec<f64>> = Vec::with_capacity(nocc);
    let mut buf = vec![0.0_f64; n * nvir];

    for i in 0..nocc {
        let mut tmp = vec![0.0_f64; nvir * n * n];

        for lam in 0..n {
            for sig in 0..n {
                for x in buf.iter_mut() {
                    *x = 0.0;
                }
                for mu in 0..n {
                    for nu in 0..n {
                        let eri_val = eri[eri_1d_index(mu, nu, lam, sig)];
                        if eri_val.abs() < 1e-15 {
                            continue;
                        }
                        // buf[mu, a] += C_vir[a][nu] * eri_val -- use daxpy
                        for a in 0..nvir {
                            buf[mu * nvir + a] += c_vir[a][nu] * eri_val;
                        }
                    }
                }

                for mu in 0..n {
                    let c_mi = c_occ[i][mu];
                    if c_mi.abs() < 1e-15 {
                        continue;
                    }
                    for a in 0..nvir {
                        tmp[a * n * n + lam * n + sig] += c_mi * buf[mu * nvir + a];
                    }
                }
            }
        }

        half_transformed.push(tmp);
    }

    // Step 2: Contract using SIMD dot products.
    // (ia|jb) = Σ_λ C_occ[j][λ] * Σ_σ C_vir[b][σ] * tmp[a, λ, σ]
    //         = Σ_λ C_occ[j][λ] * dot(C_vir[b], tmp[a, λ, :])
    let mut emp2 = 0.0_f64;

    for i in 0..nocc {
        let tmp_i = &half_transformed[i];
        for j in i..nocc {
            for a in 0..nvir {
                for b in 0..nvir {
                    // (ia|jb) = Σ_λ C_occ[j][λ] * dot(C_vir[b], tmp_i[a*n*n + λ*n .. +n])
                    let mut iajb = 0.0_f64;
                    let base_a = a * n * n;
                    for lam in 0..n {
                        let c_lj = c_occ[j][lam];
                        if c_lj.abs() < 1e-15 { continue; }
                        let row_start = base_a + lam * n;
                        iajb += c_lj * simd_utils::dot(&c_vir[b], &tmp_i[row_start..row_start + n]);
                    }

                    // (ib|ja) = Σ_λ C_occ[j][λ] * dot(C_vir[a], tmp_i[b*n*n + λ*n .. +n])
                    let mut ibja = 0.0_f64;
                    let base_b = b * n * n;
                    for lam in 0..n {
                        let c_lj = c_occ[j][lam];
                        if c_lj.abs() < 1e-15 { continue; }
                        let row_start = base_b + lam * n;
                        ibja += c_lj * simd_utils::dot(&c_vir[a], &tmp_i[row_start..row_start + n]);
                    }

                    let aa = nocc + a;
                    let bb = nocc + b;
                    let denom = epsilon[i] + epsilon[j] - epsilon[aa] - epsilon[bb];
                    let factor = if i == j { 1.0 } else { 2.0 };
                    emp2 += factor * iajb * (2.0 * iajb - ibja) / denom;
                }
            }
        }
    }

    emp2
}

/// Compute UMP2 correlation energy from converged UHF results.
///
/// `eri`: 8-fold symmetric ERI data
/// `ca`: alpha MO coefficient matrix (N*N, row-major)
/// `cb`: beta MO coefficient matrix (N*N, row-major)
/// `epsilon_a`: alpha orbital energies (N elements)
/// `epsilon_b`: beta orbital energies (N elements)
/// `nocc_a`: number of occupied alpha orbitals
/// `nocc_b`: number of occupied beta orbitals
/// `n`: number of basis functions
pub fn compute_ump2(
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
    let n2 = n * n;
    let n3 = n2 * n;

    // Full MO integral transforms for the 3 spin blocks
    let mo_aa = full_mo_transform_4c(ca, ca, ca, ca, eri, n);
    let mo_bb = full_mo_transform_4c(cb, cb, cb, cb, eri, n);
    let mo_ab = full_mo_transform_4c(ca, ca, cb, cb, eri, n);

    let mut emp2 = 0.0_f64;

    // αα contribution: Σ_{i<j}^α Σ_{a<b}^α |(ia|jb)-(ib|ja)|² / D
    for i in 0..nocc_a {
        for j in (i + 1)..nocc_a {
            for a in 0..nvir_a {
                for b in (a + 1)..nvir_a {
                    let aa = nocc_a + a;
                    let bb = nocc_a + b;
                    let iajb = mo_aa[i * n3 + aa * n2 + j * n + bb];
                    let ibja = mo_aa[i * n3 + bb * n2 + j * n + aa];
                    let antisym = iajb - ibja;
                    let denom = epsilon_a[i] + epsilon_a[j] - epsilon_a[aa] - epsilon_a[bb];
                    emp2 += antisym * antisym / denom;
                }
            }
        }
    }

    // ββ contribution: Σ_{i<j}^β Σ_{a<b}^β |(ia|jb)-(ib|ja)|² / D
    for i in 0..nocc_b {
        for j in (i + 1)..nocc_b {
            for a in 0..nvir_b {
                for b in (a + 1)..nvir_b {
                    let aa = nocc_b + a;
                    let bb = nocc_b + b;
                    let iajb = mo_bb[i * n3 + aa * n2 + j * n + bb];
                    let ibja = mo_bb[i * n3 + bb * n2 + j * n + aa];
                    let antisym = iajb - ibja;
                    let denom = epsilon_b[i] + epsilon_b[j] - epsilon_b[aa] - epsilon_b[bb];
                    emp2 += antisym * antisym / denom;
                }
            }
        }
    }

    // αβ contribution: Σ_{iα,jβ,aα,bβ} |(iα aα|jβ bβ)|² / D
    for i in 0..nocc_a {
        for j in 0..nocc_b {
            for a in 0..nvir_a {
                for b in 0..nvir_b {
                    let aa = nocc_a + a;
                    let bb = nocc_b + b;
                    let iajb = mo_ab[i * n3 + aa * n2 + j * n + bb];
                    let denom = epsilon_a[i] + epsilon_b[j] - epsilon_a[aa] - epsilon_b[bb];
                    emp2 += iajb * iajb / denom;
                }
            }
        }
    }

    emp2
}
