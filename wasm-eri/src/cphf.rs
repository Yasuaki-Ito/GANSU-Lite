/// CPHF: MO-ERI 4-index transform + Preconditioned CG solver.
///
/// transformERItoMO: AO ERIs → full MO basis (pq|rs).  O(N⁵)
/// solveCPHF: A U^x = B^x via preconditioned conjugate gradient.

use crate::simd_utils::{dot, daxpy};
use crate::integrals::eri_1d_index;

/// Four-quarter AO→MO integral transform.
/// C: MO coefficient matrix (n × n, row-major, C[mu*n+p] = C_{μp})
/// eri_sym: AO ERIs in 8-fold symmetric packed storage
/// n: number of basis functions (= number of MOs)
/// Returns: full (pq|rs) in row-major n⁴ array
pub fn transform_eri_to_mo(c: &[f64], eri_sym: &[f64], n: usize) -> Vec<f64> {
    let n2 = n * n;
    let n3 = n2 * n;
    let n4 = n3 * n;

    // Helper to read from 8-fold symmetric storage
    let eri_get = |i: usize, j: usize, k: usize, l: usize| -> f64 {
        eri_sym[eri_1d_index(i, j, k, l)]
    };

    // Quarter 1: (pν|λσ) = Σ_μ C_{μp} (μν|λσ)
    let mut half1 = vec![0.0; n4];
    for p in 0..n {
        for nu in 0..n {
            for lam in 0..n {
                for sig in 0..n {
                    let mut val = 0.0;
                    for mu in 0..n {
                        val += c[mu * n + p] * eri_get(mu, nu, lam, sig);
                    }
                    half1[((p * n + nu) * n + lam) * n + sig] = val;
                }
            }
        }
    }

    // Quarter 2: (pq|λσ) = Σ_ν C_{νq} (pν|λσ)
    let mut half2 = vec![0.0; n4];
    for p in 0..n {
        for q in 0..n {
            for lam in 0..n {
                for sig in 0..n {
                    let mut val = 0.0;
                    for nu in 0..n {
                        val += c[nu * n + q] * half1[((p * n + nu) * n + lam) * n + sig];
                    }
                    half2[((p * n + q) * n + lam) * n + sig] = val;
                }
            }
        }
    }
    drop(half1);

    // Quarter 3: (pq|rσ) = Σ_λ C_{λr} (pq|λσ)
    let mut half3 = vec![0.0; n4];
    for p in 0..n {
        for q in 0..n {
            for r in 0..n {
                for sig in 0..n {
                    let mut val = 0.0;
                    for lam in 0..n {
                        val += c[lam * n + r] * half2[((p * n + q) * n + lam) * n + sig];
                    }
                    half3[((p * n + q) * n + r) * n + sig] = val;
                }
            }
        }
    }
    drop(half2);

    // Quarter 4: (pq|rs) = Σ_σ C_{σs} (pq|rσ)
    let mut mo_eri = vec![0.0; n4];
    for p in 0..n {
        for q in 0..n {
            for r in 0..n {
                for s in 0..n {
                    let mut val = 0.0;
                    for sig in 0..n {
                        val += c[sig * n + s] * half3[((p * n + q) * n + r) * n + sig];
                    }
                    mo_eri[((p * n + q) * n + r) * n + s] = val;
                }
            }
        }
    }

    mo_eri
}

/// Solve CPHF equations: A U^x = B^x for all perturbations via preconditioned CG.
///
/// mo_eri: full MO-basis ERIs (n⁴, row-major)
/// eps: orbital energies (nmo)
/// rhs: right-hand side vectors, flat [n_pert × dim] where dim = nocc*nvir
/// nocc, nvir, nmo: orbital counts
/// n_pert: number of perturbation directions (3N)
/// tol: convergence tolerance
/// max_iter: maximum CG iterations
/// Returns: solution vectors U^x, flat [n_pert × dim]
pub fn solve_cphf(
    mo_eri: &[f64],
    eps: &[f64],
    rhs: &[f64],
    nocc: usize,
    nvir: usize,
    nmo: usize,
    n_pert: usize,
    tol: f64,
    max_iter: usize,
) -> Vec<f64> {
    let dim = nocc * nvir;

    // Precompute diagonal: ε_a - ε_i
    let mut diagonal = vec![0.0; dim];
    for i in 0..nocc {
        for a in 0..nvir {
            diagonal[i * nvir + a] = eps[nocc + a] - eps[i];
        }
    }

    let mut result = vec![0.0; n_pert * dim];
    let mut r_vec = vec![0.0; dim];
    let mut z_vec = vec![0.0; dim];
    let mut p_vec = vec![0.0; dim];
    let mut ap_vec = vec![0.0; dim];

    for pert in 0..n_pert {
        let b_offset = pert * dim;
        let x_offset = pert * dim;

        // x = 0, r = b
        for idx in 0..dim {
            result[x_offset + idx] = 0.0;
            r_vec[idx] = rhs[b_offset + idx];
        }

        // z = M⁻¹ r (preconditioner: divide by diagonal)
        for idx in 0..dim {
            let d = diagonal[idx];
            z_vec[idx] = if d.abs() > 1e-12 { r_vec[idx] / d } else { r_vec[idx] };
        }

        // p = z
        p_vec.copy_from_slice(&z_vec);

        let mut rz = dot(&r_vec, &z_vec);
        let b_norm = dot(&r_vec, &r_vec).sqrt();
        if b_norm < 1e-15 { continue; }

        for _iter in 0..max_iter {
            // Ap = A p (apply CPHF operator)
            apply_cphf_operator(&diagonal, mo_eri, &p_vec, &mut ap_vec, nocc, nvir, nmo);

            // α = rz / (p·Ap)
            let p_ap = dot(&p_vec, &ap_vec);
            let alpha = rz / if p_ap.abs() > 1e-30 { p_ap } else { 1e-30 };

            // x += α p, r -= α Ap
            daxpy(alpha, &p_vec, &mut result[x_offset..x_offset + dim]);
            daxpy(-alpha, &ap_vec, &mut r_vec);

            // Convergence check
            let r_norm = dot(&r_vec, &r_vec).sqrt();
            if r_norm / b_norm < tol {
                break;
            }

            // z = M⁻¹ r
            for idx in 0..dim {
                let d = diagonal[idx];
                z_vec[idx] = if d.abs() > 1e-12 { r_vec[idx] / d } else { r_vec[idx] };
            }

            // β = (r·z_new) / rz_old
            let rz_new = dot(&r_vec, &z_vec);
            let beta = rz_new / if rz.abs() > 1e-30 { rz } else { 1e-30 };

            // p = z + β p
            for idx in 0..dim {
                p_vec[idx] = z_vec[idx] + beta * p_vec[idx];
            }

            rz = rz_new;
        }
    }

    result
}

/// Apply CPHF operator: output = A * input
/// A_{ai,bj} = (ε_a - ε_i) δ_{ab}δ_{ij} + 4(ai|bj) - (ab|ij) - (aj|ib)
fn apply_cphf_operator(
    diagonal: &[f64],
    mo_eri: &[f64],
    input: &[f64],
    output: &mut [f64],
    nocc: usize,
    nvir: usize,
    nmo: usize,
) {
    let dim = nocc * nvir;

    // Diagonal part
    for idx in 0..dim {
        output[idx] = diagonal[idx] * input[idx];
    }

    // 2-electron part
    let nmo2 = nmo * nmo;
    let nmo3 = nmo2 * nmo;
    for i in 0..nocc {
        for a in 0..nvir {
            let ia = i * nvir + a;
            let a_mo = nocc + a;
            let mut sum = 0.0;
            for j in 0..nocc {
                for b in 0..nvir {
                    let u_bj = input[j * nvir + b];
                    if u_bj.abs() < 1e-15 { continue; }
                    let b_mo = nocc + b;

                    // (ai|bj) in Mulliken notation: mo_eri[a_mo*nmo³ + i*nmo² + b_mo*nmo + j]
                    let aibj = mo_eri[a_mo * nmo3 + i * nmo2 + b_mo * nmo + j];
                    let abij = mo_eri[a_mo * nmo3 + b_mo * nmo2 + i * nmo + j];
                    let ajib = mo_eri[a_mo * nmo3 + j * nmo2 + i * nmo + b_mo];

                    sum += (4.0 * aibj - abij - ajib) * u_bj;
                }
            }
            output[ia] += sum;
        }
    }
}
