/// RI (Resolution of Identity) Coulomb / Exchange / Fock matrix construction.
///
/// B matrix layout: B[P * N * N + mu * N + nu], size = naux * N * N
/// Density / coreH / coefficients: row-major N * N
/// Coefficients C[mu * N + i] (row-major Matrix.data from JS)

use crate::boys::boys_all;
use crate::integrals::{
    PrimShell, ShellGroup, angular_momentums, e_coefficients, group_shells, primitive_norm,
    shell_type_to_num_basis, RBuffer,
};
use crate::simd_utils;

// ---------------------------------------------------------------------------
// RI Setup: 2-center, 3-center integrals, Cholesky, B matrix
// ---------------------------------------------------------------------------

/// Compute the full RI setup: 2c integrals → Cholesky → 3c integrals → B matrix.
///
/// Input:
///   primary_shells_flat: packed PrimShell data (7 fields per shell)
///   primary_norm: CGTO normalization factors for primary basis
///   nbasis: number of primary basis functions
///   aux_shells_flat: packed PrimShell data for auxiliary basis
///   aux_norm: CGTO normalization factors for auxiliary basis
///   naux: number of auxiliary basis functions
///
/// Returns: Vec<f64> with layout [naux_effective, B_matrix...]
///   First element is naux_effective (after Cholesky dropping).
///   Remaining naux_effective * nbasis * nbasis elements are the B matrix.
pub fn compute_ri_setup(
    primary_shells: &[PrimShell],
    primary_norm: &[f64],
    nbasis: usize,
    aux_shells: &[PrimShell],
    aux_norm: &[f64],
    naux: usize,
) -> Vec<f64> {
    let primary_groups = group_shells(primary_shells);
    let aux_groups = group_shells(aux_shells);

    // 1. 2-center integrals (P|Q)
    let vpq = compute_2c_integrals(&aux_groups, aux_norm, naux);

    // 2. Cholesky decomposition with threshold dropping
    let chol_thresh = 1e-8;
    let (l_mat, aux_valid, n_dropped) = cholesky_with_drop(&vpq, naux, chol_thresh);
    let _ = n_dropped; // info only

    // 3. 3-center integrals (μν|P)
    let three = compute_3c_integrals(&primary_groups, primary_norm, nbasis, &aux_groups, aux_norm, naux);

    // 4. B = L^{-1} * three (forward substitution)
    let nn = nbasis * nbasis;
    let mut b = vec![0.0f64; naux * nn];
    for k in 0..nn {
        for i in 0..naux {
            if aux_valid[i] == 0 { continue; }
            let mut sum = three[i * nn + k];
            for j in 0..i {
                if aux_valid[j] == 0 { continue; }
                sum -= l_mat[i * naux + j] * b[j * nn + k];
            }
            b[i * nn + k] = sum / l_mat[i * naux + i];
        }
    }

    // Pack result: [naux as f64, B...]
    // We keep the full naux dimension (dropped rows are zero in B).
    let mut result = Vec::with_capacity(1 + naux * nn);
    result.push(naux as f64);
    result.extend_from_slice(&b);
    result
}

/// 2-center Coulomb integrals: (P|Q)
fn compute_2c_integrals(
    aux_groups: &[ShellGroup],
    aux_norm: &[f64],
    naux: usize,
) -> Vec<f64> {
    let mut v = vec![0.0f64; naux * naux];
    let mut r_buf = RBuffer::new();
    let mut boys_buf = Vec::new();

    for i_grp in 0..aux_groups.len() {
        for j_grp in i_grp..aux_groups.len() {
            compute_2c_shell_pair(
                &aux_groups[i_grp], &aux_groups[j_grp],
                aux_norm, &mut v, naux, i_grp == j_grp,
                &mut r_buf, &mut boys_buf,
            );
        }
    }
    v
}

#[allow(clippy::too_many_arguments)]
fn compute_2c_shell_pair(
    grp_a: &ShellGroup, grp_b: &ShellGroup,
    norm_factors: &[f64], v: &mut [f64], naux: usize, same_group: bool,
    r_buf: &mut RBuffer, boys_buf: &mut Vec<f64>,
) {
    let la = grp_a.shell_type;
    let lb = grp_b.shell_type;
    let ang_a = angular_momentums(la);
    let ang_b = angular_momentums(lb);
    let num_a = shell_type_to_num_basis(la);
    let num_b = shell_type_to_num_basis(lb);

    for ia in 0..num_a {
        let (lxa, lya, lza) = ang_a[ia];
        let mu_a = grp_a.basis_index + ia;
        let ib_start = if same_group { ia } else { 0 };

        for ib in ib_start..num_b {
            let (lxb, lyb, lzb) = ang_b[ib];
            let mu_b = grp_b.basis_index + ib;
            let mut val = 0.0;

            for prim_a in &grp_a.primitives {
                let n_a = primitive_norm(prim_a.exponent, lxa, lya, lza);
                for prim_b in &grp_b.primitives {
                    let n_b = primitive_norm(prim_b.exponent, lxb, lyb, lzb);
                    let alpha = prim_a.exponent;
                    let beta = prim_b.exponent;
                    let p = alpha;
                    let q = beta;

                    let ex_a = e_coefficients(lxa, 0, p, 0.0, 0.0);
                    let ey_a = e_coefficients(lya, 0, p, 0.0, 0.0);
                    let ez_a = e_coefficients(lza, 0, p, 0.0, 0.0);
                    let ex_b = e_coefficients(lxb, 0, q, 0.0, 0.0);
                    let ey_b = e_coefficients(lyb, 0, q, 0.0, 0.0);
                    let ez_b = e_coefficients(lzb, 0, q, 0.0, 0.0);

                    let rpqx = prim_a.x - prim_b.x;
                    let rpqy = prim_a.y - prim_b.y;
                    let rpqz = prim_a.z - prim_b.z;
                    let rpq2 = rpqx * rpqx + rpqy * rpqy + rpqz * rpqz;
                    let eta = p * q / (p + q);

                    let max_n = (lxa + lya + lza + lxb + lyb + lzb) as usize;
                    if boys_buf.len() < max_n + 1 { boys_buf.resize(max_n + 1, 0.0); }
                    boys_all(max_n, eta * rpq2, boys_buf);
                    r_buf.compute(max_n, eta, rpqx, rpqy, rpqz, boys_buf);

                    let mut sum = 0.0;
                    for t1 in 0..=lxa as usize {
                        for u1 in 0..=lya as usize {
                            for v1 in 0..=lza as usize {
                                let e1 = ex_a[t1] * ey_a[u1] * ez_a[v1];
                                if e1.abs() < 1e-15 { continue; }
                                for t2 in 0..=lxb as usize {
                                    for u2 in 0..=lyb as usize {
                                        for v2 in 0..=lzb as usize {
                                            let e2 = ex_b[t2] * ey_b[u2] * ez_b[v2];
                                            if e2.abs() < 1e-15 { continue; }
                                            let sign = if (t2 + u2 + v2) & 1 == 1 { -1.0 } else { 1.0 };
                                            sum += e1 * e2 * sign * r_buf.get(t1 + t2, u1 + u2, v1 + v2, 0);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    let prefactor = 2.0 * std::f64::consts::PI.powf(2.5) / (p * q * (p + q).sqrt());
                    val += prim_a.coefficient * prim_b.coefficient * n_a * n_b * prefactor * sum;
                }
            }

            val *= norm_factors[mu_a] * norm_factors[mu_b];
            v[mu_a * naux + mu_b] += val;
            if mu_a != mu_b {
                v[mu_b * naux + mu_a] = v[mu_a * naux + mu_b];
            }
        }
    }
}

/// Cholesky decomposition with threshold dropping.
/// Returns (L matrix, valid flags, number dropped).
fn cholesky_with_drop(vpq: &[f64], naux: usize, thresh: f64) -> (Vec<f64>, Vec<u8>, usize) {
    let mut l = vec![0.0f64; naux * naux];
    let mut valid = vec![0u8; naux];
    let mut n_dropped = 0usize;

    for i in 0..naux {
        for j in 0..i {
            if valid[j] == 0 { continue; }
            let mut sum = vpq[i * naux + j];
            for k in 0..j {
                if valid[k] == 0 { continue; }
                sum -= l[i * naux + k] * l[j * naux + k];
            }
            l[i * naux + j] = sum / l[j * naux + j];
        }
        let mut diag = vpq[i * naux + i];
        for k in 0..i {
            if valid[k] == 0 { continue; }
            diag -= l[i * naux + k] * l[i * naux + k];
        }
        if diag < thresh {
            l[i * naux + i] = 0.0;
            for j in 0..i { l[i * naux + j] = 0.0; }
            n_dropped += 1;
        } else {
            l[i * naux + i] = diag.sqrt();
            valid[i] = 1;
        }
    }

    (l, valid, n_dropped)
}

/// 3-center Coulomb integrals: (μν|P)
fn compute_3c_integrals(
    primary_groups: &[ShellGroup],
    primary_norm: &[f64],
    nbasis: usize,
    aux_groups: &[ShellGroup],
    aux_norm: &[f64],
    naux: usize,
) -> Vec<f64> {
    let mut three = vec![0.0f64; naux * nbasis * nbasis];
    let mut r_buf = RBuffer::new();
    let mut boys_buf = Vec::new();

    for i_grp in 0..primary_groups.len() {
        for j_grp in i_grp..primary_groups.len() {
            for p_grp in 0..aux_groups.len() {
                compute_3c_shell_triple(
                    &primary_groups[i_grp], &primary_groups[j_grp], &aux_groups[p_grp],
                    primary_norm, aux_norm, nbasis, &mut three, i_grp == j_grp,
                    &mut r_buf, &mut boys_buf,
                );
            }
        }
    }
    three
}

#[allow(clippy::too_many_arguments)]
fn compute_3c_shell_triple(
    grp_a: &ShellGroup, grp_b: &ShellGroup, grp_p: &ShellGroup,
    primary_norm: &[f64], aux_norm: &[f64],
    nbasis: usize, three: &mut [f64], same_ab: bool,
    r_buf: &mut RBuffer, boys_buf: &mut Vec<f64>,
) {
    let la = grp_a.shell_type;
    let lb = grp_b.shell_type;
    let lp = grp_p.shell_type;
    let ang_a = angular_momentums(la);
    let ang_b = angular_momentums(lb);
    let ang_p = angular_momentums(lp);
    let num_a = shell_type_to_num_basis(la);
    let num_b = shell_type_to_num_basis(lb);
    let num_p = shell_type_to_num_basis(lp);
    let nn = nbasis * nbasis;

    for ia in 0..num_a {
        let (lxa, lya, lza) = ang_a[ia];
        let mu_a = grp_a.basis_index + ia;
        let ib_start = if same_ab { ia } else { 0 };

        for ib in ib_start..num_b {
            let (lxb, lyb, lzb) = ang_b[ib];
            let mu_b = grp_b.basis_index + ib;

            for ip in 0..num_p {
                let (lxp, lyp, lzp) = ang_p[ip];
                let mu_p = grp_p.basis_index + ip;
                let mut val = 0.0;

                for prim_a in &grp_a.primitives {
                    let n_a = primitive_norm(prim_a.exponent, lxa, lya, lza);
                    for prim_b in &grp_b.primitives {
                        let n_b = primitive_norm(prim_b.exponent, lxb, lyb, lzb);
                        let alpha = prim_a.exponent;
                        let beta = prim_b.exponent;
                        let p = alpha + beta;
                        let mu_ab = alpha * beta / p;
                        let ab2 = (prim_a.x - prim_b.x).powi(2)
                            + (prim_a.y - prim_b.y).powi(2)
                            + (prim_a.z - prim_b.z).powi(2);
                        let kab = (-mu_ab * ab2).exp();
                        let px = (alpha * prim_a.x + beta * prim_b.x) / p;
                        let py = (alpha * prim_a.y + beta * prim_b.y) / p;
                        let pz = (alpha * prim_a.z + beta * prim_b.z) / p;

                        let ex_ab = e_coefficients(lxa, lxb, p, px - prim_a.x, px - prim_b.x);
                        let ey_ab = e_coefficients(lya, lyb, p, py - prim_a.y, py - prim_b.y);
                        let ez_ab = e_coefficients(lza, lzb, p, pz - prim_a.z, pz - prim_b.z);

                        for prim_p in &grp_p.primitives {
                            let n_p = primitive_norm(prim_p.exponent, lxp, lyp, lzp);
                            let gamma = prim_p.exponent;
                            let q = gamma;

                            let ex_p = e_coefficients(lxp, 0, q, 0.0, 0.0);
                            let ey_p = e_coefficients(lyp, 0, q, 0.0, 0.0);
                            let ez_p = e_coefficients(lzp, 0, q, 0.0, 0.0);

                            let rpqx = px - prim_p.x;
                            let rpqy = py - prim_p.y;
                            let rpqz = pz - prim_p.z;
                            let rpq2 = rpqx * rpqx + rpqy * rpqy + rpqz * rpqz;
                            let eta = p * q / (p + q);

                            let max_n = (lxa + lya + lza + lxb + lyb + lzb + lxp + lyp + lzp) as usize;
                            if boys_buf.len() < max_n + 1 { boys_buf.resize(max_n + 1, 0.0); }
                            boys_all(max_n, eta * rpq2, boys_buf);
                            r_buf.compute(max_n, eta, rpqx, rpqy, rpqz, boys_buf);

                            let mut sum = 0.0;
                            for t1 in 0..=(lxa + lxb) as usize {
                                for u1 in 0..=(lya + lyb) as usize {
                                    for v1 in 0..=(lza + lzb) as usize {
                                        let e1 = ex_ab[t1] * ey_ab[u1] * ez_ab[v1];
                                        if e1.abs() < 1e-15 { continue; }
                                        for t2 in 0..=lxp as usize {
                                            for u2 in 0..=lyp as usize {
                                                for v2 in 0..=lzp as usize {
                                                    let e2 = ex_p[t2] * ey_p[u2] * ez_p[v2];
                                                    if e2.abs() < 1e-15 { continue; }
                                                    let sign = if (t2 + u2 + v2) & 1 == 1 { -1.0 } else { 1.0 };
                                                    sum += e1 * e2 * sign * r_buf.get(t1 + t2, u1 + u2, v1 + v2, 0);
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            let prefactor = 2.0 * std::f64::consts::PI.powf(2.5) / (p * q * (p + q).sqrt());
                            val += prim_a.coefficient * prim_b.coefficient * prim_p.coefficient
                                * n_a * n_b * n_p * kab * prefactor * sum;
                        }
                    }
                }

                val *= primary_norm[mu_a] * primary_norm[mu_b] * aux_norm[mu_p];
                three[mu_p * nn + mu_a * nbasis + mu_b] += val;
                if mu_a != mu_b {
                    three[mu_p * nn + mu_b * nbasis + mu_a] += val;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Coulomb / Exchange / Fock (existing code below)
// ---------------------------------------------------------------------------

/// Compute RI Coulomb matrix J.
///   d_P = Σ_k density[k] * B[P*NN + k]
///   J[k] = Σ_P B[P*NN + k] * d_P
pub fn compute_ri_coulomb(
    b_matrix: &[f64],
    density: &[f64],
    naux: usize,
    nbasis: usize,
) -> Vec<f64> {
    let nn = nbasis * nbasis;

    // d_P = Σ_k density[k] * B[P*NN + k]
    let mut d = vec![0.0f64; naux];
    for p in 0..naux {
        let off = p * nn;
        d[p] = simd_utils::dot(&density[..nn], &b_matrix[off..off + nn]);
    }

    // J[k] = Σ_P B[P*NN + k] * d_P
    let mut j = vec![0.0f64; nn];
    for p in 0..naux {
        let dp = d[p];
        if dp.abs() < 1e-15 {
            continue;
        }
        let off = p * nn;
        simd_utils::daxpy(dp, &b_matrix[off..off + nn], &mut j[..nn]);
    }

    j
}

/// Compute RI Exchange matrix K.
///   X[P*N*nocc + mu*nocc + i] = Σ_λ B[P*NN + mu*N + λ] * C[λ*N + i]
///   K[mu*N + nu] = Σ_P Σ_i X[P*N*nocc + mu*nocc + i] * X[P*N*nocc + nu*nocc + i]
pub fn compute_ri_exchange(
    b_matrix: &[f64],
    coefficients: &[f64],
    naux: usize,
    nbasis: usize,
    nocc: usize,
) -> Vec<f64> {
    let nn = nbasis * nbasis;
    let n = nbasis;

    // X[P][mu][i] half-transform
    let x_size = naux * n * nocc;
    let mut x = vec![0.0f64; x_size];

    for p in 0..naux {
        let b_off = p * nn;
        let x_off = p * n * nocc;
        for mu in 0..n {
            let b_row = b_off + mu * n;
            let x_row = x_off + mu * nocc;
            for i in 0..nocc {
                let mut sum = 0.0;
                for lam in 0..n {
                    sum += b_matrix[b_row + lam] * coefficients[lam * n + i];
                }
                x[x_row + i] = sum;
            }
        }
    }

    // K[mu*N + nu] = Σ_P Σ_i X[P,mu,i] * X[P,nu,i]
    let mut k = vec![0.0f64; nn];
    for p in 0..naux {
        let x_off = p * n * nocc;
        for mu in 0..n {
            let x_mu = x_off + mu * nocc;
            for nu in mu..n {
                let x_nu = x_off + nu * nocc;
                let sum = simd_utils::dot(&x[x_mu..x_mu + nocc], &x[x_nu..x_nu + nocc]);
                k[mu * n + nu] += sum;
                if mu != nu {
                    k[nu * n + mu] += sum;
                }
            }
        }
    }

    k
}

/// Compute RI-RHF Fock matrix: F = H + J - K
pub fn compute_ri_fock_rhf(
    b_matrix: &[f64],
    density: &[f64],
    core_h: &[f64],
    coefficients: &[f64],
    naux: usize,
    nbasis: usize,
    nocc: usize,
) -> Vec<f64> {
    let j = compute_ri_coulomb(b_matrix, density, naux, nbasis);
    let k = compute_ri_exchange(b_matrix, coefficients, naux, nbasis, nocc);
    let nn = nbasis * nbasis;
    let mut f = vec![0.0f64; nn];
    for i in 0..nn {
        f[i] = core_h[i] + j[i] - k[i];
    }
    f
}

/// Compute RI-UHF Fock matrices: Fα = H + J(Ptot) - Kα, Fβ = H + J(Ptot) - Kβ
/// Returns concatenated [Fα..., Fβ...] (2*N*N elements)
pub fn compute_ri_fock_uhf(
    b_matrix: &[f64],
    density_total: &[f64],
    core_h: &[f64],
    coeff_alpha: &[f64],
    nocc_alpha: usize,
    coeff_beta: &[f64],
    nocc_beta: usize,
    naux: usize,
    nbasis: usize,
) -> Vec<f64> {
    let j = compute_ri_coulomb(b_matrix, density_total, naux, nbasis);
    let ka = compute_ri_exchange(b_matrix, coeff_alpha, naux, nbasis, nocc_alpha);
    let kb = compute_ri_exchange(b_matrix, coeff_beta, naux, nbasis, nocc_beta);
    let nn = nbasis * nbasis;
    let mut result = vec![0.0f64; 2 * nn];
    for i in 0..nn {
        let hj = core_h[i] + j[i];
        result[i] = hj - ka[i];
        result[nn + i] = hj - kb[i];
    }
    result
}

/// Compute RI-MP2 energy.
/// B_ia^P = Σ_μν C[μ,i] C[ν,a] B_P^{μν}
/// E_MP2 = Σ_{ij} Σ_{ab} (ia|jb)[2(ia|jb) - (ib|ja)] / (εi + εj - εa - εb)
pub fn compute_ri_mp2(
    b_matrix: &[f64],
    coefficients: &[f64],
    epsilon: &[f64],
    naux: usize,
    nocc: usize,
    nbasis: usize,
) -> f64 {
    let n = nbasis;
    let nn = n * n;
    let nvir = n - nocc;

    // MO transform: bmo[i * nvir * naux + a * naux + P]
    let mut bmo = vec![0.0f64; nocc * nvir * naux];
    let mut b_half = vec![0.0f64; n * nvir];

    for p in 0..naux {
        let b_off = p * nn;

        // Half-transform ν→a: bHalf[μ,a] = Σ_ν B_P^{μν} C[ν, nocc+a]
        b_half.fill(0.0);
        for mu in 0..n {
            for nu in 0..n {
                let b_val = b_matrix[b_off + mu * n + nu];
                if b_val.abs() < 1e-15 { continue; }
                for a in 0..nvir {
                    b_half[mu * nvir + a] += b_val * coefficients[nu * n + nocc + a];
                }
            }
        }

        // Full transform μ→i: bmo[i,a,P] = Σ_μ C[μ,i] * bHalf[μ,a]
        for mu in 0..n {
            for i in 0..nocc {
                let c_mi = coefficients[mu * n + i];
                if c_mi.abs() < 1e-15 { continue; }
                for a in 0..nvir {
                    bmo[i * nvir * naux + a * naux + p] += c_mi * b_half[mu * nvir + a];
                }
            }
        }
    }

    // Energy evaluation
    let mut emp2 = 0.0;
    for i in 0..nocc {
        for j in i..nocc {
            for a in 0..nvir {
                for b in 0..nvir {
                    let ia_off = i * nvir * naux + a * naux;
                    let jb_off = j * nvir * naux + b * naux;
                    let iajb = simd_utils::dot(&bmo[ia_off..ia_off + naux], &bmo[jb_off..jb_off + naux]);

                    let ib_off = i * nvir * naux + b * naux;
                    let ja_off = j * nvir * naux + a * naux;
                    let ibja = simd_utils::dot(&bmo[ib_off..ib_off + naux], &bmo[ja_off..ja_off + naux]);

                    let denom = epsilon[i] + epsilon[j] - epsilon[nocc + a] - epsilon[nocc + b];
                    let factor = if i == j { 1.0 } else { 2.0 };
                    emp2 += factor * iajb * (2.0 * iajb - ibja) / denom;
                }
            }
        }
    }

    emp2
}

/// Compute RI-UMP2 energy (αα + ββ + αβ contributions).
pub fn compute_ri_ump2(
    b_matrix: &[f64],
    ca: &[f64],
    cb: &[f64],
    epsilon_a: &[f64],
    epsilon_b: &[f64],
    naux: usize,
    nocc_a: usize,
    nocc_b: usize,
    nbasis: usize,
) -> f64 {
    let n = nbasis;
    let nn = n * n;
    let nvir_a = n - nocc_a;
    let nvir_b = n - nocc_b;

    // MO transform for alpha
    let bmo_a = ri_mo_transform(b_matrix, ca, naux, n, nn, nocc_a, nvir_a);
    // MO transform for beta
    let bmo_b = ri_mo_transform(b_matrix, cb, naux, n, nn, nocc_b, nvir_b);

    let mut emp2 = 0.0;

    // αα contribution
    for i in 0..nocc_a {
        for j in (i + 1)..nocc_a {
            for a in 0..nvir_a {
                for b in (a + 1)..nvir_a {
                    let (iajb, ibja) = dot_pair(&bmo_a, &bmo_a, i, j, a, b, nvir_a, nvir_a, naux);
                    let antisym = iajb - ibja;
                    let denom = epsilon_a[i] + epsilon_a[j]
                        - epsilon_a[nocc_a + a] - epsilon_a[nocc_a + b];
                    emp2 += antisym * antisym / denom;
                }
            }
        }
    }

    // ββ contribution
    for i in 0..nocc_b {
        for j in (i + 1)..nocc_b {
            for a in 0..nvir_b {
                for b in (a + 1)..nvir_b {
                    let (iajb, ibja) = dot_pair(&bmo_b, &bmo_b, i, j, a, b, nvir_b, nvir_b, naux);
                    let antisym = iajb - ibja;
                    let denom = epsilon_b[i] + epsilon_b[j]
                        - epsilon_b[nocc_b + a] - epsilon_b[nocc_b + b];
                    emp2 += antisym * antisym / denom;
                }
            }
        }
    }

    // αβ contribution
    for i in 0..nocc_a {
        for j in 0..nocc_b {
            for a in 0..nvir_a {
                for b in 0..nvir_b {
                    let ia_off = i * nvir_a * naux + a * naux;
                    let jb_off = j * nvir_b * naux + b * naux;
                    let iajb = simd_utils::dot(&bmo_a[ia_off..ia_off + naux], &bmo_b[jb_off..jb_off + naux]);
                    let denom = epsilon_a[i] + epsilon_b[j]
                        - epsilon_a[nocc_a + a] - epsilon_b[nocc_b + b];
                    emp2 += iajb * iajb / denom;
                }
            }
        }
    }

    emp2
}

/// Half+full MO transform helper
fn ri_mo_transform(
    b_matrix: &[f64],
    coeff: &[f64],
    naux: usize,
    n: usize,
    nn: usize,
    nocc: usize,
    nvir: usize,
) -> Vec<f64> {
    let mut bmo = vec![0.0f64; nocc * nvir * naux];
    let mut b_half = vec![0.0f64; n * nvir];

    for p in 0..naux {
        let b_off = p * nn;
        b_half.fill(0.0);
        for mu in 0..n {
            for nu in 0..n {
                let b_val = b_matrix[b_off + mu * n + nu];
                if b_val.abs() < 1e-15 { continue; }
                for a in 0..nvir {
                    b_half[mu * nvir + a] += b_val * coeff[nu * n + nocc + a];
                }
            }
        }
        for mu in 0..n {
            for i in 0..nocc {
                let c_mi = coeff[mu * n + i];
                if c_mi.abs() < 1e-15 { continue; }
                for a in 0..nvir {
                    bmo[i * nvir * naux + a * naux + p] += c_mi * b_half[mu * nvir + a];
                }
            }
        }
    }

    bmo
}

/// Compute (ia|jb) and (ib|ja) dot products
fn dot_pair(
    bmo1: &[f64], bmo2: &[f64],
    i: usize, j: usize, a: usize, b: usize,
    nvir1: usize, nvir2: usize, naux: usize,
) -> (f64, f64) {
    let ia_off = i * nvir1 * naux + a * naux;
    let jb_off = j * nvir2 * naux + b * naux;
    let ib_off = i * nvir1 * naux + b * naux;
    let ja_off = j * nvir2 * naux + a * naux;
    let iajb = simd_utils::dot(&bmo1[ia_off..ia_off + naux], &bmo2[jb_off..jb_off + naux]);
    let ibja = simd_utils::dot(&bmo1[ib_off..ib_off + naux], &bmo2[ja_off..ja_off + naux]);
    (iajb, ibja)
}
