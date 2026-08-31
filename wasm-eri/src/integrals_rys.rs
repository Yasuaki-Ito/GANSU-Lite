/// Rys quadrature ERI algorithm — port of JS implementation
/// Modified Chebyshev + Golub-Welsch for roots/weights, 1D VRR+HRR per direction.

use crate::boys::boys_all;
use crate::integrals::{
    angular_momentums, eri_1d_index, primitive_norm, ShellGroup,
};

/// Compute Rys roots and weights from Boys function moments
fn rys_roots_weights(n_roots: usize, t_val: f64, boys_buf: &mut Vec<f64>) -> (Vec<f64>, Vec<f64>) {
    let mut roots = vec![0.0; n_roots];
    let mut weights = vec![0.0; n_roots];
    if n_roots == 0 { return (roots, weights); }

    let max_n = 2 * n_roots - 1;
    if boys_buf.len() < max_n + 1 { boys_buf.resize(max_n + 1, 0.0); }
    boys_all(max_n, t_val, boys_buf);

    let mut alpha = vec![0.0; n_roots];
    let mut beta = vec![0.0; n_roots];

    // Modified Chebyshev
    let mut sig = vec![vec![0.0; 2 * n_roots]; n_roots + 1];
    for k in 0..2 * n_roots { sig[0][k] = boys_buf[k]; }

    alpha[0] = boys_buf[1] / boys_buf[0];
    beta[0] = boys_buf[0];

    for l in 1..n_roots {
        for k in l..2 * n_roots - l {
            sig[l][k] = sig[l - 1][k + 1] - alpha[l - 1] * sig[l - 1][k]
                - if l >= 2 { beta[l - 1] * sig[l - 2][k] } else { 0.0 };
        }
        alpha[l] = sig[l][l + 1] / sig[l][l] - sig[l - 1][l] / sig[l - 1][l - 1];
        beta[l] = sig[l][l] / sig[l - 1][l - 1];
    }

    // Golub-Welsch: QL eigenvalue decomposition
    let mut d = vec![0.0; n_roots];
    let mut e = vec![0.0; n_roots];
    let mut z = vec![0.0; n_roots];
    for i in 0..n_roots {
        d[i] = alpha[i];
        z[i] = if i == 0 { 1.0 } else { 0.0 };
        if i > 0 { e[i] = beta[i].abs().sqrt(); }
    }
    tql_eigenvectors(&mut d, &mut e, &mut z, n_roots);

    for i in 0..n_roots {
        roots[i] = d[i];
        weights[i] = beta[0] * z[i] * z[i];
    }
    (roots, weights)
}

/// QL algorithm with implicit shift for symmetric tridiagonal eigenvalues
fn tql_eigenvectors(d: &mut [f64], e: &mut [f64], z: &mut [f64], n: usize) {
    for l in 0..n {
        let mut iter = 0;
        loop {
            let mut m = l;
            while m < n - 1 {
                let dd = d[m].abs() + d[m + 1].abs();
                if e[m + 1].abs() <= 1e-15 * dd { break; }
                m += 1;
            }
            if m == l { break; }
            iter += 1;
            if iter > 100 { break; }

            let mut g = (d[l + 1] - d[l]) / (2.0 * e[l + 1]);
            let mut r = (g * g + 1.0).sqrt();
            g = d[m] - d[l] + e[l + 1] / (g + if g >= 0.0 { r } else { -r });

            let mut s = 1.0_f64;
            let mut c = 1.0_f64;
            let mut pp = 0.0_f64;
            let mut i = m as isize - 1;
            while i >= l as isize {
                let iu = i as usize;
                let f = s * e[iu + 1];
                let b = c * e[iu + 1];
                r = (f * f + g * g).sqrt();
                if iu + 2 < n { e[iu + 2] = r; }
                if r == 0.0 {
                    d[iu + 1] -= pp;
                    if m + 1 < n { e[m + 1] = 0.0; }
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[iu + 1] - pp;
                r = (d[iu] - g) * s + 2.0 * c * b;
                pp = s * r;
                d[iu + 1] = g + pp;
                g = c * r - b;
                let fz = z[iu + 1];
                z[iu + 1] = s * z[iu] + c * fz;
                z[iu] = c * z[iu] - s * fz;
                i -= 1;
            }
            d[l] -= pp;
            if l + 1 < n { e[l + 1] = g; }
            if m + 1 < n { e[m + 1] = 0.0; }
        }
    }
}

/// 1D Rys integral with VRR + HRR for one Cartesian direction
fn rys_1d(la: i32, lb: i32, lc: i32, ld: i32,
          pa: f64, qc: f64, b10: f64, b01: f64, b00: f64,
          ab: f64, cd: f64) -> f64 {
    let max_a = (la + lb) as usize;
    let max_c = (lc + ld) as usize;
    let str_c = max_c + 1;

    // VRR
    let mut vrr = vec![0.0; (max_a + 1) * str_c];
    vrr[0] = 1.0;
    for c in 0..max_c {
        vrr[c + 1] = qc * vrr[c] + if c > 0 { c as f64 * b01 * vrr[c - 1] } else { 0.0 };
    }
    for a in 0..max_a {
        for c in 0..=max_c {
            vrr[(a + 1) * str_c + c] = pa * vrr[a * str_c + c]
                + if a > 0 { a as f64 * b10 * vrr[(a - 1) * str_c + c] } else { 0.0 }
                + if c > 0 { c as f64 * b00 * vrr[a * str_c + c - 1] } else { 0.0 };
        }
    }

    // Ket HRR
    let mut ket_arr = vrr;
    let mut c_size = max_c + 1;
    for _dd in 0..ld as usize {
        let new_c_size = c_size - 1;
        let mut new_ket = vec![0.0; (max_a + 1) * new_c_size];
        for a in 0..=max_a {
            for c in 0..new_c_size {
                new_ket[a * new_c_size + c] = ket_arr[a * c_size + c + 1] + cd * ket_arr[a * c_size + c];
            }
        }
        ket_arr = new_ket;
        c_size = new_c_size;
    }

    // Bra HRR
    let lc_u = lc as usize;
    let mut bra_arr = vec![0.0; max_a + 1];
    for a in 0..=max_a {
        bra_arr[a] = ket_arr[a * (lc_u + 1) + lc_u];
    }
    for _bb in 0..lb as usize {
        let new_size = bra_arr.len() - 1;
        let mut new_bra = vec![0.0; new_size];
        for a in 0..new_size {
            new_bra[a] = bra_arr[a + 1] + ab * bra_arr[a];
        }
        bra_arr = new_bra;
    }
    bra_arr[la as usize]
}

/// Precompute primitive norms
fn precompute_norms(ang: &[(i32, i32, i32)], prims: &[crate::integrals::PrimShell]) -> Vec<Vec<f64>> {
    ang.iter()
        .map(|&(lx, ly, lz)| prims.iter().map(|ps| primitive_norm(ps.exponent, lx, ly, lz)).collect())
        .collect()
}

/// Compute shell quartet via Rys quadrature
#[allow(clippy::too_many_arguments)]
pub fn compute_shell_quartet_rys(
    grp_a: &ShellGroup, grp_b: &ShellGroup,
    grp_c: &ShellGroup, grp_d: &ShellGroup,
    norm_factors: &[f64], _num_basis: usize, eri: &mut [f64],
    same_ab: bool, same_cd: bool, same_bra_ket: bool,
    boys_buf: &mut Vec<f64>,
) {
    let la = grp_a.shell_type as i32;
    let lb = grp_b.shell_type as i32;
    let lc = grp_c.shell_type as i32;
    let ld = grp_d.shell_type as i32;
    let ang_a = angular_momentums(grp_a.shell_type);
    let ang_b = angular_momentums(grp_b.shell_type);
    let ang_c = angular_momentums(grp_c.shell_type);
    let ang_d = angular_momentums(grp_d.shell_type);
    let num_a = ang_a.len();
    let num_b = ang_b.len();
    let num_c = ang_c.len();
    let num_d = ang_d.len();
    let n_roots = ((la + lb + lc + ld + 1) as usize + 1) / 2;
    let p_norms_a = precompute_norms(ang_a, &grp_a.primitives);
    let p_norms_b = precompute_norms(ang_b, &grp_b.primitives);
    let p_norms_c = precompute_norms(ang_c, &grp_c.primitives);
    let p_norms_d = precompute_norms(ang_d, &grp_d.primitives);

    for ia in 0..num_a {
        let (lxa, lya, lza) = ang_a[ia];
        let mu_a = grp_a.basis_index + ia;
        let ib_start = if same_ab { ia } else { 0 };
        for ib in ib_start..num_b {
            let (lxb, lyb, lzb) = ang_b[ib];
            let mu_b = grp_b.basis_index + ib;
            let ij = if mu_a >= mu_b { mu_a * (mu_a + 1) / 2 + mu_b } else { mu_b * (mu_b + 1) / 2 + mu_a };
            for ic in 0..num_c {
                let (lxc, lyc, lzc) = ang_c[ic];
                let mu_c = grp_c.basis_index + ic;
                let id_start = if same_cd { ic } else { 0 };
                for id in id_start..num_d {
                    let (lxd, lyd, lzd) = ang_d[id];
                    let mu_d = grp_d.basis_index + id;
                    if same_bra_ket {
                        let kl = if mu_c >= mu_d { mu_c * (mu_c + 1) / 2 + mu_d } else { mu_d * (mu_d + 1) / 2 + mu_c };
                        if ij < kl { continue; }
                    }

                    let mut val = 0.0_f64;
                    for (p_a, prim_a) in grp_a.primitives.iter().enumerate() {
                        let alpha = prim_a.exponent;
                        let (ax, ay, az) = (prim_a.x, prim_a.y, prim_a.z);
                        for (p_b, prim_b) in grp_b.primitives.iter().enumerate() {
                            let beta_v = prim_b.exponent;
                            let (bx, by, bz) = (prim_b.x, prim_b.y, prim_b.z);
                            let p = alpha + beta_v;
                            let mu_ab = alpha * beta_v / p;
                            let kab = (-mu_ab * ((ax - bx).powi(2) + (ay - by).powi(2) + (az - bz).powi(2))).exp();
                            let px = (alpha * ax + beta_v * bx) / p;
                            let py = (alpha * ay + beta_v * by) / p;
                            let pz = (alpha * az + beta_v * bz) / p;

                            for (p_c, prim_c) in grp_c.primitives.iter().enumerate() {
                                let gamma_v = prim_c.exponent;
                                let (cx, cy, cz) = (prim_c.x, prim_c.y, prim_c.z);
                                for (p_d, prim_d) in grp_d.primitives.iter().enumerate() {
                                    let delta_v = prim_d.exponent;
                                    let (dx, dy, dz) = (prim_d.x, prim_d.y, prim_d.z);
                                    let q = gamma_v + delta_v;
                                    let mu_cd = gamma_v * delta_v / q;
                                    let kcd = (-mu_cd * ((cx - dx).powi(2) + (cy - dy).powi(2) + (cz - dz).powi(2))).exp();
                                    let qx = (gamma_v * cx + delta_v * dx) / q;
                                    let qy = (gamma_v * cy + delta_v * dy) / q;
                                    let qz = (gamma_v * cz + delta_v * dz) / q;

                                    let pqx = px - qx;
                                    let pqy = py - qy;
                                    let pqz = pz - qz;
                                    let t_val = p * q / (p + q) * (pqx * pqx + pqy * pqy + pqz * pqz);

                                    let (rys_roots, rys_weights) = rys_roots_weights(n_roots, t_val, boys_buf);
                                    let ab_x = ax - bx;
                                    let ab_y = ay - by;
                                    let ab_z = az - bz;
                                    let cd_x = cx - dx;
                                    let cd_y = cy - dy;
                                    let cd_z = cz - dz;

                                    let mut sum = 0.0_f64;
                                    for r in 0..n_roots {
                                        let t2 = rys_roots[r];
                                        let b00 = t2 / (2.0 * (p + q));
                                        let b10 = (1.0 - t2 * q / (p + q)) / (2.0 * p);
                                        let b01 = (1.0 - t2 * p / (p + q)) / (2.0 * q);
                                        let pa_x = px - ax + t2 * (qx - px) * q / (p + q);
                                        let pa_y = py - ay + t2 * (qy - py) * q / (p + q);
                                        let pa_z = pz - az + t2 * (qz - pz) * q / (p + q);
                                        let qc_x = qx - cx + t2 * (px - qx) * p / (p + q);
                                        let qc_y = qy - cy + t2 * (py - qy) * p / (p + q);
                                        let qc_z = qz - cz + t2 * (pz - qz) * p / (p + q);

                                        let ix = rys_1d(lxa, lxb, lxc, lxd, pa_x, qc_x, b10, b01, b00, ab_x, cd_x);
                                        let iy = rys_1d(lya, lyb, lyc, lyd, pa_y, qc_y, b10, b01, b00, ab_y, cd_y);
                                        let iz = rys_1d(lza, lzb, lzc, lzd, pa_z, qc_z, b10, b01, b00, ab_z, cd_z);
                                        sum += rys_weights[r] * ix * iy * iz;
                                    }

                                    val += prim_a.coefficient * prim_b.coefficient
                                        * prim_c.coefficient * prim_d.coefficient
                                        * kab * kcd
                                        * p_norms_a[ia][p_a] * p_norms_b[ib][p_b]
                                        * p_norms_c[ic][p_c] * p_norms_d[id][p_d]
                                        * 2.0 * std::f64::consts::PI.powf(2.5)
                                        / (p * q * (p + q).sqrt()) * sum;
                                }
                            }
                        }
                    }

                    let idx = eri_1d_index(mu_a, mu_b, mu_c, mu_d);
                    eri[idx] += val * norm_factors[mu_a] * norm_factors[mu_b]
                        * norm_factors[mu_c] * norm_factors[mu_d];
                }
            }
        }
    }
}
