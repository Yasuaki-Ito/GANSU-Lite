/// RHF nuclear gradient — 2-electron part.
/// Port of twoElectronGradient from gradient.ts to Rust.
///
/// Computes: grad[3*atom+d] = 0.5 Σ_μνλσ Γ_μνλσ d(μν|λσ)/dR_{atom,d}
/// Using McMurchie-Davidson E-coefficient derivatives with boundary chain rule.

use crate::integrals::{
    angular_momentums, e_coefficients, primitive_norm, group_shells, PrimShell, RBuffer,
};
use crate::boys::boys_all;

/// E-coefficient derivative: d/dA [K E^{ij}_t] / K = 2α E^{i+1,j}_t − i E^{i−1,j}_t
fn e_coeff_deriv(i: i32, j: i32, alpha: f64, p: f64, pa: f64, pb: f64) -> Vec<f64> {
    let max_t = (i + j) as usize;
    let ep1 = e_coefficients(i + 1, j, p, pa, pb);
    let em1 = if i > 0 {
        Some(e_coefficients(i - 1, j, p, pa, pb))
    } else {
        None
    };
    let mut d = vec![0.0; max_t + 1];
    for t in 0..=max_t {
        d[t] = 2.0 * alpha * if t < ep1.len() { ep1[t] } else { 0.0 }
            - if let Some(ref em) = em1 {
                if t < em.len() { i as f64 * em[t] } else { 0.0 }
            } else {
                0.0
            };
    }
    d
}

/// B-center derivative: d/dB [K E^{ij}_t] / K = 2β E^{i,j+1}_t − j E^{i,j−1}_t
fn e_coeff_deriv_b(i: i32, j: i32, beta: f64, p: f64, pa: f64, pb: f64) -> Vec<f64> {
    let max_t = (i + j) as usize;
    let ep1 = e_coefficients(i, j + 1, p, pa, pb);
    let em1 = if j > 0 {
        Some(e_coefficients(i, j - 1, p, pa, pb))
    } else {
        None
    };
    let mut d = vec![0.0; max_t + 1];
    for t in 0..=max_t {
        d[t] = 2.0 * beta * if t < ep1.len() { ep1[t] } else { 0.0 }
            - if let Some(ref em) = em1 {
                if t < em.len() { j as f64 * em[t] } else { 0.0 }
            } else {
                0.0
            };
    }
    d
}

/// Compute 2-electron gradient contribution.
/// shells_flat: [exponent, coefficient, x, y, z, shellType, basisIndex, atomIndex] × nShells
/// density_flat: n×n density matrix (row-major)
/// norms: CGTO normalization factors
/// Returns: gradient vector (3*nAtom)
pub fn compute_2e_gradient(
    shells: &[PrimShell],
    density: &[f64],
    norms: &[f64],
    n: usize,
    n_atom: usize,
) -> Vec<f64> {
    let mut grad = vec![0.0; 3 * n_atom];
    let groups = group_shells(shells);
    let n_grp = groups.len();
    let mut rbuf = RBuffer::new();

    for ig_a in 0..n_grp {
        let g_a = &groups[ig_a];
        for ig_b in 0..n_grp {
            let g_b = &groups[ig_b];
            let ig_ab = ig_a * n_grp + ig_b;
            for ig_c in 0..n_grp {
                let g_c = &groups[ig_c];
                for ig_d in 0..n_grp {
                    let g_d = &groups[ig_d];
                    let ig_cd = ig_c * n_grp + ig_d;
                    // Bra-ket symmetry
                    if ig_ab < ig_cd { continue; }
                    let bra_ket_sym = if ig_ab == ig_cd { 1.0 } else { 2.0 };

                    let ang_a = angular_momentums(g_a.shell_type);
                    let ang_b = angular_momentums(g_b.shell_type);
                    let ang_c = angular_momentums(g_c.shell_type);
                    let ang_d = angular_momentums(g_d.shell_type);

                    for (ia, &(lxa, lya, lza)) in ang_a.iter().enumerate() {
                        let mu = g_a.basis_index + ia;
                        for (ib, &(lxb, lyb, lzb)) in ang_b.iter().enumerate() {
                            let nu = g_b.basis_index + ib;
                            for (ic, &(lxc, lyc, lzc)) in ang_c.iter().enumerate() {
                                let lam = g_c.basis_index + ic;
                                for (id, &(lxd, lyd, lzd)) in ang_d.iter().enumerate() {
                                    let sig = g_d.basis_index + id;

                                    let gamma = density[mu * n + nu] * density[lam * n + sig]
                                        - 0.5 * density[mu * n + lam] * density[nu * n + sig];
                                    if gamma.abs() < 1e-15 { continue; }

                                    for p_a in &g_a.primitives {
                                        let n_a = primitive_norm(p_a.exponent, lxa, lya, lza);
                                        for p_b in &g_b.primitives {
                                            let n_b = primitive_norm(p_b.exponent, lxb, lyb, lzb);
                                            let alpha = p_a.exponent;
                                            let beta = p_b.exponent;
                                            let p = alpha + beta;
                                            let ax = p_a.x; let ay = p_a.y; let az = p_a.z;
                                            let bx = p_b.x; let by = p_b.y; let bz = p_b.z;
                                            let kab = (-(alpha * beta / p) * ((ax-bx).powi(2) + (ay-by).powi(2) + (az-bz).powi(2))).exp();
                                            let px = (alpha*ax + beta*bx) / p;
                                            let py = (alpha*ay + beta*by) / p;
                                            let pz = (alpha*az + beta*bz) / p;

                                            let ex_ab = e_coefficients(lxa, lxb, p, px-ax, px-bx);
                                            let ey_ab = e_coefficients(lya, lyb, p, py-ay, py-by);
                                            let ez_ab = e_coefficients(lza, lzb, p, pz-az, pz-bz);
                                            let dex_a = e_coeff_deriv(lxa, lxb, alpha, p, px-ax, px-bx);
                                            let dey_a = e_coeff_deriv(lya, lyb, alpha, p, py-ay, py-by);
                                            let dez_a = e_coeff_deriv(lza, lzb, alpha, p, pz-az, pz-bz);
                                            let dex_b = e_coeff_deriv_b(lxa, lxb, beta, p, px-ax, px-bx);
                                            let dey_b = e_coeff_deriv_b(lya, lyb, beta, p, py-ay, py-by);
                                            let dez_b = e_coeff_deriv_b(lza, lzb, beta, p, pz-az, pz-bz);

                                            for p_c in &g_c.primitives {
                                                let n_c = primitive_norm(p_c.exponent, lxc, lyc, lzc);
                                                for p_d in &g_d.primitives {
                                                    let n_d = primitive_norm(p_d.exponent, lxd, lyd, lzd);
                                                    let gam = p_c.exponent;
                                                    let del = p_d.exponent;
                                                    let q = gam + del;
                                                    let cx = p_c.x; let cy = p_c.y; let cz = p_c.z;
                                                    let dx = p_d.x; let dy = p_d.y; let dz = p_d.z;
                                                    let kcd = (-(gam * del / q) * ((cx-dx).powi(2) + (cy-dy).powi(2) + (cz-dz).powi(2))).exp();
                                                    let qx = (gam*cx + del*dx) / q;
                                                    let qy = (gam*cy + del*dy) / q;
                                                    let qz = (gam*cz + del*dz) / q;
                                                    let rpqx = px - qx;
                                                    let rpqy = py - qy;
                                                    let rpqz = pz - qz;
                                                    let eta = p * q / (p + q);
                                                    let max_n = (lxa+lya+lza+lxb+lyb+lzb+lxc+lyc+lzc+lxd+lyd+lzd) as usize;

                                                    let mut boys = vec![0.0; max_n + 3];
                                                    boys_all(max_n + 1, eta * (rpqx*rpqx + rpqy*rpqy + rpqz*rpqz), &mut boys);
                                                    rbuf.compute(max_n + 1, eta, rpqx, rpqy, rpqz, &boys);

                                                    let ex_cd = e_coefficients(lxc, lxd, q, qx-cx, qx-dx);
                                                    let ey_cd = e_coefficients(lyc, lyd, q, qy-cy, qy-dy);
                                                    let ez_cd = e_coefficients(lzc, lzd, q, qz-cz, qz-dz);
                                                    let dex_c = e_coeff_deriv(lxc, lxd, gam, q, qx-cx, qx-dx);
                                                    let dey_c = e_coeff_deriv(lyc, lyd, gam, q, qy-cy, qy-dy);
                                                    let dez_c = e_coeff_deriv(lzc, lzd, gam, q, qz-cz, qz-dz);

                                                    let prim_coeff = p_a.coefficient * p_b.coefficient * p_c.coefficient * p_d.coefficient
                                                        * kab * kcd * n_a * n_b * n_c * n_d
                                                        * norms[mu] * norms[nu] * norms[lam] * norms[sig];
                                                    let pre2e = 2.0 * std::f64::consts::PI.powf(2.5) / (p * q * (p + q).sqrt());

                                                    let mt1 = (lxa + lxb) as usize;
                                                    let mu1 = (lya + lyb) as usize;
                                                    let mv1 = (lza + lzb) as usize;
                                                    let mt2 = (lxc + lxd) as usize;
                                                    let mu2 = (lyc + lyd) as usize;
                                                    let mv2 = (lzc + lzd) as usize;

                                                    // Accumulate derivatives for all 4 centers × 3 directions
                                                    let mut ga = [0.0; 3];
                                                    let mut gb = [0.0; 3];
                                                    let mut gc = [0.0; 3];

                                                    for t1 in 0..=mt1 {
                                                        for u1 in 0..=mu1 {
                                                            for v1 in 0..=mv1 {
                                                                let e1 = ex_ab[t1] * ey_ab[u1] * ez_ab[v1];
                                                                for t2 in 0..=mt2 {
                                                                    for u2 in 0..=mu2 {
                                                                        for v2 in 0..=mv2 {
                                                                            let e2 = ex_cd[t2] * ey_cd[u2] * ez_cd[v2];
                                                                            let sign = if (t2 + u2 + v2) & 1 != 0 { -1.0 } else { 1.0 };
                                                                            let rv = sign * rbuf.get(t1+t2, u1+u2, v1+v2, 0);
                                                                            let e12 = e1 * e2;
                                                                            let rx = sign * rbuf.get(t1+t2+1, u1+u2, v1+v2, 0);
                                                                            let ry = sign * rbuf.get(t1+t2, u1+u2+1, v1+v2, 0);
                                                                            let rz = sign * rbuf.get(t1+t2, u1+u2, v1+v2+1, 0);

                                                                            // Center A
                                                                            ga[0] += dex_a[t1]*ey_ab[u1]*ez_ab[v1]*e2*rv;
                                                                            ga[1] += ex_ab[t1]*dey_a[u1]*ez_ab[v1]*e2*rv;
                                                                            ga[2] += ex_ab[t1]*ey_ab[u1]*dez_a[v1]*e2*rv;
                                                                            if t1 == mt1 { ga[0] += e12 * (alpha/p) * rx; }
                                                                            if u1 == mu1 { ga[1] += e12 * (alpha/p) * ry; }
                                                                            if v1 == mv1 { ga[2] += e12 * (alpha/p) * rz; }
                                                                            // Center B
                                                                            gb[0] += dex_b[t1]*ey_ab[u1]*ez_ab[v1]*e2*rv;
                                                                            gb[1] += ex_ab[t1]*dey_b[u1]*ez_ab[v1]*e2*rv;
                                                                            gb[2] += ex_ab[t1]*ey_ab[u1]*dez_b[v1]*e2*rv;
                                                                            if t1 == mt1 { gb[0] += e12 * (beta/p) * rx; }
                                                                            if u1 == mu1 { gb[1] += e12 * (beta/p) * ry; }
                                                                            if v1 == mv1 { gb[2] += e12 * (beta/p) * rz; }
                                                                            // Center C
                                                                            gc[0] += e1*(dex_c[t2]*ey_cd[u2]*ez_cd[v2])*rv;
                                                                            gc[1] += e1*(ex_cd[t2]*dey_c[u2]*ez_cd[v2])*rv;
                                                                            gc[2] += e1*(ex_cd[t2]*ey_cd[u2]*dez_c[v2])*rv;
                                                                            if t2 == mt2 { gc[0] -= e12 * (gam/q) * rx; }
                                                                            if u2 == mu2 { gc[1] -= e12 * (gam/q) * ry; }
                                                                            if v2 == mv2 { gc[2] -= e12 * (gam/q) * rz; }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }

                                                    // d/dD via 4-center TI
                                                    let w = 0.5 * bra_ket_sym * gamma * prim_coeff * pre2e;
                                                    let atom_a = p_a.atom_index as usize;
                                                    let atom_b = p_b.atom_index as usize;
                                                    let atom_c = p_c.atom_index as usize;
                                                    let atom_d = p_d.atom_index as usize;
                                                    for d in 0..3 {
                                                        grad[3*atom_a + d] += w * ga[d];
                                                        grad[3*atom_b + d] += w * gb[d];
                                                        grad[3*atom_c + d] += w * gc[d];
                                                        grad[3*atom_d + d] -= w * (ga[d] + gb[d] + gc[d]);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    grad
}
