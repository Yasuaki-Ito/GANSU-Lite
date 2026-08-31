/// McMurchie-Davidson 2-electron integrals — port of integrals2e.ts

use crate::boys::boys_all;

/// Angular momentum cartesian exponents: ANGULAR_MOMENTUMS[shell_type][basis_idx] = (lx,ly,lz)
const ANG_S: &[(i32, i32, i32)] = &[(0, 0, 0)];
const ANG_P: &[(i32, i32, i32)] = &[(1, 0, 0), (0, 1, 0), (0, 0, 1)];
const ANG_D: &[(i32, i32, i32)] = &[
    (2, 0, 0),
    (0, 2, 0),
    (0, 0, 2),
    (1, 1, 0),
    (1, 0, 1),
    (0, 1, 1),
];
const ANG_F: &[(i32, i32, i32)] = &[
    (3, 0, 0),
    (0, 3, 0),
    (0, 0, 3),
    (1, 2, 0),
    (2, 1, 0),
    (2, 0, 1),
    (1, 0, 2),
    (0, 1, 2),
    (0, 2, 1),
    (1, 1, 1),
];

const ANG_G: &[(i32, i32, i32)] = &[
    (4, 0, 0), (0, 4, 0), (0, 0, 4),
    (3, 1, 0), (3, 0, 1), (1, 3, 0), (0, 3, 1), (1, 0, 3), (0, 1, 3),
    (2, 2, 0), (2, 0, 2), (0, 2, 2),
    (2, 1, 1), (1, 2, 1), (1, 1, 2),
];

pub fn angular_momentums(shell_type: u32) -> &'static [(i32, i32, i32)] {
    match shell_type {
        0 => ANG_S,
        1 => ANG_P,
        2 => ANG_D,
        3 => ANG_F,
        4 => ANG_G,
        _ => panic!("unsupported shell type {}", shell_type),
    }
}

pub fn shell_type_to_num_basis(shell_type: u32) -> usize {
    let l = shell_type as usize;
    (l + 1) * (l + 2) / 2
}

/// Double factorial n!!
fn double_factorial(n: i32) -> f64 {
    if n <= 1 {
        return 1.0;
    }
    let mut result = 1.0_f64;
    let mut i = n;
    while i >= 2 {
        result *= i as f64;
        i -= 2;
    }
    result
}

/// Primitive Gaussian normalization factor
pub fn primitive_norm(alpha: f64, lx: i32, ly: i32, lz: i32) -> f64 {
    let l = lx + ly + lz;
    let dbl_fact = double_factorial(2 * lx - 1)
        * double_factorial(2 * ly - 1)
        * double_factorial(2 * lz - 1);
    (2.0_f64).powi(l) / dbl_fact.sqrt()
        * (2.0 / std::f64::consts::PI).powf(0.75)
        * alpha.powf((2 * l + 3) as f64 / 4.0)
}

/// E-coefficients via McMurchie-Davidson recurrence.
/// Returns array of size (i+j+1) indexed by t.
pub fn e_coefficients(i: i32, j: i32, p: f64, xpa: f64, xpb: f64) -> Vec<f64> {
    let max_t = (i + j) as usize;
    let ni = (i + 1) as usize;
    let nj = (j + 1) as usize;
    let nt = max_t + 2;
    let size = ni * nj * nt;
    let mut e_arr = vec![0.0_f64; size];

    let idx = |a: usize, b: usize, t: usize| -> usize { a * nj * nt + b * nt + t };

    e_arr[idx(0, 0, 0)] = 1.0;

    // Build up in a direction
    for a in 0..i as usize {
        for b in 0..=j as usize {
            for t in 0..=a + b + 1 {
                let mut val = xpa * if t <= a + b { e_arr[idx(a, b, t)] } else { 0.0 };
                if t > 0 {
                    val += (1.0 / (2.0 * p)) * e_arr[idx(a, b, t - 1)];
                }
                if t <= a + b {
                    val += (t + 1) as f64
                        * if t + 1 <= a + b {
                            e_arr[idx(a, b, t + 1)]
                        } else {
                            0.0
                        };
                }
                e_arr[idx(a + 1, b, t)] += val;
            }
        }
    }

    // Build up in b direction
    let iu = i as usize;
    for b in 0..j as usize {
        for t in 0..=iu + b + 1 {
            let mut val = xpb * if t <= iu + b { e_arr[idx(iu, b, t)] } else { 0.0 };
            if t > 0 {
                val += (1.0 / (2.0 * p)) * e_arr[idx(iu, b, t - 1)];
            }
            if t <= iu + b {
                val += (t + 1) as f64
                    * if t + 1 <= iu + b {
                        e_arr[idx(iu, b, t + 1)]
                    } else {
                        0.0
                    };
            }
            e_arr[idx(iu, b + 1, t)] += val;
        }
    }

    // Extract result
    let mut result = vec![0.0_f64; max_t + 1];
    for t in 0..=max_t {
        result[t] = e_arr[idx(iu, j as usize, t)];
    }
    result
}

/// R-function buffer (reusable across calls)
pub struct RBuffer {
    buf: Vec<f64>,
    dim: usize,
}

impl RBuffer {
    pub fn new() -> Self {
        Self {
            buf: Vec::new(),
            dim: 0,
        }
    }

    pub fn compute(
        &mut self,
        max_order: usize,
        eta: f64,
        pcx: f64,
        pcy: f64,
        pcz: f64,
        boys: &[f64],
    ) {
        let dim = max_order + 1;
        let d1 = dim + 1;
        let size = d1 * d1 * d1 * d1;

        if size > self.buf.len() {
            self.buf.resize(size, 0.0);
        } else {
            self.buf[..size].fill(0.0);
        }
        self.dim = dim;

        let d2 = d1 * d1;
        let d3 = d2 * d1;

        for n in 0..=max_order {
            self.buf[n * d3] = (-2.0 * eta).powi(n as i32) * boys[n];
        }

        for n in (0..max_order).rev() {
            for t in 0..=max_order - n {
                for u in 0..=max_order - n - t {
                    for v in 0..=max_order - n - t - u {
                        if t + u + v == 0 {
                            continue;
                        }
                        if v > 0 {
                            let mut val =
                                pcz * self.buf[(n + 1) * d3 + t * d2 + u * d1 + v - 1];
                            if v >= 2 {
                                val += (v - 1) as f64
                                    * self.buf[(n + 1) * d3 + t * d2 + u * d1 + v - 2];
                            }
                            self.buf[n * d3 + t * d2 + u * d1 + v] = val;
                        } else if u > 0 {
                            let mut val =
                                pcy * self.buf[(n + 1) * d3 + t * d2 + (u - 1) * d1 + v];
                            if u >= 2 {
                                val += (u - 1) as f64
                                    * self.buf[(n + 1) * d3 + t * d2 + (u - 2) * d1 + v];
                            }
                            self.buf[n * d3 + t * d2 + u * d1 + v] = val;
                        } else {
                            let mut val =
                                pcx * self.buf[(n + 1) * d3 + (t - 1) * d2 + u * d1 + v];
                            if t >= 2 {
                                val += (t - 1) as f64
                                    * self.buf[(n + 1) * d3 + (t - 2) * d2 + u * d1 + v];
                            }
                            self.buf[n * d3 + t * d2 + u * d1 + v] = val;
                        }
                    }
                }
            }
        }
    }

    #[inline]
    pub fn get(&self, t: usize, u: usize, v: usize, n: usize) -> f64 {
        let d1 = self.dim + 1;
        self.buf[n * d1 * d1 * d1 + t * d1 * d1 + u * d1 + v]
    }
}

/// ERI 1D index with 8-fold symmetry
#[inline]
pub fn eri_1d_index(i: usize, j: usize, k: usize, l: usize) -> usize {
    let (ii, jj) = if i >= j { (i, j) } else { (j, i) };
    let (kk, ll) = if k >= l { (k, l) } else { (l, k) };
    let mut ij = ii * (ii + 1) / 2 + jj;
    let mut kl = kk * (kk + 1) / 2 + ll;
    if ij < kl {
        std::mem::swap(&mut ij, &mut kl);
    }
    ij * (ij + 1) / 2 + kl
}

/// Primitive shell data packed for WASM transfer
#[derive(Clone)]
pub struct PrimShell {
    pub exponent: f64,
    pub coefficient: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub shell_type: u32,
    pub basis_index: u32,
    pub atom_index: u32,
}

/// Shell group (contracted shell)
pub struct ShellGroup {
    pub primitives: Vec<PrimShell>,
    pub basis_index: usize,
    pub shell_type: u32,
}

pub fn group_shells(shells: &[PrimShell]) -> Vec<ShellGroup> {
    let mut groups: Vec<ShellGroup> = Vec::new();
    for ps in shells {
        let should_new = groups.last().map_or(true, |g| {
            g.basis_index != ps.basis_index as usize || g.shell_type != ps.shell_type
        });
        if should_new {
            groups.push(ShellGroup {
                primitives: Vec::new(),
                basis_index: ps.basis_index as usize,
                shell_type: ps.shell_type,
            });
        }
        groups.last_mut().unwrap().primitives.push(ps.clone());
    }
    groups
}

/// Compute Schwarz bounds Q[ij] for each shell pair
fn compute_schwarz_bounds(
    shell_groups: &[ShellGroup],
    norm_factors: &[f64],
    num_basis: usize,
) -> Vec<f64> {
    let n_grp = shell_groups.len();
    let mut q = vec![0.0_f64; n_grp * n_grp];

    let eri_size = num_basis * (num_basis + 1) / 2;
    let total_size = eri_size * (eri_size + 1) / 2;
    let mut diag_eri = vec![0.0_f64; total_size];
    let mut r_buf = RBuffer::new();
    let mut boys_buf = Vec::new();

    for i_grp in 0..n_grp {
        for j_grp in i_grp..n_grp {
            diag_eri.iter_mut().for_each(|v| *v = 0.0);
            compute_shell_quartet(
                &shell_groups[i_grp],
                &shell_groups[j_grp],
                &shell_groups[i_grp],
                &shell_groups[j_grp],
                norm_factors,
                num_basis,
                &mut diag_eri,
                i_grp == j_grp,
                i_grp == j_grp,
                true,
                &mut r_buf,
                &mut boys_buf,
            );

            let num_a = shell_type_to_num_basis(shell_groups[i_grp].shell_type);
            let num_b = shell_type_to_num_basis(shell_groups[j_grp].shell_type);
            let mut max_val = 0.0_f64;
            for ia in 0..num_a {
                let mu_a = shell_groups[i_grp].basis_index + ia;
                let ib_start = if i_grp == j_grp { ia } else { 0 };
                for ib in ib_start..num_b {
                    let mu_b = shell_groups[j_grp].basis_index + ib;
                    let idx = eri_1d_index(mu_a, mu_b, mu_a, mu_b);
                    let v = diag_eri[idx].abs();
                    if v > max_val {
                        max_val = v;
                    }
                }
            }

            let qval = max_val.sqrt();
            q[i_grp * n_grp + j_grp] = qval;
            q[j_grp * n_grp + i_grp] = qval;
        }
    }
    q
}

#[allow(clippy::too_many_arguments)]
fn compute_shell_quartet(
    grp_a: &ShellGroup,
    grp_b: &ShellGroup,
    grp_c: &ShellGroup,
    grp_d: &ShellGroup,
    norm_factors: &[f64],
    _num_basis: usize,
    eri: &mut [f64],
    same_ab: bool,
    same_cd: bool,
    same_bra_ket: bool,
    r_buf: &mut RBuffer,
    boys_buf: &mut Vec<f64>,
) {
    let ang_a = angular_momentums(grp_a.shell_type);
    let ang_b = angular_momentums(grp_b.shell_type);
    let ang_c = angular_momentums(grp_c.shell_type);
    let ang_d = angular_momentums(grp_d.shell_type);
    let num_a = ang_a.len();
    let num_b = ang_b.len();
    let num_c = ang_c.len();
    let num_d = ang_d.len();

    // Pre-compute primitive norms
    let p_norms_a = precompute_prim_norms(ang_a, &grp_a.primitives);
    let p_norms_b = precompute_prim_norms(ang_b, &grp_b.primitives);
    let p_norms_c = precompute_prim_norms(ang_c, &grp_c.primitives);
    let p_norms_d = precompute_prim_norms(ang_d, &grp_d.primitives);

    for ia in 0..num_a {
        let (lxa, lya, lza) = ang_a[ia];
        let mu_a = grp_a.basis_index + ia;
        let ib_start = if same_ab { ia } else { 0 };

        for ib in ib_start..num_b {
            let (lxb, lyb, lzb) = ang_b[ib];
            let mu_b = grp_b.basis_index + ib;
            let ij = if mu_a >= mu_b {
                mu_a * (mu_a + 1) / 2 + mu_b
            } else {
                mu_b * (mu_b + 1) / 2 + mu_a
            };

            for ic in 0..num_c {
                let (lxc, lyc, lzc) = ang_c[ic];
                let mu_c = grp_c.basis_index + ic;
                let id_start = if same_cd { ic } else { 0 };

                for id in id_start..num_d {
                    let (lxd, lyd, lzd) = ang_d[id];
                    let mu_d = grp_d.basis_index + id;

                    if same_bra_ket {
                        let kl = if mu_c >= mu_d {
                            mu_c * (mu_c + 1) / 2 + mu_d
                        } else {
                            mu_d * (mu_d + 1) / 2 + mu_c
                        };
                        if ij < kl {
                            continue;
                        }
                    }

                    let mut val = 0.0_f64;

                    for (p_a, prim_a) in grp_a.primitives.iter().enumerate() {
                        for (p_b, prim_b) in grp_b.primitives.iter().enumerate() {
                            let alpha = prim_a.exponent;
                            let beta = prim_b.exponent;
                            let p = alpha + beta;
                            let mu_ab = alpha * beta / p;
                            let (ax, ay, az) = (prim_a.x, prim_a.y, prim_a.z);
                            let (bx, by, bz) = (prim_b.x, prim_b.y, prim_b.z);
                            let ab2 = (ax - bx).powi(2) + (ay - by).powi(2) + (az - bz).powi(2);
                            let kab = (-mu_ab * ab2).exp();
                            let px = (alpha * ax + beta * bx) / p;
                            let py = (alpha * ay + beta * by) / p;
                            let pz = (alpha * az + beta * bz) / p;

                            let ex_ab =
                                e_coefficients(lxa, lxb, p, px - ax, px - bx);
                            let ey_ab =
                                e_coefficients(lya, lyb, p, py - ay, py - by);
                            let ez_ab =
                                e_coefficients(lza, lzb, p, pz - az, pz - bz);

                            for (p_c, prim_c) in grp_c.primitives.iter().enumerate() {
                                for (p_d, prim_d) in grp_d.primitives.iter().enumerate() {
                                    let gamma = prim_c.exponent;
                                    let delta = prim_d.exponent;
                                    let q = gamma + delta;
                                    let mu_cd = gamma * delta / q;
                                    let (cx, cy, cz) = (prim_c.x, prim_c.y, prim_c.z);
                                    let (dx, dy, dz) = (prim_d.x, prim_d.y, prim_d.z);
                                    let cd2 = (cx - dx).powi(2)
                                        + (cy - dy).powi(2)
                                        + (cz - dz).powi(2);
                                    let kcd = (-mu_cd * cd2).exp();
                                    let qx = (gamma * cx + delta * dx) / q;
                                    let qy = (gamma * cy + delta * dy) / q;
                                    let qz = (gamma * cz + delta * dz) / q;

                                    let ex_cd = e_coefficients(lxc, lxd, q, qx - cx, qx - dx);
                                    let ey_cd = e_coefficients(lyc, lyd, q, qy - cy, qy - dy);
                                    let ez_cd = e_coefficients(lzc, lzd, q, qz - cz, qz - dz);

                                    let rpqx = px - qx;
                                    let rpqy = py - qy;
                                    let rpqz = pz - qz;
                                    let rpq2 = rpqx * rpqx + rpqy * rpqy + rpqz * rpqz;
                                    let eta = p * q / (p + q);

                                    let max_n = (lxa + lya + lza + lxb + lyb + lzb + lxc + lyc
                                        + lzc
                                        + lxd
                                        + lyd
                                        + lzd)
                                        as usize;

                                    if boys_buf.len() < max_n + 1 {
                                        boys_buf.resize(max_n + 1, 0.0);
                                    }
                                    boys_all(max_n, eta * rpq2, boys_buf);
                                    r_buf.compute(max_n, eta, rpqx, rpqy, rpqz, boys_buf);

                                    let mut sum = 0.0_f64;
                                    for t1 in 0..=(lxa + lxb) as usize {
                                        for u1 in 0..=(lya + lyb) as usize {
                                            for v1 in 0..=(lza + lzb) as usize {
                                                let e1 = ex_ab[t1] * ey_ab[u1] * ez_ab[v1];
                                                if e1.abs() < 1e-15 {
                                                    continue;
                                                }
                                                for t2 in 0..=(lxc + lxd) as usize {
                                                    for u2 in 0..=(lyc + lyd) as usize {
                                                        for v2 in 0..=(lzc + lzd) as usize {
                                                            let e2 = ex_cd[t2] * ey_cd[u2]
                                                                * ez_cd[v2];
                                                            if e2.abs() < 1e-15 {
                                                                continue;
                                                            }
                                                            let sign = if (t2 + u2 + v2) & 1 == 1
                                                            {
                                                                -1.0
                                                            } else {
                                                                1.0
                                                            };
                                                            sum += e1
                                                                * e2
                                                                * sign
                                                                * r_buf.get(
                                                                    t1 + t2,
                                                                    u1 + u2,
                                                                    v1 + v2,
                                                                    0,
                                                                );
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    let coeff = prim_a.coefficient
                                        * prim_b.coefficient
                                        * prim_c.coefficient
                                        * prim_d.coefficient
                                        * kab
                                        * kcd
                                        * p_norms_a[ia][p_a]
                                        * p_norms_b[ib][p_b]
                                        * p_norms_c[ic][p_c]
                                        * p_norms_d[id][p_d];
                                    let prefactor = 2.0
                                        * std::f64::consts::PI.powf(2.5)
                                        / (p * q * (p + q).sqrt());

                                    val += coeff * prefactor * sum;
                                }
                            }
                        }
                    }

                    val *= norm_factors[mu_a]
                        * norm_factors[mu_b]
                        * norm_factors[mu_c]
                        * norm_factors[mu_d];

                    let idx = eri_1d_index(mu_a, mu_b, mu_c, mu_d);
                    eri[idx] += val;
                }
            }
        }
    }
}

fn precompute_prim_norms(
    ang: &[(i32, i32, i32)],
    prims: &[PrimShell],
) -> Vec<Vec<f64>> {
    ang.iter()
        .map(|&(lx, ly, lz)| {
            prims
                .iter()
                .map(|ps| primitive_norm(ps.exponent, lx, ly, lz))
                .collect()
        })
        .collect()
}

/// Main entry point: compute all ERIs with 8-fold symmetry + Schwarz screening.
/// Uses MD for L<=2, OS for L>2 (auto dispatch per shell quartet).
/// Returns flat f64 array.
pub fn compute_eris(
    shells: &[PrimShell],
    norm_factors: &[f64],
    num_basis: usize,
    schwarz_threshold: f64,
) -> Vec<f64> {
    let shell_groups = group_shells(shells);
    let n_grp = shell_groups.len();

    let eri_size = num_basis * (num_basis + 1) / 2;
    let total_size = eri_size * (eri_size + 1) / 2;
    let mut eri = vec![0.0_f64; total_size];

    let schwarz = compute_schwarz_bounds(&shell_groups, norm_factors, num_basis);

    let mut r_buf = RBuffer::new();
    let mut boys_buf = Vec::new();
    let mut os_vrr_buf = crate::integrals_os::OsVrrBuf::new();

    for i_grp in 0..n_grp {
        for j_grp in i_grp..n_grp {
            let q_ij = schwarz[i_grp * n_grp + j_grp];

            for k_grp in 0..n_grp {
                for l_grp in k_grp..n_grp {
                    let ij_min = i_grp.min(j_grp);
                    let ij_max = i_grp.max(j_grp);
                    let kl_min = k_grp.min(l_grp);
                    let kl_max = k_grp.max(l_grp);
                    if ij_max < kl_max || (ij_max == kl_max && ij_min < kl_min) {
                        continue;
                    }

                    let q_kl = schwarz[k_grp * n_grp + l_grp];
                    if q_ij * q_kl < schwarz_threshold {
                        continue;
                    }

                    let l_total = shell_groups[i_grp].shell_type
                        + shell_groups[j_grp].shell_type
                        + shell_groups[k_grp].shell_type
                        + shell_groups[l_grp].shell_type;

                    if l_total <= 1 {
                        compute_shell_quartet(
                            &shell_groups[i_grp],
                            &shell_groups[j_grp],
                            &shell_groups[k_grp],
                            &shell_groups[l_grp],
                            norm_factors,
                            num_basis,
                            &mut eri,
                            i_grp == j_grp,
                            k_grp == l_grp,
                            i_grp == k_grp && j_grp == l_grp,
                            &mut r_buf,
                            &mut boys_buf,
                        );
                    } else if l_total <= 6 {
                        crate::integrals_os::compute_shell_quartet_os(
                            &shell_groups[i_grp],
                            &shell_groups[j_grp],
                            &shell_groups[k_grp],
                            &shell_groups[l_grp],
                            norm_factors,
                            num_basis,
                            &mut eri,
                            i_grp == j_grp,
                            k_grp == l_grp,
                            i_grp == k_grp && j_grp == l_grp,
                            &mut os_vrr_buf,
                            &mut boys_buf,
                        );
                    } else {
                        crate::integrals_rys::compute_shell_quartet_rys(
                            &shell_groups[i_grp],
                            &shell_groups[j_grp],
                            &shell_groups[k_grp],
                            &shell_groups[l_grp],
                            norm_factors,
                            num_basis,
                            &mut eri,
                            i_grp == j_grp,
                            k_grp == l_grp,
                            i_grp == k_grp && j_grp == l_grp,
                            &mut boys_buf,
                        );
                    }
                }
            }
        }
    }

    eri
}
