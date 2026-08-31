/// Obara-Saika ERI algorithm — 6-phase VRR + recursive HRR
/// Port of the JS implementation in integrals2e.ts

use std::collections::HashMap;
use crate::boys::boys_all;
use crate::integrals::{
    angular_momentums, eri_1d_index, primitive_norm,
    PrimShell, ShellGroup,
};

/// Precompute primitive norms for a shell group
fn precompute_norms(ang: &[(i32, i32, i32)], prims: &[PrimShell]) -> Vec<Vec<f64>> {
    ang.iter()
        .map(|&(lx, ly, lz)| {
            prims.iter().map(|ps| primitive_norm(ps.exponent, lx, ly, lz)).collect()
        })
        .collect()
}

/// VRR buffer (module-level reuse via caller)
pub struct OsVrrBuf {
    data: Vec<f64>,
    sa: usize,  // strides
    say: usize,
    saz: usize,
    scx: usize,
    scy: usize,
    scz: usize,
}

impl OsVrrBuf {
    pub fn new() -> Self {
        Self { data: Vec::new(), sa: 0, say: 0, saz: 0, scx: 0, scy: 0, scz: 0 }
    }

    /// Compute VRR for all [ax,ay,az | cx,cy,cz]^(m) needed
    pub fn compute(
        &mut self,
        lab: usize, lcd: usize,
        p: f64, q: f64,
        pax: f64, pay: f64, paz: f64,
        wpx: f64, wpy: f64, wpz: f64,
        qcx: f64, qcy: f64, qcz: f64,
        wqx: f64, wqy: f64, wqz: f64,
        boys: &[f64],
    ) {
        let rho_p = q / (p + q);
        let rho_q = p / (p + q);
        let inv2p = 0.5 / p;
        let inv2q = 0.5 / q;
        let inv2pq = 0.5 / (p + q);

        let sa = lab + 1;
        let sc = lcd + 1;
        let sm = lab + lcd + 1;
        let scz_s = sm;
        let scy_s = sc * scz_s;
        let scx_s = sc * scy_s;
        let saz_s = sc * scx_s;
        let say_s = sa * saz_s;
        let sax_s = sa * say_s;
        let total = sa * sax_s;

        if self.data.len() < total {
            self.data.resize(total, 0.0);
        } else {
            self.data[..total].fill(0.0);
        }
        self.sa = sax_s;
        self.say = say_s;
        self.saz = saz_s;
        self.scx = scx_s;
        self.scy = scy_s;
        self.scz = scz_s;

        let vi = |ax: usize, ay: usize, az: usize, cx: usize, cy: usize, cz: usize, m: usize| -> usize {
            ax * sax_s + ay * say_s + az * saz_s + cx * scx_s + cy * scy_s + cz * scz_s + m
        };

        // Base case
        for m in 0..=lab + lcd {
            self.data[vi(0, 0, 0, 0, 0, 0, m)] = boys[m];
        }

        // Phase 1: Bra x
        for ax in 0..lab {
            let max_m = lab + lcd - ax - 1;
            for m in 0..=max_m {
                let mut v = pax * self.data[vi(ax, 0, 0, 0, 0, 0, m)]
                    + wpx * self.data[vi(ax, 0, 0, 0, 0, 0, m + 1)];
                if ax > 0 {
                    v += ax as f64 * inv2p * (self.data[vi(ax - 1, 0, 0, 0, 0, 0, m)]
                        - rho_p * self.data[vi(ax - 1, 0, 0, 0, 0, 0, m + 1)]);
                }
                self.data[vi(ax + 1, 0, 0, 0, 0, 0, m)] = v;
            }
        }

        // Phase 2: Bra y
        for ax in 0..=lab {
            for ay in 0..lab.saturating_sub(ax) {
                let max_m = lab + lcd - ax - ay - 1;
                for m in 0..=max_m {
                    let mut v = pay * self.data[vi(ax, ay, 0, 0, 0, 0, m)]
                        + wpy * self.data[vi(ax, ay, 0, 0, 0, 0, m + 1)];
                    if ay > 0 {
                        v += ay as f64 * inv2p * (self.data[vi(ax, ay - 1, 0, 0, 0, 0, m)]
                            - rho_p * self.data[vi(ax, ay - 1, 0, 0, 0, 0, m + 1)]);
                    }
                    self.data[vi(ax, ay + 1, 0, 0, 0, 0, m)] = v;
                }
            }
        }

        // Phase 3: Bra z
        for ax in 0..=lab {
            for ay in 0..=lab - ax {
                for az in 0..lab.saturating_sub(ax + ay) {
                    let max_m = lab + lcd - ax - ay - az - 1;
                    for m in 0..=max_m {
                        let mut v = paz * self.data[vi(ax, ay, az, 0, 0, 0, m)]
                            + wpz * self.data[vi(ax, ay, az, 0, 0, 0, m + 1)];
                        if az > 0 {
                            v += az as f64 * inv2p * (self.data[vi(ax, ay, az - 1, 0, 0, 0, m)]
                                - rho_p * self.data[vi(ax, ay, az - 1, 0, 0, 0, m + 1)]);
                        }
                        self.data[vi(ax, ay, az + 1, 0, 0, 0, m)] = v;
                    }
                }
            }
        }

        // Phase 4: Ket x (coupling: ax)
        for cx in 0..lcd {
            for ax in 0..=lab {
                for ay in 0..=lab - ax {
                    for az in 0..=lab - ax - ay {
                        let max_m = lcd - cx - 1;
                        for m in 0..=max_m {
                            let mut v = qcx * self.data[vi(ax, ay, az, cx, 0, 0, m)]
                                + wqx * self.data[vi(ax, ay, az, cx, 0, 0, m + 1)];
                            if cx > 0 {
                                v += cx as f64 * inv2q * (self.data[vi(ax, ay, az, cx - 1, 0, 0, m)]
                                    - rho_q * self.data[vi(ax, ay, az, cx - 1, 0, 0, m + 1)]);
                            }
                            if ax > 0 {
                                v += ax as f64 * inv2pq * self.data[vi(ax - 1, ay, az, cx, 0, 0, m + 1)];
                            }
                            self.data[vi(ax, ay, az, cx + 1, 0, 0, m)] = v;
                        }
                    }
                }
            }
        }

        // Phase 5: Ket y (coupling: ay)
        for cx in 0..=lcd {
            for cy in 0..lcd.saturating_sub(cx) {
                for ax in 0..=lab {
                    for ay in 0..=lab - ax {
                        for az in 0..=lab - ax - ay {
                            let max_m = lcd - cx - cy - 1;
                            for m in 0..=max_m {
                                let mut v = qcy * self.data[vi(ax, ay, az, cx, cy, 0, m)]
                                    + wqy * self.data[vi(ax, ay, az, cx, cy, 0, m + 1)];
                                if cy > 0 {
                                    v += cy as f64 * inv2q * (self.data[vi(ax, ay, az, cx, cy - 1, 0, m)]
                                        - rho_q * self.data[vi(ax, ay, az, cx, cy - 1, 0, m + 1)]);
                                }
                                if ay > 0 {
                                    v += ay as f64 * inv2pq * self.data[vi(ax, ay - 1, az, cx, cy, 0, m + 1)];
                                }
                                self.data[vi(ax, ay, az, cx, cy + 1, 0, m)] = v;
                            }
                        }
                    }
                }
            }
        }

        // Phase 6: Ket z (coupling: az)
        for cx in 0..=lcd {
            for cy in 0..=lcd - cx {
                for cz in 0..lcd.saturating_sub(cx + cy) {
                    for ax in 0..=lab {
                        for ay in 0..=lab - ax {
                            for az in 0..=lab - ax - ay {
                                let max_m = lcd - cx - cy - cz - 1;
                                for m in 0..=max_m {
                                    let mut v = qcz * self.data[vi(ax, ay, az, cx, cy, cz, m)]
                                        + wqz * self.data[vi(ax, ay, az, cx, cy, cz, m + 1)];
                                    if cz > 0 {
                                        v += cz as f64 * inv2q * (self.data[vi(ax, ay, az, cx, cy, cz - 1, m)]
                                            - rho_q * self.data[vi(ax, ay, az, cx, cy, cz - 1, m + 1)]);
                                    }
                                    if az > 0 {
                                        v += az as f64 * inv2pq * self.data[vi(ax, ay, az - 1, cx, cy, cz, m + 1)];
                                    }
                                    self.data[vi(ax, ay, az, cx, cy, cz + 1, m)] = v;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// Read VRR value at m=0
    #[inline]
    pub fn get(&self, ax: usize, ay: usize, az: usize, cx: usize, cy: usize, cz: usize) -> f64 {
        self.data[ax * self.sa + ay * self.say + az * self.saz
            + cx * self.scx + cy * self.scy + cz * self.scz]
    }
}

/// HRR key encoder (12 components × 4 bits = 48 bits)
#[inline]
fn hrr_key(ax: u8, ay: u8, az: u8, bx: u8, by: u8, bz: u8,
           cx: u8, cy: u8, cz: u8, dx: u8, dy: u8, dz: u8) -> u64 {
    (ax as u64) << 44 | (ay as u64) << 40 | (az as u64) << 36
    | (bx as u64) << 32 | (by as u64) << 28 | (bz as u64) << 24
    | (cx as u64) << 20 | (cy as u64) << 16 | (cz as u64) << 12
    | (dx as u64) << 8 | (dy as u64) << 4 | (dz as u64)
}

/// Recursive HRR with memoization
fn hrr_eval(
    ax: i32, ay: i32, az: i32, bx: i32, by: i32, bz: i32,
    cx: i32, cy: i32, cz: i32, dx: i32, dy: i32, dz: i32,
    vrr: &OsVrrBuf,
    abx: f64, aby: f64, abz: f64, cdx: f64, cdy: f64, cdz: f64,
    memo: &mut HashMap<u64, f64>,
) -> f64 {
    let key = hrr_key(ax as u8, ay as u8, az as u8, bx as u8, by as u8, bz as u8,
                      cx as u8, cy as u8, cz as u8, dx as u8, dy as u8, dz as u8);
    if let Some(&v) = memo.get(&key) {
        return v;
    }

    let result = if bx + by + bz > 0 {
        // Bra HRR
        if bx > 0 {
            hrr_eval(ax + 1, ay, az, bx - 1, by, bz, cx, cy, cz, dx, dy, dz, vrr, abx, aby, abz, cdx, cdy, cdz, memo)
            + abx * hrr_eval(ax, ay, az, bx - 1, by, bz, cx, cy, cz, dx, dy, dz, vrr, abx, aby, abz, cdx, cdy, cdz, memo)
        } else if by > 0 {
            hrr_eval(ax, ay + 1, az, bx, by - 1, bz, cx, cy, cz, dx, dy, dz, vrr, abx, aby, abz, cdx, cdy, cdz, memo)
            + aby * hrr_eval(ax, ay, az, bx, by - 1, bz, cx, cy, cz, dx, dy, dz, vrr, abx, aby, abz, cdx, cdy, cdz, memo)
        } else {
            hrr_eval(ax, ay, az + 1, bx, by, bz - 1, cx, cy, cz, dx, dy, dz, vrr, abx, aby, abz, cdx, cdy, cdz, memo)
            + abz * hrr_eval(ax, ay, az, bx, by, bz - 1, cx, cy, cz, dx, dy, dz, vrr, abx, aby, abz, cdx, cdy, cdz, memo)
        }
    } else if dx + dy + dz > 0 {
        // Ket HRR
        if dx > 0 {
            hrr_eval(ax, ay, az, 0, 0, 0, cx + 1, cy, cz, dx - 1, dy, dz, vrr, abx, aby, abz, cdx, cdy, cdz, memo)
            + cdx * hrr_eval(ax, ay, az, 0, 0, 0, cx, cy, cz, dx - 1, dy, dz, vrr, abx, aby, abz, cdx, cdy, cdz, memo)
        } else if dy > 0 {
            hrr_eval(ax, ay, az, 0, 0, 0, cx, cy + 1, cz, dx, dy - 1, dz, vrr, abx, aby, abz, cdx, cdy, cdz, memo)
            + cdy * hrr_eval(ax, ay, az, 0, 0, 0, cx, cy, cz, dx, dy - 1, dz, vrr, abx, aby, abz, cdx, cdy, cdz, memo)
        } else {
            hrr_eval(ax, ay, az, 0, 0, 0, cx, cy, cz + 1, dx, dy, dz - 1, vrr, abx, aby, abz, cdx, cdy, cdz, memo)
            + cdz * hrr_eval(ax, ay, az, 0, 0, 0, cx, cy, cz, dx, dy, dz - 1, vrr, abx, aby, abz, cdx, cdy, cdz, memo)
        }
    } else {
        vrr.get(ax as usize, ay as usize, az as usize, cx as usize, cy as usize, cz as usize)
    };

    memo.insert(key, result);
    result
}

/// Compute shell quartet via Obara-Saika VRR + HRR
#[allow(clippy::too_many_arguments)]
pub fn compute_shell_quartet_os(
    grp_a: &ShellGroup, grp_b: &ShellGroup,
    grp_c: &ShellGroup, grp_d: &ShellGroup,
    norm_factors: &[f64], _num_basis: usize, eri: &mut [f64],
    same_ab: bool, same_cd: bool, same_bra_ket: bool,
    vrr_buf: &mut OsVrrBuf, boys_buf: &mut Vec<f64>,
) {
    let la = grp_a.shell_type as usize;
    let lb = grp_b.shell_type as usize;
    let lc = grp_c.shell_type as usize;
    let ld = grp_d.shell_type as usize;
    let lab = la + lb;
    let lcd = lc + ld;
    let ang_a = angular_momentums(grp_a.shell_type);
    let ang_b = angular_momentums(grp_b.shell_type);
    let ang_c = angular_momentums(grp_c.shell_type);
    let ang_d = angular_momentums(grp_d.shell_type);
    let num_a = ang_a.len();
    let num_b = ang_b.len();
    let num_c = ang_c.len();
    let num_d = ang_d.len();
    let p_norms_a = precompute_norms(ang_a, &grp_a.primitives);
    let p_norms_b = precompute_norms(ang_b, &grp_b.primitives);
    let p_norms_c = precompute_norms(ang_c, &grp_c.primitives);
    let p_norms_d = precompute_norms(ang_d, &grp_d.primitives);

    let mut contracted = vec![0.0_f64; num_a * num_b * num_c * num_d];

    for (p_a, prim_a) in grp_a.primitives.iter().enumerate() {
        let alpha = prim_a.exponent;
        let (ax, ay, az) = (prim_a.x, prim_a.y, prim_a.z);
        for (p_b, prim_b) in grp_b.primitives.iter().enumerate() {
            let beta = prim_b.exponent;
            let (bx, by, bz) = (prim_b.x, prim_b.y, prim_b.z);
            let p = alpha + beta;
            let mu_ab = alpha * beta / p;
            let ab2 = (ax - bx).powi(2) + (ay - by).powi(2) + (az - bz).powi(2);
            let kab = (-mu_ab * ab2).exp();
            let px = (alpha * ax + beta * bx) / p;
            let py = (alpha * ay + beta * by) / p;
            let pz = (alpha * az + beta * bz) / p;

            for (p_c, prim_c) in grp_c.primitives.iter().enumerate() {
                let gamma = prim_c.exponent;
                let (cx, cy, cz) = (prim_c.x, prim_c.y, prim_c.z);
                for (p_d, prim_d) in grp_d.primitives.iter().enumerate() {
                    let delta = prim_d.exponent;
                    let (dx, dy, dz) = (prim_d.x, prim_d.y, prim_d.z);
                    let q = gamma + delta;
                    let mu_cd = gamma * delta / q;
                    let cd2 = (cx - dx).powi(2) + (cy - dy).powi(2) + (cz - dz).powi(2);
                    let kcd = (-mu_cd * cd2).exp();
                    let qx = (gamma * cx + delta * dx) / q;
                    let qy = (gamma * cy + delta * dy) / q;
                    let qz = (gamma * cz + delta * dz) / q;

                    let eta = p * q / (p + q);
                    let pqx = px - qx;
                    let pqy = py - qy;
                    let pqz = pz - qz;
                    let t_val = eta * (pqx * pqx + pqy * pqy + pqz * pqz);

                    let max_n = lab + lcd;
                    if boys_buf.len() < max_n + 1 {
                        boys_buf.resize(max_n + 1, 0.0);
                    }
                    boys_all(max_n, t_val, boys_buf);

                    let prefactor = 2.0 * std::f64::consts::PI.powf(2.5)
                        / (p * q * (p + q).sqrt()) * kab * kcd;

                    let wx = (p * px + q * qx) / (p + q);
                    let wy = (p * py + q * qy) / (p + q);
                    let wz = (p * pz + q * qz) / (p + q);

                    vrr_buf.compute(lab, lcd, p, q,
                        px - ax, py - ay, pz - az,
                        wx - px, wy - py, wz - pz,
                        qx - cx, qy - cy, qz - cz,
                        wx - qx, wy - qy, wz - qz,
                        boys_buf);

                    let ab_x = ax - bx;
                    let ab_y = ay - by;
                    let ab_z = az - bz;
                    let cd_x = cx - dx;
                    let cd_y = cy - dy;
                    let cd_z = cz - dz;
                    let mut memo = HashMap::new();

                    for ia in 0..num_a {
                        let (lxa, lya, lza) = ang_a[ia];
                        for ib in 0..num_b {
                            let (lxb, lyb, lzb) = ang_b[ib];
                            let c_ab = prim_a.coefficient * prim_b.coefficient
                                * p_norms_a[ia][p_a] * p_norms_b[ib][p_b];
                            for ic in 0..num_c {
                                let (lxc, lyc, lzc) = ang_c[ic];
                                for id in 0..num_d {
                                    let (lxd, lyd, lzd) = ang_d[id];
                                    let c_cd = prim_c.coefficient * prim_d.coefficient
                                        * p_norms_c[ic][p_c] * p_norms_d[id][p_d];
                                    let val = hrr_eval(
                                        lxa, lya, lza, lxb, lyb, lzb,
                                        lxc, lyc, lzc, lxd, lyd, lzd,
                                        vrr_buf, ab_x, ab_y, ab_z, cd_x, cd_y, cd_z,
                                        &mut memo);
                                    contracted[(ia * num_b + ib) * num_c * num_d + ic * num_d + id]
                                        += prefactor * c_ab * c_cd * val;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Store into ERI with symmetry
    for ia in 0..num_a {
        let mu_a = grp_a.basis_index + ia;
        let ib_start = if same_ab { ia } else { 0 };
        for ib in ib_start..num_b {
            let mu_b = grp_b.basis_index + ib;
            let ij = if mu_a >= mu_b { mu_a * (mu_a + 1) / 2 + mu_b } else { mu_b * (mu_b + 1) / 2 + mu_a };
            for ic in 0..num_c {
                let mu_c = grp_c.basis_index + ic;
                let id_start = if same_cd { ic } else { 0 };
                for id in id_start..num_d {
                    let mu_d = grp_d.basis_index + id;
                    if same_bra_ket {
                        let kl = if mu_c >= mu_d { mu_c * (mu_c + 1) / 2 + mu_d } else { mu_d * (mu_d + 1) / 2 + mu_c };
                        if ij < kl { continue; }
                    }
                    let idx = eri_1d_index(mu_a, mu_b, mu_c, mu_d);
                    eri[idx] += contracted[(ia * num_b + ib) * num_c * num_d + ic * num_d + id]
                        * norm_factors[mu_a] * norm_factors[mu_b]
                        * norm_factors[mu_c] * norm_factors[mu_d];
                }
            }
        }
    }
}
