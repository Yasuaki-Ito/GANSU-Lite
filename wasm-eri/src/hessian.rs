/// RHF 2-electron skeleton Hessian — port of twoElectronHessian from analyticalHessian.ts.
///
/// Computes: H[3a+d1, 3b+d2] = Σ_μνλσ Γ d²(μν|λσ)/dR_{a,d1} dR_{b,d2}
/// Using McMurchie-Davidson E-coefficient shift formulas.
///
/// Optimised: RBuffer/boys_buf allocated once, primitive norms precomputed,
/// E-coefficient cache uses stack arrays.

use crate::boys::boys_all;
use crate::integrals::{
    angular_momentums, e_coefficients, primitive_norm, group_shells, PrimShell, RBuffer,
};

/// Max E-coefficient length for shifts up to ±2 with L≤4 (g-functions).
/// (la+2)+(lb+2)+1 = la+lb+5.  For two g-shells that's 4+4+5 = 13.
const MAX_EC_LEN: usize = 16;

/// Fixed-size E-coefficient entry to avoid heap allocation.
#[derive(Clone, Copy)]
struct Ec {
    data: [f64; MAX_EC_LEN],
    len: usize,
    valid: bool,
}

impl Default for Ec {
    fn default() -> Self {
        Self { data: [0.0; MAX_EC_LEN], len: 0, valid: false }
    }
}

impl Ec {
    fn from_vec(v: &[f64]) -> Self {
        let mut ec = Self::default();
        let n = v.len().min(MAX_EC_LEN);
        ec.data[..n].copy_from_slice(&v[..n]);
        ec.len = n;
        ec.valid = true;
        ec
    }

    #[inline(always)]
    fn get(&self, i: usize) -> f64 {
        if i < self.len { self.data[i] } else { 0.0 }
    }
}

/// Precompute primitive norms for each angular component × primitive in a shell group.
fn precompute_prim_norms(
    ang: &[(i32, i32, i32)],
    prims: &[PrimShell],
) -> Vec<Vec<f64>> {
    ang.iter()
        .map(|&(lx, ly, lz)| {
            prims.iter().map(|ps| primitive_norm(ps.exponent, lx, ly, lz)).collect()
        })
        .collect()
}

/// Compute 2-electron skeleton Hessian.
pub fn compute_2e_hessian(
    shells: &[PrimShell],
    density: &[f64],
    norms: &[f64],
    n: usize,
    n_atom: usize,
) -> Vec<f64> {
    let ndim = 3 * n_atom;
    let mut hess = vec![0.0; ndim * ndim];
    let groups = group_shells(shells);
    let n_grp = groups.len();

    // Reusable buffers — allocated ONCE
    let mut rbuf = RBuffer::new();
    let mut boys_buf = vec![0.0; 64]; // grown as needed

    // Precompute all primitive norms per shell group
    let all_pnorms: Vec<Vec<Vec<f64>>> = groups.iter().map(|g| {
        precompute_prim_norms(angular_momentums(g.shell_type), &g.primitives)
    }).collect();

    for ig_a in 0..n_grp {
        let g_a = &groups[ig_a];
        for ig_b in 0..n_grp {
            let g_b = &groups[ig_b];
            for ig_c in 0..n_grp {
                let g_c = &groups[ig_c];
                for ig_d in 0..n_grp {
                    let g_d = &groups[ig_d];

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

                                    let d_ab = density[mu * n + nu];
                                    let d_cd = density[lam * n + sig];
                                    let d_ac = density[mu * n + lam];
                                    let d_bd = density[nu * n + sig];
                                    let d_ad = density[mu * n + sig];
                                    let d_bc = density[nu * n + lam];
                                    let density_w = 0.5 * d_ab * d_cd - 0.125 * (d_ac * d_bd + d_ad * d_bc);
                                    if density_w.abs() < 1e-15 { continue; }

                                    let la = [[lxa, lya, lza], [lxb, lyb, lzb], [lxc, lyc, lzc], [lxd, lyd, lzd]];
                                    let norm_mu_nu_lam_sig = norms[mu] * norms[nu] * norms[lam] * norms[sig];

                                    for (ip_a, p_a) in g_a.primitives.iter().enumerate() {
                                        let n_a = all_pnorms[ig_a][ia][ip_a];
                                        for (ip_b, p_b) in g_b.primitives.iter().enumerate() {
                                            let n_b = all_pnorms[ig_b][ib][ip_b];
                                            let alpha = p_a.exponent;
                                            let beta = p_b.exponent;
                                            let p_val = alpha + beta;
                                            let aa = [p_a.x, p_a.y, p_a.z];
                                            let bb = [p_b.x, p_b.y, p_b.z];
                                            let kab = (-(alpha * beta / p_val)
                                                * ((aa[0]-bb[0]).powi(2) + (aa[1]-bb[1]).powi(2) + (aa[2]-bb[2]).powi(2))).exp();
                                            let px = (alpha*aa[0] + beta*bb[0]) / p_val;
                                            let py = (alpha*aa[1] + beta*bb[1]) / p_val;
                                            let pz = (alpha*aa[2] + beta*bb[2]) / p_val;
                                            let pa = [px-aa[0], py-aa[1], pz-aa[2]];
                                            let pb = [px-bb[0], py-bb[1], pz-bb[2]];
                                            let coeff_ab = p_a.coefficient * p_b.coefficient * n_a * n_b * kab;

                                            for (ip_c, p_c) in g_c.primitives.iter().enumerate() {
                                                let n_c = all_pnorms[ig_c][ic][ip_c];
                                                for (ip_d, p_d) in g_d.primitives.iter().enumerate() {
                                                    let n_d = all_pnorms[ig_d][id][ip_d];
                                                    let gam = p_c.exponent;
                                                    let del = p_d.exponent;
                                                    let q_val = gam + del;
                                                    let cc = [p_c.x, p_c.y, p_c.z];
                                                    let dd = [p_d.x, p_d.y, p_d.z];
                                                    let kcd = (-(gam * del / q_val)
                                                        * ((cc[0]-dd[0]).powi(2) + (cc[1]-dd[1]).powi(2) + (cc[2]-dd[2]).powi(2))).exp();

                                                    let coeff = coeff_ab * p_c.coefficient * p_d.coefficient
                                                        * n_c * n_d * norm_mu_nu_lam_sig;
                                                    let w = density_w * coeff;
                                                    if w.abs() < 1e-18 { continue; }

                                                    let exps = [alpha, beta, gam, del];
                                                    let atoms4 = [
                                                        p_a.atom_index as usize,
                                                        p_b.atom_index as usize,
                                                        p_c.atom_index as usize,
                                                        p_d.atom_index as usize,
                                                    ];

                                                    let eta = p_val * q_val / (p_val + q_val);
                                                    let qx = (gam*cc[0] + del*dd[0]) / q_val;
                                                    let qy = (gam*cc[1] + del*dd[1]) / q_val;
                                                    let qz = (gam*cc[2] + del*dd[2]) / q_val;
                                                    let rpq = [px-qx, py-qy, pz-qz];
                                                    let rpq2 = rpq[0]*rpq[0] + rpq[1]*rpq[1] + rpq[2]*rpq[2];

                                                    let pre2e = 2.0 * std::f64::consts::PI.powf(2.5)
                                                        / (p_val * q_val * (p_val + q_val).sqrt())
                                                        * kab * kcd;

                                                    let qc = [qx-cc[0], qy-cc[1], qz-cc[2]];
                                                    let qd = [qx-dd[0], qy-dd[1], qz-dd[2]];

                                                    let max_n_orig = (lxa+lya+lza+lxb+lyb+lzb+lxc+lyc+lzc+lxd+lyd+lzd) as usize;
                                                    let max_n = max_n_orig + 2;

                                                    // Reuse boys_buf (grow if needed)
                                                    if boys_buf.len() < max_n + 2 {
                                                        boys_buf.resize(max_n + 2, 0.0);
                                                    }
                                                    boys_all(max_n, eta * rpq2, &mut boys_buf);
                                                    rbuf.compute(max_n, eta, rpq[0], rpq[1], rpq[2], &boys_buf);

                                                    // E-coefficient cache on stack: [dir][s1+2][s2+2]
                                                    let mut bra_ec = [[[Ec::default(); 5]; 5]; 3];
                                                    let mut ket_ec = [[[Ec::default(); 5]; 5]; 3];
                                                    for d in 0..3 {
                                                        for s1 in -2i32..=2 {
                                                            for s2 in -2i32..=2 {
                                                                let l1b = la[0][d] + s1;
                                                                let l2b = la[1][d] + s2;
                                                                if l1b >= 0 && l2b >= 0 {
                                                                    let v = e_coefficients(l1b, l2b, p_val, pa[d], pb[d]);
                                                                    bra_ec[d][(s1+2) as usize][(s2+2) as usize] = Ec::from_vec(&v);
                                                                }
                                                                let l1k = la[2][d] + s1;
                                                                let l2k = la[3][d] + s2;
                                                                if l1k >= 0 && l2k >= 0 {
                                                                    let v = e_coefficients(l1k, l2k, q_val, qc[d], qd[d]);
                                                                    ket_ec[d][(s1+2) as usize][(s2+2) as usize] = Ec::from_vec(&v);
                                                                }
                                                            }
                                                        }
                                                    }

                                                    // Inline ERI with AM shifts
                                                    let eri3d = |dl: [[i32; 3]; 4]| -> f64 {
                                                        for d in 0..3 {
                                                            if la[0][d]+dl[0][d] < 0 || la[1][d]+dl[1][d] < 0
                                                                || la[2][d]+dl[2][d] < 0 || la[3][d]+dl[3][d] < 0 {
                                                                return 0.0;
                                                            }
                                                        }
                                                        let ex = &bra_ec[0][(dl[0][0]+2) as usize][(dl[1][0]+2) as usize];
                                                        let ey = &bra_ec[1][(dl[0][1]+2) as usize][(dl[1][1]+2) as usize];
                                                        let ez = &bra_ec[2][(dl[0][2]+2) as usize][(dl[1][2]+2) as usize];
                                                        let fx = &ket_ec[0][(dl[2][0]+2) as usize][(dl[3][0]+2) as usize];
                                                        let fy = &ket_ec[1][(dl[2][1]+2) as usize][(dl[3][1]+2) as usize];
                                                        let fz = &ket_ec[2][(dl[2][2]+2) as usize][(dl[3][2]+2) as usize];
                                                        if !ex.valid || !ey.valid || !ez.valid || !fx.valid || !fy.valid || !fz.valid {
                                                            return 0.0;
                                                        }
                                                        let max_bx = (la[0][0]+dl[0][0]+la[1][0]+dl[1][0]) as usize;
                                                        let max_by = (la[0][1]+dl[0][1]+la[1][1]+dl[1][1]) as usize;
                                                        let max_bz = (la[0][2]+dl[0][2]+la[1][2]+dl[1][2]) as usize;
                                                        let max_kx = (la[2][0]+dl[2][0]+la[3][0]+dl[3][0]) as usize;
                                                        let max_ky = (la[2][1]+dl[2][1]+la[3][1]+dl[3][1]) as usize;
                                                        let max_kz = (la[2][2]+dl[2][2]+la[3][2]+dl[3][2]) as usize;

                                                        if max_bx+max_by+max_bz+max_kx+max_ky+max_kz == 0 {
                                                            return ex.get(0)*ey.get(0)*ez.get(0)*fx.get(0)*fy.get(0)*fz.get(0)
                                                                * rbuf.get(0,0,0,0) * pre2e;
                                                        }

                                                        let mut sum = 0.0;
                                                        for t1 in 0..=max_bx {
                                                            for u1 in 0..=max_by {
                                                                for v1 in 0..=max_bz {
                                                                    let e1 = ex.get(t1)*ey.get(u1)*ez.get(v1);
                                                                    if e1.abs() < 1e-18 { continue; }
                                                                    for t2 in 0..=max_kx {
                                                                        for u2 in 0..=max_ky {
                                                                            for v2 in 0..=max_kz {
                                                                                let sign = if (t2+u2+v2) & 1 != 0 { -1.0 } else { 1.0 };
                                                                                sum += e1 * fx.get(t2)*fy.get(u2)*fz.get(v2)
                                                                                    * sign * rbuf.get(t1+t2, u1+u2, v1+v2, 0);
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                        sum * pre2e
                                                    };

                                                    let eri0 = eri3d([[0;3];4]);

                                                    // d²(μν|λσ)/dR for all center pairs
                                                    for c1 in 0..4usize {
                                                        for c2 in c1..4usize {
                                                            for d1 in 0..3usize {
                                                                let d2_start = if c1 == c2 { d1 } else { 0 };
                                                                for d2 in d2_start..3usize {
                                                                    let val;
                                                                    if c1 == c2 {
                                                                        let e = exps[c1];
                                                                        let l_val = la[c1][d1];
                                                                        if d1 == d2 {
                                                                            let mut sh = [[0i32;3];4]; sh[c1][d1] = 2;
                                                                            let mut sm = [[0i32;3];4]; sm[c1][d1] = -2;
                                                                            val = 4.0*e*e*eri3d(sh)
                                                                                - 2.0*e*(2*l_val+1) as f64 * eri0
                                                                                + if l_val >= 2 { (l_val*(l_val-1)) as f64 * eri3d(sm) } else { 0.0 };
                                                                        } else {
                                                                            let l2v = la[c1][d2];
                                                                            let mut sh_pp = [[0i32;3];4]; sh_pp[c1][d1]=1; sh_pp[c1][d2]+=1;
                                                                            let mut sh_mp = [[0i32;3];4]; sh_mp[c1][d1]=-1; sh_mp[c1][d2]+=1;
                                                                            let mut sh_pm = [[0i32;3];4]; sh_pm[c1][d1]=1; sh_pm[c1][d2]-=1;
                                                                            let mut sh_mm = [[0i32;3];4]; sh_mm[c1][d1]=-1; sh_mm[c1][d2]-=1;
                                                                            val = 4.0*e*e*eri3d(sh_pp)
                                                                                - if l_val > 0 { 2.0*e*l_val as f64*eri3d(sh_mp) } else { 0.0 }
                                                                                - if l2v > 0 { 2.0*e*l2v as f64*eri3d(sh_pm) } else { 0.0 }
                                                                                + if l_val > 0 && l2v > 0 { (l_val*l2v) as f64*eri3d(sh_mm) } else { 0.0 };
                                                                        }
                                                                    } else {
                                                                        let e1 = exps[c1]; let e2 = exps[c2];
                                                                        let l1v = la[c1][d1]; let l2v = la[c2][d2];
                                                                        let mut sh_pp = [[0i32;3];4]; sh_pp[c1][d1]=1; sh_pp[c2][d2]=1;
                                                                        let mut sh_mp = [[0i32;3];4]; sh_mp[c1][d1]=-1; sh_mp[c2][d2]=1;
                                                                        let mut sh_pm = [[0i32;3];4]; sh_pm[c1][d1]=1; sh_pm[c2][d2]=-1;
                                                                        let mut sh_mm = [[0i32;3];4]; sh_mm[c1][d1]=-1; sh_mm[c2][d2]=-1;
                                                                        val = 4.0*e1*e2*eri3d(sh_pp)
                                                                            - if l1v > 0 { 2.0*e2*l1v as f64*eri3d(sh_mp) } else { 0.0 }
                                                                            - if l2v > 0 { 2.0*e1*l2v as f64*eri3d(sh_pm) } else { 0.0 }
                                                                            + if l1v > 0 && l2v > 0 { (l1v*l2v) as f64*eri3d(sh_mm) } else { 0.0 };
                                                                    }

                                                                    if val.abs() < 1e-18 { continue; }
                                                                    let val_w = val * w;
                                                                    let g1 = 3*atoms4[c1]+d1;
                                                                    let g2 = 3*atoms4[c2]+d2;
                                                                    let l1_local = c1*3+d1;
                                                                    let l2_local = c2*3+d2;

                                                                    if g1 != g2 {
                                                                        hess[g1*ndim+g2] += val_w;
                                                                        hess[g2*ndim+g1] += val_w;
                                                                    } else if l1_local == l2_local {
                                                                        hess[g1*ndim+g2] += val_w;
                                                                    } else {
                                                                        hess[g1*ndim+g2] += 2.0*val_w;
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
                }
            }
        }
    }
    hess
}
