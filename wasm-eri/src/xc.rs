/// XC integration: evaluate basis functions on grid, compute density, functional, and V_xc matrix
///
/// Functional IDs: 0=SVWN, 1=BLYP, 2=PBE, 3=B3LYP

use crate::integrals::{angular_momentums, shell_type_to_num_basis, primitive_norm, PrimShell};
use crate::simd_utils;

const PI: f64 = std::f64::consts::PI;

// ─── Integer power for small exponents ───
#[inline(always)]
fn ipow(x: f64, n: i32) -> f64 {
    match n {
        0 => 1.0,
        1 => x,
        2 => x * x,
        3 => x * x * x,
        _ => x.powi(n),
    }
}

// ─── Shell norm cache ───
struct ShellCache {
    p_norms: Vec<f64>,
    r2_max: f64,
}

fn build_shell_caches(shells: &[PrimShell]) -> Vec<ShellCache> {
    let log_thresh: f64 = 69.08; // -ln(1e-30)
    shells
        .iter()
        .map(|s| {
            let ang = angular_momentums(s.shell_type);
            let n_bf = shell_type_to_num_basis(s.shell_type);
            let mut p_norms = Vec::with_capacity(n_bf);
            for bf in 0..n_bf {
                let (lx, ly, lz) = ang[bf];
                p_norms.push(primitive_norm(s.exponent, lx, ly, lz));
            }
            ShellCache {
                p_norms,
                r2_max: log_thresh / s.exponent,
            }
        })
        .collect()
}

// ─── Basis function evaluation with screening ───
struct BasisEval {
    phi: Vec<f64>,
    dphi_dx: Vec<f64>,
    dphi_dy: Vec<f64>,
    dphi_dz: Vec<f64>,
    sig_idx: Vec<usize>,
}

fn evaluate_basis_screened(
    shells: &[PrimShell],
    caches: &[ShellCache],
    norm_factors: &[f64],
    nbasis: usize,
    px: f64, py: f64, pz: f64,
    need_grad: bool,
) -> BasisEval {
    let mut phi = vec![0.0f64; nbasis];
    let mut dphi_dx = vec![0.0f64; nbasis];
    let mut dphi_dy = vec![0.0f64; nbasis];
    let mut dphi_dz = vec![0.0f64; nbasis];

    for (s, shell) in shells.iter().enumerate() {
        let cache = &caches[s];
        let dx = px - shell.x;
        let dy = py - shell.y;
        let dz = pz - shell.z;
        let r2 = dx * dx + dy * dy + dz * dz;
        if r2 > cache.r2_max {
            continue;
        }
        let exp_val = shell.coefficient * (-shell.exponent * r2).exp();

        let ang = angular_momentums(shell.shell_type);
        let n_bf = shell_type_to_num_basis(shell.shell_type);

        for bf in 0..n_bf {
            let (lx, ly, lz) = ang[bf];
            let basis_idx = shell.basis_index as usize + bf;
            let coeff = norm_factors[basis_idx] * cache.p_norms[bf] * exp_val;

            let angular = ipow(dx, lx) * ipow(dy, ly) * ipow(dz, lz);
            phi[basis_idx] += coeff * angular;

            if need_grad {
                let exp_alpha = -2.0 * shell.exponent;
                let pxl = ipow(dx, lx);
                let pyl = ipow(dy, ly);
                let pzl = ipow(dz, lz);

                let mut d_ang_dx = exp_alpha * dx * angular;
                if lx > 0 {
                    d_ang_dx += lx as f64 * ipow(dx, lx - 1) * pyl * pzl;
                }
                let mut d_ang_dy = exp_alpha * dy * angular;
                if ly > 0 {
                    d_ang_dy += ly as f64 * pxl * ipow(dy, ly - 1) * pzl;
                }
                let mut d_ang_dz = exp_alpha * dz * angular;
                if lz > 0 {
                    d_ang_dz += lz as f64 * pxl * pyl * ipow(dz, lz - 1);
                }

                dphi_dx[basis_idx] += coeff * d_ang_dx;
                dphi_dy[basis_idx] += coeff * d_ang_dy;
                dphi_dz[basis_idx] += coeff * d_ang_dz;
            }
        }
    }

    // Build significant index list
    let thresh = 1e-15f64;
    let mut sig_idx = Vec::new();
    for i in 0..nbasis {
        if phi[i].abs() > thresh
            || (need_grad
                && (dphi_dx[i].abs() > thresh
                    || dphi_dy[i].abs() > thresh
                    || dphi_dz[i].abs() > thresh))
        {
            sig_idx.push(i);
        }
    }

    BasisEval {
        phi,
        dphi_dx,
        dphi_dy,
        dphi_dz,
        sig_idx,
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// XC Functionals
// ═══════════════════════════════════════════════════════════════════════════

struct XCOutput {
    exc: f64,
    vrho_a: f64,
    vrho_b: f64,
    vgamma_aa: f64,
    vgamma_bb: f64,
    vgamma_ab: f64,
}

// ─── Slater exchange ───
const CX: f64 = 0.738_558_766_382_022_4; // 0.75 * (3/π)^(1/3)
const TWO_1_3: f64 = 1.259_921_049_894_873_2; // 2^(1/3)

fn slater_exchange(rho_a: f64, rho_b: f64) -> (f64, f64, f64) {
    let rho = rho_a + rho_b;
    if rho < 1e-20 {
        return (0.0, 0.0, 0.0);
    }
    let fa = if rho_a > 1e-20 { rho_a.powf(1.0 / 3.0) } else { 0.0 };
    let fb = if rho_b > 1e-20 { rho_b.powf(1.0 / 3.0) } else { 0.0 };
    let exc = -TWO_1_3 * CX * (rho_a * fa + rho_b * fb) / rho;
    let vra = if rho_a > 1e-20 { -(4.0 / 3.0) * TWO_1_3 * CX * fa } else { 0.0 };
    let vrb = if rho_b > 1e-20 { -(4.0 / 3.0) * TWO_1_3 * CX * fb } else { 0.0 };
    (exc, vra, vrb)
}

// ─── VWN5 correlation ───
fn vwn5_correlation(rho_a: f64, rho_b: f64) -> (f64, f64, f64) {
    let rho = rho_a + rho_b;
    if rho < 1e-20 {
        return (0.0, 0.0, 0.0);
    }

    let rs = (3.0 / (4.0 * PI * rho)).powf(1.0 / 3.0);
    let zeta = (rho_a - rho_b) / rho;

    let a0 = 0.0310907f64; let b0 = 3.72744f64; let c0 = 12.9352f64; let x00 = -0.10498f64;
    let a1 = 0.01554535f64; let b1 = 7.06042f64; let c1 = 18.0578f64; let x01 = -0.32500f64;
    let aa = -1.0 / (6.0 * PI * PI); let ba = 1.13107f64; let ca = 13.0045f64; let x0a = -0.0047584f64;

    let vwn_eps = |a: f64, b: f64, c: f64, x0: f64, x: f64| -> (f64, f64) {
        let xc = x * x + b * x + c;
        let x0c = x0 * x0 + b * x0 + c;
        let q = (4.0 * c - b * b).sqrt();
        let eps = a * (
            (x * x / xc).ln()
            + 2.0 * b / q * (q / (2.0 * x + b)).atan()
            - b * x0 / x0c * (
                ((x - x0) * (x - x0) / xc).ln()
                + 2.0 * (b + 2.0 * x0) / q * (q / (2.0 * x + b)).atan()
            )
        );
        let dxdx = 2.0 * x + b;
        let denom1 = (2.0 * x + b).powi(2) + q * q;
        let depsdx = a * (
            2.0 / x - dxdx / xc
            - 4.0 * b / (q * denom1)
            - b * x0 / x0c * (
                2.0 / (x - x0) - dxdx / xc
                - 4.0 * (b + 2.0 * x0) / (q * denom1)
            )
        );
        let deps = depsdx / (2.0 * x); // dε/drs (x = sqrt(rs))
        (eps, deps)
    };

    let sqrt_rs = rs.sqrt();
    let (ec0, dec0) = vwn_eps(a0, b0, c0, x00, sqrt_rs);
    let (ec1, dec1) = vwn_eps(a1, b1, c1, x01, sqrt_rs);
    let (ac, dac) = vwn_eps(aa, ba, ca, x0a, sqrt_rs);

    let fzeta_denom = 2.0 * (2.0f64.powf(1.0 / 3.0) - 1.0);
    let f_val = (if (1.0 + zeta) > 1e-20 { (1.0 + zeta).powf(4.0 / 3.0) } else { 0.0 })
        + (if (1.0 - zeta) > 1e-20 { (1.0 - zeta).powf(4.0 / 3.0) } else { 0.0 })
        - 2.0;
    let fz = f_val / fzeta_denom;
    let fpp0 = 4.0 / (9.0 * fzeta_denom);

    let z4 = zeta * zeta * zeta * zeta;
    let exc = ec0 + ac * fz / fpp0 * (1.0 - z4) + (ec1 - ec0) * fz * z4;

    let dec_drs = dec0 + dac * fz / fpp0 * (1.0 - z4) + (dec1 - dec0) * fz * z4;

    let dfz_a = if (1.0 + zeta) > 1e-20 { (4.0 / 3.0) * (1.0 + zeta).powf(1.0 / 3.0) } else { 0.0 };
    let dfz_b = if (1.0 - zeta) > 1e-20 { -(4.0 / 3.0) * (1.0 - zeta).powf(1.0 / 3.0) } else { 0.0 };
    let dfzdz = (dfz_a + dfz_b) / fzeta_denom;

    let dexc_dzeta = dfzdz * (ac / fpp0 * (1.0 - z4) + (ec1 - ec0) * z4)
        + fz * (-4.0 * ac / fpp0 * zeta.powi(3) + 4.0 * (ec1 - ec0) * zeta.powi(3));

    let vc_common = exc - rs / 3.0 * dec_drs;
    let vrho_a = if rho > 1e-20 { vc_common + dexc_dzeta * (1.0 - zeta) } else { 0.0 };
    let vrho_b = if rho > 1e-20 { vc_common - dexc_dzeta * (1.0 + zeta) } else { 0.0 };

    (exc, vrho_a, vrho_b)
}

// ─── B88 exchange ───
fn b88_exchange(rho_a: f64, rho_b: f64, gamma_aa: f64, gamma_bb: f64) -> XCOutput {
    let (s_exc, s_vra, s_vrb) = slater_exchange(rho_a, rho_b);
    let beta = 0.0042f64;

    let eval_b88 = |rs: f64, gs: f64| -> f64 {
        if rs < 1e-20 { return 0.0; }
        let r43 = rs.powf(4.0 / 3.0);
        let x = gs.max(0.0).sqrt() / r43;
        -beta * r43 * x * x / (1.0 + 6.0 * beta * x * x.asinh())
    };

    let exc_gga = eval_b88(rho_a, gamma_aa) + eval_b88(rho_b, gamma_bb);

    let h = 1e-6f64;
    let e0a = eval_b88(rho_a, gamma_aa);
    let e0b = eval_b88(rho_b, gamma_bb);
    let vr_a = if rho_a > 1e-20 { (eval_b88(rho_a + h, gamma_aa) - eval_b88(rho_a - h, gamma_aa)) / (2.0 * h) } else { 0.0 };
    let vr_b = if rho_b > 1e-20 { (eval_b88(rho_b + h, gamma_bb) - eval_b88(rho_b - h, gamma_bb)) / (2.0 * h) } else { 0.0 };
    let vg_aa = if rho_a > 1e-20 { (eval_b88(rho_a, gamma_aa + h) - e0a) / h } else { 0.0 };
    let vg_bb = if rho_b > 1e-20 { (eval_b88(rho_b, gamma_bb + h) - e0b) / h } else { 0.0 };

    let rho = rho_a + rho_b;
    XCOutput {
        exc: if rho > 1e-20 { (s_exc * rho + exc_gga) / rho } else { 0.0 },
        vrho_a: s_vra + vr_a,
        vrho_b: s_vrb + vr_b,
        vgamma_aa: vg_aa,
        vgamma_bb: vg_bb,
        vgamma_ab: 0.0,
    }
}

// ─── LYP correlation ───
fn lyp_correlation(rho_a: f64, rho_b: f64, gamma_aa: f64, gamma_bb: f64, gamma_ab: f64) -> XCOutput {
    let rho = rho_a + rho_b;
    if rho < 1e-20 {
        return XCOutput { exc: 0.0, vrho_a: 0.0, vrho_b: 0.0, vgamma_aa: 0.0, vgamma_bb: 0.0, vgamma_ab: 0.0 };
    }

    let a = 0.04918f64; let b = 0.132f64; let c = 0.2533f64; let d = 0.349f64;
    let cf = 0.3 * (3.0 * PI * PI).powf(2.0 / 3.0);
    let cf_scaled = 2.0f64.powf(11.0 / 3.0) * cf;

    let eval_lyp = |ra: f64, rb: f64, gaa: f64, gbb: f64, gab: f64| -> f64 {
        let r = ra + rb;
        if r < 1e-20 || ra < 0.0 || rb < 0.0 { return 0.0; }
        let r13 = r.powf(-1.0 / 3.0);
        let om = (-c * r13).exp() / (1.0 + d * r13) * r.powf(-11.0 / 3.0);
        let dl = c * r13 + d * r13 / (1.0 + d * r13);
        let g = gaa + gbb + 2.0 * gab;
        let t1 = -a * 4.0 * ra * rb / (r * (1.0 + d * r13));
        let t2 = -a * b * om * ra * rb * (
            (47.0 / 18.0 - 7.0 * dl / 18.0) * g
            + (-5.0 / 2.0 + dl / 18.0) * (gaa + gbb)
            + (11.0 - dl) / 9.0 * (ra * gaa / r + rb * gbb / r)
        );
        let t3 = -a * b * om * (
            -2.0 / 3.0 * r * r * g
            + (2.0 / 3.0 * r * r - ra * ra) * gbb
            + (2.0 / 3.0 * r * r - rb * rb) * gaa
        );
        let t4 = -a * b * om * ra * rb * cf_scaled * (ra.powf(8.0 / 3.0) + rb.powf(8.0 / 3.0));
        t1 + t2 + t3 + t4
    };

    let e0 = eval_lyp(rho_a, rho_b, gamma_aa, gamma_bb, gamma_ab);
    let exc_val = e0 / rho;

    let h = 1e-6f64;
    let vr_a = if rho_a > h {
        (eval_lyp(rho_a + h, rho_b, gamma_aa, gamma_bb, gamma_ab)
         - eval_lyp(rho_a - h, rho_b, gamma_aa, gamma_bb, gamma_ab)) / (2.0 * h)
    } else {
        (eval_lyp(rho_a + h, rho_b, gamma_aa, gamma_bb, gamma_ab) - e0) / h
    };
    let vr_b = if rho_b > h {
        (eval_lyp(rho_a, rho_b + h, gamma_aa, gamma_bb, gamma_ab)
         - eval_lyp(rho_a, rho_b - h, gamma_aa, gamma_bb, gamma_ab)) / (2.0 * h)
    } else {
        (eval_lyp(rho_a, rho_b + h, gamma_aa, gamma_bb, gamma_ab) - e0) / h
    };
    let vg_aa = (eval_lyp(rho_a, rho_b, gamma_aa + h, gamma_bb, gamma_ab) - e0) / h;
    let vg_bb = (eval_lyp(rho_a, rho_b, gamma_aa, gamma_bb + h, gamma_ab) - e0) / h;
    let vg_ab = (eval_lyp(rho_a, rho_b, gamma_aa, gamma_bb, gamma_ab + h) - e0) / h;

    XCOutput { exc: exc_val, vrho_a: vr_a, vrho_b: vr_b, vgamma_aa: vg_aa, vgamma_bb: vg_bb, vgamma_ab: vg_ab }
}

// ─── PBE exchange ───
fn pbe_exchange(rho_a: f64, rho_b: f64, gamma_aa: f64, gamma_bb: f64) -> XCOutput {
    let (s_exc, s_vra, s_vrb) = slater_exchange(rho_a, rho_b);
    let kappa = 0.8040f64;
    let mu = 0.2195149727645171f64;

    let eval_pbe_s = |r: f64, g: f64| -> f64 {
        if r < 1e-20 { return 0.0; }
        let k = (6.0 * PI * PI * r).powf(1.0 / 3.0);
        let ss = g / (4.0 * r * r * k * k);
        let f = 1.0 + kappa - kappa / (1.0 + mu * ss / kappa);
        r * (-TWO_1_3 * CX * r.powf(1.0 / 3.0)) * (f - 1.0)
    };

    let mut exc_gga = 0.0;
    let mut vr_a = 0.0; let mut vr_b = 0.0;
    let mut vg_aa = 0.0; let mut vg_bb = 0.0;
    let h = 1e-6f64;

    for &(rho_s, gamma_ss, is_a) in &[(rho_a, gamma_aa, true), (rho_b, gamma_bb, false)] {
        if rho_s < 1e-20 { continue; }
        exc_gga += eval_pbe_s(rho_s, gamma_ss);
        let vr = (eval_pbe_s(rho_s + h, gamma_ss) - eval_pbe_s(rho_s - h, gamma_ss)) / (2.0 * h);
        let vg = (eval_pbe_s(rho_s, gamma_ss + h) - eval_pbe_s(rho_s, gamma_ss - h)) / (2.0 * h);
        if is_a { vr_a = vr; vg_aa = vg; } else { vr_b = vr; vg_bb = vg; }
    }

    let rho = rho_a + rho_b;
    XCOutput {
        exc: if rho > 1e-20 { (s_exc * rho + exc_gga) / rho } else { 0.0 },
        vrho_a: s_vra + vr_a,
        vrho_b: s_vrb + vr_b,
        vgamma_aa: vg_aa,
        vgamma_bb: vg_bb,
        vgamma_ab: 0.0,
    }
}

// ─── PBE correlation ───
fn pbe_correlation(rho_a: f64, rho_b: f64, gamma_aa: f64, gamma_bb: f64, gamma_ab: f64) -> XCOutput {
    let (vwn_exc, vwn_vra, vwn_vrb) = vwn5_correlation(rho_a, rho_b);
    let rho = rho_a + rho_b;
    if rho < 1e-20 {
        return XCOutput { exc: vwn_exc, vrho_a: vwn_vra, vrho_b: vwn_vrb, vgamma_aa: 0.0, vgamma_bb: 0.0, vgamma_ab: 0.0 };
    }

    let beta_pbe = 0.06672455060314922f64;
    let gamma_pbe = 0.031090690869654895f64;

    let eval_h = |ra: f64, rb: f64, gaa: f64, gbb: f64, gab: f64| -> f64 {
        let r = ra + rb;
        if r < 1e-20 || ra < 0.0 || rb < 0.0 { return 0.0; }
        let vwn_e = vwn5_correlation(ra, rb).0;
        let z = (ra - rb) / r;
        let ph = 0.5 * (if (1.0 + z) > 1e-20 { (1.0 + z).powf(2.0 / 3.0) } else { 0.0 })
               + 0.5 * (if (1.0 - z) > 1e-20 { (1.0 - z).powf(2.0 / 3.0) } else { 0.0 });
        let ph3 = ph * ph * ph;
        let kf = (3.0 * PI * PI * r).powf(1.0 / 3.0);
        let ks = (4.0 * kf / PI).sqrt();
        let t2 = (gaa + gbb + 2.0 * gab) / (4.0 * ph * ph * ks * ks * r * r);
        let exp_val = (-vwn_e / gamma_pbe).exp();
        let a = beta_pbe / gamma_pbe / (if (exp_val - 1.0).abs() > 1e-30 { exp_val - 1.0 } else { 1e-30 });
        let at2 = a * t2;
        let h_val = gamma_pbe * ph3 * (1.0 + beta_pbe / gamma_pbe * t2 * (1.0 + at2) / (1.0 + at2 + at2 * at2)).ln();
        r * h_val // return rho*H
    };

    let h0 = eval_h(rho_a, rho_b, gamma_aa, gamma_bb, gamma_ab);
    let h_per_particle = if rho > 1e-20 { h0 / rho } else { 0.0 };

    let h = 1e-6f64;
    let vr_a = vwn_vra + if rho_a > h {
        (eval_h(rho_a + h, rho_b, gamma_aa, gamma_bb, gamma_ab) - eval_h(rho_a - h, rho_b, gamma_aa, gamma_bb, gamma_ab)) / (2.0 * h)
    } else {
        (eval_h(rho_a + h, rho_b, gamma_aa, gamma_bb, gamma_ab) - h0) / h
    };
    let vr_b = vwn_vrb + if rho_b > h {
        (eval_h(rho_a, rho_b + h, gamma_aa, gamma_bb, gamma_ab) - eval_h(rho_a, rho_b - h, gamma_aa, gamma_bb, gamma_ab)) / (2.0 * h)
    } else {
        (eval_h(rho_a, rho_b + h, gamma_aa, gamma_bb, gamma_ab) - h0) / h
    };
    let vg_aa = (eval_h(rho_a, rho_b, gamma_aa + h, gamma_bb, gamma_ab) - h0) / h;
    let vg_bb = (eval_h(rho_a, rho_b, gamma_aa, gamma_bb + h, gamma_ab) - h0) / h;
    let vg_ab = (eval_h(rho_a, rho_b, gamma_aa, gamma_bb, gamma_ab + h) - h0) / h;

    XCOutput {
        exc: vwn_exc + h_per_particle,
        vrho_a: vr_a, vrho_b: vr_b,
        vgamma_aa: vg_aa, vgamma_bb: vg_bb, vgamma_ab: vg_ab,
    }
}

// ─── Functional dispatch ───
// func_id: 0=SVWN, 1=BLYP, 2=PBE, 3=B3LYP
fn evaluate_functional(func_id: u32, rho_a: f64, rho_b: f64, gaa: f64, gbb: f64, gab: f64) -> XCOutput {
    match func_id {
        0 => { // SVWN
            let (ex_exc, ex_vra, ex_vrb) = slater_exchange(rho_a, rho_b);
            let (ec_exc, ec_vra, ec_vrb) = vwn5_correlation(rho_a, rho_b);
            XCOutput {
                exc: ex_exc + ec_exc,
                vrho_a: ex_vra + ec_vra,
                vrho_b: ex_vrb + ec_vrb,
                vgamma_aa: 0.0, vgamma_bb: 0.0, vgamma_ab: 0.0,
            }
        }
        1 => { // BLYP
            let ex = b88_exchange(rho_a, rho_b, gaa, gbb);
            let ec = lyp_correlation(rho_a, rho_b, gaa, gbb, gab);
            XCOutput {
                exc: ex.exc + ec.exc,
                vrho_a: ex.vrho_a + ec.vrho_a,
                vrho_b: ex.vrho_b + ec.vrho_b,
                vgamma_aa: ex.vgamma_aa + ec.vgamma_aa,
                vgamma_bb: ex.vgamma_bb + ec.vgamma_bb,
                vgamma_ab: ex.vgamma_ab + ec.vgamma_ab,
            }
        }
        2 => { // PBE
            let ex = pbe_exchange(rho_a, rho_b, gaa, gbb);
            let ec = pbe_correlation(rho_a, rho_b, gaa, gbb, gab);
            XCOutput {
                exc: ex.exc + ec.exc,
                vrho_a: ex.vrho_a + ec.vrho_a,
                vrho_b: ex.vrho_b + ec.vrho_b,
                vgamma_aa: ex.vgamma_aa + ec.vgamma_aa,
                vgamma_bb: ex.vgamma_bb + ec.vgamma_bb,
                vgamma_ab: ex.vgamma_ab + ec.vgamma_ab,
            }
        }
        3 => { // B3LYP
            let a0 = 0.20f64; let ax = 0.72f64; let ac = 0.81f64;
            let (s_exc, s_vra, s_vrb) = slater_exchange(rho_a, rho_b);
            let b88 = b88_exchange(rho_a, rho_b, gaa, gbb);
            let b88corr_exc = b88.exc - s_exc;
            let b88corr_vra = b88.vrho_a - s_vra;
            let b88corr_vrb = b88.vrho_b - s_vrb;
            let (vwn_exc, vwn_vra, vwn_vrb) = vwn5_correlation(rho_a, rho_b);
            let lyp = lyp_correlation(rho_a, rho_b, gaa, gbb, gab);
            XCOutput {
                exc: (1.0 - a0) * s_exc + ax * b88corr_exc + (1.0 - ac) * vwn_exc + ac * lyp.exc,
                vrho_a: (1.0 - a0) * s_vra + ax * b88corr_vra + (1.0 - ac) * vwn_vra + ac * lyp.vrho_a,
                vrho_b: (1.0 - a0) * s_vrb + ax * b88corr_vrb + (1.0 - ac) * vwn_vrb + ac * lyp.vrho_b,
                vgamma_aa: ax * b88.vgamma_aa + ac * lyp.vgamma_aa,
                vgamma_bb: ax * b88.vgamma_bb + ac * lyp.vgamma_bb,
                vgamma_ab: ac * lyp.vgamma_ab,
            }
        }
        _ => XCOutput { exc: 0.0, vrho_a: 0.0, vrho_b: 0.0, vgamma_aa: 0.0, vgamma_bb: 0.0, vgamma_ab: 0.0 },
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Main compute_xc function
// ═══════════════════════════════════════════════════════════════════════════

/// Compute XC energy and V_xc matrix.
///
/// Returns Vec<f64> with layout:
///   [exc, num_electrons, vxc_a (nbasis*nbasis), vxc_b (nbasis*nbasis if unrestricted)]
///
/// grid_flat: [x0, y0, z0, w0, x1, y1, z1, w1, ...]
/// density_b: empty slice for RKS (restricted)
/// func_id: 0=SVWN, 1=BLYP, 2=PBE, 3=B3LYP
/// need_grad: true for GGA/hybrid functionals
pub fn compute_xc(
    shells: &[PrimShell],
    norm_factors: &[f64],
    density_a: &[f64],
    density_b: &[f64],
    grid_flat: &[f64],
    nbasis: usize,
    func_id: u32,
    need_grad: bool,
) -> Vec<f64> {
    let npts = grid_flat.len() / 4;
    let nn = nbasis * nbasis;
    let unrestricted = !density_b.is_empty();

    let caches = build_shell_caches(shells);

    let mut vxc_a = vec![0.0f64; nn];
    let mut vxc_b = if unrestricted { vec![0.0f64; nn] } else { vec![] };
    let mut exc = 0.0f64;
    let mut num_elec = 0.0f64;

    // Reusable buffers for gathered basis function data (avoid per-grid-point allocation)
    let mut phi_sig = Vec::with_capacity(nbasis);
    let mut dphi_dx_sig = Vec::with_capacity(nbasis);
    let mut dphi_dy_sig = Vec::with_capacity(nbasis);
    let mut dphi_dz_sig = Vec::with_capacity(nbasis);
    let mut d_row_a = Vec::with_capacity(nbasis);
    let mut d_row_b = Vec::with_capacity(nbasis);

    for gp in 0..npts {
        let off = gp * 4;
        let px = grid_flat[off];
        let py = grid_flat[off + 1];
        let pz = grid_flat[off + 2];
        let weight = grid_flat[off + 3];

        let eval = evaluate_basis_screened(shells, &caches, norm_factors, nbasis, px, py, pz, need_grad);
        let nsig = eval.sig_idx.len();
        if nsig == 0 { continue; }

        // Gather significant basis function values into contiguous arrays
        phi_sig.clear();
        for ii in 0..nsig {
            phi_sig.push(eval.phi[eval.sig_idx[ii]]);
        }
        if need_grad {
            dphi_dx_sig.clear();
            dphi_dy_sig.clear();
            dphi_dz_sig.clear();
            for ii in 0..nsig {
                let idx = eval.sig_idx[ii];
                dphi_dx_sig.push(eval.dphi_dx[idx]);
                dphi_dy_sig.push(eval.dphi_dy[idx]);
                dphi_dz_sig.push(eval.dphi_dz[idx]);
            }
        }

        // Compute density using SIMD dot products on gathered arrays
        let mut rho_a = 0.0f64;
        let mut rho_b = 0.0f64;
        let mut grad_ax = 0.0f64; let mut grad_ay = 0.0f64; let mut grad_az = 0.0f64;
        let mut grad_bx = 0.0f64; let mut grad_by = 0.0f64; let mut grad_bz = 0.0f64;

        for ii in 0..nsig {
            let mu = eval.sig_idx[ii];
            let phi_mu = phi_sig[ii];
            let mu_off = mu * nbasis;

            // Gather density row for significant indices
            d_row_a.clear();
            for jj in 0..nsig {
                d_row_a.push(density_a[mu_off + eval.sig_idx[jj]]);
            }

            // ρ_α += φ_μ * Σ_ν D_α[μ,ν] * φ_ν  (SIMD dot)
            let dp = simd_utils::dot(&d_row_a, &phi_sig);
            rho_a += phi_mu * dp;

            if unrestricted {
                d_row_b.clear();
                for jj in 0..nsig {
                    d_row_b.push(density_b[mu_off + eval.sig_idx[jj]]);
                }
                rho_b += phi_mu * simd_utils::dot(&d_row_b, &phi_sig);
            }

            if need_grad {
                let dx_mu = dphi_dx_sig[ii];
                let dy_mu = dphi_dy_sig[ii];
                let dz_mu = dphi_dz_sig[ii];
                // ∇ρ_α += Σ_ν D[μ,ν] * (∇φ_μ * φ_ν + φ_μ * ∇φ_ν)
                //       = ∇φ_μ * dot(D_row, φ_sig) + φ_μ * dot(D_row, ∇φ_sig)
                grad_ax += dx_mu * dp + phi_mu * simd_utils::dot(&d_row_a, &dphi_dx_sig);
                grad_ay += dy_mu * dp + phi_mu * simd_utils::dot(&d_row_a, &dphi_dy_sig);
                grad_az += dz_mu * dp + phi_mu * simd_utils::dot(&d_row_a, &dphi_dz_sig);
                if unrestricted {
                    let dp_b = simd_utils::dot(&d_row_b, &phi_sig);
                    grad_bx += dx_mu * dp_b + phi_mu * simd_utils::dot(&d_row_b, &dphi_dx_sig);
                    grad_by += dy_mu * dp_b + phi_mu * simd_utils::dot(&d_row_b, &dphi_dy_sig);
                    grad_bz += dz_mu * dp_b + phi_mu * simd_utils::dot(&d_row_b, &dphi_dz_sig);
                }
            }
        }

        // RKS: split total density into α and β
        if !unrestricted {
            rho_b = rho_a * 0.5;
            rho_a *= 0.5;
            if need_grad {
                grad_bx = grad_ax * 0.5; grad_by = grad_ay * 0.5; grad_bz = grad_az * 0.5;
                grad_ax *= 0.5; grad_ay *= 0.5; grad_az *= 0.5;
            }
        }

        // Clamp to avoid negative densities from numerical noise in density matrix
        if rho_a < 0.0 { rho_a = 0.0; }
        if rho_b < 0.0 { rho_b = 0.0; }
        let rho = rho_a + rho_b;
        if rho < 1e-20 { continue; }

        let gaa = if need_grad { grad_ax * grad_ax + grad_ay * grad_ay + grad_az * grad_az } else { 0.0 };
        let gbb = if need_grad { grad_bx * grad_bx + grad_by * grad_by + grad_bz * grad_bz } else { 0.0 };
        let gab = if need_grad { grad_ax * grad_bx + grad_ay * grad_by + grad_az * grad_bz } else { 0.0 };

        num_elec += rho * weight;
        let xc_out = evaluate_functional(func_id, rho_a, rho_b, gaa, gbb, gab);
        exc += xc_out.exc * rho * weight;

        // GGA pre-compute: fα = 2*vgammaAA*∇ρα + vgammaAB*∇ρβ
        let (f_ax, f_ay, f_az, f_bx, f_by, f_bz) = if need_grad {
            let vg_aa = xc_out.vgamma_aa;
            let vg_ab = xc_out.vgamma_ab;
            let vg_bb = xc_out.vgamma_bb;
            let fax = 2.0 * vg_aa * grad_ax + vg_ab * grad_bx;
            let fay = 2.0 * vg_aa * grad_ay + vg_ab * grad_by;
            let faz = 2.0 * vg_aa * grad_az + vg_ab * grad_bz;
            let (fbx, fby, fbz) = if unrestricted {
                (2.0 * vg_bb * grad_bx + vg_ab * grad_ax,
                 2.0 * vg_bb * grad_by + vg_ab * grad_ay,
                 2.0 * vg_bb * grad_bz + vg_ab * grad_az)
            } else { (0.0, 0.0, 0.0) };
            (fax, fay, faz, fbx, fby, fbz)
        } else {
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        };

        // Build V_xc matrix using gathered arrays
        // For LDA: Vxc[μ,ν] += w * vrho * φ_μ * φ_ν
        // For GGA: also += w * (f · ∇(φ_μ*φ_ν))
        for ii in 0..nsig {
            let mu = eval.sig_idx[ii];
            let phi_mu = phi_sig[ii];

            if !need_grad {
                // LDA path: rank-1 update via daxpy on gathered array
                let scaled_a = xc_out.vrho_a * phi_mu * weight;
                // Accumulate for jj >= ii: vxc[mu, sig_idx[jj]] += scaled_a * phi_sig[jj]
                for jj in ii..nsig {
                    let nu = eval.sig_idx[jj];
                    let wa = scaled_a * phi_sig[jj];
                    vxc_a[mu * nbasis + nu] += wa;
                    if mu != nu { vxc_a[nu * nbasis + mu] += wa; }
                }
                if unrestricted {
                    let scaled_b = xc_out.vrho_b * phi_mu * weight;
                    for jj in ii..nsig {
                        let nu = eval.sig_idx[jj];
                        let wb = scaled_b * phi_sig[jj];
                        vxc_b[mu * nbasis + nu] += wb;
                        if mu != nu { vxc_b[nu * nbasis + mu] += wb; }
                    }
                }
            } else {
                // GGA path
                let dx_mu = dphi_dx_sig[ii];
                let dy_mu = dphi_dy_sig[ii];
                let dz_mu = dphi_dz_sig[ii];

                for jj in ii..nsig {
                    let nu = eval.sig_idx[jj];
                    let phi_nu = phi_sig[jj];

                    let mut val_a = xc_out.vrho_a * phi_mu * phi_nu;
                    let mut val_b = if unrestricted { xc_out.vrho_b * phi_mu * phi_nu } else { 0.0 };

                    let nabla_x = dx_mu * phi_nu + phi_mu * dphi_dx_sig[jj];
                    let nabla_y = dy_mu * phi_nu + phi_mu * dphi_dy_sig[jj];
                    let nabla_z = dz_mu * phi_nu + phi_mu * dphi_dz_sig[jj];
                    val_a += f_ax * nabla_x + f_ay * nabla_y + f_az * nabla_z;
                    if unrestricted {
                        val_b += f_bx * nabla_x + f_by * nabla_y + f_bz * nabla_z;
                    }

                    let wa = val_a * weight;
                    vxc_a[mu * nbasis + nu] += wa;
                    if mu != nu { vxc_a[nu * nbasis + mu] += wa; }

                    if unrestricted {
                        let wb = val_b * weight;
                        vxc_b[mu * nbasis + nu] += wb;
                        if mu != nu { vxc_b[nu * nbasis + mu] += wb; }
                    }
                }
            }
        }
    }

    // Pack result: [exc, num_elec, vxc_a..., vxc_b...]
    let mut result = Vec::with_capacity(2 + nn + if unrestricted { nn } else { 0 });
    result.push(exc);
    result.push(num_elec);
    result.extend_from_slice(&vxc_a);
    if unrestricted {
        result.extend_from_slice(&vxc_b);
    }
    result
}
