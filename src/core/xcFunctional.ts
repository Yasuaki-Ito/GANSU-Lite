/** Exchange-correlation functionals for DFT */

// Input/output for XC evaluation at a single grid point
export interface XCInput {
  rhoA: number;      // alpha density
  rhoB: number;      // beta density
  gammaAA?: number;  // |∇ρα|²  (GGA only)
  gammaBB?: number;  // |∇ρβ|²  (GGA only)
  gammaAB?: number;  // ∇ρα·∇ρβ (GGA only)
  tauA?: number;     // (1/2) Σ_i^occ |∇φ_iα|² (meta-GGA)
  tauB?: number;     // (1/2) Σ_i^occ |∇φ_iβ|² (meta-GGA)
}

export interface XCOutput {
  exc: number;       // energy density ε_xc(r)
  vrhoA: number;     // ∂ε/∂ρα
  vrhoB: number;     // ∂ε/∂ρβ
  vgammaAA?: number; // ∂ε/∂γAA (GGA)
  vgammaBB?: number; // ∂ε/∂γBB (GGA)
  vgammaAB?: number; // ∂ε/∂γAB (GGA)
  vtauA?: number;    // ∂ε/∂τα (meta-GGA, GKS treatment)
  vtauB?: number;    // ∂ε/∂τβ (meta-GGA, GKS treatment)
}

export type FunctionalType = 'LDA' | 'GGA' | 'meta-GGA' | 'hybrid';

/** Range-separation parameters for RSH functionals (e.g. CAM-B3LYP, ωB97X-D).
 *  HF-exchange fraction at distance r:  α + β · erf(ω r) */
export interface RangeSeparation {
  alpha: number;  // HF fraction at r → 0
  beta: number;   // additional HF at r → ∞ (so full HF at infinity is α + β)
  omega: number;  // range-separation parameter (Bohr⁻¹)
}

export interface XCFunctional {
  name: string;
  type: FunctionalType;
  exactExchangeFraction: number; // = α (HF fraction at r=0)
  rangeSeparation?: RangeSeparation; // present for RSH functionals
  evaluate(inp: XCInput): XCOutput;
}

// ── LDA: Slater exchange ──
// ε_x = -C_x * ρ^(1/3), C_x = (3/4)(3/π)^(1/3)
const CX = 0.75 * Math.pow(3 / Math.PI, 1 / 3); // ≈ 0.7386

function slaterExchange(rhoA: number, rhoB: number): XCOutput {
  const rho = rhoA + rhoB;
  if (rho < 1e-20) return { exc: 0, vrhoA: 0, vrhoB: 0 };
  // Spin-polarized: E_x = 2^(1/3) * C_x * (ρα^(4/3) + ρβ^(4/3))
  const fA = rhoA > 1e-20 ? Math.pow(rhoA, 1 / 3) : 0;
  const fB = rhoB > 1e-20 ? Math.pow(rhoB, 1 / 3) : 0;
  // Note: negative rho values are clamped at the integration level (xcIntegration.ts)
  const twoThird = Math.pow(2, 1 / 3);
  const exc = -twoThird * CX * (rhoA * fA + rhoB * fB) / rho;
  const vrhoA = rhoA > 1e-20 ? -(4 / 3) * twoThird * CX * fA : 0;
  const vrhoB = rhoB > 1e-20 ? -(4 / 3) * twoThird * CX * fB : 0;
  return { exc, vrhoA, vrhoB };
}

// ── LDA: VWN5 correlation ──
// Vosko-Wilk-Nusair parameterization V (1980)
function vwn5Correlation(rhoA: number, rhoB: number): XCOutput {
  const rho = rhoA + rhoB;
  if (rho < 1e-20) return { exc: 0, vrhoA: 0, vrhoB: 0 };

  const rs = Math.pow(3 / (4 * Math.PI * rho), 1 / 3);
  const zeta = (rhoA - rhoB) / rho;

  // VWN5 parameters (Hartree): A values are Ry/2 from original paper
  const A0 = 0.0310907, b0 = 3.72744, c0 = 12.9352, x00 = -0.10498;
  const A1 = 0.01554535, b1 = 7.06042, c1 = 18.0578, x01 = -0.32500;
  // Alpha_c parameters
  const Aa = -1 / (6 * Math.PI * Math.PI), ba = 1.13107, ca = 13.0045, x0a = -0.0047584;

  function vwnEps(A: number, b: number, c: number, x0: number, sqrtrs: number): { eps: number; deps: number } {
    const x = sqrtrs;
    const X = x * x + b * x + c;
    const X0 = x0 * x0 + b * x0 + c;
    const Q = Math.sqrt(4 * c - b * b);

    const eps = A * (
      Math.log(x * x / X) +
      2 * b / Q * Math.atan(Q / (2 * x + b)) -
      b * x0 / X0 * (
        Math.log((x - x0) * (x - x0) / X) +
        2 * (b + 2 * x0) / Q * Math.atan(Q / (2 * x + b))
      )
    );

    // Derivative: dε/d(rs) = dε/dx * dx/drs = dε/dx * 1/(2*sqrt(rs))
    //
    // For atan(Q/(2x+b)) the chain rule gives:
    //   d/dx[atan(Q/(2x+b))] = -2Q / ((2x+b)² + Q²) = -2Q / denom1
    // So d/dx[(2b/Q)·atan(...)] = (2b/Q)·(-2Q/denom1) = -4b/denom1  (the Q's cancel)
    // (Earlier code mistakenly carried an extra 1/Q factor, giving V_c errors of ~1e-3.)
    const dXdx = 2 * x + b;
    const denom1 = (2 * x + b) * (2 * x + b) + Q * Q;
    const depsdx = A * (
      2 / x - dXdx / X -
      4 * b / denom1 -
      b * x0 / X0 * (
        2 / (x - x0) - dXdx / X -
        4 * (b + 2 * x0) / denom1
      )
    );

    const deps = depsdx / (2 * sqrtrs); // dε/drs
    return { eps, deps };
  }

  const sqrtrs = Math.sqrt(rs);
  const { eps: ec0, deps: dec0 } = vwnEps(A0, b0, c0, x00, sqrtrs);
  const { eps: ec1, deps: dec1 } = vwnEps(A1, b1, c1, x01, sqrtrs);
  const { eps: ac, deps: dac } = vwnEps(Aa, ba, ca, x0a, sqrtrs);

  // f(ζ) spin-interpolation
  const fzeta_denom = 2 * (Math.pow(2, 1 / 3) - 1);
  const f = ((1 + zeta) > 1e-20 ? Math.pow(1 + zeta, 4 / 3) : 0) +
            ((1 - zeta) > 1e-20 ? Math.pow(1 - zeta, 4 / 3) : 0) - 2;
  const fz = f / fzeta_denom;
  const fpp0 = 4 / (9 * fzeta_denom); // f''(0)

  const exc = ec0 + ac * fz / fpp0 * (1 - zeta * zeta * zeta * zeta) +
              (ec1 - ec0) * fz * zeta * zeta * zeta * zeta;

  // dε_c/drs
  const dec_drs = dec0 + dac * fz / fpp0 * (1 - zeta * zeta * zeta * zeta) +
                  (dec1 - dec0) * fz * zeta * zeta * zeta * zeta;

  // df/dζ
  const dfz_A = (1 + zeta) > 1e-20 ? (4 / 3) * Math.pow(1 + zeta, 1 / 3) : 0;
  const dfz_B = (1 - zeta) > 1e-20 ? -(4 / 3) * Math.pow(1 - zeta, 1 / 3) : 0;
  const dfzdz = (dfz_A + dfz_B) / fzeta_denom;

  const dexc_dzeta = dfzdz * (ac / fpp0 * (1 - zeta * zeta * zeta * zeta) +
                     (ec1 - ec0) * zeta * zeta * zeta * zeta) +
                     fz * (-4 * ac / fpp0 * zeta * zeta * zeta +
                     4 * (ec1 - ec0) * zeta * zeta * zeta);

  // V_c = ε_c - rs/3 * dε_c/drs ± (dε_c/dζ) correction
  // dζ/dρα = (1 - ζ)/ρ, dζ/dρβ = -(1 + ζ)/ρ
  const vc_common = exc - rs / 3 * dec_drs;
  const vrhoA = vc_common + dexc_dzeta * (1 - zeta) / (rho > 1e-20 ? 1 : 0);
  const vrhoB = vc_common - dexc_dzeta * (1 + zeta) / (rho > 1e-20 ? 1 : 0);

  return { exc, vrhoA: rho > 1e-20 ? vrhoA : 0, vrhoB: rho > 1e-20 ? vrhoB : 0 };
}

// ── GGA: B88 exchange ──
// ΔE_x^σ = -β ∫ ρ_σ^(4/3) x² / (1 + 6β x arcsinh(x)) dr
// where x = |∇ρ_σ| / ρ_σ^(4/3)
function b88Exchange(rhoA: number, rhoB: number,
                     gammaAA: number, gammaBB: number): XCOutput {
  const slater = slaterExchange(rhoA, rhoB);
  const beta = 0.0042;

  // B88 GGA correction energy density (per volume) for one spin
  const evalB88 = (rs: number, gs: number): number => {
    if (rs < 1e-20) return 0;
    const r43 = Math.pow(rs, 4 / 3);
    const x = Math.sqrt(Math.max(0, gs)) / r43;
    return -beta * r43 * x * x / (1 + 6 * beta * x * Math.asinh(x));
  };

  const excGGA = evalB88(rhoA, gammaAA) + evalB88(rhoB, gammaBB);

  // Finite-difference derivatives
  const h = 1e-6;
  const e0A = evalB88(rhoA, gammaAA);
  const e0B = evalB88(rhoB, gammaBB);
  const vrA = rhoA > 1e-20 ? (evalB88(rhoA + h, gammaAA) - evalB88(rhoA - h, gammaAA)) / (2 * h) : 0;
  const vrB = rhoB > 1e-20 ? (evalB88(rhoB + h, gammaBB) - evalB88(rhoB - h, gammaBB)) / (2 * h) : 0;
  const vgAA = rhoA > 1e-20 ? (evalB88(rhoA, gammaAA + h) - e0A) / h : 0;
  const vgBB = rhoB > 1e-20 ? (evalB88(rhoB, gammaBB + h) - e0B) / h : 0;

  const rho = rhoA + rhoB;
  return {
    exc: rho > 1e-20 ? (slater.exc * rho + excGGA) / rho : 0,
    vrhoA: slater.vrhoA + vrA,
    vrhoB: slater.vrhoB + vrB,
    vgammaAA: vgAA,
    vgammaBB: vgBB,
    vgammaAB: 0,
  };
}

// ── GGA: LYP correlation ──
function lypCorrelation(rhoA: number, rhoB: number,
                        gammaAA: number, gammaBB: number, gammaAB: number): XCOutput {
  const rho = rhoA + rhoB;
  if (rho < 1e-20) return { exc: 0, vrhoA: 0, vrhoB: 0, vgammaAA: 0, vgammaBB: 0, vgammaAB: 0 };

  const a = 0.04918, b = 0.132, c = 0.2533, d = 0.349;
  const CF = 0.3 * Math.pow(3 * Math.PI * Math.PI, 2 / 3);

  const rho13 = Math.pow(rho, -1 / 3);
  const omega = Math.exp(-c * rho13) / (1 + d * rho13) * Math.pow(rho, -11 / 3);
  const delta = c * rho13 + d * rho13 / (1 + d * rho13);
  const gamma = gammaAA + gammaBB + 2 * gammaAB;

  const term1 = -a * 4 * rhoA * rhoB / (rho * (1 + d * rho13));
  const term2 = -a * b * omega * rhoA * rhoB * (
    (47 / 18 - 7 * delta / 18) * gamma +
    (-5 / 2 + delta / 18) * (gammaAA + gammaBB) +
    (11 - delta) / 9 * (rhoA * gammaAA / rho + rhoB * gammaBB / rho)
  );
  const term3 = -a * b * omega * (
    -2 / 3 * rho * rho * gamma +
    (2 / 3 * rho * rho - rhoA * rhoA) * gammaBB +
    (2 / 3 * rho * rho - rhoB * rhoB) * gammaAA
  );
  const CF_scaled = Math.pow(2, 11 / 3) * CF; // spin-scaled Thomas-Fermi
  const term4 = -a * b * omega * rhoA * rhoB * (
    CF_scaled * (Math.pow(Math.max(0, rhoA), 8 / 3) + Math.pow(Math.max(0, rhoB), 8 / 3))
  );

  const exc_val = (term1 + term2 + term3 + term4) / rho;

  // For LYP derivatives, use finite differences for simplicity in Phase 1
  const h = 1e-6;
  const evalLYP = (rA: number, rB: number, gAA: number, gBB: number, gAB: number): number => {
    const r = rA + rB;
    if (r < 1e-20 || rA < 0 || rB < 0) return 0;
    const r13 = Math.pow(r, -1 / 3);
    const om = Math.exp(-c * r13) / (1 + d * r13) * Math.pow(r, -11 / 3);
    const dl = c * r13 + d * r13 / (1 + d * r13);
    const g = gAA + gBB + 2 * gAB;
    const t1 = -a * 4 * rA * rB / (r * (1 + d * r13));
    const t2 = -a * b * om * rA * rB * (
      (47 / 18 - 7 * dl / 18) * g +
      (-5 / 2 + dl / 18) * (gAA + gBB) +
      (11 - dl) / 9 * (rA * gAA / r + rB * gBB / r)
    );
    const t3 = -a * b * om * (
      -2 / 3 * r * r * g +
      (2 / 3 * r * r - rA * rA) * gBB +
      (2 / 3 * r * r - rB * rB) * gAA
    );
    const cfs = Math.pow(2, 11 / 3) * CF;
    const t4 = -a * b * om * rA * rB * cfs * (Math.pow(rA, 8 / 3) + Math.pow(rB, 8 / 3));
    return t1 + t2 + t3 + t4;
  };

  const e0 = evalLYP(rhoA, rhoB, gammaAA, gammaBB, gammaAB);
  // Use one-sided difference when spin density is too small for central difference
  const vrA = rhoA > h
    ? (evalLYP(rhoA + h, rhoB, gammaAA, gammaBB, gammaAB) -
       evalLYP(rhoA - h, rhoB, gammaAA, gammaBB, gammaAB)) / (2 * h)
    : (evalLYP(rhoA + h, rhoB, gammaAA, gammaBB, gammaAB) - e0) / h;
  const vrB = rhoB > h
    ? (evalLYP(rhoA, rhoB + h, gammaAA, gammaBB, gammaAB) -
       evalLYP(rhoA, rhoB - h, gammaAA, gammaBB, gammaAB)) / (2 * h)
    : (evalLYP(rhoA, rhoB + h, gammaAA, gammaBB, gammaAB) - e0) / h;
  const vgAA = (evalLYP(rhoA, rhoB, gammaAA + h, gammaBB, gammaAB) - e0) / h;
  const vgBB = (evalLYP(rhoA, rhoB, gammaAA, gammaBB + h, gammaAB) - e0) / h;
  const vgAB = (evalLYP(rhoA, rhoB, gammaAA, gammaBB, gammaAB + h) - e0) / h;

  return { exc: exc_val, vrhoA: vrA, vrhoB: vrB, vgammaAA: vgAA, vgammaBB: vgBB, vgammaAB: vgAB };
}

// ── PBE exchange ──
function pbeExchange(rhoA: number, rhoB: number,
                     gammaAA: number, gammaBB: number): XCOutput {
  const slater = slaterExchange(rhoA, rhoB);
  const kappa = 0.8040, mu = 0.2195149727645171;

  let excGGA = 0, vrA = 0, vrB = 0, vgAA = 0, vgBB = 0;

  for (const [rho_s, gamma_ss, isA] of [
    [rhoA, gammaAA, true], [rhoB, gammaBB, false],
  ] as [number, number, boolean][]) {
    if (rho_s < 1e-20) continue;
    const kf = Math.pow(6 * Math.PI * Math.PI * rho_s, 1 / 3);
    const s2 = gamma_ss / (4 * rho_s * rho_s * kf * kf);
    const Fx = 1 + kappa - kappa / (1 + mu * s2 / kappa);
    const ex_unif = -Math.pow(2, 1 / 3) * CX * Math.pow(rho_s, 1 / 3);
    const exc_s = rho_s * ex_unif * (Fx - 1); // GGA correction only
    excGGA += exc_s;

    // Simplified finite-diff derivative for Phase 1
    const h = 1e-6;
    const evalPBE_s = (r: number, g: number): number => {
      if (r < 1e-20) return 0;
      const k = Math.pow(6 * Math.PI * Math.PI * r, 1 / 3);
      const ss = g / (4 * r * r * k * k);
      const F = 1 + kappa - kappa / (1 + mu * ss / kappa);
      return r * (-Math.pow(2, 1 / 3) * CX * Math.pow(r, 1 / 3)) * (F - 1);
    };
    const vr_s = (evalPBE_s(rho_s + h, gamma_ss) - evalPBE_s(rho_s - h, gamma_ss)) / (2 * h);
    const vg_s = (evalPBE_s(rho_s, gamma_ss + h) - evalPBE_s(rho_s, gamma_ss - h)) / (2 * h);

    if (isA) { vrA = vr_s; vgAA = vg_s; }
    else { vrB = vr_s; vgBB = vg_s; }
  }

  const rho = rhoA + rhoB;
  return {
    exc: rho > 1e-20 ? (slater.exc * rho + excGGA) / rho : 0,
    vrhoA: slater.vrhoA + vrA,
    vrhoB: slater.vrhoB + vrB,
    vgammaAA: vgAA,
    vgammaBB: vgBB,
    vgammaAB: 0,
  };
}

// ── PBE correlation ──
function pbeCorrelation(rhoA: number, rhoB: number,
                        gammaAA: number, gammaBB: number, gammaAB: number): XCOutput {
  // PBE correlation on top of VWN5
  const vwn = vwn5Correlation(rhoA, rhoB);
  const rho = rhoA + rhoB;
  if (rho < 1e-20) return vwn;

  const beta_PBE = 0.06672455060314922;
  const gamma_PBE = 0.031090690869654895;

  const zeta = (rhoA - rhoB) / rho;
  const phi = 0.5 * ((1 + zeta) > 1e-20 ? Math.pow(1 + zeta, 2 / 3) : 0) +
              0.5 * ((1 - zeta) > 1e-20 ? Math.pow(1 - zeta, 2 / 3) : 0);
  const phi3 = phi * phi * phi;
  const kF = Math.pow(3 * Math.PI * Math.PI * rho, 1 / 3);
  const ks = Math.sqrt(4 * kF / Math.PI);
  const t2 = (gammaAA + gammaBB + 2 * gammaAB) / (4 * phi * phi * ks * ks * rho * rho);

  const A = beta_PBE / gamma_PBE / (Math.exp(-vwn.exc / gamma_PBE) - 1 || 1e-30);
  const At2 = A * t2;
  const H = gamma_PBE * phi3 * Math.log(1 + beta_PBE / gamma_PBE * t2 * (1 + At2) / (1 + At2 + At2 * At2));

  // Use finite differences for derivatives
  // evalH returns rho*H (energy per volume) so finite-diff gives correct KS potential
  const h = 1e-6;
  const evalH = (rA: number, rB: number, gAA: number, gBB: number, gAB: number): number => {
    const r = rA + rB;
    if (r < 1e-20 || rA < 0 || rB < 0) return 0;
    const vwnE = vwn5Correlation(rA, rB).exc;
    const z = (rA - rB) / r;
    const ph = 0.5 * ((1 + z) > 1e-20 ? Math.pow(1 + z, 2 / 3) : 0) +
               0.5 * ((1 - z) > 1e-20 ? Math.pow(1 - z, 2 / 3) : 0);
    const ph3 = ph * ph * ph;
    const kF_ = Math.pow(3 * Math.PI * Math.PI * r, 1 / 3);
    const ks_ = Math.sqrt(4 * kF_ / Math.PI);
    const t2_ = (gAA + gBB + 2 * gAB) / (4 * ph * ph * ks_ * ks_ * r * r);
    const A_ = beta_PBE / gamma_PBE / (Math.exp(-vwnE / gamma_PBE) - 1 || 1e-30);
    const At2_ = A_ * t2_;
    const H_ = gamma_PBE * ph3 * Math.log(1 + beta_PBE / gamma_PBE * t2_ * (1 + At2_) / (1 + At2_ + At2_ * At2_));
    return r * H_; // return rho*H (energy per volume)
  };

  const H0 = evalH(rhoA, rhoB, gammaAA, gammaBB, gammaAB);
  const vrA = vwn.vrhoA + (rhoA > h
    ? (evalH(rhoA + h, rhoB, gammaAA, gammaBB, gammaAB) -
       evalH(rhoA - h, rhoB, gammaAA, gammaBB, gammaAB)) / (2 * h)
    : (evalH(rhoA + h, rhoB, gammaAA, gammaBB, gammaAB) - H0) / h);
  const vrB = vwn.vrhoB + (rhoB > h
    ? (evalH(rhoA, rhoB + h, gammaAA, gammaBB, gammaAB) -
       evalH(rhoA, rhoB - h, gammaAA, gammaBB, gammaAB)) / (2 * h)
    : (evalH(rhoA, rhoB + h, gammaAA, gammaBB, gammaAB) - H0) / h);
  const vgAA = (evalH(rhoA, rhoB, gammaAA + h, gammaBB, gammaAB) - H0) / h;
  const vgBB = (evalH(rhoA, rhoB, gammaAA, gammaBB + h, gammaAB) - H0) / h;
  const vgAB = (evalH(rhoA, rhoB, gammaAA, gammaBB, gammaAB + h) - H0) / h;

  return {
    exc: vwn.exc + H,
    vrhoA: vrA, vrhoB: vrB,
    vgammaAA: vgAA, vgammaBB: vgBB, vgammaAB: vgAB,
  };
}

// ── Functional constructors ──

function makeSVWN(): XCFunctional {
  return {
    name: 'SVWN',
    type: 'LDA',
    exactExchangeFraction: 0,
    evaluate(inp: XCInput): XCOutput {
      const ex = slaterExchange(inp.rhoA, inp.rhoB);
      const ec = vwn5Correlation(inp.rhoA, inp.rhoB);
      return {
        exc: ex.exc + ec.exc,
        vrhoA: ex.vrhoA + ec.vrhoA,
        vrhoB: ex.vrhoB + ec.vrhoB,
      };
    },
  };
}

function makeBLYP(): XCFunctional {
  return {
    name: 'BLYP',
    type: 'GGA',
    exactExchangeFraction: 0,
    evaluate(inp: XCInput): XCOutput {
      const ex = b88Exchange(inp.rhoA, inp.rhoB, inp.gammaAA ?? 0, inp.gammaBB ?? 0);
      const ec = lypCorrelation(inp.rhoA, inp.rhoB, inp.gammaAA ?? 0, inp.gammaBB ?? 0, inp.gammaAB ?? 0);
      return {
        exc: ex.exc + ec.exc,
        vrhoA: ex.vrhoA + ec.vrhoA,
        vrhoB: ex.vrhoB + ec.vrhoB,
        vgammaAA: (ex.vgammaAA ?? 0) + (ec.vgammaAA ?? 0),
        vgammaBB: (ex.vgammaBB ?? 0) + (ec.vgammaBB ?? 0),
        vgammaAB: (ex.vgammaAB ?? 0) + (ec.vgammaAB ?? 0),
      };
    },
  };
}

function makePBE(): XCFunctional {
  return {
    name: 'PBE',
    type: 'GGA',
    exactExchangeFraction: 0,
    evaluate(inp: XCInput): XCOutput {
      const ex = pbeExchange(inp.rhoA, inp.rhoB, inp.gammaAA ?? 0, inp.gammaBB ?? 0);
      const ec = pbeCorrelation(inp.rhoA, inp.rhoB, inp.gammaAA ?? 0, inp.gammaBB ?? 0, inp.gammaAB ?? 0);
      return {
        exc: ex.exc + ec.exc,
        vrhoA: ex.vrhoA + ec.vrhoA,
        vrhoB: ex.vrhoB + ec.vrhoB,
        vgammaAA: (ex.vgammaAA ?? 0) + (ec.vgammaAA ?? 0),
        vgammaBB: (ex.vgammaBB ?? 0) + (ec.vgammaBB ?? 0),
        vgammaAB: (ex.vgammaAB ?? 0) + (ec.vgammaAB ?? 0),
      };
    },
  };
}

function makePBE0(): XCFunctional {
  // E_xc^PBE0 = 0.25 E_x^HF + 0.75 E_x^PBE + E_c^PBE
  return {
    name: 'PBE0',
    type: 'hybrid',
    exactExchangeFraction: 0.25,
    evaluate(inp: XCInput): XCOutput {
      const ax = 0.75; // PBE exchange fraction (0.25 HF replaces the rest)
      const gAA = inp.gammaAA ?? 0, gBB = inp.gammaBB ?? 0, gAB = inp.gammaAB ?? 0;
      const ex = pbeExchange(inp.rhoA, inp.rhoB, gAA, gBB);
      const ec = pbeCorrelation(inp.rhoA, inp.rhoB, gAA, gBB, gAB);
      return {
        exc: ax * ex.exc + ec.exc,
        vrhoA: ax * ex.vrhoA + ec.vrhoA,
        vrhoB: ax * ex.vrhoB + ec.vrhoB,
        vgammaAA: ax * (ex.vgammaAA ?? 0) + (ec.vgammaAA ?? 0),
        vgammaBB: ax * (ex.vgammaBB ?? 0) + (ec.vgammaBB ?? 0),
        vgammaAB: ax * (ex.vgammaAB ?? 0) + (ec.vgammaAB ?? 0),
      };
    },
  };
}

// ── TPSS meta-GGA (Tao, Perdew, Staroverov, Scuseria, PRL 91, 146401, 2003) ──

/** TPSS exchange contribution from one spin channel (energy density per volume).
 *  Direct spin convention: e_x^σ = -2^(1/3) C_x ρ_σ^(4/3) · F_x(p_pol, α_pol)
 *  where F_x arguments are computed from the polarised reference (2ρ_σ). */
function tpssExchangeEnergyDensity(rhoSig: number, gammaSig: number, tauSig: number): number {
  if (rhoSig < 1e-20) return 0;
  // Clamp physical inputs: |∇ρ|² and τ are non-negative by definition.
  // FD-perturbed γ in V_xc evaluation can otherwise go slightly negative and
  // propagate as NaN through tauW = sigma/(8ρ) → z = min(1, sigma/0) → -∞ in t3/t5.
  if (gammaSig < 0) gammaSig = 0;
  if (tauSig < 0) tauSig = 0;
  // Polarised reference for F_x evaluation
  const rhoPol = 2 * rhoSig;
  const sigmaPol = 4 * gammaSig;        // |∇(2ρ)|²
  const tauPol = 2 * tauSig;            // 2τ_σ for polarised reference

  // Reduced gradient p = s²
  const denom_p = 4 * Math.pow(3 * Math.PI * Math.PI, 2 / 3) * Math.pow(rhoPol, 8 / 3);
  let p = sigmaPol / denom_p; if (p < 0) p = 0;
  const tauUnif = 0.3 * Math.pow(3 * Math.PI * Math.PI, 2 / 3) * Math.pow(rhoPol, 5 / 3);
  const tauW = sigmaPol / (8 * rhoPol);
  const tauPos = Math.max(tauPol, tauW);
  // z = τ_W/τ_pos ∈ [0, 1]. When tauPos=0 (no kinetic density at all, e.g. when
  // FD-perturbed γ is clamped to 0 and τ is unset), 0/0 is NaN — convention: 0.
  const z = tauPos > 1e-30 ? Math.min(1.0, tauW / tauPos) : 0;
  const alpha = Math.max(0, (tauPos - tauW) / Math.max(tauUnif, 1e-30));

  // TPSS constants
  const KAPPA = 0.804, MU = 0.21951, b = 0.40, c = 1.59096, e = 1.537;

  const safe1 = 1 + b * alpha * (alpha - 1);
  const sqrt_safe = Math.sqrt(Math.max(safe1, 1e-30));
  const qb = (9 / 20) * (alpha - 1) / sqrt_safe + (2 / 3) * p;

  const z2 = z * z;
  const onePz2 = 1 + z2;
  const t1 = (10 / 81 + c * z2 / (onePz2 * onePz2)) * p;
  const t2 = (146 / 2025) * qb * qb;
  // t3: -(73/405) q̃_b · sqrt[(3/5)² z² + (1/2) p²]
  const t3 = -(73 / 405) * qb * Math.sqrt(0.36 * z2 + 0.5 * p * p);
  const t4 = (1 / KAPPA) * Math.pow(10 / 81, 2) * p * p;
  // t5: 2√e · (10/81) · (3/5)² · z²   (note: (3/5)² is a plain number, not (3π²)^(2/3))
  const t5 = 2 * Math.sqrt(e) * (10 / 81) * 0.36 * z2;
  const t6 = e * MU * p * p * p;
  const denomX = (1 + Math.sqrt(e) * p);
  const xx = (t1 + t2 + t3 + t4 + t5 + t6) / (denomX * denomX);

  const Fx = 1 + KAPPA - KAPPA / (1 + xx / KAPPA);

  // Direct spin formula matching existing slaterExchange convention:
  //   e_x[ρ_α, ρ_β] = -2^(1/3) C_x (ρ_α^(4/3) + ρ_β^(4/3))   (LDA, F_x=1)
  // → per-spin contribution = -2^(1/3) C_x ρ_σ^(4/3) · F_x_σ
  return -Math.pow(2, 1 / 3) * CX * Math.pow(rhoSig, 4 / 3) * Fx;
}

/** TPSS exchange (sums spin contributions, returns total ε_x at point and density-weighted form). */
function tpssExchange(rhoA: number, rhoB: number,
  gAA: number, gBB: number,
  tauA: number, tauB: number): { excTimesRho: number } {
  const eA = rhoA > 1e-20 ? tpssExchangeEnergyDensity(rhoA, gAA, tauA) : 0;
  const eB = rhoB > 1e-20 ? tpssExchangeEnergyDensity(rhoB, gBB, tauB) : 0;
  return { excTimesRho: eA + eB };
}

/** TPSS correlation (full PKZB form, Tao-Perdew-Staroverov-Scuseria 2003, Eq. 12-14).
 *
 *   ε_c^PKZB = ε_c^PBE · (1 + C(ζ, ξ)·R²)
 *            - (1 + C(ζ, ξ))·R²·Σ_σ (ρ_σ/ρ) · ε_c^PBE(ρ_σ, 0, ∇ρ_σ)
 *   ε_c^TPSS = ε_c^PKZB · (1 + d · ε_c^PKZB · R³)
 *
 *   R = τ_W/τ,  τ_W = |∇ρ|²/(8ρ),  τ = τ_α + τ_β
 *   C(ζ, ξ) = C(ζ, 0) / [1 + ξ²·((1+ζ)^{-4/3} + (1-ζ)^{-4/3})/2]^4
 *   C(ζ, 0) = 0.53 + 0.87 ζ² + 0.50 ζ⁴ + 2.26 ζ⁶
 *   ξ = |∇ζ|/(2·(3π²ρ)^{1/3})
 *   d = 2.8
 *
 *   The "self-interaction correction" piece (subtracts spin-isolated PBE) suppresses
 *   correlation in single-orbital regions where R → 1. Returns the per-volume
 *   energy density ε_c·ρ. */
function tpssCorrelationEnergyDensity(rhoA: number, rhoB: number,
  gAA: number, gBB: number, gAB: number,
  tauA: number, tauB: number): number {
  const rho = rhoA + rhoB;
  if (rho < 1e-20) return 0;
  // Clamp non-negative inputs (see tpssExchangeEnergyDensity comment)
  if (gAA < 0) gAA = 0;
  if (gBB < 0) gBB = 0;
  if (tauA < 0) tauA = 0;
  if (tauB < 0) tauB = 0;
  // gAB can legitimately be negative (∇ρ_α · ∇ρ_β); leave it alone.
  const tau = tauA + tauB;
  const gammaTotal = gAA + 2 * gAB + gBB;  // |∇ρ|²
  const tauW = gammaTotal / (8 * rho);
  // R = τ_W / τ; clip to [0, 1] (physical bound; numerical noise may push slightly above)
  let R = tau > 1e-20 ? tauW / tau : 0;
  if (R > 1) R = 1;
  if (R < 0) R = 0;

  // Full system PBE correlation per particle
  const pbeFull = pbeCorrelation(rhoA, rhoB, gAA, gBB, gAB).exc;

  // Spin-isolated PBE: each spin treated as fully spin-polarized (ρ_σ, 0)
  // Note: PBE c gives correlation energy per particle of the spin-isolated system.
  // The "max" in the original TPSS paper is for a different formulation;
  // here we use the standard spin-isolated form per-spin.
  const pbeA = rhoA > 1e-20 ? pbeCorrelation(rhoA, 0, gAA, 0, 0).exc : 0;
  const pbeB = rhoB > 1e-20 ? pbeCorrelation(0, rhoB, 0, gBB, 0).exc : 0;
  const sumSig = (rhoA / rho) * pbeA + (rhoB / rho) * pbeB;

  // Spin polarization ζ and its reduced gradient ξ
  const zeta = (rhoA - rhoB) / rho;
  // |∇ζ|² = |∇(ρ_α - ρ_β)/ρ - (ρ_α - ρ_β)·∇ρ/ρ²|²
  // For the closed-shell case (ζ=0) this is well-defined; for the polarised case
  // we approximate using the spin-channel gradients:
  //   ∇ρ_σ = (1±ζ)/2 · ∇ρ + ρ/2 · (±)∇ζ  → ∇ζ = 2·(∇ρ_α - ρ_α/ρ·∇ρ)/ρ_total
  // Simpler: use |∇ρ_α - ρ_α/ρ·∇ρ_total|² · (2/ρ)² for |∇ζ|²
  // |∇ρ_α|² = gAA, |∇ρ_β|² = gBB, ∇ρ_α·∇ρ_β = gAB
  // |∇ζ|² = |2/ρ · (∇ρ_α - (ρ_α/ρ)·∇ρ_total)|²
  //       = (4/ρ²) · [|∇ρ_α|² - 2(ρ_α/ρ)·∇ρ_α·∇ρ_total + (ρ_α/ρ)²·|∇ρ_total|²]
  const xa = rhoA / rho;
  const dotAtot = gAA + gAB;     // ∇ρ_α · ∇ρ_total
  const gradZeta2 = (4 / (rho * rho)) * (gAA - 2 * xa * dotAtot + xa * xa * gammaTotal);
  const kF = Math.pow(3 * Math.PI * Math.PI * rho, 1 / 3);
  const xi2 = gradZeta2 / Math.max(4 * kF * kF, 1e-30);

  // C(ζ, 0) — TPSS Eq. 13
  const z2 = zeta * zeta;
  const C0 = 0.53 + 0.87 * z2 + 0.50 * z2 * z2 + 2.26 * z2 * z2 * z2;
  // Spin-symmetry suppression factor for |ζ| ≈ 1 (singular; clip)
  const oneP = 1 + zeta, oneM = 1 - zeta;
  const oneP43 = oneP > 1e-10 ? Math.pow(oneP, -4 / 3) : 1e10;
  const oneM43 = oneM > 1e-10 ? Math.pow(oneM, -4 / 3) : 1e10;
  const denomC4 = Math.pow(1 + xi2 * (oneP43 + oneM43) / 2, 4);
  const C = C0 / Math.max(denomC4, 1e-30);

  const R2 = R * R;
  const R3 = R2 * R;

  // PKZB part
  const pkzb = pbeFull * (1 + C * R2) - (1 + C) * R2 * sumSig;

  // TPSS revision: multiply by (1 + d·PKZB·R³)
  const D_TPSS = 2.8;
  const tpss = pkzb * (1 + D_TPSS * pkzb * R3);

  return tpss * rho;  // per-volume energy density
}

/** Numerical-derivative helper: returns dF/dx via central FD. */
function fdDeriv(f: () => number, x: number, setX: (v: number) => void, eps: number): number {
  const h = Math.max(eps, Math.abs(x) * 1e-6);
  setX(x + h); const fp = f();
  setX(x - h); const fm = f();
  setX(x);
  return (fp - fm) / (2 * h);
}

function makeTPSS(): XCFunctional {
  return {
    name: 'TPSS',
    type: 'meta-GGA',
    exactExchangeFraction: 0,
    evaluate(inp: XCInput): XCOutput {
      const rho = inp.rhoA + inp.rhoB;
      if (rho < 1e-20) {
        return { exc: 0, vrhoA: 0, vrhoB: 0, vgammaAA: 0, vgammaBB: 0, vgammaAB: 0, vtauA: 0, vtauB: 0 };
      }
      const gAA = inp.gammaAA ?? 0, gBB = inp.gammaBB ?? 0, gAB = inp.gammaAB ?? 0;
      const tA = inp.tauA ?? 0, tB = inp.tauB ?? 0;

      // Local mutable copies for FD
      let rhoA = inp.rhoA, rhoB = inp.rhoB;
      let gammaAA = gAA, gammaBB = gBB, gammaAB = gAB;
      let tauA = tA, tauB = tB;

      // Energy-times-rho (so that ε_xc returned per-electron)
      const ETR = (): number => {
        const ex = tpssExchange(rhoA, rhoB, gammaAA, gammaBB, tauA, tauB).excTimesRho;
        const ec = tpssCorrelationEnergyDensity(rhoA, rhoB, gammaAA, gammaBB, gammaAB, tauA, tauB);
        return ex + ec;
      };

      const etr0 = ETR();
      const exc = etr0 / rho; // per-electron

      // Numerical derivatives w.r.t. each variable.
      // V = ∂(ε·ρ)/∂var = derivative of total energy density.
      //
      // Use a density-relative FD step. fdDeriv internally takes max(eps, |x|·1e-6),
      // so a smaller "eps" floor lets it scale better at low densities. Roundoff
      // floor at 1e-12 (~ sqrt(double precision)).
      const epsRho   = Math.max(1e-12, rho * 1e-6);
      const epsGamma = Math.max(1e-12, Math.max(gammaAA, gammaBB, Math.abs(gammaAB)) * 1e-6);
      const epsTau   = Math.max(1e-12, Math.max(tauA, tauB) * 1e-6);
      const vrhoA    = fdDeriv(ETR, rhoA,    v => { rhoA    = v; }, epsRho);
      const vrhoB    = fdDeriv(ETR, rhoB,    v => { rhoB    = v; }, epsRho);
      const vgammaAA = fdDeriv(ETR, gammaAA, v => { gammaAA = v; }, epsGamma);
      const vgammaBB = fdDeriv(ETR, gammaBB, v => { gammaBB = v; }, epsGamma);
      const vgammaAB = fdDeriv(ETR, gammaAB, v => { gammaAB = v; }, epsGamma);
      const vtauA    = fdDeriv(ETR, tauA,    v => { tauA    = v; }, epsTau);
      const vtauB    = fdDeriv(ETR, tauB,    v => { tauB    = v; }, epsTau);

      // Clamp: at low-density basis-set artifact points the reduced gradient
      // s = √γ/(2k_F ρ) can become unphysically large (~150), making the
      // TPSS F_x formula nearly singular and FD-derivative magnitudes blow up
      // (10^10+ for vγ). Clamp to a generous physical bound to avoid SCF
      // divergence — these points contribute weight·V·D ≈ tiny anyway.
      const VMAX = 1e3;
      const safe = (v: number) => Math.abs(v) > VMAX ? Math.sign(v) * VMAX : v;
      return {
        exc,
        vrhoA: safe(vrhoA), vrhoB: safe(vrhoB),
        vgammaAA: safe(vgammaAA), vgammaBB: safe(vgammaBB), vgammaAB: safe(vgammaAB),
        vtauA: safe(vtauA), vtauB: safe(vtauB),
      };
    },
  };
}

function makeB3LYP(): XCFunctional {
  // E_xc = (1-a0)*E_x^Slater + a0*E_x^HF + ax*ΔE_x^B88 + (1-ac)*E_c^VWN + ac*E_c^LYP
  // a0 = 0.20, ax = 0.72, ac = 0.81
  return {
    name: 'B3LYP',
    type: 'hybrid',
    exactExchangeFraction: 0.20,
    evaluate(inp: XCInput): XCOutput {
      const a0 = 0.20, ax = 0.72, ac = 0.81;
      const gAA = inp.gammaAA ?? 0, gBB = inp.gammaBB ?? 0, gAB = inp.gammaAB ?? 0;

      const slater = slaterExchange(inp.rhoA, inp.rhoB);
      const b88 = b88Exchange(inp.rhoA, inp.rhoB, gAA, gBB);
      // B88 GGA correction = b88 total - slater
      const b88corr_exc = b88.exc - slater.exc;
      const b88corr_vrhoA = b88.vrhoA - slater.vrhoA;
      const b88corr_vrhoB = b88.vrhoB - slater.vrhoB;

      const vwn = vwn5Correlation(inp.rhoA, inp.rhoB);
      const lyp = lypCorrelation(inp.rhoA, inp.rhoB, gAA, gBB, gAB);

      return {
        exc: (1 - a0) * slater.exc + ax * b88corr_exc + (1 - ac) * vwn.exc + ac * lyp.exc,
        vrhoA: (1 - a0) * slater.vrhoA + ax * b88corr_vrhoA + (1 - ac) * vwn.vrhoA + ac * lyp.vrhoA,
        vrhoB: (1 - a0) * slater.vrhoB + ax * b88corr_vrhoB + (1 - ac) * vwn.vrhoB + ac * lyp.vrhoB,
        vgammaAA: ax * (b88.vgammaAA ?? 0) + ac * (lyp.vgammaAA ?? 0),
        vgammaBB: ax * (b88.vgammaBB ?? 0) + ac * (lyp.vgammaBB ?? 0),
        vgammaAB: ac * (lyp.vgammaAB ?? 0),
      };
    },
  };
}

// ── Range-Separated Hybrids ──

/** Iikura-Tsuneda-Yanai-Hirao (2001) short-range LDA exchange attenuation factor.
 *  F^SR(a) gives the fraction of LDA exchange in the SHORT-range region only.
 *  a = ω / (2 k_F),   k_F = (3π²ρ)^(1/3) */
function lsdaSRfactor(rho: number, omega: number): number {
  if (rho < 1e-20 || omega <= 0) return 1.0; // degenerate to full
  const kF = Math.pow(3 * Math.PI * Math.PI * rho, 1 / 3);
  const a = omega / (2 * kF);
  if (a > 50) return 0; // a→∞ limit: SR part vanishes
  // Iikura formula: F^SR(a) = 1 - (8/3) a [√π·erf(1/(2a)) + (2a - 4a³) exp(-1/(4a²)) - 3a + 4a³]
  const erf_inv2a = lsdaErf(1 / (2 * a));
  const exp_inv4a2 = Math.exp(-1 / (4 * a * a));
  const inner = Math.sqrt(Math.PI) * erf_inv2a + (2 * a - 4 * a * a * a) * exp_inv4a2 - 3 * a + 4 * a * a * a;
  let F = 1 - (8 / 3) * a * inner;
  if (F < 0) F = 0;
  if (F > 1) F = 1;
  return F;
}

/** erf via Abramowitz-Stegun 7.1.26 (max error 1.5e-7), local copy to avoid import */
function lsdaErf(x: number): number {
  const sign = x < 0 ? -1 : 1;
  const t = 1 / (1 + 0.3275911 * Math.abs(x));
  const y = 1 - ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
  return sign * y;
}

/** Short-range Slater exchange (per-electron). */
function slaterExchangeSR(rhoA: number, rhoB: number, omega: number): XCOutput {
  const slater = slaterExchange(rhoA, rhoB);
  const rho = rhoA + rhoB;
  if (rho < 1e-20 || omega <= 0) return slater;
  // Use total density for the attenuation factor (simplification).
  const F = lsdaSRfactor(rho, omega);
  return {
    exc: slater.exc * F,
    vrhoA: slater.vrhoA * F,
    vrhoB: slater.vrhoB * F,
  };
}

/** Short-range B88 GGA exchange (lite implementation).
 *  Applies the LDA SR attenuation factor F^SR_LDA(ρ, ω) uniformly to the full B88
 *  exchange (energy density and all derivatives). This is the "lite" approach
 *  used in several simplified RSH implementations — it captures the main range
 *  separation effect on the GGA correction without requiring the full
 *  Henderson-Janesko-Scuseria / Akinaga-Ten-no s-dependent attenuation.
 *
 *  Strictly speaking the SR factor depends on both ρ and the reduced gradient
 *  s = |∇ρ|/(2 k_F ρ); using only F^SR_LDA underestimates SR attenuation in
 *  high-gradient regions (atomic cores, density tails). Adequate for pedagogical
 *  use in the textbook RSH chapter.
 */
function b88ExchangeSR(
  rhoA: number, rhoB: number, gammaAA: number, gammaBB: number, omega: number,
): XCOutput {
  const b88 = b88Exchange(rhoA, rhoB, gammaAA, gammaBB);
  const rho = rhoA + rhoB;
  if (rho < 1e-20 || omega <= 0) return b88;
  const F = lsdaSRfactor(rho, omega);
  return {
    exc: b88.exc * F,
    vrhoA: b88.vrhoA * F,
    vrhoB: b88.vrhoB * F,
    vgammaAA: (b88.vgammaAA ?? 0) * F,
    vgammaBB: (b88.vgammaBB ?? 0) * F,
  };
}



/** CAM-B3LYP (Yanai-Tew-Handy 2004): α=0.19, β=0.46, ω=0.33.
 *  Implementation: HF range-split correctly + SR-LDA attenuation of LDA exchange
 *  + SR-B88-lite (LDA SR factor applied uniformly to B88 GGA correction).
 *  DFT-exchange weights chosen so that HF + DFT operator weights sum to 1 across all r. */
function makeCAMB3LYP(): XCFunctional {
  const ALPHA = 0.19, BETA = 0.46, OMEGA = 0.33;
  return {
    name: 'CAM-B3LYP',
    type: 'hybrid',
    exactExchangeFraction: ALPHA,
    rangeSeparation: { alpha: ALPHA, beta: BETA, omega: OMEGA },
    evaluate(inp: XCInput): XCOutput {
      const ax = 0.72, ac = 0.81;
      const gAA = inp.gammaAA ?? 0, gBB = inp.gammaBB ?? 0, gAB = inp.gammaAB ?? 0;
      // Apply the same range separation to BOTH Slater (LDA exchange) and B88 GGA.
      // Total DFT exchange = (1-α) · X^SR + (1-α-β) · X^LR
      //                    = X^full · (1-α-β) + X^SR · β
      const slaterFull = slaterExchange(inp.rhoA, inp.rhoB);
      const slaterSR = slaterExchangeSR(inp.rhoA, inp.rhoB, OMEGA);
      const wSR = (1 - ALPHA);
      const wFullMinusSR = (1 - ALPHA - BETA);
      const slaterTotal_exc = wSR * slaterSR.exc + wFullMinusSR * (slaterFull.exc - slaterSR.exc);
      const slaterTotal_vrhoA = wSR * slaterSR.vrhoA + wFullMinusSR * (slaterFull.vrhoA - slaterSR.vrhoA);
      const slaterTotal_vrhoB = wSR * slaterSR.vrhoB + wFullMinusSR * (slaterFull.vrhoB - slaterSR.vrhoB);

      // SR-B88: same range-separation treatment as Slater.
      const b88Full = b88Exchange(inp.rhoA, inp.rhoB, gAA, gBB);
      const b88SR = b88ExchangeSR(inp.rhoA, inp.rhoB, gAA, gBB, OMEGA);
      const b88corrFull_exc = b88Full.exc - slaterFull.exc;
      const b88corrSR_exc = b88SR.exc - slaterSR.exc;
      const b88corrFull_vrhoA = b88Full.vrhoA - slaterFull.vrhoA;
      const b88corrSR_vrhoA = b88SR.vrhoA - slaterSR.vrhoA;
      const b88corrFull_vrhoB = b88Full.vrhoB - slaterFull.vrhoB;
      const b88corrSR_vrhoB = b88SR.vrhoB - slaterSR.vrhoB;
      const b88corrFull_vgAA = (b88Full.vgammaAA ?? 0);
      const b88corrSR_vgAA = (b88SR.vgammaAA ?? 0);  // slater has no γ derivative
      const b88corrFull_vgBB = (b88Full.vgammaBB ?? 0);
      const b88corrSR_vgBB = (b88SR.vgammaBB ?? 0);
      const b88corrTotal_exc = wSR * b88corrSR_exc + wFullMinusSR * (b88corrFull_exc - b88corrSR_exc);
      const b88corrTotal_vrhoA = wSR * b88corrSR_vrhoA + wFullMinusSR * (b88corrFull_vrhoA - b88corrSR_vrhoA);
      const b88corrTotal_vrhoB = wSR * b88corrSR_vrhoB + wFullMinusSR * (b88corrFull_vrhoB - b88corrSR_vrhoB);
      const b88corrTotal_vgAA = wSR * b88corrSR_vgAA + wFullMinusSR * (b88corrFull_vgAA - b88corrSR_vgAA);
      const b88corrTotal_vgBB = wSR * b88corrSR_vgBB + wFullMinusSR * (b88corrFull_vgBB - b88corrSR_vgBB);

      const vwn = vwn5Correlation(inp.rhoA, inp.rhoB);
      const lyp = lypCorrelation(inp.rhoA, inp.rhoB, gAA, gBB, gAB);

      return {
        exc: slaterTotal_exc + ax * b88corrTotal_exc + (1 - ac) * vwn.exc + ac * lyp.exc,
        vrhoA: slaterTotal_vrhoA + ax * b88corrTotal_vrhoA + (1 - ac) * vwn.vrhoA + ac * lyp.vrhoA,
        vrhoB: slaterTotal_vrhoB + ax * b88corrTotal_vrhoB + (1 - ac) * vwn.vrhoB + ac * lyp.vrhoB,
        vgammaAA: ax * b88corrTotal_vgAA + ac * (lyp.vgammaAA ?? 0),
        vgammaBB: ax * b88corrTotal_vgBB + ac * (lyp.vgammaBB ?? 0),
        vgammaAB: ac * (lyp.vgammaAB ?? 0),
      };
    },
  };
}

/** ωB97X-D (Chai-Head-Gordon 2008): α=0.157706, β=0.842294, ω=0.30.
 *  Long-range corrected (α+β=1, full HF at infinity).
 *  Implementation: SR-LDA exchange + SR-B88-lite GGA correction + LYP correlation.
 *  True ωB97X uses re-optimised ωB97 short-range exchange (different functional form). */
function makeWB97XD(): XCFunctional {
  const ALPHA = 0.157706, BETA = 0.842294, OMEGA = 0.30;
  return {
    name: 'ωB97X-D',
    type: 'hybrid',
    exactExchangeFraction: ALPHA,
    rangeSeparation: { alpha: ALPHA, beta: BETA, omega: OMEGA },
    evaluate(inp: XCInput): XCOutput {
      const gAA = inp.gammaAA ?? 0, gBB = inp.gammaBB ?? 0, gAB = inp.gammaAB ?? 0;
      const slaterFull = slaterExchange(inp.rhoA, inp.rhoB);
      const slaterSR = slaterExchangeSR(inp.rhoA, inp.rhoB, OMEGA);
      const wSR = (1 - ALPHA);
      const wFullMinusSR = (1 - ALPHA - BETA);
      const slaterTotal_exc = wSR * slaterSR.exc + wFullMinusSR * (slaterFull.exc - slaterSR.exc);
      const slaterTotal_vrhoA = wSR * slaterSR.vrhoA + wFullMinusSR * (slaterFull.vrhoA - slaterSR.vrhoA);
      const slaterTotal_vrhoB = wSR * slaterSR.vrhoB + wFullMinusSR * (slaterFull.vrhoB - slaterSR.vrhoB);

      // SR-B88: same range-separation treatment as Slater (with ax = 1.0 for ωB97X-D's pure GGA flavour).
      const b88Full = b88Exchange(inp.rhoA, inp.rhoB, gAA, gBB);
      const b88SR = b88ExchangeSR(inp.rhoA, inp.rhoB, gAA, gBB, OMEGA);
      const b88corrFull_exc = b88Full.exc - slaterFull.exc;
      const b88corrSR_exc = b88SR.exc - slaterSR.exc;
      const b88corrFull_vrhoA = b88Full.vrhoA - slaterFull.vrhoA;
      const b88corrSR_vrhoA = b88SR.vrhoA - slaterSR.vrhoA;
      const b88corrFull_vrhoB = b88Full.vrhoB - slaterFull.vrhoB;
      const b88corrSR_vrhoB = b88SR.vrhoB - slaterSR.vrhoB;
      const b88corrFull_vgAA = (b88Full.vgammaAA ?? 0);
      const b88corrSR_vgAA = (b88SR.vgammaAA ?? 0);
      const b88corrFull_vgBB = (b88Full.vgammaBB ?? 0);
      const b88corrSR_vgBB = (b88SR.vgammaBB ?? 0);
      const b88corr_exc = wSR * b88corrSR_exc + wFullMinusSR * (b88corrFull_exc - b88corrSR_exc);
      const b88corr_vrhoA = wSR * b88corrSR_vrhoA + wFullMinusSR * (b88corrFull_vrhoA - b88corrSR_vrhoA);
      const b88corr_vrhoB = wSR * b88corrSR_vrhoB + wFullMinusSR * (b88corrFull_vrhoB - b88corrSR_vrhoB);
      const b88corr_vgAA = wSR * b88corrSR_vgAA + wFullMinusSR * (b88corrFull_vgAA - b88corrSR_vgAA);
      const b88corr_vgBB = wSR * b88corrSR_vgBB + wFullMinusSR * (b88corrFull_vgBB - b88corrSR_vgBB);

      const lyp = lypCorrelation(inp.rhoA, inp.rhoB, gAA, gBB, gAB);
      return {
        exc: slaterTotal_exc + b88corr_exc + lyp.exc,
        vrhoA: slaterTotal_vrhoA + b88corr_vrhoA + lyp.vrhoA,
        vrhoB: slaterTotal_vrhoB + b88corr_vrhoB + lyp.vrhoB,
        vgammaAA: b88corr_vgAA + (lyp.vgammaAA ?? 0),
        vgammaBB: b88corr_vgBB + (lyp.vgammaBB ?? 0),
        vgammaAB: (lyp.vgammaAB ?? 0),
      };
    },
  };
}

export type FunctionalName = 'SVWN' | 'BLYP' | 'PBE' | 'PBE0' | 'TPSS' | 'B3LYP'
  | 'CAM-B3LYP' | 'ωB97X-D';

const FUNCTIONAL_REGISTRY: Record<FunctionalName, () => XCFunctional> = {
  SVWN: makeSVWN,
  BLYP: makeBLYP,
  PBE: makePBE,
  PBE0: makePBE0,
  TPSS: makeTPSS,
  B3LYP: makeB3LYP,
  'CAM-B3LYP': makeCAMB3LYP,
  'ωB97X-D': makeWB97XD,
};

export function createFunctional(name: FunctionalName): XCFunctional {
  return FUNCTIONAL_REGISTRY[name]();
}

export const FUNCTIONAL_NAMES: FunctionalName[] = ['SVWN', 'BLYP', 'PBE', 'PBE0', 'TPSS', 'B3LYP', 'CAM-B3LYP', 'ωB97X-D'];
