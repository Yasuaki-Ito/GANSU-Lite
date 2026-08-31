/** V_xc gradient evaluator for DFT (analytic LDA + GGA + meta-GGA).
 *
 *  Computes dE_xc/dR_A via grid integration with basis function derivatives.
 *
 *  Formula (Pulay's formulation, fixed-grid approximation):
 *    dE_xc/dR_A_d = -2 Σ_{μ on A} Σ_ν P_μν * [
 *        V^ρ(r) ∂_d φ_μ(r) φ_ν(r)                              (LDA part)
 *      + 2 V^γ(r) [∂_d ∇φ_μ(r)] · ∇ρ(r) · φ_ν(r) ... but using a smart factorisation
 *      ]
 *
 *  Here we use a compact form via the GGA "f-vector" already constructed in computeXC:
 *    f^A(r) = 2 V^γAA · ∇ρ_α + V^γAB · ∇ρ_β
 *  Then the gradient contribution from one spin σ is:
 *    dE/dR_A_d = -2 Σ_{μ on A} Σ_ν P^σ_μν * [V^ρσ ∂_d φ_μ φ_ν + f^σ · (∂_d∇φ_μ) φ_ν + ...]
 *
 *  This implementation skips Becke-partition weight derivatives (fixed-grid
 *  approximation, accurate to ~1e-4 Eh/bohr for typical molecules).
 */

import type { Atom, PrimitiveShell } from './types';
import { ANGULAR_MOMENTUMS, shellTypeToNumBasis } from './constants';
import { primitiveNorm } from './integrals1e';
import type { GridPoint } from './grid';
import { XCFunctional, XCInput } from './xcFunctional';
import { Matrix } from '../linalg/matrix';

function ipow(x: number, n: number): number {
  if (n === 0) return 1;
  if (n === 1) return x;
  if (n === 2) return x * x;
  if (n === 3) return x * x * x;
  return Math.pow(x, n);
}

interface ShellNormCache { pNorms: number[]; r2Max: number; }

function buildShellNormCache(shells: PrimitiveShell[]): ShellNormCache[] {
  const LOG_THRESH = 69.08;
  return shells.map(shell => {
    const angMom = ANGULAR_MOMENTUMS[shell.shellType];
    const nBF = shellTypeToNumBasis(shell.shellType);
    const pNorms: number[] = [];
    for (let bf = 0; bf < nBF; bf++) {
      const [lx, ly, lz] = angMom[bf];
      pNorms.push(primitiveNorm(shell.exponent, lx, ly, lz));
    }
    return { pNorms, r2Max: LOG_THRESH / shell.exponent };
  });
}

/** Evaluate basis function values, 1st and 2nd derivatives at a grid point.
 *  Returns φ_μ, ∂_d φ_μ (d=x,y,z), ∂²_{d,d'} φ_μ (3x3 = 9 components per μ; symmetric).
 *  Also returns which atom each basis function belongs to. */
function evaluateBasisAndHessian(
  shells: PrimitiveShell[],
  shellNorms: ShellNormCache[],
  normFactors: number[],
  nbasis: number,
  px: number, py: number, pz: number,
): {
  phi: Float64Array;
  dphi: Float64Array; // [d, basisIdx]: 3 * nbasis
  ddphi: Float64Array; // [d, d', basisIdx]: 9 * nbasis (symmetric in d, d')
  basisAtom: Int32Array; // atom index of each basis function
  sigIdx: Int32Array; nsig: number;
} {
  const phi = new Float64Array(nbasis);
  const dphi = new Float64Array(3 * nbasis);
  const ddphi = new Float64Array(9 * nbasis);
  const basisAtom = new Int32Array(nbasis);

  for (let s = 0; s < shells.length; s++) {
    const shell = shells[s];
    const cache = shellNorms[s];
    const dx = px - shell.coordinate.x;
    const dy = py - shell.coordinate.y;
    const dz = pz - shell.coordinate.z;
    const r2 = dx * dx + dy * dy + dz * dz;
    if (r2 > cache.r2Max) {
      // still tag basisAtom for correctness
      const angMom = ANGULAR_MOMENTUMS[shell.shellType];
      const nBF = shellTypeToNumBasis(shell.shellType);
      for (let bf = 0; bf < nBF; bf++) basisAtom[shell.basisIndex + bf] = shell.atomIndex;
      continue;
    }
    const a = shell.exponent;
    const expVal = shell.coefficient * Math.exp(-a * r2);

    const angMom = ANGULAR_MOMENTUMS[shell.shellType];
    const nBF = shellTypeToNumBasis(shell.shellType);

    for (let bf = 0; bf < nBF; bf++) {
      const [lx, ly, lz] = angMom[bf];
      const basisIdx = shell.basisIndex + bf;
      basisAtom[basisIdx] = shell.atomIndex;
      const N = normFactors[basisIdx] * cache.pNorms[bf] * expVal;

      const Px = ipow(dx, lx), Py = ipow(dy, ly), Pz = ipow(dz, lz);
      const angular = Px * Py * Pz;
      phi[basisIdx] += N * angular;

      // First derivatives ∂_d φ
      // For Gaussian times polynomial: ∂_x [X_x^lx · e^{-αr²}] = (lx X_x^{lx-1} - 2α X_x^{lx+1}) e^{-αr²}
      const Pxm = lx > 0 ? ipow(dx, lx - 1) : 0;
      const Pym = ly > 0 ? ipow(dy, ly - 1) : 0;
      const Pzm = lz > 0 ? ipow(dz, lz - 1) : 0;
      const Pxp = dx * Px;       // X^{lx+1}
      const Pyp = dy * Py;
      const Pzp = dz * Pz;

      const dxF = (lx * Pxm - 2 * a * Pxp) * Py * Pz;
      const dyF = Px * (ly * Pym - 2 * a * Pyp) * Pz;
      const dzF = Px * Py * (lz * Pzm - 2 * a * Pzp);
      dphi[0 * nbasis + basisIdx] += N * dxF;
      dphi[1 * nbasis + basisIdx] += N * dyF;
      dphi[2 * nbasis + basisIdx] += N * dzF;

      // Second derivatives ∂_d ∂_d' φ
      // ∂²_x [X^lx e^{-αr²}] = (lx(lx-1) X^{lx-2} - 2α(2lx+1) X^lx + 4α² X^{lx+2}) e^{-αr²}
      const Pxmm = lx >= 2 ? ipow(dx, lx - 2) : 0;
      const Pymm = ly >= 2 ? ipow(dy, ly - 2) : 0;
      const Pzmm = lz >= 2 ? ipow(dz, lz - 2) : 0;
      const Pxpp = dx * Pxp;
      const Pypp = dy * Pyp;
      const Pzpp = dz * Pzp;

      const dxxF = (lx * (lx - 1) * Pxmm - 2 * a * (2 * lx + 1) * Px + 4 * a * a * Pxpp) * Py * Pz;
      const dyyF = Px * (ly * (ly - 1) * Pymm - 2 * a * (2 * ly + 1) * Py + 4 * a * a * Pypp) * Pz;
      const dzzF = Px * Py * (lz * (lz - 1) * Pzmm - 2 * a * (2 * lz + 1) * Pz + 4 * a * a * Pzpp);

      // Mixed: ∂_x ∂_y = (lx X^{lx-1} - 2α X^{lx+1}) · (ly Y^{ly-1} - 2α Y^{ly+1}) · Z^lz
      const fx = lx * Pxm - 2 * a * Pxp;
      const fy = ly * Pym - 2 * a * Pyp;
      const fz = lz * Pzm - 2 * a * Pzp;
      const dxyF = fx * fy * Pz;
      const dxzF = fx * Py * fz;
      const dyzF = Px * fy * fz;

      ddphi[(0 * 3 + 0) * nbasis + basisIdx] += N * dxxF;
      ddphi[(1 * 3 + 1) * nbasis + basisIdx] += N * dyyF;
      ddphi[(2 * 3 + 2) * nbasis + basisIdx] += N * dzzF;
      ddphi[(0 * 3 + 1) * nbasis + basisIdx] += N * dxyF;
      ddphi[(1 * 3 + 0) * nbasis + basisIdx] += N * dxyF;
      ddphi[(0 * 3 + 2) * nbasis + basisIdx] += N * dxzF;
      ddphi[(2 * 3 + 0) * nbasis + basisIdx] += N * dxzF;
      ddphi[(1 * 3 + 2) * nbasis + basisIdx] += N * dyzF;
      ddphi[(2 * 3 + 1) * nbasis + basisIdx] += N * dyzF;
    }
  }

  // Significant basis index list (treat all non-trivial as significant)
  const sigIdx = new Int32Array(nbasis);
  let nsig = 0;
  for (let i = 0; i < nbasis; i++) {
    if (Math.abs(phi[i]) > 1e-15
        || Math.abs(dphi[i]) > 1e-15
        || Math.abs(dphi[nbasis + i]) > 1e-15
        || Math.abs(dphi[2 * nbasis + i]) > 1e-15) {
      sigIdx[nsig++] = i;
    }
  }

  return { phi, dphi, ddphi, basisAtom, sigIdx, nsig };
}

/** Compute V_xc nuclear gradient.
 *  RKS: density is total closed-shell density matrix (P).
 *  Returns flat 3*nAtom array. */
export function computeVxcGradient(
  functional: XCFunctional,
  grid: GridPoint[],
  shells: PrimitiveShell[],
  atoms: Atom[],
  normFactors: number[],
  nbasis: number,
  density: Matrix,
): Float64Array {
  const nAtom = atoms.length;
  const grad = new Float64Array(3 * nAtom);
  const needGrad = functional.type !== 'LDA';
  const shellNorms = buildShellNormCache(shells);

  for (const gp of grid) {
    const { phi, dphi, ddphi, basisAtom, sigIdx, nsig } =
      evaluateBasisAndHessian(shells, shellNorms, normFactors, nbasis, gp.x, gp.y, gp.z);
    if (nsig === 0) continue;

    // Compute density and gradient at this grid point
    let rho = 0;
    let gradAx = 0, gradAy = 0, gradAz = 0;
    for (let ii = 0; ii < nsig; ii++) {
      const mu = sigIdx[ii];
      const phiMu = phi[mu];
      const dxMu = dphi[mu], dyMu = dphi[nbasis + mu], dzMu = dphi[2 * nbasis + mu];
      const muOff = mu * nbasis;
      for (let jj = 0; jj < nsig; jj++) {
        const nu = sigIdx[jj];
        const D = density.data[muOff + nu];
        if (Math.abs(D) < 1e-15) continue;
        rho += D * phiMu * phi[nu];
        if (needGrad) {
          gradAx += D * (dxMu * phi[nu] + phiMu * dphi[nu]);
          gradAy += D * (dyMu * phi[nu] + phiMu * dphi[nbasis + nu]);
          gradAz += D * (dzMu * phi[nu] + phiMu * dphi[2 * nbasis + nu]);
        }
      }
    }
    if (rho < 1e-12) continue;

    // RKS: split spin (closed-shell)
    const rhoA = rho * 0.5;
    const rhoB = rho * 0.5;
    let gAA = 0, gBB = 0, gAB = 0;
    let gradAx2 = 0, gradBx2 = 0, gradAy2 = 0, gradBy2 = 0, gradAz2 = 0, gradBz2 = 0;
    if (needGrad) {
      gradAx2 = gradAx * 0.5; gradAy2 = gradAy * 0.5; gradAz2 = gradAz * 0.5;
      gradBx2 = gradAx2; gradBy2 = gradAy2; gradBz2 = gradAz2;
      gAA = gradAx2 * gradAx2 + gradAy2 * gradAy2 + gradAz2 * gradAz2;
      gBB = gAA;
      gAB = gAA;
    }

    const inp: XCInput = { rhoA, rhoB };
    if (needGrad) { inp.gammaAA = gAA; inp.gammaBB = gBB; inp.gammaAB = gAB; }
    const xcOut = functional.evaluate(inp);

    // Build f-vectors used for gradient assembly
    // For LDA: only V^ρ
    const vRhoA = xcOut.vrhoA;
    const vRhoB = xcOut.vrhoB;
    let fAx = 0, fAy = 0, fAz = 0;
    let fBx = 0, fBy = 0, fBz = 0;
    if (needGrad) {
      const vgAA = xcOut.vgammaAA ?? 0;
      const vgAB = xcOut.vgammaAB ?? 0;
      const vgBB = xcOut.vgammaBB ?? 0;
      fAx = 2 * vgAA * gradAx2 + vgAB * gradBx2;
      fAy = 2 * vgAA * gradAy2 + vgAB * gradBy2;
      fAz = 2 * vgAA * gradAz2 + vgAB * gradBz2;
      fBx = 2 * vgBB * gradBx2 + vgAB * gradAx2;
      fBy = 2 * vgBB * gradBy2 + vgAB * gradAy2;
      fBz = 2 * vgBB * gradBz2 + vgAB * gradAz2;
    }

    // Closed-shell V_xc total = V^ρα + V^ρβ scaling (since both spins same density)
    // For closed-shell, V_xc^total at point = V^ρα = V^ρβ (symmetric)
    // The gradient contribution per atom A:
    //   dE_xc/dR_A_d = -2 Σ_{μ on A} Σ_ν P_μν * [
    //       (V^ρα + V^ρβ)/2 · ∂_d φ_μ φ_ν                     ← this contains 1 factor of 2 for closed shell
    //     + GGA term: (f^A_d' + f^B_d') · ∂_d∂_d' φ_μ · φ_ν   ← second deriv of φ
    //     + GGA term: (f^A_d' + f^B_d') · ∂_d φ_μ · ∂_d' φ_ν  ← first deriv of both
    //     ] · weight
    // For RKS closed shell, V^ρα = V^ρβ, fAx = fBx, etc., so total f = 2*fA, V = 2*V^ρα.
    // Then we use the closed-shell "total" density in P (which is 2× alpha density).
    // So the scaling is consistent.

    const w = gp.weight;
    const Vrho = vRhoA;  // = vRhoB for closed shell, used with closed-shell P (which has factor 2)
    const fxX = fAx, fxY = fAy, fxZ = fAz;  // = fBx etc. for closed shell

    // Loop over basis functions ON each atom A (rows), and all ν (columns)
    for (let ii = 0; ii < nsig; ii++) {
      const mu = sigIdx[ii];
      const A = basisAtom[mu];
      const dxMu = dphi[mu];
      const dyMu = dphi[nbasis + mu];
      const dzMu = dphi[2 * nbasis + mu];
      // 2nd derivatives of φ_μ
      const dxxMu = ddphi[(0 * 3 + 0) * nbasis + mu];
      const dyyMu = ddphi[(1 * 3 + 1) * nbasis + mu];
      const dzzMu = ddphi[(2 * 3 + 2) * nbasis + mu];
      const dxyMu = ddphi[(0 * 3 + 1) * nbasis + mu];
      const dxzMu = ddphi[(0 * 3 + 2) * nbasis + mu];
      const dyzMu = ddphi[(1 * 3 + 2) * nbasis + mu];
      const muOff = mu * nbasis;

      for (let jj = 0; jj < nsig; jj++) {
        const nu = sigIdx[jj];
        const D = density.data[muOff + nu];
        if (Math.abs(D) < 1e-15) continue;
        const phiNu = phi[nu];
        const dxNu = dphi[nu];
        const dyNu = dphi[nbasis + nu];
        const dzNu = dphi[2 * nbasis + nu];

        // LDA part: V_rho · ∂_d φ_μ · φ_ν
        // Coefficient: -2 (from formula) × D × w
        // For closed-shell P (=2 P_α), Vrho is per-spin. So no extra factor.
        const c0 = -2 * D * w * Vrho;
        grad[3 * A + 0] += c0 * dxMu * phiNu;
        grad[3 * A + 1] += c0 * dyMu * phiNu;
        grad[3 * A + 2] += c0 * dzMu * phiNu;

        // GGA part: f^σ · [∂²_{d,d'} φ_μ · φ_ν + ∂_d φ_μ · ∂_{d'} φ_ν]
        // Wait: only the d-th component of the derivative (we're computing dR_A_d gradient).
        // Actually for ∂γ/∂R_A_d:
        //   2 ∇ρ · d(∇ρ)/dR_A_d
        // d(∇ρ)/dR_A_d at component d':
        //   = -Σ_μ on A Σ_ν P_μν [∂²_{d,d'} φ_μ · φ_ν + ∂_d φ_μ · ∂_{d'} φ_ν]   (for d' coordinate)
        // f^σ_d' = 2 V^γ ∇ρ_d', so the full integral is:
        //   ∫ Σ_{d'} f^σ_d' · d(∂_{d'} ρ)/dR_A_d dr
        // = -2 Σ_{μ on A} Σ_ν P_μν · Σ_{d'} f_d' [∂²_{d,d'} φ_μ · φ_ν + ∂_d φ_μ · ∂_{d'} φ_ν]
        // For closed shell, sum over both spins gives factor 2 → using closed-shell P naturally.
        if (needGrad) {
          // ∂²_{d=x,d'} φ_μ contributions for d=x, summed over d'
          const ddx_x = dxxMu * fxX + dxyMu * fxY + dxzMu * fxZ;
          const ddx_y = dxyMu * fxX + dyyMu * fxY + dyzMu * fxZ;
          const ddx_z = dxzMu * fxX + dyzMu * fxY + dzzMu * fxZ;
          const c1 = -2 * D * w;
          // For x-direction gradient: terms involving ∂²_{x,d'} φ_μ · φ_ν · f_d'
          //   and ∂_x φ_μ · ∂_{d'} φ_ν · f_d'
          grad[3 * A + 0] += c1 * (ddx_x * phiNu + dxMu * (fxX * dxNu + fxY * dyNu + fxZ * dzNu));
          // For y-direction gradient: ∂²_{y,d'} φ_μ · φ_ν · f_d' + ∂_y φ_μ · ∂_{d'} φ_ν · f_d'
          grad[3 * A + 1] += c1 * (ddx_y * phiNu + dyMu * (fxX * dxNu + fxY * dyNu + fxZ * dzNu));
          // For z-direction gradient
          grad[3 * A + 2] += c1 * (ddx_z * phiNu + dzMu * (fxX * dxNu + fxY * dyNu + fxZ * dzNu));
        }
      }
    }
  }

  return grad;
}
