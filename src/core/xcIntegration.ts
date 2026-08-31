/** Evaluate density on grid, compute V_xc matrix and E_xc energy */

import { PrimitiveShell } from './types';
import { ANGULAR_MOMENTUMS, shellTypeToNumBasis } from './constants';
import { primitiveNorm } from './integrals1e';
import { GridPoint } from './grid';
import { XCFunctional, XCInput } from './xcFunctional';
import { isXCWasmAvailable, computeXCWasm, packShells } from './eriWasm';

/** Map functional name to WASM func_id */
const FUNC_ID_MAP: Record<string, number> = {
  SVWN: 0, BLYP: 1, PBE: 2, B3LYP: 3,
};

/** Integer power for small exponents (avoids Math.pow overhead) */
function ipow(x: number, n: number): number {
  if (n === 0) return 1;
  if (n === 1) return x;
  if (n === 2) return x * x;
  if (n === 3) return x * x * x;
  return Math.pow(x, n);
}

/** Pre-computed primitive normalization factors and screening data per shell */
interface ShellNormCache {
  pNorms: number[];  // primitiveNorm for each basis function in shell
  r2Max: number;     // distance threshold: skip if r² > r2Max (exp(-α*r²) < 1e-30)
}

function buildShellNormCache(shells: PrimitiveShell[]): ShellNormCache[] {
  const LOG_THRESH = 69.08; // -ln(1e-30) ≈ 69.08
  return shells.map(shell => {
    const angMom = ANGULAR_MOMENTUMS[shell.shellType];
    const nBF = shellTypeToNumBasis(shell.shellType);
    const pNorms: number[] = [];
    for (let bf = 0; bf < nBF; bf++) {
      const [lx, ly, lz] = angMom[bf];
      pNorms.push(primitiveNorm(shell.exponent, lx, ly, lz));
    }
    const r2Max = LOG_THRESH / shell.exponent;
    return { pNorms, r2Max };
  });
}

/**
 * Evaluate basis functions at a grid point, returning values and the list of
 * significant (non-negligible) basis function indices.
 */
export { buildShellNormCache, evaluateBasisFunctionsScreened };

function evaluateBasisFunctionsScreened(
  shells: PrimitiveShell[],
  shellNorms: ShellNormCache[],
  normFactors: number[],
  nbasis: number,
  px: number, py: number, pz: number,
  needGradient: boolean,
): {
  phi: Float64Array; dphiDx: Float64Array; dphiDy: Float64Array; dphiDz: Float64Array;
  sigIdx: Int32Array; nsig: number;
} {
  const phi = new Float64Array(nbasis);
  const dphiDx = new Float64Array(nbasis);
  const dphiDy = new Float64Array(nbasis);
  const dphiDz = new Float64Array(nbasis);

  for (let s = 0; s < shells.length; s++) {
    const shell = shells[s];
    const cache = shellNorms[s];
    const dx = px - shell.coordinate.x;
    const dy = py - shell.coordinate.y;
    const dz = pz - shell.coordinate.z;
    const r2 = dx * dx + dy * dy + dz * dz;
    if (r2 > cache.r2Max) continue;  // pre-screen: skip exp() entirely
    const expVal = shell.coefficient * Math.exp(-shell.exponent * r2);

    const angMom = ANGULAR_MOMENTUMS[shell.shellType];
    const nBF = shellTypeToNumBasis(shell.shellType);

    for (let bf = 0; bf < nBF; bf++) {
      const [lx, ly, lz] = angMom[bf];
      const basisIdx = shell.basisIndex + bf;
      const coeff = normFactors[basisIdx] * cache.pNorms[bf] * expVal;

      const angular = ipow(dx, lx) * ipow(dy, ly) * ipow(dz, lz);
      phi[basisIdx] += coeff * angular;

      if (needGradient) {
        const expAlpha = -2 * shell.exponent;
        const pxl = ipow(dx, lx);
        const pyl = ipow(dy, ly);
        const pzl = ipow(dz, lz);

        let dAngDx = expAlpha * dx * angular;
        if (lx > 0) dAngDx += lx * ipow(dx, lx - 1) * pyl * pzl;

        let dAngDy = expAlpha * dy * angular;
        if (ly > 0) dAngDy += ly * pxl * ipow(dy, ly - 1) * pzl;

        let dAngDz = expAlpha * dz * angular;
        if (lz > 0) dAngDz += lz * pxl * pyl * ipow(dz, lz - 1);

        dphiDx[basisIdx] += coeff * dAngDx;
        dphiDy[basisIdx] += coeff * dAngDy;
        dphiDz[basisIdx] += coeff * dAngDz;
      }
    }
  }

  // Build list of significant basis functions
  const sigIdx = new Int32Array(nbasis);
  let nsig = 0;
  const thresh = 1e-15;
  for (let i = 0; i < nbasis; i++) {
    if (Math.abs(phi[i]) > thresh ||
        (needGradient && (
          Math.abs(dphiDx[i]) > thresh ||
          Math.abs(dphiDy[i]) > thresh ||
          Math.abs(dphiDz[i]) > thresh))) {
      sigIdx[nsig++] = i;
    }
  }

  return { phi, dphiDx, dphiDy, dphiDz, sigIdx, nsig };
}

export interface XCResult {
  exc: number;
  vxcA: Float64Array;
  vxcB?: Float64Array;
  numElectrons: number;
}

/**
 * Compute XC energy and V_xc matrix.
 * RKS: densityA = total density matrix, densityB = null.
 * UKS: densityA = alpha density, densityB = beta density.
 */
export function computeXC(
  functional: XCFunctional,
  grid: GridPoint[],
  shells: PrimitiveShell[],
  normFactors: number[],
  nbasis: number,
  densityA: Float64Array,
  densityB: Float64Array | null = null,
): XCResult {
  // Try WASM path first (WASM does not yet support meta-GGA)
  const funcId = FUNC_ID_MAP[functional.name];
  if (functional.type !== 'meta-GGA' && isXCWasmAvailable() && funcId !== undefined) {
    const shellsFlat = packShells(shells);
    const normF64 = new Float64Array(normFactors);
    const dbFlat = densityB ?? new Float64Array(0);
    // Pack grid into flat array [x0,y0,z0,w0, x1,y1,z1,w1, ...]
    const gridFlat = new Float64Array(grid.length * 4);
    for (let i = 0; i < grid.length; i++) {
      const off = i * 4;
      gridFlat[off] = grid[i].x;
      gridFlat[off + 1] = grid[i].y;
      gridFlat[off + 2] = grid[i].z;
      gridFlat[off + 3] = grid[i].weight;
    }
    const needGrad = functional.type !== 'LDA';
    const result = computeXCWasm(shellsFlat, normF64, densityA, dbFlat, gridFlat, nbasis, funcId, needGrad);
    if (result) {
      return result;
    }
  }

  const needGrad = functional.type !== 'LDA';
  const needTau = functional.type === 'meta-GGA';
  const vxcA = new Float64Array(nbasis * nbasis);
  let vxcB: Float64Array | undefined;
  if (densityB) vxcB = new Float64Array(nbasis * nbasis);

  // Pre-compute shell normalization factors
  const shellNorms = buildShellNormCache(shells);

  let exc = 0;
  let numElec = 0;

  for (const gp of grid) {
    const { phi, dphiDx, dphiDy, dphiDz, sigIdx, nsig } =
      evaluateBasisFunctionsScreened(shells, shellNorms, normFactors, nbasis,
                                      gp.x, gp.y, gp.z, needGrad);

    if (nsig === 0) continue;

    // Compute density using only significant basis functions
    let rhoA = 0, rhoB = 0;
    let gradAx = 0, gradAy = 0, gradAz = 0;
    let gradBx = 0, gradBy = 0, gradBz = 0;
    let tauA = 0, tauB = 0;

    for (let ii = 0; ii < nsig; ii++) {
      const mu = sigIdx[ii];
      const phiMu = phi[mu];
      const dxMu = dphiDx[mu], dyMu = dphiDy[mu], dzMu = dphiDz[mu];
      const muOff = mu * nbasis;

      for (let jj = 0; jj < nsig; jj++) {
        const nu = sigIdx[jj];
        const dA = densityA[muOff + nu];
        const pp = phiMu * phi[nu];
        rhoA += dA * pp;

        if (densityB) {
          rhoB += densityB[muOff + nu] * pp;
        }

        if (needGrad) {
          const dpx = dxMu * phi[nu] + phiMu * dphiDx[nu];
          const dpy = dyMu * phi[nu] + phiMu * dphiDy[nu];
          const dpz = dzMu * phi[nu] + phiMu * dphiDz[nu];
          gradAx += dA * dpx;
          gradAy += dA * dpy;
          gradAz += dA * dpz;
          if (needTau) {
            // ∇φ_μ · ∇φ_ν, τ accumulator before factor of 1/2
            const dpdot = dxMu * dphiDx[nu] + dyMu * dphiDy[nu] + dzMu * dphiDz[nu];
            tauA += dA * dpdot;
            if (densityB) tauB += densityB[muOff + nu] * dpdot;
          }
          if (densityB) {
            const dB = densityB[muOff + nu];
            gradBx += dB * dpx;
            gradBy += dB * dpy;
            gradBz += dB * dpz;
          }
        }
      }
    }

    // RKS: densityA is total density, split into α and β
    if (!densityB) {
      rhoB = rhoA * 0.5;
      rhoA *= 0.5;
      if (needGrad) {
        gradBx = gradAx * 0.5; gradBy = gradAy * 0.5; gradBz = gradAz * 0.5;
        gradAx *= 0.5; gradAy *= 0.5; gradAz *= 0.5;
      }
      if (needTau) { tauB = tauA * 0.5; tauA *= 0.5; }
    }
    // Convert τ-accumulator (Σ P_μν ∇φ_μ·∇φ_ν = 2τ_spin) to τ definition.
    if (needTau) { tauA *= 0.5; tauB *= 0.5; }

    // Clamp to avoid negative densities from numerical noise in the density matrix
    if (rhoA < 0) rhoA = 0;
    if (rhoB < 0) rhoB = 0;
    const rho = rhoA + rhoB;
    if (rho < 1e-20) continue;

    const inp: XCInput = { rhoA, rhoB };
    if (needGrad) {
      inp.gammaAA = gradAx * gradAx + gradAy * gradAy + gradAz * gradAz;
      inp.gammaBB = gradBx * gradBx + gradBy * gradBy + gradBz * gradBz;
      inp.gammaAB = gradAx * gradBx + gradAy * gradBy + gradAz * gradBz;
    }
    if (needTau) { inp.tauA = tauA; inp.tauB = tauB; }

    numElec += rho * gp.weight;
    const xcOut = functional.evaluate(inp);
    exc += xcOut.exc * rho * gp.weight;

    // Pre-compute GGA vector: fα = 2*vgammaAA*∇ρα + vgammaAB*∇ρβ
    let fAx = 0, fAy = 0, fAz = 0;
    let fBx = 0, fBy = 0, fBz = 0;
    if (needGrad) {
      const vgAA = xcOut.vgammaAA ?? 0;
      const vgAB = xcOut.vgammaAB ?? 0;
      const vgBB = xcOut.vgammaBB ?? 0;
      fAx = 2 * vgAA * gradAx + vgAB * gradBx;
      fAy = 2 * vgAA * gradAy + vgAB * gradBy;
      fAz = 2 * vgAA * gradAz + vgAB * gradBz;
      if (densityB) {
        fBx = 2 * vgBB * gradBx + vgAB * gradAx;
        fBy = 2 * vgBB * gradBy + vgAB * gradAy;
        fBz = 2 * vgBB * gradBz + vgAB * gradAz;
      }
    }

    // Build V_xc matrix using only significant basis functions
    for (let ii = 0; ii < nsig; ii++) {
      const mu = sigIdx[ii];
      const phiMu = phi[mu];
      const dxMu = dphiDx[mu], dyMu = dphiDy[mu], dzMu = dphiDz[mu];

      for (let jj = ii; jj < nsig; jj++) {
        const nu = sigIdx[jj];
        const phiNu = phi[nu];

        let valA = xcOut.vrhoA * phiMu * phiNu;
        let valB = densityB ? xcOut.vrhoB * phiMu * phiNu : 0;

        if (needGrad) {
          const nabla_x = dxMu * phiNu + phiMu * dphiDx[nu];
          const nabla_y = dyMu * phiNu + phiMu * dphiDy[nu];
          const nabla_z = dzMu * phiNu + phiMu * dphiDz[nu];
          valA += fAx * nabla_x + fAy * nabla_y + fAz * nabla_z;
          if (densityB) {
            valB += fBx * nabla_x + fBy * nabla_y + fBz * nabla_z;
          }
        }

        if (needTau) {
          // GKS contribution: F_xc^μν += (1/2) v_τ ∇φ_μ · ∇φ_ν
          const dpdot = dxMu * dphiDx[nu] + dyMu * dphiDy[nu] + dzMu * dphiDz[nu];
          valA += 0.5 * (xcOut.vtauA ?? 0) * dpdot;
          if (densityB) valB += 0.5 * (xcOut.vtauB ?? 0) * dpdot;
        }

        const wA = valA * gp.weight;
        vxcA[mu * nbasis + nu] += wA;
        if (mu !== nu) vxcA[nu * nbasis + mu] += wA;

        if (vxcB) {
          const wB = valB * gp.weight;
          vxcB[mu * nbasis + nu] += wB;
          if (mu !== nu) vxcB[nu * nbasis + mu] += wB;
        }
      }
    }
  }

  return { exc, vxcA, vxcB, numElectrons: numElec };
}
