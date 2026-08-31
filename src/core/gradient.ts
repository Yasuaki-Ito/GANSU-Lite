/** Analytical RHF nuclear gradients via McMurchie-Davidson derivatives.
 *
 *  dE/dR_A = dVnn/dR_A
 *          + Σ_μν D_μν d(T+V)_μν/dR_A
 *          - Σ_μν W_μν dS_μν/dR_A
 *          + Σ_μν (dG_μν/dR_A) D_μν
 *
 *  Key derivative formula (Helgaker 9.3.30):
 *    d/dA_x [K_AB E^{ij}_t] = K_AB [2α E^{i+1,j}_t − i E^{i-1,j}_t]
 */

import type { Atom, PrimitiveShell } from './types';
import { ANGULAR_MOMENTUMS, shellTypeToNumBasis } from './constants';
import { eCoefficients, primitiveNorm } from './integrals1e';
import { boysAll } from './boys';
import { Matrix } from '../linalg/matrix';
import type { XCFunctional } from './xcFunctional';
import type { GridPoint } from './grid';
import { computeVxcGradient } from './xcGradient';

// ─── W-matrix ─────────────────────────────────────────────────────────────────

/** W_ij = 2 Σ_{k∈occ} ε_k C_ik C_jk (energy-weighted density matrix) */
export function computeWMatrix(
  C: Matrix, eps: Float64Array, nocc: number, n: number,
): Matrix {
  const W = new Matrix(n, n);
  for (let i = 0; i < n; i++)
    for (let j = 0; j <= i; j++) {
      let v = 0;
      for (let k = 0; k < nocc; k++) v += 2 * eps[k] * C.get(i, k) * C.get(j, k);
      W.set(i, j, v);
      W.set(j, i, v);
    }
  return W;
}

// ─── Nuclear repulsion gradient ───────────────────────────────────────────────

export function nuclearRepulsionGradient(atoms: Atom[]): Float64Array {
  const grad = new Float64Array(3 * atoms.length);
  for (let a = 0; a < atoms.length; a++)
    for (let b = 0; b < atoms.length; b++) {
      if (a === b) continue;
      const dx = atoms[a].coordinate.x - atoms[b].coordinate.x;
      const dy = atoms[a].coordinate.y - atoms[b].coordinate.y;
      const dz = atoms[a].coordinate.z - atoms[b].coordinate.z;
      const r3 = (dx * dx + dy * dy + dz * dz) ** 1.5;
      const f = atoms[a].atomicNumber * atoms[b].atomicNumber / r3;
      grad[3 * a] -= f * dx;
      grad[3 * a + 1] -= f * dy;
      grad[3 * a + 2] -= f * dz;
    }
  return grad;
}

// ─── Primitive shell grouping ─────────────────────────────────────────────────

interface ShellGrp {
  prims: PrimitiveShell[];
  bi: number;  // basisIndex
  st: number;  // shellType
}

function groupShells(shells: PrimitiveShell[]): ShellGrp[] {
  const g: ShellGrp[] = [];
  let cur: ShellGrp | null = null;
  for (const ps of shells) {
    if (!cur || cur.bi !== ps.basisIndex || cur.st !== ps.shellType) {
      cur = { prims: [], bi: ps.basisIndex, st: ps.shellType };
      g.push(cur);
    }
    cur.prims.push(ps);
  }
  return g;
}

// ─── E-coefficient derivative helpers ─────────────────────────────────────────

/** d/dA_x [K E^{ij}_t] / K = 2α E^{i+1,j}_t − i E^{i-1,j}_t
 *  Returns array of derivatives for t = 0..i+j */
function eCoeffDeriv(
  i: number, j: number, alpha: number, p: number, PA: number, PB: number,
): Float64Array {
  const maxT = i + j;
  const Ep1 = eCoefficients(i + 1, j, p, PA, PB);
  const Em1 = i > 0 ? eCoefficients(i - 1, j, p, PA, PB) : null;
  const d = new Float64Array(maxT + 1);
  for (let t = 0; t <= maxT; t++) {
    d[t] = 2 * alpha * (t < Ep1.length ? Ep1[t] : 0)
      - (Em1 && t < Em1.length ? i * Em1[t] : 0);
  }
  return d;
}

/** d/dB_x [K E^{ij}_t] / K = 2β E^{i,j+1}_t − j E^{i,j-1}_t */
function eCoeffDerivB(
  i: number, j: number, beta: number, p: number, PA: number, PB: number,
): Float64Array {
  const maxT = i + j;
  const Ep1 = eCoefficients(i, j + 1, p, PA, PB);
  const Em1 = j > 0 ? eCoefficients(i, j - 1, p, PA, PB) : null;
  const d = new Float64Array(maxT + 1);
  for (let t = 0; t <= maxT; t++) {
    d[t] = 2 * beta * (t < Ep1.length ? Ep1[t] : 0)
      - (Em1 && t < Em1.length ? j * Em1[t] : 0);
  }
  return d;
}

// ─── R-function (reusable buffer) ─────────────────────────────────────────────

let _rBuf = new Float64Array(0);
let _rDim = 0;

function computeR(maxN: number, eta: number, PCx: number, PCy: number, PCz: number, boys: Float64Array) {
  const dim = maxN + 1, d1 = dim + 1, size = d1 ** 4;
  if (size > _rBuf.length) _rBuf = new Float64Array(size);
  else _rBuf.fill(0, 0, size);
  _rDim = dim;
  const d2 = d1 * d1, d3 = d2 * d1;
  for (let n = 0; n <= maxN; n++) _rBuf[n * d3] = (-2 * eta) ** n * boys[n];
  for (let n = maxN - 1; n >= 0; n--)
    for (let t = 0; t <= maxN - n; t++)
      for (let u = 0; u <= maxN - n - t; u++)
        for (let v = 0; v <= maxN - n - t - u; v++) {
          if (t + u + v === 0) continue;
          if (v > 0) {
            let val = PCz * _rBuf[(n + 1) * d3 + t * d2 + u * d1 + v - 1];
            if (v >= 2) val += (v - 1) * _rBuf[(n + 1) * d3 + t * d2 + u * d1 + v - 2];
            _rBuf[n * d3 + t * d2 + u * d1 + v] = val;
          } else if (u > 0) {
            let val = PCy * _rBuf[(n + 1) * d3 + t * d2 + (u - 1) * d1];
            if (u >= 2) val += (u - 1) * _rBuf[(n + 1) * d3 + t * d2 + (u - 2) * d1];
            _rBuf[n * d3 + t * d2 + u * d1] = val;
          } else {
            let val = PCx * _rBuf[(n + 1) * d3 + (t - 1) * d2];
            if (t >= 2) val += (t - 1) * _rBuf[(n + 1) * d3 + (t - 2) * d2];
            _rBuf[n * d3 + t * d2] = val;
          }
        }
}

function R(t: number, u: number, v: number, n: number): number {
  const d1 = _rDim + 1;
  return _rBuf[n * d1 * d1 * d1 + t * d1 * d1 + u * d1 + v];
}

// ─── 1-electron gradient: overlap (Pulay), kinetic, nuclear attraction ────────

/** Compute 1-electron gradient contributions:
 *  grad += Σ_μν D_μν (dT_μν/dR + dV_μν/dR) − Σ_μν W_μν dS_μν/dR
 *  Uses μ≥ν loop with explicit translational invariance (standard approach). */
export function oneElectronGradient(
  shells: PrimitiveShell[], atoms: Atom[], norms: number[], n: number,
  D: Matrix, W: Matrix,
  components = { pulay: true, kinetic: true, nucAttr: true },
): Float64Array {
  const grad = new Float64Array(3 * atoms.length);
  const groups = groupShells(shells);

  for (let igA = 0; igA < groups.length; igA++) {
    const gA = groups[igA];
    for (let igB = igA; igB < groups.length; igB++) {
      const gB = groups[igB];
      const angA = ANGULAR_MOMENTUMS[gA.st], angB = ANGULAR_MOMENTUMS[gB.st];
      for (let ia = 0; ia < angA.length; ia++) {
        const [lxa, lya, lza] = angA[ia];
        const mu = gA.bi + ia;
        const ibStart = (igA === igB) ? ia : 0;
        for (let ib = ibStart; ib < angB.length; ib++) {
          const [lxb, lyb, lzb] = angB[ib];
          const nu = gB.bi + ib;

          // Symmetry factor: 2 for off-diagonal, 1 for diagonal
          const sym = (mu === nu) ? 1.0 : 2.0;
          const dij = D.get(mu, nu) * sym;
          const wij = W.get(mu, nu) * sym;

          for (const pA of gA.prims) {
            const nA = primitiveNorm(pA.exponent, lxa, lya, lza);
            for (const pB of gB.prims) {
              const nB = primitiveNorm(pB.exponent, lxb, lyb, lzb);
              const alpha = pA.exponent, beta = pB.exponent, p = alpha + beta;
              const Ax = pA.coordinate.x, Ay = pA.coordinate.y, Az = pA.coordinate.z;
              const Bx = pB.coordinate.x, By = pB.coordinate.y, Bz = pB.coordinate.z;
              const Kab = Math.exp(-alpha * beta / p * ((Ax - Bx) ** 2 + (Ay - By) ** 2 + (Az - Bz) ** 2));
              const Px = (alpha * Ax + beta * Bx) / p;
              const Py = (alpha * Ay + beta * By) / p;
              const Pz = (alpha * Az + beta * Bz) / p;
              const coeff = pA.coefficient * pB.coefficient * Kab * nA * nB * norms[mu] * norms[nu];
              const prefS = (Math.PI / p) ** 1.5;

              // E-coefficients and their A-derivatives
              const Ex = eCoefficients(lxa, lxb, p, Px - Ax, Px - Bx);
              const Ey = eCoefficients(lya, lyb, p, Py - Ay, Py - By);
              const Ez = eCoefficients(lza, lzb, p, Pz - Az, Pz - Bz);
              const dEx = eCoeffDeriv(lxa, lxb, alpha, p, Px - Ax, Px - Bx);
              const dEy = eCoeffDeriv(lya, lyb, alpha, p, Py - Ay, Py - By);
              const dEz = eCoeffDeriv(lza, lzb, alpha, p, Pz - Az, Pz - Bz);

              const atomA = pA.atomIndex, atomB = pB.atomIndex;

              // ── Overlap Pulay: -W_μν dS_μν/dA, dS/dB = -dS/dA ──
              if (components.pulay) {
                const dSx = prefS * dEx[0] * Ey[0] * Ez[0];
                const dSy = prefS * Ex[0] * dEy[0] * Ez[0];
                const dSz = prefS * Ex[0] * Ey[0] * dEz[0];
                const ws = -coeff * wij;
                grad[3 * atomA]     += ws * dSx;
                grad[3 * atomA + 1] += ws * dSy;
                grad[3 * atomA + 2] += ws * dSz;
                grad[3 * atomB]     -= ws * dSx;
                grad[3 * atomB + 1] -= ws * dSy;
                grad[3 * atomB + 2] -= ws * dSz;
              }

              // ── Kinetic: D_μν dT_μν/dA, dT/dB = -dT/dA ──
              if (components.kinetic) {
                const dTx = kineticDeriv(lxa, lya, lza, lxb, lyb, lzb, beta, p, alpha,
                  Px - Ax, Py - Ay, Pz - Az, Px - Bx, Py - By, Pz - Bz, 0);
                const dTy = kineticDeriv(lxa, lya, lza, lxb, lyb, lzb, beta, p, alpha,
                  Px - Ax, Py - Ay, Pz - Az, Px - Bx, Py - By, Pz - Bz, 1);
                const dTz = kineticDeriv(lxa, lya, lza, lxb, lyb, lzb, beta, p, alpha,
                  Px - Ax, Py - Ay, Pz - Az, Px - Bx, Py - By, Pz - Bz, 2);
                const wd = coeff * dij * prefS;
                grad[3 * atomA]     += wd * dTx;
                grad[3 * atomA + 1] += wd * dTy;
                grad[3 * atomA + 2] += wd * dTz;
                grad[3 * atomB]     -= wd * dTx;
                grad[3 * atomB + 1] -= wd * dTy;
                grad[3 * atomB + 2] -= wd * dTz;
              }

              // ── Nuclear attraction: D_μν dV_μν/dA ──
              if (components.nucAttr) {
              // dV/dC = HF (Hellmann-Feynman force on nucleus)
              const prefV = 2 * Math.PI / p;
              for (let iC = 0; iC < atoms.length; iC++) {
                const Z = atoms[iC].atomicNumber;
                const Cx = atoms[iC].coordinate.x, Cy = atoms[iC].coordinate.y, Cz = atoms[iC].coordinate.z;
                const PCx = Px - Cx, PCy = Py - Cy, PCz = Pz - Cz;
                const maxN = lxa + lya + lza + lxb + lyb + lzb;
                const boys = boysAll(maxN + 1, p * (PCx * PCx + PCy * PCy + PCz * PCz));
                computeR(maxN + 1, p, PCx, PCy, PCz, boys);

                let gAx = 0, gAy = 0, gAz = 0;
                let gCx = 0, gCy = 0, gCz = 0;
                const maxT = lxa + lxb, maxU = lya + lyb, maxV = lza + lzb;
                for (let t = 0; t <= maxT; t++)
                  for (let u = 0; u <= maxU; u++)
                    for (let v = 0; v <= maxV; v++) {
                      const e = Ex[t] * Ey[u] * Ez[v];
                      const Rv = R(t, u, v, 0);
                      // Pulay A: D_t = 2α E^{i+1,j}_t − i E^{i-1,j}_t (Helgaker expansion)
                      gAx += dEx[t] * Ey[u] * Ez[v] * Rv;
                      gAy += Ex[t] * dEy[u] * Ez[v] * Rv;
                      gAz += Ex[t] * Ey[u] * dEz[v] * Rv;
                      // Chain rule boundary: (α/p) R_{shifted} only at max index per direction
                      if (t === maxT) gAx += e * (alpha / p) * R(t + 1, u, v, 0);
                      if (u === maxU) gAy += e * (alpha / p) * R(t, u + 1, v, 0);
                      if (v === maxV) gAz += e * (alpha / p) * R(t, u, v + 1, 0);
                      // Hellmann-Feynman on C (no boundary issue — pure R derivative)
                      gCx -= e * R(t + 1, u, v, 0);
                      gCy -= e * R(t, u + 1, v, 0);
                      gCz -= e * R(t, u, v + 1, 0);
                    }

                const wv = -Z * coeff * dij * prefV;
                // Pulay on atom A
                grad[3 * atomA]     += wv * gAx;
                grad[3 * atomA + 1] += wv * gAy;
                grad[3 * atomA + 2] += wv * gAz;
                // Hellmann-Feynman on nuclear center C
                grad[3 * iC]     += wv * gCx;
                grad[3 * iC + 1] += wv * gCy;
                grad[3 * iC + 2] += wv * gCz;
                // 3-center translational invariance: dV/dB = -(dV/dA + dV/dC)
                grad[3 * atomB]     -= wv * (gAx + gCx);
                grad[3 * atomB + 1] -= wv * (gAy + gCy);
                grad[3 * atomB + 2] -= wv * (gAz + gCz);
              }
              } // end nucAttr
            }
          }
        }
      }
    }
  }
  return grad;
}

/** Kinetic energy derivative in direction dir (0=x,1=y,2=z) w.r.t. center A.
 *  Returns the derivative of (π/p)^{-3/2} T_{ij} (caller multiplies by prefS). */
function kineticDeriv(
  lxa: number, lya: number, lza: number, lxb: number, lyb: number, lzb: number,
  beta: number, p: number, alpha: number,
  PAx: number, PAy: number, PAz: number, PBx: number, PBy: number, PBz: number,
  dir: number,
): number {
  // Helper: E^{a,b}_0 overlap coeff
  function E0(ax: number, bx: number, PA: number, PB: number): number {
    if (ax < 0 || bx < 0) return 0;
    return eCoefficients(ax, bx, p, PA, PB)[0];
  }
  // Helper: dE^{a,b}_0/dA
  function dE0(ax: number, bx: number, PA: number, PB: number): number {
    if (bx < 0) return 0;
    return 2 * alpha * E0(ax + 1, bx, PA, PB) - (ax > 0 ? ax * E0(ax - 1, bx, PA, PB) : 0);
  }

  const PA = [PAx, PAy, PAz], PB = [PBx, PBy, PBz];
  const la = [lxa, lya, lza], lb = [lxb, lyb, lzb];

  // 3D overlap product and its derivative in direction dir
  function S3d(bShift: number[], useDerivDir: boolean): number {
    let val = 1;
    for (let d = 0; d < 3; d++) {
      const a = la[d], b = lb[d] + bShift[d];
      if (d === dir && useDerivDir) val *= dE0(a, b, PA[d], PB[d]);
      else val *= E0(a, b, PA[d], PB[d]);
    }
    return val;
  }

  // T = Σ_{d=x,y,z} [β(2lb_d+1)S − 2β²S(lb_d+2) − lb_d(lb_d-1)/2 S(lb_d-2)]
  // dT/dA_dir = same with S → dS/dA_dir for direction `dir`
  let result = 0;
  for (let d = 0; d < 3; d++) {
    const shift = [0, 0, 0];
    const l = lb[d];
    // β(2l+1) S
    result += beta * (2 * l + 1) * S3d(shift, true);
    // -2β² S(l+2)
    shift[d] = 2;
    result -= 2 * beta * beta * S3d(shift, true);
    shift[d] = 0;
    // -l(l-1)/2 S(l-2)
    if (l >= 2) {
      shift[d] = -2;
      result -= 0.5 * l * (l - 1) * S3d(shift, true);
      shift[d] = 0;
    }
  }
  return result;
}

// ─── 2-electron gradient ──────────────────────────────────────────────────────

/** Compute 2-electron gradient:
 *  dE_2e/dR = 0.5 Σ_μνλσ Γ_μνλσ d(μν|λσ)/dR
 *  where Γ_μνλσ = D_μν D_λσ − 0.5 D_μλ D_νσ (for RHF)
 *
 *  Uses bra-ket symmetry (AB|CD) = (CD|AB) to halve computation.
 *  Computes explicit derivatives for all 4 centers (A, B, C, D). */
export function twoElectronGradient(
  shells: PrimitiveShell[], norms: number[], n: number, D: Matrix,
  hfFraction: number = 1.0,
): Float64Array {
  let nAtom = 0;
  for (const ps of shells) if (ps.atomIndex >= nAtom) nAtom = ps.atomIndex + 1;
  const grad = new Float64Array(3 * nAtom);
  const groups = groupShells(shells);
  const nGrp = groups.length;

  for (let igA = 0; igA < nGrp; igA++) {
    const gA = groups[igA];
    for (let igB = 0; igB < nGrp; igB++) {
      const gB = groups[igB];
      const igAB = igA * nGrp + igB;
      for (let igC = 0; igC < nGrp; igC++) {
        const gC = groups[igC];
        for (let igD = 0; igD < nGrp; igD++) {
          const gD = groups[igD];
          const igCD = igC * nGrp + igD;
          // Bra-ket symmetry: only compute AB >= CD (upper triangle)
          if (igAB < igCD) continue;
          const braKetSym = igAB === igCD ? 1.0 : 2.0;
          const angA = ANGULAR_MOMENTUMS[gA.st], angB = ANGULAR_MOMENTUMS[gB.st];
          const angC = ANGULAR_MOMENTUMS[gC.st], angD = ANGULAR_MOMENTUMS[gD.st];

          for (let ia = 0; ia < angA.length; ia++) {
            const [lxa, lya, lza] = angA[ia]; const mu = gA.bi + ia;
            for (let ib = 0; ib < angB.length; ib++) {
              const [lxb, lyb, lzb] = angB[ib]; const nu = gB.bi + ib;
              for (let ic = 0; ic < angC.length; ic++) {
                const [lxc, lyc, lzc] = angC[ic]; const lam = gC.bi + ic;
                for (let id = 0; id < angD.length; id++) {
                  const [lxd, lyd, lzd] = angD[id]; const sig = gD.bi + id;

                  const gamma = D.get(mu, nu) * D.get(lam, sig)
                    - hfFraction * 0.5 * D.get(mu, lam) * D.get(nu, sig);
                  if (Math.abs(gamma) < 1e-15) continue;

                  for (const pA of gA.prims) {
                    const nA = primitiveNorm(pA.exponent, lxa, lya, lza);
                    for (const pB of gB.prims) {
                      const nB = primitiveNorm(pB.exponent, lxb, lyb, lzb);
                      const alpha = pA.exponent, beta = pB.exponent, p = alpha + beta;
                      const Ax = pA.coordinate.x, Ay = pA.coordinate.y, Az = pA.coordinate.z;
                      const Bx = pB.coordinate.x, By = pB.coordinate.y, Bz = pB.coordinate.z;
                      const Kab = Math.exp(-alpha * beta / p * ((Ax - Bx) ** 2 + (Ay - By) ** 2 + (Az - Bz) ** 2));
                      const Px = (alpha * Ax + beta * Bx) / p;
                      const Py = (alpha * Ay + beta * By) / p;
                      const Pz = (alpha * Az + beta * Bz) / p;

                      const ExAB = eCoefficients(lxa, lxb, p, Px - Ax, Px - Bx);
                      const EyAB = eCoefficients(lya, lyb, p, Py - Ay, Py - By);
                      const EzAB = eCoefficients(lza, lzb, p, Pz - Az, Pz - Bz);
                      // Bra A-derivatives
                      const dExA = eCoeffDeriv(lxa, lxb, alpha, p, Px - Ax, Px - Bx);
                      const dEyA = eCoeffDeriv(lya, lyb, alpha, p, Py - Ay, Py - By);
                      const dEzA = eCoeffDeriv(lza, lzb, alpha, p, Pz - Az, Pz - Bz);
                      // Bra B-derivatives
                      const dExB = eCoeffDerivB(lxa, lxb, beta, p, Px - Ax, Px - Bx);
                      const dEyB = eCoeffDerivB(lya, lyb, beta, p, Py - Ay, Py - By);
                      const dEzB = eCoeffDerivB(lza, lzb, beta, p, Pz - Az, Pz - Bz);

                      for (const pC of gC.prims) {
                        const nC = primitiveNorm(pC.exponent, lxc, lyc, lzc);
                        for (const pD of gD.prims) {
                          const nD = primitiveNorm(pD.exponent, lxd, lyd, lzd);
                          const gam = pC.exponent, del = pD.exponent, q = gam + del;
                          const Cx = pC.coordinate.x, Cy = pC.coordinate.y, Cz = pC.coordinate.z;
                          const DDx = pD.coordinate.x, DDy = pD.coordinate.y, DDz = pD.coordinate.z;
                          const Kcd = Math.exp(-gam * del / q * ((Cx - DDx) ** 2 + (Cy - DDy) ** 2 + (Cz - DDz) ** 2));
                          const Qx = (gam * Cx + del * DDx) / q;
                          const Qy = (gam * Cy + del * DDy) / q;
                          const Qz = (gam * Cz + del * DDz) / q;
                          const RPQx = Px - Qx, RPQy = Py - Qy, RPQz = Pz - Qz;
                          const eta = p * q / (p + q);
                          const maxN = lxa + lya + lza + lxb + lyb + lzb + lxc + lyc + lzc + lxd + lyd + lzd;

                          const boys = boysAll(maxN + 1, eta * (RPQx * RPQx + RPQy * RPQy + RPQz * RPQz));
                          computeR(maxN + 1, eta, RPQx, RPQy, RPQz, boys);

                          const ExCD = eCoefficients(lxc, lxd, q, Qx - Cx, Qx - DDx);
                          const EyCD = eCoefficients(lyc, lyd, q, Qy - Cy, Qy - DDy);
                          const EzCD = eCoefficients(lzc, lzd, q, Qz - Cz, Qz - DDz);
                          // Ket C-derivatives
                          const dExC = eCoeffDeriv(lxc, lxd, gam, q, Qx - Cx, Qx - DDx);
                          const dEyC = eCoeffDeriv(lyc, lyd, gam, q, Qy - Cy, Qy - DDy);
                          const dEzC = eCoeffDeriv(lzc, lzd, gam, q, Qz - Cz, Qz - DDz);

                          const primCoeff = pA.coefficient * pB.coefficient * pC.coefficient * pD.coefficient
                            * Kab * Kcd * nA * nB * nC * nD * norms[mu] * norms[nu] * norms[lam] * norms[sig];
                          const pre2e = 2 * Math.PI ** 2.5 / (p * q * Math.sqrt(p + q));

                          // Accumulate derivatives for all 4 centers
                          let gAx = 0, gAy = 0, gAz = 0;
                          let gBx = 0, gBy = 0, gBz = 0;
                          let gCx = 0, gCy = 0, gCz = 0;
                          const mT1 = lxa + lxb, mU1 = lya + lyb, mV1 = lza + lzb;
                          const mT2 = lxc + lxd, mU2 = lyc + lyd, mV2 = lzc + lzd;

                          for (let t1 = 0; t1 <= mT1; t1++)
                            for (let u1 = 0; u1 <= mU1; u1++)
                              for (let v1 = 0; v1 <= mV1; v1++) {
                                const e1 = ExAB[t1] * EyAB[u1] * EzAB[v1];
                                for (let t2 = 0; t2 <= mT2; t2++)
                                  for (let u2 = 0; u2 <= mU2; u2++)
                                    for (let v2 = 0; v2 <= mV2; v2++) {
                                      const e2 = ExCD[t2] * EyCD[u2] * EzCD[v2];
                                      const sign = ((t2 + u2 + v2) & 1) ? -1 : 1;
                                      const Rv = sign * R(t1 + t2, u1 + u2, v1 + v2, 0);
                                      const e12 = e1 * e2;

                                      // d/dA: Helgaker expansion + boundary chain rule through P
                                      gAx += (dExA[t1] * EyAB[u1] * EzAB[v1]) * e2 * Rv;
                                      gAy += (ExAB[t1] * dEyA[u1] * EzAB[v1]) * e2 * Rv;
                                      gAz += (ExAB[t1] * EyAB[u1] * dEzA[v1]) * e2 * Rv;
                                      if (t1 === mT1) gAx += e12 * (alpha / p) * sign * R(t1 + t2 + 1, u1 + u2, v1 + v2, 0);
                                      if (u1 === mU1) gAy += e12 * (alpha / p) * sign * R(t1 + t2, u1 + u2 + 1, v1 + v2, 0);
                                      if (v1 === mV1) gAz += e12 * (alpha / p) * sign * R(t1 + t2, u1 + u2, v1 + v2 + 1, 0);

                                      // d/dB: same with β
                                      gBx += (dExB[t1] * EyAB[u1] * EzAB[v1]) * e2 * Rv;
                                      gBy += (ExAB[t1] * dEyB[u1] * EzAB[v1]) * e2 * Rv;
                                      gBz += (ExAB[t1] * EyAB[u1] * dEzB[v1]) * e2 * Rv;
                                      if (t1 === mT1) gBx += e12 * (beta / p) * sign * R(t1 + t2 + 1, u1 + u2, v1 + v2, 0);
                                      if (u1 === mU1) gBy += e12 * (beta / p) * sign * R(t1 + t2, u1 + u2 + 1, v1 + v2, 0);
                                      if (v1 === mV1) gBz += e12 * (beta / p) * sign * R(t1 + t2, u1 + u2, v1 + v2 + 1, 0);

                                      // d/dC: Helgaker ket expansion + boundary chain rule through Q
                                      gCx += e1 * (dExC[t2] * EyCD[u2] * EzCD[v2]) * Rv;
                                      gCy += e1 * (ExCD[t2] * dEyC[u2] * EzCD[v2]) * Rv;
                                      gCz += e1 * (ExCD[t2] * EyCD[u2] * dEzC[v2]) * Rv;
                                      if (t2 === mT2) gCx -= e12 * (gam / q) * sign * R(t1 + t2 + 1, u1 + u2, v1 + v2, 0);
                                      if (u2 === mU2) gCy -= e12 * (gam / q) * sign * R(t1 + t2, u1 + u2 + 1, v1 + v2, 0);
                                      if (v2 === mV2) gCz -= e12 * (gam / q) * sign * R(t1 + t2, u1 + u2, v1 + v2 + 1, 0);
                                    }
                              }

                          // d/dD via 4-center TI: d/dD = -(d/dA + d/dB + d/dC)
                          const w = 0.5 * braKetSym * gamma * primCoeff * pre2e;
                          const atomA = pA.atomIndex, atomB = pB.atomIndex;
                          const atomC = pC.atomIndex, atomD = pD.atomIndex;
                          grad[3 * atomA]     += w * gAx;
                          grad[3 * atomA + 1] += w * gAy;
                          grad[3 * atomA + 2] += w * gAz;
                          grad[3 * atomB]     += w * gBx;
                          grad[3 * atomB + 1] += w * gBy;
                          grad[3 * atomB + 2] += w * gBz;
                          grad[3 * atomC]     += w * gCx;
                          grad[3 * atomC + 1] += w * gCy;
                          grad[3 * atomC + 2] += w * gCz;
                          grad[3 * atomD]     -= w * (gAx + gBx + gCx);
                          grad[3 * atomD + 1] -= w * (gAy + gBy + gCy);
                          grad[3 * atomD + 2] -= w * (gAz + gBz + gCz);
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
  return grad;
}

// ─── Skeleton Fock derivative matrix dG/dR ───────────────────────────────────

/** Compute dG_μν/dR for ALL perturbation directions at once.
 *  Returns array of 3*nAtom matrices (one per perturbation direction).
 *  Single pass over shell quartets — much faster than calling per-direction. */
export function computeAllSkeletonDGMatrices(
  shells: PrimitiveShell[], norms: number[], n: number, D: Matrix, nAtom: number,
): Matrix[] {
  const ndim = 3 * nAtom;
  const dGs: Matrix[] = [];
  for (let i = 0; i < ndim; i++) dGs.push(new Matrix(n, n));
  const groups = groupShells(shells);

  for (const gA of groups)
    for (const gB of groups)
      for (const gC of groups)
        for (const gD of groups) {
          const angA = ANGULAR_MOMENTUMS[gA.st], angB = ANGULAR_MOMENTUMS[gB.st];
          const angC = ANGULAR_MOMENTUMS[gC.st], angD = ANGULAR_MOMENTUMS[gD.st];

          for (let ia = 0; ia < angA.length; ia++) {
            const [lxa, lya, lza] = angA[ia]; const mu = gA.bi + ia;
            for (let ib = 0; ib < angB.length; ib++) {
              const [lxb, lyb, lzb] = angB[ib]; const nu = gB.bi + ib;
              for (let ic = 0; ic < angC.length; ic++) {
                const [lxc, lyc, lzc] = angC[ic]; const lam = gC.bi + ic;
                for (let id = 0; id < angD.length; id++) {
                  const [lxd, lyd, lzd] = angD[id]; const sig = gD.bi + id;

                  // For dG: each quartet contributes to Coulomb and Exchange
                  const Dls = D.get(lam, sig);
                  const Dns = D.get(nu, sig);
                  if (Math.abs(Dls) < 1e-15 && Math.abs(Dns) < 1e-15) continue;

                  for (const pA of gA.prims) {
                    const nA = primitiveNorm(pA.exponent, lxa, lya, lza);
                    for (const pB of gB.prims) {
                      const nB = primitiveNorm(pB.exponent, lxb, lyb, lzb);
                      const alpha = pA.exponent, beta = pB.exponent, p = alpha + beta;
                      const Ax = pA.coordinate.x, Ay = pA.coordinate.y, Az = pA.coordinate.z;
                      const Bx = pB.coordinate.x, By = pB.coordinate.y, Bz = pB.coordinate.z;
                      const Kab = Math.exp(-alpha * beta / p * ((Ax-Bx)**2 + (Ay-By)**2 + (Az-Bz)**2));
                      const Px = (alpha*Ax + beta*Bx)/p, Py = (alpha*Ay + beta*By)/p, Pz = (alpha*Az + beta*Bz)/p;

                      const ExAB = eCoefficients(lxa, lxb, p, Px-Ax, Px-Bx);
                      const EyAB = eCoefficients(lya, lyb, p, Py-Ay, Py-By);
                      const EzAB = eCoefficients(lza, lzb, p, Pz-Az, Pz-Bz);
                      const dExA = eCoeffDeriv(lxa, lxb, alpha, p, Px-Ax, Px-Bx);
                      const dEyA = eCoeffDeriv(lya, lyb, alpha, p, Py-Ay, Py-By);
                      const dEzA = eCoeffDeriv(lza, lzb, alpha, p, Pz-Az, Pz-Bz);
                      const dExB = eCoeffDerivB(lxa, lxb, beta, p, Px-Ax, Px-Bx);
                      const dEyB = eCoeffDerivB(lya, lyb, beta, p, Py-Ay, Py-By);
                      const dEzB = eCoeffDerivB(lza, lzb, beta, p, Pz-Az, Pz-Bz);

                      for (const pC of gC.prims) {
                        const nC = primitiveNorm(pC.exponent, lxc, lyc, lzc);
                        for (const pD of gD.prims) {
                          const nD = primitiveNorm(pD.exponent, lxd, lyd, lzd);
                          const gam = pC.exponent, del = pD.exponent, q = gam + del;
                          const Cx = pC.coordinate.x, Cy = pC.coordinate.y, Cz = pC.coordinate.z;
                          const DDx = pD.coordinate.x, DDy = pD.coordinate.y, DDz = pD.coordinate.z;
                          const Kcd = Math.exp(-gam * del / q * ((Cx-DDx)**2 + (Cy-DDy)**2 + (Cz-DDz)**2));
                          const Qx = (gam*Cx + del*DDx)/q, Qy = (gam*Cy + del*DDy)/q, Qz = (gam*Cz + del*DDz)/q;
                          const RPQx = Px-Qx, RPQy = Py-Qy, RPQz = Pz-Qz;
                          const eta = p * q / (p + q);
                          const maxN = lxa+lya+lza+lxb+lyb+lzb+lxc+lyc+lzc+lxd+lyd+lzd;

                          const boys = boysAll(maxN + 1, eta * (RPQx*RPQx + RPQy*RPQy + RPQz*RPQz));
                          computeR(maxN + 1, eta, RPQx, RPQy, RPQz, boys);

                          const ExCD = eCoefficients(lxc, lxd, q, Qx-Cx, Qx-DDx);
                          const EyCD = eCoefficients(lyc, lyd, q, Qy-Cy, Qy-DDy);
                          const EzCD = eCoefficients(lzc, lzd, q, Qz-Cz, Qz-DDz);
                          const dExC = eCoeffDeriv(lxc, lxd, gam, q, Qx-Cx, Qx-DDx);
                          const dEyC = eCoeffDeriv(lyc, lyd, gam, q, Qy-Cy, Qy-DDy);
                          const dEzC = eCoeffDeriv(lzc, lzd, gam, q, Qz-Cz, Qz-DDz);

                          const primCoeff = pA.coefficient * pB.coefficient * pC.coefficient * pD.coefficient
                            * Kab * Kcd * nA * nB * nC * nD * norms[mu] * norms[nu] * norms[lam] * norms[sig];
                          const pre2e = 2 * Math.PI ** 2.5 / (p * q * Math.sqrt(p + q));

                          // Compute d(μν|λσ)/dR for ALL 4 centers × 3 directions simultaneously
                          const mT1 = lxa+lxb, mU1 = lya+lyb, mV1 = lza+lzb;
                          const mT2 = lxc+lxd, mU2 = lyc+lyd, mV2 = lzc+lzd;
                          // derivs[center][dir]: 4 centers × 3 directions
                          const derivs = [[0,0,0],[0,0,0],[0,0,0],[0,0,0]];

                          const dExD = eCoeffDerivB(lxc,lxd,del,q,Qx-Cx,Qx-DDx);
                          const dEyD = eCoeffDerivB(lyc,lyd,del,q,Qy-Cy,Qy-DDy);
                          const dEzD = eCoeffDerivB(lzc,lzd,del,q,Qz-Cz,Qz-DDz);

                          for (let t1 = 0; t1 <= mT1; t1++)
                            for (let u1 = 0; u1 <= mU1; u1++)
                              for (let v1 = 0; v1 <= mV1; v1++) {
                                const e1 = ExAB[t1]*EyAB[u1]*EzAB[v1];
                                for (let t2 = 0; t2 <= mT2; t2++)
                                  for (let u2 = 0; u2 <= mU2; u2++)
                                    for (let v2 = 0; v2 <= mV2; v2++) {
                                      const e2 = ExCD[t2]*EyCD[u2]*EzCD[v2];
                                      const sign = ((t2+u2+v2)&1)?-1:1;
                                      const Rv = sign*R(t1+t2,u1+u2,v1+v2,0);
                                      const e12 = e1*e2;
                                      const Rx = sign*R(t1+t2+1,u1+u2,v1+v2,0);
                                      const Ry = sign*R(t1+t2,u1+u2+1,v1+v2,0);
                                      const Rz = sign*R(t1+t2,u1+u2,v1+v2+1,0);

                                      // Center A (x,y,z)
                                      derivs[0][0] += (dExA[t1]*EyAB[u1]*EzAB[v1])*e2*Rv;
                                      derivs[0][1] += (ExAB[t1]*dEyA[u1]*EzAB[v1])*e2*Rv;
                                      derivs[0][2] += (ExAB[t1]*EyAB[u1]*dEzA[v1])*e2*Rv;
                                      if(t1===mT1) derivs[0][0] += e12*(alpha/p)*Rx;
                                      if(u1===mU1) derivs[0][1] += e12*(alpha/p)*Ry;
                                      if(v1===mV1) derivs[0][2] += e12*(alpha/p)*Rz;
                                      // Center B
                                      derivs[1][0] += (dExB[t1]*EyAB[u1]*EzAB[v1])*e2*Rv;
                                      derivs[1][1] += (ExAB[t1]*dEyB[u1]*EzAB[v1])*e2*Rv;
                                      derivs[1][2] += (ExAB[t1]*EyAB[u1]*dEzB[v1])*e2*Rv;
                                      if(t1===mT1) derivs[1][0] += e12*(beta/p)*Rx;
                                      if(u1===mU1) derivs[1][1] += e12*(beta/p)*Ry;
                                      if(v1===mV1) derivs[1][2] += e12*(beta/p)*Rz;
                                      // Center C
                                      derivs[2][0] += e1*(dExC[t2]*EyCD[u2]*EzCD[v2])*Rv;
                                      derivs[2][1] += e1*(ExCD[t2]*dEyC[u2]*EzCD[v2])*Rv;
                                      derivs[2][2] += e1*(ExCD[t2]*EyCD[u2]*dEzC[v2])*Rv;
                                      if(t2===mT2) derivs[2][0] -= e12*(gam/q)*Rx;
                                      if(u2===mU2) derivs[2][1] -= e12*(gam/q)*Ry;
                                      if(v2===mV2) derivs[2][2] -= e12*(gam/q)*Rz;
                                      // Center D
                                      derivs[3][0] += e1*(dExD[t2]*EyCD[u2]*EzCD[v2])*Rv;
                                      derivs[3][1] += e1*(ExCD[t2]*dEyD[u2]*EzCD[v2])*Rv;
                                      derivs[3][2] += e1*(ExCD[t2]*EyCD[u2]*dEzD[v2])*Rv;
                                      if(t2===mT2) derivs[3][0] -= e12*(del/q)*Rx;
                                      if(u2===mU2) derivs[3][1] -= e12*(del/q)*Ry;
                                      if(v2===mV2) derivs[3][2] -= e12*(del/q)*Rz;
                                    }
                              }

                          const coef = primCoeff * pre2e;
                          const atomIndices = [pA.atomIndex, pB.atomIndex, pC.atomIndex, pD.atomIndex];
                          // Accumulate to dG matrices for all perturbation directions
                          for (let ci = 0; ci < 4; ci++) {
                            const ai = atomIndices[ci];
                            for (let dd = 0; dd < 3; dd++) {
                              const dv = derivs[ci][dd] * coef;
                              if (Math.abs(dv) < 1e-18) continue;
                              const idx = 3 * ai + dd;
                              // Coulomb
                              dGs[idx].set(mu, nu, dGs[idx].get(mu, nu) + Dls * dv);
                              // Exchange
                              dGs[idx].set(mu, lam, dGs[idx].get(mu, lam) - 0.5 * Dns * dv);
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

  return dGs;
}

/** Compute dG for a single perturbation direction (wrapper). */
export function computeSkeletonDGMatrix(
  shells: PrimitiveShell[], norms: number[], n: number, D: Matrix,
  pertAtom: number, pertDir: number,
): Matrix {
  let nAtom = 0;
  for (const ps of shells) if (ps.atomIndex >= nAtom) nAtom = ps.atomIndex + 1;
  const all = computeAllSkeletonDGMatrices(shells, norms, n, D, nAtom);
  return all[3 * pertAtom + pertDir];
}

// ─── Naive 1-electron gradient (debug version) ──────────────────────────────

/** Naive 1e gradient: full μ,ν loop, explicit A/B derivatives, NO translational
 *  invariance. For debugging — should give identical results to oneElectronGradient. */
export function oneElectronGradientNaive(
  shells: PrimitiveShell[], atoms: Atom[], norms: number[], _n: number,
  D: Matrix, W: Matrix,
): Float64Array {
  const grad = new Float64Array(3 * atoms.length);
  const groups = groupShells(shells);

  for (const gA of groups)
    for (const gB of groups) {
      const angA = ANGULAR_MOMENTUMS[gA.st], angB = ANGULAR_MOMENTUMS[gB.st];
      for (let ia = 0; ia < angA.length; ia++) {
        const [lxa, lya, lza] = angA[ia];
        const mu = gA.bi + ia;
        for (let ib = 0; ib < angB.length; ib++) {
          const [lxb, lyb, lzb] = angB[ib];
          const nu = gB.bi + ib;
          const dij = D.get(mu, nu);
          const wij = W.get(mu, nu);

          for (const pA of gA.prims) {
            const nA = primitiveNorm(pA.exponent, lxa, lya, lza);
            for (const pB of gB.prims) {
              const nB = primitiveNorm(pB.exponent, lxb, lyb, lzb);
              const alpha = pA.exponent, beta = pB.exponent, p = alpha + beta;
              const Ax = pA.coordinate.x, Ay = pA.coordinate.y, Az = pA.coordinate.z;
              const Bx = pB.coordinate.x, By = pB.coordinate.y, Bz = pB.coordinate.z;
              const Kab = Math.exp(-alpha * beta / p * ((Ax - Bx) ** 2 + (Ay - By) ** 2 + (Az - Bz) ** 2));
              const Px = (alpha * Ax + beta * Bx) / p;
              const Py = (alpha * Ay + beta * By) / p;
              const Pz = (alpha * Az + beta * Bz) / p;
              const coeff = pA.coefficient * pB.coefficient * Kab * nA * nB * norms[mu] * norms[nu];
              const prefS = (Math.PI / p) ** 1.5;

              const Ex = eCoefficients(lxa, lxb, p, Px - Ax, Px - Bx);
              const Ey = eCoefficients(lya, lyb, p, Py - Ay, Py - By);
              const Ez = eCoefficients(lza, lzb, p, Pz - Az, Pz - Bz);
              const dExA = eCoeffDeriv(lxa, lxb, alpha, p, Px - Ax, Px - Bx);
              const dEyA = eCoeffDeriv(lya, lyb, alpha, p, Py - Ay, Py - By);
              const dEzA = eCoeffDeriv(lza, lzb, alpha, p, Pz - Az, Pz - Bz);
              const dExB = eCoeffDerivB(lxa, lxb, beta, p, Px - Ax, Px - Bx);
              const dEyB = eCoeffDerivB(lya, lyb, beta, p, Py - Ay, Py - By);
              const dEzB = eCoeffDerivB(lza, lzb, beta, p, Pz - Az, Pz - Bz);

              const atomA = pA.atomIndex, atomB = pB.atomIndex;

              // ── Overlap Pulay: -W dS/dR ──
              {
                const ws = -coeff * wij;
                // dS/dA
                grad[3 * atomA]     += ws * prefS * dExA[0] * Ey[0] * Ez[0];
                grad[3 * atomA + 1] += ws * prefS * Ex[0] * dEyA[0] * Ez[0];
                grad[3 * atomA + 2] += ws * prefS * Ex[0] * Ey[0] * dEzA[0];
                // dS/dB (explicit, no TI)
                grad[3 * atomB]     += ws * prefS * dExB[0] * Ey[0] * Ez[0];
                grad[3 * atomB + 1] += ws * prefS * Ex[0] * dEyB[0] * Ez[0];
                grad[3 * atomB + 2] += ws * prefS * Ex[0] * Ey[0] * dEzB[0];
              }

              // ── Kinetic: D dT/dR ──
              {
                const wd = coeff * dij * prefS;
                const dTAx = kineticDeriv(lxa, lya, lza, lxb, lyb, lzb, beta, p, alpha,
                  Px - Ax, Py - Ay, Pz - Az, Px - Bx, Py - By, Pz - Bz, 0);
                const dTAy = kineticDeriv(lxa, lya, lza, lxb, lyb, lzb, beta, p, alpha,
                  Px - Ax, Py - Ay, Pz - Az, Px - Bx, Py - By, Pz - Bz, 1);
                const dTAz = kineticDeriv(lxa, lya, lza, lxb, lyb, lzb, beta, p, alpha,
                  Px - Ax, Py - Ay, Pz - Az, Px - Bx, Py - By, Pz - Bz, 2);
                // dT/dA
                grad[3 * atomA]     += wd * dTAx;
                grad[3 * atomA + 1] += wd * dTAy;
                grad[3 * atomA + 2] += wd * dTAz;
                // dT/dB = -dT/dA (2-center TI, always valid)
                grad[3 * atomB]     -= wd * dTAx;
                grad[3 * atomB + 1] -= wd * dTAy;
                grad[3 * atomB + 2] -= wd * dTAz;
              }

              // ── Nuclear attraction: D dV/dR ──
              {
                const prefV = 2 * Math.PI / p;
                for (let iC = 0; iC < atoms.length; iC++) {
                  const Z = atoms[iC].atomicNumber;
                  const Cx = atoms[iC].coordinate.x, Cy = atoms[iC].coordinate.y, Cz = atoms[iC].coordinate.z;
                  const PCx = Px - Cx, PCy = Py - Cy, PCz = Pz - Cz;
                  const maxN = lxa + lya + lza + lxb + lyb + lzb;
                  const boys = boysAll(maxN + 1, p * (PCx * PCx + PCy * PCy + PCz * PCz));
                  computeR(maxN + 1, p, PCx, PCy, PCz, boys);

                  let gAx = 0, gAy = 0, gAz = 0;
                  let gBx = 0, gBy = 0, gBz = 0;
                  let gCx = 0, gCy = 0, gCz = 0;
                  const maxT = lxa + lxb, maxU = lya + lyb, maxV = lza + lzb;
                  for (let t = 0; t <= maxT; t++)
                    for (let u = 0; u <= maxU; u++)
                      for (let v = 0; v <= maxV; v++) {
                        const e = Ex[t] * Ey[u] * Ez[v];
                        const Rv = R(t, u, v, 0);
                        // Pulay A: Helgaker expansion + boundary chain rule
                        gAx += dExA[t] * Ey[u] * Ez[v] * Rv;
                        gAy += Ex[t] * dEyA[u] * Ez[v] * Rv;
                        gAz += Ex[t] * Ey[u] * dEzA[v] * Rv;
                        if (t === maxT) gAx += e * (alpha / p) * R(t + 1, u, v, 0);
                        if (u === maxU) gAy += e * (alpha / p) * R(t, u + 1, v, 0);
                        if (v === maxV) gAz += e * (alpha / p) * R(t, u, v + 1, 0);
                        // Pulay B: same structure with β
                        gBx += dExB[t] * Ey[u] * Ez[v] * Rv;
                        gBy += Ex[t] * dEyB[u] * Ez[v] * Rv;
                        gBz += Ex[t] * Ey[u] * dEzB[v] * Rv;
                        if (t === maxT) gBx += e * (beta / p) * R(t + 1, u, v, 0);
                        if (u === maxU) gBy += e * (beta / p) * R(t, u + 1, v, 0);
                        if (v === maxV) gBz += e * (beta / p) * R(t, u, v + 1, 0);
                        // HF on C (pure R derivative, no boundary issue)
                        gCx -= e * R(t + 1, u, v, 0);
                        gCy -= e * R(t, u + 1, v, 0);
                        gCz -= e * R(t, u, v + 1, 0);
                      }

                  const wv = -Z * coeff * dij * prefV;
                  grad[3 * atomA]     += wv * gAx;
                  grad[3 * atomA + 1] += wv * gAy;
                  grad[3 * atomA + 2] += wv * gAz;
                  grad[3 * atomB]     += wv * gBx;
                  grad[3 * atomB + 1] += wv * gBy;
                  grad[3 * atomB + 2] += wv * gBz;
                  grad[3 * iC]        += wv * gCx;
                  grad[3 * iC + 1]    += wv * gCy;
                  grad[3 * iC + 2]    += wv * gCz;
                }
              }
            }
          }
        }
      }
    }
  return grad;
}

// ─── Full RHF gradient ────────────────────────────────────────────────────────

export interface DFTGradientContext {
  functional: XCFunctional;
  grid: GridPoint[];
}

export function computeRHFGradient(
  shells: PrimitiveShell[], atoms: Atom[], norms: number[], n: number, nocc: number,
  D: Matrix, C: Matrix, eps: Float64Array,
  components = { nuc: true, pulay: true, kinetic: true, nucAttr: true, twoEl: true },
  dftCtx?: DFTGradientContext,
): { total: Float64Array; nuc: Float64Array; oneEl: Float64Array; twoEl: Float64Array; xc?: Float64Array } {
  const nAt = atoms.length;
  const total = new Float64Array(3 * nAt);

  const nuc = components.nuc ? nuclearRepulsionGradient(atoms) : new Float64Array(3 * nAt);
  const W = computeWMatrix(C, eps, nocc, n);
  const oneEl = (components.pulay || components.kinetic || components.nucAttr)
    ? oneElectronGradient(shells, atoms, norms, n, D, W, components)
    : new Float64Array(3 * nAt);
  // For DFT (KS), exchange in the 2e term is scaled by exact-exchange fraction
  // (hfFraction=0 for pure DFT skips Coulomb-only formation of K). For RSH, the
  // long-range exchange portion (β K_LR) is currently not included — deferred.
  const hfFrac = dftCtx ? dftCtx.functional.exactExchangeFraction : 1.0;
  const twoEl = components.twoEl ? twoElectronGradient(shells, norms, n, D, hfFrac) : new Float64Array(3 * nAt);

  // V_xc gradient (DFT only)
  let xc: Float64Array | undefined;
  if (dftCtx) {
    xc = computeVxcGradient(dftCtx.functional, dftCtx.grid, shells, atoms, norms, n, D);
  }

  // Debug: check translational invariance (warn only for HF-only components;
  // xc has expected violation from skipped Becke weight derivatives)
  for (const [name, g] of [['nuc', nuc], ['1e', oneEl], ['2e', twoEl]] as const) {
    let sx = 0, sy = 0, sz = 0;
    for (let a = 0; a < nAt; a++) { sx += g[3*a]; sy += g[3*a+1]; sz += g[3*a+2]; }
    if (Math.abs(sx) + Math.abs(sy) + Math.abs(sz) > 1e-6)
      console.warn(`${name} translational invariance violation: [${sx.toExponential(3)}, ${sy.toExponential(3)}, ${sz.toExponential(3)}]`);
  }

  for (let i = 0; i < total.length; i++) {
    total[i] = nuc[i] + oneEl[i] + twoEl[i];
    if (xc) total[i] += xc[i];
  }
  return { total, nuc, oneEl, twoEl, xc };
}
