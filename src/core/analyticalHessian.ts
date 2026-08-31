/** Analytical RHF Hessian — ported from GANSU C++ (GPU).
 *
 *  Components:
 *    d²E/dR²  = d²Vnn/dR²                    (nuclear repulsion)
 *             + Σ_μν D_μν d²H_μν/dR²         (1-electron: S, T, V)
 *             - Σ_μν W_μν d²S_μν/dR²         (overlap Pulay)
 *             + Σ_μνλσ Γ_μνλσ d²(μν|λσ)/dR²  (2-electron)
 *             + CPHF response terms
 *
 *  Uses McMurchie-Davidson E-coefficient shift formulas for derivatives.
 *  Overlap 1D integral: S(l1,l2) = E^{l1,l2}_0(p, PA, PB)
 */

import type { Atom, PrimitiveShell } from './types';
import { ANGULAR_MOMENTUMS } from './constants';
import { eCoefficients, primitiveNorm } from './integrals1e';
import { boysAll } from './boys';
import type { ERIStored } from './eri';
import { Matrix } from '../linalg/matrix';
import { computeWMatrix, computeAllSkeletonDGMatrices } from './gradient';
import { transformERItoMO, solveCPHF } from './cphf';
import { computeOneElectronIntegrals } from './integrals1e';
import {
  isWasmAvailable, compute2eHessianWasm, transformERItoMOWasm,
  solveCPHFWasm, packShellsExt,
} from './eriWasm';

// ── Overlap 1D integral ──────────────────────────────────────────────

/** S(l1, l2, α, β, AB) = E^{l1,l2}_0 * K * (π/p)^{1/2} per direction.
 *  Here we return just E^{l1,l2}_0 (caller handles prefactor). */
function S1d(l1: number, l2: number, p: number, PA: number, PB: number): number {
  if (l1 < 0 || l2 < 0) return 0;
  return eCoefficients(l1, l2, p, PA, PB)[0];
}

// ── Kinetic 1D integral ──────────────────────────────────────────────

/** T(l1, l2) = β(2l2+1)S(l1,l2) - 2β²S(l1,l2+2) - l2(l2-1)/2 S(l1,l2-2) */
function T1d(l1: number, l2: number, beta: number, p: number, PA: number, PB: number): number {
  return beta * (2 * l2 + 1) * S1d(l1, l2, p, PA, PB)
    - 2 * beta * beta * S1d(l1, l2 + 2, p, PA, PB)
    - (l2 >= 2 ? 0.5 * l2 * (l2 - 1) * S1d(l1, l2 - 2, p, PA, PB) : 0);
}

// ── Shift formulas for derivatives ───────────────────────────────────

/** First derivative of S(l1,l2) w.r.t. center A:
 *  dS/dA = 2α·S(l1+1,l2) - l1·S(l1-1,l2) */
function dS_A(l1: number, l2: number, alpha: number, p: number, PA: number, PB: number): number {
  return 2 * alpha * S1d(l1 + 1, l2, p, PA, PB)
    - (l1 > 0 ? l1 * S1d(l1 - 1, l2, p, PA, PB) : 0);
}

/** Second derivative of S(l1,l2) w.r.t. center A (same direction):
 *  d²S/dA² = 4α²·S(l1+2,l2) - 2α(2l1+1)·S(l1,l2) + l1(l1-1)·S(l1-2,l2) */
function d2S_AA(l1: number, l2: number, alpha: number, p: number, PA: number, PB: number): number {
  return 4 * alpha * alpha * S1d(l1 + 2, l2, p, PA, PB)
    - 2 * alpha * (2 * l1 + 1) * S1d(l1, l2, p, PA, PB)
    + (l1 >= 2 ? l1 * (l1 - 1) * S1d(l1 - 2, l2, p, PA, PB) : 0);
}

/** First derivative of T(l1,l2) w.r.t. center A:
 *  Same shift formula applied to T instead of S */
function dT_A(l1: number, l2: number, alpha: number, beta: number, p: number, PA: number, PB: number): number {
  return 2 * alpha * T1d(l1 + 1, l2, beta, p, PA, PB)
    - (l1 > 0 ? l1 * T1d(l1 - 1, l2, beta, p, PA, PB) : 0);
}

/** Second derivative of T(l1,l2) w.r.t. center A (same direction) */
function d2T_AA(l1: number, l2: number, alpha: number, beta: number, p: number, PA: number, PB: number): number {
  return 4 * alpha * alpha * T1d(l1 + 2, l2, beta, p, PA, PB)
    - 2 * alpha * (2 * l1 + 1) * T1d(l1, l2, beta, p, PA, PB)
    + (l1 >= 2 ? l1 * (l1 - 1) * T1d(l1 - 2, l2, beta, p, PA, PB) : 0);
}

// ── Helper: shell grouping ───────────────────────────────────────────

interface ShellGrp { prims: PrimitiveShell[]; bi: number; st: number; }

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

// ── Nuclear repulsion Hessian ────────────────────────────────────────

function nuclearRepulsionHessian(atoms: Atom[]): Float64Array {
  const nAt = atoms.length, ndim = 3 * nAt;
  const H = new Float64Array(ndim * ndim);
  for (let a = 0; a < nAt; a++)
    for (let b = a + 1; b < nAt; b++) {
      const Za = atoms[a].atomicNumber, Zb = atoms[b].atomicNumber;
      const dx = atoms[a].coordinate.x - atoms[b].coordinate.x;
      const dy = atoms[a].coordinate.y - atoms[b].coordinate.y;
      const dz = atoms[a].coordinate.z - atoms[b].coordinate.z;
      const r2 = dx * dx + dy * dy + dz * dz;
      const r = Math.sqrt(r2);
      const r5 = r2 * r2 * r;
      const f = Za * Zb / r5;
      const dd = [dx, dy, dz];

      for (let d1 = 0; d1 < 3; d1++)
        for (let d2 = 0; d2 < 3; d2++) {
          const val = f * (3 * dd[d1] * dd[d2] - (d1 === d2 ? r2 : 0));
          // AB block (off-diagonal atoms): negative
          H[(3 * a + d1) * ndim + (3 * b + d2)] -= val;
          H[(3 * b + d2) * ndim + (3 * a + d1)] -= val;
          // AA and BB blocks: positive
          H[(3 * a + d1) * ndim + (3 * a + d2)] += val;
          H[(3 * b + d1) * ndim + (3 * b + d2)] += val;
        }
    }
  return H;
}

// ── Overlap & Kinetic Hessian (1-electron skeleton) ──────────────────

/** Compute skeleton 1e Hessian: -W·d²S + D·d²T + D·d²V
 *  For overlap and kinetic only (no nuclear attraction — that needs R function). */
function oneElectronHessianSK(
  shells: PrimitiveShell[], norms: number[], n: number,
  D: Matrix, W: Matrix, nAtom: number,
): Float64Array {
  const ndim = 3 * nAtom;
  const H = new Float64Array(ndim * ndim);
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
          const sym = (mu === nu) ? 1.0 : 2.0;
          const dij = D.get(mu, nu) * sym;
          const wij = -W.get(mu, nu) * sym; // note: -W for overlap Pulay

          for (const pA of gA.prims) {
            const nA = primitiveNorm(pA.exponent, lxa, lya, lza);
            for (const pB of gB.prims) {
              const nB = primitiveNorm(pB.exponent, lxb, lyb, lzb);
              const alpha = pA.exponent, beta = pB.exponent, p = alpha + beta;
              const Ax = pA.coordinate.x, Ay = pA.coordinate.y, Az = pA.coordinate.z;
              const Bx = pB.coordinate.x, By = pB.coordinate.y, Bz = pB.coordinate.z;
              const Kab = Math.exp(-alpha * beta / p * ((Ax-Bx)**2 + (Ay-By)**2 + (Az-Bz)**2));
              const Px = (alpha * Ax + beta * Bx) / p;
              const Py = (alpha * Ay + beta * By) / p;
              const Pz = (alpha * Az + beta * Bz) / p;
              const coeff = pA.coefficient * pB.coefficient * Kab * nA * nB * norms[mu] * norms[nu];
              const prefS = (Math.PI / p) ** 1.5;

              const PAx = Px - Ax, PAy = Py - Ay, PAz = Pz - Az;
              const PBx = Px - Bx, PBy = Py - By, PBz = Pz - Bz;

              const la = [lxa, lya, lza], lb = [lxb, lyb, lzb];
              const PA = [PAx, PAy, PAz], PB = [PBx, PBy, PBz];

              // Base 1D integrals in each direction
              const Sx = [S1d(la[0], lb[0], p, PA[0], PB[0]),
                          S1d(la[1], lb[1], p, PA[1], PB[1]),
                          S1d(la[2], lb[2], p, PA[2], PB[2])];

              // First derivatives dS/dA per direction
              const dSA = [dS_A(la[0], lb[0], alpha, p, PA[0], PB[0]),
                           dS_A(la[1], lb[1], alpha, p, PA[1], PB[1]),
                           dS_A(la[2], lb[2], alpha, p, PA[2], PB[2])];

              // Second derivatives d²S/dA² per direction
              const d2SA = [d2S_AA(la[0], lb[0], alpha, p, PA[0], PB[0]),
                            d2S_AA(la[1], lb[1], alpha, p, PA[1], PB[1]),
                            d2S_AA(la[2], lb[2], alpha, p, PA[2], PB[2])];

              // Kinetic 1D
              const Tx = [T1d(la[0], lb[0], beta, p, PA[0], PB[0]),
                          T1d(la[1], lb[1], beta, p, PA[1], PB[1]),
                          T1d(la[2], lb[2], beta, p, PA[2], PB[2])];
              const dTA = [dT_A(la[0], lb[0], alpha, beta, p, PA[0], PB[0]),
                           dT_A(la[1], lb[1], alpha, beta, p, PA[1], PB[1]),
                           dT_A(la[2], lb[2], alpha, beta, p, PA[2], PB[2])];
              const d2TA = [d2T_AA(la[0], lb[0], alpha, beta, p, PA[0], PB[0]),
                            d2T_AA(la[1], lb[1], alpha, beta, p, PA[1], PB[1]),
                            d2T_AA(la[2], lb[2], alpha, beta, p, PA[2], PB[2])];

              const iA = pA.atomIndex, iB = pB.atomIndex;

              // AA block (same direction and cross direction)
              const hess_AA = new Float64Array(9);
              for (let d1 = 0; d1 < 3; d1++) {
                for (let d2 = d1; d2 < 3; d2++) {
                  let valS: number, valT: number;
                  if (d1 === d2) {
                    // Same direction: d²S uses d2SA, others use Sx
                    valS = d2SA[d1];
                    valT = d2TA[d1];
                    for (let k = 0; k < 3; k++) {
                      if (k !== d1) { valS *= Sx[k]; valT = valT * Sx[k] + (k < d1 ? 0 : 0); }
                    }
                    // Actually need full product: d²T/dA²_d1 = d2TA[d1]*Sx[other1]*Sx[other2] for S part
                    // and for T: d2TA[d1]*Sx[other1]*Sx[other2] + ... cross terms
                    // Full 3D kinetic: T_3d = Tx*Sy*Sz + Sx*Ty*Sz + Sx*Sy*Tz
                    // d²T_3d/dA_d1² = d2TA[d1]*Sy*Sz + Sx*d2TA[d1_of_T]*Sz + ...
                    // This is getting complex. Let me use a cleaner approach.
                  }
                  // Skip complex inline, use direct product below
                }
              }

              // Direct computation of d²(S_3D)/dA_d1 dA_d2 and d²(T_3D)/dA_d1 dA_d2
              // S_3D = Sx * Sy * Sz
              // T_3D = Tx*Sy*Sz + Sx*Ty*Sz + Sx*Sy*Tz
              for (let d1 = 0; d1 < 3; d1++) {
                for (let d2 = d1; d2 < 3; d2++) {
                  let valS = 0, valT = 0;

                  if (d1 === d2) {
                    // d²S_3D/dA_d² = d2SA[d] * Π(Sx[k], k≠d)
                    let prodOther = 1;
                    for (let k = 0; k < 3; k++) if (k !== d1) prodOther *= Sx[k];
                    valS = d2SA[d1] * prodOther;

                    // d²T_3D/dA_d² = Σ_comp d²(Tcomp)/dA_d²
                    // Tcomp = Tx*Sy*Sz (comp=0), Sx*Ty*Sz (comp=1), Sx*Sy*Tz (comp=2)
                    for (let comp = 0; comp < 3; comp++) {
                      let term = 1;
                      for (let k = 0; k < 3; k++) {
                        if (k === comp) {
                          term *= (k === d1) ? d2TA[k] : Tx[k];
                        } else {
                          term *= (k === d1) ? d2SA[k] : Sx[k];
                        }
                      }
                      valT += term;
                    }
                  } else {
                    // d²S_3D/dA_d1 dA_d2 = dSA[d1] * dSA[d2] * Sx[d3]  (d3 = other)
                    const d3 = 3 - d1 - d2; // 0+1+2=3
                    valS = dSA[d1] * dSA[d2] * Sx[d3];

                    // d²T_3D/dA_d1 dA_d2 = Σ_comp (product rule)
                    for (let comp = 0; comp < 3; comp++) {
                      let term = 1;
                      for (let k = 0; k < 3; k++) {
                        const isT = (k === comp);
                        const isDeriv1 = (k === d1);
                        const isDeriv2 = (k === d2);
                        if (isDeriv1 && isDeriv2) {
                          // both derivatives on same direction (only when d1==d2, handled above)
                          term *= isT ? d2TA[k] : d2SA[k];
                        } else if (isDeriv1) {
                          term *= isT ? dTA[k] : dSA[k];
                        } else if (isDeriv2) {
                          term *= isT ? dTA[k] : dSA[k];
                        } else {
                          term *= isT ? Tx[k] : Sx[k];
                        }
                      }
                      valT += term;
                    }
                  }

                  const w_total = coeff * prefS * (wij * valS + dij * valT);
                  if (w_total === 0) continue;

                  // Write to Hessian using translational invariance
                  const gA1 = 3 * iA + d1, gA2 = 3 * iA + d2;
                  const gB1 = 3 * iB + d1, gB2 = 3 * iB + d2;

                  H[gA1 * ndim + gA2] += w_total;
                  H[gB1 * ndim + gB2] += w_total;
                  H[gA1 * ndim + gB2] -= w_total;
                  H[gB1 * ndim + gA2] -= w_total;

                  if (d1 !== d2) {
                    H[gA2 * ndim + gA1] += w_total;
                    H[gB2 * ndim + gB1] += w_total;
                    H[gA2 * ndim + gB1] -= w_total;
                    H[gB2 * ndim + gA1] -= w_total;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return H;
}

// ── Nuclear attraction Hessian ────────────────────────────────────────

/** Compute the full 3D nuclear attraction integral V(la,lb) for given angular momenta.
 *  Uses McMurchie-Davidson: V = -Z (2π/p) K Σ_{tuv} E^x_t E^y_u E^z_v R_{tuv}^0 */
function computeNAI3D(
  la: number[], lb: number[], alpha: number, beta: number,
  A: number[], B: number[], C: number[],
): number {
  const p = alpha + beta;
  const P = [(alpha * A[0] + beta * B[0]) / p, (alpha * A[1] + beta * B[1]) / p, (alpha * A[2] + beta * B[2]) / p];
  const PA = [P[0] - A[0], P[1] - A[1], P[2] - A[2]];
  const PB = [P[0] - B[0], P[1] - B[1], P[2] - B[2]];
  const PC = [P[0] - C[0], P[1] - C[1], P[2] - C[2]];

  const Ex = eCoefficients(la[0], lb[0], p, PA[0], PB[0]);
  const Ey = eCoefficients(la[1], lb[1], p, PA[1], PB[1]);
  const Ez = eCoefficients(la[2], lb[2], p, PA[2], PB[2]);

  const maxN = la[0] + lb[0] + la[1] + lb[1] + la[2] + lb[2];
  const PC2 = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];
  const boys = boysAll(maxN, p * PC2);

  // R function via downward recurrence
  const dim = maxN + 1, d1 = dim + 1, size = d1 ** 4;
  const Rbuf = new Float64Array(size);
  const d2 = d1 * d1, d3 = d2 * d1;
  for (let nn = 0; nn <= maxN; nn++) Rbuf[nn * d3] = (-2 * p) ** nn * boys[nn];
  for (let nn = maxN - 1; nn >= 0; nn--)
    for (let t = 0; t <= maxN - nn; t++)
      for (let u = 0; u <= maxN - nn - t; u++)
        for (let v = 0; v <= maxN - nn - t - u; v++) {
          if (t + u + v === 0) continue;
          if (v > 0) { let val = PC[2] * Rbuf[(nn+1)*d3+t*d2+u*d1+v-1]; if (v >= 2) val += (v-1)*Rbuf[(nn+1)*d3+t*d2+u*d1+v-2]; Rbuf[nn*d3+t*d2+u*d1+v] = val; }
          else if (u > 0) { let val = PC[1] * Rbuf[(nn+1)*d3+t*d2+(u-1)*d1]; if (u >= 2) val += (u-1)*Rbuf[(nn+1)*d3+t*d2+(u-2)*d1]; Rbuf[nn*d3+t*d2+u*d1] = val; }
          else { let val = PC[0] * Rbuf[(nn+1)*d3+(t-1)*d2]; if (t >= 2) val += (t-1)*Rbuf[(nn+1)*d3+(t-2)*d2]; Rbuf[nn*d3+t*d2] = val; }
        }
  const R = (t: number, u: number, v: number) => Rbuf[t * d2 + u * d1 + v]; // n=0

  const prefV = 2 * Math.PI / p;
  let sum = 0;
  for (let t = 0; t <= la[0] + lb[0]; t++)
    for (let u = 0; u <= la[1] + lb[1]; u++)
      for (let v = 0; v <= la[2] + lb[2]; v++)
        sum += Ex[t] * Ey[u] * Ez[v] * R(t, u, v);

  const Kab = Math.exp(-alpha * beta / p * ((A[0]-B[0])**2 + (A[1]-B[1])**2 + (A[2]-B[2])**2));
  return prefV * Kab * sum;
}

/** Nuclear attraction Hessian: d²(D·V)/dR²
 *  Uses angular momentum shift formula on full NAI. */
function nuclearAttractionHessian(
  shells: PrimitiveShell[], atoms: Atom[], norms: number[], n: number,
  D: Matrix, nAtom: number,
): Float64Array {
  const ndim = 3 * nAtom;
  const H = new Float64Array(ndim * ndim);
  const groups = groupShells(shells);

  for (let igA = 0; igA < groups.length; igA++) {
    const gA = groups[igA];
    for (let igB = igA; igB < groups.length; igB++) {
      const gB = groups[igB];
      const angA = ANGULAR_MOMENTUMS[gA.st], angB = ANGULAR_MOMENTUMS[gB.st];

      for (let ia = 0; ia < angA.length; ia++) {
        const la = angA[ia];
        const mu = gA.bi + ia;
        const ibStart = (igA === igB) ? ia : 0;

        for (let ib = ibStart; ib < angB.length; ib++) {
          const lb = angB[ib];
          const nu = gB.bi + ib;
          const sym = (mu === nu) ? 1.0 : 2.0;

          for (const pA of gA.prims) {
            const nA = primitiveNorm(pA.exponent, la[0], la[1], la[2]);
            for (const pB of gB.prims) {
              const nB = primitiveNorm(pB.exponent, lb[0], lb[1], lb[2]);
              const alpha = pA.exponent, beta = pB.exponent;
              const A = [pA.coordinate.x, pA.coordinate.y, pA.coordinate.z];
              const B = [pB.coordinate.x, pB.coordinate.y, pB.coordinate.z];
              const coeff = pA.coefficient * pB.coefficient * nA * nB * norms[mu] * norms[nu];
              const iA = pA.atomIndex, iB = pB.atomIndex;

              for (let iC = 0; iC < nAtom; iC++) {
                const Z = atoms[iC].atomicNumber;
                if (Z === 0) continue;
                const C = [atoms[iC].coordinate.x, atoms[iC].coordinate.y, atoms[iC].coordinate.z];

                const w = -Z * coeff * sym * D.get(mu, nu);
                if (Math.abs(w) < 1e-18) continue;

                // Use shift formula: d/dA_d V = 2α V(la[d]+1,...) - la[d] V(la[d]-1,...)
                // Helper to compute V with shifted angular momentum
                const V = (dlaA: number[], dlaB: number[]) => {
                  const sla = [la[0]+dlaA[0], la[1]+dlaA[1], la[2]+dlaA[2]];
                  const slb = [lb[0]+dlaB[0], lb[1]+dlaB[1], lb[2]+dlaB[2]];
                  if (sla[0] < 0 || sla[1] < 0 || sla[2] < 0) return 0;
                  if (slb[0] < 0 || slb[1] < 0 || slb[2] < 0) return 0;
                  return computeNAI3D(sla, slb, alpha, beta, A, B, C);
                };

                const V0 = V([0,0,0], [0,0,0]);

                // Compute AA, BB, AB blocks [3×3 each]
                const hAA = new Float64Array(9), hBB = new Float64Array(9), hAB = new Float64Array(9);

                for (let d1 = 0; d1 < 3; d1++) {
                  const sh1A = [0,0,0], sh1Am = [0,0,0], sh1B = [0,0,0], sh1Bm = [0,0,0];
                  sh1A[d1] = 1; sh1Am[d1] = -1; sh1B[d1] = 1; sh1Bm[d1] = -1;

                  for (let d2 = 0; d2 < 3; d2++) {
                    const sh2A = [0,0,0], sh2Am = [0,0,0], sh2B = [0,0,0], sh2Bm = [0,0,0];
                    sh2A[d2] = 1; sh2Am[d2] = -1; sh2B[d2] = 1; sh2Bm[d2] = -1;

                    if (d1 === d2) {
                      const d = d1;
                      // d²V/dA_d² = 4α² V(la[d]+2) - 2α(2la[d]+1) V0 + la[d](la[d]-1) V(la[d]-2)
                      const sh2 = [0,0,0]; sh2[d] = 2;
                      const shm2 = [0,0,0]; shm2[d] = -2;
                      hAA[d1*3+d2] = 4*alpha*alpha * V(sh2,[0,0,0])
                                    - 2*alpha*(2*la[d]+1) * V0
                                    + (la[d] >= 2 ? la[d]*(la[d]-1) * V(shm2,[0,0,0]) : 0);

                      hBB[d1*3+d2] = 4*beta*beta * V([0,0,0],sh2)
                                    - 2*beta*(2*lb[d]+1) * V0
                                    + (lb[d] >= 2 ? lb[d]*(lb[d]-1) * V([0,0,0],shm2) : 0);

                      // d²V/dA_d dB_d = [dA_d][dB_d] applied sequentially
                      const dA_V = (shB: number[]) => 2*alpha * V(sh1A, shB) - (la[d] > 0 ? la[d] * V(sh1Am, shB) : 0);
                      hAB[d1*3+d2] = 2*beta * dA_V(sh2B) - (lb[d] > 0 ? lb[d] * dA_V(sh2Bm) : 0);
                    } else {
                      // d²V/dA_d1 dA_d2: product of two shifts on A
                      // = [2α V(la[d1]+1) - la[d1] V(la[d1]-1)] shifted in d2
                      const sh12 = [0,0,0]; sh12[d1] = 1; sh12[d2] = 1;
                      const shm12 = [0,0,0]; shm12[d1] = -1; shm12[d2] = 1;
                      const sh1m2 = [0,0,0]; sh1m2[d1] = 1; sh1m2[d2] = -1;
                      const shm1m2 = [0,0,0]; shm1m2[d1] = -1; shm1m2[d2] = -1;

                      hAA[d1*3+d2] = 4*alpha*alpha * V(sh12,[0,0,0])
                                    - 2*alpha*(la[d1] > 0 ? la[d1] : 0) * V(shm12,[0,0,0])
                                    - 2*alpha*(la[d2] > 0 ? la[d2] : 0) * V(sh1m2,[0,0,0])
                                    + (la[d1] > 0 && la[d2] > 0 ? la[d1]*la[d2] * V(shm1m2,[0,0,0]) : 0);

                      const sh12B = [0,0,0]; sh12B[d1] = 1; sh12B[d2] = 1;
                      const shm12B = [0,0,0]; shm12B[d1] = -1; shm12B[d2] = 1;
                      const sh1m2B = [0,0,0]; sh1m2B[d1] = 1; sh1m2B[d2] = -1;
                      const shm1m2B = [0,0,0]; shm1m2B[d1] = -1; shm1m2B[d2] = -1;

                      hBB[d1*3+d2] = 4*beta*beta * V([0,0,0],sh12B)
                                    - 2*beta*(lb[d1] > 0 ? lb[d1] : 0) * V([0,0,0],shm12B)
                                    - 2*beta*(lb[d2] > 0 ? lb[d2] : 0) * V([0,0,0],sh1m2B)
                                    + (lb[d1] > 0 && lb[d2] > 0 ? lb[d1]*lb[d2] * V([0,0,0],shm1m2B) : 0);

                      // AB cross: dA_d1 * dB_d2
                      const dA_V2 = (shB: number[]) => 2*alpha * V(sh1A, shB) - (la[d1] > 0 ? la[d1] * V(sh1Am, shB) : 0);
                      hAB[d1*3+d2] = 2*beta * dA_V2(sh2B) - (lb[d2] > 0 ? lb[d2] * dA_V2(sh2Bm) : 0);
                    }
                  }
                }

                // Write to global Hessian: AA, BB, AB blocks + TI for C
                for (let d1 = 0; d1 < 3; d1++) {
                  for (let d2 = 0; d2 < 3; d2++) {
                    const vAA = w * hAA[d1*3+d2];
                    const vBB = w * hBB[d1*3+d2];
                    const vAB = w * hAB[d1*3+d2];
                    const vBA = w * hAB[d2*3+d1]; // transpose

                    H[(3*iA+d1)*ndim+(3*iA+d2)] += vAA;
                    H[(3*iB+d1)*ndim+(3*iB+d2)] += vBB;
                    H[(3*iA+d1)*ndim+(3*iB+d2)] += vAB;
                    H[(3*iB+d1)*ndim+(3*iA+d2)] += vBA;

                    // AC = -(AA + AB), CA = transpose
                    const vAC = -(vAA + vAB);
                    H[(3*iA+d1)*ndim+(3*iC+d2)] += vAC;
                    H[(3*iC+d2)*ndim+(3*iA+d1)] += vAC;

                    // BC = -(BA + BB), CB = transpose
                    const vBC = -(vBA + vBB);
                    H[(3*iB+d1)*ndim+(3*iC+d2)] += vBC;
                    H[(3*iC+d2)*ndim+(3*iB+d1)] += vBC;

                    // CC = AA + AB + BA + BB
                    H[(3*iC+d1)*ndim+(3*iC+d2)] += vAA + vAB + vBA + vBB;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return H;
}

// ── Single primitive ERI via McMurchie-Davidson ──────────────────────

/** Compute a single primitive ERI (la lb | lc ld) with given angular momenta.
 *  Returns the bare integral without contraction coefficients/norms. */
function computePrimERI(
  la: number[], lb: number[], lc: number[], ld: number[],
  alpha: number, beta: number, gamma: number, delta: number,
  A: number[], B: number[], C: number[], DD: number[],
): number {
  const p = alpha + beta, q = gamma + delta, eta = p * q / (p + q);
  const P = [(alpha*A[0]+beta*B[0])/p, (alpha*A[1]+beta*B[1])/p, (alpha*A[2]+beta*B[2])/p];
  const Q = [(gamma*C[0]+delta*DD[0])/q, (gamma*C[1]+delta*DD[1])/q, (gamma*C[2]+delta*DD[2])/q];
  const PA = [P[0]-A[0], P[1]-A[1], P[2]-A[2]], PB = [P[0]-B[0], P[1]-B[1], P[2]-B[2]];
  const QC = [Q[0]-C[0], Q[1]-C[1], Q[2]-C[2]], QD = [Q[0]-DD[0], Q[1]-DD[1], Q[2]-DD[2]];
  const RPQ = [P[0]-Q[0], P[1]-Q[1], P[2]-Q[2]];

  const Kab = Math.exp(-alpha*beta/p*((A[0]-B[0])**2+(A[1]-B[1])**2+(A[2]-B[2])**2));
  const Kcd = Math.exp(-gamma*delta/q*((C[0]-DD[0])**2+(C[1]-DD[1])**2+(C[2]-DD[2])**2));

  const ExAB = eCoefficients(la[0], lb[0], p, PA[0], PB[0]);
  const EyAB = eCoefficients(la[1], lb[1], p, PA[1], PB[1]);
  const EzAB = eCoefficients(la[2], lb[2], p, PA[2], PB[2]);
  const ExCD = eCoefficients(lc[0], ld[0], q, QC[0], QD[0]);
  const EyCD = eCoefficients(lc[1], ld[1], q, QC[1], QD[1]);
  const EzCD = eCoefficients(lc[2], ld[2], q, QC[2], QD[2]);

  const maxN = la[0]+lb[0]+lc[0]+ld[0]+la[1]+lb[1]+lc[1]+ld[1]+la[2]+lb[2]+lc[2]+ld[2];
  const RPQ2 = RPQ[0]*RPQ[0]+RPQ[1]*RPQ[1]+RPQ[2]*RPQ[2];
  const boys = boysAll(maxN, eta * RPQ2);

  // R function
  const dim = maxN + 1, d1 = dim + 1, size = d1 ** 4;
  const Rbuf = new Float64Array(size);
  const d2 = d1 * d1, d3 = d2 * d1;
  for (let nn = 0; nn <= maxN; nn++) Rbuf[nn * d3] = (-2 * eta) ** nn * boys[nn];
  for (let nn = maxN - 1; nn >= 0; nn--)
    for (let t = 0; t <= maxN - nn; t++)
      for (let u = 0; u <= maxN - nn - t; u++)
        for (let v = 0; v <= maxN - nn - t - u; v++) {
          if (t + u + v === 0) continue;
          if (v > 0) { let val = RPQ[2]*Rbuf[(nn+1)*d3+t*d2+u*d1+v-1]; if (v>=2) val+=(v-1)*Rbuf[(nn+1)*d3+t*d2+u*d1+v-2]; Rbuf[nn*d3+t*d2+u*d1+v]=val; }
          else if (u > 0) { let val = RPQ[1]*Rbuf[(nn+1)*d3+t*d2+(u-1)*d1]; if (u>=2) val+=(u-1)*Rbuf[(nn+1)*d3+t*d2+(u-2)*d1]; Rbuf[nn*d3+t*d2+u*d1]=val; }
          else { let val = RPQ[0]*Rbuf[(nn+1)*d3+(t-1)*d2]; if (t>=2) val+=(t-1)*Rbuf[(nn+1)*d3+(t-2)*d2]; Rbuf[nn*d3+t*d2]=val; }
        }

  const prefERI = 2 * Math.PI ** 2.5 / (p * q * Math.sqrt(p + q));
  let sum = 0;
  for (let t1 = 0; t1 <= la[0]+lb[0]; t1++)
    for (let u1 = 0; u1 <= la[1]+lb[1]; u1++)
      for (let v1 = 0; v1 <= la[2]+lb[2]; v1++) {
        const e1 = ExAB[t1] * EyAB[u1] * EzAB[v1];
        for (let t2 = 0; t2 <= lc[0]+ld[0]; t2++)
          for (let u2 = 0; u2 <= lc[1]+ld[1]; u2++)
            for (let v2 = 0; v2 <= lc[2]+ld[2]; v2++) {
              const e2 = ExCD[t2] * EyCD[u2] * EzCD[v2];
              const sign = ((t2+u2+v2) & 1) ? -1 : 1;
              sum += e1 * e2 * sign * Rbuf[(t1+t2)*d2+(u1+u2)*d1+(v1+v2)];
            }
      }

  return prefERI * Kab * Kcd * sum;
}

// ── 2-electron Hessian via angular momentum shift ────────────────────

/** Compute 2-electron skeleton Hessian: Σ Γ d²(μν|λσ)/dR²
 *  Uses angular momentum shift formula on full ERIs. */
function twoElectronHessian(
  shells: PrimitiveShell[], norms: number[], n: number,
  D: Matrix, nAtom: number,
): Float64Array {
  const ndim = 3 * nAtom;
  const H = new Float64Array(ndim * ndim);
  const groups = groupShells(shells);
  const nGrp = groups.length;

  // Loop over unique shell quartets (bra-ket symmetry)
  // Full loop — no symmetry optimization for correctness
  for (let igA = 0; igA < nGrp; igA++) {
    const gA = groups[igA];
    for (let igB = 0; igB < nGrp; igB++) {
      const gB = groups[igB];
      for (let igC = 0; igC < nGrp; igC++) {
        const gC = groups[igC];
        for (let igD = 0; igD < nGrp; igD++) {
          const gD = groups[igD];

          const angA = ANGULAR_MOMENTUMS[gA.st], angB = ANGULAR_MOMENTUMS[gB.st];
          const angC = ANGULAR_MOMENTUMS[gC.st], angD = ANGULAR_MOMENTUMS[gD.st];

          for (let ia = 0; ia < angA.length; ia++) {
            const laA = angA[ia]; const mu = gA.bi + ia;
            for (let ib = 0; ib < angB.length; ib++) {
              const laB = angB[ib]; const nu = gB.bi + ib;
              for (let ic = 0; ic < angC.length; ic++) {
                const laC = angC[ic]; const lam = gC.bi + ic;
                for (let id = 0; id < angD.length; id++) {
                  const laD = angD[id]; const sig = gD.bi + id;

                  // Density weighting (GANSU form for full loop, no shell symmetry)
                  const D_ab = D.get(mu, nu), D_cd = D.get(lam, sig);
                  const D_ac = D.get(mu, lam), D_bd = D.get(nu, sig);
                  const D_ad = D.get(mu, sig), D_bc = D.get(nu, lam);
                  const density_w = 0.5 * D_ab * D_cd - 0.125 * (D_ac * D_bd + D_ad * D_bc);
                  if (Math.abs(density_w) < 1e-15) continue;

                  for (const pA of gA.prims) {
                    for (const pB of gB.prims) {
                      for (const pC of gC.prims) {
                        for (const pD of gD.prims) {
                          const alpha = pA.exponent, beta = pB.exponent;
                          const gam = pC.exponent, del = pD.exponent;
                          const nA = primitiveNorm(alpha, laA[0], laA[1], laA[2]);
                          const nB = primitiveNorm(beta, laB[0], laB[1], laB[2]);
                          const nC = primitiveNorm(gam, laC[0], laC[1], laC[2]);
                          const nD = primitiveNorm(del, laD[0], laD[1], laD[2]);
                          const coeff = pA.coefficient * pB.coefficient * pC.coefficient * pD.coefficient
                            * nA * nB * nC * nD * norms[mu] * norms[nu] * norms[lam] * norms[sig];
                          const w = density_w * coeff;
                          if (Math.abs(w) < 1e-18) continue;

                          const AA = [pA.coordinate.x, pA.coordinate.y, pA.coordinate.z];
                          const BB = [pB.coordinate.x, pB.coordinate.y, pB.coordinate.z];
                          const CC = [pC.coordinate.x, pC.coordinate.y, pC.coordinate.z];
                          const DD = [pD.coordinate.x, pD.coordinate.y, pD.coordinate.z];
                          const atomA = pA.atomIndex, atomB = pB.atomIndex;
                          const atomC = pC.atomIndex, atomD = pD.atomIndex;

                          // ── Fast 2e Hessian via 1D integral tables ──
                          // Compute E-coefficients and R function ONCE (extended by +2 for shifts)
                          const p = alpha + beta, q = gam + del, eta = p * q / (p + q);
                          const Px = (alpha*AA[0]+beta*BB[0])/p, Py = (alpha*AA[1]+beta*BB[1])/p, Pz = (alpha*AA[2]+beta*BB[2])/p;
                          const Qx = (gam*CC[0]+del*DD[0])/q, Qy = (gam*CC[1]+del*DD[1])/q, Qz = (gam*CC[2]+del*DD[2])/q;
                          const RPQ = [Px-Qx, Py-Qy, Pz-Qz];
                          const Kab = Math.exp(-alpha*beta/p*((AA[0]-BB[0])**2+(AA[1]-BB[1])**2+(AA[2]-BB[2])**2));
                          const Kcd = Math.exp(-gam*del/q*((CC[0]-DD[0])**2+(CC[1]-DD[1])**2+(CC[2]-DD[2])**2));
                          const pre2e = 2 * Math.PI ** 2.5 / (p * q * Math.sqrt(p + q)) * Kab * Kcd;

                          const la = [laA, laB, laC, laD]; // [4][3]
                          const exps = [alpha, beta, gam, del];
                          const atoms4 = [atomA, atomB, atomC, atomD];

                          // Compute 1D base integral I_d(l1,l2,l3,l4) for each direction d
                          // Using E-coefficients + R function contraction
                          // Extend by +2 for second derivatives
                          const PA = [Px-AA[0],Py-AA[1],Pz-AA[2]], PB = [Px-BB[0],Py-BB[1],Pz-BB[2]];
                          const QC = [Qx-CC[0],Qy-CC[1],Qz-CC[2]], QD = [Qx-DD[0],Qy-DD[1],Qy-DD[2]];
                          // Fix: QD z-component
                          QD[2] = Qz - DD[2];

                          // For shift formula, we need I with each center's AM shifted by -2..+2
                          // Max AM per direction: la[c][d] + 2 for any center c
                          // Total: l1+l2+l3+l4 up to original + 4 (two shifts of +2)
                          const maxN_orig = la[0][0]+la[1][0]+la[2][0]+la[3][0]+la[0][1]+la[1][1]+la[2][1]+la[3][1]+la[0][2]+la[1][2]+la[2][2]+la[3][2];
                          const maxN = maxN_orig + 2; // room for shift: max +2 on one center at a time

                          const RPQ2 = RPQ[0]*RPQ[0]+RPQ[1]*RPQ[1]+RPQ[2]*RPQ[2];
                          const boysV = boysAll(maxN, eta * RPQ2);

                          // R function (extended)
                          const dim2 = maxN + 1, dd1 = dim2 + 1, sz = dd1 ** 4;
                          const Rbuf = new Float64Array(sz);
                          const dd2 = dd1 * dd1, dd3 = dd2 * dd1;
                          for (let nn = 0; nn <= maxN; nn++) Rbuf[nn * dd3] = (-2 * eta) ** nn * boysV[nn];
                          for (let nn = maxN - 1; nn >= 0; nn--)
                            for (let t = 0; t <= maxN - nn; t++)
                              for (let u = 0; u <= maxN - nn - t; u++)
                                for (let v = 0; v <= maxN - nn - t - u; v++) {
                                  if (t + u + v === 0) continue;
                                  if (v > 0) { let vv = RPQ[2]*Rbuf[(nn+1)*dd3+t*dd2+u*dd1+v-1]; if (v>=2) vv+=(v-1)*Rbuf[(nn+1)*dd3+t*dd2+u*dd1+v-2]; Rbuf[nn*dd3+t*dd2+u*dd1+v]=vv; }
                                  else if (u > 0) { let vv = RPQ[1]*Rbuf[(nn+1)*dd3+t*dd2+(u-1)*dd1]; if (u>=2) vv+=(u-1)*Rbuf[(nn+1)*dd3+t*dd2+(u-2)*dd1]; Rbuf[nn*dd3+t*dd2+u*dd1]=vv; }
                                  else { let vv = RPQ[0]*Rbuf[(nn+1)*dd3+(t-1)*dd2]; if (t>=2) vv+=(t-1)*Rbuf[(nn+1)*dd3+(t-2)*dd2]; Rbuf[nn*dd3+t*dd2]=vv; }
                                }
                          const RR = (t: number, u: number, v: number) => Rbuf[t * dd2 + u * dd1 + v];

                          // ── Pre-compute all E-coefficient tables needed ──
                          // For shifts -2..+2 on each center in each direction
                          // braEcache[dir][l1_shift+2][l2_shift+2] = E^{la[0][d]+s1, la[1][d]+s2}
                          // s1, s2 ∈ {-2,-1,0,+1,+2}, indexed as s+2 = 0..4
                          const braEc: (Float64Array|null)[][][] = [];
                          const ketEc: (Float64Array|null)[][][] = [];
                          for (let d = 0; d < 3; d++) {
                            const b: (Float64Array|null)[][] = [];
                            const k: (Float64Array|null)[][] = [];
                            for (let s1 = -2; s1 <= 2; s1++) {
                              const brow: (Float64Array|null)[] = [];
                              const krow: (Float64Array|null)[] = [];
                              for (let s2 = -2; s2 <= 2; s2++) {
                                const l1b = la[0][d]+s1, l2b = la[1][d]+s2;
                                brow.push(l1b >= 0 && l2b >= 0 ? eCoefficients(l1b, l2b, p, PA[d], PB[d]) : null);
                                const l1k = la[2][d]+s1, l2k = la[3][d]+s2;
                                krow.push(l1k >= 0 && l2k >= 0 ? eCoefficients(l1k, l2k, q, QC[d], QD[d]) : null);
                              }
                              b.push(brow);
                              k.push(krow);
                            }
                            braEc.push(b);
                            ketEc.push(k);
                          }

                          // Fast eri3D using cached E-coefficients
                          function eri3D(dl: number[][]): number {
                            // dl[4][3]: AM shifts for each center × direction
                            // Map to cache indices
                            const sA = dl[0], sB = dl[1], sC = dl[2], sD = dl[3];
                            // Check bounds
                            for (let d = 0; d < 3; d++) {
                              if (la[0][d]+sA[d] < 0 || la[1][d]+sB[d] < 0 || la[2][d]+sC[d] < 0 || la[3][d]+sD[d] < 0) return 0;
                            }
                            const Ex = braEc[0][sA[0]+2][sB[0]+2];
                            const Ey = braEc[1][sA[1]+2][sB[1]+2];
                            const Ez = braEc[2][sA[2]+2][sB[2]+2];
                            const Fx = ketEc[0][sC[0]+2][sD[0]+2];
                            const Fy = ketEc[1][sC[1]+2][sD[1]+2];
                            const Fz = ketEc[2][sC[2]+2][sD[2]+2];
                            if (!Ex||!Ey||!Ez||!Fx||!Fy||!Fz) return 0;
                            const maxBx = la[0][0]+sA[0]+la[1][0]+sB[0];
                            const maxBy = la[0][1]+sA[1]+la[1][1]+sB[1];
                            const maxBz = la[0][2]+sA[2]+la[1][2]+sB[2];
                            const maxKx = la[2][0]+sC[0]+la[3][0]+sD[0];
                            const maxKy = la[2][1]+sC[1]+la[3][1]+sD[1];
                            const maxKz = la[2][2]+sC[2]+la[3][2]+sD[2];
                            // Fast path: all-s case (maxB=maxK=0)
                            if (maxBx+maxBy+maxBz+maxKx+maxKy+maxKz === 0) {
                              return Ex[0]*Ey[0]*Ez[0]*Fx[0]*Fy[0]*Fz[0]*RR(0,0,0) * pre2e;
                            }
                            let sum = 0;
                            for (let t1=0;t1<=maxBx;t1++)
                              for (let u1=0;u1<=maxBy;u1++)
                                for (let v1=0;v1<=maxBz;v1++) {
                                  const e1 = Ex[t1]*Ey[u1]*Ez[v1];
                                  if (Math.abs(e1) < 1e-18) continue;
                                  for (let t2=0;t2<=maxKx;t2++)
                                    for (let u2=0;u2<=maxKy;u2++)
                                      for (let v2=0;v2<=maxKz;v2++) {
                                        const sign=((t2+u2+v2)&1)?-1:1;
                                        sum += e1*Fx[t2]*Fy[u2]*Fz[v2]*sign*RR(t1+t2,u1+u2,v1+v2);
                                      }
                                }
                            return sum * pre2e;
                          }

                          const Z4 = [[0,0,0],[0,0,0],[0,0,0],[0,0,0]];
                          const eri0 = eri3D(Z4);

                          // Reusable shift arrays (avoid allocation in hot loop)
                          const sh = [[0,0,0],[0,0,0],[0,0,0],[0,0,0]];
                          function resetSh() { for (let i=0;i<4;i++) sh[i][0]=sh[i][1]=sh[i][2]=0; }
                          function eriSh(c1: number, d1: number, s1: number, c2?: number, d2?: number, s2?: number): number {
                            resetSh();
                            sh[c1][d1] = s1;
                            if (c2 !== undefined) sh[c2!][d2!] = s2!;
                            return eri3D(sh);
                          }

                          for (let c1 = 0; c1 < 4; c1++) {
                            for (let c2 = c1; c2 < 4; c2++) {
                              for (let d1 = 0; d1 < 3; d1++) {
                                for (let d2 = (c1 === c2 ? d1 : 0); d2 < 3; d2++) {
                                  let val: number;

                                  if (c1 === c2) {
                                    const e = exps[c1], l = la[c1][d1];
                                    if (d1 === d2) {
                                      val = 4*e*e*eriSh(c1,d1,2) - 2*e*(2*l+1)*eri0 + (l>=2?l*(l-1)*eriSh(c1,d1,-2):0);
                                    } else {
                                      const l2 = la[c1][d2];
                                      // Same center, different dirs: need combined shifts
                                      resetSh(); sh[c1][d1]=1; sh[c1][d2]=1; const v12=eri3D(sh);
                                      resetSh(); sh[c1][d1]=-1; sh[c1][d2]=1; const vm12=eri3D(sh);
                                      resetSh(); sh[c1][d1]=1; sh[c1][d2]=-1; const v1m2=eri3D(sh);
                                      resetSh(); sh[c1][d1]=-1; sh[c1][d2]=-1; const vm1m2=eri3D(sh);
                                      val = 4*e*e*v12 - (l>0?2*e*l*vm12:0) - (l2>0?2*e*l2*v1m2:0) + (l>0&&l2>0?l*l2*vm1m2:0);
                                    }
                                  } else {
                                    const e1=exps[c1], e2=exps[c2], l1v=la[c1][d1], l2v=la[c2][d2];
                                    const vpp=eriSh(c1,d1,1,c2,d2,1);
                                    const vmp=l1v>0?eriSh(c1,d1,-1,c2,d2,1):0;
                                    const vpm=l2v>0?eriSh(c1,d1,1,c2,d2,-1):0;
                                    const vmm=l1v>0&&l2v>0?eriSh(c1,d1,-1,c2,d2,-1):0;
                                    val = 4*e1*e2*vpp - 2*e2*l1v*vmp - 2*e1*l2v*vpm + l1v*l2v*vmm;
                                  }

                                  if (Math.abs(val) < 1e-18) continue;
                                  val *= w;

                                  const g1 = 3 * atoms4[c1] + d1;
                                  const g2 = 3 * atoms4[c2] + d2;
                                  const l1 = c1 * 3 + d1, l2 = c2 * 3 + d2;

                                  if (g1 !== g2) {
                                    H[g1 * ndim + g2] += val;
                                    H[g2 * ndim + g1] += val;
                                  } else if (l1 === l2) {
                                    H[g1 * ndim + g2] += val;
                                  } else {
                                    // Different local centers map to same global index
                                    H[g1 * ndim + g2] += 2 * val;
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
  return H;
}

// ── Skeleton dG matrix (2e Fock derivative) via shift formula ─────────

/** Compute dG_μν/dR for a single perturbation direction and ADD to dF matrix.
 *  dG_μν = Σ_λσ D_λσ [d(μν|λσ)/dR - 0.5 d(μλ|νσ)/dR]
 *  Uses angular momentum shift formula on primitive ERIs. */
function computeSkeletonDG(
  shells: PrimitiveShell[], norms: number[], n: number,
  D: Matrix, pertAtom: number, pertDir: number, dF: Matrix,
): void {
  const groups = groupShells(shells);

  // For each shell quartet (μν|λσ), compute d(μν|λσ)/dR
  // dR acts on all functions centered on pertAtom
  for (const gA of groups)
    for (const gB of groups)
      for (const gC of groups)
        for (const gD of groups) {
          const angA = ANGULAR_MOMENTUMS[gA.st], angB = ANGULAR_MOMENTUMS[gB.st];
          const angC = ANGULAR_MOMENTUMS[gC.st], angD = ANGULAR_MOMENTUMS[gD.st];

          for (let ia = 0; ia < angA.length; ia++) {
            const laA = angA[ia]; const mu = gA.bi + ia;
            for (let ib = 0; ib < angB.length; ib++) {
              const laB = angB[ib]; const nu = gB.bi + ib;
              for (let ic = 0; ic < angC.length; ic++) {
                const laC = angC[ic]; const lam = gC.bi + ic;
                for (let id = 0; id < angD.length; id++) {
                  const laD = angD[id]; const sig = gD.bi + id;

                  // Weight: Coulomb D_λσ for (μν|λσ), exchange D_λσ for (μλ|νσ)
                  const Dls = D.get(lam, sig);
                  const Dns = D.get(nu, sig); // for exchange: D_νσ in (μλ|νσ)
                  if (Math.abs(Dls) < 1e-15 && Math.abs(Dns) < 1e-15) continue;

                  for (const pA of gA.prims) {
                    for (const pB of gB.prims) {
                      for (const pC of gC.prims) {
                        for (const pD of gD.prims) {
                          const alpha = pA.exponent, beta = pB.exponent;
                          const gam = pC.exponent, del = pD.exponent;
                          const nA = primitiveNorm(alpha, laA[0], laA[1], laA[2]);
                          const nB = primitiveNorm(beta, laB[0], laB[1], laB[2]);
                          const nC = primitiveNorm(gam, laC[0], laC[1], laC[2]);
                          const nD = primitiveNorm(del, laD[0], laD[1], laD[2]);
                          const coeff = pA.coefficient * pB.coefficient * pC.coefficient * pD.coefficient
                            * nA * nB * nC * nD * norms[mu] * norms[nu] * norms[lam] * norms[sig];
                          if (Math.abs(coeff) < 1e-18) continue;

                          const AA = [pA.coordinate.x, pA.coordinate.y, pA.coordinate.z];
                          const BB = [pB.coordinate.x, pB.coordinate.y, pB.coordinate.z];
                          const CC = [pC.coordinate.x, pC.coordinate.y, pC.coordinate.z];
                          const DD2 = [pD.coordinate.x, pD.coordinate.y, pD.coordinate.z];

                          // Compute d(μν|λσ)/dR_pertAtom_pertDir via shift formula
                          // Sum shift over all 4 centers that are on pertAtom
                          const d = pertDir;
                          const Z = [0,0,0];
                          let dERI = 0;

                          const centers = [
                            { exp: alpha, la: laA, atom: pA.atomIndex },
                            { exp: beta, la: laB, atom: pB.atomIndex },
                            { exp: gam, la: laC, atom: pC.atomIndex },
                            { exp: del, la: laD, atom: pD.atomIndex },
                          ];

                          for (let ci = 0; ci < 4; ci++) {
                            if (centers[ci].atom !== pertAtom) continue;
                            const s1 = [0,0,0]; s1[d] = 1;
                            const sm1 = [0,0,0]; sm1[d] = -1;
                            const da = [Z,Z,Z,Z].map(x => [...x]);
                            const dm = [Z,Z,Z,Z].map(x => [...x]);
                            da[ci] = s1; dm[ci] = sm1;

                            const eri_p = computePrimERI(
                              [laA[0]+da[0][0],laA[1]+da[0][1],laA[2]+da[0][2]],
                              [laB[0]+da[1][0],laB[1]+da[1][1],laB[2]+da[1][2]],
                              [laC[0]+da[2][0],laC[1]+da[2][1],laC[2]+da[2][2]],
                              [laD[0]+da[3][0],laD[1]+da[3][1],laD[2]+da[3][2]],
                              alpha, beta, gam, del, AA, BB, CC, DD2);
                            const eri_m = (centers[ci].la[d] > 0) ? computePrimERI(
                              [laA[0]+dm[0][0],laA[1]+dm[0][1],laA[2]+dm[0][2]],
                              [laB[0]+dm[1][0],laB[1]+dm[1][1],laB[2]+dm[1][2]],
                              [laC[0]+dm[2][0],laC[1]+dm[2][1],laC[2]+dm[2][2]],
                              [laD[0]+dm[3][0],laD[1]+dm[3][1],laD[2]+dm[3][2]],
                              alpha, beta, gam, del, AA, BB, CC, DD2) : 0;

                            dERI += 2 * centers[ci].exp * eri_p - centers[ci].la[d] * eri_m;
                          }

                          // dG_μν += D_λσ * d(μν|λσ)/dR * coeff
                          dF.set(mu, nu, dF.get(mu, nu) + Dls * dERI * coeff);
                          // Exchange: -0.5 D_λσ * d(μλ|νσ)/dR
                          // This requires (μλ|νσ) derivative — different index ordering
                          // Skip exchange for now, add later for full accuracy
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

// ── Contracted ERI from shell data ───────────────────────────────────

/** Compute a contracted ERI (μν|λσ) from primitive shells. */
function computeContractedERI(
  mu: number, nu: number, lam: number, sig: number,
  shells: PrimitiveShell[], norms: number[], n: number,
): number {
  // Find primitives for each basis function
  const groups = groupShells(shells);
  // Map basis index to group+component
  function findGroup(bi: number): { grp: ShellGrp; comp: number } | null {
    for (const g of groups) {
      const ang = ANGULAR_MOMENTUMS[g.st];
      for (let c = 0; c < ang.length; c++) {
        if (g.bi + c === bi) return { grp: g, comp: c };
      }
    }
    return null;
  }

  const gA = findGroup(mu), gB = findGroup(nu), gC = findGroup(lam), gD = findGroup(sig);
  if (!gA || !gB || !gC || !gD) return 0;

  const angA = ANGULAR_MOMENTUMS[gA.grp.st][gA.comp];
  const angB = ANGULAR_MOMENTUMS[gB.grp.st][gB.comp];
  const angC = ANGULAR_MOMENTUMS[gC.grp.st][gC.comp];
  const angD = ANGULAR_MOMENTUMS[gD.grp.st][gD.comp];

  let sum = 0;
  for (const pA of gA.grp.prims) {
    const nA = primitiveNorm(pA.exponent, angA[0], angA[1], angA[2]);
    for (const pB of gB.grp.prims) {
      const nB = primitiveNorm(pB.exponent, angB[0], angB[1], angB[2]);
      for (const pC of gC.grp.prims) {
        const nC = primitiveNorm(pC.exponent, angC[0], angC[1], angC[2]);
        for (const pD of gD.grp.prims) {
          const nD = primitiveNorm(pD.exponent, angD[0], angD[1], angD[2]);
          const val = computePrimERI(
            angA, angB, angC, angD,
            pA.exponent, pB.exponent, pC.exponent, pD.exponent,
            [pA.coordinate.x, pA.coordinate.y, pA.coordinate.z],
            [pB.coordinate.x, pB.coordinate.y, pB.coordinate.z],
            [pC.coordinate.x, pC.coordinate.y, pC.coordinate.z],
            [pD.coordinate.x, pD.coordinate.y, pD.coordinate.z],
          );
          sum += val * pA.coefficient * pB.coefficient * pC.coefficient * pD.coefficient
               * nA * nB * nC * nD * norms[mu] * norms[nu] * norms[lam] * norms[sig];
        }
      }
    }
  }
  return sum;
}

// ── Assemble analytical Hessian (skeleton only, no CPHF yet) ─────────

export function computeAnalyticalHessianSkeleton(
  shells: PrimitiveShell[], atoms: Atom[], norms: number[], n: number,
  nocc: number, D: Matrix, C: Matrix, eps: Float64Array,
): Float64Array {
  const nAt = atoms.length, ndim = 3 * nAt;
  const W = computeWMatrix(C, eps, nocc, n);

  // Nuclear repulsion
  const Hnn = nuclearRepulsionHessian(atoms);

  // 1-electron (overlap + kinetic) — skeleton part
  const H1e = oneElectronHessianSK(shells, norms, n, D, W, nAt);

  // Nuclear attraction
  const Hna = nuclearAttractionHessian(shells, atoms, norms, n, D, nAt);

  // 2-electron — try WASM first
  let H2e: Float64Array;
  const wasmH2e = isWasmAvailable() ? compute2eHessianWasm(
    packShellsExt(shells),
    new Float64Array(D.data),
    new Float64Array(norms),
    n, nAt,
  ) : null;
  if (wasmH2e) {
    H2e = wasmH2e;
  } else {
    H2e = twoElectronHessian(shells, norms, n, D, nAt);
  }

  // Combine
  const H = new Float64Array(ndim * ndim);
  for (let i = 0; i < ndim * ndim; i++) H[i] = Hnn[i] + H1e[i] + Hna[i] + H2e[i];

  return H;
}

/** Full analytical Hessian including CPHF response terms.
 *  H_full = H_skeleton + Σ_x Σ_y U^x · (response contributions)
 *
 *  The CPHF response adds: Σ_x (dD/dR_x)^T F^y = 4 Σ_{ai} U^x_{ai} B^y_{ai}
 *  where B^y is the skeleton Fock derivative in MO basis. */
export async function computeFullAnalyticalHessian(
  shells: PrimitiveShell[], atoms: Atom[], norms: number[], n: number,
  nocc: number, D: Matrix, C: Matrix, eps: Float64Array,
  eri: ERIStored,
  onProgress?: (msg: string) => void,
): Promise<Float64Array> {
  const nAt = atoms.length, ndim = 3 * nAt;
  const nvir = n - nocc;

  const yield_ = () => new Promise<void>(r => setTimeout(r, 0));

  onProgress?.('Skeleton Hessian...');
  await yield_();
  const Hskel = computeAnalyticalHessianSkeleton(shells, atoms, norms, n, nocc, D, C, eps);

  if (nvir === 0) return Hskel; // no virtual orbitals → no CPHF needed

  // ── CPHF ──
  onProgress?.('MO-ERI transform...');
  await yield_();
  let moERI: Float64Array;
  const wasmMOERI = isWasmAvailable() ? transformERItoMOWasm(
    new Float64Array(C.data),
    new Float64Array(eri.data),
    n,
  ) : null;
  if (wasmMOERI) {
    moERI = wasmMOERI;
  } else {
    moERI = transformERItoMO(C, eri, n);
  }

  // Build CPHF RHS for each perturbation direction
  // B^x_{ai} = F^x_{ai} - ε_i S^x_{ai}  (in MO basis)
  // F^x = dH/dR_x + dG/dR_x (skeleton Fock derivative)
  // We compute these via FD of AO integrals (practical for small basis)
  onProgress?.('CPHF RHS...');
  await yield_();
  const dim = nocc * nvir;
  const rhs = new Float64Array(ndim * dim);
  const h = 5e-6;

  // Pre-compute all dG matrices in one pass (avoids 3N separate quartet loops)
  onProgress?.('Skeleton dG matrices...');
  await yield_();
  const allDG = computeAllSkeletonDGMatrices(shells, norms, n, D, nAt);

  for (let pert = 0; pert < ndim; pert++) {
    const atomIdx = Math.floor(pert / 3), dirIdx = pert % 3;

    // Displaced geometries
    function shiftAtoms(delta: number): Atom[] {
      return atoms.map((a, i) => {
        if (i !== atomIdx) return a;
        const c = { ...a.coordinate };
        if (dirIdx === 0) c.x += delta;
        else if (dirIdx === 1) c.y += delta;
        else c.z += delta;
        return { ...a, coordinate: c };
      });
    }
    function shiftShells(delta: number): PrimitiveShell[] {
      return shells.map(ps => {
        if (ps.atomIndex !== atomIdx) return ps;
        const c = { ...ps.coordinate };
        if (dirIdx === 0) c.x += delta;
        else if (dirIdx === 1) c.y += delta;
        else c.z += delta;
        return { ...ps, coordinate: c };
      });
    }

    const atomsP = shiftAtoms(h), atomsM = shiftAtoms(-h);
    const shellsP = shiftShells(h), shellsM = shiftShells(-h);

    // 1e integrals at ±h
    const intP = computeOneElectronIntegrals(shellsP, atomsP, norms, n);
    const intM = computeOneElectronIntegrals(shellsM, atomsM, norms, n);

    // dS/dR and dF/dR (skeleton Fock = H + G with fixed D)
    // dH/dR via FD
    const dS = new Matrix(n, n);
    const dF = new Matrix(n, n);
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++) {
        dS.set(i, j, (intP.overlap.get(i, j) - intM.overlap.get(i, j)) / (2 * h));
        dF.set(i, j, (intP.coreHamiltonian.get(i, j) - intM.coreHamiltonian.get(i, j)) / (2 * h));
      }

    // Add dG/dR from pre-computed matrices
    const dGmat = allDG[pert];
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++)
        dF.set(i, j, dF.get(i, j) + dGmat.get(i, j));

    // Transform to MO basis: dF_pq = Σ_μν C_μp dF_μν C_νq
    // dS_pq = Σ_μν C_μp dS_μν C_νq
    const dF_mo = new Matrix(n, n);
    const dS_mo = new Matrix(n, n);
    for (let p = 0; p < n; p++)
      for (let q = 0; q < n; q++) {
        let fv = 0, sv = 0;
        for (let mu = 0; mu < n; mu++)
          for (let nu = 0; nu < n; nu++) {
            const cc = C.get(mu, p) * C.get(nu, q);
            fv += cc * dF.get(mu, nu);
            sv += cc * dS.get(mu, nu);
          }
        dF_mo.set(p, q, fv);
        dS_mo.set(p, q, sv);
      }

    // B^x_{ia} = dF_{ai} - ε_i dS_{ai}  (note: a=virtual, i=occupied)
    for (let i = 0; i < nocc; i++)
      for (let a = 0; a < nvir; a++) {
        const ai = i * nvir + a;
        const aMO = nocc + a;
        rhs[pert * dim + ai] = -(dF_mo.get(aMO, i) - eps[i] * dS_mo.get(aMO, i));
      }
  }

  // Step 2b: occ-occ density response correction to CPHF RHS
  // D_oo = -2 C s1oo C^T → G(D_oo) → project to MO vir-occ → subtract from RHS
  {
    // Compute s1oo for each perturbation (FD of S in MO occ-occ block)
    for (let pert = 0; pert < ndim; pert++) {
      const atomIdx = Math.floor(pert / 3), dirIdx = pert % 3;
      function shiftShellsRHS(delta: number): PrimitiveShell[] {
        return shells.map(ps => {
          if (ps.atomIndex !== atomIdx) return ps;
          const c = { ...ps.coordinate };
          if (dirIdx === 0) c.x += delta; else if (dirIdx === 1) c.y += delta; else c.z += delta;
          return { ...ps, coordinate: c };
        });
      }
      function shiftAtomsRHS(delta: number): Atom[] {
        return atoms.map((a, i) => {
          if (i !== atomIdx) return a;
          const c = { ...a.coordinate };
          if (dirIdx === 0) c.x += delta; else if (dirIdx === 1) c.y += delta; else c.z += delta;
          return { ...a, coordinate: c };
        });
      }
      const sp = shiftShellsRHS(h), sm = shiftShellsRHS(-h);
      const ap = shiftAtomsRHS(h), am = shiftAtomsRHS(-h);
      const ip = computeOneElectronIntegrals(sp, ap, norms, n);
      const im = computeOneElectronIntegrals(sm, am, norms, n);

      // s1oo[i,j] = C^T dS C (occ-occ block)
      const s1oo_pert = new Float64Array(nocc * nocc);
      for (let i = 0; i < nocc; i++)
        for (let j = 0; j < nocc; j++) {
          let v = 0;
          for (let mu = 0; mu < n; mu++)
            for (let nu = 0; nu < n; nu++)
              v += C.get(mu, i) * ((ip.overlap.get(mu, nu) - im.overlap.get(mu, nu)) / (2 * h)) * C.get(nu, j);
          s1oo_pert[i * nocc + j] = v;
        }

      // D_oo = -2 C_occ s1oo C_occ^T
      const D_oo = new Matrix(n, n);
      for (let mu = 0; mu < n; mu++)
        for (let nu = 0; nu < n; nu++) {
          let v = 0;
          for (let i = 0; i < nocc; i++)
            for (let j = 0; j < nocc; j++)
              v += C.get(mu, i) * s1oo_pert[i * nocc + j] * C.get(nu, j);
          D_oo.set(mu, nu, -2 * v);
        }

      // G(D_oo) = Σ_λσ D_oo_λσ [(μν|λσ) - 0.5(μλ|νσ)]
      const G_oo = new Matrix(n, n);
      for (let mu = 0; mu < n; mu++)
        for (let nu = 0; nu < n; nu++) {
          let v = 0;
          for (let lam = 0; lam < n; lam++)
            for (let sig = 0; sig < n; sig++) {
              const d = D_oo.get(lam, sig);
              if (Math.abs(d) < 1e-15) continue;
              v += d * (eri.get(mu, nu, lam, sig) - 0.5 * eri.get(mu, lam, nu, sig));
            }
          G_oo.set(mu, nu, v);
        }

      // G_oo_MO = C^T G_oo C, then subtract vir-occ block from RHS
      for (let i = 0; i < nocc; i++)
        for (let a = 0; a < nvir; a++) {
          const aMO = nocc + a;
          let v = 0;
          for (let mu = 0; mu < n; mu++)
            for (let nu = 0; nu < n; nu++)
              v += C.get(mu, aMO) * G_oo.get(mu, nu) * C.get(nu, i);
          rhs[pert * dim + i * nvir + a] -= v;
        }
    }
  }

  // Solve CPHF — try WASM first
  onProgress?.('CPHF solving...');
  await yield_();
  let U: Float64Array;
  const wasmU = isWasmAvailable() ? solveCPHFWasm(
    moERI, new Float64Array(eps), rhs, nocc, nvir, n, ndim, 1e-8, 200,
  ) : null;
  if (wasmU) {
    onProgress?.(`CPHF: ${ndim} perturbations solved (WASM)`);
    U = wasmU;
  } else {
    U = solveCPHF(moERI, eps, rhs, nocc, nvir, n, ndim, 1e-8, 200, onProgress);
  }

  // ── Build response quantities from CPHF solution ──
  onProgress?.('Response densities...');
  await yield_();

  // For each perturbation, build:
  //   dm1[pert]  = response density in AO basis
  //   dm1e[pert] = energy-weighted response density
  //   mo_e1[pert] = orbital energy response (occ×occ)
  // h1ao[pert] = dF_AO (skeleton Fock derivative, already computed as dF above)
  // s1ao[pert] = dS_AO (already computed)

  // We stored dS and dF matrices per perturbation in AO basis.
  // Need to recompute them or store. Recompute:
  const h1ao: Matrix[] = [], s1ao: Matrix[] = [];
  const s1oo: Float64Array[] = []; // S^(1) in occ-occ MO block
  const dm1_all: Matrix[] = [];
  const dm1e_all: Matrix[] = [];
  const mo_e1_all: Float64Array[] = [];

  for (let pert = 0; pert < ndim; pert++) {
    const atomIdx = Math.floor(pert / 3), dirIdx = pert % 3;
    function shiftShells2(delta: number): PrimitiveShell[] {
      return shells.map(ps => {
        if (ps.atomIndex !== atomIdx) return ps;
        const c = { ...ps.coordinate };
        if (dirIdx === 0) c.x += delta; else if (dirIdx === 1) c.y += delta; else c.z += delta;
        return { ...ps, coordinate: c };
      });
    }
    function shiftAtoms2(delta: number): Atom[] {
      return atoms.map((a, i) => {
        if (i !== atomIdx) return a;
        const c = { ...a.coordinate };
        if (dirIdx === 0) c.x += delta; else if (dirIdx === 1) c.y += delta; else c.z += delta;
        return { ...a, coordinate: c };
      });
    }
    const sp = shiftShells2(h), sm = shiftShells2(-h);
    const ap = shiftAtoms2(h), am = shiftAtoms2(-h);
    const ip = computeOneElectronIntegrals(sp, ap, norms, n);
    const im = computeOneElectronIntegrals(sm, am, norms, n);

    const dS_ao = new Matrix(n, n), dH_ao = new Matrix(n, n);
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++) {
        dS_ao.set(i, j, (ip.overlap.get(i, j) - im.overlap.get(i, j)) / (2 * h));
        dH_ao.set(i, j, (ip.coreHamiltonian.get(i, j) - im.coreHamiltonian.get(i, j)) / (2 * h));
      }
    // Add dG from pre-computed matrices
    const dGm = allDG[pert];
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++)
        dH_ao.set(i, j, dH_ao.get(i, j) + dGm.get(i, j));
    h1ao.push(dH_ao); // h1ao = dH + dG = skeleton Fock derivative
    s1ao.push(dS_ao);

    // S^(1) in occ-occ MO block: s1oo[i,j] = C^T dS C restricted to occ×occ
    const s1oo_pert = new Float64Array(nocc * nocc);
    for (let i = 0; i < nocc; i++)
      for (let j = 0; j < nocc; j++) {
        let v = 0;
        for (let mu = 0; mu < n; mu++)
          for (let nu = 0; nu < n; nu++)
            v += C.get(mu, i) * dS_ao.get(mu, nu) * C.get(nu, j);
        s1oo_pert[i * nocc + j] = v;
      }
    s1oo.push(s1oo_pert);

    // Build mo1[p,i] = U[a,i] (vir-occ block) + (-0.5 S^(1)[j,i]) (occ-occ block)
    const mo1 = new Float64Array(n * nocc); // mo1[p*nocc+i]
    for (let i = 0; i < nocc; i++)
      for (let a = 0; a < nvir; a++)
        mo1[(nocc + a) * nocc + i] = U[pert * dim + i * nvir + a];
    for (let j = 0; j < nocc; j++)
      for (let i = 0; i < nocc; i++)
        mo1[j * nocc + i] = -0.5 * s1oo_pert[j * nocc + i];

    // temp[μ,i] = Σ_p C[μ,p] mo1[p,i]
    const temp = new Float64Array(n * nocc);
    for (let mu = 0; mu < n; mu++)
      for (let i = 0; i < nocc; i++) {
        let v = 0;
        for (let p = 0; p < n; p++) v += C.get(mu, p) * mo1[p * nocc + i];
        temp[mu * nocc + i] = v;
      }

    // dm1[μ,ν] = Σ_i temp[μ,i] C[ν,i]
    const dm1 = new Matrix(n, n);
    for (let mu = 0; mu < n; mu++)
      for (let nu = 0; nu < n; nu++) {
        let v = 0;
        for (let i = 0; i < nocc; i++) v += temp[mu * nocc + i] * C.get(nu, i);
        dm1.set(mu, nu, v);
      }
    dm1_all.push(dm1);

    // dm1e[μ,ν] = Σ_i temp[μ,i] ε_i C[ν,i]
    const dm1e = new Matrix(n, n);
    for (let mu = 0; mu < n; mu++)
      for (let nu = 0; nu < n; nu++) {
        let v = 0;
        for (let i = 0; i < nocc; i++) v += temp[mu * nocc + i] * eps[i] * C.get(nu, i);
        dm1e.set(mu, nu, v);
      }
    dm1e_all.push(dm1e);

    // Compute vhf1 = G(D1) where D1 = 2(dm1 + dm1^T)
    // Then F_tot_MO = C^T (h1ao + vhf1) C
    // mo_e1[i,j] = F_tot_MO[i,j] - 0.5*(ε_i+ε_j)*s1oo[i,j]
    // For vhf1, we need G(D1) = Σ_λσ D1_λσ [(μν|λσ) - 0.5(μλ|νσ)]
    // Use stored AO ERIs
    const D1 = new Matrix(n, n);
    for (let mu = 0; mu < n; mu++)
      for (let nu = 0; nu < n; nu++)
        D1.set(mu, nu, 2 * (dm1.get(mu, nu) + dm1.get(nu, mu)));

    const vhf1 = new Matrix(n, n);
    for (let mu = 0; mu < n; mu++)
      for (let nu = 0; nu < n; nu++) {
        let v = 0;
        for (let lam = 0; lam < n; lam++)
          for (let sig = 0; sig < n; sig++) {
            const d1 = D1.get(lam, sig);
            if (Math.abs(d1) < 1e-15) continue;
            v += d1 * (eri.get(mu, nu, lam, sig) - 0.5 * eri.get(mu, lam, nu, sig));
          }
        vhf1.set(mu, nu, v);
      }

    // F_tot = h1ao + vhf1
    const F_tot_mo = new Float64Array(n * n);
    for (let p = 0; p < n; p++)
      for (let q = 0; q < n; q++) {
        let v = 0;
        for (let mu = 0; mu < n; mu++)
          for (let nu = 0; nu < n; nu++)
            v += C.get(mu, p) * (h1ao[pert].get(mu, nu) + vhf1.get(mu, nu)) * C.get(nu, q);
        F_tot_mo[p * n + q] = v;
      }

    const mo_e1 = new Float64Array(nocc * nocc);
    for (let i = 0; i < nocc; i++)
      for (let j = 0; j < nocc; j++)
        mo_e1[i * nocc + j] = F_tot_mo[i * n + j] - 0.5 * (eps[i] + eps[j]) * s1oo_pert[i * nocc + j];
    mo_e1_all.push(mo_e1);
  }

  // ── Assemble response Hessian ──
  onProgress?.('Assembling response Hessian...');
  await yield_();
  for (let x = 0; x < ndim; x++) {
    for (let y = 0; y < ndim; y++) {
      let term1 = 0, term2 = 0, term3 = 0;
      for (let mu = 0; mu < n; mu++)
        for (let nu = 0; nu < n; nu++) {
          term1 += h1ao[x].get(mu, nu) * dm1_all[y].get(mu, nu);
          term2 += s1ao[x].get(mu, nu) * dm1e_all[y].get(mu, nu);
        }
      for (let k = 0; k < nocc * nocc; k++)
        term3 += s1oo[x][k] * mo_e1_all[y][k];

      Hskel[x * ndim + y] += 4 * term1 - 4 * term2 - 2 * term3;
    }
  }

  // Symmetrize
  for (let i = 0; i < ndim; i++)
    for (let j = i + 1; j < ndim; j++) {
      const avg = 0.5 * (Hskel[i * ndim + j] + Hskel[j * ndim + i]);
      Hskel[i * ndim + j] = avg;
      Hskel[j * ndim + i] = avg;
    }

  onProgress?.('Analytical Hessian complete');
  return Hskel;
}
