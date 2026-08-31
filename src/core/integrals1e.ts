/** 1-Electron integrals via McMurchie-Davidson — overlap, kinetic, nuclear attraction */

import type { Atom, PrimitiveShell } from './types';
import { ANGULAR_MOMENTUMS, shellTypeToNumBasis } from './constants';
import { Matrix } from '../linalg/matrix';
import { boysAll } from './boys';

/** Double factorial: n!! = n * (n-2) * (n-4) * ... (1 or 2) */
function doubleFactorial(n: number): number {
  if (n <= 1) return 1;
  let result = 1;
  for (let i = n; i >= 2; i -= 2) result *= i;
  return result;
}

/** Primitive Gaussian normalization factor N(α, lx, ly, lz)
 *  N = 2^L / sqrt((2lx-1)!! * (2ly-1)!! * (2lz-1)!!) * (2/π)^{3/4} * α^{(2L+3)/4} */
export function primitiveNorm(alpha: number, lx: number, ly: number, lz: number): number {
  const L = lx + ly + lz;
  const dblFact = doubleFactorial(2 * lx - 1) * doubleFactorial(2 * ly - 1) * doubleFactorial(2 * lz - 1);
  return Math.pow(2, L) / Math.sqrt(dblFact)
    * Math.pow(2.0 / Math.PI, 0.75)
    * Math.pow(alpha, (2 * L + 3) / 4.0);
}

/** Compute overlap and core Hamiltonian (kinetic + nuclear attraction) matrices */
export function computeOneElectronIntegrals(
  primitiveShells: PrimitiveShell[],
  atoms: Atom[],
  normFactors: number[],
  numBasis: number,
): { overlap: Matrix; coreHamiltonian: Matrix; kinetic: Matrix } {
  const S = new Matrix(numBasis, numBasis);
  const T = new Matrix(numBasis, numBasis);
  const V = new Matrix(numBasis, numBasis);

  // Group primitive shells by (basisIndex, shellType)
  // Each contracted shell = all primitives with same basisIndex
  const shellGroups = groupPrimitiveShells(primitiveShells);

  for (let iGrp = 0; iGrp < shellGroups.length; iGrp++) {
    for (let jGrp = iGrp; jGrp < shellGroups.length; jGrp++) {
      const grpA = shellGroups[iGrp];
      const grpB = shellGroups[jGrp];
      computeShellPairIntegrals(grpA, grpB, atoms, normFactors, S, T, V, iGrp === jGrp);
    }
  }

  // Core Hamiltonian = T + V
  const H = new Matrix(numBasis, numBasis);
  for (let k = 0; k < numBasis * numBasis; k++) {
    H.data[k] = T.data[k] + V.data[k];
  }

  return { overlap: S, coreHamiltonian: H, kinetic: T };
}

interface ShellGroup {
  primitives: PrimitiveShell[];
  basisIndex: number;
  shellType: number;
}

function groupPrimitiveShells(shells: PrimitiveShell[]): ShellGroup[] {
  const groups: ShellGroup[] = [];
  let current: ShellGroup | null = null;

  for (const ps of shells) {
    if (!current || current.basisIndex !== ps.basisIndex || current.shellType !== ps.shellType) {
      current = { primitives: [], basisIndex: ps.basisIndex, shellType: ps.shellType };
      groups.push(current);
    }
    current.primitives.push(ps);
  }

  return groups;
}

function computeShellPairIntegrals(
  grpA: ShellGroup,
  grpB: ShellGroup,
  atoms: Atom[],
  normFactors: number[],
  S: Matrix,
  T: Matrix,
  V: Matrix,
  sameGroup: boolean,
) {
  const la = grpA.shellType;
  const lb = grpB.shellType;
  const angA = ANGULAR_MOMENTUMS[la];
  const angB = ANGULAR_MOMENTUMS[lb];
  const numA = shellTypeToNumBasis(la);
  const numB = shellTypeToNumBasis(lb);

  for (let ia = 0; ia < numA; ia++) {
    const [lxa, lya, lza] = angA[ia];
    const muA = grpA.basisIndex + ia;
    const normA = normFactors[muA];

    // When same group, only iterate ib >= ia to avoid double counting
    const ibStart = sameGroup ? ia : 0;
    for (let ib = ibStart; ib < numB; ib++) {
      const [lxb, lyb, lzb] = angB[ib];
      const muB = grpB.basisIndex + ib;
      const normB = normFactors[muB];

      let sVal = 0, tVal = 0, vVal = 0;

      // Contract over primitive pairs
      for (const primA of grpA.primitives) {
        const pNormA = primitiveNorm(primA.exponent, lxa, lya, lza);
        for (const primB of grpB.primitives) {
          const pNormB = primitiveNorm(primB.exponent, lxb, lyb, lzb);
          const alpha = primA.exponent;
          const beta = primB.exponent;
          const p = alpha + beta;
          const mu = alpha * beta / p;
          const Ax = primA.coordinate.x, Ay = primA.coordinate.y, Az = primA.coordinate.z;
          const Bx = primB.coordinate.x, By = primB.coordinate.y, Bz = primB.coordinate.z;
          const AB2 = (Ax - Bx) ** 2 + (Ay - By) ** 2 + (Az - Bz) ** 2;
          const Kab = Math.exp(-mu * AB2);
          const coeff = primA.coefficient * primB.coefficient * Kab * pNormA * pNormB;

          const Px = (alpha * Ax + beta * Bx) / p;
          const Py = (alpha * Ay + beta * By) / p;
          const Pz = (alpha * Az + beta * Bz) / p;

          // E-coefficients for each Cartesian direction
          const Ex = eCoefficients(lxa, lxb, p, Px - Ax, Px - Bx);
          const Ey = eCoefficients(lya, lyb, p, Py - Ay, Py - By);
          const Ez = eCoefficients(lza, lzb, p, Pz - Az, Pz - Bz);

          // Overlap: S = (π/p)^{3/2} * Ex[0] * Ey[0] * Ez[0]
          const prefS = Math.pow(Math.PI / p, 1.5);
          sVal += coeff * prefS * Ex[0] * Ey[0] * Ez[0];

          // Kinetic energy: T = -1/2 ∇²
          // T_x direction: β(2lb+1)S(la,lb) - 2β²S(la,lb+2) - lb(lb-1)/2 S(la,lb-2)
          tVal += coeff * prefS * kineticContribution(
            alpha, beta, p,
            lxa, lya, lza, lxb, lyb, lzb,
            Px - Ax, Py - Ay, Pz - Az,
            Px - Bx, Py - By, Pz - Bz,
          );

          // Nuclear attraction: V = -Σ_C Z_C * (2π/p) * Σ E * R
          const prefV = 2 * Math.PI / p;
          for (const atom of atoms) {
            const PCx = Px - atom.coordinate.x;
            const PCy = Py - atom.coordinate.y;
            const PCz = Pz - atom.coordinate.z;
            const RPC2 = PCx * PCx + PCy * PCy + PCz * PCz;

            const maxN = lxa + lya + lza + lxb + lyb + lzb;
            const boys = boysAll(maxN, p * RPC2);
            const R = rFunction(maxN, p, PCx, PCy, PCz, boys);

            let nucSum = 0;
            for (let t = 0; t <= lxa + lxb; t++) {
              for (let u = 0; u <= lya + lyb; u++) {
                for (let v = 0; v <= lza + lzb; v++) {
                  nucSum += Ex[t] * Ey[u] * Ez[v] * rGet(R, t, u, v, 0, maxN);
                }
              }
            }
            vVal += coeff * (-atom.atomicNumber) * prefV * nucSum;
          }
        }
      }

      S.set(muA, muB, S.get(muA, muB) + sVal * normA * normB);
      T.set(muA, muB, T.get(muA, muB) + tVal * normA * normB);
      V.set(muA, muB, V.get(muA, muB) + vVal * normA * normB);

      if (muA !== muB) {
        S.set(muB, muA, S.get(muA, muB));
        T.set(muB, muA, T.get(muA, muB));
        V.set(muB, muA, V.get(muA, muB));
      }
    }
  }
}

/** E-coefficients E^{i,j}_t via McMurchie-Davidson recurrence.
 *  Returns array of size (i+j+1) indexed by t from 0 to i+j. */
export function eCoefficients(
  i: number, j: number, p: number, XPA: number, XPB: number,
): Float64Array {
  const maxT = i + j;
  // E[a][b][t] with 0<=a<=i, 0<=b<=j, 0<=t<=a+b
  const size = (i + 1) * (j + 1) * (maxT + 2);
  const E = new Float64Array(size);
  const idx = (a: number, b: number, t: number) =>
    a * (j + 1) * (maxT + 2) + b * (maxT + 2) + t;

  E[idx(0, 0, 0)] = 1.0;

  // Build up in a direction
  for (let a = 0; a < i; a++) {
    for (let b = 0; b <= j; b++) {
      for (let t = 0; t <= a + b + 1; t++) {
        let val = XPA * (t <= a + b ? E[idx(a, b, t)] : 0);
        if (t > 0) val += (1 / (2 * p)) * E[idx(a, b, t - 1)];
        if (t <= a + b) val += (t + 1) * (t + 1 <= a + b ? E[idx(a, b, t + 1)] : 0);
        E[idx(a + 1, b, t)] += val;
      }
    }
  }

  // Build up in b direction
  for (let b = 0; b < j; b++) {
    for (let t = 0; t <= i + b + 1; t++) {
      let val = XPB * (t <= i + b ? E[idx(i, b, t)] : 0);
      if (t > 0) val += (1 / (2 * p)) * E[idx(i, b, t - 1)];
      if (t <= i + b) val += (t + 1) * (t + 1 <= i + b ? E[idx(i, b, t + 1)] : 0);
      E[idx(i, b + 1, t)] += val;
    }
  }

  // Extract E^{i,j}_t for t = 0..i+j
  const result = new Float64Array(maxT + 1);
  for (let t = 0; t <= maxT; t++) {
    result[t] = E[idx(i, j, t)];
  }
  return result;
}

/** Kinetic energy contribution using overlap integrals with shifted angular momenta */
function kineticContribution(
  _alpha: number, beta: number, p: number,
  lxa: number, lya: number, lza: number,
  lxb: number, lyb: number, lzb: number,
  XPA_x: number, XPA_y: number, XPA_z: number,
  XPB_x: number, XPB_y: number, XPB_z: number,
): number {
  // Note: caller already multiplies by (π/p)^{3/2}, so we don't include it here

  // Helper: overlap integral for given angular momenta
  function overlapComponent(la: number, lb: number, XPA: number, XPB: number): number {
    if (la < 0 || lb < 0) return 0;
    const E = eCoefficients(la, lb, p, XPA, XPB);
    return E[0];
  }

  function S3(aX: number, aY: number, aZ: number, bX: number, bY: number, bZ: number): number {
    if (aX < 0 || aY < 0 || aZ < 0 || bX < 0 || bY < 0 || bZ < 0) return 0;
    return overlapComponent(aX, bX, XPA_x, XPB_x)
      * overlapComponent(aY, bY, XPA_y, XPB_y)
      * overlapComponent(aZ, bZ, XPA_z, XPB_z);
  }

  // T = -1/2 * (d²/dx² + d²/dy² + d²/dz²)
  // For each direction q with angular momentum lb_q:
  // contribution = beta*(2*lb_q+1)*S(la,lb) - 2*beta^2*S(la,lb_q+2) - lb_q*(lb_q-1)/2*S(la,lb_q-2)
  let T = 0;

  // x-direction
  T += beta * (2 * lxb + 1) * S3(lxa, lya, lza, lxb, lyb, lzb)
    - 2 * beta * beta * S3(lxa, lya, lza, lxb + 2, lyb, lzb)
    - 0.5 * lxb * (lxb - 1) * S3(lxa, lya, lza, lxb - 2, lyb, lzb);

  // y-direction
  T += beta * (2 * lyb + 1) * S3(lxa, lya, lza, lxb, lyb, lzb)
    - 2 * beta * beta * S3(lxa, lya, lza, lxb, lyb + 2, lzb)
    - 0.5 * lyb * (lyb - 1) * S3(lxa, lya, lza, lxb, lyb - 2, lzb);

  // z-direction
  T += beta * (2 * lzb + 1) * S3(lxa, lya, lza, lxb, lyb, lzb)
    - 2 * beta * beta * S3(lxa, lya, lza, lxb, lyb, lzb + 2)
    - 0.5 * lzb * (lzb - 1) * S3(lxa, lya, lza, lxb, lyb, lzb - 2);

  return T;
}

/** R-function (Hermite Coulomb integrals) via downward recurrence from Boys function.
 *  Returns a flat array of R^n_{t,u,v} values. */
function rFunction(
  maxOrder: number, p: number,
  PCx: number, PCy: number, PCz: number,
  boys: Float64Array,
): Float64Array {
  // R[n][t][u][v] - store as flat array
  const dim = maxOrder + 1;
  const size = (dim + 1) * (dim + 1) * (dim + 1) * (dim + 1); // n, t, u, v
  const R = new Float64Array(size);

  const rIdx = (n: number, t: number, u: number, v: number) =>
    n * (dim + 1) * (dim + 1) * (dim + 1) + t * (dim + 1) * (dim + 1) + u * (dim + 1) + v;

  // Base: R^n_{0,0,0} = (-2p)^n * F_n(p * |PC|²)
  for (let n = 0; n <= maxOrder; n++) {
    R[rIdx(n, 0, 0, 0)] = Math.pow(-2 * p, n) * boys[n];
  }

  // Recurrence: R^n_{t+1,u,v} = t * R^{n+1}_{t-1,u,v} + XPC * R^{n+1}_{t,u,v}
  for (let n = maxOrder - 1; n >= 0; n--) {
    for (let t = 0; t <= maxOrder - n; t++) {
      for (let u = 0; u <= maxOrder - n - t; u++) {
        for (let v = 0; v <= maxOrder - n - t - u; v++) {
          if (t + u + v === 0) continue; // already set

          if (v > 0) {
            let val = PCz * R[rIdx(n + 1, t, u, v - 1)];
            if (v >= 2) val += (v - 1) * R[rIdx(n + 1, t, u, v - 2)];
            R[rIdx(n, t, u, v)] = val;
          } else if (u > 0) {
            let val = PCy * R[rIdx(n + 1, t, u - 1, v)];
            if (u >= 2) val += (u - 1) * R[rIdx(n + 1, t, u - 2, v)];
            R[rIdx(n, t, u, v)] = val;
          } else {
            let val = PCx * R[rIdx(n + 1, t - 1, u, v)];
            if (t >= 2) val += (t - 1) * R[rIdx(n + 1, t - 2, u, v)];
            R[rIdx(n, t, u, v)] = val;
          }
        }
      }
    }
  }

  return R;
}

function rGet(R: Float64Array, t: number, u: number, v: number, n: number, maxOrder: number): number {
  const dim = maxOrder + 1;
  return R[n * (dim + 1) * (dim + 1) * (dim + 1) + t * (dim + 1) * (dim + 1) + u * (dim + 1) + v];
}
