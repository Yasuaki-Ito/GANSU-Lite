/** Dipole integrals via McMurchie-Davidson.
 *
 *  <μ|r_q|ν> for q = x, y, z
 *  Using identity: r_q = (r_q - B_q) + B_q
 *  so the integral becomes: E^{la,lb+1}_0 + B_q * E^{la,lb}_0 (times overlap in other directions)
 */

import type { PrimitiveShell } from './types';
import { ANGULAR_MOMENTUMS, shellTypeToNumBasis } from './constants';
import { Matrix } from '../linalg/matrix';
import { eCoefficients, primitiveNorm } from './integrals1e';

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

/** Compute dipole integral matrices Dx, Dy, Dz in AO basis (Bohr coordinates). */
export function computeDipoleIntegrals(
  primitiveShells: PrimitiveShell[],
  normFactors: number[],
  numBasis: number,
): { Dx: Matrix; Dy: Matrix; Dz: Matrix } {
  const Dx = new Matrix(numBasis, numBasis);
  const Dy = new Matrix(numBasis, numBasis);
  const Dz = new Matrix(numBasis, numBasis);

  const shellGroups = groupPrimitiveShells(primitiveShells);

  for (let iGrp = 0; iGrp < shellGroups.length; iGrp++) {
    for (let jGrp = iGrp; jGrp < shellGroups.length; jGrp++) {
      const grpA = shellGroups[iGrp];
      const grpB = shellGroups[jGrp];
      computeShellPairDipole(grpA, grpB, normFactors, Dx, Dy, Dz, iGrp === jGrp);
    }
  }

  return { Dx, Dy, Dz };
}

function computeShellPairDipole(
  grpA: ShellGroup,
  grpB: ShellGroup,
  normFactors: number[],
  Dx: Matrix,
  Dy: Matrix,
  Dz: Matrix,
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

    const ibStart = sameGroup ? ia : 0;
    for (let ib = ibStart; ib < numB; ib++) {
      const [lxb, lyb, lzb] = angB[ib];
      const muB = grpB.basisIndex + ib;
      const normB = normFactors[muB];

      let dx = 0, dy = 0, dz = 0;

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

          const Ex = eCoefficients(lxa, lxb, p, Px - Ax, Px - Bx);
          const Ey = eCoefficients(lya, lyb, p, Py - Ay, Py - By);
          const Ez = eCoefficients(lza, lzb, p, Pz - Az, Pz - Bz);

          const prefS = Math.pow(Math.PI / p, 1.5);

          // <μ|x|ν> = (π/p)^{3/2} * [E^{lxa, lxb+1}_0 + Bx * E^{lxa, lxb}_0] * Ey_0 * Ez_0
          const ExUp = eCoefficients(lxa, lxb + 1, p, Px - Ax, Px - Bx);
          dx += coeff * prefS * (ExUp[0] + Bx * Ex[0]) * Ey[0] * Ez[0];

          // <μ|y|ν>
          const EyUp = eCoefficients(lya, lyb + 1, p, Py - Ay, Py - By);
          dy += coeff * prefS * Ex[0] * (EyUp[0] + By * Ey[0]) * Ez[0];

          // <μ|z|ν>
          const EzUp = eCoefficients(lza, lzb + 1, p, Pz - Az, Pz - Bz);
          dz += coeff * prefS * Ex[0] * Ey[0] * (EzUp[0] + Bz * Ez[0]);
        }
      }

      const valX = dx * normA * normB;
      const valY = dy * normA * normB;
      const valZ = dz * normA * normB;

      Dx.set(muA, muB, Dx.get(muA, muB) + valX);
      Dy.set(muA, muB, Dy.get(muA, muB) + valY);
      Dz.set(muA, muB, Dz.get(muA, muB) + valZ);

      if (muA !== muB) {
        Dx.set(muB, muA, Dx.get(muA, muB));
        Dy.set(muB, muA, Dy.get(muA, muB));
        Dz.set(muB, muA, Dz.get(muA, muB));
      }
    }
  }
}
