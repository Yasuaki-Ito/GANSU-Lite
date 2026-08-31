/** SAD (Superposition of Atomic Densities) initial guess for SCF.
 *
 *  Builds a diagonal density matrix from atomic ground-state electron
 *  configurations. Each basis function gets a fractional occupation based
 *  on its angular momentum type and the atom's electron count.
 *
 *  For split-valence basis sets, electrons of each angular momentum type
 *  are distributed equally among all contracted gaussians of that type.
 */

import type { Atom, BasisRange } from './types';
import { BasisSet } from './basisSet';
import { atomicNumberToElementName, shellNameToType, shellTypeToNumBasis } from './constants';
import { Matrix } from '../linalg/matrix';

/** Atomic ground-state electron counts by angular momentum.
 *  s: total s-electrons, p: total p-electrons, d: total d-electrons */
const ATOMIC_ELECTRONS: Record<number, { s: number; p: number; d: number }> = {
  1:  { s: 1, p: 0, d: 0 },  // H:  1s¹
  2:  { s: 2, p: 0, d: 0 },  // He: 1s²
  3:  { s: 3, p: 0, d: 0 },  // Li: [He] 2s¹
  4:  { s: 4, p: 0, d: 0 },  // Be: [He] 2s²
  5:  { s: 4, p: 1, d: 0 },  // B:  [He] 2s² 2p¹
  6:  { s: 4, p: 2, d: 0 },  // C:  [He] 2s² 2p²
  7:  { s: 4, p: 3, d: 0 },  // N:  [He] 2s² 2p³
  8:  { s: 4, p: 4, d: 0 },  // O:  [He] 2s² 2p⁴
  9:  { s: 4, p: 5, d: 0 },  // F:  [He] 2s² 2p⁵
  10: { s: 4, p: 6, d: 0 },  // Ne: [He] 2s² 2p⁶
  11: { s: 5, p: 6, d: 0 },  // Na: [Ne] 3s¹
  12: { s: 6, p: 6, d: 0 },  // Mg: [Ne] 3s²
  13: { s: 6, p: 7, d: 0 },  // Al: [Ne] 3s² 3p¹
  14: { s: 6, p: 8, d: 0 },  // Si: [Ne] 3s² 3p²
  15: { s: 6, p: 9, d: 0 },  // P:  [Ne] 3s² 3p³
  16: { s: 6, p: 10, d: 0 }, // S:  [Ne] 3s² 3p⁴
  17: { s: 6, p: 11, d: 0 }, // Cl: [Ne] 3s² 3p⁵
  18: { s: 6, p: 12, d: 0 }, // Ar: [Ne] 3s² 3p⁶
  19: { s: 7, p: 12, d: 0 }, // K:  [Ar] 4s¹
  20: { s: 8, p: 12, d: 0 }, // Ca: [Ar] 4s²
  21: { s: 8, p: 12, d: 1 }, // Sc: [Ar] 3d¹ 4s²
  22: { s: 8, p: 12, d: 2 }, // Ti: [Ar] 3d² 4s²
  23: { s: 8, p: 12, d: 3 }, // V:  [Ar] 3d³ 4s²
  24: { s: 7, p: 12, d: 5 }, // Cr: [Ar] 3d⁵ 4s¹
  25: { s: 8, p: 12, d: 5 }, // Mn: [Ar] 3d⁵ 4s²
  26: { s: 8, p: 12, d: 6 }, // Fe: [Ar] 3d⁶ 4s²
  27: { s: 8, p: 12, d: 7 }, // Co: [Ar] 3d⁷ 4s²
  28: { s: 8, p: 12, d: 8 }, // Ni: [Ar] 3d⁸ 4s²
  29: { s: 7, p: 12, d: 10 }, // Cu: [Ar] 3d¹⁰ 4s¹
  30: { s: 8, p: 12, d: 10 }, // Zn: [Ar] 3d¹⁰ 4s²
  31: { s: 8, p: 13, d: 10 }, // Ga
  32: { s: 8, p: 14, d: 10 }, // Ge
  33: { s: 8, p: 15, d: 10 }, // As
  34: { s: 8, p: 16, d: 10 }, // Se
  35: { s: 8, p: 17, d: 10 }, // Br
  36: { s: 8, p: 18, d: 10 }, // Kr
};

/** Build SOAD initial density matrix.
 *
 *  For each atom, distributes electrons equally among basis functions
 *  of the same angular momentum type:
 *    occupation_per_s_function = (total_s_electrons) / (num_s_contracted_gaussians)
 *    occupation_per_p_function = (total_p_electrons) / (num_p_contracted_gaussians × 3)
 *    etc.
 *
 *  Returns a diagonal matrix (total density P = Pα + Pβ).
 */
export function computeSOADDensity(
  atoms: Atom[],
  basisSet: BasisSet,
  numBasis: number,
  atomToBasisRange: BasisRange[],
): Matrix {
  const P = new Matrix(numBasis, numBasis);

  for (let ai = 0; ai < atoms.length; ai++) {
    const Z = atoms[ai].atomicNumber;
    const config = ATOMIC_ELECTRONS[Z];
    if (!config) continue; // unsupported element

    const elementName = atomicNumberToElementName(Z);
    const ebs = basisSet.get(elementName);

    // Count contracted gaussians by shell type
    let ns = 0, np = 0, nd = 0;
    for (const cg of ebs.contractedGausses) {
      const t = cg.type.toUpperCase();
      if (t === 'S') ns++;
      else if (t === 'P') np++;
      else if (t === 'D') nd++;
    }

    // Occupation per basis function
    const occS = ns > 0 ? config.s / ns : 0;
    const occP = np > 0 ? config.p / (np * 3) : 0;
    const occD = nd > 0 ? config.d / (nd * 6) : 0;

    // Assign diagonal elements in the same order as Molecular constructor
    let mu = atomToBasisRange[ai].startIndex;
    for (const cg of ebs.contractedGausses) {
      const shellType = shellNameToType(cg.type);
      const nfuncs = shellTypeToNumBasis(shellType);
      let occ: number;
      if (shellType === 0) occ = occS;
      else if (shellType === 1) occ = occP;
      else if (shellType === 2) occ = occD;
      else occ = 0;

      for (let k = 0; k < nfuncs; k++) {
        P.set(mu + k, mu + k, occ);
      }
      mu += nfuncs;
    }
  }

  return P;
}
