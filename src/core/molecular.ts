/** Molecular class — port of GANSU molecular.hpp */

import type { Atom, PrimitiveShell, ShellTypeInfo, BasisRange } from './types';
import { atomicNumberToElementName, shellNameToType, shellTypeToNumBasis } from './constants';
import { BasisSet } from './basisSet';
import { Matrix } from '../linalg/matrix';
import { buildCartesianToSpherical, type ContractedShellLayout } from './sphericalBasis';

export class Molecular {
  readonly atoms: Atom[];
  readonly primitiveShells: PrimitiveShell[];
  readonly shellTypeInfos: ShellTypeInfo[];
  readonly atomToBasisRange: BasisRange[];
  readonly numBasis: number;
  readonly cgtoNormalizationFactors: number[];
  readonly numElectrons: number;
  readonly numAlphaSpins: number;
  readonly numBetaSpins: number;
  /** Contracted shells in basis order — what a Cartesian→spherical map is built from. */
  readonly contractedShells: ContractedShellLayout[];
  /**
   * Cartesian → spherical map T (numBasis × numSpherical) when the basis set is
   * used with pure d/f functions; null for a Cartesian basis or one without d/f.
   * numBasis always counts Cartesian functions, because every integral is
   * computed over them; the spherical restriction is applied in the SCF.
   */
  readonly sphericalTransform: Matrix | null;
  /** Functions in the variational basis: the spherical count if pure, else numBasis. */
  readonly numSpherical: number;

  constructor(atoms: Atom[], basisSet: BasisSet, charge = 0, betaToAlpha = 0) {
    this.atoms = atoms;
    if (atoms.length === 0) throw new Error('No atoms given');

    // Build primitive shells and normalization factors
    const primitiveShells: PrimitiveShell[] = [];
    const normFactors: number[] = [];
    const basisRanges: BasisRange[] = [];
    const contracted: ContractedShellLayout[] = [];
    let basisIndex = 0;

    for (let ai = 0; ai < atoms.length; ai++) {
      const atom = atoms[ai];
      const startIndex = basisIndex;
      const elementName = atomicNumberToElementName(atom.atomicNumber);
      const ebs = basisSet.get(elementName);

      for (const cg of ebs.contractedGausses) {
        const shellType = shellNameToType(cg.type);
        contracted.push({ basisIndex, shellType });

        for (const prim of cg.primitives) {
          primitiveShells.push({
            exponent: prim.exponent,
            coefficient: prim.coefficient,
            coordinate: { ...atom.coordinate },
            shellType,
            basisIndex,
            atomIndex: ai,
          });
        }

        basisIndex += shellTypeToNumBasis(shellType);
        normFactors.push(...cg.getNormalizationFactors());
      }

      basisRanges.push({ startIndex, endIndex: basisIndex });
    }

    this.numBasis = basisIndex;
    this.cgtoNormalizationFactors = normFactors;
    this.atomToBasisRange = basisRanges;
    this.contractedShells = contracted;

    const sph = basisSet.pure ? buildCartesianToSpherical(contracted, basisIndex) : null;
    this.sphericalTransform = sph?.T ?? null;
    this.numSpherical = sph?.numSpherical ?? basisIndex;

    if (this.numBasis !== normFactors.length) {
      throw new Error(`Basis count mismatch: numBasis=${this.numBasis}, normFactors=${normFactors.length}`);
    }

    // Sort primitive shells by shell type, then by basis index
    primitiveShells.sort((a, b) =>
      a.shellType !== b.shellType ? a.shellType - b.shellType : a.basisIndex - b.basisIndex
    );

    // Compute shell type infos
    const maxShellType = primitiveShells[primitiveShells.length - 1].shellType;
    const infos: ShellTypeInfo[] = [];
    for (let t = 0; t <= maxShellType; t++) {
      infos.push({ count: 0, startIndex: 0 });
    }
    for (const ps of primitiveShells) {
      infos[ps.shellType].count++;
    }
    infos[0].startIndex = 0;
    for (let t = 1; t < infos.length; t++) {
      infos[t].startIndex = infos[t - 1].startIndex + infos[t - 1].count;
    }

    this.primitiveShells = primitiveShells;
    this.shellTypeInfos = infos;

    // Electron count
    let totalZ = 0;
    for (const atom of atoms) totalZ += atom.atomicNumber;
    this.numElectrons = totalZ - charge;
    if (this.numElectrons < 1) throw new Error('Number of electrons is less than 1');

    this.numAlphaSpins = Math.ceil(this.numElectrons / 2) + betaToAlpha;
    this.numBetaSpins = Math.floor(this.numElectrons / 2) - betaToAlpha;
    if (this.numBetaSpins < 0) throw new Error('Number of beta-spin electrons is negative');
  }
}
