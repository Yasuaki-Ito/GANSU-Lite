/** HF/KS Builder factory — port of GANSU builder.hpp */

import { HF } from './hf';
import { RHF } from './rhf';
import { UHF } from './uhf';
import { ROHF } from './rohf';
import { Molecular } from './molecular';
import { createFunctional, type FunctionalName } from './xcFunctional';
import { buildMolecularGrid, type GridLevel } from './grid';
export type HFMethod = 'RHF' | 'UHF' | 'ROHF';

export interface DFTConfig {
  functional: FunctionalName;
  gridLevel?: GridLevel;
}

function buildFromMol(mol: Molecular, Constructor: new (...args: ConstructorParameters<typeof RHF>) => HF): HF {
  return new Constructor(
    mol.numBasis,
    mol.numElectrons,
    mol.numAlphaSpins,
    mol.numBetaSpins,
    mol.atoms,
    mol.primitiveShells,
    mol.shellTypeInfos,
    mol.atomToBasisRange,
    mol.cgtoNormalizationFactors,
  );
}

export function buildHF(mol: Molecular, method: HFMethod = 'RHF', dft?: DFTConfig): HF {
  let hf: HF;
  switch (method) {
    case 'RHF':
      hf = buildFromMol(mol, RHF);
      break;
    case 'UHF':
      hf = buildFromMol(mol, UHF);
      break;
    case 'ROHF':
      hf = buildFromMol(mol, ROHF);
      break;
    default:
      throw new Error(`Unknown HF method: ${method}`);
  }

  if (dft) {
    const functional = createFunctional(dft.functional);
    const targetLevel = dft.gridLevel ?? 'medium';
    const grid = buildMolecularGrid(mol.atoms, targetLevel);
    // Build coarse grid for initial SCF iterations if target is not already coarse
    const coarseGrid = targetLevel !== 'coarse'
      ? buildMolecularGrid(mol.atoms, 'coarse')
      : undefined;
    hf.setDFT(functional, grid, coarseGrid);
  }

  return hf;
}
