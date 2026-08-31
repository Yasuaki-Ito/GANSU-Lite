/**
 * Reusable theory selector for demo pages (Walsh, PES Scan, Charges, Basis Set).
 * Provides a compact `<select>` + a one-call HF/DFT builder that handles RI-J
 * for pure DFT functionals.
 *
 * Usage in a demo page:
 *
 *   1. Add the select to the page UI:
 *        <label>Theory: ${theorySelectHTML('theory-sel')}</label>
 *
 *   2. When running SCF, replace
 *        const hf = buildHF(mol, 'RHF') as RHF;
 *      with
 *        const choice = (root.querySelector<HTMLSelectElement>('#theory-sel')!).value as TheoryChoice;
 *        const hf = await buildHFOrDFT(mol, basisSet, choice) as RHF;
 */

import type { DFTConfig, HFMethod } from '../core/builder';
import type { FunctionalName } from '../core/xcFunctional';
import { Molecular } from '../core/molecular';
import { BasisSet } from '../core/basisSet';
import { RIData, generateAutoAuxBasis } from '../core/ri';
import { atomicNumberToElementName } from '../core/constants';
import { buildHF } from '../core/builder';
import type { HF } from '../core/hf';

export type TheoryChoice = 'HF' | FunctionalName;

/** Functionals to expose in the dropdown. SVWN/BLYP/PBE = pure DFT;
 *  B3LYP/PBE0/TPSS = hybrid; CAM-B3LYP/ωB97X-D omitted (RSH-lite simplification
 *  is less interesting for the demo pages). */
const THEORY_OPTIONS: TheoryChoice[] = ['HF', 'SVWN', 'BLYP', 'PBE', 'PBE0', 'B3LYP', 'TPSS'];

const THEORY_TIPS: Record<TheoryChoice, string> = {
  'HF': 'Hartree-Fock (no electron correlation)',
  'SVWN': 'LDA: Slater exchange + VWN5 correlation',
  'BLYP': 'GGA: Becke88 exchange + LYP correlation',
  'PBE': 'GGA: Perdew-Burke-Ernzerhof exchange-correlation',
  'PBE0': 'Hybrid GGA: 25% HF + 75% PBE exchange + PBE correlation',
  'B3LYP': 'Hybrid GGA: 20% HF + 8% Slater + 72% B88 + 19% VWN5 + 81% LYP',
  'TPSS': 'meta-GGA: τ-dependent exchange + PKZB correlation',
  'CAM-B3LYP': 'Range-separated hybrid (RSH-lite simplification)',
  'ωB97X-D': 'Long-range corrected hybrid (RSH-lite simplification)',
};

/** Display label for a theory choice. HF → "HF"; functional → "DFT/<name>". */
export function theoryDisplayLabel(choice: TheoryChoice): string {
  return choice === 'HF' ? 'HF' : `DFT/${choice}`;
}

/** Restricted/unrestricted variant prefix combined with theory.
 *  e.g. "RHF" / "UHF" for HF, "RKS/B3LYP" / "UKS/B3LYP" for DFT. */
export function theoryMethodLabel(choice: TheoryChoice, restricted: boolean): string {
  if (choice === 'HF') return restricted ? 'RHF' : 'UHF';
  return restricted ? `RKS/${choice}` : `UKS/${choice}`;
}

/** Render a `<select>` element for choosing HF or a DFT functional.
 *
 *  Optional `exclude` list lets pages skip slow functionals. TPSS uses
 *  numerical FD for V_xc (~5-10× slower than pure DFT or standard hybrids),
 *  so iterative pages (Geom Opt, Freq Analysis, PES Scan) can omit it. */
export function theorySelectHTML(
  idAttr: string,
  current: TheoryChoice = 'HF',
  extraStyle: string = '',
  exclude: TheoryChoice[] = [],
): string {
  const opts = THEORY_OPTIONS
    .filter(o => !exclude.includes(o))
    .map(o => `<option value="${o}" title="${THEORY_TIPS[o] ?? ''}"${o === current ? ' selected' : ''}>${theoryDisplayLabel(o)}</option>`)
    .join('');
  return `<select id="${idAttr}" class="theory-select" style="${extraStyle}">${opts}</select>`;
}

/** Iterative pages exclude TPSS — its FD-based V_xc is too slow for many SCF runs. */
export const HEAVY_ITERATIVE_EXCLUDE: TheoryChoice[] = ['TPSS'];

/** Build an HF or DFT (KS) instance from a theory choice. For pure DFT
 *  (exact-exchange = 0) it auto-generates a Coulomb-fitting auxiliary basis
 *  and switches to RI-J — same path as the worker uses. */
export async function buildHFOrDFT(
  mol: Molecular, basisSet: BasisSet, theory: TheoryChoice,
  method: HFMethod = 'RHF',
  gridLevel: NonNullable<DFTConfig['gridLevel']> = 'medium',
): Promise<HF> {
  const dftConfig: DFTConfig | undefined = (theory === 'HF')
    ? undefined
    : { functional: theory, gridLevel };
  const hf = buildHF(mol, method, dftConfig);
  if (hf.xcFunctional && hf.xcFunctional.exactExchangeFraction === 0) {
    const elementNames = [...new Set(mol.atoms.map(a => atomicNumberToElementName(a.atomicNumber)))];
    const auxBasis = generateAutoAuxBasis(basisSet, elementNames);
    const ri = await RIData.build(
      mol.primitiveShells, mol.cgtoNormalizationFactors, mol.numBasis, mol.atoms, auxBasis,
    );
    (hf as unknown as { setRIData: (r: RIData) => void }).setRIData(ri);
  }
  return hf;
}

/** Compact label with the dropdown, for header rows. */
export function theoryControlHTML(idAttr: string, current: TheoryChoice = 'HF'): string {
  return `<label class="theory-ctrl">Theory ${theorySelectHTML(idAttr, current)}</label>`;
}
