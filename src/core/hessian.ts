/** Numerical Hessian (via finite differences of analytic gradient)
 *  and harmonic vibrational frequency analysis. */

import type { Atom, PrimitiveShell } from './types';
import { Molecular } from './molecular';
import { BasisSet } from './basisSet';
import { buildHF, type DFTConfig } from './builder';
import { RHF } from './rhf';
import { computeRHFGradient } from './gradient';
import { computeDipoleMoment } from './properties';
import { computeFullAnalyticalHessian } from './analyticalHessian';
import type { ERIStored } from './eri';
import { Matrix } from '../linalg/matrix';
import { atomicNumberToElementName } from './constants';
import { RIData, generateAutoAuxBasis } from './ri';

// ── Atomic masses (amu, most abundant isotope) ───────────────────────

const ATOMIC_MASSES: Record<number, number> = {
  1: 1.00794, 2: 4.00260, 3: 6.941, 4: 9.01218, 5: 10.811,
  6: 12.011, 7: 14.007, 8: 15.999, 9: 18.998, 10: 20.180,
  11: 22.990, 12: 24.305, 13: 26.982, 14: 28.086, 15: 30.974,
  16: 32.065, 17: 35.453, 18: 39.948,
};

export function getAtomicMass(Z: number): number {
  return ATOMIC_MASSES[Z] ?? Z * 2;
}

// ── Constants ────────────────────────────────────────────────────────

const BOHR_TO_M = 5.29177210903e-11;
const HARTREE_TO_J = 4.3597447222071e-18;
const AMU_TO_KG = 1.66053906660e-27;
const SPEED_OF_LIGHT = 2.99792458e10; // cm/s
const TWOPI = 2 * Math.PI;

// ── Numerical Hessian ────────────────────────────────────────────────

export interface ThermoData {
  temperature: number;     // K
  pressure: number;        // atm
  E_elec: number;         // Hartree (electronic energy, input)
  ZPE: number;            // kcal/mol
  E_tot: number;          // Hartree (E_elec + thermal corrections)
  H: number;              // Hartree (enthalpy)
  G: number;              // Hartree (Gibbs free energy)
  S_tot: number;          // cal/(mol·K) total entropy
  S_trans: number;        // translational entropy
  S_rot: number;          // rotational entropy
  S_vib: number;          // vibrational entropy
  S_elec: number;         // electronic entropy
  Cv_tot: number;         // cal/(mol·K) heat capacity
  thermalCorr: number;    // Hartree (thermal correction to energy)
  enthalpyCorr: number;   // Hartree (thermal correction to enthalpy)
  gibbsCorr: number;      // Hartree (thermal correction to Gibbs)
}

export interface HessianResult {
  hessian: Float64Array;       // 3N × 3N flat
  n3: number;                  // 3N
  frequencies: number[];       // cm⁻¹ (sorted, negative = imaginary)
  modes: Float64Array[];       // eigenvectors in mass-weighted coords
  zpe: number;                 // zero-point energy in Hartree
  masses: number[];            // amu per atom
  intensities: number[];       // IR intensities (km/mol)
  thermo?: ThermoData;         // thermodynamic quantities
}

/** Compute the Hessian by central finite differences of the analytic gradient.
 *  For N atoms, requires 6N SCF+gradient evaluations.
 *  When dftConfig is provided, uses RKS gradient (analytic). */
export async function computeNumericalHessian(
  atomicNumbers: number[],
  refCoords: Float64Array,
  basis: BasisSet,
  charge: number,
  h: number = 1e-3,
  onProgress?: (msg: string) => void,
  dftConfig?: DFTConfig,
): Promise<HessianResult> {
  const nAtom = atomicNumbers.length;
  const n3 = 3 * nAtom;

  /** Run SCF + gradient + dipole at given coordinates */
  async function evalGradDipole(coords: Float64Array): Promise<{ grad: Float64Array; dipole: [number, number, number] }> {
    const atoms: Atom[] = atomicNumbers.map((z, i) => ({
      atomicNumber: z,
      coordinate: { x: coords[3 * i], y: coords[3 * i + 1], z: coords[3 * i + 2] },
      atomIndex: i,
    }));
    const mol = new Molecular(atoms, basis, charge);
    const rhf = buildHF(mol, 'RHF', dftConfig) as RHF;
    // For pure DFT, set up RI-J auxiliary basis (matches main worker path)
    if (rhf.xcFunctional && rhf.xcFunctional.exactExchangeFraction === 0) {
      const elementNames = [...new Set(mol.atoms.map(a => atomicNumberToElementName(a.atomicNumber)))];
      const auxBasis = generateAutoAuxBasis(basis, elementNames);
      const riData = await RIData.build(
        mol.primitiveShells, mol.cgtoNormalizationFactors, mol.numBasis, mol.atoms, auxBasis,
      );
      rhf.setRIData(riData);
    }
    await rhf.solve({ eriBackend: 'auto' });
    const nocc = mol.numAlphaSpins;
    const dftCtx = (dftConfig && rhf.xcFunctional && rhf.grid)
      ? { functional: rhf.xcFunctional, grid: rhf.grid }
      : undefined;
    const gResult = computeRHFGradient(
      mol.primitiveShells, mol.atoms, mol.cgtoNormalizationFactors,
      mol.numBasis, nocc, rhf.density, rhf.coefficients, rhf.orbitalEnergies,
      undefined, dftCtx,
    );
    const dip = computeDipoleMoment(rhf.density, mol.atoms, mol.primitiveShells, mol.cgtoNormalizationFactors, mol.numBasis);
    return { grad: gResult.total, dipole: [dip.x, dip.y, dip.z] };
  }

  const hessian = new Float64Array(n3 * n3);
  // dμ/dR: dipole derivative matrix (3 × 3N)
  const dDipole = new Float64Array(3 * n3); // [dμx/dR0, dμy/dR0, dμz/dR0, dμx/dR1, ...]

  for (let i = 0; i < n3; i++) {
    onProgress?.(`Displacing ${i + 1}/${n3}`);
    // Yield to UI so progress updates are painted
    await new Promise<void>(r => setTimeout(r, 0));

    const coordsP = Float64Array.from(refCoords);
    const coordsM = Float64Array.from(refCoords);
    coordsP[i] += h;
    coordsM[i] -= h;

    const [resP, resM] = await Promise.all([evalGradDipole(coordsP), evalGradDipole(coordsM)]);

    for (let j = 0; j < n3; j++) {
      hessian[i * n3 + j] = (resP.grad[j] - resM.grad[j]) / (2 * h);
    }
    // Dipole derivatives
    dDipole[3 * i]     = (resP.dipole[0] - resM.dipole[0]) / (2 * h);
    dDipole[3 * i + 1] = (resP.dipole[1] - resM.dipole[1]) / (2 * h);
    dDipole[3 * i + 2] = (resP.dipole[2] - resM.dipole[2]) / (2 * h);
  }

  // Symmetrize
  for (let i = 0; i < n3; i++)
    for (let j = i + 1; j < n3; j++) {
      const avg = 0.5 * (hessian[i * n3 + j] + hessian[j * n3 + i]);
      hessian[i * n3 + j] = avg;
      hessian[j * n3 + i] = avg;
    }

  // ── Vibrational analysis ──
  const masses = atomicNumbers.map(getAtomicMass);
  const { frequencies, modes, zpe, sortedEigvecs } = diagonalizeMassWeighted(hessian, n3, masses, refCoords);

  // ── IR intensities ──
  // I_k ∝ |∂μ/∂Q_k|² where ∂μ_α/∂Q_k = Σ_i (∂μ_α/∂R_i) * V_{ik} / √m_i
  // Use the SAME eigenvectors as the frequency calculation
  const intensities: number[] = [];
  const AU2_TO_KMMOL = 974.9; // (e/√amu)² → km/mol

  for (let ki = 0; ki < frequencies.length; ki++) {
    if (Math.abs(frequencies[ki]) < 50) { intensities.push(0); continue; }
    let dmu_x = 0, dmu_y = 0, dmu_z = 0;
    for (let i = 0; i < n3; i++) {
      const Vik = sortedEigvecs[ki * n3 + i]; // ki-th sorted eigenvector, i-th component
      const invSqrtM = 1 / Math.sqrt(masses[Math.floor(i / 3)]);
      dmu_x += dDipole[3 * i]     * Vik * invSqrtM;
      dmu_y += dDipole[3 * i + 1] * Vik * invSqrtM;
      dmu_z += dDipole[3 * i + 2] * Vik * invSqrtM;
    }
    intensities.push((dmu_x * dmu_x + dmu_y * dmu_y + dmu_z * dmu_z) * AU2_TO_KMMOL);
  }

  // ── Thermodynamic quantities ──
  const thermo = computeThermodynamics(
    frequencies, masses, refCoords, atomicNumbers.length, 0, // E_elec=0 placeholder, caller fills in
    298.15, 1.0,
  );

  return { hessian, n3, frequencies, modes, zpe, masses, intensities, thermo };
}

/** Compute Hessian analytically (skeleton + CPHF).
 *  Faster than numerical for 5+ atoms, same accuracy. */
export async function computeAnalyticalHessianFull(
  atomicNumbers: number[],
  refCoords: Float64Array,
  basis: BasisSet,
  charge: number,
  onProgress?: (msg: string) => void,
): Promise<HessianResult> {
  const nAtom = atomicNumbers.length;
  const n3 = 3 * nAtom;

  // Run SCF for Hessian reference
  onProgress?.('SCF for reference...');
  await new Promise<void>(r => setTimeout(r, 0));
  const atoms: Atom[] = atomicNumbers.map((z, i) => ({
    atomicNumber: z,
    coordinate: { x: refCoords[3*i], y: refCoords[3*i+1], z: refCoords[3*i+2] },
    atomIndex: i,
  }));
  const mol = new Molecular(atoms, basis, charge);
  const rhf = buildHF(mol, 'RHF') as RHF;
  const energy = await rhf.solve({ eriBackend: 'auto' });
  const nocc = mol.numAlphaSpins, n = mol.numBasis;

  // Analytical Hessian
  const hessian = await computeFullAnalyticalHessian(
    mol.primitiveShells, mol.atoms, mol.cgtoNormalizationFactors, n,
    nocc, rhf.density, rhf.coefficients, rhf.orbitalEnergies, rhf.eriStore,
    onProgress,
  );

  // Vibrational analysis
  const masses = atomicNumbers.map(getAtomicMass);
  const { frequencies, modes, zpe, sortedEigvecs } = diagonalizeMassWeighted(hessian, n3, masses, refCoords);

  // IR intensities via dipole FD (still needed — analytical dipole derivs require CPHF dipole response)
  // Use a lighter approach: FD of dipole only (no SCF reconvergence, just rebuild integrals)
  onProgress?.('IR intensities...');
  await new Promise<void>(r => setTimeout(r, 0));
  const h = 5e-4;
  const dDipole = new Float64Array(3 * n3);
  for (let i = 0; i < n3; i++) {
    const coordsP = Float64Array.from(refCoords); coordsP[i] += h;
    const coordsM = Float64Array.from(refCoords); coordsM[i] -= h;
    const atomsP = atomicNumbers.map((z, j) => ({ atomicNumber: z, coordinate: { x: coordsP[3*j], y: coordsP[3*j+1], z: coordsP[3*j+2] }, atomIndex: j }));
    const atomsM = atomicNumbers.map((z, j) => ({ atomicNumber: z, coordinate: { x: coordsM[3*j], y: coordsM[3*j+1], z: coordsM[3*j+2] }, atomIndex: j }));
    const molP = new Molecular(atomsP, basis, charge), molM = new Molecular(atomsM, basis, charge);
    // Use reference density (fixed D) for dipole at displaced geometry
    const dipP = computeDipoleMoment(rhf.density, molP.atoms, molP.primitiveShells, molP.cgtoNormalizationFactors, n);
    const dipM = computeDipoleMoment(rhf.density, molM.atoms, molM.primitiveShells, molM.cgtoNormalizationFactors, n);
    dDipole[3*i]   = (dipP.x - dipM.x) / (2*h);
    dDipole[3*i+1] = (dipP.y - dipM.y) / (2*h);
    dDipole[3*i+2] = (dipP.z - dipM.z) / (2*h);
  }

  const intensities: number[] = [];
  const AU2_TO_KMMOL = 974.9;
  for (let ki = 0; ki < frequencies.length; ki++) {
    if (Math.abs(frequencies[ki]) < 50) { intensities.push(0); continue; }
    let dmu_x = 0, dmu_y = 0, dmu_z = 0;
    for (let i = 0; i < n3; i++) {
      const Vik = sortedEigvecs[ki * n3 + i];
      const invSqrtM = 1 / Math.sqrt(masses[Math.floor(i / 3)]);
      dmu_x += dDipole[3*i]   * Vik * invSqrtM;
      dmu_y += dDipole[3*i+1] * Vik * invSqrtM;
      dmu_z += dDipole[3*i+2] * Vik * invSqrtM;
    }
    intensities.push((dmu_x*dmu_x + dmu_y*dmu_y + dmu_z*dmu_z) * AU2_TO_KMMOL);
  }

  const thermo = computeThermodynamics(frequencies, masses, refCoords, nAtom, energy, 298.15, 1.0);
  return { hessian, n3, frequencies, modes, zpe, masses, intensities, thermo };
}

/** Auto-select numerical or analytical Hessian based on molecule/basis size.
 *  CPHF overhead (MO-ERI O(N⁵), RHS build, CG) makes analytical slower
 *  for small basis sets. Use nbasis as criterion:
 *    nbasis < 20: numerical (SCF is cheap, 6N evaluations still fast)
 *    nbasis >= 20: analytical (CPHF amortises better)
 *  When dftConfig is provided, force numerical (DFT analytical Hessian needs CPKS+XC kernel — not implemented). */
export async function computeHessianAuto(
  atomicNumbers: number[],
  refCoords: Float64Array,
  basis: BasisSet,
  charge: number,
  h: number = 5e-4,
  onProgress?: (msg: string) => void,
  dftConfig?: DFTConfig,
): Promise<HessianResult> {
  const nAtom = atomicNumbers.length;
  // Estimate nbasis to choose method
  const atoms: Atom[] = atomicNumbers.map((z, i) => ({
    atomicNumber: z,
    coordinate: { x: refCoords[3*i], y: refCoords[3*i+1], z: refCoords[3*i+2] },
    atomIndex: i,
  }));
  const mol = new Molecular(atoms, basis, charge);
  const nbasis = mol.numBasis;

  if (dftConfig) {
    onProgress?.(`Numerical KS (${nAtom} atoms, ${nbasis} basis, ${dftConfig.functional})`);
    return computeNumericalHessian(atomicNumbers, refCoords, basis, charge, h, onProgress, dftConfig);
  } else if (nbasis >= 20) {
    onProgress?.(`Analytical (${nAtom} atoms, ${nbasis} basis)`);
    return computeAnalyticalHessianFull(atomicNumbers, refCoords, basis, charge, onProgress);
  } else {
    onProgress?.(`Numerical (${nAtom} atoms, ${nbasis} basis)`);
    return computeNumericalHessian(atomicNumbers, refCoords, basis, charge, h, onProgress);
  }
}

// ── Mass-weighted Hessian diagonalization ─────────────────────────────

/** Build orthonormal basis for the 6 (or 5 for linear) translation/rotation
 *  modes in mass-weighted coordinates. Returns rows of the basis as Float64Arrays. */
function buildTRBasis(coords: Float64Array | undefined, masses: number[], n3: number): Float64Array[] {
  if (!coords) return [];
  const nAt = masses.length;
  // COM (mass-weighted)
  let M = 0, cx = 0, cy = 0, cz = 0;
  for (let i = 0; i < nAt; i++) { M += masses[i]; cx += masses[i] * coords[3*i]; cy += masses[i] * coords[3*i+1]; cz += masses[i] * coords[3*i+2]; }
  cx /= M; cy /= M; cz /= M;
  // 3 translations and 3 rotations in mass-weighted coords
  // (mass-weighted vector v_i has component v_{3i+d} that already includes √m_i)
  const raw: Float64Array[] = [];
  for (let d = 0; d < 3; d++) {
    const v = new Float64Array(n3);
    for (let i = 0; i < nAt; i++) v[3*i+d] = Math.sqrt(masses[i]);
    raw.push(v);
  }
  for (let axis = 0; axis < 3; axis++) {
    const v = new Float64Array(n3);
    for (let i = 0; i < nAt; i++) {
      const sq = Math.sqrt(masses[i]);
      const rx = coords[3*i] - cx, ry = coords[3*i+1] - cy, rz = coords[3*i+2] - cz;
      if (axis === 0) { v[3*i+1] = -rz * sq; v[3*i+2] = ry * sq; }
      else if (axis === 1) { v[3*i] = rz * sq; v[3*i+2] = -rx * sq; }
      else { v[3*i] = -ry * sq; v[3*i+1] = rx * sq; }
    }
    raw.push(v);
  }
  // Gram-Schmidt
  const orth: Float64Array[] = [];
  for (const r of raw) {
    const v = Float64Array.from(r);
    for (const b of orth) { let d = 0; for (let i = 0; i < n3; i++) d += v[i] * b[i]; for (let i = 0; i < n3; i++) v[i] -= d * b[i]; }
    let nm = 0; for (let i = 0; i < n3; i++) nm += v[i] * v[i]; nm = Math.sqrt(nm);
    if (nm > 1e-10) { for (let i = 0; i < n3; i++) v[i] /= nm; orth.push(v); }
  }
  return orth;
}

function diagonalizeMassWeighted(
  hessian: Float64Array, n3: number, masses: number[], coords?: Float64Array,
): { frequencies: number[]; modes: Float64Array[]; zpe: number; sortedEigvecs: Float64Array } {
  const nAtom = masses.length;

  // Build mass-weighted Hessian: F_ij = H_ij / sqrt(m_i * m_j)
  const mwH = new Float64Array(n3 * n3);
  for (let i = 0; i < n3; i++) {
    const mi = masses[Math.floor(i / 3)];
    for (let j = 0; j < n3; j++) {
      const mj = masses[Math.floor(j / 3)];
      mwH[i * n3 + j] = hessian[i * n3 + j] / Math.sqrt(mi * mj);
    }
  }

  // Project out translations and rotations in mass-weighted coords:
  //   H' = (1 - P) H (1 - P)
  // where P = Σ_k |t_k><t_k| projects onto the TR subspace.
  // This filters spurious low-frequency modes from gradient/grid noise in DFT.
  const trBasis = buildTRBasis(coords, masses, n3);
  if (trBasis.length > 0) {
    // Subtract Σ_k (t_k t_k^T H + H t_k t_k^T - t_k (t_k^T H t_k) t_k^T)
    // Equivalent: build (I-P) and apply twice.
    // Apply P from left: (PH)_ij = Σ_k t_ki Σ_l t_kl H_lj
    // Build (I-P)H first, then ((I-P)H)(I-P).
    // Step A: H -= Σ_k t_k (t_k^T H)
    const tkH: Float64Array[] = [];
    for (const t of trBasis) {
      const v = new Float64Array(n3);
      for (let j = 0; j < n3; j++) { let s = 0; for (let i = 0; i < n3; i++) s += t[i] * mwH[i*n3+j]; v[j] = s; }
      tkH.push(v);
    }
    for (let k = 0; k < trBasis.length; k++) {
      const t = trBasis[k], r = tkH[k];
      for (let i = 0; i < n3; i++) for (let j = 0; j < n3; j++) mwH[i*n3+j] -= t[i] * r[j];
    }
    // Step B: H -= Σ_k (H t_k) t_k^T
    for (const t of trBasis) {
      const v = new Float64Array(n3);
      for (let i = 0; i < n3; i++) { let s = 0; for (let j = 0; j < n3; j++) s += mwH[i*n3+j] * t[j]; v[i] = s; }
      for (let i = 0; i < n3; i++) for (let j = 0; j < n3; j++) mwH[i*n3+j] -= v[i] * t[j];
    }
    // Optional: explicitly zero out the matrix on the TR subspace to ensure
    // diagonalization gives near-zero eigenvalues for those modes.
    // (already zero from the projection; rounding is small)
  }

  // Diagonalize using Jacobi method
  const { eigenvalues, eigenvectors } = jacobiEigen(mwH, n3);

  // Convert eigenvalues to frequencies in cm⁻¹
  // ω² = λ (in atomic units: Eh/(bohr² · amu))
  // Convert: Eh/bohr² → J/m², amu → kg
  // ν(cm⁻¹) = (1/(2πc)) * sqrt(λ * Eh_to_J / (bohr_to_m² * amu_to_kg))
  const convFactor = HARTREE_TO_J / (BOHR_TO_M * BOHR_TO_M * AMU_TO_KG);
  const frequencies: number[] = [];
  const modes: Float64Array[] = [];

  // Sort eigenvalues by absolute value
  const indices = Array.from({ length: n3 }, (_, i) => i);
  indices.sort((a, b) => eigenvalues[a] - eigenvalues[b]);

  for (const idx of indices) {
    const lambda = eigenvalues[idx];
    let freq: number;
    if (lambda >= 0) {
      freq = Math.sqrt(lambda * convFactor) / (TWOPI * SPEED_OF_LIGHT);
    } else {
      freq = -Math.sqrt(-lambda * convFactor) / (TWOPI * SPEED_OF_LIGHT);
    }
    frequencies.push(freq);

    // Extract eigenvector and un-mass-weight
    const mode = new Float64Array(n3);
    for (let i = 0; i < n3; i++) {
      mode[i] = eigenvectors[i * n3 + idx] / Math.sqrt(masses[Math.floor(i / 3)]);
    }
    // Normalize
    let norm = 0;
    for (let i = 0; i < n3; i++) norm += mode[i] * mode[i];
    norm = Math.sqrt(norm);
    if (norm > 1e-14) for (let i = 0; i < n3; i++) mode[i] /= norm;
    modes.push(mode);
  }

  // Build sorted eigenvector matrix for IR intensity calculation
  // sortedEigvecs[k * n3 + i] = V_{i, indices[k]} (i-th component of k-th sorted eigenvector)
  const sortedEigvecs = new Float64Array(n3 * n3);
  for (let k = 0; k < indices.length; k++) {
    const idx = indices[k];
    for (let i = 0; i < n3; i++) {
      sortedEigvecs[k * n3 + i] = eigenvectors[i * n3 + idx];
    }
  }

  // Zero-point energy: ZPE = 0.5 * Σ ℏω (for real frequencies)
  const CM_TO_HARTREE = 4.556335252767e-6;
  let zpe = 0;
  for (const f of frequencies) {
    if (f > 10) zpe += 0.5 * f * CM_TO_HARTREE;
  }

  return { frequencies, modes, zpe, sortedEigvecs };
}

// ── Jacobi eigenvalue algorithm for symmetric matrices ───────────────

// ── Thermodynamic quantities from partition functions ──────────────────

const KB = 1.380649e-23;           // J/K (Boltzmann)
const HBAR_SI = 1.054571817e-34;   // J·s
const H_PLANCK = 6.62607015e-34;   // J·s
const NA = 6.02214076e23;          // /mol
const R_GAS = KB * NA;             // J/(mol·K) = 8.314
const CAL_PER_J = 0.239006;       // cal/J
const HARTREE_PER_J = 1 / HARTREE_TO_J;
const ATM_TO_PA = 101325;

/** Compute thermodynamic quantities at given T, P using ideal gas / rigid rotor / harmonic oscillator. */
export function computeThermodynamics(
  frequencies: number[], masses: number[], coords: Float64Array, nAtom: number,
  E_elec: number, T: number, P: number,
): ThermoData {
  const kT = KB * T;
  const RT = R_GAS * T; // J/mol

  // Total mass
  let totalMass = 0;
  for (const m of masses) totalMass += m;
  const totalMass_kg = totalMass * AMU_TO_KG;

  // ── Translational contributions ──
  // q_trans = (2πmkT/h²)^{3/2} × V/N  (V/N = kT/P for ideal gas)
  const lambda_trans = Math.sqrt(H_PLANCK * H_PLANCK / (TWOPI * totalMass_kg * kT)); // thermal de Broglie
  const q_trans = (kT / (ATM_TO_PA * P)) / (lambda_trans * lambda_trans * lambda_trans);
  const S_trans = R_GAS * (Math.log(q_trans) + 5 / 2) * CAL_PER_J; // cal/(mol·K)  Sackur-Tetrode
  const E_trans = 1.5 * RT; // J/mol
  const Cv_trans = 1.5 * R_GAS * CAL_PER_J; // cal/(mol·K)

  // ── Rotational contributions (rigid rotor) ──
  // Compute principal moments of inertia
  let cx = 0, cy = 0, cz = 0;
  for (let i = 0; i < nAtom; i++) {
    const m = masses[i] * AMU_TO_KG;
    cx += m * coords[3 * i]; cy += m * coords[3 * i + 1]; cz += m * coords[3 * i + 2];
  }
  cx /= totalMass_kg; cy /= totalMass_kg; cz /= totalMass_kg;

  // Inertia tensor (in kg·m²), coords in bohr → convert to m
  const I = [0, 0, 0, 0, 0, 0, 0, 0, 0]; // 3×3
  for (let i = 0; i < nAtom; i++) {
    const m = masses[i] * AMU_TO_KG;
    const x = (coords[3 * i] - cx) * BOHR_TO_M;
    const y = (coords[3 * i + 1] - cy) * BOHR_TO_M;
    const z = (coords[3 * i + 2] - cz) * BOHR_TO_M;
    I[0] += m * (y * y + z * z); I[4] += m * (x * x + z * z); I[8] += m * (x * x + y * y);
    I[1] -= m * x * y; I[2] -= m * x * z; I[5] -= m * y * z;
  }
  I[3] = I[1]; I[6] = I[2]; I[7] = I[5];

  // Diagonalize inertia tensor
  const { eigenvalues: Iprinc } = jacobiEigen(new Float64Array(I), 3);
  const Ia = Math.max(Iprinc[0], 1e-50), Ib = Math.max(Iprinc[1], 1e-50), Ic = Math.max(Iprinc[2], 1e-50);

  // Detect linearity: one principal moment ≈ 0
  const isLinear = nAtom >= 2 && (Ia < 1e-46 || Ib < 1e-46 || Ic < 1e-46);
  const sigma = 1; // rotational symmetry number (simplified; 1 for asymmetric)

  let S_rot = 0, E_rot = 0, Cv_rot = 0;
  if (nAtom === 1) {
    // Atom: no rotation
    S_rot = 0; E_rot = 0; Cv_rot = 0;
  } else if (isLinear) {
    // Linear: q_rot = 8π²IkT/(σh²)
    const Ilin = Math.max(Ia, Ib, Ic); // largest non-zero moment
    const theta_rot = H_PLANCK * H_PLANCK / (8 * Math.PI * Math.PI * Ilin * KB);
    S_rot = R_GAS * (Math.log(T / (sigma * theta_rot)) + 1) * CAL_PER_J;
    E_rot = RT; // J/mol
    Cv_rot = R_GAS * CAL_PER_J;
  } else {
    // Nonlinear: q_rot = √π/(σ) × (8π²kT/h²)^{3/2} × √(Ia×Ib×Ic)
    const theta_A = H_PLANCK * H_PLANCK / (8 * Math.PI * Math.PI * Ia * KB);
    const theta_B = H_PLANCK * H_PLANCK / (8 * Math.PI * Math.PI * Ib * KB);
    const theta_C = H_PLANCK * H_PLANCK / (8 * Math.PI * Math.PI * Ic * KB);
    const q_rot = Math.sqrt(Math.PI) / sigma * Math.pow(T, 1.5) / Math.sqrt(theta_A * theta_B * theta_C);
    S_rot = R_GAS * (Math.log(q_rot) + 1.5) * CAL_PER_J;
    E_rot = 1.5 * RT;
    Cv_rot = 1.5 * R_GAS * CAL_PER_J;
  }

  // ── Vibrational contributions (harmonic oscillator) ──
  let S_vib = 0, E_vib = 0, Cv_vib = 0, zpe_J = 0;
  const vibFreqs = frequencies.filter(f => f > 50); // real vibrational modes only

  for (const nu_cm of vibFreqs) {
    const nu_Hz = nu_cm * SPEED_OF_LIGHT; // cm⁻¹ → Hz (c in cm/s)
    const x = H_PLANCK * nu_Hz / kT; // hν/kT

    zpe_J += 0.5 * H_PLANCK * nu_Hz * NA; // J/mol

    if (x < 500) { // avoid overflow
      const ex = Math.exp(x);
      const exm1 = ex - 1;
      E_vib += RT * x / exm1; // without ZPE
      S_vib += R_GAS * (x / exm1 - Math.log(1 - 1 / ex)) * CAL_PER_J;
      Cv_vib += R_GAS * x * x * ex / (exm1 * exm1) * CAL_PER_J;
    }
  }

  // ── Electronic contribution ──
  const S_elec = 0; // closed-shell singlet: ln(1) = 0
  const E_elec_thermal = 0;
  const Cv_elec = 0;

  // ── Totals ──
  const S_tot = S_trans + S_rot + S_vib + S_elec;
  const Cv_tot = Cv_trans + Cv_rot + Cv_vib + Cv_elec;

  // Thermal correction to energy (in Hartree)
  const ZPE_hartree = zpe_J / NA / HARTREE_TO_J;
  const thermalCorr = (E_trans + E_rot + E_vib) / NA / HARTREE_TO_J + ZPE_hartree;
  const enthalpyCorr = thermalCorr + RT / NA / HARTREE_TO_J; // + kT per molecule
  const TS = T * S_tot / CAL_PER_J / 1000 / NA / HARTREE_TO_J; // T×S in Hartree
  const gibbsCorr = enthalpyCorr - TS;

  return {
    temperature: T,
    pressure: P,
    E_elec,
    ZPE: ZPE_hartree * 627.509, // kcal/mol
    E_tot: E_elec + thermalCorr,
    H: E_elec + enthalpyCorr,
    G: E_elec + gibbsCorr,
    S_tot, S_trans, S_rot, S_vib, S_elec,
    Cv_tot,
    thermalCorr,
    enthalpyCorr,
    gibbsCorr,
  };
}

export function jacobiEigen(A: Float64Array, n: number): { eigenvalues: Float64Array; eigenvectors: Float64Array } {
  // Work on a copy
  const a = Float64Array.from(A);
  const v = new Float64Array(n * n);
  for (let i = 0; i < n; i++) v[i * n + i] = 1; // identity

  const maxIter = 100 * n * n;
  for (let sweep = 0; sweep < maxIter; sweep++) {
    // Find largest off-diagonal element
    let maxVal = 0, p = 0, q = 1;
    for (let i = 0; i < n; i++)
      for (let j = i + 1; j < n; j++) {
        const val = Math.abs(a[i * n + j]);
        if (val > maxVal) { maxVal = val; p = i; q = j; }
      }
    if (maxVal < 1e-14) break;

    // Compute rotation
    const app = a[p * n + p], aqq = a[q * n + q], apq = a[p * n + q];
    const theta = 0.5 * Math.atan2(2 * apq, app - aqq);
    const c = Math.cos(theta), s = Math.sin(theta);

    // Update A
    for (let i = 0; i < n; i++) {
      if (i === p || i === q) continue;
      const aip = a[i * n + p], aiq = a[i * n + q];
      a[i * n + p] = a[p * n + i] = c * aip + s * aiq;
      a[i * n + q] = a[q * n + i] = -s * aip + c * aiq;
    }
    a[p * n + p] = c * c * app + 2 * s * c * apq + s * s * aqq;
    a[q * n + q] = s * s * app - 2 * s * c * apq + c * c * aqq;
    a[p * n + q] = a[q * n + p] = 0;

    // Update eigenvectors
    for (let i = 0; i < n; i++) {
      const vip = v[i * n + p], viq = v[i * n + q];
      v[i * n + p] = c * vip + s * viq;
      v[i * n + q] = -s * vip + c * viq;
    }
  }

  const eigenvalues = new Float64Array(n);
  for (let i = 0; i < n; i++) eigenvalues[i] = a[i * n + i];

  return { eigenvalues, eigenvectors: v };
}
