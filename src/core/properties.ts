/** Post-SCF property calculations */

import type { Atom, BasisRange } from './types';
import { Matrix, createFloatArray, type FloatArray } from '../linalg/matrix';
import { matmul, traceProduct } from '../linalg/matmul';
import { jacobiEigen } from '../linalg/eigendecomposition';
import { computeDipoleIntegrals } from './integralsDipole';
import type { PrimitiveShell } from './types';

/** Mulliken population analysis.
 *  q_A = Z_A - Σ_{μ∈A} (P·S)_{μμ}
 */
export function computeMullikenCharges(
  densityMatrix: Matrix,
  overlapMatrix: Matrix,
  atoms: Atom[],
  atomToBasisRange: BasisRange[],
): FloatArray {
  const PS = matmul(densityMatrix, overlapMatrix);
  const charges = createFloatArray(atoms.length);

  for (let a = 0; a < atoms.length; a++) {
    const { startIndex, endIndex } = atomToBasisRange[a];
    let pop = 0;
    for (let mu = startIndex; mu < endIndex; mu++) {
      pop += PS.get(mu, mu);
    }
    charges[a] = atoms[a].atomicNumber - pop;
  }

  return charges;
}

export interface DipoleMoment {
  x: number;  // a.u.
  y: number;
  z: number;
  total: number;  // a.u.
  debye: number;  // Debye
}

const AU_TO_DEBYE = 2.5417464;

/** Compute electric dipole moment: μ = Σ_A Z_A·R_A - Tr(P·D)
 *  where D_x, D_y, D_z are the dipole integral matrices. */
export function computeDipoleMoment(
  densityMatrix: Matrix,
  atoms: Atom[],
  primitiveShells: PrimitiveShell[],
  normFactors: number[],
  numBasis: number,
): DipoleMoment {
  // Nuclear contribution
  let nucX = 0, nucY = 0, nucZ = 0;
  for (const atom of atoms) {
    nucX += atom.atomicNumber * atom.coordinate.x;
    nucY += atom.atomicNumber * atom.coordinate.y;
    nucZ += atom.atomicNumber * atom.coordinate.z;
  }

  // Electronic contribution: -Tr(P·D)
  const { Dx, Dy, Dz } = computeDipoleIntegrals(primitiveShells, normFactors, numBasis);

  let elecX = 0, elecY = 0, elecZ = 0;
  const n = numBasis;
  for (let i = 0; i < n; i++) {
    for (let k = 0; k < n; k++) {
      const Pik = densityMatrix.get(i, k);
      elecX += Pik * Dx.get(k, i);
      elecY += Pik * Dy.get(k, i);
      elecZ += Pik * Dz.get(k, i);
    }
  }

  const x = nucX - elecX;
  const y = nucY - elecY;
  const z = nucZ - elecZ;
  const total = Math.sqrt(x * x + y * y + z * z);

  return { x, y, z, total, debye: total * AU_TO_DEBYE };
}

// ---------------------------------------------------------------------------
// <S²> expectation value (spin contamination diagnostic for UHF)
// ---------------------------------------------------------------------------

export interface S2Result {
  exact: number;    // S(S+1) where S = (Nα-Nβ)/2
  computed: number; // <S²> from density matrices
}

/** Compute <S²> = S_exact(S_exact+1) + Nβ - Tr(Pα·S·Pβ·S) */
export function computeS2(
  densityAlpha: Matrix, densityBeta: Matrix,
  overlapMatrix: Matrix,
  numAlpha: number, numBeta: number,
): S2Result {
  const S = (numAlpha - numBeta) / 2;
  const exact = S * (S + 1);

  // Tr(Pα·S·Pβ·S)
  const PaS = matmul(densityAlpha, overlapMatrix);
  const PbS = matmul(densityBeta, overlapMatrix);
  const computed = exact + numBeta - traceProduct(PaS, PbS);

  return { exact, computed };
}

// ---------------------------------------------------------------------------
// Löwdin population analysis
// ---------------------------------------------------------------------------

/** Löwdin charges: q_A = Z_A - Σ_{μ∈A} (S^{1/2}·P·S^{1/2})_{μμ} */
export function computeLowdinCharges(
  densityMatrix: Matrix,
  overlapMatrix: Matrix,
  atoms: Atom[],
  atomToBasisRange: BasisRange[],
): FloatArray {
  const n = overlapMatrix.rows;

  // Eigendecompose S: S = U·diag(s)·U^T
  const { eigenvalues, eigenvectors: U } = jacobiEigen(overlapMatrix);

  // Build S^{1/2} = U·diag(√s)·U^T
  const Shalf = new Matrix(n, n);
  for (let i = 0; i < n; i++) {
    const sqrtS = eigenvalues[i] > 1e-10 ? Math.sqrt(eigenvalues[i]) : 0;
    for (let mu = 0; mu < n; mu++) {
      for (let nu = 0; nu < n; nu++) {
        Shalf.add(mu, nu, U.get(mu, i) * sqrtS * U.get(nu, i));
      }
    }
  }

  // P_Löwdin = S^{1/2}·P·S^{1/2}
  const SP = matmul(Shalf, densityMatrix);
  const SPS = matmul(SP, Shalf);

  const charges = createFloatArray(atoms.length);
  for (let a = 0; a < atoms.length; a++) {
    const { startIndex, endIndex } = atomToBasisRange[a];
    let pop = 0;
    for (let mu = startIndex; mu < endIndex; mu++) {
      pop += SPS.get(mu, mu);
    }
    charges[a] = atoms[a].atomicNumber - pop;
  }

  return charges;
}

// ---------------------------------------------------------------------------
// Wiberg/Mayer bond order
// ---------------------------------------------------------------------------

export interface BondOrder {
  i: number;  // atom index
  j: number;  // atom index
  order: number;
}

/** Mayer bond order: B_{AB} = 2·Σ_{μ∈A,ν∈B} [(PαS)_{μν}(PαS)_{νμ} + (PβS)_{μν}(PβS)_{νμ}]
 *  Equivalent to Σ [(PS)_{μν}(PS)_{νμ} + (PsS)_{μν}(PsS)_{νμ}] where P=Pα+Pβ, Ps=Pα-Pβ.
 *  For RHF, pass Pα = Pβ = P_total/2. */
export function computeWibergBondOrder(
  densityAlpha: Matrix, densityBeta: Matrix,
  overlapMatrix: Matrix,
  atoms: Atom[],
  atomToBasisRange: BasisRange[],
  threshold = 0.1,
): BondOrder[] {
  const PaS = matmul(densityAlpha, overlapMatrix);
  const PbS = matmul(densityBeta, overlapMatrix);

  const bonds: BondOrder[] = [];
  for (let a = 0; a < atoms.length; a++) {
    const ra = atomToBasisRange[a];
    for (let b = a + 1; b < atoms.length; b++) {
      const rb = atomToBasisRange[b];
      let order = 0;
      for (let mu = ra.startIndex; mu < ra.endIndex; mu++) {
        for (let nu = rb.startIndex; nu < rb.endIndex; nu++) {
          order += PaS.get(mu, nu) * PaS.get(nu, mu)
                 + PbS.get(mu, nu) * PbS.get(nu, mu);
        }
      }
      order *= 2; // Factor of 2 from Mayer's general formula
      if (order >= threshold) {
        bonds.push({ i: a, j: b, order });
      }
    }
  }

  return bonds;
}

// ---------------------------------------------------------------------------
// Energy component breakdown
// ---------------------------------------------------------------------------

export interface EnergyComponents {
  oneElectron: number;
  twoElectron: number;
  kinetic: number;
  nuclearAttraction: number;
}

/** Compute energy breakdown: 1e/2e/kinetic/nuclear-attraction.
 *  For RHF: pass Fα = Fβ = F, Pα = Pβ = P/2.
 *  For UHF/ROHF: pass actual spin-separated matrices. */
export function computeEnergyComponents(
  densityTotal: Matrix,
  densityAlpha: Matrix, densityBeta: Matrix,
  coreHamiltonian: Matrix, kineticMatrix: Matrix,
  fockAlpha: Matrix, fockBeta: Matrix,
): EnergyComponents {
  const oneElectron = traceProduct(densityTotal, coreHamiltonian);
  const kinetic = traceProduct(densityTotal, kineticMatrix);
  const nuclearAttraction = oneElectron - kinetic;

  // E_2e = 0.5 * (Tr(Pα·(Fα-H)) + Tr(Pβ·(Fβ-H)))
  const n = coreHamiltonian.rows;
  const FaMinusH = new Matrix(n, n);
  const FbMinusH = new Matrix(n, n);
  for (let k = 0; k < n * n; k++) {
    FaMinusH.data[k] = fockAlpha.data[k] - coreHamiltonian.data[k];
    FbMinusH.data[k] = fockBeta.data[k] - coreHamiltonian.data[k];
  }
  const twoElectron = 0.5 * (traceProduct(densityAlpha, FaMinusH)
                            + traceProduct(densityBeta, FbMinusH));

  return { oneElectron, twoElectron, kinetic, nuclearAttraction };
}
