/** Resolution of Identity (RI / density fitting) approximation.
 *
 *  (μν|λσ) ≈ Σ_P B_P^{μν} B_P^{λσ}
 *  where B_P^{μν} = Σ_Q (P|Q)^{-1/2} (Q|μν)
 */

import type { Atom, PrimitiveShell } from './types';
import { ANGULAR_MOMENTUMS, shellTypeToNumBasis, shellNameToType } from './constants';
import { primitiveNorm } from './integrals1e';
import { boysAll } from './boys';
import { osVRR, hrrEval } from './integrals2e';
import { Matrix } from '../linalg/matrix';
// jacobiEigenAsync no longer needed — using Cholesky decomposition instead
import { BasisSet, ContractedGauss, ElementBasisSet } from './basisSet';
import { atomicNumberToElementName } from './constants';
import {
  isRIWasmAvailable,
  computeRISetupWasm,
  computeRIFockRhfWasm,
  computeRIFockUhfWasm,
  packShells,
} from './eriWasm';

// ---------------------------------------------------------------------------
// Auto auxiliary basis generation (product basis)
// ---------------------------------------------------------------------------

/**
 * Generate an auxiliary basis set automatically from a primary basis set.
 * For each element, pairs of primitive exponents (α_i, l_a) and (α_j, l_b)
 * produce auxiliary functions with exponent α_i + α_j and angular momenta
 * |l_a - l_b| to l_a + l_b.  Each auxiliary function is uncontracted (single
 * primitive with coefficient 1.0).  Near-duplicate exponents within the same
 * angular momentum are merged.
 */
export function generateAutoAuxBasis(primaryBasis: BasisSet, elementNames: string[]): BasisSet {
  const SHELL_NAMES = ['S', 'P', 'D', 'F', 'G', 'H', 'I'];
  const auxBasis = new BasisSet();

  for (const elem of new Set(elementNames)) {
    const ebs = primaryBasis.get(elem);

    // Collect all (exponent, shellType) pairs across primitives
    const primInfo: { exp: number; l: number }[] = [];
    for (const cg of ebs.contractedGausses) {
      const l = shellNameToType(cg.type);
      for (const p of cg.primitives) {
        primInfo.push({ exp: p.exponent, l });
      }
    }

    // Generate product exponents for each target angular momentum
    // Cap aux angular momentum at primary l_max + 1
    const primaryLmax = primInfo.reduce((mx, p) => Math.max(mx, p.l), 0);
    const auxLmax = Math.min(primaryLmax + 1, SHELL_NAMES.length - 1);
    const auxExponents = new Map<number, number[]>();

    for (let i = 0; i < primInfo.length; i++) {
      for (let j = i; j < primInfo.length; j++) {
        const gamma = primInfo[i].exp + primInfo[j].exp;
        const lMin = Math.abs(primInfo[i].l - primInfo[j].l);
        // +1 accounts for angular momentum from displacement between centers
        const lMax = Math.min(primInfo[i].l + primInfo[j].l + 1, auxLmax);
        for (let l = lMin; l <= lMax; l++) {
          let arr = auxExponents.get(l);
          if (!arr) { arr = []; auxExponents.set(l, arr); }
          arr.push(gamma);
        }
      }
    }

    // Deduplicate: keep exponents spaced by at least factor of 3, cap per l
    const MAX_PER_L = 5;
    const auxEbs = new ElementBasisSet();
    auxEbs.elementName = elem;

    for (const [l, exps] of auxExponents) {
      exps.sort((a, b) => b - a); // descending
      const unique: number[] = [];
      for (const e of exps) {
        if (unique.length === 0 || e < unique[unique.length - 1] * 0.33) {
          unique.push(e);
        }
      }
      // If still too many, subsample evenly in log-space
      let selected = unique;
      if (selected.length > MAX_PER_L) {
        const step = (selected.length - 1) / (MAX_PER_L - 1);
        selected = [];
        for (let k = 0; k < MAX_PER_L; k++) {
          selected.push(unique[Math.round(k * step)]);
        }
      }
      for (const exp of selected) {
        const cg = new ContractedGauss(SHELL_NAMES[l]);
        cg.addPrimitive(exp, 1.0);
        auxEbs.addContractedGauss(cg);
      }
    }

    auxBasis.add(auxEbs);
  }

  return auxBasis;
}

// ---------------------------------------------------------------------------
// Auxiliary basis shell construction
// ---------------------------------------------------------------------------

interface AuxShellGroup {
  primitives: PrimitiveShell[];
  basisIndex: number;
  shellType: number;
}

export function buildAuxShells(
  atoms: Atom[], auxBasis: BasisSet,
): { shells: PrimitiveShell[]; normFactors: number[]; naux: number } {
  const shells: PrimitiveShell[] = [];
  const normFactors: number[] = [];
  let basisIndex = 0;

  for (let ai = 0; ai < atoms.length; ai++) {
    const atom = atoms[ai];
    const elementName = atomicNumberToElementName(atom.atomicNumber);
    const ebs = auxBasis.get(elementName);

    for (const cg of ebs.contractedGausses) {
      const st = shellNameToType(cg.type);
      for (const prim of cg.primitives) {
        shells.push({
          exponent: prim.exponent,
          coefficient: prim.coefficient,
          coordinate: { ...atom.coordinate },
          shellType: st,
          basisIndex,
          atomIndex: ai,
        });
      }
      normFactors.push(...cg.getNormalizationFactors());
      basisIndex += shellTypeToNumBasis(st);
    }
  }

  // Sort by shellType then basisIndex (same as Molecular)
  shells.sort((a, b) =>
    a.shellType !== b.shellType ? a.shellType - b.shellType : a.basisIndex - b.basisIndex
  );

  return { shells, normFactors, naux: basisIndex };
}

function groupShells(shells: PrimitiveShell[]): AuxShellGroup[] {
  const groups: AuxShellGroup[] = [];
  let current: AuxShellGroup | null = null;
  for (const ps of shells) {
    if (!current || current.basisIndex !== ps.basisIndex || current.shellType !== ps.shellType) {
      current = { primitives: [], basisIndex: ps.basisIndex, shellType: ps.shellType };
      groups.push(current);
    }
    current.primitives.push(ps);
  }
  return groups;
}

// ---------------------------------------------------------------------------
// 2-center Coulomb integrals: (P|Q) for auxiliary basis
// ---------------------------------------------------------------------------

async function compute2CenterIntegrals(
  auxGroups: AuxShellGroup[],
  auxNorm: number[],
  naux: number,
): Promise<Matrix> {
  const V = new Matrix(naux, naux);

  for (let iGrp = 0; iGrp < auxGroups.length; iGrp++) {
    for (let jGrp = iGrp; jGrp < auxGroups.length; jGrp++) {
      compute2CShellPair(auxGroups[iGrp], auxGroups[jGrp], auxNorm, V, iGrp === jGrp);
    }
    // Yield after each outer shell group
    await yieldToUI();
  }

  return V;
}

function compute2CShellPair(
  grpA: AuxShellGroup, grpB: AuxShellGroup,
  normFactors: number[], V: Matrix, sameGroup: boolean,
) {
  // 2-center (P|Q) treated as (P,0|Q,0) — single uncontracted functions
  const la = grpA.shellType, lb = grpB.shellType;
  const angA = ANGULAR_MOMENTUMS[la], angB = ANGULAR_MOMENTUMS[lb];
  const numA = shellTypeToNumBasis(la), numB = shellTypeToNumBasis(lb);

  // Precompute primitive norms
  const pNormsA: Float64Array[] = [];
  for (let ia = 0; ia < numA; ia++) {
    const [lx, ly, lz] = angA[ia];
    const norms = new Float64Array(grpA.primitives.length);
    for (let pp = 0; pp < grpA.primitives.length; pp++)
      norms[pp] = primitiveNorm(grpA.primitives[pp].exponent, lx, ly, lz);
    pNormsA.push(norms);
  }
  const pNormsB: Float64Array[] = [];
  for (let ib = 0; ib < numB; ib++) {
    const [lx, ly, lz] = angB[ib];
    const norms = new Float64Array(grpB.primitives.length);
    for (let pp = 0; pp < grpB.primitives.length; pp++)
      norms[pp] = primitiveNorm(grpB.primitives[pp].exponent, lx, ly, lz);
    pNormsB.push(norms);
  }

  // Accumulate contracted values
  const contracted = new Float64Array(numA * numB);

  for (let pA = 0; pA < grpA.primitives.length; pA++) {
    const primA = grpA.primitives[pA];
    const alpha = primA.exponent;
    const Ax = primA.coordinate.x, Ay = primA.coordinate.y, Az = primA.coordinate.z;
    // Single function: p = alpha, P = A
    const p = alpha;

    for (let pB = 0; pB < grpB.primitives.length; pB++) {
      const primB = grpB.primitives[pB];
      const beta = primB.exponent;
      const Bx = primB.coordinate.x, By = primB.coordinate.y, Bz = primB.coordinate.z;
      const q = beta;

      const eta = p * q / (p + q);
      const PQx = Ax - Bx, PQy = Ay - By, PQz = Az - Bz;
      const T = eta * (PQx * PQx + PQy * PQy + PQz * PQz);
      const boys = boysAll(la + lb, T);
      const prefactor = 2 * Math.pow(Math.PI, 2.5) / (p * q * Math.sqrt(p + q));

      // W = (p*A + q*B) / (p+q)
      const Wx = (p * Ax + q * Bx) / (p + q);
      const Wy = (p * Ay + q * By) / (p + q);
      const Wz = (p * Az + q * Bz) / (p + q);

      // VRR: (P,0|Q,0) → la on bra, lb on ket, no HRR needed (lb_shell = ld_shell = 0)
      const { V: vrrV, vi } = osVRR(la, lb, p, q,
        0, 0, 0,  // PA = P - A = 0 (single function on A)
        Wx - Ax, Wy - Ay, Wz - Az,
        0, 0, 0,  // QC = Q - B = 0 (single function on B)
        Wx - Bx, Wy - By, Wz - Bz,
        boys);

      // No HRR needed (lb = ld = 0), read VRR directly
      for (let ia = 0; ia < numA; ia++) {
        const [lxa, lya, lza] = angA[ia];
        const cA = primA.coefficient * pNormsA[ia][pA];
        for (let ib = 0; ib < numB; ib++) {
          const [lxb, lyb, lzb] = angB[ib];
          const cB = primB.coefficient * pNormsB[ib][pB];
          contracted[ia * numB + ib] += prefactor * cA * cB
            * vrrV[vi(lxa, lya, lza, lxb, lyb, lzb, 0)];
        }
      }
    }
  }

  // Store into matrix
  for (let ia = 0; ia < numA; ia++) {
    const muA = grpA.basisIndex + ia;
    const ibStart = sameGroup ? ia : 0;
    for (let ib = ibStart; ib < numB; ib++) {
      const muB = grpB.basisIndex + ib;
      const val = contracted[ia * numB + ib] * normFactors[muA] * normFactors[muB];
      V.set(muA, muB, V.get(muA, muB) + val);
      if (muA !== muB) V.set(muB, muA, V.get(muA, muB));
    }
  }
}

// ---------------------------------------------------------------------------
// 3-center Coulomb integrals: (μν|P)
// Stored as three[P * N * N + mu * N + nu], symmetric in μν
// ---------------------------------------------------------------------------

function compute3CenterIntegrals(
  primaryGroups: AuxShellGroup[],
  primaryNorm: number[],
  nbasis: number,
  auxGroups: AuxShellGroup[],
  auxNorm: number[],
  naux: number,
): Float64Array {
  const three = new Float64Array(naux * nbasis * nbasis);

  for (let iGrp = 0; iGrp < primaryGroups.length; iGrp++) {
    for (let jGrp = iGrp; jGrp < primaryGroups.length; jGrp++) {
      for (let pGrp = 0; pGrp < auxGroups.length; pGrp++) {
        compute3CShellTriple(
          primaryGroups[iGrp], primaryGroups[jGrp], auxGroups[pGrp],
          primaryNorm, auxNorm, nbasis, three, iGrp === jGrp,
        );
      }
    }
  }

  return three;
}

function compute3CShellTriple(
  grpA: AuxShellGroup, grpB: AuxShellGroup, grpP: AuxShellGroup,
  primaryNorm: number[], auxNorm: number[],
  nbasis: number, three: Float64Array, sameAB: boolean,
) {
  // 3-center (μν|P) treated as (μν|P,0) — bra is pair (A,B), ket is single P
  const la = grpA.shellType, lb = grpB.shellType, lp = grpP.shellType;
  const angA = ANGULAR_MOMENTUMS[la], angB = ANGULAR_MOMENTUMS[lb], angP = ANGULAR_MOMENTUMS[lp];
  const numA = shellTypeToNumBasis(la), numB = shellTypeToNumBasis(lb), numP = shellTypeToNumBasis(lp);
  const NN = nbasis * nbasis;
  const LAB = la + lb;

  // Precompute primitive norms
  const pNormsA: Float64Array[] = [];
  for (let ia = 0; ia < numA; ia++) {
    const [lx, ly, lz] = angA[ia];
    const norms = new Float64Array(grpA.primitives.length);
    for (let pp = 0; pp < grpA.primitives.length; pp++)
      norms[pp] = primitiveNorm(grpA.primitives[pp].exponent, lx, ly, lz);
    pNormsA.push(norms);
  }
  const pNormsB: Float64Array[] = [];
  for (let ib = 0; ib < numB; ib++) {
    const [lx, ly, lz] = angB[ib];
    const norms = new Float64Array(grpB.primitives.length);
    for (let pp = 0; pp < grpB.primitives.length; pp++)
      norms[pp] = primitiveNorm(grpB.primitives[pp].exponent, lx, ly, lz);
    pNormsB.push(norms);
  }
  const pNormsP: Float64Array[] = [];
  for (let ip = 0; ip < numP; ip++) {
    const [lx, ly, lz] = angP[ip];
    const norms = new Float64Array(grpP.primitives.length);
    for (let pp = 0; pp < grpP.primitives.length; pp++)
      norms[pp] = primitiveNorm(grpP.primitives[pp].exponent, lx, ly, lz);
    pNormsP.push(norms);
  }

  // Accumulated contracted integrals [ia * numB * numP + ib * numP + ip]
  const contracted = new Float64Array(numA * numB * numP);

  for (let pA = 0; pA < grpA.primitives.length; pA++) {
    const primA = grpA.primitives[pA];
    const alpha = primA.exponent;
    const Ax = primA.coordinate.x, Ay = primA.coordinate.y, Az = primA.coordinate.z;

    for (let pB = 0; pB < grpB.primitives.length; pB++) {
      const primB = grpB.primitives[pB];
      const beta = primB.exponent;
      const Bx = primB.coordinate.x, By = primB.coordinate.y, Bz = primB.coordinate.z;
      const p = alpha + beta;
      const muAB = alpha * beta / p;
      const AB2 = (Ax - Bx) ** 2 + (Ay - By) ** 2 + (Az - Bz) ** 2;
      const Kab = Math.exp(-muAB * AB2);
      const Px = (alpha * Ax + beta * Bx) / p;
      const Py = (alpha * Ay + beta * By) / p;
      const Pz = (alpha * Az + beta * Bz) / p;

      for (let pP = 0; pP < grpP.primitives.length; pP++) {
        const primP = grpP.primitives[pP];
        const gamma = primP.exponent;
        const q = gamma; // single function
        const Cx = primP.coordinate.x, Cy = primP.coordinate.y, Cz = primP.coordinate.z;

        const eta = p * q / (p + q);
        const PQx = Px - Cx, PQy = Py - Cy, PQz = Pz - Cz;
        const T = eta * (PQx * PQx + PQy * PQy + PQz * PQz);
        const boys = boysAll(LAB + lp, T);
        const prefactor = 2 * Math.pow(Math.PI, 2.5) / (p * q * Math.sqrt(p + q)) * Kab;

        const Wx = (p * Px + q * Cx) / (p + q);
        const Wy = (p * Py + q * Cy) / (p + q);
        const Wz = (p * Pz + q * Cz) / (p + q);

        // VRR: bra LAB = la+lb, ket LCD = lp (D shell = s, ld=0)
        const { V: vrrV, vi } = osVRR(LAB, lp, p, q,
          Px - Ax, Py - Ay, Pz - Az,
          Wx - Px, Wy - Py, Wz - Pz,
          0, 0, 0,  // QC = Q - C = 0 (single function on C)
          Wx - Cx, Wy - Cy, Wz - Cz,
          boys);

        // HRR: transfer lb from A to B (ket has ld=0, no ket HRR needed)
        const ABx = Ax - Bx, ABy = Ay - By, ABz = Az - Bz;
        const memo = new Map<number, number>();

        for (let ia = 0; ia < numA; ia++) {
          const [lxa, lya, lza] = angA[ia];
          const cA = primA.coefficient * pNormsA[ia][pA];
          for (let ib = 0; ib < numB; ib++) {
            const [lxb, lyb, lzb] = angB[ib];
            const cB = primB.coefficient * pNormsB[ib][pB];
            for (let ip = 0; ip < numP; ip++) {
              const [lxp, lyp, lzp] = angP[ip];
              const cP = primP.coefficient * pNormsP[ip][pP];
              // HRR for bra: (lxa,lya,lza, lxb,lyb,lzb | lxp,lyp,lzp, 0,0,0)
              const val = hrrEval(
                lxa, lya, lza, lxb, lyb, lzb,
                lxp, lyp, lzp, 0, 0, 0,
                vrrV, vi, ABx, ABy, ABz, 0, 0, 0, memo);
              contracted[(ia * numB + ib) * numP + ip] += prefactor * cA * cB * cP * val;
            }
          }
        }
      }
    }
  }

  // Store into three-center array
  for (let ia = 0; ia < numA; ia++) {
    const muA = grpA.basisIndex + ia;
    const ibStart = sameAB ? ia : 0;
    for (let ib = ibStart; ib < numB; ib++) {
      const muB = grpB.basisIndex + ib;
      for (let ip = 0; ip < numP; ip++) {
        const muP = grpP.basisIndex + ip;
        const val = contracted[(ia * numB + ib) * numP + ip]
          * primaryNorm[muA] * primaryNorm[muB] * auxNorm[muP];
        three[muP * NN + muA * nbasis + muB] += val;
        if (muA !== muB) {
          three[muP * NN + muB * nbasis + muA] += val;
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// RIData: precomputed B matrix and Fock build methods
// ---------------------------------------------------------------------------

/** Yield to browser event loop to prevent timeout dialogs */
function yieldToUI(): Promise<void> {
  return new Promise<void>(r => setTimeout(r, 0));
}

export class RIData {
  /** B[P * N * N + mu * N + nu] = fitting coefficient */
  readonly bMatrix: Float64Array;
  readonly naux: number;
  readonly nbasis: number;

  private constructor(bMatrix: Float64Array, naux: number, nbasis: number) {
    this.bMatrix = bMatrix;
    this.naux = naux;
    this.nbasis = nbasis;
  }

  static async build(
    primaryShells: PrimitiveShell[],
    primaryNorm: number[],
    nbasis: number,
    atoms: Atom[],
    auxBasis: BasisSet,
    onProgress?: (msg: string) => void,
    options?: { forceJS?: boolean },
  ): Promise<RIData> {
    const { shells: auxShells, normFactors: auxNorm, naux } = buildAuxShells(atoms, auxBasis);

    // Try WASM path
    if (!options?.forceJS && isRIWasmAvailable()) {
      onProgress?.('Computing RI setup via WASM...');
      const primaryShellsFlat = packShells(primaryShells);
      const auxShellsFlat = packShells(auxShells);
      const result = computeRISetupWasm(
        primaryShellsFlat, new Float64Array(primaryNorm), nbasis,
        auxShellsFlat, new Float64Array(auxNorm), naux,
      );
      if (result) {
        return new RIData(result.bMatrix, result.naux, nbasis);
      }
    }

    // JS fallback
    const primaryGroups = groupShells(primaryShells);
    const auxGroups = groupShells(auxShells);

    // 1. (P|Q) 2-center integrals
    onProgress?.('Computing 2-center integrals (P|Q)...');
    const Vpq = await compute2CenterIntegrals(auxGroups, auxNorm, naux);
    await yieldToUI();

    // 2. Cholesky decomposition: V = L * L^T
    //    Drop near-linearly-dependent aux functions (small diagonal → zero row)
    onProgress?.('Cholesky decomposition of (P|Q)...');
    const CHOL_THRESH = 1e-8;
    const L = new Float64Array(naux * naux);
    const auxValid = new Uint8Array(naux); // 1 = kept, 0 = dropped
    let nDropped = 0;
    for (let i = 0; i < naux; i++) {
      // Off-diagonal: L[i,j] for j < i
      for (let j = 0; j < i; j++) {
        if (!auxValid[j]) continue; // skip dropped columns
        let sum = Vpq.get(i, j);
        for (let k = 0; k < j; k++) {
          if (!auxValid[k]) continue;
          sum -= L[i * naux + k] * L[j * naux + k];
        }
        L[i * naux + j] = sum / L[j * naux + j];
      }
      // Diagonal: L[i,i]
      let diag = Vpq.get(i, i);
      for (let k = 0; k < i; k++) {
        if (!auxValid[k]) continue;
        diag -= L[i * naux + k] * L[i * naux + k];
      }
      if (diag < CHOL_THRESH) {
        // Near-singular: drop this auxiliary function
        L[i * naux + i] = 0;
        // Zero out off-diagonal elements of this row
        for (let j = 0; j < i; j++) L[i * naux + j] = 0;
        nDropped++;
      } else {
        L[i * naux + i] = Math.sqrt(diag);
        auxValid[i] = 1;
      }
    }
    if (nDropped > 0) {
      onProgress?.(`Dropped ${nDropped}/${naux} near-dependent aux functions`);
    }
    await yieldToUI();

    // 3. 3-center integrals (μν|Q) — yield between shell groups
    onProgress?.('Computing 3-center integrals (μν|P)...');
    const three = new Float64Array(naux * nbasis * nbasis);
    let shellCount = 0;
    const totalShellPairs = primaryGroups.length * (primaryGroups.length + 1) / 2 * auxGroups.length;
    for (let iGrp = 0; iGrp < primaryGroups.length; iGrp++) {
      for (let jGrp = iGrp; jGrp < primaryGroups.length; jGrp++) {
        for (let pGrp = 0; pGrp < auxGroups.length; pGrp++) {
          compute3CShellTriple(
            primaryGroups[iGrp], primaryGroups[jGrp], auxGroups[pGrp],
            primaryNorm, auxNorm, nbasis, three, iGrp === jGrp,
          );
          shellCount++;
        }
        // Yield every outer pair to keep browser responsive
        const pct = Math.round(100 * shellCount / totalShellPairs);
        onProgress?.(`Computing 3-center integrals... ${pct}%`);
        await yieldToUI();
      }
    }
    await yieldToUI();

    // 4. B = L^{-1} * three  (forward substitution for each μν pair)
    //    B_P^{μν} = [L^{-1} * three]_{P,μν}
    //    so that Σ_P B_P^{μν} B_P^{λσ} = threeᵀ V^{-1} three
    onProgress?.('Building B matrix (L^{-1} × three)...');
    const NN = nbasis * nbasis;
    const B = new Float64Array(naux * NN);
    // Forward substitution: for each column k of three, solve L * x = three[:,k]
    // Skip dropped auxiliary functions (auxValid[i] === 0)
    for (let k = 0; k < NN; k++) {
      for (let i = 0; i < naux; i++) {
        if (!auxValid[i]) continue; // dropped — B row stays zero
        let sum = three[i * NN + k];
        for (let j = 0; j < i; j++) {
          if (!auxValid[j]) continue;
          sum -= L[i * naux + j] * B[j * NN + k];
        }
        B[i * NN + k] = sum / L[i * naux + i];
      }
    }
    await yieldToUI();

    return new RIData(B, naux, nbasis);
  }

  /** Coulomb matrix: J[μν] = Σ_P B_P^{μν} * d_P
   *  where d_P = Σ_λσ P[λσ] * B_P^{λσ} */
  buildCoulomb(density: Matrix): Matrix {
    const N = this.nbasis, naux = this.naux, NN = N * N;
    const B = this.bMatrix;

    // d_P = Σ_λσ P[λσ] * B_P^{λσ}
    const d = new Float64Array(naux);
    for (let P = 0; P < naux; P++) {
      let sum = 0;
      const off = P * NN;
      for (let k = 0; k < NN; k++) {
        sum += density.data[k] * B[off + k];
      }
      d[P] = sum;
    }

    // J[μν] = Σ_P B_P^{μν} * d_P
    const J = new Matrix(N, N);
    for (let P = 0; P < naux; P++) {
      if (Math.abs(d[P]) < 1e-15) continue;
      const off = P * NN;
      for (let k = 0; k < NN; k++) {
        J.data[k] += B[off + k] * d[P];
      }
    }

    return J;
  }

  /** Check if RI WASM acceleration is available */
  get wasmAvailable(): boolean {
    return isRIWasmAvailable();
  }

  /** Exchange matrix for RHF: K[μν] = Σ_P Σ_i X_P^{μi} X_P^{νi}
   *  where X_P^{μi} = Σ_λ C[λ,i] * B_P^{μλ} */
  buildExchange(coefficients: Matrix, nocc: number): Matrix {
    const N = this.nbasis, naux = this.naux;
    const B = this.bMatrix;

    // X[P][mu][i] = Σ_λ B_P^{μλ} * C[λ,i]
    // Store as X[P * N * nocc + mu * nocc + i]
    const X = new Float64Array(naux * N * nocc);
    for (let P = 0; P < naux; P++) {
      const Boff = P * N * N;
      const Xoff = P * N * nocc;
      for (let mu = 0; mu < N; mu++) {
        for (let i = 0; i < nocc; i++) {
          let sum = 0;
          for (let lam = 0; lam < N; lam++) {
            sum += B[Boff + mu * N + lam] * coefficients.get(lam, i);
          }
          X[Xoff + mu * nocc + i] = sum;
        }
      }
    }

    // K[μν] = Σ_P Σ_i X_P^{μi} * X_P^{νi}
    const K = new Matrix(N, N);
    for (let P = 0; P < naux; P++) {
      const Xoff = P * N * nocc;
      for (let mu = 0; mu < N; mu++) {
        for (let nu = mu; nu < N; nu++) {
          let sum = 0;
          for (let i = 0; i < nocc; i++) {
            sum += X[Xoff + mu * nocc + i] * X[Xoff + nu * nocc + i];
          }
          K.data[mu * N + nu] += sum;
          if (mu !== nu) K.data[nu * N + mu] += sum;
        }
      }
    }

    return K;
  }

  /** Build RHF Fock matrix: F = H + J - K
   *  J uses density P (which includes factor 2 for closed shell).
   *  K = Σ_i X^{μi} X^{νi} uses coefficients directly (no factor 2),
   *  so no 0.5 multiplier is needed (unlike conventional where K_conv = Σ P*(μk|νl) has factor 2). */
  buildFockRHF(density: Matrix, coreH: Matrix, coefficients: Matrix, nocc: number): Matrix {
    const N = this.nbasis;

    // Try WASM path
    if (isRIWasmAvailable()) {
      const result = computeRIFockRhfWasm(
        this.bMatrix, density.data, coreH.data, coefficients.data,
        this.naux, N, nocc,
      );
      if (result) {
        const F = new Matrix(N, N);
        F.data.set(result);
        return F;
      }
    }

    // JS fallback
    const J = this.buildCoulomb(density);
    const K = this.buildExchange(coefficients, nocc);
    const F = new Matrix(N, N);
    for (let k = 0; k < N * N; k++) {
      F.data[k] = coreH.data[k] + J.data[k] - K.data[k];
    }
    return F;
  }

  /** Build UHF/ROHF Fock matrices: Fα = H + J(Ptot) - Kα, Fβ = H + J(Ptot) - Kβ
   *  J uses total density Ptot = Pα + Pβ.
   *  Kα = Σ_i X^{μi}_α X^{νi}_α uses alpha coefficients directly.
   *  Kβ = Σ_i X^{μi}_β X^{νi}_β uses beta coefficients directly. */
  buildFockUHF(
    densityTotal: Matrix, coreH: Matrix,
    coeffAlpha: Matrix, noccAlpha: number,
    coeffBeta: Matrix, noccBeta: number,
  ): { fockAlpha: Matrix; fockBeta: Matrix } {
    const N = this.nbasis;
    const nn = N * N;

    // Try WASM path
    if (isRIWasmAvailable()) {
      const result = computeRIFockUhfWasm(
        this.bMatrix, densityTotal.data, coreH.data,
        coeffAlpha.data, noccAlpha,
        coeffBeta.data, noccBeta,
        this.naux, N,
      );
      if (result) {
        const Fa = new Matrix(N, N);
        const Fb = new Matrix(N, N);
        Fa.data.set(result.fockAlpha);
        Fb.data.set(result.fockBeta);
        return { fockAlpha: Fa, fockBeta: Fb };
      }
    }

    // JS fallback
    const J = this.buildCoulomb(densityTotal);
    const Ka = this.buildExchange(coeffAlpha, noccAlpha);
    const Kb = this.buildExchange(coeffBeta, noccBeta);
    const Fa = new Matrix(N, N);
    const Fb = new Matrix(N, N);
    for (let k = 0; k < nn; k++) {
      const hj = coreH.data[k] + J.data[k];
      Fa.data[k] = hj - Ka.data[k];
      Fb.data[k] = hj - Kb.data[k];
    }
    return { fockAlpha: Fa, fockBeta: Fb };
  }

  /** Get B tensor for MO transformation (used by RI-MP2).
   *  Returns B_ia^P = Σ_μν C[μ,i] C[ν,a] B_P^{μν}
   *  Stored as bmo[i * nvir * naux + a * naux + P] */
  transformToMO(C: Matrix, nocc: number): Float64Array {
    const N = this.nbasis, naux = this.naux;
    const nvir = N - nocc;
    const B = this.bMatrix;

    // Step 1: half-transform B_P^{μa} = Σ_ν B_P^{μν} C[ν, a+nocc]
    // Then: B_ia^P = Σ_μ C[μ,i] * B_P^{μa}
    const bmo = new Float64Array(nocc * nvir * naux);

    // Temporary for half-transform: bHalf[mu * nvir + a] per P
    const bHalf = new Float64Array(N * nvir);

    for (let P = 0; P < naux; P++) {
      const Boff = P * N * N;

      // Half-transform ν→a: bHalf[μ,a] = Σ_ν B_P^{μν} C[ν,a+nocc]
      bHalf.fill(0);
      for (let mu = 0; mu < N; mu++) {
        for (let nu = 0; nu < N; nu++) {
          const bVal = B[Boff + mu * N + nu];
          if (Math.abs(bVal) < 1e-15) continue;
          for (let a = 0; a < nvir; a++) {
            bHalf[mu * nvir + a] += bVal * C.get(nu, nocc + a);
          }
        }
      }

      // Full transform μ→i: bmo[i,a,P] = Σ_μ C[μ,i] * bHalf[μ,a]
      for (let mu = 0; mu < N; mu++) {
        for (let i = 0; i < nocc; i++) {
          const Cmi = C.get(mu, i);
          if (Math.abs(Cmi) < 1e-15) continue;
          for (let a = 0; a < nvir; a++) {
            bmo[i * nvir * naux + a * naux + P] += Cmi * bHalf[mu * nvir + a];
          }
        }
      }
    }

    return bmo;
  }
}
