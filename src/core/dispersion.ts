/**
 * DFT-D2 dispersion correction (Grimme 2006).
 *
 * Reference:
 *   S. Grimme, J. Comput. Chem. 27, 1787 (2006). DOI 10.1002/jcc.20495.
 *
 * Energy formula:
 *   E_disp = -s6 Σ_{A<B} C6_AB / R_AB^6 · f_damp(R_AB)
 *   C6_AB  = sqrt(C6_A · C6_B)
 *   f_damp = 1 / (1 + exp(-d · (R_AB / R_r - 1)))     (Fermi-type damping)
 *   R_r    = R0_A + R0_B                              (vdW radii sum)
 *
 * Parameters (Z = 1–36) below are taken from Table 1 of the original paper.
 * D2 is the simplest of the Grimme family; D3 / D3BJ improve accuracy via
 * coordination-number-dependent C6 but are not implemented here yet.
 */

import type { Atom } from './types';

/** Conversion factor: 1 J·nm⁶/mol = 17.34528 Eh·Bohr⁶ (atomic units). */
const C6_JNM6MOL_TO_AU = 17.34528;

/** Bohr / Angstrom conversion (Å → Bohr). */
const ANG_TO_BOHR = 1.8897259886;

/** Damping steepness coefficient (Grimme 2006). */
const D2_DAMP_D = 20.0;

/** Per-element parameters: C6 in J·nm⁶/mol, R0 in Å.
 *  Z = 1 (H) … 36 (Kr). Index = atomic number (0 unused). */
const D2_C6_JNM6MOL: number[] = [
  0,
  /* 1  H  */ 0.14,  /* 2  He */ 0.08,
  /* 3  Li */ 1.61,  /* 4  Be */ 1.61,
  /* 5  B  */ 3.13,  /* 6  C  */ 1.75,  /* 7  N  */ 1.23,
  /* 8  O  */ 0.70,  /* 9  F  */ 0.75,  /* 10 Ne */ 0.63,
  /* 11 Na */ 5.71,  /* 12 Mg */ 5.71,
  /* 13 Al */ 10.79, /* 14 Si */ 9.23,  /* 15 P  */ 7.84,
  /* 16 S  */ 5.57,  /* 17 Cl */ 5.07,  /* 18 Ar */ 4.61,
  /* 19 K  */ 10.80, /* 20 Ca */ 10.80,
  /* 21 Sc */ 10.80, /* 22 Ti */ 10.80, /* 23 V  */ 10.80,
  /* 24 Cr */ 10.80, /* 25 Mn */ 10.80, /* 26 Fe */ 10.80,
  /* 27 Co */ 10.80, /* 28 Ni */ 10.80, /* 29 Cu */ 10.80, /* 30 Zn */ 10.80,
  /* 31 Ga */ 16.99, /* 32 Ge */ 17.10, /* 33 As */ 16.37,
  /* 34 Se */ 12.64, /* 35 Br */ 12.47, /* 36 Kr */ 12.01,
];

const D2_R0_ANG: number[] = [
  0,
  /* 1  H  */ 1.001, /* 2  He */ 1.012,
  /* 3  Li */ 0.825, /* 4  Be */ 1.408,
  /* 5  B  */ 1.485, /* 6  C  */ 1.452, /* 7  N  */ 1.397,
  /* 8  O  */ 1.342, /* 9  F  */ 1.287, /* 10 Ne */ 1.243,
  /* 11 Na */ 1.144, /* 12 Mg */ 1.364,
  /* 13 Al */ 1.639, /* 14 Si */ 1.716, /* 15 P  */ 1.705,
  /* 16 S  */ 1.683, /* 17 Cl */ 1.639, /* 18 Ar */ 1.595,
  /* 19 K  */ 1.485, /* 20 Ca */ 1.474,
  /* 21 Sc */ 1.562, /* 22 Ti */ 1.562, /* 23 V  */ 1.562,
  /* 24 Cr */ 1.562, /* 25 Mn */ 1.562, /* 26 Fe */ 1.562,
  /* 27 Co */ 1.562, /* 28 Ni */ 1.562, /* 29 Cu */ 1.562, /* 30 Zn */ 1.562,
  /* 31 Ga */ 1.650, /* 32 Ge */ 1.727, /* 33 As */ 1.760,
  /* 34 Se */ 1.771, /* 35 Br */ 1.749, /* 36 Kr */ 1.727,
];

/** s6 scaling factor per functional (Grimme 2006).
 *  Functionals not listed default to 1.0 (acceptable approximation). */
const D2_S6: Record<string, number> = {
  'BLYP':  1.20,
  'BP86':  1.05,
  'B3LYP': 1.05,
  'PBE':   0.75,
  'PBE0':  0.60,
  'TPSS':  1.00,
  // No official D2 value for SVWN (LDA) — fall back to 1.0 default.
};

export interface D2Result {
  energy: number;       // Total dispersion energy in Hartree
  pairs: number;        // Number of atom pairs included
  s6: number;           // Scaling factor used
}

/** Compute Grimme D2 dispersion correction for a molecule.
 *  @param atoms atoms with .atomicNumber and .coordinate (in Bohr)
 *  @param functional name used to pick s6 (case-insensitive)
 *  @returns dispersion energy in Hartree (negative for attractive interactions)
 */
export function computeD2Dispersion(
  atoms: Atom[],
  functional?: string,
): D2Result {
  const key = (functional ?? '').toUpperCase();
  const s6 = D2_S6[key] ?? 1.0;

  let energy = 0;
  let pairs = 0;
  for (let a = 0; a < atoms.length; a++) {
    const Za = atoms[a].atomicNumber;
    if (Za < 1 || Za > 36) continue;
    const C6a = D2_C6_JNM6MOL[Za];
    const R0a = D2_R0_ANG[Za];
    if (C6a === 0) continue;
    for (let b = a + 1; b < atoms.length; b++) {
      const Zb = atoms[b].atomicNumber;
      if (Zb < 1 || Zb > 36) continue;
      const C6b = D2_C6_JNM6MOL[Zb];
      const R0b = D2_R0_ANG[Zb];
      if (C6b === 0) continue;

      const dx = atoms[a].coordinate.x - atoms[b].coordinate.x;
      const dy = atoms[a].coordinate.y - atoms[b].coordinate.y;
      const dz = atoms[a].coordinate.z - atoms[b].coordinate.z;
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz); // Bohr (coords are in Bohr)
      if (r < 1e-6) continue;

      // C6 combination rule: geometric mean
      const C6ab_jnm6 = Math.sqrt(C6a * C6b);
      const C6ab = C6ab_jnm6 * C6_JNM6MOL_TO_AU; // Eh·Bohr⁶

      // R0 combination: sum of vdW radii (in Bohr for damping comparison)
      const Rr = (R0a + R0b) * ANG_TO_BOHR;

      const f_damp = 1.0 / (1.0 + Math.exp(-D2_DAMP_D * (r / Rr - 1.0)));
      energy -= s6 * f_damp * C6ab / Math.pow(r, 6);
      pairs++;
    }
  }
  return { energy, pairs, s6 };
}

/** Whether the functional has an official D2 s6 parameter. */
export function isD2Parameterized(functional: string): boolean {
  return D2_S6[functional.toUpperCase()] !== undefined;
}

// ─────────────────────────────────────────────────────────────────────
// D3(BJ)-lite — Becke-Johnson rational damping with C6+C8 terms.
//
// Compared to full Grimme D3, this implementation:
//   - keeps the D2 C6 table (no coordination-number-dependent C6 interpolation)
//   - uses the BJ damping form 1/(r^n + (a1·R_BJ + a2)^n)
//   - includes the C8 term with C8_AB = 3·sqrt(Q_A·Q_B)·C6_AB
//
// This captures the *key new physics of D3-BJ* (rational damping, no overshoot,
// C8 contribution) while staying compact. For real production work users should
// use a code with the full ~32k-entry D3 reference table.
//
// Reference for BJ damping:
//   Grimme, Ehrlich, Goerigk, J. Comput. Chem. 32, 1456 (2011).
// ─────────────────────────────────────────────────────────────────────

/** Q_A values from BJ theory: Q_A = sqrt(Z_A) · <r²>_A / <r⁴>_A^(1/2).
 *  Approximate values per element (Z = 1–18). Higher Z falls back to ~5. */
const D3_Q: number[] = [
  0,
  /* H  */ 0.45,
  /* He */ 0.92,
  /* Li */ 4.27,  /* Be */ 4.21,
  /* B  */ 4.07,  /* C  */ 3.40,  /* N  */ 2.78,
  /* O  */ 2.18,  /* F  */ 1.79,  /* Ne */ 1.49,
  /* Na */ 7.26,  /* Mg */ 8.27,
  /* Al */ 8.94,  /* Si */ 8.62,  /* P  */ 7.80,
  /* S  */ 6.94,  /* Cl */ 6.16,  /* Ar */ 5.51,
];

/** Functional-specific BJ parameters (s6, s8, a1, a2).
 *  From Grimme et al., 2011. Common values for popular functionals. */
const D3BJ_PARAMS: Record<string, { s6: number; s8: number; a1: number; a2: number }> = {
  'B3LYP': { s6: 1.000, s8: 1.9889, a1: 0.3981, a2: 4.4211 },
  'BLYP':  { s6: 1.000, s8: 2.6996, a1: 0.4298, a2: 4.2359 },
  'PBE':   { s6: 1.000, s8: 0.7875, a1: 0.4289, a2: 4.4407 },
  'PBE0':  { s6: 1.000, s8: 1.2177, a1: 0.4145, a2: 4.8593 },
  'TPSS':  { s6: 1.000, s8: 1.9435, a1: 0.4535, a2: 4.4752 },
};

export interface D3BJResult {
  energy: number;
  pairs: number;
  s6: number;
  s8: number;
  a1: number;
  a2: number;
}

/** Compute D3(BJ)-lite dispersion correction.
 *  Uses D2's C6 table but BJ rational damping plus C8 term. */
export function computeD3BJDispersion(
  atoms: Atom[],
  functional?: string,
): D3BJResult {
  const key = (functional ?? '').toUpperCase();
  const p = D3BJ_PARAMS[key] ?? { s6: 1.0, s8: 1.0, a1: 0.4, a2: 5.0 };

  let energy = 0;
  let pairs = 0;
  for (let a = 0; a < atoms.length; a++) {
    const Za = atoms[a].atomicNumber;
    if (Za < 1 || Za >= D2_C6_JNM6MOL.length) continue;
    const C6a = D2_C6_JNM6MOL[Za];
    if (C6a === 0) continue;
    const Qa = D3_Q[Za] ?? 5.0;
    for (let b = a + 1; b < atoms.length; b++) {
      const Zb = atoms[b].atomicNumber;
      if (Zb < 1 || Zb >= D2_C6_JNM6MOL.length) continue;
      const C6b = D2_C6_JNM6MOL[Zb];
      if (C6b === 0) continue;
      const Qb = D3_Q[Zb] ?? 5.0;

      const dx = atoms[a].coordinate.x - atoms[b].coordinate.x;
      const dy = atoms[a].coordinate.y - atoms[b].coordinate.y;
      const dz = atoms[a].coordinate.z - atoms[b].coordinate.z;
      const r2 = dx * dx + dy * dy + dz * dz;
      if (r2 < 1e-12) continue;
      const r = Math.sqrt(r2);

      // C6 in atomic units (Eh·Bohr⁶)
      const C6ab = Math.sqrt(C6a * C6b) * C6_JNM6MOL_TO_AU;
      // C8 = 3·sqrt(Q_A·Q_B)·C6_AB (BJ recipe; coefficient 3 for s/p atoms)
      const C8ab = 3.0 * Math.sqrt(Qa * Qb) * C6ab;
      // BJ critical radius
      const R_BJ = Math.sqrt(C8ab / C6ab);
      const damp = p.a1 * R_BJ + p.a2;
      const damp6 = Math.pow(damp, 6);
      const damp8 = Math.pow(damp, 8);
      const r6 = r2 * r2 * r2;
      const r8 = r6 * r2;

      energy -= p.s6 * C6ab / (r6 + damp6);
      energy -= p.s8 * C8ab / (r8 + damp8);
      pairs++;
    }
  }
  return { energy, pairs, ...p };
}

/** Whether the functional has an official D3-BJ parameter set. */
export function isD3BJParameterized(functional: string): boolean {
  return D3BJ_PARAMS[functional.toUpperCase()] !== undefined;
}

export type DispersionMethod = 'none' | 'd2' | 'd3bj';

export interface DispersionResult {
  method: DispersionMethod;
  energy: number;
  s6?: number;
  s8?: number;
  a1?: number;
  a2?: number;
  pairs: number;
}

/** Top-level dispatcher: compute dispersion correction by method name. */
export function computeDispersion(
  atoms: Atom[],
  method: DispersionMethod,
  functional?: string,
): DispersionResult {
  if (method === 'd2') {
    const r = computeD2Dispersion(atoms, functional);
    return { method: 'd2', energy: r.energy, s6: r.s6, pairs: r.pairs };
  } else if (method === 'd3bj') {
    const r = computeD3BJDispersion(atoms, functional);
    return { method: 'd3bj', energy: r.energy, s6: r.s6, s8: r.s8, a1: r.a1, a2: r.a2, pairs: r.pairs };
  }
  return { method: 'none', energy: 0, pairs: 0 };
}
