/** Physical constants, element maps, angular momentums — mirroring GANSU types.hpp + parameters.h */

// Unit conversions
export const ANGSTROM_TO_BOHR = 1.8897259886;
export const BOHR_TO_ANGSTROM = 0.529177210903;

export function angstromToBohr(a: number): number {
  return a * ANGSTROM_TO_BOHR;
}

export function bohrToAngstrom(b: number): number {
  return b * BOHR_TO_ANGSTROM;
}

// Boys function parameters (from parameters.h)
export const A_TR = 0.352905920120321;
export const B_TR = 0.015532762923351;
export const LUT_N_RANGE = 24;
export const LUT_XI_RANGE = 32.0;
export const LUT_XI_INTERVAL = 0.03125;
export const LUT_NUM_XI = 1024;
export const LUT_K_MAX = 5;

// Normalization constants
export const TWO_TIMES_PI_TO_THE_2_POINT_5_TH_POWER = 34.986836655249725;

// Element name -> atomic number
export const ELEMENT_NAME_TO_ATOMIC_NUMBER: Record<string, number> = {
  H: 1, He: 2, Li: 3, Be: 4, B: 5, C: 6, N: 7, O: 8, F: 9, Ne: 10,
  Na: 11, Mg: 12, Al: 13, Si: 14, P: 15, S: 16, Cl: 17, Ar: 18, K: 19, Ca: 20,
  Sc: 21, Ti: 22, V: 23, Cr: 24, Mn: 25, Fe: 26, Co: 27, Ni: 28, Cu: 29, Zn: 30,
  Ga: 31, Ge: 32, As: 33, Se: 34, Br: 35, Kr: 36, Rb: 37, Sr: 38, Y: 39, Zr: 40,
  Nb: 41, Mo: 42, Tc: 43, Ru: 44, Rh: 45, Pd: 46, Ag: 47, Cd: 48, In: 49, Sn: 50,
  Sb: 51, Te: 52, I: 53, Xe: 54, Cs: 55, Ba: 56, La: 57, Ce: 58, Pr: 59, Nd: 60,
  Pm: 61, Sm: 62, Eu: 63, Gd: 64, Tb: 65, Dy: 66, Ho: 67, Er: 68, Tm: 69, Yb: 70,
  Lu: 71, Hf: 72, Ta: 73, W: 74, Re: 75, Os: 76, Ir: 77, Pt: 78, Au: 79, Hg: 80,
  Tl: 81, Pb: 82, Bi: 83, Po: 84, At: 85, Rn: 86, Fr: 87, Ra: 88, Ac: 89, Th: 90,
  Pa: 91, U: 92, Np: 93, Pu: 94, Am: 95, Cm: 96, Bk: 97, Cf: 98, Es: 99, Fm: 100,
  Md: 101, No: 102, Lr: 103, Rf: 104, Db: 105, Sg: 106, Bh: 107, Hs: 108, Mt: 109, Ds: 110,
  Rg: 111, Cn: 112, Nh: 113, Fl: 114, Mc: 115, Lv: 116, Ts: 117, Og: 118,
};

// Ground-state spin multiplicity (2S+1) by atomic number — Hund's rule on the
// neutral atom electronic configuration. Z = 1–18 covered.
export const ATOMIC_GROUND_MULT: Record<number, number> = {
  1: 2,  2: 1,  3: 2,  4: 1,  5: 2,  6: 3,  7: 4,  8: 3,  9: 2, 10: 1,
  11: 2, 12: 1, 13: 2, 14: 3, 15: 4, 16: 3, 17: 2, 18: 1,
};

// Atomic number -> element name (1-indexed)
export const ATOMIC_NUMBER_TO_ELEMENT_NAME: string[] = [
  'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
  'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
  'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
  'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
  'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
  'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
  'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
  'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
  'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
  'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
  'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
  'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og',
];

export function elementNameToAtomicNumber(name: string): number {
  const n = ELEMENT_NAME_TO_ATOMIC_NUMBER[name];
  if (n === undefined) throw new Error(`Unknown element: ${name}`);
  return n;
}

export function atomicNumberToElementName(z: number): string {
  return ATOMIC_NUMBER_TO_ELEMENT_NAME[z - 1];
}

// Shell name <-> shell type
export const SHELL_NAME_TO_TYPE: Record<string, number> = {
  s: 0, p: 1, d: 2, f: 3, g: 4, h: 5, i: 6,
  S: 0, P: 1, D: 2, F: 3, G: 4, H: 5, I: 6,
};

export const SHELL_TYPE_TO_NAME = ['s', 'p', 'd', 'f', 'g', 'h', 'i'];

export function shellNameToType(name: string): number {
  const t = SHELL_NAME_TO_TYPE[name];
  if (t === undefined) throw new Error(`Unknown shell name: ${name}`);
  return t;
}

/** Number of Cartesian basis functions for a given shell type: (l+1)(l+2)/2 */
export function shellTypeToNumBasis(shellType: number): number {
  return ((shellType + 1) * (shellType + 2)) / 2;
}

/** Cartesian exponents (lx, ly, lz) for each basis function of each shell type */
export const ANGULAR_MOMENTUMS: number[][][] = [
  // s
  [[0, 0, 0]],
  // p
  [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
  // d
  [[2, 0, 0], [0, 2, 0], [0, 0, 2], [1, 1, 0], [1, 0, 1], [0, 1, 1]],
  // f
  [[3, 0, 0], [0, 3, 0], [0, 0, 3], [1, 2, 0], [2, 1, 0], [2, 0, 1], [1, 0, 2], [0, 1, 2], [0, 2, 1], [1, 1, 1]],
  // g
  [
    [4, 0, 0], [0, 4, 0], [0, 0, 4],
    [3, 1, 0], [3, 0, 1], [1, 3, 0], [0, 3, 1], [1, 0, 3], [0, 1, 3],
    [2, 2, 0], [2, 0, 2], [0, 2, 2],
    [2, 1, 1], [1, 2, 1], [1, 1, 2],
  ],
];
