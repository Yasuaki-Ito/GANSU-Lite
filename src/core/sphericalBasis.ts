/**
 * Cartesian → real-solid-harmonic transformation for d and f shells.
 *
 * Integrals are always computed over Cartesian Gaussians (6 d, 10 f components).
 * A spherical basis (5 d, 7 f) is the subspace spanned by the solid harmonics
 * r^l Y_lm, so it is reached with a constant matrix T (n_cart × n_sph): any
 * spherical-basis quantity is T^T (Cartesian quantity) T. T depends only on angles,
 * not on geometry, which is why gradients need no extra terms.
 *
 * GANSU normalises every Cartesian component individually (primitiveNorm carries
 * 1/sqrt((2a-1)!!(2b-1)!!(2c-1)!!)), unlike codes that give all components the
 * normalisation of x^l. The coefficients below are therefore built for that
 * convention, and each column is normalised with the analytic one-centre overlap
 * of the shell rather than taken from a table written for another code.
 *
 * Columns are ordered m = -l … +l, as in PySCF/libcint.
 */

import { Matrix } from '../linalg/matrix';
import { ANGULAR_MOMENTUMS } from './constants';

/** One contracted shell's place in the Cartesian basis. */
export interface ContractedShellLayout {
  /** First Cartesian basis-function index of the shell. */
  basisIndex: number;
  /** 0 = s, 1 = p, 2 = d, 3 = f, … */
  shellType: number;
}

/** (n-1)!! for even n ≥ 0, i.e. the angular integral factor; 0 for odd n. */
function evenMomentFactor(n: number): number {
  if (n % 2 !== 0) return 0;
  let r = 1;
  for (let k = n - 1; k > 1; k -= 2) r *= k;
  return r;
}

/** Monomial key "abc" → coefficient, for one solid harmonic. */
type Monomials = Record<string, number>;

/**
 * Real solid harmonics as polynomials in x, y, z (unnormalised), m = -l … +l.
 * Only d and f are tabulated: no basis shipped with GANSU has g functions, and a
 * silent wrong answer is worse than refusing.
 */
const SOLID_HARMONICS: Record<number, Monomials[]> = {
  2: [
    { '110': 1 },                                   // m=-2  xy
    { '011': 1 },                                   // m=-1  yz
    { '002': 2, '200': -1, '020': -1 },            // m= 0  2z² - x² - y²
    { '101': 1 },                                   // m=+1  xz
    { '200': 1, '020': -1 },                       // m=+2  x² - y²
  ],
  3: [
    { '210': 3, '030': -1 },                       // m=-3  3x²y - y³
    { '111': 1 },                                   // m=-2  xyz
    { '012': 4, '210': -1, '030': -1 },            // m=-1  4yz² - x²y - y³
    { '003': 2, '201': -3, '021': -3 },            // m= 0  2z³ - 3x²z - 3y²z
    { '102': 4, '300': -1, '120': -1 },            // m=+1  4xz² - x³ - xy²
    { '201': 1, '021': -1 },                       // m=+2  x²z - y²z
    { '300': 1, '120': -3 },                       // m=+3  x³ - 3xy²
  ],
};

/**
 * Columns of T for one shell, expressed over GANSU's individually-normalised
 * Cartesian components (ANGULAR_MOMENTUMS order). Each column has unit norm.
 */
export function shellSphericalColumns(l: number): number[][] {
  const comps = ANGULAR_MOMENTUMS[l];
  if (l < 2) {
    // s and p: Cartesian and spherical coincide.
    return comps.map((_, k) => comps.map((__, j) => (j === k ? 1 : 0)));
  }
  const harmonics = SOLID_HARMONICS[l];
  if (!harmonics) {
    throw new Error(`Spherical harmonics are not implemented for l = ${l}; use a Cartesian basis.`);
  }
  const key = (c: number[]) => `${c[0]}${c[1]}${c[2]}`;
  // Angular overlap of monomials x^a y^b z^c with x^a' y^b' z^c' (common radial
  // factor dropped — it is identical for every component of a shell).
  const G = (p: number[], q: number[]) =>
    evenMomentFactor(p[0] + q[0]) * evenMomentFactor(p[1] + q[1]) * evenMomentFactor(p[2] + q[2]);

  return harmonics.map(h => {
    const s = comps.map(c => h[key(c)] ?? 0);
    let norm2 = 0;
    for (let i = 0; i < comps.length; i++) {
      for (let j = 0; j < comps.length; j++) norm2 += s[i] * s[j] * G(comps[i], comps[j]);
    }
    // monomial = sqrt(G(c,c)) × (normalised component), up to the shared radial constant
    return comps.map((c, i) => s[i] * Math.sqrt(G(c, c)) / Math.sqrt(norm2));
  });
}

/** Number of spherical functions in a shell of angular momentum l. */
export const numSphericalInShell = (l: number) => 2 * l + 1;

/**
 * Build T (n_cart × n_sph) for a whole molecule. Returns null when the basis has
 * no shell with l ≥ 2, in which case Cartesian and spherical are the same basis.
 */
export function buildCartesianToSpherical(
  shells: ContractedShellLayout[],
  numCartesian: number,
): { T: Matrix; numSpherical: number } | null {
  if (!shells.some(s => s.shellType >= 2)) return null;
  const numSpherical = shells.reduce((acc, s) => acc + numSphericalInShell(s.shellType), 0);
  const T = new Matrix(numCartesian, numSpherical);
  let col = 0;
  for (const s of shells) {
    const cols = shellSphericalColumns(s.shellType);
    for (const v of cols) {
      for (let k = 0; k < v.length; k++) if (v[k] !== 0) T.set(s.basisIndex + k, col, v[k]);
      col++;
    }
  }
  return { T, numSpherical };
}

/**
 * Whether a basis set is conventionally used with pure (spherical) d/f functions.
 * Dunning (cc-pVXZ, aug-) and Karlsruhe (def2-) sets were optimised as spherical;
 * Pople and STO sets are conventionally Cartesian (6d), as in Gaussian's defaults.
 */
export function conventionallySpherical(basisName: string): boolean {
  const n = basisName.trim().toLowerCase();
  return /^(aug-)?(cc-p|pc-)/.test(n) || n.startsWith('def2-') || n.startsWith('aug-cc-');
}
