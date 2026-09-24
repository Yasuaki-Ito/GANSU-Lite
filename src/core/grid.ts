/** Molecular integration grid for DFT: product angular + Mura-Knowles radial + Becke partitioning */

import { Atom } from './types';

export interface GridPoint {
  x: number;
  y: number;
  z: number;
  weight: number;
}

// Bragg-Slater radii in Bohr (used for Becke partitioning size adjustment)
const BRAGG_SLATER_RADII: Record<number, number> = {
  1: 0.661, 2: 0.567,
  3: 2.646, 4: 1.890,
  5: 1.606, 6: 1.417, 7: 1.228, 8: 1.134, 9: 1.039, 10: 0.945,
  11: 3.402, 12: 2.646,
  13: 2.268, 14: 2.079, 15: 1.890, 16: 1.890, 17: 1.701, 18: 1.606,
  19: 4.158, 20: 3.402,
  21: 3.024, 22: 2.646, 23: 2.457, 24: 2.457, 25: 2.457, 26: 2.457,
  27: 2.268, 28: 2.268, 29: 2.268, 30: 2.457,
  31: 2.362, 32: 2.268, 33: 2.173, 34: 2.079, 35: 2.079, 36: 1.985,
};

function getBraggSlaterRadius(z: number): number {
  return BRAGG_SLATER_RADII[z] ?? 2.0;
}

/** Mura-Knowles radial scaling parameter R (Bohr).
 *  Mura & Knowles (J. Chem. Phys. 104, 9848 (1996)) use a fixed R for every
 *  element: 7.0 Bohr, lowered to 5.0 for the alkali / alkaline-earth metals
 *  (Li, Be, Na, Mg, K, Ca) whose valence density is very diffuse. PySCF uses
 *  the same convention.
 *
 *  Earlier versions used the Bragg-Slater radius as R (with a 5.0 special case
 *  for H). That under-samples the outer region of compact atoms: Ar has
 *  BS = 1.6 Bohr, so a 50-point grid placed only 2-3 radial shells beyond
 *  3 Bohr. In Ar2 at 3.8 Å the Becke cell extends to 3.6 Bohr, ρ was
 *  mis-integrated by 2.7e-2 electrons and the B3LYP curve showed a spurious
 *  -18 kJ/mol minimum; even the 75-point "fine" grid was 7 kJ/mol off.
 *  With R = 7.0 the same 50-point grid integrates ρ to 3e-4 and the curve is
 *  converged (-0.5 kJ/mol vs a 75-point reference at -0.02).
 *
 *  Hydrogen keeps R = 5.0 (the value this code has always used for H). H2 is
 *  already grid-converged with it, and R = 7.0 for H exposes a sensitivity of the
 *  PBE *potential* (not the energy) to the outer radial shells: H2/STO-3G PBE
 *  HOMO moves by +0.5 eV on the 50-point grid while E_xc is unchanged to 1e-8.
 *  That is a PBE-potential issue to be fixed separately; until then R = 5.0
 *  keeps the validated H behaviour. */
function getMuraKnowlesR(z: number): number {
  if (z === 1) return 5.0;
  return (z === 3 || z === 4 || z === 11 || z === 12 || z === 19 || z === 20) ? 5.0 : 7.0;
}

export type GridLevel = 'coarse' | 'medium' | 'fine';

interface GridSpec {
  nrad: number;
  ntheta: number;  // Gauss-Legendre points for θ
  /** Multiplier on the Mura-Knowles radial scale R (experiments only). */
  rScale?: number;
  /** Replace the per-element radial scale R(Z) entirely (experiments only). */
  rFunc?: (z: number) => number;
}

const GRID_SPECS: Record<GridLevel, GridSpec> = {
  coarse:  { nrad: 20, ntheta: 7 },   // ~98 angular pts/atom, ~1960 total/atom
  medium:  { nrad: 50, ntheta: 11 },   // ~242 angular pts/atom, ~12100 total/atom
  fine:    { nrad: 75, ntheta: 15 },   // ~450 angular pts/atom, ~33750 total/atom
};

// ── Gauss-Legendre quadrature nodes and weights on [-1, 1] ──

function gaussLegendre(n: number): { nodes: Float64Array; weights: Float64Array } {
  const nodes = new Float64Array(n);
  const weights = new Float64Array(n);

  for (let i = 0; i < n; i++) {
    // Initial guess
    let x = Math.cos(Math.PI * (i + 0.75) / (n + 0.5));

    // Newton iterations to find root of P_n(x)
    for (let iter = 0; iter < 100; iter++) {
      let p0 = 1, p1 = x;
      for (let j = 2; j <= n; j++) {
        const p2 = ((2 * j - 1) * x * p1 - (j - 1) * p0) / j;
        p0 = p1;
        p1 = p2;
      }
      // p1 = P_n(x), derivative: n*(x*P_n - P_{n-1})/(x²-1)
      const dp = n * (x * p1 - p0) / (x * x - 1);
      const dx = p1 / dp;
      x -= dx;
      if (Math.abs(dx) < 1e-15) break;
    }

    nodes[i] = x;
    // Weight: 2 / ((1-x²) * [P'_n(x)]²)
    let p0 = 1, p1 = x;
    for (let j = 2; j <= n; j++) {
      const p2 = ((2 * j - 1) * x * p1 - (j - 1) * p0) / j;
      p0 = p1;
      p1 = p2;
    }
    const dp = n * (x * p1 - p0) / (x * x - 1);
    weights[i] = 2 / ((1 - x * x) * dp * dp);
  }

  return { nodes, weights };
}

// Product angular grid: Gauss-Legendre for θ, uniform for φ
// Returns points on the unit sphere with weights summing to 4π
function productAngularGrid(ntheta: number): { x: number; y: number; z: number; w: number }[] {
  const nphi = 2 * ntheta;
  const gl = gaussLegendre(ntheta);
  const pts: { x: number; y: number; z: number; w: number }[] = [];

  for (let i = 0; i < ntheta; i++) {
    const cosTheta = gl.nodes[i];
    const sinTheta = Math.sqrt(1 - cosTheta * cosTheta);
    // GL weight for cos(θ) integration × 2π/nphi for φ integration
    const wTheta = gl.weights[i];

    for (let j = 0; j < nphi; j++) {
      const phi = 2 * Math.PI * j / nphi;
      const wPhi = 2 * Math.PI / nphi;
      pts.push({
        x: sinTheta * Math.cos(phi),
        y: sinTheta * Math.sin(phi),
        z: cosTheta,
        w: wTheta * wPhi,  // sum = 4π
      });
    }
  }
  return pts;
}

// ── Mura-Knowles radial grid ──
// r_i = -R * ln(1 - x³), x = (i+0.5)/nrad
function muraKnowlesRadial(nrad: number, braggR: number): { r: number; w: number }[] {
  const points: { r: number; w: number }[] = [];
  for (let i = 0; i < nrad; i++) {
    const x = (i + 0.5) / nrad;
    const x3 = x * x * x;
    const r = -braggR * Math.log(1 - x3);
    const dr = 3 * braggR * x * x / (1 - x3) / nrad;
    points.push({ r, w: dr * r * r }); // r² dr Jacobian
  }
  return points;
}

// ── Becke partitioning ──
function beckeStep(mu: number): number {
  let s = mu;
  s = 1.5 * s - 0.5 * s * s * s;
  s = 1.5 * s - 0.5 * s * s * s;
  s = 1.5 * s - 0.5 * s * s * s;
  return s;
}

function computeBeckeWeights(
  atoms: Atom[],
  px: number, py: number, pz: number,
  distCache: Float64Array,  // natom × natom inter-atomic distances
): Float64Array {
  const natom = atoms.length;
  if (natom === 1) return new Float64Array([1.0]);

  const raw = new Float64Array(natom);
  raw.fill(1.0);

  // Distances from grid point to each atom
  const dist = new Float64Array(natom);
  for (let i = 0; i < natom; i++) {
    const ai = atoms[i].coordinate;
    const dx = px - ai.x, dy = py - ai.y, dz = pz - ai.z;
    dist[i] = Math.sqrt(dx * dx + dy * dy + dz * dz);
  }

  for (let i = 0; i < natom; i++) {
    for (let j = 0; j < natom; j++) {
      if (i === j) continue;
      const dij = distCache[i * natom + j];
      if (dij < 1e-14) continue;
      let mu = (dist[i] - dist[j]) / dij;

      // Bragg-Slater size adjustment
      const chi = getBraggSlaterRadius(atoms[i].atomicNumber) /
                  getBraggSlaterRadius(atoms[j].atomicNumber);
      const u = (chi - 1) / (chi + 1);
      const a = u / (u * u - 1);
      const aij = Math.max(-0.5, Math.min(0.5, a));
      mu += aij * (1 - mu * mu);

      raw[i] *= 0.5 * (1 - beckeStep(mu));
    }
  }

  let sum = 0;
  for (let i = 0; i < natom; i++) sum += raw[i];
  if (sum > 1e-30) {
    for (let i = 0; i < natom; i++) raw[i] /= sum;
  }
  return raw;
}

/** Build molecular integration grid */
export function buildMolecularGrid(
  atoms: Atom[],
  level: GridLevel = 'medium',
  overrideSpec?: Partial<GridSpec>,
): GridPoint[] {
  const spec: GridSpec = { ...GRID_SPECS[level] };
  if (overrideSpec) {
    for (const [k, v] of Object.entries(overrideSpec)) {
      if (v !== undefined) (spec as unknown as Record<string, unknown>)[k] = v;
    }
  }
  const angularGrid = productAngularGrid(spec.ntheta);
  const natom = atoms.length;

  // Pre-compute inter-atomic distances
  const distCache = new Float64Array(natom * natom);
  for (let i = 0; i < natom; i++) {
    const ai = atoms[i].coordinate;
    for (let j = i + 1; j < natom; j++) {
      const aj = atoms[j].coordinate;
      const dx = ai.x - aj.x, dy = ai.y - aj.y, dz = ai.z - aj.z;
      const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
      distCache[i * natom + j] = d;
      distCache[j * natom + i] = d;
    }
  }

  const points: GridPoint[] = [];

  for (let iAtom = 0; iAtom < natom; iAtom++) {
    const atom = atoms[iAtom];
    const cx = atom.coordinate.x, cy = atom.coordinate.y, cz = atom.coordinate.z;
    const radialR = (spec.rFunc ?? getMuraKnowlesR)(atom.atomicNumber) * (spec.rScale ?? 1);
    const radialGrid = muraKnowlesRadial(spec.nrad, radialR);

    for (const rp of radialGrid) {
      for (const ap of angularGrid) {
        const x = cx + rp.r * ap.x;
        const y = cy + rp.r * ap.y;
        const z = cz + rp.r * ap.z;

        const beckeW = computeBeckeWeights(atoms, x, y, z, distCache);
        const totalWeight = rp.w * ap.w * beckeW[iAtom];

        if (totalWeight > 1e-30) {
          points.push({ x, y, z, weight: totalWeight });
        }
      }
    }
  }

  return points;
}
