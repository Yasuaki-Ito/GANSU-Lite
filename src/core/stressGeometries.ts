/**
 * Deterministic size-ramp geometries for the cross-device stress test.
 *
 * The stress test measures *where each device stops working*, so every device
 * must be fed byte-identical geometries. Everything here is generated from the
 * cluster index alone — no RNG, no external files — so `waterCluster(6)` yields
 * exactly the same XYZ on a desktop and on an iPhone.
 */

/** Experimental water monomer: r(OH) = 0.9572 Å, angle(HOH) = 104.52°. */
const R_OH = 0.9572;
const ANG_HOH = 104.52 * Math.PI / 180;

/** Simple-cubic lattice spacing (Å) between monomer centres (O atoms). */
const LATTICE_A = 3.0;

/** Monomer in its local frame: O at the origin, C2 axis along +z. */
const MONOMER: ReadonlyArray<readonly [string, number, number, number]> = [
  ['O', 0, 0, 0],
  ['H', R_OH * Math.sin(ANG_HOH / 2), 0, R_OH * Math.cos(ANG_HOH / 2)],
  ['H', -R_OH * Math.sin(ANG_HOH / 2), 0, R_OH * Math.cos(ANG_HOH / 2)],
];

/**
 * The n lattice sites closest to the origin, giving a compact (roughly
 * spherical) cluster. Ties are broken on (k, j, i) so the ordering is total
 * and platform-independent.
 */
function compactLatticeSites(n: number): Array<[number, number, number]> {
  const reach = Math.ceil(Math.cbrt(n)) + 1;
  const sites: Array<[number, number, number]> = [];
  for (let i = -reach; i <= reach; i++)
    for (let j = -reach; j <= reach; j++)
      for (let k = -reach; k <= reach; k++)
        sites.push([i, j, k]);
  sites.sort((a, b) => {
    const da = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
    const db = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
    if (da !== db) return da - db;
    if (a[2] !== b[2]) return a[2] - b[2];
    if (a[1] !== b[1]) return a[1] - b[1];
    return a[0] - b[0];
  });
  return sites.slice(0, n);
}

/**
 * Orientation of monomer `i`, as ZYZ Euler angles taken from a golden-angle
 * sequence. This spreads the monomer dipoles quasi-uniformly over orientations
 * without any randomness, so no two neighbours end up in a pathological
 * head-to-head arrangement.
 */
function orientation(i: number): [number, number, number] {
  const GOLDEN = Math.PI * (3 - Math.sqrt(5));  // ≈ 2.39996 rad
  return [i * GOLDEN, i * 1.2345678, i * 0.7654321];
}

function rotateZYZ(
  v: readonly [number, number, number],
  [a, b, c]: [number, number, number],
): [number, number, number] {
  const rz = (p: [number, number, number], t: number): [number, number, number] =>
    [p[0] * Math.cos(t) - p[1] * Math.sin(t), p[0] * Math.sin(t) + p[1] * Math.cos(t), p[2]];
  const ry = (p: [number, number, number], t: number): [number, number, number] =>
    [p[0] * Math.cos(t) + p[2] * Math.sin(t), p[1], -p[0] * Math.sin(t) + p[2] * Math.cos(t)];
  return rz(ry(rz([v[0], v[1], v[2]], c), b), a);
}

/** XYZ text for the water cluster (H2O)_n. */
export function waterCluster(n: number): string {
  const sites = compactLatticeSites(n);
  const lines: string[] = [];
  for (let i = 0; i < n; i++) {
    const [li, lj, lk] = sites[i];
    const angles = orientation(i);
    const ox = li * LATTICE_A, oy = lj * LATTICE_A, oz = lk * LATTICE_A;
    for (const [sym, x, y, z] of MONOMER) {
      const [rx, ry_, rz_] = rotateZYZ([x, y, z], angles);
      lines.push(
        `${sym.padEnd(2)} ${(ox + rx).toFixed(6).padStart(12)} ` +
        `${(oy + ry_).toFixed(6).padStart(12)} ${(oz + rz_).toFixed(6).padStart(12)}`,
      );
    }
  }
  return `${n * 3}\n(H2O)${n} — cubic lattice a=${LATTICE_A} A, golden-angle orientations\n${lines.join('\n')}`;
}

// ── n-alkane series (alternative ramp: C_n H_{2n+2}) ──

const R_CC = 1.526, R_CH = 1.09;
const TETRA = 109.47 * Math.PI / 180;

/** XYZ text for the all-anti n-alkane C_nH_{2n+2}. */
export function nAlkane(n: number): string {
  // All-anti backbone in the xz-plane: a planar zig-zag along x.
  const half = TETRA / 2;
  const dx = R_CC * Math.sin(half);
  const dz = R_CC * Math.cos(half);
  const C: Array<[number, number, number]> = [];
  for (let i = 0; i < n; i++) C.push([i * dx, 0, (i % 2 === 0 ? 0 : dz)]);

  const lines: string[] = [];
  for (const c of C) lines.push(fmtAtom('C', c));

  // Two H per carbon, perpendicular to the backbone plane (±y),
  // plus one extra H on each terminal carbon along the chain axis.
  for (let i = 0; i < n; i++) {
    const [cx, cy, cz] = C[i];
    const up = i % 2 === 0 ? 1 : -1;
    const hz = cz + up * R_CH * Math.cos(half);
    const hy = R_CH * Math.sin(half);
    lines.push(fmtAtom('H', [cx, cy + hy, hz]));
    lines.push(fmtAtom('H', [cx, cy - hy, hz]));
  }
  // Terminal hydrogens extend the zig-zag past each end carbon.
  const first = C[0], last = C[n - 1];
  const dirFirst: [number, number, number] = [-Math.sin(half), 0, (n > 1 && C[1][2] > first[2]) ? -Math.cos(half) : Math.cos(half)];
  const dirLast: [number, number, number] = [Math.sin(half), 0, (n > 1 && C[n - 2][2] > last[2]) ? -Math.cos(half) : Math.cos(half)];
  lines.push(fmtAtom('H', [first[0] + R_CH * dirFirst[0], first[1], first[2] + R_CH * dirFirst[2]]));
  lines.push(fmtAtom('H', [last[0] + R_CH * dirLast[0], last[1], last[2] + R_CH * dirLast[2]]));

  return `${n + 2 * n + 2}\nC${n}H${2 * n + 2} all-anti\n${lines.join('\n')}`;
}

function fmtAtom(sym: string, [x, y, z]: [number, number, number]): string {
  return `${sym.padEnd(2)} ${x.toFixed(6).padStart(12)} ${y.toFixed(6).padStart(12)} ${z.toFixed(6).padStart(12)}`;
}

// ── Series registry ──

export interface StressSeries {
  id: string;
  label: string;
  /** Human-readable name of member `n`, e.g. "(H2O)6". */
  name: (n: number) => string;
  /** Atom count of member `n`. */
  natoms: (n: number) => number;
  xyz: (n: number) => string;
  /** Default ladder of `n` values, finest near the expected breaking point. */
  ladder: number[];
}

export const STRESS_SERIES: StressSeries[] = [
  {
    id: 'water',
    label: 'Water cluster (H₂O)ₙ',
    name: (n) => `(H2O)${n}`,
    natoms: (n) => 3 * n,
    xyz: waterCluster,
    ladder: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 24, 28, 32],
  },
  {
    id: 'alkane',
    label: 'n-Alkane CₙH₂ₙ₊₂',
    name: (n) => `C${n}H${2 * n + 2}`,
    natoms: (n) => 3 * n + 2,
    xyz: nAlkane,
    ladder: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20],
  },
];
