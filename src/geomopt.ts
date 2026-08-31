/**
 * GANSU Lite — Geometry Optimization demo page.
 * Algorithms: SD, CG (FR/PR/HS/DY), Quasi-Newton (BFGS/DFP/SR1), GDIIS
 * With Armijo backtracking line search and translation/rotation projection.
 */

import { parseXYZ } from './core/parseXYZ';
import { atomicNumberToElementName } from './core/constants';
import { BasisSet } from './core/basisSet';
import { Molecular } from './core/molecular';
import { theorySelectHTML, buildHFOrDFT, theoryDisplayLabel, HEAVY_ITERATIVE_EXCLUDE, type TheoryChoice } from './ui/theoryControls';
import { RHF } from './core/rhf';
import { computeRHFGradient } from './core/gradient';
// Hessian import kept for potential future use
// import { computeNumericalHessian, jacobiEigen as jacobiEigenLocal } from './core/hessian';
import { initTheme, toggleTheme, getThemeColors } from './ui/theme';
import { t, initLang, toggleLang } from './ui/i18n';
import { initWasm, getActiveBackend } from './core/eriWasm';
import { renderHeader } from './ui/nav';
import { createGeomOptViewer, type GeomOptViewer } from './ui/moleculeViewer3D';

// ── Scenario definitions ─────────────────────────────────────────────

interface OptScenario {
  id: string; labelKey: string; descKey: string;
  charge: number; mult: number; xyz: string; basisName: string;
}

const B = 1.8897259886;

function xyz(atoms: [string, number, number, number][]): string {
  return `${atoms.length}\ngeomopt\n${atoms.map(([el, x, y, z]) => `${el}  ${x.toFixed(6)}  ${y.toFixed(6)}  ${z.toFixed(6)}`).join('\n')}`;
}

const SCENARIOS: OptScenario[] = [
  { id: 'h2', labelKey: 'gopt.scenH2', descKey: 'gopt.descH2', charge: 0, mult: 1, basisName: 'sto-3g',
    xyz: xyz([['H', 0, 0, 0], ['H', 0, 0, 2.0 * B]]) },
  { id: 'hf', labelKey: 'gopt.scenHF', descKey: 'gopt.descHF', charge: 0, mult: 1, basisName: 'sto-3g',
    xyz: xyz([['H', 0, 0, 0], ['F', 0, 0, 2.2 * B]]) },
  { id: 'lih', labelKey: 'gopt.scenLiH', descKey: 'gopt.descLiH', charge: 0, mult: 1, basisName: 'sto-3g',
    xyz: xyz([['Li', 0, 0, 0], ['H', 0, 0, 2.0 * B]]) },
  { id: 'h2o', labelKey: 'gopt.scenH2O', descKey: 'gopt.descH2O', charge: 0, mult: 1, basisName: 'sto-3g',
    xyz: xyz([['O', 0, 0, 0], ['H', 1.1*B, 0, 0], ['H', 0, 1.1*B, 0]]) },
  { id: 'nh3', labelKey: 'gopt.scenNH3', descKey: 'gopt.descNH3', charge: 0, mult: 1, basisName: 'sto-3g',
    xyz: (() => { const R=1.05*B, th=95*Math.PI/180, z=R*Math.cos(th), r=R*Math.sin(th);
      return xyz([['N',0,0,0],['H',r,0,z],['H',-r*.5,r*Math.sqrt(3)/2,z],['H',-r*.5,-r*Math.sqrt(3)/2,z]]); })() },
  { id: 'beh2', labelKey: 'gopt.scenBeH2', descKey: 'gopt.descBeH2', charge: 0, mult: 1, basisName: 'sto-3g',
    xyz: xyz([['Be', 0, 0, 0], ['H', 1.4*B, 0, 1.0*B], ['H', -1.4*B, 0, 1.0*B]]) },
  { id: 'hcn', labelKey: 'gopt.scenHCN', descKey: 'gopt.descHCN', charge: 0, mult: 1, basisName: 'sto-3g',
    xyz: xyz([['H', 0, 0, 0], ['C', 0, 0, 1.2*B], ['N', 0, 0, 1.2*B+1.3*B]]) },
  { id: 'ch4', labelKey: 'gopt.scenCH4', descKey: 'gopt.descCH4', charge: 0, mult: 1, basisName: 'sto-3g',
    xyz: (() => { const r=1.15*B, a=100*Math.PI/180, s=Math.sin(a), c=Math.cos(a);
      return xyz([['C',0,0,0],['H',r*s,0,r*c],['H',-r*s*.5,r*s*Math.sqrt(3)/2,r*c],
        ['H',-r*s*.5,-r*s*Math.sqrt(3)/2,r*c],['H',0,0,-r]]); })() },
  { id: 'bh3', labelKey: 'gopt.scenBH3', descKey: 'gopt.descBH3', charge: 0, mult: 1, basisName: 'sto-3g',
    xyz: (() => { const r=1.25*B;
      return xyz([['B',0,0,0],['H',r,0,.5*B],['H',-r*.5,r*Math.sqrt(3)/2,.5*B],['H',-r*.5,-r*Math.sqrt(3)/2,.5*B]]); })() },
];

// ── Algorithm types ──────────────────────────────────────────────────

type AlgoId = 'sd' | 'cg-fr' | 'cg-pr' | 'cg-hs' | 'cg-dy' | 'bfgs' | 'dfp' | 'sr1' | 'gdiis';

const ALGO_GROUPS: { label: string; items: { id: AlgoId; name: string }[] }[] = [
  { label: 'Steepest Descent', items: [{ id: 'sd', name: 'SD' }] },
  { label: 'Conjugate Gradient', items: [
    { id: 'cg-fr', name: 'CG-FR' }, { id: 'cg-pr', name: 'CG-PR' },
    { id: 'cg-hs', name: 'CG-HS' }, { id: 'cg-dy', name: 'CG-DY' }] },
  { label: 'Quasi-Newton', items: [
    { id: 'bfgs', name: 'BFGS' }, { id: 'dfp', name: 'DFP' }, { id: 'sr1', name: 'SR1' }] },
  { label: 'GDIIS', items: [{ id: 'gdiis', name: 'GDIIS' }] },
];

// ── Types & State ────────────────────────────────────────────────────

interface OptStep {
  energy: number; maxForce: number; rmsForce: number;
  coords: number[]; gradient: number[];
}

let selectedScen = SCENARIOS[0];
let selectedAlgo: AlgoId = 'bfgs';
let steps: OptStep[] = [];
let running = false, stopRequested = false, done = false, converged = false;
let elapsedMs = 0, maxIter = 50, selectedIdx = 0;
let theoryChoice: TheoryChoice = 'HF';
const forceThreshold = 4.5e-4;
let geomViewer: GeomOptViewer | null = null;
let eriOpt: 'wasm' | 'auto' = 'auto';

// ── Basis loading ────────────────────────────────────────────────────

const basisCache = new Map<string, BasisSet>();
async function loadBasis(name: string): Promise<BasisSet> {
  const c = basisCache.get(name);
  if (c) return c;
  const url = `${import.meta.env.BASE_URL}basis/${name.toLowerCase()}.gbs`;
  const bs = BasisSet.fromGBS(await (await fetch(url)).text());
  basisCache.set(name, bs);
  return bs;
}

// ── Translation & Rotation Projection ────────────────────────────────

/** Project out translation and rotation from gradient (6 DOF, or 5 for linear). */
function projectTR(grad: Float64Array, coords: Float64Array, nAtom: number): void {
  const n = 3 * nAtom;
  // Build translation vectors (3)
  const tvecs: Float64Array[] = [];
  for (let d = 0; d < 3; d++) {
    const v = new Float64Array(n);
    for (let i = 0; i < nAtom; i++) v[3 * i + d] = 1;
    tvecs.push(v);
  }
  // Build rotation vectors (3): R_x cross r_i, etc.
  // Center of mass (equal weight for simplicity)
  let cx = 0, cy = 0, cz = 0;
  for (let i = 0; i < nAtom; i++) { cx += coords[3*i]; cy += coords[3*i+1]; cz += coords[3*i+2]; }
  cx /= nAtom; cy /= nAtom; cz /= nAtom;

  const rvecs: Float64Array[] = [];
  for (let axis = 0; axis < 3; axis++) {
    const v = new Float64Array(n);
    for (let i = 0; i < nAtom; i++) {
      const rx = coords[3*i] - cx, ry = coords[3*i+1] - cy, rz = coords[3*i+2] - cz;
      // cross product of unit axis vector with r
      if (axis === 0) { v[3*i+1] = -rz; v[3*i+2] = ry; }       // x-axis: (0,−z,y)
      else if (axis === 1) { v[3*i] = rz; v[3*i+2] = -rx; }     // y-axis: (z,0,−x)
      else { v[3*i] = -ry; v[3*i+1] = rx; }                      // z-axis: (−y,x,0)
    }
    rvecs.push(v);
  }

  // Gram-Schmidt orthonormalize all 6 vectors, skip near-zero (linear molecule)
  const basis: Float64Array[] = [];
  for (const raw of [...tvecs, ...rvecs]) {
    const v = new Float64Array(raw);
    for (const b of basis) {
      let dot = 0;
      for (let i = 0; i < n; i++) dot += v[i] * b[i];
      for (let i = 0; i < n; i++) v[i] -= dot * b[i];
    }
    let norm = 0;
    for (let i = 0; i < n; i++) norm += v[i] * v[i];
    norm = Math.sqrt(norm);
    if (norm > 1e-10) {
      for (let i = 0; i < n; i++) v[i] /= norm;
      basis.push(v);
    }
  }

  // Project gradient: g = g - Σ (g·b) b
  for (const b of basis) {
    let dot = 0;
    for (let i = 0; i < n; i++) dot += grad[i] * b[i];
    for (let i = 0; i < n; i++) grad[i] -= dot * b[i];
  }
}

// ── Armijo backtracking line search ──────────────────────────────────

async function armijoLineSearch(
  coords: Float64Array, dir: Float64Array, energy: number, grad: Float64Array,
  atomicNumbers: number[], basis: BasisSet, charge: number,
  initialStep: number,
  theory: TheoryChoice,
): Promise<{ coords: Float64Array; energy: number; grad: Float64Array; step: number }> {
  const n = coords.length;
  const c1 = 1e-4; // sufficient decrease parameter
  let dirGrad = 0;
  for (let i = 0; i < n; i++) dirGrad += dir[i] * grad[i];
  if (dirGrad >= 0) dirGrad = -1e-8; // safety: ensure descent

  let alpha = initialStep;
  for (let trial = 0; trial < 6; trial++) {
    const newCoords = new Float64Array(n);
    for (let i = 0; i < n; i++) newCoords[i] = coords[i] + alpha * dir[i];

    const atoms = atomicNumbers.map((z, i) => ({
      atomicNumber: z,
      coordinate: { x: newCoords[3*i], y: newCoords[3*i+1], z: newCoords[3*i+2] },
      atomIndex: i,
    }));
    const mol = new Molecular(atoms, basis, charge);
    const rhf = await buildHFOrDFT(mol, basis, theory, 'RHF') as RHF;
    const newEnergy = await rhf.solve({ eriBackend: eriOpt });

    if (newEnergy <= energy + c1 * alpha * dirGrad) {
      // Armijo condition satisfied — compute gradient at accepted point
      const nocc = mol.numAlphaSpins;
      const dftCtx = (rhf.xcFunctional && rhf.grid) ? { functional: rhf.xcFunctional, grid: rhf.grid } : undefined;
      const gr = computeRHFGradient(mol.primitiveShells, mol.atoms, mol.cgtoNormalizationFactors,
        mol.numBasis, nocc, rhf.density, rhf.coefficients, rhf.orbitalEnergies, undefined, dftCtx);
      projectTR(gr.total, newCoords, atomicNumbers.length);
      return { coords: newCoords, energy: newEnergy, grad: gr.total, step: alpha };
    }
    alpha *= 0.5;
    await new Promise<void>(r => setTimeout(r, 0));
  }
  // Fallback: accept last trial
  const newCoords = new Float64Array(n);
  for (let i = 0; i < n; i++) newCoords[i] = coords[i] + alpha * dir[i];
  return { coords: newCoords, energy, grad, step: alpha };
}

// ── Conjugate Gradient beta formulas ─────────────────────────────────

function cgBeta(id: AlgoId, g: Float64Array, gPrev: Float64Array, dPrev: Float64Array): number {
  const n = g.length;
  let gg = 0, gpgp = 0, gDg = 0, gpd = 0, gd = 0;
  for (let i = 0; i < n; i++) {
    const dy = g[i] - gPrev[i];
    gg += g[i] * g[i];
    gpgp += gPrev[i] * gPrev[i];
    gDg += g[i] * dy;
    gpd += dy * dPrev[i];
    gd += g[i] * dPrev[i];
  }
  if (gpgp < 1e-30) return 0;
  switch (id) {
    case 'cg-fr': return gg / gpgp;                                  // Fletcher-Reeves
    case 'cg-pr': return Math.max(0, gDg / gpgp);                   // Polak-Ribière+
    case 'cg-hs': return gpd > 1e-30 ? gDg / gpd : 0;              // Hestenes-Stiefel
    case 'cg-dy': return gpd > 1e-30 ? gg / gpd : 0;               // Dai-Yuan
    default: return 0;
  }
}

// ── Quasi-Newton Hessian updates ─────────────────────────────────────

function updateHessInv(id: AlgoId, H: Float64Array, s: Float64Array, y: Float64Array, n: number): Float64Array {
  let sy = 0;
  for (let i = 0; i < n; i++) sy += s[i] * y[i];

  if (id === 'sr1') {
    // SR1: H += (s - Hy)(s - Hy)^T / (s - Hy)^T y
    const sHy = new Float64Array(n);
    for (let i = 0; i < n; i++) { let v = s[i]; for (let j = 0; j < n; j++) v -= H[i*n+j]*y[j]; sHy[i] = v; }
    let den = 0;
    for (let i = 0; i < n; i++) den += sHy[i] * y[i];
    if (Math.abs(den) < 1e-12 * Math.sqrt(sy > 0 ? sy : 1)) return H; // skip if denominator too small
    const newH = new Float64Array(H);
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++)
        newH[i*n+j] += sHy[i] * sHy[j] / den;
    return newH;
  }

  if (sy < 1e-14) return H;
  const rho = 1 / sy;

  const Hy = new Float64Array(n);
  for (let i = 0; i < n; i++) { let v = 0; for (let j = 0; j < n; j++) v += H[i*n+j]*y[j]; Hy[i] = v; }
  let yHy = 0;
  for (let i = 0; i < n; i++) yHy += y[i] * Hy[i];

  const newH = new Float64Array(n * n);

  if (id === 'bfgs') {
    // BFGS: (I - ρsy^T)H(I - ρys^T) + ρss^T
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++)
        newH[i*n+j] = H[i*n+j] - rho*(s[i]*Hy[j]+Hy[i]*s[j]) + rho*(rho*yHy+1)*s[i]*s[j];
  } else {
    // DFP: H - (Hy)(Hy)^T/(yHy) + ss^T/(sy)
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++)
        newH[i*n+j] = H[i*n+j] - Hy[i]*Hy[j]/(yHy||1) + rho*s[i]*s[j];
  }
  return newH;
}

// ── GDIIS ────────────────────────────────────────────────────────────

interface GDIISEntry { coords: Float64Array; grad: Float64Array; }
const gdiisHistory: GDIISEntry[] = [];
const GDIIS_MAX = 6;

function gdiisExtrapolate(n: number): Float64Array | null {
  const m = gdiisHistory.length;
  if (m < 2) return null;

  // Build B matrix: B_ij = e_i · e_j (error vectors = gradients)
  const B = new Float64Array((m + 1) * (m + 1));
  for (let i = 0; i < m; i++)
    for (let j = 0; j < m; j++) {
      let dot = 0;
      for (let k = 0; k < n; k++) dot += gdiisHistory[i].grad[k] * gdiisHistory[j].grad[k];
      B[i * (m + 1) + j] = dot;
    }
  // Lagrange constraint row/col
  for (let i = 0; i < m; i++) { B[m * (m + 1) + i] = -1; B[i * (m + 1) + m] = -1; }
  B[m * (m + 1) + m] = 0;

  // Solve B c = [0,...,0,-1] by Gaussian elimination
  const rhs = new Float64Array(m + 1);
  rhs[m] = -1;

  // Augmented matrix
  const dim = m + 1;
  const A = new Float64Array(dim * (dim + 1));
  for (let i = 0; i < dim; i++) {
    for (let j = 0; j < dim; j++) A[i * (dim + 1) + j] = B[i * dim + j];
    A[i * (dim + 1) + dim] = rhs[i];
  }

  // Partial pivoting Gauss elimination
  for (let col = 0; col < dim; col++) {
    let maxRow = col, maxVal = Math.abs(A[col * (dim + 1) + col]);
    for (let row = col + 1; row < dim; row++) {
      const v = Math.abs(A[row * (dim + 1) + col]);
      if (v > maxVal) { maxVal = v; maxRow = row; }
    }
    if (maxVal < 1e-14) return null;
    if (maxRow !== col) {
      for (let j = 0; j <= dim; j++) {
        const tmp = A[col * (dim + 1) + j];
        A[col * (dim + 1) + j] = A[maxRow * (dim + 1) + j];
        A[maxRow * (dim + 1) + j] = tmp;
      }
    }
    const pivot = A[col * (dim + 1) + col];
    for (let j = col; j <= dim; j++) A[col * (dim + 1) + j] /= pivot;
    for (let row = 0; row < dim; row++) {
      if (row === col) continue;
      const f = A[row * (dim + 1) + col];
      for (let j = col; j <= dim; j++) A[row * (dim + 1) + j] -= f * A[col * (dim + 1) + j];
    }
  }

  const c = new Float64Array(m);
  for (let i = 0; i < m; i++) c[i] = A[i * (dim + 1) + dim];

  // Extrapolated coordinates
  const result = new Float64Array(n);
  for (let i = 0; i < m; i++)
    for (let k = 0; k < n; k++)
      result[k] += c[i] * gdiisHistory[i].coords[k];
  return result;
}

// ── Optimization loop ────────────────────────────────────────────────

async function runOptimization(): Promise<void> {
  running = true; stopRequested = false; steps = []; done = false; converged = false; selectedIdx = 0;
  gdiisHistory.length = 0;
  render();

  const scen = selectedScen;
  const startTime = performance.now();
  const base = (import.meta as any).env?.BASE_URL || '/';
  await initWasm(base);
  const _backend = getActiveBackend();
  eriOpt = _backend !== 'js' ? 'wasm' : 'auto';
  const basis = await loadBasis(scen.basisName);
  const initAtoms = parseXYZ(scen.xyz);
  const nAtom = initAtoms.length;
  const n = 3 * nAtom;
  let coords: Float64Array<ArrayBufferLike> = new Float64Array(n);
  for (let i = 0; i < nAtom; i++) {
    coords[3*i] = initAtoms[i].coordinate.x;
    coords[3*i+1] = initAtoms[i].coordinate.y;
    coords[3*i+2] = initAtoms[i].coordinate.z;
  }
  const atomicNumbers = initAtoms.map(a => a.atomicNumber);

  let stepSize = 0.3;
  let prevGrad: Float64Array | null = null;
  let prevDir: Float64Array | null = null;
  let hessInv: Float64Array | null = null;
  let hessScaled = false;
  const isCG = selectedAlgo.startsWith('cg-');
  const isQN = selectedAlgo === 'bfgs' || selectedAlgo === 'dfp' || selectedAlgo === 'sr1';
  const isGDIIS = selectedAlgo === 'gdiis';

  for (let iter = 0; iter < maxIter; iter++) {
    if (stopRequested) break;

    // SCF + gradient
    const atoms = atomicNumbers.map((z, i) => ({
      atomicNumber: z,
      coordinate: { x: coords[3*i], y: coords[3*i+1], z: coords[3*i+2] },
      atomIndex: i,
    }));
    const mol = new Molecular(atoms, basis, scen.charge);
    const rhf = await buildHFOrDFT(mol, basis, theoryChoice, 'RHF') as RHF;
    const energy = await rhf.solve({ eriBackend: eriOpt });
    const nocc = mol.numAlphaSpins;
    const dftCtx = (rhf.xcFunctional && rhf.grid) ? { functional: rhf.xcFunctional, grid: rhf.grid } : undefined;
    const gResult = computeRHFGradient(mol.primitiveShells, mol.atoms, mol.cgtoNormalizationFactors,
      mol.numBasis, nocc, rhf.density, rhf.coefficients, rhf.orbitalEnergies, undefined, dftCtx);
    const grad = gResult.total;

    // Project out translation & rotation
    projectTR(grad, coords as Float64Array, nAtom);

    let maxF = 0, rmsSum = 0;
    for (let i = 0; i < n; i++) { maxF = Math.max(maxF, Math.abs(grad[i])); rmsSum += grad[i] * grad[i]; }

    steps.push({ energy, maxForce: maxF, rmsForce: Math.sqrt(rmsSum / n), coords: Array.from(coords), gradient: Array.from(grad) });
    selectedIdx = steps.length - 1;
    updateProgress(iter + 1, maxF);
    renderGraphs(); setupSlider(); renderMolVis();

    if (maxF < forceThreshold) { converged = true; break; }

    // Energy went up? Revert to previous geometry and retry with smaller step
    // (skip for Newton and tuned QN which have their own trust region)
    if (steps.length >= 2 && energy > steps[steps.length - 2].energy && (!isQN || !hessScaled)) {
      const prev = steps[steps.length - 2];
      for (let i = 0; i < n; i++) coords[i] = prev.coords[i];
      stepSize *= 0.5;
      prevDir = null;
      await new Promise<void>(r => setTimeout(r, 0));
      continue;
    }

    // ── Compute search direction ──
    let dir = new Float64Array(n);

    if (selectedAlgo === 'sd') {
      for (let i = 0; i < n; i++) dir[i] = -grad[i];

    } else if (isCG) {
      let beta = 0;
      if (prevGrad && prevDir) beta = cgBeta(selectedAlgo, grad, prevGrad, prevDir);
      for (let i = 0; i < n; i++) dir[i] = -grad[i] + beta * (prevDir ? prevDir[i] : 0);
      let dg = 0;
      for (let i = 0; i < n; i++) dg += dir[i] * grad[i];
      if (dg > 0) for (let i = 0; i < n; i++) dir[i] = -grad[i];

    } else if (isQN) {
      if (!hessInv) {
        hessInv = new Float64Array(n * n);
        for (let i = 0; i < n; i++) hessInv[i * n + i] = 1.0;
      }
      if (prevGrad && steps.length >= 2) {
        const pc = steps[steps.length - 2].coords;
        const s = new Float64Array(n), y = new Float64Array(n);
        for (let i = 0; i < n; i++) { s[i] = (coords as Float64Array)[i] - pc[i]; y[i] = grad[i] - prevGrad[i]; }
        // grad is already TR-projected, so y = grad - prevGrad is also projected.
        // s may have TR components from numerical noise; project it too.
        projectTR(s, coords as Float64Array, nAtom);
        let sy = 0, yy = 0, ss = 0;
        for (let i = 0; i < n; i++) { sy += s[i] * y[i]; yy += y[i] * y[i]; ss += s[i] * s[i]; }
        if (sy > 1e-14 && ss > 1e-20) {
          // Scale initial Hessian once: H0 = (s·y)/(y·y) * I
          if (!hessScaled && yy > 1e-14) {
            const gamma = sy / yy;
            for (let i = 0; i < n; i++) hessInv![i * n + i] = gamma;
            hessScaled = true;
          }
          hessInv = updateHessInv(selectedAlgo, hessInv!, s, y, n);
        }
      }
      // dir = -H⁻¹ g  (grad is already TR-projected)
      for (let i = 0; i < n; i++) {
        let v = 0;
        for (let j = 0; j < n; j++) v += hessInv![i * n + j] * (-grad[j]);
        dir[i] = v;
      }
      // Ensure descent; if not, reset Hessian and use SD
      let dg = 0;
      for (let i = 0; i < n; i++) dg += dir[i] * grad[i];
      if (dg > 0) {
        hessInv = new Float64Array(n * n);
        for (let i = 0; i < n; i++) hessInv[i * n + i] = 1.0;
        hessScaled = false;
        for (let i = 0; i < n; i++) dir[i] = -grad[i];
      }

    } else if (isGDIIS) {
      gdiisHistory.push({ coords: Float64Array.from(coords as Float64Array), grad: Float64Array.from(grad) });
      if (gdiisHistory.length > GDIIS_MAX) gdiisHistory.shift();
      const extrap = gdiisExtrapolate(n);
      if (extrap) {
        for (let i = 0; i < n; i++) dir[i] = extrap[i] - (coords as Float64Array)[i];
        const dNorm = Math.sqrt(dir.reduce((s, v) => s + v * v, 0));
        let dg = 0;
        for (let i = 0; i < n; i++) dg += dir[i] * grad[i];
        if (dNorm < 1e-12 || dg > 0) {
          for (let i = 0; i < n; i++) dir[i] = -grad[i];
        } else {
          // GDIIS: use extrapolated point directly
          const newCoords = new Float64Array(n);
          for (let i = 0; i < n; i++) newCoords[i] = (coords as Float64Array)[i] + dir[i];
          coords = newCoords;
          prevGrad = Float64Array.from(grad);
          prevDir = dir;
          await new Promise<void>(r => setTimeout(r, 0));
          continue;
        }
      } else {
        for (let i = 0; i < n; i++) dir[i] = -grad[i];
      }
    }

    // ── Take step ──
    const trustRadius = 0.5; // bohr max step
    const norm = Math.sqrt(dir.reduce((s, v) => s + v * v, 0));

    if (isQN && hessScaled) {
      // QN with tuned Hessian: dir has natural scale; limit max step
      const scale = norm > trustRadius ? trustRadius / norm : 1.0;
      const newCoords = new Float64Array(n);
      for (let i = 0; i < n; i++) newCoords[i] = (coords as Float64Array)[i] + scale * dir[i];
      coords = newCoords;
    } else {
      // SD, CG, GDIIS, or QN before first Hessian update: normalize and use stepSize
      const scale = Math.min(stepSize, trustRadius) / (norm || 1);
      const newCoords = new Float64Array(n);
      for (let i = 0; i < n; i++) newCoords[i] = (coords as Float64Array)[i] + scale * dir[i];
      coords = newCoords;
    }

    // Grow step size gradually if energy decreased
    if (!isQN || !hessScaled) {
      stepSize = Math.min(stepSize * 1.2, trustRadius);
    }

    prevGrad = Float64Array.from(grad);
    prevDir = dir;
    await new Promise<void>(r => setTimeout(r, 0));
  }

  elapsedMs = performance.now() - startTime;
  running = false; done = true;
  render();
}

// ── Styles ───────────────────────────────────────────────────────────

let stylesInjected = false;
function injectStyles(): void {
  if (stylesInjected) return;
  stylesInjected = true;
  const style = document.createElement('style');
  style.textContent = `
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body { background: var(--color-bg); color: var(--color-text);
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; min-height: 100vh; }
    .opt-page { max-width: 980px; margin: 0 auto; padding: 16px 20px; }
    .opt-content { display: flex; gap: 20px; }
    .opt-panel { background: var(--color-surface); border: 1px solid var(--color-border); border-radius: 10px; padding: 16px 18px; }
    .opt-controls { flex: 0 0 280px; }
    .opt-graph-panel { flex: 1; min-width: 0; display: flex; flex-direction: column; align-items: center; padding-top: 8px; }
    .opt-controls h2 {
      font-size: 0.72rem; font-weight: 600; color: var(--color-text-secondary);
      text-transform: uppercase; letter-spacing: 0.04em; margin: 10px 0 5px; }
    .opt-controls h2:first-child { margin-top: 0; }
    .opt-scenario-grid { display: flex; flex-direction: column; gap: 5px; }
    .opt-scenario-card {
      padding: 7px 12px; border: 2px solid var(--color-border); border-radius: 8px;
      cursor: pointer; transition: all 0.15s; display: flex; align-items: baseline; gap: 8px; }
    .opt-scenario-card:hover { background: var(--color-surface-alt); border-color: var(--color-text-dim); }
    .opt-scenario-card.selected { border-color: var(--color-accent); background: var(--color-surface-alt); box-shadow: 0 0 0 1px var(--color-accent); }
    .opt-scenario-card strong { font-size: 0.85rem; white-space: nowrap; }
    .conv-desc { font-size: 0.7rem; color: var(--color-text-dim); }

    .gopt-algo-group { margin-bottom: 4px; }
    .gopt-algo-label { font-size: 0.65rem; color: var(--color-text-dim); margin-bottom: 2px; }
    .gopt-algo-row { display: flex; gap: 3px; margin-bottom: 3px; }
    .gopt-algo-btn {
      flex: 1; padding: 4px 2px; border: 1px solid var(--color-border); border-radius: 5px;
      background: none; color: var(--color-text); font-size: 0.7rem; cursor: pointer;
      transition: all 0.12s; text-align: center; }
    .gopt-algo-btn:hover { background: var(--color-surface-alt); }
    .gopt-algo-btn.active { border-color: var(--color-accent); background: var(--color-surface-alt); font-weight: 600; box-shadow: 0 0 0 1px var(--color-accent); }

    .opt-run-btn {
      width: 100%; margin-top: 10px; padding: 10px; border: none; border-radius: 8px;
      font-size: 0.85rem; font-weight: 600; cursor: pointer;
      background: var(--color-accent); color: var(--color-accent-on); transition: background 0.15s; }
    .opt-run-btn:hover:not([disabled]) { background: var(--color-accent-hover); }
    .opt-run-btn[disabled] { opacity: 0.6; cursor: not-allowed; }
    .opt-stop-btn {
      width: 100%; margin-top: 6px; padding: 8px; border: 1px solid var(--color-error); border-radius: 8px;
      font-size: 0.8rem; font-weight: 600; cursor: pointer; background: none; color: var(--color-error); }

    .opt-progress { height: 6px; background: var(--color-progress-bg, #e0e4ea); border-radius: 3px; margin-top: 10px; overflow: hidden; }
    .opt-progress-bar { height: 100%; background: var(--color-accent); border-radius: 3px; transition: width 0.2s ease; }
    .opt-progress-text { font-size: 0.72rem; color: var(--color-text-dim); text-align: center; margin-top: 4px; }

    .opt-summary { margin-top: 10px; padding: 10px; background: var(--color-surface-alt); border-radius: 8px; border: 1px solid var(--color-border); }
    .opt-summary h3 { font-size: 0.8rem; margin-bottom: 6px; }
    .opt-summary table { width: 100%; font-size: 0.72rem; border-collapse: collapse; }
    .opt-summary td { padding: 2px 4px; }

    #graph-container { width: 100%; text-align: center; }
    #mol-vis { width: 100%; max-width: 400px; margin: 6px auto 0; }
    .opt-hint { color: var(--color-text-dim); font-size: 0.85rem; padding: 60px 20px; }
    .gopt-slider { width: 100%; max-width: 500px; margin: 6px auto 0; display: flex; align-items: center; gap: 10px; }
    .gopt-slider input[type=range] { flex: 1; cursor: pointer; }
    .gopt-slider-label { font-size: 0.78rem; font-weight: 600; min-width: 90px; text-align: center; }
    .gopt-step-info { font-size: 0.72rem; color: var(--color-text-dim); text-align: center; margin-top: 2px; }
    .gopt-xyz-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 8px; margin-top: 10px; }
    .gopt-xyz-block { background: var(--color-surface-alt); border: 1px solid var(--color-border); border-radius: 6px; padding: 6px 8px; overflow: hidden; }
    .gopt-xyz-label { font-size: 0.65rem; color: var(--color-text-dim); text-transform: uppercase; letter-spacing: 0.04em; margin-bottom: 3px; font-weight: 600; }
    .gopt-xyz-pre { font-family: ui-monospace, 'SF Mono', Consolas, monospace; font-size: 0.65rem; line-height: 1.35; color: var(--color-text); white-space: pre; overflow-x: auto; margin: 0; }
    @media (max-width: 700px) { .gopt-xyz-grid { grid-template-columns: 1fr; } }
    @media (max-width: 700px) { .opt-content { flex-direction: column; } .opt-controls { flex: none; } }
  `;
  document.head.appendChild(style);
}

// ── Main render ──────────────────────────────────────────────────────

const root = document.getElementById('app')!;

function render(): void {
  geomViewer = null; // DOM is about to be rebuilt
  root.innerHTML = `
    <div class="opt-page">
      ${renderHeader('geomopt')}
      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${t('gopt.molecule')}</h2>
          <div class="opt-scenario-grid" id="scen-grid"></div>

          <h2>${t('gopt.algorithm')}</h2>
          <div id="algo-area"></div>

          <div class="theory-row" style="margin:8px 0;display:flex;align-items:center;gap:8px;font-size:0.9rem;">
            <span>Theory:</span>${theorySelectHTML('theory-sel', theoryChoice, '', HEAVY_ITERATIVE_EXCLUDE)}
          </div>

          <div style="margin:6px 0 0;">
            <label style="font-size:0.72rem;">${t('gopt.maxIter')}
              <input type="number" id="max-iter" min="5" max="200" value="${maxIter}" style="width:50px;margin-left:4px;">
            </label>
          </div>

          <button id="run-btn" class="opt-run-btn" ${running ? 'disabled' : ''}>
            ${running ? t('gopt.running') : t('gopt.run')}
          </button>
          ${running ? `<button id="stop-btn" class="opt-stop-btn">${t('gopt.stop')}</button>` : ''}
          <div id="progress-area"></div>
        </div>
        <div class="opt-panel opt-graph-panel">
          <div id="graph-container">
            ${!done && !running && steps.length === 0 ? `<p class="opt-hint">${t('gopt.waiting')}</p>` : ''}
          </div>
          <div id="slider-area"></div>
          <div id="mol-vis"></div>
          <div id="summary-area"></div>
        </div>
      </div>
    </div>`;

  injectStyles();

  // Scenario cards
  const grid = root.querySelector('#scen-grid')!;
  for (const scen of SCENARIOS) {
    const card = document.createElement('div');
    card.className = 'opt-scenario-card' + (scen.id === selectedScen.id ? ' selected' : '');
    card.innerHTML = `<strong>${t(scen.labelKey)}</strong><span class="conv-desc">${t(scen.descKey)}</span>`;
    card.addEventListener('click', () => {
      if (running) return;
      selectedScen = scen; steps = []; done = false; converged = false; selectedIdx = 0; geomViewer = null;
      render();
    });
    grid.appendChild(card);
  }

  // Algorithm grouped buttons
  const algoArea = root.querySelector('#algo-area')!;
  for (const group of ALGO_GROUPS) {
    const div = document.createElement('div');
    div.className = 'gopt-algo-group';
    div.innerHTML = `<div class="gopt-algo-label">${group.label}</div><div class="gopt-algo-row">${
      group.items.map(it =>
        `<button class="gopt-algo-btn${it.id === selectedAlgo ? ' active' : ''}" data-algo="${it.id}">${it.name}</button>`
      ).join('')
    }</div>`;
    algoArea.appendChild(div);
  }
  for (const btn of root.querySelectorAll<HTMLButtonElement>('.gopt-algo-btn')) {
    btn.addEventListener('click', () => {
      if (running) return;
      selectedAlgo = btn.dataset.algo as AlgoId;
      for (const b of root.querySelectorAll('.gopt-algo-btn')) b.classList.remove('active');
      btn.classList.add('active');
    });
  }

  // Events
  root.querySelector('#nav-theme')!.addEventListener('click', () => { toggleTheme(); render(); });
  root.querySelector('#nav-lang')!.addEventListener('click', () => { toggleLang(); render(); });
  const theorySel = root.querySelector<HTMLSelectElement>('#theory-sel');
  if (theorySel) theorySel.addEventListener('change', () => { theoryChoice = theorySel.value as TheoryChoice; });
  root.querySelector('#run-btn')!.addEventListener('click', () => {
    if (!running) { maxIter = +(root.querySelector<HTMLInputElement>('#max-iter')?.value ?? 50); runOptimization(); }
  });
  root.querySelector('#stop-btn')?.addEventListener('click', () => { stopRequested = true; });

  if (steps.length > 0) { renderGraphs(); renderSummary(); setupSlider(); renderMolVis(); }
}

// ── Progress / Summary ───────────────────────────────────────────────

function updateProgress(step: number, maxF: number, detail?: string): void {
  const area = root.querySelector('#progress-area');
  if (!area) return;
  const text = detail ?? `Step ${step} | Max|F| = ${maxF.toExponential(2)}`;
  area.innerHTML = `
    <div class="opt-progress"><div class="opt-progress-bar" style="width:${Math.min(100,step/maxIter*100).toFixed(0)}%"></div></div>
    <div class="opt-progress-text">${text}</div>`;
}

function renderSummary(): void {
  const el = root.querySelector('#summary-area');
  if (!el || steps.length === 0) return;
  const first = steps[0], last = steps[steps.length - 1];
  const sc = converged ? 'var(--color-accent)' : 'var(--color-error)';
  const st = converged ? `\u2705 ${t('gopt.converged')}` : done ? `\u26a0\ufe0f ${t('gopt.notConverged')}` : t('gopt.running');
  const algoItem = ALGO_GROUPS.flatMap(g => g.items).find(it => it.id === selectedAlgo);
  const d = (k: string, v: string) => `<tr><td style="color:var(--color-text-dim)">${k}</td><td>${v}</td></tr>`;

  // Format XYZ block (Bohr → Angstrom)
  const BOHR_TO_ANG = 0.529177210903;
  const atomNums = parseXYZ(selectedScen.xyz).map(a => a.atomicNumber);
  const fmtXYZ = (coords: number[]): string => {
    const lines: string[] = [];
    for (let i = 0; i < atomNums.length; i++) {
      const sym = atomicNumberToElementName(atomNums[i]);
      const x = (coords[3*i]   * BOHR_TO_ANG).toFixed(6).padStart(11);
      const y = (coords[3*i+1] * BOHR_TO_ANG).toFixed(6).padStart(11);
      const z = (coords[3*i+2] * BOHR_TO_ANG).toFixed(6).padStart(11);
      lines.push(`${sym.padEnd(2)}${x}${y}${z}`);
    }
    return lines.join('\n');
  };

  el.innerHTML = `<div class="opt-summary"><h3 style="color:${sc}">${st}</h3><table>
    ${d('Theory', theoryDisplayLabel(theoryChoice))}
    ${d(t('gopt.algorithm'), algoItem?.name ?? selectedAlgo)}
    ${d(t('gopt.steps'), String(steps.length))}
    ${d(t('gopt.initEnergy'), first.energy.toFixed(8) + ' Eh')}
    ${d(t('gopt.finalEnergy'), last.energy.toFixed(8) + ' Eh')}
    ${d(t('gopt.energyChange'), (last.energy - first.energy).toFixed(8) + ' Eh')}
    ${d(t('gopt.maxForce'), last.maxForce.toExponential(4) + ' Eh/bohr')}
    ${d(t('gopt.threshold'), forceThreshold.toExponential(1) + ' Eh/bohr')}
    ${done ? d(t('gopt.time'), (elapsedMs / 1000).toFixed(1) + ' s') : ''}
  </table>
  <div class="gopt-xyz-grid">
    <div class="gopt-xyz-block">
      <div class="gopt-xyz-label">Initial (\u00c5)</div>
      <pre class="gopt-xyz-pre">${fmtXYZ(Array.from(first.coords))}</pre>
    </div>
    <div class="gopt-xyz-block">
      <div class="gopt-xyz-label">Final (\u00c5)</div>
      <pre class="gopt-xyz-pre">${fmtXYZ(Array.from(last.coords))}</pre>
    </div>
  </div>
  </div>`;
}

// ── Slider ───────────────────────────────────────────────────────────

function setupSlider(): void {
  const area = root.querySelector('#slider-area');
  if (!area || steps.length < 2) return;
  const s = steps[selectedIdx];
  area.innerHTML = `
    <div class="gopt-slider">
      <span class="gopt-slider-label">Step ${selectedIdx+1}/${steps.length}</span>
      <input type="range" id="step-slider" min="0" max="${steps.length-1}" value="${selectedIdx}" step="1">
    </div>
    <div class="gopt-step-info">E = ${s.energy.toFixed(8)} Eh &nbsp; Max|F| = ${s.maxForce.toExponential(3)}</div>`;
  root.querySelector<HTMLInputElement>('#step-slider')!.addEventListener('input', (ev) => {
    selectedIdx = +(ev.target as HTMLInputElement).value;
    renderGraphs(); renderMolVis();
    const ss = steps[selectedIdx];
    area.querySelector('.gopt-slider-label')!.textContent = `Step ${selectedIdx+1}/${steps.length}`;
    area.querySelector('.gopt-step-info')!.innerHTML = `E = ${ss.energy.toFixed(8)} Eh &nbsp; Max|F| = ${ss.maxForce.toExponential(3)}`;
  });
}

// ── Graphs ───────────────────────────────────────────────────────────

function renderGraphs(): void {
  const container = root.querySelector('#graph-container');
  if (!container || steps.length < 1) return;
  const W = 500, H = 155, PAD = { l: 70, r: 20, t: 24, b: 28 };
  const c = getThemeColors();

  function makeSVG(title: string, yLabel: string, values: number[], useLog: boolean): string {
    const nn = values.length, xMin = 1, xMax = Math.max(nn, 2);
    let yMin: number, yMax: number;
    if (useLog) {
      const lv = values.map(v => Math.log10(Math.max(v, 1e-15)));
      yMin = Math.min(...lv); yMax = Math.max(...lv);
      if (yMax - yMin < 0.5) { yMin -= 0.25; yMax += 0.25; }
    } else {
      yMin = Math.min(...values); yMax = Math.max(...values);
      const m = (yMax - yMin) * 0.1 || 0.001; yMin -= m; yMax += m;
    }
    const pw = W-PAD.l-PAD.r, ph = H-PAD.t-PAD.b;
    const sx = (x: number) => PAD.l+(x-xMin)/(xMax-xMin)*pw;
    const sy = (y: number) => PAD.t+(1-(y-yMin)/(yMax-yMin))*ph;

    let svg = `<svg viewBox="0 0 ${W} ${H}" style="width:100%;max-width:${W}px;background:${c.surface};border-radius:6px;margin-bottom:4px;">`;
    svg += `<text x="${W/2}" y="14" text-anchor="middle" font-size="10" fill="${c.titleSvg}" font-weight="600">${title}</text>`;
    svg += `<line x1="${PAD.l}" y1="${PAD.t}" x2="${PAD.l}" y2="${H-PAD.b}" stroke="${c.grid}"/>`;
    svg += `<line x1="${PAD.l}" y1="${H-PAD.b}" x2="${W-PAD.r}" y2="${H-PAD.b}" stroke="${c.grid}"/>`;
    svg += `<text x="12" y="${H/2}" text-anchor="middle" font-size="7" fill="${c.dim}" transform="rotate(-90,12,${H/2})">${yLabel}</text>`;
    svg += `<text x="${PAD.l+pw/2}" y="${H-3}" text-anchor="middle" font-size="7" fill="${c.dim}">${t('gopt.xStep')}</text>`;
    for (let i = 0; i <= 4; i++) {
      const yv = yMin+(yMax-yMin)*i/4, py = sy(yv);
      svg += `<line x1="${PAD.l-3}" y1="${py}" x2="${W-PAD.r}" y2="${py}" stroke="${c.grid}" stroke-dasharray="2,3"/>`;
      svg += `<text x="${PAD.l-5}" y="${py+3}" text-anchor="end" font-size="6.5" fill="${c.dim}">${useLog?`1e${yv.toFixed(1)}`:yv.toPrecision(6)}</text>`;
    }
    if (useLog) {
      const thY = sy(Math.log10(forceThreshold));
      if (thY >= PAD.t && thY <= H-PAD.b)
        svg += `<line x1="${PAD.l}" y1="${thY}" x2="${W-PAD.r}" y2="${thY}" stroke="#ef4444" stroke-dasharray="4,3"/><text x="${W-PAD.r-2}" y="${thY-3}" text-anchor="end" font-size="7" fill="#ef4444">threshold</text>`;
    }
    if (nn >= 1) {
      const pts = values.map((v, i) => `${sx(i+1).toFixed(1)},${sy(useLog?Math.log10(Math.max(v,1e-15)):v).toFixed(1)}`);
      svg += `<polyline points="${pts.join(' ')}" fill="none" stroke="${c.occupied}" stroke-width="1.5"/>`;
      for (let i = 0; i < nn; i++) {
        const yv = useLog?Math.log10(Math.max(values[i],1e-15)):values[i];
        const sel = i === selectedIdx;
        svg += `<circle cx="${sx(i+1).toFixed(1)}" cy="${sy(yv).toFixed(1)}" r="${sel?5:2.5}" fill="${sel?c.accent:c.occupied}" ${sel?`stroke="${c.titleSvg}" stroke-width="1.5"`:''} />`;
      }
    }
    svg += '</svg>';
    return svg;
  }
  container.innerHTML =
    makeSVG(t('gopt.graphEnergy'), t('gopt.yEnergy'), steps.map(s => s.energy), false) +
    makeSVG(t('gopt.graphGrad'), t('gopt.yForce'), steps.map(s => s.maxForce), true);
}

// ── 3D Molecule visualization with gradient arrows ───────────────────

function renderMolVis(): void {
  const el = root.querySelector<HTMLDivElement>('#mol-vis');
  if (!el || steps.length === 0) return;
  const step = steps[selectedIdx];
  const atomNums = parseXYZ(selectedScen.xyz).map(a => a.atomicNumber);
  if (geomViewer) {
    geomViewer.update(step.coords, step.gradient, forceThreshold);
  } else {
    geomViewer = createGeomOptViewer(el, atomNums, step.coords, step.gradient, forceThreshold);
  }
}

// ── Init ─────────────────────────────────────────────────────────────

initTheme(); initLang(); render();
