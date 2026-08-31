/**
 * Cross-device performance benchmark for the paper.
 *
 * Spec:
 * - 4 fixed workloads (H2O HF/6-31G(d,p), H2O geom opt, CO2 Hessian, H2O B3LYP/6-31G(d,p))
 * - 4 runs per workload: 1 warmup discarded, 3 timed → median
 * - Auto-select best available backend (WASM-SIMD > WASM > JS)
 * - Records device info (UA, RAM, cores), backend, initial load time
 * - JSON + CSV export (clipboard + download)
 */

import { parseXYZ } from './core/parseXYZ';
import { BasisSet } from './core/basisSet';
import { Molecular } from './core/molecular';
import { buildHF, type DFTConfig } from './core/builder';
import { RHF } from './core/rhf';
import { computeRHFGradient } from './core/gradient';
import { computeHessianAuto } from './core/hessian';
import { initWasm, isWasmAvailable, getActiveBackend } from './core/eriWasm';
import { computeERIs } from './core/integrals2e';
import { ERIStored } from './core/eri';
import { Matrix } from './linalg/matrix';
import type { EriBackend } from './core/hf';
import { gatherDeviceInfo, getInitialLoadMs, type DeviceInfo } from './core/deviceInfo';
import { stressPanelHTML, wireStressPanel } from './ui/stressPanel';

// ── Constants ──
const RUNS_TOTAL = 4;
const RUNS_DISCARD = 1;
const CONV_TOL = 1e-7;
const ANG_TO_BOHR = 1.8897259886;

// ── Fixed test geometries ──
const H2O_XYZ = `3
H2O standard
O   0.000000   0.000000   0.117790
H   0.000000   0.755450  -0.471160
H   0.000000  -0.755450  -0.471160`;

const H2O_OFF_XYZ = `3
H2O displaced
O   0.000000   0.000000   0.000000
H   0.000000   0.850000   0.620000
H   0.000000  -0.850000   0.620000`;

const CO2_XYZ = `3
CO2 linear
C   0.000000   0.000000   0.000000
O   0.000000   0.000000   1.160000
O   0.000000   0.000000  -1.160000`;

// ── Basis cache ──
const basisCache = new Map<string, BasisSet>();
async function loadBasis(name: string): Promise<BasisSet> {
  if (basisCache.has(name)) return basisCache.get(name)!;
  const url = `${import.meta.env.BASE_URL}basis/${name}.gbs`;
  const r = await fetch(url);
  if (!r.ok) throw new Error(`Failed to load basis ${name}: ${r.status}`);
  const bs = BasisSet.fromGBS(await r.text());
  basisCache.set(name, bs);
  return bs;
}

// ── Workloads ──
interface Workload {
  id: string;
  label: string;
  run: (backend: EriBackend) => Promise<void>;
}

/** H2O HF/6-31G(d,p) single point — basic SCF */
async function workH2O_HF_631gdp(backend: EriBackend): Promise<void> {
  const basis = await loadBasis('6-31g(d,p)');
  const atoms = parseXYZ(H2O_XYZ);
  const mol = new Molecular(atoms, basis, 0);
  const hf = buildHF(mol, 'RHF');
  await hf.solve({ eriBackend: backend });
}

/** H2O geom opt HF/STO-3G — iterative SCF from displaced start */
async function workH2O_GeomOpt(backend: EriBackend): Promise<void> {
  const basis = await loadBasis('sto-3g');
  const atoms = parseXYZ(H2O_OFF_XYZ);
  let coords = new Float64Array(atoms.length * 3);
  for (let i = 0; i < atoms.length; i++) {
    coords[3*i]   = atoms[i].coordinate.x;
    coords[3*i+1] = atoms[i].coordinate.y;
    coords[3*i+2] = atoms[i].coordinate.z;
  }
  const STEP = 0.4;
  const MAX_ITER = 25;
  for (let iter = 0; iter < MAX_ITER; iter++) {
    const at = atoms.map((a, i) => ({
      ...a,
      coordinate: { x: coords[3*i], y: coords[3*i+1], z: coords[3*i+2] },
    }));
    const mol = new Molecular(at, basis, 0);
    const hf = buildHF(mol, 'RHF') as RHF;
    await hf.solve({ eriBackend: backend });
    const grad = computeRHFGradient(
      mol.primitiveShells, mol.atoms, mol.cgtoNormalizationFactors,
      mol.numBasis, mol.numAlphaSpins, hf.density, hf.coefficients, hf.orbitalEnergies,
    ).total;
    let maxF = 0;
    for (let i = 0; i < grad.length; i++) if (Math.abs(grad[i]) > maxF) maxF = Math.abs(grad[i]);
    if (maxF < 5e-4) break;
    for (let i = 0; i < grad.length; i++) coords[i] -= STEP * grad[i];
  }
}

/** CO2 Hessian HF/STO-3G — 6N grad evals (heaviest derivative task) */
async function workCO2_Hess(backend: EriBackend): Promise<void> {
  const basis = await loadBasis('sto-3g');
  const atoms = parseXYZ(CO2_XYZ);
  const refCoords = new Float64Array(atoms.length * 3);
  for (let i = 0; i < atoms.length; i++) {
    refCoords[3*i]   = atoms[i].coordinate.x * ANG_TO_BOHR;
    refCoords[3*i+1] = atoms[i].coordinate.y * ANG_TO_BOHR;
    refCoords[3*i+2] = atoms[i].coordinate.z * ANG_TO_BOHR;
  }
  await computeHessianAuto(
    atoms.map(a => a.atomicNumber), refCoords, basis, 0, 5e-4,
    undefined, undefined,
  );
}

/** H2O B3LYP/6-31G(d,p) single point — DFT reality check */
async function workH2O_B3LYP_631gdp(backend: EriBackend): Promise<void> {
  const basis = await loadBasis('6-31g(d,p)');
  const atoms = parseXYZ(H2O_XYZ);
  const mol = new Molecular(atoms, basis, 0);
  const dft: DFTConfig = { functional: 'B3LYP', gridLevel: 'medium' };
  const hf = buildHF(mol, 'RHF', dft);
  await hf.solve({ eriBackend: backend });
}

const WORKLOADS: Workload[] = [
  { id: 'h2o_hf_631gdp',      label: 'H2O HF/6-31G(d,p)',       run: workH2O_HF_631gdp },
  { id: 'h2o_geomopt_sto3g',  label: 'H2O geom opt HF/STO-3G',  run: workH2O_GeomOpt },
  { id: 'co2_hess_sto3g',     label: 'CO2 Hessian HF/STO-3G',   run: workCO2_Hess },
  { id: 'h2o_b3lyp_631gdp',   label: 'H2O B3LYP/6-31G(d,p)',    run: workH2O_B3LYP_631gdp },
];

// ── Median helper ──
function median(xs: number[]): number {
  const s = [...xs].sort((a, b) => a - b);
  const m = s.length;
  if (m === 0) return NaN;
  return m % 2 === 1 ? s[(m - 1) / 2] : (s[m/2 - 1] + s[m/2]) / 2;
}

// ── Device info: see ./core/deviceInfo (shared with the stress panel) ──

// ── Direct WASM binary loader (for backend comparison only).
//    Loads a specific wasm_eri binary and exposes its ERI + Fock functions.
//    Used to compare JS / WASM-baseline / WASM-SIMD on the same device. ──

interface DirectWasm {
  computeERIs: (shellsFlat: Float64Array, normFactors: Float64Array, n: number, threshold: number) => Float64Array;
  computeFockRhf: (eri: Float64Array, dens: Float64Array, coreH: Float64Array, n: number) => Float64Array;
}

async function loadWasmBinaryDirect(url: string): Promise<DirectWasm> {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`Cannot fetch ${url}`);
  const bytes = await r.arrayBuffer();
  const mod = await WebAssembly.compile(bytes);
  let instance: WebAssembly.Instance;
  const imports: WebAssembly.Imports = {
    './wasm_eri_bg.js': {
      __wbindgen_init_externref_table: () => {
        const table = (instance.exports as unknown as { __wbindgen_externrefs: WebAssembly.Table }).__wbindgen_externrefs;
        const offset = table.grow(4);
        table.set(0, undefined);
        table.set(offset + 0, undefined);
        table.set(offset + 1, null);
        table.set(offset + 2, true);
        table.set(offset + 3, false);
      },
    },
  };
  instance = await WebAssembly.instantiate(mod, imports);
  const exp = instance.exports as unknown as Record<string, Function> & { memory: WebAssembly.Memory };
  if (exp.__wbindgen_start) exp.__wbindgen_start();
  const memory = exp.memory;
  const malloc = exp.__wbindgen_malloc as (n: number, a: number) => number;
  const free = exp.__wbindgen_free as (p: number, n: number, a: number) => void;
  const writeF64 = (data: ArrayLike<number>): [number, number] => {
    const ptr = malloc(data.length * 8, 8) >>> 0;
    new Float64Array(memory.buffer).set(data as never, ptr / 8);
    return [ptr, data.length];
  };
  const readF64 = (result: [number, number]): Float64Array => {
    const [ptr, len] = result;
    const out = new Float64Array(memory.buffer).slice(ptr / 8, ptr / 8 + len);
    free(ptr, len * 8, 8);
    return out;
  };
  return {
    computeERIs: (shellsFlat, normFactors, n, threshold) => {
      const [sp, sl] = writeF64(shellsFlat);
      const [np, nl] = writeF64(normFactors);
      return readF64((exp.compute_eris_wasm as Function)(sp, sl, np, nl, n, threshold));
    },
    computeFockRhf: (eri, dens, coreH, n) => {
      const [ep, el] = writeF64(eri);
      const [dp, dl] = writeF64(dens);
      const [cp, cl] = writeF64(coreH);
      return readF64((exp.compute_fock_rhf as Function)(ep, el, dp, dl, cp, cl, n));
    },
  };
}

// Lazily-loaded WASM binaries for backend comparison
let _wasmBase: DirectWasm | null | undefined;  // undefined = not tried, null = failed
let _wasmSimd: DirectWasm | null | undefined;
async function getWasmBase(): Promise<DirectWasm | null> {
  if (_wasmBase !== undefined) return _wasmBase;
  try { _wasmBase = await loadWasmBinaryDirect(`${import.meta.env.BASE_URL}wasm/wasm_eri_bg.wasm`); }
  catch { _wasmBase = null; }
  return _wasmBase;
}
async function getWasmSimd(): Promise<DirectWasm | null> {
  if (_wasmSimd !== undefined) return _wasmSimd;
  try { _wasmSimd = await loadWasmBinaryDirect(`${import.meta.env.BASE_URL}wasm/wasm_eri_simd_bg.wasm`); }
  catch { _wasmSimd = null; }
  return _wasmSimd;
}

// ── Backend selection ──
async function selectBestBackend(): Promise<{ backend: EriBackend; label: string }> {
  // initWasm needs the deploy base URL to fetch /GANSU-Lite/wasm/... on GitHub Pages.
  await initWasm(import.meta.env.BASE_URL);
  if (isWasmAvailable()) {
    const active = getActiveBackend();
    if (active === 'wasm-simd') return { backend: 'wasm', label: 'WASM-SIMD' };
    return { backend: 'wasm', label: 'WASM' };
  }
  return { backend: 'js', label: 'JS' };
}

// ── Pure-JS Fock matrix construction (RHF), 4-index loop ──
//     Used to measure the TypeScript baseline in the backend comparison.
function eri4Idx(i: number, j: number, k: number, l: number): number {
  let ii = i, jj = j; if (ii < jj) { const t = ii; ii = jj; jj = t; }
  let kk = k, ll = l; if (kk < ll) { const t = kk; kk = ll; ll = t; }
  let ij = ii * (ii + 1) / 2 + jj, kl = kk * (kk + 1) / 2 + ll;
  if (ij < kl) { const t = ij; ij = kl; kl = t; }
  return ij * (ij + 1) / 2 + kl;
}
function computeFockJS(eriData: Float64Array, density: Float64Array, coreH: Float64Array, n: number): Float64Array {
  const F = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = i; j < n; j++) {
      let G = 0;
      for (let k = 0; k < n; k++) {
        for (let l = 0; l < n; l++) {
          const P = density[k * n + l];
          if (Math.abs(P) < 1e-10) continue;
          G += P * (eriData[eri4Idx(i, j, k, l)] - 0.5 * eriData[eri4Idx(i, k, j, l)]);
        }
      }
      const v = coreH[i * n + j] + G;
      F[i * n + j] = v;
      F[j * n + i] = v;
    }
  }
  return F;
}

// ── Manual SCF loop for backend comparison ──
//     Reuses HF infrastructure for setup but invokes the chosen Fock function
//     directly. Returns ERI time and SCF-loop time separately. */
type FockFn = (eri: Float64Array, density: Float64Array, coreH: Float64Array, n: number) => Float64Array | Promise<Float64Array>;

interface ManualScfResult { eriMs: number; scfMs: number; iters: number; energy: number; }

async function runManualScf(
  mol: Molecular, eriData: Float64Array, eriPreparMs: number, fockFn: FockFn,
): Promise<ManualScfResult> {
  const n = mol.numBasis;
  const eri = new ERIStored(eriData, n);
  const hf = buildHF(mol, 'RHF') as RHF;
  // Bypass eri precompute by injecting; setup integrals + transform manually.
  (hf as unknown as { computeNuclearRepulsionEnergy: () => number }).computeNuclearRepulsionEnergy();
  (hf as unknown as { computeCoreHamiltonianMatrix: () => void }).computeCoreHamiltonianMatrix();
  (hf as unknown as { computeTransformMatrix: () => void }).computeTransformMatrix();
  (hf as unknown as { eri: ERIStored }).eri = eri;
  (hf as unknown as { guessInitialFockMatrix: () => void }).guessInitialFockMatrix();
  const coreHData = ((hf as unknown as { coreHamiltonianMatrix: Matrix }).coreHamiltonianMatrix).data as Float64Array;

  let prevE = 0, iters = 0, energy = 0;
  const t0 = performance.now();
  for (let iter = 0; iter < 200; iter++) {
    (hf as unknown as { computeCoefficientMatrix: () => void }).computeCoefficientMatrix();
    (hf as unknown as { computeDensityMatrix: () => void }).computeDensityMatrix();
    const densityData = ((hf as unknown as { densityMatrix: Matrix }).densityMatrix).data as Float64Array;
    const fockData = await fockFn(eriData, densityData, coreHData, n);
    (hf as unknown as { fockMatrix: Matrix }).fockMatrix = new Matrix(n, n, fockData);
    energy = (hf as unknown as { computeEnergy: () => number }).computeEnergy();
    if (iter > 0 && Math.abs(energy - prevE) < 1e-8) { iters = iter + 1; break; }
    (hf as unknown as { updateFockMatrix: () => void }).updateFockMatrix();
    prevE = energy;
  }
  return { eriMs: eriPreparMs, scfMs: performance.now() - t0, iters, energy };
}

// ── Backend comparison cases ──
interface CompareCase { molecule: string; xyzFile: string; basis: string; }
const COMPARE_CASES: CompareCase[] = [
  { molecule: 'H2O',     xyzFile: 'H2O.xyz',     basis: '6-31g(d,p)' },
  { molecule: 'CO2',     xyzFile: 'CO2.xyz',     basis: '6-31g(d,p)' },
  { molecule: 'Benzene', xyzFile: 'Benzene.xyz', basis: '6-31g(d,p)' },
];
const COMPARE_RUNS = 3;  // median of 3, no warmup discard (matches paper spec)

interface BackendTiming {
  eriMedianMs: number;
  scfMedianMs: number;
  totalMedianMs: number;
  iters: number;
  eriRuns: number[];
  scfRuns: number[];
  error?: string;
}
interface CompareRow {
  molecule: string;
  basis: string;
  nbasis: number;
  js: BackendTiming;
  wasmBase: BackendTiming;
  wasmSimd: BackendTiming;
}

async function loadCaseMolecule(c: CompareCase): Promise<{ mol: Molecular; xyz: string }> {
  const basis = await loadBasis(c.basis);
  const xyzUrl = `${import.meta.env.BASE_URL}xyz/${c.xyzFile}`;
  const r = await fetch(xyzUrl);
  if (!r.ok) throw new Error(`Cannot load ${xyzUrl}`);
  const xyz = await r.text();
  const atoms = parseXYZ(xyz);
  return { mol: new Molecular(atoms, basis, 0), xyz };
}

async function runOneCompareCase(c: CompareCase, base: DirectWasm | null, simd: DirectWasm | null, status: (s: string) => void): Promise<CompareRow> {
  const { mol } = await loadCaseMolecule(c);
  const n = mol.numBasis;
  const shells = mol.primitiveShells;
  const norms = mol.cgtoNormalizationFactors;
  const shellsFlat = (await import('./core/eriWasm')).packShells(shells);
  const normsF64 = new Float64Array(norms);

  async function bench(label: string, eriFn: () => Promise<Float64Array> | Float64Array, fockFn: FockFn): Promise<BackendTiming> {
    const eriTimes: number[] = [], scfTimes: number[] = [];
    let iters = 0;
    let error: string | undefined;
    try {
      for (let i = 0; i < COMPARE_RUNS; i++) {
        status(`  ${label} run ${i + 1}/${COMPARE_RUNS}…`);
        const t0 = performance.now();
        const eriData = await eriFn();
        const eriMs = performance.now() - t0;
        const r = await runManualScf(mol, eriData, eriMs, fockFn);
        eriTimes.push(r.eriMs);
        scfTimes.push(r.scfMs);
        if (i === 0) iters = r.iters;
        await new Promise(res => setTimeout(res, 30));
      }
    } catch (e) {
      error = (e instanceof Error) ? e.message : String(e);
    }
    const em = eriTimes.length ? median(eriTimes) : NaN;
    const sm = scfTimes.length ? median(scfTimes) : NaN;
    return {
      eriMedianMs: em, scfMedianMs: sm, totalMedianMs: em + sm,
      iters, eriRuns: eriTimes, scfRuns: scfTimes,
      ...(error ? { error } : {}),
    };
  }

  const jsT = await bench('JS',
    () => computeERIs(shells, norms, n, 1e-10, undefined, 'js'),
    (e, d, c, nn) => computeFockJS(e, d, c, nn));
  const baseT = base
    ? await bench('WASM baseline',
        () => base.computeERIs(shellsFlat, normsF64, n, 1e-10),
        (e, d, c, nn) => base.computeFockRhf(e, d, c, nn))
    : { eriMedianMs: NaN, scfMedianMs: NaN, totalMedianMs: NaN, iters: 0, eriRuns: [], scfRuns: [], error: 'WASM baseline not available' };
  const simdT = simd
    ? await bench('WASM+SIMD',
        () => simd.computeERIs(shellsFlat, normsF64, n, 1e-10),
        (e, d, c, nn) => simd.computeFockRhf(e, d, c, nn))
    : { eriMedianMs: NaN, scfMedianMs: NaN, totalMedianMs: NaN, iters: 0, eriRuns: [], scfRuns: [], error: 'SIMD not available' };

  return { molecule: c.molecule, basis: c.basis, nbasis: n, js: jsT, wasmBase: baseT, wasmSimd: simdT };
}

let compareResults: CompareRow[] | null = null;

async function runBackendCompare(): Promise<void> {
  const btn = document.getElementById('cb-run-compare') as HTMLButtonElement | null;
  const stopBtn = document.getElementById('cb-stop-compare') as HTMLButtonElement | null;
  if (btn) btn.disabled = true;
  if (stopBtn) stopBtn.disabled = false;
  const status = (s: string) => { const el = document.getElementById('cb-compare-status'); if (el) el.textContent = s; };
  status('Loading WASM binaries…');
  const base = await getWasmBase();
  const simd = await getWasmSimd();
  const rows: CompareRow[] = [];
  for (let i = 0; i < COMPARE_CASES.length; i++) {
    const c = COMPARE_CASES[i];
    status(`[${i + 1}/${COMPARE_CASES.length}] ${c.molecule}/${c.basis}…`);
    try {
      rows.push(await runOneCompareCase(c, base, simd, status));
      renderCompareTable(rows);
    } catch (e) {
      const msg = (e instanceof Error) ? e.message : String(e);
      status(`  ERROR: ${msg}`);
    }
  }
  compareResults = rows;
  status('Done.');
  if (btn) btn.disabled = false;
  if (stopBtn) stopBtn.disabled = true;
}

function renderCompareTable(rows: CompareRow[]): void {
  const section = document.getElementById('cb-compare-section') as HTMLElement | null;
  if (!section) return;
  section.style.display = '';
  const tbody = document.getElementById('cb-compare-body') as HTMLTableSectionElement | null;
  if (!tbody) return;
  tbody.innerHTML = '';
  const fmt = (n: number) => isFinite(n) ? n.toFixed(0) : '—';
  const sp = (base: number, fast: number) => (isFinite(base) && isFinite(fast) && fast > 0) ? `${(base / fast).toFixed(1)}×` : '—';
  for (const r of rows) {
    const tr = document.createElement('tr');
    tr.innerHTML = `
      <td>${r.molecule}</td><td>${r.basis}</td><td class="num">${r.nbasis}</td><td class="num">${r.js.iters || '—'}</td>
      <td class="num">${fmt(r.js.eriMedianMs)}</td>
      <td class="num">${fmt(r.js.scfMedianMs)}</td>
      <td class="num"><b>${fmt(r.js.totalMedianMs)}</b></td>
      <td class="num">${fmt(r.wasmBase.eriMedianMs)}</td>
      <td class="num">${fmt(r.wasmBase.scfMedianMs)}</td>
      <td class="num"><b>${fmt(r.wasmBase.totalMedianMs)}</b></td>
      <td class="num">${sp(r.js.totalMedianMs, r.wasmBase.totalMedianMs)}</td>
      <td class="num">${fmt(r.wasmSimd.eriMedianMs)}</td>
      <td class="num">${fmt(r.wasmSimd.scfMedianMs)}</td>
      <td class="num"><b>${fmt(r.wasmSimd.totalMedianMs)}</b></td>
      <td class="num">${sp(r.js.totalMedianMs, r.wasmSimd.totalMedianMs)}</td>
    `;
    tbody.appendChild(tr);
  }
}

// ── Result type ──
interface WorkloadResult {
  id: string;
  label: string;
  warmupMs: number;
  runs: number[];      // 3 timed runs (after discard)
  medianMs: number;
  error?: string;
}

interface BenchOutput {
  device: DeviceInfo;
  backend: string;
  initialLoadMs: number | null;
  workloads: WorkloadResult[];
  toolUrl: string;
}

// ── UI ──
const root = document.getElementById('app')!;
let allResults: BenchOutput | null = null;
let detectedBackendLabel = 'detecting…';

function render(): void {
  const dev = gatherDeviceInfo();
  const loadMs = getInitialLoadMs();
  root.innerHTML = `
    <div class="cb-page">
      <header class="cb-header">
        <h1>GANSU Lite — Cross-device Benchmark</h1>
        <p class="cb-subtitle">Reproducible end-to-end timing for the paper. Median of 3 timed runs (+ 1 warmup discarded).</p>
      </header>

      <section class="cb-panel">
        <h2>Device & Backend</h2>
        <table class="cb-info">
          <tr><th>User-Agent</th><td>${dev.ua}</td></tr>
          <tr><th>Platform</th><td>${dev.platform}</td></tr>
          <tr><th>Device memory</th><td>${dev.deviceMemoryGB ? `${dev.deviceMemoryGB} GB` : 'unknown'}</td></tr>
          <tr><th>Cores (navigator)</th><td>${dev.hardwareConcurrency}</td></tr>
          <tr><th>Screen</th><td>${dev.screen}</td></tr>
          <tr><th>Initial load</th><td>${loadMs != null ? `${loadMs} ms` : 'unavailable'}</td></tr>
          <tr><th>Backend</th><td id="cb-backend">detecting…</td></tr>
        </table>
      </section>

      <section class="cb-panel">
        <h2>Workloads</h2>
        <p class="cb-note">All fixed: SCF tol = 1e-7, DIIS, default initial guess. Single Web Worker thread (no parallelism — single-core perf + memory bandwidth).</p>
        <ul class="cb-workloads">
          ${WORKLOADS.map(w => `<li id="cb-row-${w.id}"><span class="cb-wl-label">${w.label}</span><span class="cb-wl-status">pending</span></li>`).join('')}
        </ul>
        <button id="cb-run">Run all workloads</button>
        <button id="cb-stop" disabled>Stop</button>
      </section>

      <section class="cb-panel">
        <h2>Backend comparison (reference device)</h2>
        <p class="cb-note">
          SCF wall-clock time for H₂O, CO₂, and benzene at 6-31G(d,p), broken down by kernel
          (ERI computation + SCF iteration loop), comparing TypeScript / baseline WebAssembly /
          WebAssembly +simd128 on this device. Medians over 3 runs each.
        </p>
        <div style="margin:6px 0 8px">
          <button id="cb-run-compare">Run backend comparison</button>
          <button id="cb-stop-compare" disabled>Stop</button>
          <span id="cb-compare-status" style="margin-left:10px;font-size:0.85rem;color:var(--color-text-dim);font-family:'Cascadia Code',monospace;"></span>
        </div>
        <div id="cb-compare-section" style="display:none;overflow-x:auto;">
          <table class="cb-results-table" style="font-size:0.78rem;">
            <thead><tr>
              <th rowspan="2">Molecule</th><th rowspan="2">Basis</th>
              <th rowspan="2" class="num">N</th><th rowspan="2" class="num">Iters</th>
              <th colspan="3" style="text-align:center;border-bottom:1px solid var(--color-border)">TypeScript (JS)</th>
              <th colspan="4" style="text-align:center;border-bottom:1px solid var(--color-border)">Baseline WASM</th>
              <th colspan="4" style="text-align:center;border-bottom:1px solid var(--color-border)">WASM +simd128</th>
            </tr><tr>
              <th class="num">ERI</th><th class="num">SCF</th><th class="num">Total</th>
              <th class="num">ERI</th><th class="num">SCF</th><th class="num">Total</th><th class="num">↑vs JS</th>
              <th class="num">ERI</th><th class="num">SCF</th><th class="num">Total</th><th class="num">↑vs JS</th>
            </tr></thead>
            <tbody id="cb-compare-body"></tbody>
          </table>
          <p class="cb-note">All times in ms. ↑vs JS = TypeScript total / WASM total.</p>
          <div style="margin-top:8px">
            <button id="cb-copy-compare-json">Copy JSON</button>
            <button id="cb-download-compare-json">Download JSON</button>
            <button id="cb-download-compare-csv">Download CSV</button>
          </div>
        </div>
      </section>

      ${stressPanelHTML()}

      <section class="cb-panel" id="cb-results" style="display:none">
        <h2>Results</h2>
        <table class="cb-results-table">
          <thead><tr>
            <th>Workload</th>
            <th class="num">Warmup (ms)</th>
            <th class="num">Run 1 (ms)</th>
            <th class="num">Run 2 (ms)</th>
            <th class="num">Run 3 (ms)</th>
            <th class="num">Median (ms)</th>
          </tr></thead>
          <tbody id="cb-results-body"></tbody>
        </table>
        <div class="cb-export">
          <button id="cb-copy-json">Copy JSON</button>
          <button id="cb-download-json">Download JSON</button>
          <button id="cb-download-csv">Download CSV</button>
        </div>
        <details class="cb-json">
          <summary>Raw JSON</summary>
          <pre id="cb-json-pre"></pre>
        </details>
      </section>
    </div>
  `;
  injectStyles();

  document.getElementById('cb-run')!.addEventListener('click', runAll);
  document.getElementById('cb-stop')!.addEventListener('click', () => { stopRequested = true; });
  document.getElementById('cb-run-compare')?.addEventListener('click', runBackendCompare);
  document.getElementById('cb-copy-compare-json')?.addEventListener('click', () => copyText(JSON.stringify(compareResults, null, 2)));
  document.getElementById('cb-download-compare-json')?.addEventListener('click', () => download('backend-comparison.json', JSON.stringify(compareResults, null, 2), 'application/json'));
  document.getElementById('cb-download-compare-csv')?.addEventListener('click', () => download('backend-comparison.csv', compareToCSV(compareResults), 'text/csv'));
  document.getElementById('cb-copy-json')?.addEventListener('click', copyJson);
  document.getElementById('cb-download-json')?.addEventListener('click', () => download('benchmark.json', JSON.stringify(allResults, null, 2), 'application/json'));
  document.getElementById('cb-download-csv')?.addEventListener('click', () => download('benchmark.csv', toCSV(allResults), 'text/csv'));

  wireStressPanel(() => detectedBackendLabel);

  selectBestBackend().then(({ label }) => {
    detectedBackendLabel = label;
    const el = document.getElementById('cb-backend');
    if (el) el.textContent = label;
  });
}

let stopRequested = false;

function setStatus(id: string, text: string): void {
  const row = document.getElementById(`cb-row-${id}`);
  if (row) {
    const status = row.querySelector('.cb-wl-status');
    if (status) status.textContent = text;
  }
}

async function runAll(): Promise<void> {
  stopRequested = false;
  const runBtn = document.getElementById('cb-run') as HTMLButtonElement;
  const stopBtn = document.getElementById('cb-stop') as HTMLButtonElement;
  runBtn.disabled = true;
  stopBtn.disabled = false;

  const { backend, label: backendLabel } = await selectBestBackend();
  const device = gatherDeviceInfo();
  const initialLoadMs = getInitialLoadMs();
  const workloadResults: WorkloadResult[] = [];

  for (const w of WORKLOADS) {
    if (stopRequested) { setStatus(w.id, 'stopped'); continue; }
    const times: number[] = [];
    let errMsg: string | undefined;
    try {
      for (let i = 0; i < RUNS_TOTAL; i++) {
        if (stopRequested) break;
        setStatus(w.id, `run ${i + 1}/${RUNS_TOTAL}…`);
        const t0 = performance.now();
        await w.run(backend);
        const dt = performance.now() - t0;
        times.push(dt);
        setStatus(w.id, `run ${i + 1}/${RUNS_TOTAL}: ${dt.toFixed(0)} ms`);
        // Yield to UI between runs
        await new Promise(r => setTimeout(r, 50));
      }
    } catch (e) {
      errMsg = (e instanceof Error) ? e.message : String(e);
      setStatus(w.id, `ERROR: ${errMsg}`);
    }
    const measured = times.slice(RUNS_DISCARD);
    const med = measured.length > 0 ? median(measured) : NaN;
    workloadResults.push({
      id: w.id, label: w.label,
      warmupMs: times[0] ?? NaN,
      runs: measured,
      medianMs: med,
      ...(errMsg ? { error: errMsg } : {}),
    });
    if (!errMsg && measured.length > 0) {
      setStatus(w.id, `median ${med.toFixed(0)} ms`);
    }
  }

  allResults = {
    device,
    backend: backendLabel,
    initialLoadMs,
    workloads: workloadResults,
    toolUrl: location.href,
  };

  renderResults(allResults);
  runBtn.disabled = false;
  stopBtn.disabled = true;
}

function renderResults(out: BenchOutput): void {
  const section = document.getElementById('cb-results') as HTMLElement;
  section.style.display = '';
  const body = document.getElementById('cb-results-body') as HTMLTableSectionElement;
  body.innerHTML = '';
  for (const w of out.workloads) {
    const tr = document.createElement('tr');
    tr.innerHTML = `
      <td>${w.label}${w.error ? ` <span class="cb-err">(${w.error})</span>` : ''}</td>
      <td class="num">${isFinite(w.warmupMs) ? w.warmupMs.toFixed(0) : '—'}</td>
      <td class="num">${w.runs[0] != null ? w.runs[0].toFixed(0) : '—'}</td>
      <td class="num">${w.runs[1] != null ? w.runs[1].toFixed(0) : '—'}</td>
      <td class="num">${w.runs[2] != null ? w.runs[2].toFixed(0) : '—'}</td>
      <td class="num"><b>${isFinite(w.medianMs) ? w.medianMs.toFixed(0) : '—'}</b></td>
    `;
    body.appendChild(tr);
  }
  (document.getElementById('cb-json-pre') as HTMLPreElement).textContent = JSON.stringify(out, null, 2);
}

function copyJson(): void {
  if (!allResults) return;
  copyText(JSON.stringify(allResults, null, 2));
}

function copyText(s: string): void {
  navigator.clipboard.writeText(s)
    .then(() => alert('Copied to clipboard'))
    .catch(() => alert('Copy failed — open browser console'));
}

function compareToCSV(rows: CompareRow[] | null): string {
  if (!rows || rows.length === 0) return '';
  const head = [
    'molecule', 'basis', 'nbasis', 'iters',
    'js_eri_ms', 'js_scf_ms', 'js_total_ms',
    'wasm_eri_ms', 'wasm_scf_ms', 'wasm_total_ms', 'speedup_wasm_vs_js',
    'simd_eri_ms', 'simd_scf_ms', 'simd_total_ms', 'speedup_simd_vs_js',
  ];
  const fmt = (n: number) => isFinite(n) ? n.toFixed(2) : '';
  const sp = (b: number, f: number) => isFinite(b) && isFinite(f) && f > 0 ? (b / f).toFixed(2) : '';
  const data = rows.map(r => [
    r.molecule, r.basis, String(r.nbasis), String(r.js.iters),
    fmt(r.js.eriMedianMs), fmt(r.js.scfMedianMs), fmt(r.js.totalMedianMs),
    fmt(r.wasmBase.eriMedianMs), fmt(r.wasmBase.scfMedianMs), fmt(r.wasmBase.totalMedianMs), sp(r.js.totalMedianMs, r.wasmBase.totalMedianMs),
    fmt(r.wasmSimd.eriMedianMs), fmt(r.wasmSimd.scfMedianMs), fmt(r.wasmSimd.totalMedianMs), sp(r.js.totalMedianMs, r.wasmSimd.totalMedianMs),
  ].map(v => `"${v}"`).join(','));
  return head.join(',') + '\n' + data.join('\n') + '\n';
}

function download(name: string, content: string, mime: string): void {
  const blob = new Blob([content], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function toCSV(out: BenchOutput | null): string {
  if (!out) return '';
  const head = ['workload_id', 'workload_label', 'backend', 'warmup_ms', 'run1_ms', 'run2_ms', 'run3_ms', 'median_ms', 'error'];
  const rows = out.workloads.map(w => [
    w.id, w.label, out.backend,
    isFinite(w.warmupMs) ? w.warmupMs.toFixed(2) : '',
    w.runs[0]?.toFixed(2) ?? '',
    w.runs[1]?.toFixed(2) ?? '',
    w.runs[2]?.toFixed(2) ?? '',
    isFinite(w.medianMs) ? w.medianMs.toFixed(2) : '',
    w.error ?? '',
  ].map(v => `"${String(v).replace(/"/g, '""')}"`).join(','));
  const meta = [
    `# device.ua,"${out.device.ua.replace(/"/g, '""')}"`,
    `# device.platform,"${out.device.platform}"`,
    `# device.deviceMemoryGB,"${out.device.deviceMemoryGB ?? ''}"`,
    `# device.hardwareConcurrency,"${out.device.hardwareConcurrency}"`,
    `# device.screen,"${out.device.screen}"`,
    `# backend,"${out.backend}"`,
    `# initialLoadMs,"${out.initialLoadMs ?? ''}"`,
    `# timestamp,"${out.device.timestamp}"`,
  ];
  return meta.join('\n') + '\n' + head.join(',') + '\n' + rows.join('\n') + '\n';
}

function injectStyles(): void {
  if (document.getElementById('cb-styles')) return;
  const s = document.createElement('style');
  s.id = 'cb-styles';
  s.textContent = `
    .cb-page { max-width: 900px; margin: 24px auto; padding: 0 16px; color: var(--color-text); }
    .cb-header h1 { margin: 0 0 4px; font-size: 1.6rem; }
    .cb-subtitle { color: var(--color-text-dim); font-size: 0.92rem; margin-bottom: 24px; }
    .cb-panel { background: var(--color-surface); border: 1px solid var(--color-border); border-radius: 8px; padding: 16px 20px; margin-bottom: 16px; }
    .cb-panel h2 { font-size: 1.1rem; margin: 0 0 12px; color: var(--color-accent); }
    .cb-note { font-size: 0.82rem; color: var(--color-text-dim); margin: 0 0 12px; }
    table.cb-info { width: 100%; font-size: 0.85rem; }
    table.cb-info th { text-align: left; padding: 4px 12px 4px 0; color: var(--color-text-dim); font-weight: 500; vertical-align: top; width: 160px; }
    table.cb-info td { padding: 4px 0; font-family: 'Cascadia Code', monospace; font-size: 0.82rem; word-break: break-all; }
    .cb-workloads { list-style: none; padding: 0; margin: 0 0 16px; }
    .cb-workloads li { display: flex; justify-content: space-between; padding: 6px 0; border-bottom: 1px solid var(--color-border); font-size: 0.9rem; }
    .cb-wl-label { font-weight: 500; }
    .cb-wl-status { color: var(--color-text-dim); font-family: 'Cascadia Code', monospace; font-size: 0.82rem; }
    .cb-panel button { background: var(--color-accent); color: #fff; border: none; padding: 8px 18px; border-radius: 6px; cursor: pointer; font-size: 0.9rem; margin-right: 8px; }
    .cb-panel button:disabled { opacity: 0.4; cursor: default; }
    .cb-panel button:hover:not([disabled]) { opacity: 0.85; }
    .cb-results-table { width: 100%; font-size: 0.85rem; border-collapse: collapse; margin-bottom: 12px; }
    .cb-results-table th, .cb-results-table td { padding: 6px 10px; border-bottom: 1px solid var(--color-border); }
    .cb-results-table th { background: var(--color-input); text-align: left; }
    .cb-results-table td.num, .cb-results-table th.num { text-align: right; font-family: 'Cascadia Code', monospace; }
    .cb-export { margin-bottom: 12px; }
    .cb-json pre { background: var(--color-input); padding: 12px; border-radius: 6px; font-size: 0.78rem; overflow-x: auto; max-height: 400px; }
    .cb-err { color: var(--color-error, #e05050); font-size: 0.75rem; }
    .cb-stress-controls { display: flex; flex-wrap: wrap; gap: 12px 18px; align-items: flex-end; margin-bottom: 4px; }
    .cb-stress-controls label { display: flex; flex-direction: column; gap: 4px; font-size: 0.78rem; color: var(--color-text-dim); }
    .cb-stress-controls select, .cb-stress-controls input { background: var(--color-input); color: var(--color-text); border: 1px solid var(--color-border); border-radius: 5px; padding: 5px 8px; font-size: 0.85rem; }
    .cb-stress-controls input[type="text"] { min-width: 260px; font-family: 'Cascadia Code', monospace; font-size: 0.78rem; }
    .cb-stress-row { background: var(--color-input); padding: 10px 12px; border-radius: 6px; font-size: 0.75rem; overflow-x: auto; white-space: pre; margin: 0 0 10px; }
  `;
  document.head.appendChild(s);
}

render();
