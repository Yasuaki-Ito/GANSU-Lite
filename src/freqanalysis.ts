/**
 * GANSU Lite — Vibrational Analysis demo page.
 * Runs SCF → analytic gradient → numerical Hessian → frequencies, IR, thermochemistry.
 * Preset molecules only (safe computation times).
 */

import { parseXYZ } from './core/parseXYZ';
import { BasisSet } from './core/basisSet';
import { Molecular } from './core/molecular';
import type { DFTConfig } from './core/builder';
import { theorySelectHTML, buildHFOrDFT, theoryDisplayLabel, HEAVY_ITERATIVE_EXCLUDE, type TheoryChoice } from './ui/theoryControls';
import { RHF } from './core/rhf';
import { computeHessianAuto, computeThermodynamics, getAtomicMass, type HessianResult, type ThermoData } from './core/hessian';
import { computeRHFGradient } from './core/gradient';
import { initTheme, toggleTheme, getThemeColors } from './ui/theme';
import { t, initLang, toggleLang } from './ui/i18n';
import { renderHeader } from './ui/nav';
import { atomicNumberToElementName } from './core/constants';
import { ProgressTracker, type StepDef } from './ui/progressTracker';
import { initWasm, getActiveBackend } from './core/eriWasm';

// ── Scenarios (pre-optimized geometries, guaranteed fast) ────────────

interface FreqScenario {
  id: string;
  label: string;
  desc: string;
  charge: number;
  atomicNumbers: number[];
  coords: number[];       // bohr, flat [x0,y0,z0,...]
  estTime: string;        // display estimate
}

const B = 1.8897259886;
const SCENARIOS: FreqScenario[] = [
  // ── Diatomics ──
  { id: 'h2', label: 'H\u2082', desc: 'R\u2091=0.712\u00c5, homonuclear',
    charge: 0, atomicNumbers: [1, 1], estTime: '<1s',
    coords: [0, 0, 0, 0, 0, 0.7122 * B] },
  { id: 'lih', label: 'LiH', desc: 'R\u2091=1.596\u00c5',
    charge: 0, atomicNumbers: [3, 1], estTime: '~1s',
    coords: [0, 0, 0, 0, 0, 1.596 * B] },
  { id: 'hf', label: 'HF', desc: 'R\u2091=0.956\u00c5',
    charge: 0, atomicNumbers: [1, 9], estTime: '~2s',
    coords: [0, 0, 0, 0, 0, 0.9555 * B] },
  { id: 'n2', label: 'N\u2082', desc: 'R\u2091=1.098\u00c5, homonuclear',
    charge: 0, atomicNumbers: [7, 7], estTime: '~6s',
    coords: [0, 0, 0, 0, 0, 1.098 * B] },
  { id: 'co', label: 'CO', desc: 'R\u2091=1.128\u00c5',
    charge: 0, atomicNumbers: [6, 8], estTime: '~6s',
    coords: [0, 0, 0, 0, 0, 1.128 * B] },
  // ── Triatomics ──
  { id: 'h2o', label: 'H\u2082O', desc: 'bent, \u03b8=104.5\u00b0',
    charge: 0, atomicNumbers: [8, 1, 1], estTime: '~6s',
    coords: [0, 0, 0.127*B, 0, 0.758*B, -0.509*B, 0, -0.758*B, -0.509*B] },
  { id: 'hcn', label: 'HCN', desc: 'linear',
    charge: 0, atomicNumbers: [1, 6, 7], estTime: '~12s',
    coords: [0, 0, 0, 0, 0, 1.066*B, 0, 0, (1.066+1.156)*B] },
  { id: 'co2', label: 'CO\u2082', desc: 'linear D\u221e\u2095',
    charge: 0, atomicNumbers: [8, 6, 8], estTime: '~12s',
    coords: [0, 0, -1.162*B, 0, 0, 0, 0, 0, 1.162*B] },
  // ── 4+ atoms ──
  { id: 'bh3', label: 'BH\u2083', desc: 'planar D\u2083\u2095',
    charge: 0, atomicNumbers: [5, 1, 1, 1], estTime: '~8s',
    coords: (() => { const r = 1.19*B;
      return [0,0,0, r,0,0, -r*.5,r*Math.sqrt(3)/2,0, -r*.5,-r*Math.sqrt(3)/2,0]; })() },
  { id: 'c2h2', label: 'C\u2082H\u2082', desc: 'linear, \u03c0 bonds',
    charge: 0, atomicNumbers: [1, 6, 6, 1], estTime: '~16s',
    coords: [0,0,0, 0,0,1.06*B, 0,0,(1.06+1.20)*B, 0,0,(1.06+1.20+1.06)*B] },
  { id: 'nh3', label: 'NH\u2083', desc: 'pyramidal C\u2083\u1d65',
    charge: 0, atomicNumbers: [7, 1, 1, 1], estTime: '~15s',
    coords: (() => {
      const R = 1.012*B, th = 107*Math.PI/180, z = R*Math.cos(th), r = R*Math.sin(th);
      return [0,0,0, r,0,z, -r*.5,r*Math.sqrt(3)/2,z, -r*.5,-r*Math.sqrt(3)/2,z];
    })() },
  { id: 'ch4', label: 'CH\u2084', desc: 'tetrahedral T\u1d48',
    charge: 0, atomicNumbers: [6, 1, 1, 1, 1], estTime: '~25s',
    coords: (() => { const a = 1.089*B/Math.sqrt(3);
      return [0,0,0, a,a,a, a,-a,-a, -a,a,-a, -a,-a,a]; })() },
];

// ── State ────────────────────────────────────────────────────────────

let selectedScen = SCENARIOS[0];
let optimizeGeom = true;
let running = false;
let stopRequested = false;
let theoryChoice: TheoryChoice = 'HF';
let result: { hess: HessianResult; thermo: ThermoData; energy: number } | null = null;
let elapsed = '';

// ── Basis loading ────────────────────────────────────────────────────

let basisSTO: BasisSet | null = null;
async function getBasis(): Promise<BasisSet> {
  if (basisSTO) return basisSTO;
  const url = `${import.meta.env.BASE_URL}basis/sto-3g.gbs`;
  basisSTO = BasisSet.fromGBS(await (await fetch(url)).text());
  return basisSTO;
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
    .opt-controls { flex: 0 0 260px; }
    .opt-results-panel { flex: 1; min-width: 0; }
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
    .opt-scenario-card.disabled { opacity: 0.4; cursor: not-allowed; }
    .opt-scenario-card.disabled:hover { background: transparent; border-color: var(--color-border); }
    .opt-scenario-card strong { font-size: 0.85rem; white-space: nowrap; }
    .conv-desc { font-size: 0.7rem; color: var(--color-text-dim); }
    .est-time { font-size: 0.65rem; color: var(--color-text-dim); margin-left: auto; }
    .opt-run-btn {
      width: 100%; margin-top: 12px; padding: 10px; border: none; border-radius: 8px;
      font-size: 0.85rem; font-weight: 600; cursor: pointer;
      background: var(--color-accent); color: var(--color-accent-on); transition: background 0.15s; }
    .opt-run-btn:hover:not([disabled]) { background: var(--color-accent-hover); }
    .opt-run-btn[disabled] { opacity: 0.6; cursor: not-allowed; }
    .opt-progress { height: 6px; background: var(--color-progress-bg, #e0e4ea); border-radius: 3px; margin-top: 10px; overflow: hidden; }
    .opt-progress-bar { height: 100%; background: var(--color-accent); border-radius: 3px; transition: width 0.2s ease; }
    .opt-progress-text { font-size: 0.72rem; color: var(--color-text-dim); text-align: center; margin-top: 4px; }
    .opt-hint { color: var(--color-text-dim); font-size: 0.85rem; padding: 40px 20px; text-align: center; }
    .freq-table { width: 100%; border-collapse: collapse; font-size: 0.78rem; margin-bottom: 10px; }
    .freq-table th, .freq-table td { padding: 3px 6px; text-align: left; border-bottom: 1px solid var(--color-border); }
    .freq-table th { font-size: 0.7rem; color: var(--color-text-dim); font-weight: 500; }
    .thermo-table { width: 100%; border-collapse: collapse; font-size: 0.78rem; }
    .thermo-table td { padding: 2px 6px; }
    .thermo-table td:first-child { color: var(--color-text-dim); }
    .section-title { font-size: 0.82rem; font-weight: 600; margin: 12px 0 6px; color: var(--color-text); }
    .section-title:first-child { margin-top: 0; }
    @media (max-width: 700px) { .opt-content { flex-direction: column; } .opt-controls { flex: none; } }
  `;
  document.head.appendChild(style);
}

// ── Main render ──────────────────────────────────────────────────────

const root = document.getElementById('app')!;

function render(): void {
  root.innerHTML = `
    <div class="opt-page">
      ${renderHeader('freqanalysis')}
      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${t('gopt.molecule')}</h2>
          <div class="opt-scenario-grid" id="scen-grid"></div>
          <div style="margin-top:8px;font-size:0.72rem;color:var(--color-text-dim);">STO-3G basis</div>
          <label style="display:flex;align-items:center;gap:6px;margin-top:8px;font-size:0.78rem;cursor:pointer;">
            <input type="checkbox" id="opt-geom-check" ${optimizeGeom ? 'checked' : ''} ${running ? 'disabled' : ''}>
            Optimize geometry first
          </label>
          <div class="theory-row" style="margin:8px 0;display:flex;align-items:center;gap:8px;font-size:0.9rem;">
            <span>Theory:</span>${theorySelectHTML('theory-sel', theoryChoice, '', HEAVY_ITERATIVE_EXCLUDE)}
          </div>
          <button id="run-btn" class="opt-run-btn" ${running ? 'disabled' : ''}>
            ${running ? 'Computing...' : 'Run Analysis'}
          </button>
          <div id="progress-area"></div>
        </div>
        <div class="opt-panel opt-results-panel">
          <div id="results-area">
            ${!result && !running ? '<p class="opt-hint">Select a molecule and click Run Analysis</p>' : ''}
          </div>
        </div>
      </div>
    </div>`;

  injectStyles();

  // Scenario cards
  // For DFT, scenarios with ≥4 atoms are disabled (6N grad evals × DFT ~5min+).
  const dftDisabled = (n: number) => theoryChoice !== 'HF' && n >= 4;
  const grid = root.querySelector('#scen-grid')!;
  for (const scen of SCENARIOS) {
    const card = document.createElement('div');
    const isDisabled = dftDisabled(scen.atomicNumbers.length);
    card.className = 'opt-scenario-card'
      + (scen.id === selectedScen.id ? ' selected' : '')
      + (isDisabled ? ' disabled' : '');
    const tip = isDisabled ? ` title="Disabled for DFT — ${scen.atomicNumbers.length} atoms × 6N grad evals × DFT is too slow. Use HF or pick a smaller molecule."` : '';
    card.innerHTML = `<strong>${scen.label}</strong><span class="conv-desc">${scen.desc}</span><span class="est-time">${scen.estTime}</span>`;
    if (tip) card.setAttribute('title', tip.slice(8, -1));
    card.addEventListener('click', () => {
      if (running || isDisabled) return;
      selectedScen = scen; result = null;
      render();
    });
    grid.appendChild(card);
  }

  root.querySelector('#nav-theme')!.addEventListener('click', () => { toggleTheme(); render(); });
  root.querySelector('#nav-lang')!.addEventListener('click', () => { toggleLang(); render(); });
  const theorySel = root.querySelector<HTMLSelectElement>('#theory-sel');
  if (theorySel) theorySel.addEventListener('change', () => {
    theoryChoice = theorySel.value as TheoryChoice;
    // If currently selected scenario is now disabled, fall back to first 3-atom-or-less scenario
    if (theoryChoice !== 'HF' && selectedScen.atomicNumbers.length >= 4) {
      const fallback = SCENARIOS.find(s => s.atomicNumbers.length < 4);
      if (fallback) { selectedScen = fallback; result = null; }
    }
    render();
  });
  root.querySelector('#run-btn')!.addEventListener('click', () => { if (!running) runAnalysis(); });
  root.querySelector<HTMLInputElement>('#opt-geom-check')!.addEventListener('change', (e) => {
    optimizeGeom = (e.target as HTMLInputElement).checked;
  });

  if (result) renderResults();
}

// ── Computation ──────────────────────────────────────────────────────

async function runAnalysis(): Promise<void> {
  running = true;
  result = null;
  render();

  const stepDefs: StepDef[] = [
    { id: 'setup', label: t('progress.setup') || 'Molecular Setup' },
    ...(optimizeGeom ? [{ id: 'geomopt', label: 'Geometry Optimization' }] : []),
    { id: 'scf', label: t('progress.scf') || 'SCF (RHF)' },
    { id: 'hessian', label: t('progress.hessian') || 'Hessian' },
    { id: 'freq', label: 'Vibrational Analysis' },
    { id: 'thermo', label: 'Thermochemistry' },
  ];

  const tracker = new ProgressTracker(stepDefs);

  try {
    const scen = selectedScen;
    const t0 = performance.now();

    // Setup
    tracker.startStep('setup');
    const base = (import.meta as any).env?.BASE_URL || '/';
    await initWasm(base);
    const wasmBackend = getActiveBackend();
    tracker.setBackend(wasmBackend === 'wasm-simd' || wasmBackend === 'wasm' ? wasmBackend : 'js');
    const basis = await getBasis();
    const eriBackend = wasmBackend !== 'js' ? 'wasm' as const : 'auto' as const;
    const buildAtoms = (coords: Float64Array | number[]) => scen.atomicNumbers.map((z, i) => ({
      atomicNumber: z,
      coordinate: { x: coords[3*i], y: coords[3*i+1], z: coords[3*i+2] },
      atomIndex: i,
    }));
    let curCoords = new Float64Array(scen.coords);
    tracker.completeStep('setup');

    // Geometry optimization (steepest descent + Armijo backtracking)
    if (optimizeGeom) {
      tracker.startStep('geomopt');
      const tol = 1e-3; // gradient threshold (Eh/bohr)
      const maxIter = 30;
      let prevE = NaN;
      let curE = NaN;
      for (let iter = 0; iter < maxIter; iter++) {
        const m = new Molecular(buildAtoms(curCoords), basis, scen.charge);
        const r = await buildHFOrDFT(m, basis, theoryChoice, 'RHF') as RHF;
        curE = await r.solve({ eriBackend });
        const dftCtx = (r.xcFunctional && r.grid) ? { functional: r.xcFunctional, grid: r.grid } : undefined;
        const g = computeRHFGradient(
          m.primitiveShells, m.atoms, m.cgtoNormalizationFactors, m.numBasis,
          m.numAlphaSpins, r.density, r.coefficients, r.orbitalEnergies, undefined, dftCtx,
        ).total;
        let maxF = 0;
        for (let i = 0; i < g.length; i++) maxF = Math.max(maxF, Math.abs(g[i]));
        tracker.updateStep('geomopt', `iter ${iter+1}  E=${curE.toFixed(6)}  |F|=${maxF.toExponential(2)}`);
        await new Promise<void>(r => setTimeout(r, 0));
        if (maxF < tol) break;
        // Armijo backtracking
        let alpha = 0.5;
        const dirGrad = -g.reduce((s, v) => s + v*v, 0);
        let accepted = false;
        for (let trial = 0; trial < 5; trial++) {
          const trialCoords = new Float64Array(curCoords.length);
          for (let i = 0; i < curCoords.length; i++) trialCoords[i] = curCoords[i] - alpha * g[i];
          const mt = new Molecular(buildAtoms(trialCoords), basis, scen.charge);
          const rt = await buildHFOrDFT(mt, basis, theoryChoice, 'RHF') as RHF;
          const eT = await rt.solve({ eriBackend });
          if (eT <= curE + 1e-4 * alpha * dirGrad) {
            curCoords = trialCoords;
            accepted = true;
            break;
          }
          alpha *= 0.5;
        }
        if (!accepted) break;
        if (!isNaN(prevE) && Math.abs(curE - prevE) < 1e-8) break;
        prevE = curE;
      }
      tracker.completeStep('geomopt', `E=${curE.toFixed(6)} Eh`);
    }

    const mol = new Molecular(buildAtoms(curCoords), basis, scen.charge);
    const rhf = await buildHFOrDFT(mol, basis, theoryChoice, 'RHF') as RHF;

    // SCF (final, at optimised geometry)
    tracker.startStep('scf');
    const energy = await rhf.solve({
      eriBackend,
      onIteration: (_iter, _e, de) => {
        tracker.updateStep('scf', `dE = ${de.toExponential(2)}`);
      },
    });
    tracker.completeStep('scf', `E = ${energy.toFixed(8)} Eh`);

    // Hessian (numerical via FD of analytic gradient; pass dftConfig if DFT chosen)
    tracker.startStep('hessian');
    const refCoords = curCoords;
    const hessDftConfig: DFTConfig | undefined = (theoryChoice === 'HF')
      ? undefined
      : { functional: theoryChoice, gridLevel: 'medium' };
    const hess = await computeHessianAuto(
      scen.atomicNumbers, refCoords, basis, scen.charge, 5e-4,
      (msg) => tracker.updateStep('hessian', msg),
      hessDftConfig,
    );
    tracker.completeStep('hessian');

    // Vibrational analysis
    tracker.startStep('freq');
    // frequencies already computed inside computeHessianAuto
    const nReal = hess.frequencies.filter(f => Math.abs(f) > 50).length;
    tracker.completeStep('freq', `${nReal} modes`);

    // Thermochemistry
    tracker.startStep('thermo');
    const thermo = computeThermodynamics(
      hess.frequencies, hess.masses, refCoords, scen.atomicNumbers.length,
      energy, 298.15, 1.0,
    );
    tracker.completeStep('thermo', `ZPE = ${thermo.ZPE.toFixed(2)} kcal/mol`);

    elapsed = ((performance.now() - t0) / 1000).toFixed(1);
    result = { hess, thermo, energy };

    setTimeout(() => tracker.close(), 800);
  } catch (e) {
    tracker.failRemaining();
    console.error(e);
    setTimeout(() => tracker.close(), 2000);
  }

  running = false;
  render();
}

// ── Results ──────────────────────────────────────────────────────────

function renderResults(): void {
  if (!result) return;
  const area = root.querySelector('#results-area');
  if (!area) return;

  const { hess, thermo: th, energy } = result;
  const c = getThemeColors();
  let html = '';

  // ── Summary ──
  html += `<div class="section-title">${selectedScen.label} / ${theoryDisplayLabel(theoryChoice)} / STO-3G &nbsp; <span style="font-weight:400;font-size:0.72rem;color:var(--color-text-dim);">(${elapsed}s)</span></div>`;
  html += `<div style="font-size:0.78rem;margin-bottom:10px;">E<sub>elec</sub> = ${energy.toFixed(8)} Eh</div>`;

  // ── Frequency table ──
  html += `<div class="section-title">Vibrational Frequencies</div>`;
  html += '<table class="freq-table"><tr><th>#</th><th>Freq (cm\u207b\u00b9)</th><th>IR (km/mol)</th><th>Type</th></tr>';
  let modeIdx = 1;
  for (let i = 0; i < hess.frequencies.length; i++) {
    const f = hess.frequencies[i];
    if (Math.abs(f) < 50) continue;
    const ir = hess.intensities[i];
    const type = f < 0 ? 'imaginary' : 'real';
    html += `<tr${f < 0 ? ' style="color:var(--color-error)"' : ''}>`;
    html += `<td>${modeIdx++}</td><td>${Math.abs(f).toFixed(1)}</td><td>${ir.toFixed(1)}</td><td>${type}</td></tr>`;
  }
  html += '</table>';

  // ── IR Spectrum SVG ──
  const vibFreqs: number[] = [], vibIR: number[] = [];
  for (let i = 0; i < hess.frequencies.length; i++) {
    if (hess.frequencies[i] > 50 && hess.intensities[i] > 0.01) {
      vibFreqs.push(hess.frequencies[i]);
      vibIR.push(hess.intensities[i]);
    }
  }
  if (vibFreqs.length > 0) {
    const W = 520, H = 180, PAD = { l: 45, r: 15, t: 24, b: 36 };
    const pw = W - PAD.l - PAD.r, ph = H - PAD.t - PAD.b;
    const xMin = 400, xMax = Math.max(5000, ...vibFreqs) + 200;
    const gamma = 30;
    const nPts = 300;
    const spectrum = new Float64Array(nPts);
    for (let p = 0; p < nPts; p++) {
      const nu = xMin + (xMax - xMin) * p / (nPts - 1);
      for (let k = 0; k < vibFreqs.length; k++) {
        const d = nu - vibFreqs[k];
        spectrum[p] += vibIR[k] * gamma * gamma / (d * d + gamma * gamma);
      }
    }
    const maxAbs = Math.max(...spectrum) || 1;
    const sx = (nu: number) => PAD.l + (nu - xMin) / (xMax - xMin) * pw;
    const sy = (a: number) => PAD.t + (a / maxAbs) * ph;

    let svg = `<svg viewBox="0 0 ${W} ${H}" style="width:100%;max-width:${W}px;background:${c.surface};border:1px solid var(--color-border);border-radius:6px;margin:8px 0;">`;
    svg += `<text x="${W/2}" y="15" text-anchor="middle" font-size="11" fill="${c.titleSvg}" font-weight="600">IR Spectrum</text>`;
    svg += `<line x1="${PAD.l}" y1="${PAD.t}" x2="${PAD.l}" y2="${H-PAD.b}" stroke="var(--color-border)"/>`;
    svg += `<line x1="${PAD.l}" y1="${H-PAD.b}" x2="${W-PAD.r}" y2="${H-PAD.b}" stroke="var(--color-border)"/>`;
    svg += `<text x="${PAD.l+pw/2}" y="${H-5}" text-anchor="middle" font-size="9" fill="${c.dim}">Wavenumber (cm\u207b\u00b9)</text>`;
    svg += `<text x="8" y="${PAD.t+ph/2}" text-anchor="middle" font-size="8" fill="${c.dim}" transform="rotate(-90,8,${PAD.t+ph/2})">Absorbance</text>`;
    for (let tick = 500; tick <= xMax; tick += 500) {
      if (tick < xMin) continue;
      const x = sx(tick);
      svg += `<line x1="${x}" y1="${H-PAD.b}" x2="${x}" y2="${H-PAD.b+3}" stroke="var(--color-border)"/>`;
      svg += `<text x="${x}" y="${H-PAD.b+13}" text-anchor="middle" font-size="7" fill="${c.dim}">${tick}</text>`;
    }
    const pts = Array.from(spectrum).map((a, p) => {
      const nu = xMin + (xMax - xMin) * p / (nPts - 1);
      return `${sx(nu).toFixed(1)},${sy(a).toFixed(1)}`;
    });
    svg += `<polygon points="${sx(xMin).toFixed(1)},${PAD.t} ${pts.join(' ')} ${sx(xMax).toFixed(1)},${PAD.t}" fill="${c.accent}" opacity="0.15"/>`;
    svg += `<polyline points="${pts.join(' ')}" fill="none" stroke="${c.accent}" stroke-width="1.5"/>`;
    for (let k = 0; k < vibFreqs.length; k++) {
      if (vibIR[k] < maxAbs * 0.05) continue;
      const x = sx(vibFreqs[k]), y = sy(vibIR[k]);
      svg += `<text x="${x}" y="${Math.max(y - 4, PAD.t + 10)}" text-anchor="middle" font-size="7.5" fill="${c.titleSvg}">${vibFreqs[k].toFixed(0)}</text>`;
    }
    svg += '</svg>';
    html += svg;
  } else {
    html += `<div style="font-size:0.78rem;color:var(--color-text-dim);margin:8px 0;padding:10px;background:var(--color-surface-alt);border-radius:6px;">No IR-active modes — homonuclear molecules (H\u2082, N\u2082, ...) have no dipole change during vibration, so infrared absorption is forbidden by symmetry.</div>`;
  }

  // ── Thermochemistry ──
  html += `<div class="section-title">Thermochemistry (${th.temperature.toFixed(1)} K, ${th.pressure.toFixed(0)} atm)</div>`;
  html += '<table class="thermo-table">';
  const row = (label: string, val: string) => `<tr><td>${label}</td><td>${val}</td></tr>`;
  html += row('Zero-Point Energy', `${th.ZPE.toFixed(2)} kcal/mol`);
  html += row('Thermal Correction to Energy', `${th.thermalCorr.toFixed(6)} Eh`);
  html += row('Thermal Correction to Enthalpy', `${th.enthalpyCorr.toFixed(6)} Eh`);
  html += row('Thermal Correction to Gibbs', `${th.gibbsCorr.toFixed(6)} Eh`);
  html += '<tr><td colspan="2" style="border-top:1px solid var(--color-border);height:6px;"></td></tr>';
  html += row('E + Thermal', `${th.E_tot.toFixed(6)} Eh`);
  html += row('H (Enthalpy)', `${th.H.toFixed(6)} Eh`);
  html += row('G (Gibbs Free Energy)', `${th.G.toFixed(6)} Eh`);
  html += '<tr><td colspan="2" style="border-top:1px solid var(--color-border);height:6px;"></td></tr>';
  html += row('S (Total)', `${th.S_tot.toFixed(2)} cal/(mol\u00b7K)`);
  html += row('&nbsp;&nbsp;Translational', `${th.S_trans.toFixed(2)}`);
  html += row('&nbsp;&nbsp;Rotational', `${th.S_rot.toFixed(2)}`);
  html += row('&nbsp;&nbsp;Vibrational', `${th.S_vib.toFixed(2)}`);
  html += row('C\u1d65 (Heat Capacity)', `${th.Cv_tot.toFixed(2)} cal/(mol\u00b7K)`);
  html += '</table>';

  area.innerHTML = html;
}

// ── Init ─────────────────────────────────────────────────────────────

initTheme(); initLang(); render();
