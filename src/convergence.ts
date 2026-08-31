/**
 * GANSU Lite — Basis Set Convergence demo page.
 * Computes RHF energy for a fixed molecule across multiple basis sets,
 * plotting energy vs basis set to visualize convergence toward CBS limit.
 */

import { parseXYZ } from './core/parseXYZ';
import { BasisSet } from './core/basisSet';
import { Molecular } from './core/molecular';
import { theorySelectHTML, buildHFOrDFT, HEAVY_ITERATIVE_EXCLUDE, type TheoryChoice } from './ui/theoryControls';
import { RHF } from './core/rhf';
import { initTheme, toggleTheme, isDark, getThemeColors } from './ui/theme';
import { t, initLang, toggleLang } from './ui/i18n';
import { renderHeader } from './ui/nav';

// ── Molecule definitions ─────────────────────────────────────────────

interface MolDef {
  id: string;
  labelKey: string;
  descKey: string;
  charge: number;
  mult: number;
  xyz: string;
}

const MOLECULES: MolDef[] = [
  {
    id: 'h2', labelKey: 'conv.scenH2', descKey: 'conv.descH2',
    charge: 0, mult: 1,
    xyz: '2\nH2\nH  0.0  0.0  0.0\nH  0.0  0.0  0.740000',
  },
  {
    id: 'hf', labelKey: 'conv.scenHF', descKey: 'conv.descHF',
    charge: 0, mult: 1,
    xyz: '2\nHF\nH  0.0  0.0  0.0\nF  0.0  0.0  0.917000',
  },
  {
    id: 'h2o', labelKey: 'conv.scenH2O', descKey: 'conv.descH2O',
    charge: 0, mult: 1,
    xyz: (() => {
      const R = 0.96, ang = 104 * Math.PI / 180, h = ang / 2;
      const hx = R * Math.sin(h), hz = R * Math.cos(h);
      return `3\nH2O\nO  0.0  0.0  0.0\nH  ${hx.toFixed(6)}  0.0  ${hz.toFixed(6)}\nH  ${(-hx).toFixed(6)}  0.0  ${hz.toFixed(6)}`;
    })(),
  },
  {
    id: 'lih', labelKey: 'conv.scenLiH', descKey: 'conv.descLiH',
    charge: 0, mult: 1,
    xyz: '2\nLiH\nLi  0.0  0.0  0.0\nH   0.0  0.0  1.596000',
  },
  {
    id: 'n2', labelKey: 'conv.scenN2', descKey: 'conv.descN2',
    charge: 0, mult: 1,
    xyz: '2\nN2\nN  0.0  0.0  0.0\nN  0.0  0.0  1.098000',
  },
  {
    id: 'nh3', labelKey: 'conv.scenNH3', descKey: 'conv.descNH3',
    charge: 0, mult: 1,
    xyz: (() => {
      const RNH = 1.012, h = 0.381; // pyramidal height
      const r = Math.sqrt(RNH * RNH - h * h);
      const s3 = Math.sqrt(3) / 2;
      return [
        '4', 'NH3',
        `N   0.000000  0.000000  ${h.toFixed(6)}`,
        `H   ${r.toFixed(6)}  0.000000  0.000000`,
        `H   ${(-r / 2).toFixed(6)}  ${(r * s3).toFixed(6)}  0.000000`,
        `H   ${(-r / 2).toFixed(6)}  ${(-r * s3).toFixed(6)}  0.000000`,
      ].join('\n');
    })(),
  },
];

const BASIS_LIST = ['STO-3G', 'STO-6G', '3-21G', '6-31G', '6-31G(d,p)', 'cc-pVDZ', 'def2-SVP', 'def2-TZVP'];
/** Triple-zeta basis sets are skipped for DFT (each is 30s-2min × molecule × DFT). */
const DFT_SKIP_BASIS = new Set(['cc-pVTZ', 'def2-TZVP']);

// ── Types ────────────────────────────────────────────────────────────

interface ConvPoint {
  basis: string;
  nbasis: number;
  energy: number;
  timeMs: number;
  error?: string;
}

// ── State ────────────────────────────────────────────────────────────

let selectedMol = MOLECULES[0];
let results: ConvPoint[] = [];
let running = false;
let stopRequested = false;
let done = false;
let elapsedMs = 0;
let theoryChoice: TheoryChoice = 'HF';

// ── Basis loading ────────────────────────────────────────────────────

const basisCache = new Map<string, BasisSet>();

async function loadBasis(name: string): Promise<BasisSet> {
  const cached = basisCache.get(name);
  if (cached) return cached;
  const url = `${import.meta.env.BASE_URL}basis/${name.toLowerCase()}.gbs`;
  const resp = await fetch(url);
  if (!resp.ok) throw new Error(`Failed to load basis set: ${name}`);
  const text = await resp.text();
  const bs = BasisSet.fromGBS(text);
  basisCache.set(name, bs);
  return bs;
}

// ── Main render ──────────────────────────────────────────────────────

const root = document.getElementById('app')!;

function render(): void {
  root.innerHTML = `
    <div class="opt-page">
      ${renderHeader('convergence')}

      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${t('conv.molecule')}</h2>
          <div class="opt-scenario-grid" id="mol-grid"></div>

          <div class="theory-row" style="margin:8px 0;display:flex;align-items:center;gap:8px;font-size:0.9rem;">
            <span>Theory:</span>${theorySelectHTML('theory-sel', theoryChoice, '', HEAVY_ITERATIVE_EXCLUDE)}
          </div>

          <button id="run-btn" class="opt-run-btn" ${running ? 'disabled' : ''}>
            ${running ? t('conv.running') : t('conv.run')}
          </button>
          ${running ? `<button id="stop-btn" class="opt-stop-btn">${t('conv.stop')}</button>` : ''}

          <div id="progress-area"></div>
          <div id="result-table"></div>
        </div>

        <div class="opt-panel opt-graph-panel">
          <div id="graph-container">
            ${!done && !running
              ? `<p class="opt-hint">${t('conv.waiting')}</p>`
              : ''}
          </div>
        </div>
      </div>
    </div>`;

  injectStyles();

  // Molecule cards
  const grid = root.querySelector('#mol-grid')!;
  const row = document.createElement('div');
  row.className = 'opt-category-row';
  for (const mol of MOLECULES) {
    const card = document.createElement('div');
    card.className = 'opt-scenario-card' + (mol.id === selectedMol.id ? ' selected' : '');
    card.innerHTML = `<strong>${t(mol.labelKey)}</strong><span class="conv-desc">${t(mol.descKey)}</span>`;
    card.addEventListener('click', () => {
      if (running) return;
      selectedMol = mol;
      results = [];
      done = false;
      render();
    });
    row.appendChild(card);
  }
  grid.appendChild(row);

  // Event listeners
  root.querySelector('#nav-theme')!.addEventListener('click', () => { toggleTheme(); render(); });
  root.querySelector('#nav-lang')!.addEventListener('click', () => { toggleLang(); render(); });
  const theorySel = root.querySelector<HTMLSelectElement>('#theory-sel');
  if (theorySel) theorySel.addEventListener('change', () => { theoryChoice = theorySel.value as TheoryChoice; });
  root.querySelector('#run-btn')!.addEventListener('click', () => { if (!running) runConvergence(); });
  root.querySelector('#stop-btn')?.addEventListener('click', () => { stopRequested = true; });

  // Render graph & table if data exists
  if (results.length > 0) {
    renderGraph();
    renderTable();
  }
}

// ── Computation ──────────────────────────────────────────────────────

async function runConvergence(): Promise<void> {
  running = true;
  stopRequested = false;
  results = [];
  done = false;
  render();

  const betaToAlpha = Math.floor((selectedMol.mult - 1) / 2);
  const startTime = performance.now();
  const atoms = parseXYZ(selectedMol.xyz);

  const isDFT = theoryChoice !== 'HF';
  for (let i = 0; i < BASIS_LIST.length; i++) {
    if (stopRequested) break;

    const basisName = BASIS_LIST[i];
    // Skip slow combos: triple-zeta basis with DFT can take 1+ min/molecule
    if (isDFT && DFT_SKIP_BASIS.has(basisName)) {
      results.push({ basis: basisName, nbasis: 0, energy: NaN, timeMs: 0, error: 'skipped (DFT × TZ too slow)' });
      continue;
    }
    updateProgress(basisName, i, BASIS_LIST.length);

    try {
      const t0 = performance.now();
      const basis = await loadBasis(basisName);
      const mol = new Molecular(atoms, basis, selectedMol.charge, betaToAlpha);
      const rhf = await buildHFOrDFT(mol, basis, theoryChoice) as RHF;
      const energy = await rhf.solve({ eriBackend: 'js' });
      const timeMs = performance.now() - t0;

      results.push({
        basis: basisName,
        nbasis: mol.numBasis,
        energy,
        timeMs,
      });

      renderGraph();
      renderTable();
    } catch (err) {
      const msg = err instanceof Error ? err.message : String(err);
      console.error(`Error computing ${basisName}:`, err);
      results.push({ basis: basisName, nbasis: 0, energy: NaN, timeMs: 0, error: msg });
    }

    await new Promise<void>(r => setTimeout(r, 0));
  }

  elapsedMs = performance.now() - startTime;
  running = false;
  done = true;
  render();
}

// ── Progress ─────────────────────────────────────────────────────────

function updateProgress(basisName: string, current: number, total: number): void {
  const area = root.querySelector('#progress-area');
  if (!area) return;
  const pct = (current / total * 100).toFixed(0);
  area.innerHTML = `
    <div class="opt-progress"><div class="opt-progress-bar" style="width:${pct}%"></div></div>
    <div class="opt-progress-text">${basisName}... (${current}/${total})</div>`;
}

// ── Results table ────────────────────────────────────────────────────

function renderTable(): void {
  const el = root.querySelector('#result-table');
  if (!el || results.length === 0) return;

  const valid = results.filter(r => isFinite(r.energy));
  const best = valid.length > 0 ? valid[valid.length - 1] : null;
  const totalTime = done ? ` (${(elapsedMs / 1000).toFixed(1)}s)` : '';

  let html = `<div class="opt-summary">
    <h3>${done ? t('conv.done') + totalTime : t('conv.running')}</h3>
    <table>
      <tr>
        <th>${t('conv.colBasis')}</th>
        <th style="text-align:right">M</th>
        <th style="text-align:right">${t('conv.colEnergy')}</th>
        <th style="text-align:right">${t('conv.colTime')}</th>
      </tr>`;

  for (const r of results) {
    const isBest = best && r.basis === best.basis;
    const tStr = r.timeMs < 1000 ? `${r.timeMs.toFixed(0)} ms` : `${(r.timeMs / 1000).toFixed(1)} s`;
    if (r.error) {
      html += `<tr style="color:var(--color-error)">
        <td>${r.basis}</td>
        <td></td>
        <td colspan="2" style="font-size:0.7rem">${r.error}</td>
      </tr>`;
    } else {
      html += `<tr${isBest ? ' style="font-weight:600"' : ''}>
        <td>${r.basis}</td>
        <td style="text-align:right">${r.nbasis}</td>
        <td style="text-align:right;font-family:monospace">${r.energy.toFixed(6)}</td>
        <td style="text-align:right">${tStr}</td>
      </tr>`;
    }
  }

  html += '</table></div>';
  el.innerHTML = html;
}

// ── SVG Graph ────────────────────────────────────────────────────────

function renderGraph(): void {
  const container = root.querySelector('#graph-container');
  if (!container || results.length === 0) return;

  const tc = getThemeColors();
  const width = 520, height = 360;
  const ml = 80, mr = 24, mt = 36, mb = 44;
  const pw = width - ml - mr, ph = height - mt - mb;

  const valid = results.filter(r => isFinite(r.energy));
  if (valid.length === 0) return;

  const energies = valid.map(r => r.energy);
  const yMin = Math.min(...energies), yMax = Math.max(...energies);
  const yPad = Math.max((yMax - yMin) * 0.15, 0.005);
  const yLow = yMin - yPad, yHigh = yMax + yPad;

  const toX = (i: number) => ml + (i / (BASIS_LIST.length - 1)) * pw;
  const toY = (y: number) => mt + ph - ((y - yLow) / (yHigh - yLow)) * ph;

  let svg = `<svg width="${width}" height="${height}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" style="max-width:100%;">`;

  // Background
  svg += `<rect x="${ml}" y="${mt}" width="${pw}" height="${ph}" fill="${tc.surface}" rx="2"/>`;

  // Grid lines (5 divisions)
  for (let i = 0; i <= 5; i++) {
    const yVal = yLow + (yHigh - yLow) * i / 5;
    const yy = toY(yVal);
    svg += `<line x1="${ml}" y1="${yy}" x2="${ml + pw}" y2="${yy}" stroke="${tc.grid}" stroke-width="0.5"/>`;
    svg += `<text x="${ml - 6}" y="${yy + 3}" text-anchor="end" font-size="9" fill="${tc.dim}">${yVal.toFixed(4)}</text>`;
  }

  // Axes
  svg += `<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt + ph}" stroke="${tc.axis}" stroke-width="1"/>`;
  svg += `<line x1="${ml}" y1="${mt + ph}" x2="${ml + pw}" y2="${mt + ph}" stroke="${tc.axis}" stroke-width="1"/>`;

  // X-axis labels (basis set names)
  for (let i = 0; i < BASIS_LIST.length; i++) {
    const xx = toX(i);
    svg += `<line x1="${xx}" y1="${mt + ph}" x2="${xx}" y2="${mt + ph + 4}" stroke="${tc.axis}" stroke-width="1"/>`;
    svg += `<text x="${xx}" y="${mt + ph + 16}" text-anchor="middle" font-size="8.5" fill="${tc.dim}">${BASIS_LIST[i]}</text>`;
  }

  // Accent color for series
  const seriesColor = isDark() ? '#00d4ff' : '#0077cc';

  // Line connecting points
  if (valid.length >= 2) {
    let path = '';
    for (let i = 0; i < valid.length; i++) {
      const idx = BASIS_LIST.indexOf(valid[i].basis);
      const px = toX(idx), py = toY(valid[i].energy);
      path += i === 0 ? `M${px.toFixed(1)},${py.toFixed(1)}` : ` L${px.toFixed(1)},${py.toFixed(1)}`;
    }
    svg += `<path d="${path}" fill="none" stroke="${seriesColor}" stroke-width="2"/>`;
  }

  // Data points + nbasis labels
  for (const r of valid) {
    const idx = BASIS_LIST.indexOf(r.basis);
    const px = toX(idx), py = toY(r.energy);
    svg += `<circle cx="${px}" cy="${py}" r="4" fill="${seriesColor}"/>`;
    // nbasis label
    svg += `<text x="${px}" y="${py - 8}" text-anchor="middle" font-size="8" fill="${tc.dim}">${r.nbasis}</text>`;
  }

  // Best estimate dashed line (from last computed point)
  if (valid.length >= 2) {
    const bestE = valid[valid.length - 1].energy;
    const by = toY(bestE);
    svg += `<line x1="${ml}" y1="${by}" x2="${ml + pw}" y2="${by}" stroke="${seriesColor}" stroke-width="0.8" stroke-dasharray="4,3" opacity="0.5"/>`;
  }

  // Title & axis labels
  svg += `<text x="${ml + pw / 2}" y="20" text-anchor="middle" font-size="12" font-weight="600" fill="${tc.titleSvg}">${t('conv.graphTitle')}</text>`;
  svg += `<text x="${ml + pw / 2}" y="${height - 4}" text-anchor="middle" font-size="10" fill="${tc.dim}">${t('conv.xBasis')}</text>`;
  svg += `<text x="14" y="${mt + ph / 2}" text-anchor="middle" font-size="10" fill="${tc.dim}" transform="rotate(-90,14,${mt + ph / 2})">${t('conv.yEnergy')}</text>`;

  svg += '</svg>';
  container.innerHTML = svg;
}

// ── Styles ───────────────────────────────────────────────────────────

let stylesInjected = false;

function injectStyles(): void {
  if (stylesInjected) return;
  stylesInjected = true;

  const style = document.createElement('style');
  style.textContent = `
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      background: var(--color-bg);
      color: var(--color-text);
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
      min-height: 100vh;
    }
    .opt-page { max-width: 960px; margin: 0 auto; padding: 16px 20px; }

    .opt-header {
      display: flex; justify-content: space-between; align-items: center;
      margin-bottom: 20px; padding-bottom: 12px;
      border-bottom: 1px solid var(--color-border);
    }
    .opt-header-left { display: flex; flex-direction: column; gap: 2px; }
    .opt-back-link { font-size: 0.75rem; color: var(--color-link); text-decoration: none; }
    .opt-back-link:hover { text-decoration: underline; }
    .opt-title { font-size: 1rem; font-weight: 600; }
    .opt-header-right { display: flex; gap: 6px; }
    .opt-header-right button {
      background: none; border: 1px solid var(--color-border); border-radius: 6px;
      padding: 4px 8px; cursor: pointer; color: var(--color-text); font-size: 0.8rem;
    }
    .opt-header-right button:hover { background: var(--color-surface-alt); }

    .opt-content { display: flex; gap: 20px; }
    .opt-panel {
      background: var(--color-surface); border: 1px solid var(--color-border);
      border-radius: 10px; padding: 16px 18px;
    }
    .opt-controls { flex: 0 0 300px; }
    .opt-graph-panel { flex: 1; min-width: 0; display: flex; flex-direction: column; align-items: center; padding-top: 8px; }
    .opt-controls h2 {
      font-size: 0.78rem; font-weight: 600; color: var(--color-text-secondary);
      text-transform: uppercase; letter-spacing: 0.04em;
      margin: 14px 0 8px;
    }
    .opt-controls h2:first-child { margin-top: 0; }

    .opt-scenario-grid { display: flex; flex-direction: column; gap: 4px; }
    .opt-category-row { display: flex; flex-wrap: wrap; gap: 5px; }
    .opt-scenario-card {
      padding: 5px 10px; border: 1px solid var(--color-border); border-radius: 6px;
      cursor: pointer; transition: all 0.15s; flex: 0 0 auto;
    }
    .opt-scenario-card:hover { background: var(--color-surface-alt); }
    .opt-scenario-card.selected {
      border-color: var(--color-accent); background: var(--color-surface-alt);
      box-shadow: 0 0 0 1px var(--color-accent);
    }
    .opt-scenario-card strong { display: block; font-size: 0.78rem; white-space: nowrap; }
    .conv-desc { display: block !important; font-size: 0.68rem; color: var(--color-text-dim); }

    .opt-run-btn {
      width: 100%; margin-top: 14px; padding: 10px; border: none; border-radius: 8px;
      font-size: 0.85rem; font-weight: 600; cursor: pointer;
      background: var(--color-accent); color: var(--color-accent-on);
      transition: background 0.15s;
    }
    .opt-run-btn:hover:not([disabled]) { background: var(--color-accent-hover); }
    .opt-run-btn[disabled] { opacity: 0.6; cursor: not-allowed; }
    .opt-stop-btn {
      width: 100%; margin-top: 6px; padding: 8px; border: 1px solid var(--color-error); border-radius: 8px;
      font-size: 0.8rem; font-weight: 600; cursor: pointer;
      background: none; color: var(--color-error);
    }

    .opt-progress {
      height: 6px; background: var(--color-progress-bg, #e0e4ea);
      border-radius: 3px; margin-top: 12px; overflow: hidden;
    }
    .opt-progress-bar {
      height: 100%; background: var(--color-accent); border-radius: 3px;
      transition: width 0.2s ease;
    }
    .opt-progress-text {
      font-size: 0.72rem; color: var(--color-text-dim); text-align: center; margin-top: 4px;
    }

    .opt-summary {
      margin-top: 14px; padding: 12px; background: var(--color-surface-alt);
      border-radius: 8px; border: 1px solid var(--color-border);
    }
    .opt-summary h3 {
      font-size: 0.82rem; color: var(--color-converged); margin-bottom: 8px;
    }
    .opt-summary table { width: 100%; font-size: 0.78rem; border-collapse: collapse; }
    .opt-summary th {
      font-size: 0.68rem; color: var(--color-text-dim); text-align: left;
      padding: 2px 4px; font-weight: 500;
    }
    .opt-summary td { padding: 3px 4px; }

    #graph-container { width: 100%; text-align: center; }
    .opt-hint { color: var(--color-text-dim); font-size: 0.85rem; padding: 60px 20px; }

    @media (max-width: 700px) {
      .opt-content { flex-direction: column; }
      .opt-controls { flex: none; }
    }
  `;
  document.head.appendChild(style);
}

// ── Init ─────────────────────────────────────────────────────────────

initLang();
initTheme();
render();
