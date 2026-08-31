/**
 * GANSU Lite — DFT Functional Ladder demo page.
 * Compares HF → SVWN (LDA) → BLYP (GGA) → PBE (GGA) → B3LYP (Hybrid)
 * for the same molecule, visualising XC energy contributions as a staircase chart.
 */

import { parseXYZ } from './core/parseXYZ';
import { BasisSet } from './core/basisSet';
import { Molecular } from './core/molecular';
import { buildHF, type DFTConfig } from './core/builder';
import type { FunctionalName } from './core/xcFunctional';
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
    id: 'h2', labelKey: 'ladder.scenH2', descKey: 'ladder.descH2',
    charge: 0, mult: 1,
    xyz: '2\nH2\nH  0.0  0.0  0.0\nH  0.0  0.0  0.740000',
  },
  {
    id: 'lih', labelKey: 'ladder.scenLiH', descKey: 'ladder.descLiH',
    charge: 0, mult: 1,
    xyz: '2\nLiH\nLi  0.0  0.0  0.0\nH   0.0  0.0  1.596000',
  },
  {
    id: 'beh2', labelKey: 'ladder.scenBeH2', descKey: 'ladder.descBeH2',
    charge: 0, mult: 1,
    xyz: '3\nBeH2\nBe  0.0  0.0  0.0\nH  0.0  0.0  1.334000\nH  0.0  0.0  -1.334000',
  },
  {
    id: 'bh3', labelKey: 'ladder.scenBH3', descKey: 'ladder.descBH3',
    charge: 0, mult: 1,
    xyz: (() => {
      const R = 1.19;
      const angles = [0, 120, 240].map(a => a * Math.PI / 180);
      const hs = angles.map(a => `H  ${(R * Math.cos(a)).toFixed(6)}  ${(R * Math.sin(a)).toFixed(6)}  0.0`);
      return `4\nBH3\nB  0.0  0.0  0.0\n${hs.join('\n')}`;
    })(),
  },
  {
    id: 'hf', labelKey: 'ladder.scenHF', descKey: 'ladder.descHF',
    charge: 0, mult: 1,
    xyz: '2\nHF\nH  0.0  0.0  0.0\nF  0.0  0.0  0.917000',
  },
  {
    id: 'h2o', labelKey: 'ladder.scenH2O', descKey: 'ladder.descH2O',
    charge: 0, mult: 1,
    xyz: (() => {
      const R = 0.96, ang = 104 * Math.PI / 180, h = ang / 2;
      const hx = R * Math.sin(h), hz = R * Math.cos(h);
      return `3\nH2O\nO  0.0  0.0  0.0\nH  ${hx.toFixed(6)}  0.0  ${hz.toFixed(6)}\nH  ${(-hx).toFixed(6)}  0.0  ${hz.toFixed(6)}`;
    })(),
  },
  {
    id: 'nh3', labelKey: 'ladder.scenNH3', descKey: 'ladder.descNH3',
    charge: 0, mult: 1,
    xyz: (() => {
      const RNH = 1.012, h = 0.381;
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
  {
    id: 'ch4', labelKey: 'ladder.scenCH4', descKey: 'ladder.descCH4',
    charge: 0, mult: 1,
    xyz: (() => {
      const R = 1.089;
      const a = R / Math.sqrt(3);
      return [
        '5', 'CH4',
        'C   0.000000  0.000000  0.000000',
        `H   ${a.toFixed(6)}  ${a.toFixed(6)}  ${a.toFixed(6)}`,
        `H   ${a.toFixed(6)}  ${(-a).toFixed(6)}  ${(-a).toFixed(6)}`,
        `H   ${(-a).toFixed(6)}  ${a.toFixed(6)}  ${(-a).toFixed(6)}`,
        `H   ${(-a).toFixed(6)}  ${(-a).toFixed(6)}  ${a.toFixed(6)}`,
      ].join('\n');
    })(),
  },
  {
    id: 'co', labelKey: 'ladder.scenCO', descKey: 'ladder.descCO',
    charge: 0, mult: 1,
    xyz: '2\nCO\nC  0.0  0.0  0.0\nO  0.0  0.0  1.128000',
  },
  {
    id: 'n2', labelKey: 'ladder.scenN2', descKey: 'ladder.descN2',
    charge: 0, mult: 1,
    xyz: '2\nN2\nN  0.0  0.0  0.0\nN  0.0  0.0  1.098000',
  },
  {
    id: 'f2', labelKey: 'ladder.scenF2', descKey: 'ladder.descF2',
    charge: 0, mult: 1,
    xyz: '2\nF2\nF  0.0  0.0  0.0\nF  0.0  0.0  1.412000',
  },
];

// ── Types ────────────────────────────────────────────────────────────

interface MethodDef {
  id: string;
  label: string;
  rung: string;          // Jacob's Ladder rung
  functional?: FunctionalName;
}

const METHODS: MethodDef[] = [
  { id: 'HF',    label: 'HF',    rung: 'Hartree\u2013Fock' },
  { id: 'SVWN',  label: 'SVWN',  rung: 'LDA',   functional: 'SVWN' },
  { id: 'BLYP',  label: 'BLYP',  rung: 'GGA',   functional: 'BLYP' },
  { id: 'PBE',   label: 'PBE',   rung: 'GGA',   functional: 'PBE' },
  { id: 'B3LYP', label: 'B3LYP', rung: 'Hybrid', functional: 'B3LYP' },
];

interface DFTResult {
  method: MethodDef;
  totalEnergy: number;
  xcEnergy: number;     // difference from HF (XC + exchange correction)
  timeMs: number;
  error?: string;
}

// ── State ────────────────────────────────────────────────────────────

let selectedMol = MOLECULES[5]; // default H2O
let results: DFTResult[] = [];
let running = false;
let stopRequested = false;
let done = false;
let currentStep = '';

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

// ── Method colors ────────────────────────────────────────────────────

const METHOD_COLORS_LIGHT: Record<string, string> = {
  HF:    '#888888',
  SVWN:  '#1f77b4',
  BLYP:  '#2ca02c',
  PBE:   '#ff7f0e',
  B3LYP: '#d62728',
};
const METHOD_COLORS_DARK: Record<string, string> = {
  HF:    '#aaaaaa',
  SVWN:  '#4dabf7',
  BLYP:  '#51cf66',
  PBE:   '#ffa94d',
  B3LYP: '#ff6b6b',
};

function methodColor(id: string): string {
  return isDark() ? (METHOD_COLORS_DARK[id] ?? '#ccc') : (METHOD_COLORS_LIGHT[id] ?? '#666');
}

// ── Main render ──────────────────────────────────────────────────────

const root = document.getElementById('app')!;

function render(): void {
  root.innerHTML = `
    <div class="opt-page">
      ${renderHeader('dft')}

      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${t('dft.molecule')}</h2>
          <div class="opt-scenario-grid" id="scen-grid"></div>

          <div class="ladder-basis-info">
            <span style="color:var(--color-text-dim);font-size:0.72rem">${t('dft.basis')}: STO-3G &middot; Grid: medium</span>
          </div>

          <button id="run-btn" class="opt-run-btn" ${running ? 'disabled' : ''}>
            ${running ? t('dft.running') : t('dft.run')}
          </button>
          ${running ? `<button id="stop-btn" class="opt-stop-btn">${t('dft.stop')}</button>` : ''}

          <div id="progress-area"></div>
          <div id="summary-area"></div>
        </div>

        <div class="opt-panel opt-graph-panel">
          <div id="graph-container">
            ${!done && !running
              ? `<p class="opt-hint">${t('dft.waiting')}</p>`
              : ''}
          </div>
        </div>
      </div>
    </div>`;

  injectStyles();

  // Scenario cards
  const grid = root.querySelector('#scen-grid')!;
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
  root.querySelector('#run-btn')!.addEventListener('click', () => { if (!running) runLadder(); });
  root.querySelector('#stop-btn')?.addEventListener('click', () => { stopRequested = true; });

  // Render charts if data exists
  if (results.length > 0) {
    renderGraph();
    renderSummary();
  }
}

// ── Computation ──────────────────────────────────────────────────────

async function runLadder(): Promise<void> {
  running = true;
  stopRequested = false;
  results = [];
  done = false;
  render();

  const mol = selectedMol;

  try {
    const basis = await loadBasis('STO-3G');
    const atoms = parseXYZ(mol.xyz);
    const betaToAlpha = Math.floor((mol.mult - 1) / 2);
    const molecular = new Molecular(atoms, basis, mol.charge, betaToAlpha);

    let hfEnergy = 0;

    for (let i = 0; i < METHODS.length; i++) {
      if (stopRequested) break;

      const meth = METHODS[i];
      currentStep = meth.label;
      updateProgress(i);

      const t0 = performance.now();

      let energy: number;
      if (meth.functional) {
        const dftConfig: DFTConfig = { functional: meth.functional, gridLevel: 'medium' };
        const hf = buildHF(molecular, 'RHF', dftConfig);
        energy = await hf.solve({
          eriBackend: 'js',
          onProgress: (msg) => updateStepInfo(msg),
        });
      } else {
        // Pure HF
        const hf = buildHF(molecular, 'RHF');
        energy = await hf.solve({
          eriBackend: 'js',
          onProgress: (msg) => updateStepInfo(msg),
        });
        hfEnergy = energy;
      }
      const timeMs = performance.now() - t0;

      results.push({
        method: meth,
        totalEnergy: energy,
        xcEnergy: energy - hfEnergy,
        timeMs,
      });

      renderGraph();
      renderSummary();
      await new Promise<void>(r => setTimeout(r, 0));
    }
  } catch (err) {
    console.error('DFT ladder computation error:', err);
  }

  finish();
}

function finish(): void {
  running = false;
  done = true;
  currentStep = '';
  render();
}

// ── Progress ─────────────────────────────────────────────────────────

function updateProgress(stepIdx: number): void {
  const area = root.querySelector('#progress-area');
  if (!area) return;
  const pct = (stepIdx / METHODS.length * 100).toFixed(0);
  area.innerHTML = `
    <div class="opt-progress"><div class="opt-progress-bar" style="width:${pct}%"></div></div>
    <div class="opt-progress-text">${currentStep} (${stepIdx + 1}/${METHODS.length})</div>`;
}

function updateStepInfo(msg: string): void {
  const area = root.querySelector('.opt-progress-text');
  if (!area) return;
  // Shorten verbose messages for compact display
  let short = msg;
  if (msg.includes('Converged')) short = msg;
  else if (msg.includes('ERIs')) short = 'ERI ' + (msg.match(/(\d+)%/)?.[0] ?? '');
  else if (msg.includes('nuclear')) short = 'V_nn';
  else if (msg.includes('one-electron')) short = '1e integrals';
  else if (msg.includes('two-electron')) short = 'ERI';
  else if (msg.includes('transformation')) short = 'S^{-1/2}';
  else if (msg.includes('SAD')) short = 'SAD guess';
  else if (msg.includes('GWH')) short = 'GWH guess';
  else if (msg.includes('initial Fock')) short = 'H_core guess';
  else if (msg.includes('grid')) short = 'grid switch';
  else if (msg.includes('SCF did not')) short = 'not converged';
  area.textContent = `${currentStep}: ${short}`;
}

// ── Summary ──────────────────────────────────────────────────────────

function renderSummary(): void {
  const el = root.querySelector('#summary-area');
  if (!el || results.length === 0) return;

  const totalTime = results.reduce((s, r) => s + r.timeMs, 0);

  let html = `<div class="opt-summary">
    <h3>${done ? t('dft.done') : t('dft.running')}</h3>
    <table>
      <tr>
        <th>${t('dft.colMethod')}</th>
        <th>${t('dft.colRung')}</th>
        <th>${t('dft.colEnergy')}</th>
        <th>${t('dft.colDiff')}</th>
        <th>${t('dft.colTime')}</th>
      </tr>`;

  for (const r of results) {
    const eStr = r.totalEnergy.toFixed(6);
    const dStr = r.method.id === 'HF' ? '\u2014' : `${(r.xcEnergy * 1000).toFixed(2)} mEh`;
    const tStr = r.timeMs < 1000 ? `${r.timeMs.toFixed(0)}ms` : `${(r.timeMs / 1000).toFixed(1)}s`;
    const color = methodColor(r.method.id);
    html += `<tr>
      <td><span style="color:${color};font-weight:600">${r.method.label}</span></td>
      <td style="font-size:0.65rem;color:var(--color-text-dim)">${r.method.rung}</td>
      <td style="font-family:monospace">${eStr}</td>
      <td style="font-family:monospace">${dStr}</td>
      <td style="text-align:right">${tStr}</td>
    </tr>`;
  }

  html += `</table>
    ${done ? `<div style="margin-top:6px;font-size:0.72rem;color:var(--color-text-dim)">${t('dft.totalTime')}: ${(totalTime / 1000).toFixed(1)}s</div>` : ''}
  </div>`;
  el.innerHTML = html;
}

// ── SVG Charts ───────────────────────────────────────────────────────

function renderGraph(): void {
  const container = root.querySelector('#graph-container');
  if (!container || results.length === 0) return;

  container.innerHTML = renderEnergyLevelChart() + renderXCBarChart();
}

/** Chart 1: Energy level staircase — horizontal lines at each method's total energy */
function renderEnergyLevelChart(): string {
  const tc = getThemeColors();

  const width = 560, height = 280;
  const ml = 82, mr = 20, mt = 36, mb = 44;
  const pw = width - ml - mr, ph = height - mt - mb;

  const energies = results.map(r => r.totalEnergy);
  const eMin = Math.min(...energies);
  const eMax = Math.max(...energies);
  const ePad = Math.max((eMax - eMin) * 0.15, 0.002);
  const eLow = eMin - ePad;
  const eHigh = eMax + ePad;

  const toY = (e: number) => mt + ph - ((e - eLow) / (eHigh - eLow)) * ph;

  const n = results.length;
  const barW = pw / METHODS.length;

  let svg = `<svg width="${width}" height="${height}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" style="max-width:100%;">`;

  // Background
  svg += `<rect x="${ml}" y="${mt}" width="${pw}" height="${ph}" fill="${tc.surface}" rx="2"/>`;

  // Y grid (5 divisions)
  for (let i = 0; i <= 5; i++) {
    const yVal = eLow + (eHigh - eLow) * i / 5;
    const yy = toY(yVal);
    svg += `<line x1="${ml}" y1="${yy}" x2="${ml + pw}" y2="${yy}" stroke="${tc.grid}" stroke-width="0.5"/>`;
    svg += `<text x="${ml - 6}" y="${yy + 3}" text-anchor="end" font-size="9" fill="${tc.dim}">${yVal.toFixed(4)}</text>`;
  }

  // Axes
  svg += `<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt + ph}" stroke="${tc.axis}" stroke-width="1"/>`;
  svg += `<line x1="${ml}" y1="${mt + ph}" x2="${ml + pw}" y2="${mt + ph}" stroke="${tc.axis}" stroke-width="1"/>`;

  // Draw bars and energy levels for each method
  for (let i = 0; i < n; i++) {
    const r = results[i];
    const color = methodColor(r.method.id);
    const cx = ml + barW * i + barW / 2;
    const y = toY(r.totalEnergy);

    // Filled bar from HF energy to this method's energy (only for DFT methods)
    if (i > 0) {
      const yHF = toY(results[0].totalEnergy);
      const barX = ml + barW * i + barW * 0.15;
      const bw = barW * 0.7;
      const yTop = Math.min(y, yHF);
      const barH = Math.abs(yHF - y);
      svg += `<rect x="${barX}" y="${yTop}" width="${bw}" height="${barH}" fill="${color}" opacity="0.25" rx="2"/>`;
    }

    // Horizontal energy level line
    const lineX1 = ml + barW * i + barW * 0.1;
    const lineX2 = ml + barW * i + barW * 0.9;
    svg += `<line x1="${lineX1}" y1="${y}" x2="${lineX2}" y2="${y}" stroke="${color}" stroke-width="3" stroke-linecap="round"/>`;

    // Energy annotation
    svg += `<text x="${cx}" y="${y - 8}" text-anchor="middle" font-size="8" fill="${color}" font-weight="600">${r.totalEnergy.toFixed(4)}</text>`;

    // Method label on x-axis
    svg += `<text x="${cx}" y="${mt + ph + 14}" text-anchor="middle" font-size="10" fill="${tc.dim}" font-weight="600">${r.method.label}</text>`;
    // Rung sub-label
    svg += `<text x="${cx}" y="${mt + ph + 26}" text-anchor="middle" font-size="7" fill="${tc.dim}" opacity="0.7">${r.method.rung}</text>`;

    // Connecting dashed line to next level
    if (i > 0) {
      const prevY = toY(results[i - 1].totalEnergy);
      const prevCx = ml + barW * (i - 1) + barW / 2;
      svg += `<line x1="${prevCx}" y1="${prevY}" x2="${cx}" y2="${y}" stroke="${tc.grid}" stroke-width="1" stroke-dasharray="4,3"/>`;
    }
  }

  // Title & Y-axis label
  svg += `<text x="${ml + pw / 2}" y="20" text-anchor="middle" font-size="12" font-weight="600" fill="${tc.titleSvg}">${t('dft.graphTitle')}</text>`;
  svg += `<text x="14" y="${mt + ph / 2}" text-anchor="middle" font-size="10" fill="${tc.dim}" transform="rotate(-90,14,${mt + ph / 2})">${t('dft.yEnergy')}</text>`;

  svg += '</svg>';
  return svg;
}

/** Chart 2: XC energy difference from HF — bar chart */
function renderXCBarChart(): string {
  const tc = getThemeColors();

  // Only show when we have DFT results
  const dftResults = results.filter(r => r.method.id !== 'HF');
  if (dftResults.length < 1) return '';

  const width = 560, height = 180;
  const ml = 82, mr = 20, mt = 28, mb = 40;
  const pw = width - ml - mr, ph = height - mt - mb;

  const n = dftResults.length;
  const barW = pw / n;

  // XC differences in mEh
  const diffs = dftResults.map(r => r.xcEnergy * 1000);
  const maxAbs = Math.max(Math.abs(Math.min(...diffs)), Math.abs(Math.max(...diffs)), 1);
  const scale = maxAbs * 1.2;

  const toY = (mEh: number) => mt + ph / 2 - (mEh / scale) * (ph / 2);

  let svg = `<svg width="${width}" height="${height}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" style="max-width:100%;margin-top:8px;">`;

  // Background
  svg += `<rect x="${ml}" y="${mt}" width="${pw}" height="${ph}" fill="${tc.surface}" rx="2"/>`;

  // Y grid
  const nTicks = 4;
  for (let i = -nTicks; i <= nTicks; i++) {
    const val = scale * i / nTicks;
    const yy = toY(val);
    if (yy < mt || yy > mt + ph) continue;
    svg += `<line x1="${ml}" y1="${yy}" x2="${ml + pw}" y2="${yy}" stroke="${tc.grid}" stroke-width="${i === 0 ? 1 : 0.5}"/>`;
    svg += `<text x="${ml - 6}" y="${yy + 3}" text-anchor="end" font-size="8" fill="${tc.dim}">${val.toFixed(0)}</text>`;
  }

  // Zero line
  const y0 = toY(0);
  svg += `<line x1="${ml}" y1="${y0}" x2="${ml + pw}" y2="${y0}" stroke="${tc.axis}" stroke-width="1"/>`;

  // Axes
  svg += `<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt + ph}" stroke="${tc.axis}" stroke-width="1"/>`;

  // Bars
  for (let i = 0; i < n; i++) {
    const r = dftResults[i];
    const color = methodColor(r.method.id);
    const cx = ml + barW * i + barW / 2;
    const barX = ml + barW * i + barW * 0.2;
    const bw = barW * 0.6;
    const y = toY(diffs[i]);
    const barH = Math.abs(y - y0);
    const barTop = Math.min(y, y0);

    svg += `<rect x="${barX}" y="${barTop}" width="${bw}" height="${barH}" fill="${color}" opacity="0.75" rx="2"/>`;

    // Value label
    const labelY = diffs[i] < 0 ? y + 12 : y - 5;
    svg += `<text x="${cx}" y="${labelY}" text-anchor="middle" font-size="9" fill="${color}" font-weight="700">${diffs[i].toFixed(1)}</text>`;

    // Method label
    svg += `<text x="${cx}" y="${mt + ph + 14}" text-anchor="middle" font-size="10" fill="${tc.dim}" font-weight="600">${r.method.label}</text>`;
  }

  // Title & Y-axis
  svg += `<text x="${ml + pw / 2}" y="16" text-anchor="middle" font-size="11" font-weight="600" fill="${tc.titleSvg}">${t('dft.xcTitle')}</text>`;
  svg += `<text x="14" y="${mt + ph / 2}" text-anchor="middle" font-size="9" fill="${tc.dim}" transform="rotate(-90,14,${mt + ph / 2})">${t('dft.yDiff')}</text>`;

  svg += '</svg>';
  return svg;
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
    .opt-controls { flex: 0 0 340px; min-width: 0; }
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

    .ladder-basis-info { margin-top: 10px; }

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
      overflow: hidden; text-overflow: ellipsis; white-space: nowrap;
    }

    .opt-summary {
      margin-top: 14px; padding: 12px; background: var(--color-surface-alt);
      border-radius: 8px; border: 1px solid var(--color-border);
    }
    .opt-summary h3 {
      font-size: 0.82rem; color: var(--color-converged); margin-bottom: 8px;
    }
    .opt-summary table { width: 100%; font-size: 0.72rem; border-collapse: collapse; }
    .opt-summary th {
      font-size: 0.65rem; color: var(--color-text-dim); text-align: left;
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
