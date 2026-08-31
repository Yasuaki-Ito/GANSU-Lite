/**
 * GANSU Lite — Walsh Diagram demo page.
 * Plots MO orbital energies vs bond angle for AH₂ molecules,
 * visualising why H₂O is bent and BeH₂ is linear.
 */

import { parseXYZ } from './core/parseXYZ';
import { BasisSet } from './core/basisSet';
import { Molecular } from './core/molecular';
import { theorySelectHTML, buildHFOrDFT, HEAVY_ITERATIVE_EXCLUDE, type TheoryChoice } from './ui/theoryControls';
import { RHF } from './core/rhf';
import { initTheme, toggleTheme, isDark, getThemeColors } from './ui/theme';
import { t, initLang, toggleLang } from './ui/i18n';
import { renderHeader } from './ui/nav';

// ── Scenario definitions ─────────────────────────────────────────────

interface WalshScenario {
  id: string;
  labelKey: string;
  descKey: string;
  charge: number;
  mult: number;
  angleMin: number;
  angleMax: number;
  steps: number;
  nCore: number;       // core MOs to skip in orbital plot
  buildXYZ: (angleDeg: number) => string;
}

function buildAH2(center: string, R: number, angleDeg: number): string {
  const rad = angleDeg * Math.PI / 180;
  const h = rad / 2;
  const hx = R * Math.sin(h);
  const hz = R * Math.cos(h);
  return `3\n${center}H2\n${center}  0.0  0.0  0.0\nH  ${hx.toFixed(6)}  0.0  ${hz.toFixed(6)}\nH  ${(-hx).toFixed(6)}  0.0  ${hz.toFixed(6)}`;
}

const SCENARIOS: WalshScenario[] = [
  {
    id: 'h2o',
    labelKey: 'walsh.scenH2O', descKey: 'walsh.descH2O',
    charge: 0, mult: 1,
    angleMin: 80, angleMax: 180, steps: 21,
    nCore: 1,
    buildXYZ: (a) => buildAH2('O', 0.96, a),
  },
  {
    id: 'beh2',
    labelKey: 'walsh.scenBeH2', descKey: 'walsh.descBeH2',
    charge: 0, mult: 1,
    angleMin: 80, angleMax: 180, steps: 21,
    nCore: 1,
    buildXYZ: (a) => buildAH2('Be', 1.334, a),
  },
  {
    id: 'ch2',
    labelKey: 'walsh.scenCH2', descKey: 'walsh.descCH2',
    charge: 0, mult: 1,
    angleMin: 80, angleMax: 180, steps: 21,
    nCore: 1,
    buildXYZ: (a) => buildAH2('C', 1.11, a),
  },
];

// ── Types ────────────────────────────────────────────────────────────

interface WalshPoint {
  angle: number;
  orbitalEnergies: number[];
  totalEnergy: number;
}

// ── State ────────────────────────────────────────────────────────────

let selectedScen = SCENARIOS[0];
let data: WalshPoint[] = [];
let nOccupied = 0;
let nBasis = 0;
let running = false;
let stopRequested = false;
let done = false;
let elapsedMs = 0;
let theoryChoice: TheoryChoice = 'HF';
let selectedIdx = 0;

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
      ${renderHeader('walsh')}

      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${t('walsh.molecule')}</h2>
          <div class="opt-scenario-grid" id="scen-grid"></div>

          <div class="theory-row" style="margin:12px 0 8px;display:flex;align-items:center;gap:8px;font-size:0.9rem;">
            <span>Theory:</span>${theorySelectHTML('theory-sel', theoryChoice, '', HEAVY_ITERATIVE_EXCLUDE)}
          </div>

          <button id="run-btn" class="opt-run-btn" ${running ? 'disabled' : ''}>
            ${running ? t('walsh.running') : t('walsh.run')}
          </button>
          ${running ? `<button id="stop-btn" class="opt-stop-btn">${t('walsh.stop')}</button>` : ''}

          <div id="progress-area"></div>
          <div id="summary-area"></div>
        </div>

        <div class="opt-panel opt-graph-panel">
          <div id="graph-container">
            ${!done && !running
              ? `<p class="opt-hint">${t('walsh.waiting')}</p>`
              : ''}
          </div>
          <div id="slider-area"></div>
          <div id="mol-vis"></div>
        </div>
      </div>
    </div>`;

  injectStyles();

  // Scenario cards
  const grid = root.querySelector('#scen-grid')!;
  const row = document.createElement('div');
  row.className = 'opt-category-row';
  for (const scen of SCENARIOS) {
    const card = document.createElement('div');
    card.className = 'opt-scenario-card' + (scen.id === selectedScen.id ? ' selected' : '');
    card.innerHTML = `<strong>${t(scen.labelKey)}</strong><span class="conv-desc">${t(scen.descKey)}</span>`;
    card.addEventListener('click', () => {
      if (running) return;
      selectedScen = scen;
      data = [];
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
  root.querySelector('#run-btn')!.addEventListener('click', () => { if (!running) runWalsh(); });
  root.querySelector('#stop-btn')?.addEventListener('click', () => { stopRequested = true; });

  // Render graph & summary if data exists
  if (data.length > 0) {
    renderGraph();
    renderSummary();
    if (done) {
      setupSlider();
      renderMolVis();
    }
  }
}

// ── Computation ──────────────────────────────────────────────────────

async function runWalsh(): Promise<void> {
  running = true;
  stopRequested = false;
  data = [];
  done = false;
  nOccupied = 0;
  nBasis = 0;
  render();

  const scen = selectedScen;
  const startTime = performance.now();
  const basis = await loadBasis('STO-3G');
  const angleStep = (scen.angleMax - scen.angleMin) / (scen.steps - 1);

  for (let i = 0; i < scen.steps; i++) {
    if (stopRequested) break;

    const angle = scen.angleMin + i * angleStep;
    updateProgress(angle, i, scen.steps);

    try {
      const xyzStr = scen.buildXYZ(angle);
      const atoms = parseXYZ(xyzStr);
      const betaToAlpha = Math.floor((scen.mult - 1) / 2);
      const mol = new Molecular(atoms, basis, scen.charge, betaToAlpha);

      if (i === 0) {
        nBasis = mol.numBasis;
        nOccupied = mol.numAlphaSpins;
      }

      const rhf = await buildHFOrDFT(mol, basis, theoryChoice) as RHF;
      const totalEnergy = await rhf.solve({ eriBackend: 'js' });
      const orbEnergies = Array.from(rhf.orbitalEnergies);

      data.push({ angle, orbitalEnergies: orbEnergies, totalEnergy });

      renderGraph();
    } catch (err) {
      console.error(`Error at angle ${angle}°:`, err);
    }

    await new Promise<void>(r => setTimeout(r, 0));
  }

  elapsedMs = performance.now() - startTime;
  running = false;
  done = true;
  render();
}

// ── Progress ─────────────────────────────────────────────────────────

function updateProgress(angle: number, current: number, total: number): void {
  const area = root.querySelector('#progress-area');
  if (!area) return;
  const pct = (current / total * 100).toFixed(0);
  area.innerHTML = `
    <div class="opt-progress"><div class="opt-progress-bar" style="width:${pct}%"></div></div>
    <div class="opt-progress-text">${angle.toFixed(0)}\u00b0 (${current + 1}/${total})</div>`;
}

// ── Summary ──────────────────────────────────────────────────────────

function renderSummary(): void {
  const el = root.querySelector('#summary-area');
  if (!el || data.length === 0) return;

  let minIdx = 0;
  for (let i = 1; i < data.length; i++) {
    if (data[i].totalEnergy < data[minIdx].totalEnergy) minIdx = i;
  }
  const nValence = nBasis - selectedScen.nCore;
  const nOccVal = nOccupied - selectedScen.nCore;
  const totalTime = (elapsedMs / 1000).toFixed(1);

  let html = `<div class="opt-summary">
    <h3>${done ? t('walsh.done') : t('walsh.running')}</h3>
    <table>
      <tr><td style="color:var(--color-text-dim)">${t('walsh.basis')}</td><td>STO-3G (M=${nBasis})</td></tr>
      <tr><td style="color:var(--color-text-dim)">${t('walsh.nOcc')}</td><td>${nOccVal} / ${nValence}</td></tr>
      <tr><td style="color:var(--color-text-dim)">${t('walsh.minAngle')}</td><td>${data[minIdx].angle.toFixed(1)}\u00b0</td></tr>
      <tr><td style="color:var(--color-text-dim)">${t('walsh.minEnergy')}</td><td style="font-family:monospace">${data[minIdx].totalEnergy.toFixed(6)} Eh</td></tr>
      ${done ? `<tr><td style="color:var(--color-text-dim)">${t('walsh.time')}</td><td>${totalTime}s</td></tr>` : ''}
    </table>
  </div>`;
  el.innerHTML = html;
}

// ── SVG Graph ────────────────────────────────────────────────────────

// MO color palettes (light / dark)
const MO_COLORS_LIGHT = [
  '#d62728', '#1f77b4', '#2ca02c', '#9467bd', '#ff7f0e',
  '#8c564b', '#e377c2', '#17becf', '#bcbd22', '#7f7f7f',
];
const MO_COLORS_DARK = [
  '#ff6b6b', '#4dabf7', '#51cf66', '#cc5de8', '#ffa94d',
  '#d4a373', '#f783ac', '#38d9f5', '#d9e363', '#adb5bd',
];

function renderGraph(): void {
  const container = root.querySelector('#graph-container');
  if (!container || data.length === 0) return;

  const tc = getThemeColors();
  const colors = isDark() ? MO_COLORS_DARK : MO_COLORS_LIGHT;
  const scen = selectedScen;
  const nCore = scen.nCore;
  const nTotal = nBasis;

  // ── Walsh diagram (orbital energies) ──
  const width = 560, height = 380;
  const ml = 72, mr = 52, mt = 36, mb = 44;
  const pw = width - ml - mr, ph = height - mt - mb;

  // Collect valence orbital energies to determine Y range
  const allEnergies: number[] = [];
  for (const pt of data) {
    for (let m = nCore; m < pt.orbitalEnergies.length; m++) {
      allEnergies.push(pt.orbitalEnergies[m]);
    }
  }
  if (allEnergies.length === 0) return;

  const yMin = Math.min(...allEnergies);
  const yMax = Math.max(...allEnergies);
  const yPad = Math.max((yMax - yMin) * 0.08, 0.05);
  const yLow = yMin - yPad;
  const yHigh = yMax + yPad;

  const toX = (angle: number) => ml + ((angle - scen.angleMin) / (scen.angleMax - scen.angleMin)) * pw;
  const toY = (e: number) => mt + ph - ((e - yLow) / (yHigh - yLow)) * ph;

  let svg = `<svg width="${width}" height="${height}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" style="max-width:100%;">`;

  // Background
  svg += `<rect x="${ml}" y="${mt}" width="${pw}" height="${ph}" fill="${tc.surface}" rx="2"/>`;

  // Y grid (5 divisions)
  for (let i = 0; i <= 5; i++) {
    const yVal = yLow + (yHigh - yLow) * i / 5;
    const yy = toY(yVal);
    svg += `<line x1="${ml}" y1="${yy}" x2="${ml + pw}" y2="${yy}" stroke="${tc.grid}" stroke-width="0.5"/>`;
    svg += `<text x="${ml - 6}" y="${yy + 3}" text-anchor="end" font-size="9" fill="${tc.dim}">${yVal.toFixed(2)}</text>`;
  }

  // X grid + tick labels
  for (let angle = scen.angleMin; angle <= scen.angleMax + 0.1; angle += 20) {
    const a = Math.round(angle);
    const xx = toX(a);
    svg += `<line x1="${xx}" y1="${mt}" x2="${xx}" y2="${mt + ph}" stroke="${tc.grid}" stroke-width="0.3"/>`;
    svg += `<line x1="${xx}" y1="${mt + ph}" x2="${xx}" y2="${mt + ph + 4}" stroke="${tc.axis}" stroke-width="1"/>`;
    svg += `<text x="${xx}" y="${mt + ph + 16}" text-anchor="middle" font-size="9" fill="${tc.dim}">${a}\u00b0</text>`;
  }

  // Axes
  svg += `<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt + ph}" stroke="${tc.axis}" stroke-width="1"/>`;
  svg += `<line x1="${ml}" y1="${mt + ph}" x2="${ml + pw}" y2="${mt + ph}" stroke="${tc.axis}" stroke-width="1"/>`;

  // Draw MO lines
  for (let m = nCore; m < nTotal; m++) {
    const moIdx = m - nCore;
    const isOcc = m < nOccupied;
    const color = colors[moIdx % colors.length];
    const strokeW = isOcc ? 2 : 1.2;
    const dash = isOcc ? '' : ' stroke-dasharray="6,3"';

    // Build path
    let path = '';
    for (let i = 0; i < data.length; i++) {
      if (m >= data[i].orbitalEnergies.length) continue;
      const px = toX(data[i].angle);
      const py = toY(data[i].orbitalEnergies[m]);
      path += i === 0 ? `M${px.toFixed(1)},${py.toFixed(1)}` : ` L${px.toFixed(1)},${py.toFixed(1)}`;
    }
    if (path) {
      svg += `<path d="${path}" fill="none" stroke="${color}" stroke-width="${strokeW}"${dash}/>`;
    }

    // Label at right edge
    if (data.length > 0) {
      const lastPt = data[data.length - 1];
      if (m < lastPt.orbitalEnergies.length) {
        const lx = toX(lastPt.angle) + 4;
        const ly = toY(lastPt.orbitalEnergies[m]) + 3;
        const label = isOcc ? `${moIdx + 1}` : `${moIdx + 1}*`;
        svg += `<text x="${lx}" y="${ly}" font-size="8" fill="${color}" font-weight="${isOcc ? '600' : '400'}">${label}</text>`;
      }
    }
  }

  // Legend (occupied = solid, virtual = dashed)
  const legY = mt + 10;
  svg += `<line x1="${ml + 8}" y1="${legY}" x2="${ml + 28}" y2="${legY}" stroke="${tc.dim}" stroke-width="2"/>`;
  svg += `<text x="${ml + 32}" y="${legY + 3}" font-size="7.5" fill="${tc.dim}">${t('walsh.occupied')}</text>`;
  svg += `<line x1="${ml + 8}" y1="${legY + 14}" x2="${ml + 28}" y2="${legY + 14}" stroke="${tc.dim}" stroke-width="1.2" stroke-dasharray="6,3"/>`;
  svg += `<text x="${ml + 32}" y="${legY + 17}" font-size="7.5" fill="${tc.dim}">${t('walsh.virtual')}</text>`;

  // Highlight selected angle (vertical line + circles on each MO)
  if (done && selectedIdx >= 0 && selectedIdx < data.length) {
    const selPt = data[selectedIdx];
    const selX = toX(selPt.angle);
    const hlColor = isDark() ? '#ffd700' : '#cc8800';
    svg += `<line x1="${selX}" y1="${mt}" x2="${selX}" y2="${mt + ph}" stroke="${hlColor}" stroke-width="1" stroke-dasharray="3,2" opacity="0.6"/>`;
    for (let m = nCore; m < nTotal; m++) {
      if (m >= selPt.orbitalEnergies.length) continue;
      const moIdx = m - nCore;
      const color = colors[moIdx % colors.length];
      const py = toY(selPt.orbitalEnergies[m]);
      svg += `<circle cx="${selX}" cy="${py}" r="4.5" fill="${color}" stroke="#fff" stroke-width="1"/>`;
    }
  }

  // Title & axis labels
  svg += `<text x="${ml + pw / 2}" y="20" text-anchor="middle" font-size="12" font-weight="600" fill="${tc.titleSvg}">${t('walsh.graphTitle')}</text>`;
  svg += `<text x="${ml + pw / 2}" y="${height - 4}" text-anchor="middle" font-size="10" fill="${tc.dim}">${t('walsh.xAngle')}</text>`;
  svg += `<text x="14" y="${mt + ph / 2}" text-anchor="middle" font-size="10" fill="${tc.dim}" transform="rotate(-90,14,${mt + ph / 2})">${t('walsh.yOrbEnergy')}</text>`;

  svg += '</svg>';

  // ── Total energy subplot ──
  const h2 = 140;
  const mt2 = 16, mb2 = 36;
  const ph2 = h2 - mt2 - mb2;

  const totalEnergies = data.map(d => d.totalEnergy);
  const teMin = Math.min(...totalEnergies);
  const teMax = Math.max(...totalEnergies);
  const tePad = Math.max((teMax - teMin) * 0.15, 0.002);
  const teLow = teMin - tePad;
  const teHigh = teMax + tePad;

  const toY2 = (e: number) => mt2 + ph2 - ((e - teLow) / (teHigh - teLow)) * ph2;

  let svg2 = `<svg width="${width}" height="${h2}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${h2}" style="max-width:100%;margin-top:8px;">`;

  // Background
  svg2 += `<rect x="${ml}" y="${mt2}" width="${pw}" height="${ph2}" fill="${tc.surface}" rx="2"/>`;

  // Y grid (3 divisions)
  for (let i = 0; i <= 3; i++) {
    const yVal = teLow + (teHigh - teLow) * i / 3;
    const yy = toY2(yVal);
    svg2 += `<line x1="${ml}" y1="${yy}" x2="${ml + pw}" y2="${yy}" stroke="${tc.grid}" stroke-width="0.5"/>`;
    svg2 += `<text x="${ml - 6}" y="${yy + 3}" text-anchor="end" font-size="8" fill="${tc.dim}">${yVal.toFixed(4)}</text>`;
  }

  // Axes
  svg2 += `<line x1="${ml}" y1="${mt2}" x2="${ml}" y2="${mt2 + ph2}" stroke="${tc.axis}" stroke-width="1"/>`;
  svg2 += `<line x1="${ml}" y1="${mt2 + ph2}" x2="${ml + pw}" y2="${mt2 + ph2}" stroke="${tc.axis}" stroke-width="1"/>`;

  // X tick labels
  for (let angle = scen.angleMin; angle <= scen.angleMax + 0.1; angle += 20) {
    const a = Math.round(angle);
    const xx = toX(a);
    svg2 += `<line x1="${xx}" y1="${mt2 + ph2}" x2="${xx}" y2="${mt2 + ph2 + 4}" stroke="${tc.axis}" stroke-width="1"/>`;
    svg2 += `<text x="${xx}" y="${mt2 + ph2 + 14}" text-anchor="middle" font-size="8" fill="${tc.dim}">${a}\u00b0</text>`;
  }

  // Total energy line
  const teColor = isDark() ? '#00d4ff' : '#0077cc';
  if (data.length >= 2) {
    let path = '';
    for (let i = 0; i < data.length; i++) {
      const px = toX(data[i].angle);
      const py = toY2(data[i].totalEnergy);
      path += i === 0 ? `M${px.toFixed(1)},${py.toFixed(1)}` : ` L${px.toFixed(1)},${py.toFixed(1)}`;
    }
    svg2 += `<path d="${path}" fill="none" stroke="${teColor}" stroke-width="2"/>`;
  }
  for (const pt of data) {
    svg2 += `<circle cx="${toX(pt.angle).toFixed(1)}" cy="${toY2(pt.totalEnergy).toFixed(1)}" r="2.5" fill="${teColor}"/>`;
  }

  // Minimum marker
  if (data.length > 0) {
    let minIdx = 0;
    for (let i = 1; i < data.length; i++) {
      if (data[i].totalEnergy < data[minIdx].totalEnergy) minIdx = i;
    }
    const mx = toX(data[minIdx].angle);
    const my = toY2(data[minIdx].totalEnergy);
    svg2 += `<circle cx="${mx}" cy="${my}" r="5" fill="none" stroke="${teColor}" stroke-width="1.5"/>`;
    svg2 += `<text x="${mx}" y="${my - 8}" text-anchor="middle" font-size="8" fill="${teColor}">${data[minIdx].angle.toFixed(0)}\u00b0</text>`;
  }

  // Highlight selected angle on total energy plot
  if (done && selectedIdx >= 0 && selectedIdx < data.length) {
    const selPt = data[selectedIdx];
    const selX = toX(selPt.angle);
    const selY = toY2(selPt.totalEnergy);
    const hlColor = isDark() ? '#ffd700' : '#cc8800';
    svg2 += `<line x1="${selX}" y1="${mt2}" x2="${selX}" y2="${mt2 + ph2}" stroke="${hlColor}" stroke-width="1" stroke-dasharray="3,2" opacity="0.6"/>`;
    svg2 += `<circle cx="${selX}" cy="${selY}" r="5" fill="${teColor}" stroke="${hlColor}" stroke-width="2"/>`;
  }

  // Labels
  svg2 += `<text x="${ml + pw / 2}" y="10" text-anchor="middle" font-size="10" font-weight="600" fill="${tc.titleSvg}">${t('walsh.totalEnergy')}</text>`;
  svg2 += `<text x="${ml + pw / 2}" y="${h2 - 2}" text-anchor="middle" font-size="9" fill="${tc.dim}">${t('walsh.xAngle')}</text>`;
  svg2 += `<text x="14" y="${mt2 + ph2 / 2}" text-anchor="middle" font-size="9" fill="${tc.dim}" transform="rotate(-90,14,${mt2 + ph2 / 2})">${t('walsh.yTotalEnergy')}</text>`;

  svg2 += '</svg>';

  container.innerHTML = svg + svg2;
}

// ── Slider ───────────────────────────────────────────────────────────

function setupSlider(): void {
  const area = root.querySelector('#slider-area');
  if (!area || data.length === 0) return;
  if (selectedIdx >= data.length) selectedIdx = 0;

  area.innerHTML = `
    <div class="walsh-slider">
      <input type="range" id="angle-slider" min="0" max="${data.length - 1}" value="${selectedIdx}" />
      <div class="walsh-slider-label">${data[selectedIdx].angle.toFixed(0)}\u00b0</div>
    </div>`;

  root.querySelector('#angle-slider')?.addEventListener('input', (e) => {
    selectedIdx = parseInt((e.target as HTMLInputElement).value);
    const label = root.querySelector('.walsh-slider-label');
    if (label) label.textContent = `${data[selectedIdx].angle.toFixed(0)}\u00b0`;
    renderGraph();
    renderMolVis();
  });
}

// ── Molecule SVG ─────────────────────────────────────────────────────

const ELEM_COLORS: Record<string, string> = {
  H: '#999', Be: '#6C0', C: '#555', O: '#F22',
};

function renderMolVis(): void {
  const el = root.querySelector('#mol-vis');
  if (!el || !done || data.length === 0) return;

  const tc = getThemeColors();
  const scen = selectedScen;
  const angle = data[selectedIdx].angle;

  const W = 220, H = 180;
  const margin = 30;

  // Parse atom coordinates
  const xyzStr = scen.buildXYZ(angle);
  const lines = xyzStr.split('\n');
  const atoms: { sym: string; x: number; z: number }[] = [];
  for (let i = 2; i < lines.length; i++) {
    const parts = lines[i].trim().split(/\s+/);
    if (parts.length >= 4) {
      atoms.push({ sym: parts[0], x: parseFloat(parts[1]), z: parseFloat(parts[3]) });
    }
  }

  // Fixed scale from extreme angles so molecule doesn't jump around
  let uMin = 0, uMax = 0, vMin = 0, vMax = 0;
  for (const a of [scen.angleMin, scen.angleMax]) {
    const str = scen.buildXYZ(a);
    const ls = str.split('\n');
    for (let j = 2; j < ls.length; j++) {
      const p = ls[j].trim().split(/\s+/);
      const ux = parseFloat(p[1]), vz = parseFloat(p[3]);
      uMin = Math.min(uMin, ux); uMax = Math.max(uMax, ux);
      vMin = Math.min(vMin, vz); vMax = Math.max(vMax, vz);
    }
  }
  const uRange = uMax - uMin || 0.01;
  const vRange = vMax - vMin || 0.01;
  const scale = Math.min((W - 2 * margin) / uRange, (H - 2 * margin) / vRange);
  const uMid = (uMin + uMax) / 2;
  const vMid = (vMin + vMax) / 2;

  const sx = (u: number) => W / 2 + (u - uMid) * scale;
  const sy = (v: number) => H / 2 - (v - vMid) * scale;

  let svg = `<svg width="${W}" height="${H}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" style="max-width:100%;">`;

  // Bonds
  for (let i = 0; i < atoms.length; i++) {
    for (let j = i + 1; j < atoms.length; j++) {
      const dx = atoms[j].x - atoms[i].x, dz = atoms[j].z - atoms[i].z;
      if (Math.sqrt(dx * dx + dz * dz) < 2.5) {
        svg += `<line x1="${sx(atoms[i].x)}" y1="${sy(atoms[i].z)}" x2="${sx(atoms[j].x)}" y2="${sy(atoms[j].z)}" stroke="${tc.grid}" stroke-width="4" stroke-linecap="round"/>`;
      }
    }
  }

  // Atoms
  for (const a of atoms) {
    const ax = sx(a.x), ay = sy(a.z);
    const r = a.sym === 'H' ? 14 : 18;
    svg += `<circle cx="${ax}" cy="${ay}" r="${r}" fill="${ELEM_COLORS[a.sym] || '#888'}" stroke="${tc.axis}" stroke-width="1.5"/>`;
    svg += `<text x="${ax}" y="${ay + 4}" text-anchor="middle" font-size="11" font-weight="bold" fill="#fff">${a.sym}</text>`;
  }

  // Angle arc annotation
  if (atoms.length >= 3) {
    const cx = sx(atoms[0].x), cy = sy(atoms[0].z);
    const arcR = 28;
    // Angles from center to each H in SVG space
    const a1 = Math.atan2(sy(atoms[1].z) - cy, sx(atoms[1].x) - cx);
    const a2 = Math.atan2(sy(atoms[2].z) - cy, sx(atoms[2].x) - cx);
    const ex1 = cx + arcR * Math.cos(a1), ey1 = cy + arcR * Math.sin(a1);
    const ex2 = cx + arcR * Math.cos(a2), ey2 = cy + arcR * Math.sin(a2);

    const ac = isDark() ? '#ffd700' : '#cc8800';
    svg += `<path d="M ${ex1.toFixed(1)} ${ey1.toFixed(1)} A ${arcR} ${arcR} 0 0 0 ${ex2.toFixed(1)} ${ey2.toFixed(1)}" fill="none" stroke="${ac}" stroke-width="1.5"/>`;

    // Label at arc midpoint
    const aMid = (a1 + a2) / 2;
    const lx = cx + (arcR + 14) * Math.cos(aMid);
    const ly = cy + (arcR + 14) * Math.sin(aMid);
    svg += `<text x="${lx.toFixed(1)}" y="${ly.toFixed(1)}" text-anchor="middle" font-size="11" font-weight="600" fill="${ac}">${angle.toFixed(0)}\u00b0</text>`;
  }

  svg += '</svg>';
  el.innerHTML = svg;
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
    #mol-vis { width: 100%; text-align: center; }
    .opt-hint { color: var(--color-text-dim); font-size: 0.85rem; padding: 60px 20px; }

    .walsh-slider {
      width: 100%; max-width: 560px; margin: 8px auto 0;
      display: flex; align-items: center; gap: 10px;
    }
    .walsh-slider input[type=range] { flex: 1; cursor: pointer; }
    .walsh-slider-label {
      font-size: 0.85rem; font-weight: 600; min-width: 40px; text-align: center;
      color: var(--color-text);
    }

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
