/**
 * GANSU Lite — Mulliken Charge Map demo page.
 * Computes Mulliken charges for various molecules and displays
 * color-mapped molecular structures to visualise electronegativity effects.
 */

import { parseXYZ } from './core/parseXYZ';
import { BasisSet } from './core/basisSet';
import { Molecular } from './core/molecular';
import { theorySelectHTML, buildHFOrDFT, type TheoryChoice } from './ui/theoryControls';
import { RHF } from './core/rhf';
import { computeMullikenCharges, computeDipoleMoment } from './core/properties';
import type { DipoleMoment } from './core/properties';
import { atomicNumberToElementName, BOHR_TO_ANGSTROM } from './core/constants';
import { initTheme, toggleTheme, getThemeColors } from './ui/theme';
import { t, initLang, toggleLang } from './ui/i18n';
import { renderHeader } from './ui/nav';
import { renderChargeViewer } from './ui/moleculeViewer3D';

// ── Molecule definitions ─────────────────────────────────────────────

interface MolDef {
  id: string;
  labelKey: string;
  descKey: string;
  category: 'homo' | 'hetero' | 'poly' | 'ion';
  charge: number;
  mult: number;
  xyz: string;
}

const MOLECULES: MolDef[] = [
  // ── Homonuclear (zero charge difference) ──
  {
    id: 'h2', labelKey: 'chg.scenH2', descKey: 'chg.descH2',
    category: 'homo', charge: 0, mult: 1,
    xyz: '2\nH2\nH  0.0  0.0  0.0\nH  0.0  0.0  0.740000',
  },
  {
    id: 'n2', labelKey: 'chg.scenN2', descKey: 'chg.descN2',
    category: 'homo', charge: 0, mult: 1,
    xyz: '2\nN2\nN  0.0  0.0  0.0\nN  0.0  0.0  1.098000',
  },
  {
    id: 'f2', labelKey: 'chg.scenF2', descKey: 'chg.descF2',
    category: 'homo', charge: 0, mult: 1,
    xyz: '2\nF2\nF  0.0  0.0  0.0\nF  0.0  0.0  1.412000',
  },

  // ── Heteronuclear diatomic ──
  {
    id: 'lih', labelKey: 'chg.scenLiH', descKey: 'chg.descLiH',
    category: 'hetero', charge: 0, mult: 1,
    xyz: '2\nLiH\nLi  0.0  0.0  0.0\nH   0.0  0.0  1.596000',
  },
  {
    id: 'lif', labelKey: 'chg.scenLiF', descKey: 'chg.descLiF',
    category: 'hetero', charge: 0, mult: 1,
    xyz: '2\nLiF\nLi  0.0  0.0  0.0\nF   0.0  0.0  1.564000',
  },
  {
    id: 'hf', labelKey: 'chg.scenHF', descKey: 'chg.descHF',
    category: 'hetero', charge: 0, mult: 1,
    xyz: '2\nHF\nH  0.0  0.0  0.0\nF  0.0  0.0  0.917000',
  },
  {
    id: 'co', labelKey: 'chg.scenCO', descKey: 'chg.descCO',
    category: 'hetero', charge: 0, mult: 1,
    xyz: '2\nCO\nC  0.0  0.0  0.0\nO  0.0  0.0  1.128000',
  },

  // ── Polyatomic ──
  {
    id: 'h2o', labelKey: 'chg.scenH2O', descKey: 'chg.descH2O',
    category: 'poly', charge: 0, mult: 1,
    xyz: (() => {
      const R = 0.96, ang = 104 * Math.PI / 180, h = ang / 2;
      const hx = R * Math.sin(h), hz = R * Math.cos(h);
      return `3\nH2O\nO  0.0  0.0  0.0\nH  ${hx.toFixed(6)}  0.0  ${hz.toFixed(6)}\nH  ${(-hx).toFixed(6)}  0.0  ${hz.toFixed(6)}`;
    })(),
  },
  {
    id: 'nh3', labelKey: 'chg.scenNH3', descKey: 'chg.descNH3',
    category: 'poly', charge: 0, mult: 1,
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
    id: 'ch4', labelKey: 'chg.scenCH4', descKey: 'chg.descCH4',
    category: 'poly', charge: 0, mult: 1,
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
    id: 'bh3', labelKey: 'chg.scenBH3', descKey: 'chg.descBH3',
    category: 'poly', charge: 0, mult: 1,
    xyz: (() => {
      const R = 1.19;
      const angles = [0, 120, 240].map(a => a * Math.PI / 180);
      const hs = angles.map(a => `H  ${(R * Math.cos(a)).toFixed(6)}  ${(R * Math.sin(a)).toFixed(6)}  0.0`);
      return `4\nBH3\nB  0.0  0.0  0.0\n${hs.join('\n')}`;
    })(),
  },

  // ── Ionic ──
  {
    id: 'heh', labelKey: 'chg.scenHeH', descKey: 'chg.descHeH',
    category: 'ion', charge: 1, mult: 1,
    xyz: '2\nHeH+\nHe  0.0  0.0  0.0\nH   0.0  0.0  0.774000',
  },
  {
    id: 'bh4', labelKey: 'chg.scenBH4', descKey: 'chg.descBH4',
    category: 'ion', charge: -1, mult: 1,
    xyz: (() => {
      const R = 1.255;
      const a = R / Math.sqrt(3);
      return [
        '5', 'BH4-',
        'B   0.000000  0.000000  0.000000',
        `H   ${a.toFixed(6)}  ${a.toFixed(6)}  ${a.toFixed(6)}`,
        `H   ${a.toFixed(6)}  ${(-a).toFixed(6)}  ${(-a).toFixed(6)}`,
        `H   ${(-a).toFixed(6)}  ${a.toFixed(6)}  ${(-a).toFixed(6)}`,
        `H   ${(-a).toFixed(6)}  ${(-a).toFixed(6)}  ${a.toFixed(6)}`,
      ].join('\n');
    })(),
  },
];

const CATEGORY_LABELS: Record<string, string> = {
  homo: 'chg.catHomo',
  hetero: 'chg.catHetero',
  poly: 'chg.catPoly',
  ion: 'chg.catIon',
};

// ── Types ────────────────────────────────────────────────────────────

interface ChargeResult {
  molId: string;
  charges: number[];
  dipole: DipoleMoment;
  atomSymbols: string[];
  atomCoords: { x: number; y: number; z: number }[];  // Ångström
  timeMs: number;
}

// ── State ────────────────────────────────────────────────────────────

let selectedMol = MOLECULES[0];
let results = new Map<string, ChargeResult>();
let running = false;
let stopRequested = false;
let runAllMode = false;
let runAllProgress = 0;
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
  const hasResult = results.has(selectedMol.id);

  root.innerHTML = `
    <div class="opt-page">
      ${renderHeader('charges')}

      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${t('chg.molecule')}</h2>
          <div class="opt-scenario-grid" id="scen-grid"></div>

          <div class="ladder-basis-info">
            <span style="color:var(--color-text-dim);font-size:0.72rem">${t('chg.basis')}: STO-3G</span>
          </div>

          <div class="theory-row" style="margin:8px 0;display:flex;align-items:center;gap:8px;font-size:0.9rem;">
            <span>Theory:</span>${theorySelectHTML('theory-sel', theoryChoice)}
          </div>

          <button id="run-btn" class="opt-run-btn" ${running ? 'disabled' : ''}>
            ${running ? t('chg.running') : (hasResult ? t('chg.rerun') : t('chg.run'))}
          </button>
          <button id="run-all-btn" class="opt-run-all-btn" ${running ? 'disabled' : ''}>
            ${t('chg.runAll')}
          </button>
          ${running ? `<button id="stop-btn" class="opt-stop-btn">${t('chg.stop')}</button>` : ''}

          <div id="progress-area"></div>
          <div id="summary-area"></div>
        </div>

        <div class="opt-panel opt-graph-panel">
          <div id="mol-vis">
            ${!hasResult && !running
              ? `<p class="opt-hint">${t('chg.waiting')}</p>`
              : ''}
          </div>
          <div id="compare-chart"></div>
        </div>
      </div>
    </div>`;

  injectStyles();

  // Scenario cards grouped by category
  const grid = root.querySelector('#scen-grid')!;
  const categories = ['homo', 'hetero', 'poly', 'ion'];
  for (const cat of categories) {
    const mols = MOLECULES.filter(m => m.category === cat);
    if (mols.length === 0) continue;

    const catLabel = document.createElement('div');
    catLabel.className = 'chg-cat-label';
    catLabel.textContent = t(CATEGORY_LABELS[cat]);
    grid.appendChild(catLabel);

    const row = document.createElement('div');
    row.className = 'opt-category-row';
    for (const mol of mols) {
      const card = document.createElement('div');
      const hasR = results.has(mol.id);
      card.className = 'opt-scenario-card'
        + (mol.id === selectedMol.id ? ' selected' : '')
        + (hasR ? ' computed' : '');
      card.innerHTML = `<strong>${t(mol.labelKey)}</strong><span class="conv-desc">${t(mol.descKey)}</span>`;
      card.addEventListener('click', () => {
        if (running) return;
        selectedMol = mol;
        render();
      });
      row.appendChild(card);
    }
    grid.appendChild(row);
  }

  // Event listeners
  root.querySelector('#nav-theme')!.addEventListener('click', () => { toggleTheme(); render(); });
  root.querySelector('#nav-lang')!.addEventListener('click', () => { toggleLang(); render(); });
  const theorySel = root.querySelector<HTMLSelectElement>('#theory-sel');
  if (theorySel) theorySel.addEventListener('change', () => { theoryChoice = theorySel.value as TheoryChoice; });
  root.querySelector('#run-btn')!.addEventListener('click', () => { if (!running) runSingle(selectedMol); });
  root.querySelector('#run-all-btn')!.addEventListener('click', () => { if (!running) runAll(); });
  root.querySelector('#stop-btn')?.addEventListener('click', () => { stopRequested = true; });

  if (hasResult) {
    renderMolVis(results.get(selectedMol.id)!);
    renderSummary(results.get(selectedMol.id)!);
  }

  if (results.size > 1) {
    renderCompareChart();
  }
}

// ── Computation ──────────────────────────────────────────────────────

async function computeMol(mol: MolDef): Promise<ChargeResult> {
  const basis = await loadBasis('STO-3G');
  const atoms = parseXYZ(mol.xyz);
  const betaToAlpha = Math.floor((mol.mult - 1) / 2);
  const molecular = new Molecular(atoms, basis, mol.charge, betaToAlpha);

  const t0 = performance.now();
  const rhf = await buildHFOrDFT(molecular, basis, theoryChoice) as RHF;
  await rhf.solve({ eriBackend: 'js' });
  const timeMs = performance.now() - t0;

  const charges = computeMullikenCharges(rhf.density, rhf.overlap, molecular.atoms, molecular.atomToBasisRange);
  const dipole = computeDipoleMoment(rhf.density, molecular.atoms, molecular.primitiveShells, molecular.cgtoNormalizationFactors, molecular.numBasis);

  const atomSymbols: string[] = [];
  const atomCoords: { x: number; y: number; z: number }[] = [];
  for (const a of molecular.atoms) {
    atomSymbols.push(atomicNumberToElementName(a.atomicNumber));
    atomCoords.push({
      x: a.coordinate.x * BOHR_TO_ANGSTROM,
      y: a.coordinate.y * BOHR_TO_ANGSTROM,
      z: a.coordinate.z * BOHR_TO_ANGSTROM,
    });
  }

  return { molId: mol.id, charges: Array.from(charges), dipole, atomSymbols, atomCoords, timeMs };
}

async function runSingle(mol: MolDef): Promise<void> {
  running = true;
  runAllMode = false;
  render();

  try {
    const result = await computeMol(mol);
    results.set(mol.id, result);
  } catch (err) {
    console.error('Charge computation error:', err);
  }

  running = false;
  render();
}

async function runAll(): Promise<void> {
  running = true;
  runAllMode = true;
  stopRequested = false;
  runAllProgress = 0;
  render();

  for (let i = 0; i < MOLECULES.length; i++) {
    if (stopRequested) break;
    const mol = MOLECULES[i];

    runAllProgress = i;
    updateProgress(i, MOLECULES.length, t(mol.labelKey));

    try {
      const result = await computeMol(mol);
      results.set(mol.id, result);
    } catch (err) {
      console.error(`Error computing ${mol.id}:`, err);
    }

    await new Promise<void>(r => setTimeout(r, 0));
  }

  running = false;
  runAllMode = false;
  render();
}

// ── Progress ─────────────────────────────────────────────────────────

function updateProgress(current: number, total: number, label: string): void {
  const area = root.querySelector('#progress-area');
  if (!area) return;
  const pct = (current / total * 100).toFixed(0);
  area.innerHTML = `
    <div class="opt-progress"><div class="opt-progress-bar" style="width:${pct}%"></div></div>
    <div class="opt-progress-text">${label} (${current + 1}/${total})</div>`;
}

// ── Summary table ────────────────────────────────────────────────────

function renderSummary(result: ChargeResult): void {
  const el = root.querySelector('#summary-area');
  if (!el) return;

  let html = `<div class="opt-summary">
    <h3>${t('chg.done')}</h3>
    <table>
      <tr>
        <th>${t('chg.colAtom')}</th>
        <th>${t('chg.colCharge')}</th>
      </tr>`;

  for (let i = 0; i < result.atomSymbols.length; i++) {
    const q = result.charges[i];
    const color = chargeColor(q);
    html += `<tr>
      <td>${result.atomSymbols[i]}${i + 1}</td>
      <td style="font-family:monospace;color:${color};font-weight:600">${q >= 0 ? '+' : ''}${q.toFixed(4)}</td>
    </tr>`;
  }

  html += `</table>
    <div style="margin-top:8px;font-size:0.75rem;">
      <span style="color:var(--color-text-dim)">${t('chg.dipole')}:</span>
      <strong>${result.dipole.debye.toFixed(3)} D</strong>
    </div>
    <div style="font-size:0.68rem;color:var(--color-text-dim);margin-top:2px">
      ${t('chg.time')}: ${result.timeMs < 1000 ? result.timeMs.toFixed(0) + 'ms' : (result.timeMs / 1000).toFixed(1) + 's'}
    </div>
  </div>`;
  el.innerHTML = html;
}

// ── Charge color mapping ─────────────────────────────────────────────

/** Map charge to color: blue(+) → white(0) → red(-) */
function chargeColor(q: number): string {
  const clamp = Math.max(-0.6, Math.min(0.6, q));
  const t = clamp / 0.6; // -1 to +1

  if (t > 0) {
    // Positive → blue
    const r = Math.round(255 * (1 - t));
    const g = Math.round(255 * (1 - t));
    return `rgb(${r},${g},255)`;
  } else {
    // Negative → red
    const s = -t;
    const g = Math.round(255 * (1 - s));
    const b = Math.round(255 * (1 - s));
    return `rgb(255,${g},${b})`;
  }
}

/** Render molecule with Three.js charge viewer. */
function renderMolVis(result: ChargeResult): void {
  const el = root.querySelector('#mol-vis') as HTMLElement | null;
  if (!el) return;
  renderChargeViewer(el, result);
}

// ── Comparison bar chart ─────────────────────────────────────────────

function renderCompareChart(): void {
  const container = root.querySelector('#compare-chart');
  if (!container || results.size < 2) return;

  const tc = getThemeColors();

  // Collect heavy-atom charges (non-H, first atom of each molecule)
  interface ChargePt { label: string; charge: number; }
  const pts: ChargePt[] = [];

  for (const mol of MOLECULES) {
    const r = results.get(mol.id);
    if (!r) continue;
    // Find the heaviest atom (most interesting charge)
    let heavyIdx = 0;
    for (let i = 1; i < r.atomSymbols.length; i++) {
      if (r.atomSymbols[i] !== 'H' && r.atomSymbols[heavyIdx] === 'H') heavyIdx = i;
    }
    const sym = r.atomSymbols[heavyIdx];
    pts.push({ label: `${t(mol.labelKey)} (${sym})`, charge: r.charges[heavyIdx] });
  }

  if (pts.length < 2) return;

  const width = 560, height = 200;
  const ml = 82, mr = 20, mt = 28, mb = 60;
  const pw = width - ml - mr, ph = height - mt - mb;

  const maxQ = Math.max(0.1, ...pts.map(p => Math.abs(p.charge))) * 1.2;
  const barW = pw / pts.length;

  const toY = (q: number) => mt + ph / 2 - (q / maxQ) * (ph / 2);

  let svg = `<svg width="${width}" height="${height}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" style="max-width:100%;margin-top:12px;">`;

  svg += `<rect x="${ml}" y="${mt}" width="${pw}" height="${ph}" fill="${tc.surface}" rx="2"/>`;

  // Zero line
  const y0 = toY(0);
  svg += `<line x1="${ml}" y1="${y0}" x2="${ml + pw}" y2="${y0}" stroke="${tc.axis}" stroke-width="1"/>`;

  // Y grid
  const gridSteps = [maxQ * 0.5, -maxQ * 0.5, maxQ, -maxQ];
  for (const q of gridSteps) {
    const yy = toY(q);
    if (yy >= mt && yy <= mt + ph) {
      svg += `<line x1="${ml}" y1="${yy}" x2="${ml + pw}" y2="${yy}" stroke="${tc.grid}" stroke-width="0.5"/>`;
      svg += `<text x="${ml - 6}" y="${yy + 3}" text-anchor="end" font-size="8" fill="${tc.dim}">${q >= 0 ? '+' : ''}${q.toFixed(2)}</text>`;
    }
  }
  svg += `<text x="${ml - 6}" y="${y0 + 3}" text-anchor="end" font-size="8" fill="${tc.dim}">0</text>`;

  // Axes
  svg += `<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt + ph}" stroke="${tc.axis}" stroke-width="1"/>`;

  // Bars
  for (let i = 0; i < pts.length; i++) {
    const p = pts[i];
    const cx = ml + barW * i + barW / 2;
    const bx = ml + barW * i + barW * 0.2;
    const bw = barW * 0.6;
    const y = toY(p.charge);
    const fill = chargeColor(p.charge);

    if (p.charge >= 0) {
      svg += `<rect x="${bx}" y="${y}" width="${bw}" height="${y0 - y}" fill="${fill}" opacity="0.85" rx="2"/>`;
    } else {
      svg += `<rect x="${bx}" y="${y0}" width="${bw}" height="${y - y0}" fill="${fill}" opacity="0.85" rx="2"/>`;
    }

    // Value
    const vy = p.charge >= 0 ? y - 5 : y + 12;
    svg += `<text x="${cx}" y="${vy}" text-anchor="middle" font-size="8" fill="${chargeColor(p.charge)}" font-weight="700">${p.charge >= 0 ? '+' : ''}${p.charge.toFixed(3)}</text>`;

    // Label
    svg += `<text x="${cx}" y="${mt + ph + 12}" text-anchor="middle" font-size="7.5" fill="${tc.dim}" transform="rotate(-30,${cx},${mt + ph + 12})">${p.label}</text>`;
  }

  // Title
  svg += `<text x="${ml + pw / 2}" y="16" text-anchor="middle" font-size="11" font-weight="600" fill="${tc.titleSvg}">${t('chg.compareTitle')}</text>`;
  svg += `<text x="14" y="${mt + ph / 2}" text-anchor="middle" font-size="9" fill="${tc.dim}" transform="rotate(-90,14,${mt + ph / 2})">${t('chg.yCharge')}</text>`;

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
    .opt-page { max-width: 1000px; margin: 0 auto; padding: 16px 20px; }

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
    .opt-category-row { display: flex; flex-wrap: wrap; gap: 5px; margin-bottom: 6px; }
    .opt-scenario-card {
      padding: 5px 10px; border: 1px solid var(--color-border); border-radius: 6px;
      cursor: pointer; transition: all 0.15s; flex: 0 0 auto;
    }
    .opt-scenario-card:hover { background: var(--color-surface-alt); }
    .opt-scenario-card.selected {
      border-color: var(--color-accent); background: var(--color-surface-alt);
      box-shadow: 0 0 0 1px var(--color-accent);
    }
    .opt-scenario-card.computed { border-left: 3px solid var(--color-converged); }
    .opt-scenario-card strong { display: block; font-size: 0.78rem; white-space: nowrap; }
    .conv-desc { display: block !important; font-size: 0.68rem; color: var(--color-text-dim); }

    .chg-cat-label {
      font-size: 0.68rem; font-weight: 600; color: var(--color-text-dim);
      text-transform: uppercase; letter-spacing: 0.04em;
      margin-top: 4px;
    }

    .ladder-basis-info { margin-top: 10px; }

    .opt-run-btn {
      width: 100%; margin-top: 14px; padding: 10px; border: none; border-radius: 8px;
      font-size: 0.85rem; font-weight: 600; cursor: pointer;
      background: var(--color-accent); color: var(--color-accent-on);
      transition: background 0.15s;
    }
    .opt-run-btn:hover:not([disabled]) { background: var(--color-accent-hover); }
    .opt-run-btn[disabled] { opacity: 0.6; cursor: not-allowed; }
    .opt-run-all-btn {
      width: 100%; margin-top: 6px; padding: 8px; border: 1px solid var(--color-accent); border-radius: 8px;
      font-size: 0.8rem; font-weight: 600; cursor: pointer;
      background: none; color: var(--color-accent);
    }
    .opt-run-all-btn:hover:not([disabled]) { background: var(--color-surface-alt); }
    .opt-run-all-btn[disabled] { opacity: 0.6; cursor: not-allowed; }
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

    #mol-vis { width: 100%; text-align: center; }
    #compare-chart { width: 100%; text-align: center; }
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
