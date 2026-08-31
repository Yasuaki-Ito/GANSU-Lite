/**
 * GANSU Lite — Potential Energy Scan demo page.
 * Scans a structural parameter (bond distance or angle) and plots energy vs parameter.
 * Supports RHF/UHF comparison with density matrix reuse for faster convergence.
 */

import { parseXYZ } from './core/parseXYZ';
import { BasisSet } from './core/basisSet';
import { Molecular } from './core/molecular';
import { theorySelectHTML, buildHFOrDFT, theoryMethodLabel, HEAVY_ITERATIVE_EXCLUDE, type TheoryChoice } from './ui/theoryControls';
import { RHF } from './core/rhf';
import { UHF } from './core/uhf';
import { Matrix } from './linalg/matrix';
import { initTheme, toggleTheme, isDark, getThemeColors } from './ui/theme';
import { t, initLang, toggleLang } from './ui/i18n';
import { renderHeader } from './ui/nav';
import { showDualPeriodicTable, COVALENT_RADII } from './ui/periodicTable';

// ── Scenario definitions ─────────────────────────────────────────────

type ParamType = 'bond' | 'angle' | 'height' | 'dihedral';
type ScenarioCategory = 'dissociation' | 'geometry';

interface Scenario {
  id: string;
  category: ScenarioCategory;
  labelKey: string;
  descKey: string;
  paramType: ParamType;
  defaultMin: number;
  defaultMax: number;
  defaultSteps: number;
  defaultCharge: number;
  defaultMult: number;
  basisOptions?: string[];
  paramAtomPair?: [number, number];
  generateXYZ: (param: number) => string;
}

const SCENARIOS: Scenario[] = [
  // ── Bond Dissociation ──
  {
    id: 'h2', category: 'dissociation',
    labelKey: 'opt.scenH2', descKey: 'opt.descH2',
    paramType: 'bond', defaultMin: 0.4, defaultMax: 3.0, defaultSteps: 20,
    defaultCharge: 0, defaultMult: 1,
    basisOptions: ['3-21G', '6-31G'],
    generateXYZ: (r) => `2\nH2 R=${r.toFixed(3)}\nH  0.0  0.0  0.0\nH  0.0  0.0  ${r.toFixed(6)}`,
  },
  {
    id: 'hf', category: 'dissociation',
    labelKey: 'opt.scenHF', descKey: 'opt.descHF',
    paramType: 'bond', defaultMin: 0.5, defaultMax: 3.0, defaultSteps: 20,
    defaultCharge: 0, defaultMult: 1,
    generateXYZ: (r) => `2\nHF R=${r.toFixed(3)}\nH  0.0  0.0  0.0\nF  0.0  0.0  ${r.toFixed(6)}`,
  },
  {
    id: 'lih', category: 'dissociation',
    labelKey: 'opt.scenLiH', descKey: 'opt.descLiH',
    paramType: 'bond', defaultMin: 0.8, defaultMax: 4.0, defaultSteps: 20,
    defaultCharge: 0, defaultMult: 1,
    generateXYZ: (r) => `2\nLiH R=${r.toFixed(3)}\nLi  0.0  0.0  0.0\nH   0.0  0.0  ${r.toFixed(6)}`,
  },
  {
    id: 'n2', category: 'dissociation',
    labelKey: 'opt.scenN2', descKey: 'opt.descN2',
    paramType: 'bond', defaultMin: 0.8, defaultMax: 3.0, defaultSteps: 20,
    defaultCharge: 0, defaultMult: 1,
    generateXYZ: (r) => `2\nN2 R=${r.toFixed(3)}\nN  0.0  0.0  0.0\nN  0.0  0.0  ${r.toFixed(6)}`,
  },
  {
    id: 'f2', category: 'dissociation',
    labelKey: 'opt.scenF2', descKey: 'opt.descF2',
    paramType: 'bond', defaultMin: 0.8, defaultMax: 3.5, defaultSteps: 20,
    defaultCharge: 0, defaultMult: 1,
    generateXYZ: (r) => `2\nF2 R=${r.toFixed(3)}\nF  0.0  0.0  0.0\nF  0.0  0.0  ${r.toFixed(6)}`,
  },
  {
    id: 'heh+', category: 'dissociation',
    labelKey: 'opt.scenHeH', descKey: 'opt.descHeH',
    paramType: 'bond', defaultMin: 0.5, defaultMax: 3.0, defaultSteps: 20,
    defaultCharge: 1, defaultMult: 1,
    generateXYZ: (r) => `2\nHeH+ R=${r.toFixed(3)}\nHe  0.0  0.0  0.0\nH   0.0  0.0  ${r.toFixed(6)}`,
  },
  {
    id: 'li2', category: 'dissociation',
    labelKey: 'opt.scenLi2', descKey: 'opt.descLi2',
    paramType: 'bond', defaultMin: 1.5, defaultMax: 5.0, defaultSteps: 20,
    defaultCharge: 0, defaultMult: 1,
    generateXYZ: (r) => `2\nLi2 R=${r.toFixed(3)}\nLi  0.0  0.0  0.0\nLi  0.0  0.0  ${r.toFixed(6)}`,
  },
  {
    id: 'he2', category: 'dissociation',
    labelKey: 'opt.scenHe2', descKey: 'opt.descHe2',
    paramType: 'bond', defaultMin: 1.0, defaultMax: 5.0, defaultSteps: 20,
    defaultCharge: 0, defaultMult: 1,
    generateXYZ: (r) => `2\nHe2 R=${r.toFixed(3)}\nHe  0.0  0.0  0.0\nHe  0.0  0.0  ${r.toFixed(6)}`,
  },
  {
    id: 'c2h2', category: 'dissociation',
    labelKey: 'opt.scenC2H2', descKey: 'opt.descC2H2',
    paramType: 'bond', defaultMin: 0.9, defaultMax: 2.5, defaultSteps: 20,
    defaultCharge: 0, defaultMult: 1,
    basisOptions: ['STO-3G'],
    paramAtomPair: [1, 2],
    generateXYZ: (r) => {
      const rCH = 1.06;
      return `4\nC2H2 CC=${r.toFixed(3)}\nH  0.0  0.0  ${(-rCH).toFixed(6)}\nC  0.0  0.0  0.0\nC  0.0  0.0  ${r.toFixed(6)}\nH  0.0  0.0  ${(r + rCH).toFixed(6)}`;
    },
  },
  {
    id: 'c2h4', category: 'dissociation',
    labelKey: 'opt.scenC2H4', descKey: 'opt.descC2H4',
    paramType: 'bond', defaultMin: 1.0, defaultMax: 2.5, defaultSteps: 20,
    defaultCharge: 0, defaultMult: 1,
    basisOptions: ['STO-3G'],
    generateXYZ: (r) => {
      const rCH = 1.08, ang = 121.7 * Math.PI / 180;
      const hx = rCH * Math.sin(ang - Math.PI / 2) * -1;
      const hz = rCH * Math.cos(ang - Math.PI / 2);
      // hx ≈ 0.935, hz ≈ 0.54  (H-C-C angle 121.7°)
      const Hx = Math.abs(hx);
      return [
        '6', `C2H4 CC=${r.toFixed(3)}`,
        `C   0.000000  0.000000  0.000000`,
        `C   0.000000  0.000000  ${r.toFixed(6)}`,
        `H   ${Hx.toFixed(6)}  0.000000  ${(-hz).toFixed(6)}`,
        `H   ${(-Hx).toFixed(6)}  0.000000  ${(-hz).toFixed(6)}`,
        `H   ${Hx.toFixed(6)}  0.000000  ${(r + hz).toFixed(6)}`,
        `H   ${(-Hx).toFixed(6)}  0.000000  ${(r + hz).toFixed(6)}`,
      ].join('\n');
    },
  },
  {
    id: 'h2o_bond', category: 'dissociation',
    labelKey: 'opt.scenH2Obond', descKey: 'opt.descH2Obond',
    paramType: 'bond', defaultMin: 0.5, defaultMax: 3.0, defaultSteps: 20,
    defaultCharge: 0, defaultMult: 1,
    basisOptions: ['3-21G', '6-31G'],
    paramAtomPair: [0, 2],
    generateXYZ: (r) => {
      // O at origin, H₁ fixed at 0.96 Å (half-angle 52°), H₂ stretched
      const s = Math.sin(52 * Math.PI / 180), c = Math.cos(52 * Math.PI / 180);
      return [
        '3', `H2O R=${r.toFixed(3)}`,
        `O   0.000000  0.000000  0.000000`,
        `H   ${(0.96 * s).toFixed(6)}  0.000000  ${(0.96 * c).toFixed(6)}`,
        `H   ${(-r * s).toFixed(6)}  0.000000  ${(r * c).toFixed(6)}`,
      ].join('\n');
    },
  },
  // ── Geometry Scan ──
  {
    id: 'h2o', category: 'geometry',
    labelKey: 'opt.scenH2O', descKey: 'opt.descH2O',
    paramType: 'angle', defaultMin: 80, defaultMax: 180, defaultSteps: 20,
    defaultCharge: 0, defaultMult: 1,
    generateXYZ: (angleDeg) => {
      const R = 0.96;
      const rad = angleDeg * Math.PI / 180;
      const halfAngle = rad / 2;
      const hx = R * Math.sin(halfAngle);
      const hz = R * Math.cos(halfAngle);
      return `3\nH2O angle=${angleDeg.toFixed(1)}\nO  0.0  0.0  0.0\nH  ${hx.toFixed(6)}  0.0  ${hz.toFixed(6)}\nH  ${(-hx).toFixed(6)}  0.0  ${hz.toFixed(6)}`;
    },
  },
  {
    id: 'beh2', category: 'geometry',
    labelKey: 'opt.scenBeH2', descKey: 'opt.descBeH2',
    paramType: 'angle', defaultMin: 90, defaultMax: 270, defaultSteps: 20,
    defaultCharge: 0, defaultMult: 1,
    generateXYZ: (angleDeg) => {
      const R = 1.33;
      const rad = angleDeg * Math.PI / 180;
      const halfAngle = rad / 2;
      const hx = R * Math.sin(halfAngle);
      const hz = R * Math.cos(halfAngle);
      return `3\nBeH2 angle=${angleDeg.toFixed(1)}\nBe  0.0  0.0  0.0\nH   ${hx.toFixed(6)}  0.0  ${hz.toFixed(6)}\nH   ${(-hx).toFixed(6)}  0.0  ${hz.toFixed(6)}`;
    },
  },
  {
    id: 'nh3', category: 'geometry',
    labelKey: 'opt.scenNH3', descKey: 'opt.descNH3',
    paramType: 'height', defaultMin: 0, defaultMax: 0.5, defaultSteps: 20,
    defaultCharge: 0, defaultMult: 1,
    generateXYZ: (h) => {
      // h = signed height of N above H₃ plane (Å). 0 = planar, ±h = pyramidal.
      const RNH = 1.012;
      const r = Math.sqrt(Math.max(0, RNH * RNH - h * h)); // H radius in xy-plane
      const s3 = Math.sqrt(3) / 2;
      return [
        '4', `NH3 h=${h.toFixed(3)}`,
        `N   0.000000  0.000000  ${h.toFixed(6)}`,
        `H   ${r.toFixed(6)}  0.000000  0.000000`,
        `H   ${(-r / 2).toFixed(6)}  ${(r * s3).toFixed(6)}  0.000000`,
        `H   ${(-r / 2).toFixed(6)}  ${(-r * s3).toFixed(6)}  0.000000`,
      ].join('\n');
    },
  },
  {
    id: 'ch2', category: 'geometry',
    labelKey: 'opt.scenCH2', descKey: 'opt.descCH2',
    paramType: 'angle', defaultMin: 90, defaultMax: 180, defaultSteps: 20,
    defaultCharge: 0, defaultMult: 1,
    generateXYZ: (angleDeg) => {
      const R = 1.08;
      const rad = angleDeg * Math.PI / 180;
      const half = rad / 2;
      const hx = R * Math.sin(half);
      const hz = R * Math.cos(half);
      return `3\nCH2 angle=${angleDeg.toFixed(1)}\nC  0.0  0.0  0.0\nH  ${hx.toFixed(6)}  0.0  ${hz.toFixed(6)}\nH  ${(-hx).toFixed(6)}  0.0  ${hz.toFixed(6)}`;
    },
  },
  {
    id: 'fhf', category: 'geometry',
    labelKey: 'opt.scenFHF', descKey: 'opt.descFHF',
    paramType: 'bond', defaultMin: 0.8, defaultMax: 1.5, defaultSteps: 20,
    defaultCharge: -1, defaultMult: 1,
    basisOptions: ['STO-3G'],
    generateXYZ: (r) => {
      // F-F distance fixed at 2.3 Å, scan F₁-H distance
      const dFF = 2.3;
      return [
        '3', `FHF- FH=${r.toFixed(3)}`,
        `F   0.000000  0.000000  0.000000`,
        `H   0.000000  0.000000  ${r.toFixed(6)}`,
        `F   0.000000  0.000000  ${dFF.toFixed(6)}`,
      ].join('\n');
    },
  },
  {
    id: 'h3plus', category: 'geometry',
    labelKey: 'opt.scenH3plus', descKey: 'opt.descH3plus',
    paramType: 'bond', defaultMin: 0.5, defaultMax: 2.5, defaultSteps: 20,
    defaultCharge: 1, defaultMult: 1,
    generateXYZ: (r) => {
      // Equilateral triangle, all H-H = r
      const s3 = Math.sqrt(3) / 2;
      return [
        '3', `H3+ R=${r.toFixed(3)}`,
        `H   0.000000  0.000000  0.000000`,
        `H   ${r.toFixed(6)}  0.000000  0.000000`,
        `H   ${(r / 2).toFixed(6)}  ${(r * s3).toFixed(6)}  0.000000`,
      ].join('\n');
    },
  },
  {
    id: 'c2h6', category: 'geometry',
    labelKey: 'opt.scenC2H6', descKey: 'opt.descC2H6',
    paramType: 'dihedral', defaultMin: 0, defaultMax: 120, defaultSteps: 12,
    defaultCharge: 0, defaultMult: 1,
    basisOptions: ['STO-3G'],
    generateXYZ: (phi) => {
      // C-C along z, scan H-C-C-H dihedral φ (0°=eclipsed, 60°=staggered)
      const rCC = 1.54, rCH = 1.09;
      const hcc = 109.47 * Math.PI / 180; // H-C-C angle
      const zH = rCH * Math.cos(Math.PI - hcc);  // H z-offset from its C
      const rp = rCH * Math.sin(Math.PI - hcc);   // H perpendicular distance
      const lines = ['8', `C2H6 dih=${phi.toFixed(1)}`];
      lines.push(`C   0.000000  0.000000  0.000000`);
      lines.push(`C   0.000000  0.000000  ${rCC.toFixed(6)}`);
      const phiRad = phi * Math.PI / 180;
      for (let k = 0; k < 3; k++) {
        const a = k * 2 * Math.PI / 3;
        lines.push(`H   ${(rp * Math.cos(a)).toFixed(6)}  ${(rp * Math.sin(a)).toFixed(6)}  ${(-zH).toFixed(6)}`);
      }
      for (let k = 0; k < 3; k++) {
        const a = k * 2 * Math.PI / 3 + phiRad;
        lines.push(`H   ${(rp * Math.cos(a)).toFixed(6)}  ${(rp * Math.sin(a)).toFixed(6)}  ${(rCC + zH).toFixed(6)}`);
      }
      return lines.join('\n');
    },
  },
];

const BASIS_OPTIONS = ['STO-3G', '3-21G', '6-31G'];

// ── Types ─────────────────────────────────────────────────────────────

interface ScanPoint {
  param: number;
  rhf: number;
  uhf?: number;
}

type SeriesKey = 'rhf' | 'uhf';

interface SeriesInfo {
  key: SeriesKey;
  label: string;
  color: string;
}

// ── State ─────────────────────────────────────────────────────────────

let selectedScenario = SCENARIOS[0];
let selectedBasis = 'STO-3G';
let enableRHF = true;
let enableUHF = false;
let charge = 0;
let multiplicity = 1;
let theoryChoice: TheoryChoice = 'HF';
let paramMin = selectedScenario.defaultMin;
let paramMax = selectedScenario.defaultMax;
let paramSteps = selectedScenario.defaultSteps;
let running = false;
let stopRequested = false;
let scanResults: ScanPoint[] = [];
let scanDone = false;
let scanElapsedMs = 0;
let selectedScanIdx = 0;

// ── Custom diatomic state ────────────────────────────────────────────
let customAtomA: { symbol: string; z: number } | null = null;
let customAtomB: { symbol: string; z: number } | null = null;

function buildCustomScenario(a: { symbol: string; z: number }, b: { symbol: string; z: number }): Scenario {
  const rA = COVALENT_RADII[a.z] ?? 1.0;
  const rB = COVALENT_RADII[b.z] ?? 1.0;
  const rSum = rA + rB;
  const totalElectrons = a.z + b.z;
  const label = a.symbol === b.symbol ? `${a.symbol}\u2082` : `${a.symbol}${b.symbol}`;
  return {
    id: `custom_${a.symbol}_${b.symbol}`,
    category: 'dissociation',
    labelKey: '',
    descKey: '',
    paramType: 'bond',
    defaultMin: Math.max(0.4, Math.round(rSum * 0.5 * 10) / 10),
    defaultMax: Math.min(5.0, Math.round(rSum * 3.0 * 10) / 10),
    defaultSteps: 20,
    defaultCharge: 0,
    defaultMult: totalElectrons % 2 === 0 ? 1 : 2,
    basisOptions: ['STO-3G'],
    generateXYZ: (r: number) =>
      `2\n${label} R=${r.toFixed(3)}\n${a.symbol}  0.0  0.0  0.0\n${b.symbol}  0.0  0.0  ${r.toFixed(6)}`,
  };
}

// ── Basis set cache ──────────────────────────────────────────────────

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

// ── Series colors ────────────────────────────────────────────────────

function getSeriesColor(key: SeriesKey): string {
  if (key === 'rhf') return isDark() ? '#00d4ff' : '#0077cc';
  return isDark() ? '#ff8844' : '#cc4400';
}

function getActiveSeries(): SeriesInfo[] {
  const series: SeriesInfo[] = [];
  if (enableRHF) series.push({ key: 'rhf', label: theoryMethodLabel(theoryChoice, true), color: getSeriesColor('rhf') });
  if (enableUHF) series.push({ key: 'uhf', label: theoryMethodLabel(theoryChoice, false), color: getSeriesColor('uhf') });
  return series;
}

// ── Main render ──────────────────────────────────────────────────────

const root = document.getElementById('app')!;

function render(): void {
  const unit = (selectedScenario.paramType === 'angle' || selectedScenario.paramType === 'dihedral') ? '°' : 'Å';

  root.innerHTML = `
    <div class="opt-page">
      ${renderHeader('optimize')}

      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${t('opt.scenario')}</h2>
          <div class="opt-scenario-grid" id="scenario-grid"></div>

          <h2>${t('opt.basis')}</h2>
          <div class="opt-basis-row" id="basis-row"></div>

          <h2>${t('opt.method')}</h2>
          <div class="opt-method-row">
            <button id="btn-rhf" class="opt-method-btn${enableRHF ? ' active' : ''}" style="--method-color:${getSeriesColor('rhf')}">${theoryMethodLabel(theoryChoice, true)}</button>
            ${selectedScenario.category !== 'geometry' ? `<button id="btn-uhf" class="opt-method-btn${enableUHF ? ' active' : ''}" style="--method-color:${getSeriesColor('uhf')}">${theoryMethodLabel(theoryChoice, false)}</button>` : ''}
          </div>

          <div class="opt-charge-row">
            <label>
              <span>${t('opt.charge')}</span>
              <input id="inp-charge" type="number" min="-3" max="3" value="${charge}" />
            </label>
            <label>
              <span>${t('opt.mult')}</span>
              <input id="inp-mult" type="number" min="1" max="5" value="${multiplicity}" />
            </label>
          </div>

          <h2>${t('opt.param')}</h2>
          <div class="opt-param-form">
            <label>
              <span>${t('opt.min')}</span>
              <input id="param-min" type="number" step="0.1" value="${paramMin}" />
              <span class="opt-unit">${unit}</span>
            </label>
            <label>
              <span>${t('opt.max')}</span>
              <input id="param-max" type="number" step="0.1" value="${paramMax}" />
              <span class="opt-unit">${unit}</span>
            </label>
            <label>
              <span>${t('opt.steps')}</span>
              <input id="param-steps" type="number" min="3" max="50" value="${paramSteps}" />
            </label>
          </div>

          <div class="theory-row" style="margin:8px 0;display:flex;align-items:center;gap:8px;font-size:0.9rem;">
            <span>Theory:</span>${theorySelectHTML('theory-sel', theoryChoice, '', HEAVY_ITERATIVE_EXCLUDE)}
          </div>

          <button id="run-btn" class="opt-run-btn" ${running ? 'disabled' : ''}>
            ${running ? t('opt.running') : t('opt.run')}
          </button>
          ${running ? `<button id="stop-btn" class="opt-stop-btn">${t('opt.stop')}</button>` : ''}

          <div id="progress-area"></div>
          <div id="result-summary"></div>
        </div>

        <div class="opt-panel opt-graph-panel">
          <div id="mol-vis"></div>
          <div id="graph-container">
            ${!scanDone && !running
              ? `<p class="opt-hint">${t('opt.waiting')}</p>`
              : ''}
          </div>
          <div id="scan-slider-area"></div>
        </div>
      </div>
    </div>`;

  injectStyles();

  // Scenario cards grouped by category
  const grid = root.querySelector('#scenario-grid')!;
  const categories: { key: ScenarioCategory; labelKey: string }[] = [
    { key: 'dissociation', labelKey: 'opt.catDissociation' },
    { key: 'geometry', labelKey: 'opt.catGeometry' },
  ];
  for (const cat of categories) {
    const items = SCENARIOS.filter(s => s.category === cat.key);
    if (items.length === 0) continue;
    const header = document.createElement('div');
    header.className = 'opt-category-header';
    header.textContent = t(cat.labelKey);
    grid.appendChild(header);
    const row = document.createElement('div');
    row.className = 'opt-category-row';
    for (const sc of items) {
      const card = document.createElement('div');
      card.className = 'opt-scenario-card' + (sc.id === selectedScenario.id ? ' selected' : '');
      card.innerHTML = `<strong>${t(sc.labelKey)}</strong><span>${t(sc.descKey)}</span>`;
      card.addEventListener('click', () => {
        if (running) return;
        selectedScenario = sc;
        paramMin = sc.defaultMin;
        paramMax = sc.defaultMax;
        paramSteps = sc.defaultSteps;
        charge = sc.defaultCharge;
        multiplicity = sc.defaultMult;
        if (sc.category === 'geometry') { enableUHF = false; enableRHF = true; }
        customAtomA = null;
        customAtomB = null;
        scanResults = [];
        scanDone = false;
        render();
      });
      row.appendChild(card);
    }
    grid.appendChild(row);
  }

  // Custom diatomic card
  {
    const header = document.createElement('div');
    header.className = 'opt-category-header';
    header.textContent = t('opt.catCustom');
    grid.appendChild(header);

    const row = document.createElement('div');
    row.className = 'opt-category-row';

    const isCustom = selectedScenario.id.startsWith('custom_');
    const card = document.createElement('div');
    card.className = 'opt-scenario-card' + (isCustom ? ' selected' : '');
    const label = customAtomA && customAtomB
      ? `${customAtomA.symbol} + ${customAtomB.symbol}`
      : t('opt.customSelect');
    card.innerHTML = `<strong>${label}</strong><span>${t('opt.catCustom')}</span>`;
    card.addEventListener('click', () => {
      if (running) return;
      showDualPeriodicTable((symA, zA, symB, zB) => {
        customAtomA = { symbol: symA, z: zA };
        customAtomB = { symbol: symB, z: zB };
        applyCustomScenario();
      }, customAtomA, customAtomB);
    });
    row.appendChild(card);
    grid.appendChild(row);
  }

  function applyCustomScenario() {
    if (!customAtomA || !customAtomB) return;
    const sc = buildCustomScenario(customAtomA, customAtomB);
    selectedScenario = sc;
    paramMin = sc.defaultMin;
    paramMax = sc.defaultMax;
    paramSteps = sc.defaultSteps;
    charge = sc.defaultCharge;
    multiplicity = sc.defaultMult;
    // Odd electrons → UHF only; even → both
    const totalE = customAtomA.z + customAtomB.z - charge;
    if (totalE % 2 !== 0) {
      enableRHF = false;
      enableUHF = true;
    } else {
      enableRHF = true;
      enableUHF = true;
    }
    scanResults = [];
    scanDone = false;
    render();
  }

  // Basis buttons
  const basisRow = root.querySelector('#basis-row')!;
  const basisList = selectedScenario.basisOptions ?? BASIS_OPTIONS;
  if (!basisList.includes(selectedBasis)) {
    selectedBasis = basisList[0];
  }
  for (const bs of basisList) {
    const btn = document.createElement('button');
    btn.className = 'opt-basis-btn' + (bs === selectedBasis ? ' selected' : '');
    btn.textContent = bs;
    btn.addEventListener('click', () => {
      if (running) return;
      selectedBasis = bs;
      scanResults = [];
      scanDone = false;
      render();
    });
    basisRow.appendChild(btn);
  }

  // Method toggle buttons
  root.querySelector('#btn-rhf')!.addEventListener('click', () => {
    if (running) return;
    if (enableRHF && !enableUHF) return; // must keep at least one
    enableRHF = !enableRHF;
    scanResults = []; scanDone = false; render();
  });
  root.querySelector('#btn-uhf')?.addEventListener('click', () => {
    if (running) return;
    if (enableUHF && !enableRHF) return;
    enableUHF = !enableUHF;
    scanResults = []; scanDone = false; render();
  });

  // Charge / Multiplicity
  root.querySelector<HTMLInputElement>('#inp-charge')!.addEventListener('change', (e) => {
    charge = parseInt((e.target as HTMLInputElement).value, 10);
  });
  root.querySelector<HTMLInputElement>('#inp-mult')!.addEventListener('change', (e) => {
    multiplicity = parseInt((e.target as HTMLInputElement).value, 10);
  });

  // Wire up events
  root.querySelector('#nav-theme')!.addEventListener('click', () => { toggleTheme(); render(); });
  root.querySelector('#nav-lang')!.addEventListener('click', () => { toggleLang(); render(); });
  const theorySel = root.querySelector<HTMLSelectElement>('#theory-sel');
  if (theorySel) theorySel.addEventListener('change', () => {
    theoryChoice = theorySel.value as TheoryChoice;
    render();  // refresh RHF/UHF button labels with new theory prefix
  });
  root.querySelector('#run-btn')!.addEventListener('click', runScan);
  root.querySelector('#stop-btn')?.addEventListener('click', () => { stopRequested = true; });
  root.querySelector<HTMLInputElement>('#param-min')!.addEventListener('change', (e) => {
    paramMin = parseFloat((e.target as HTMLInputElement).value);
  });
  root.querySelector<HTMLInputElement>('#param-max')!.addEventListener('change', (e) => {
    paramMax = parseFloat((e.target as HTMLInputElement).value);
  });
  root.querySelector<HTMLInputElement>('#param-steps')!.addEventListener('change', (e) => {
    const el = e.target as HTMLInputElement;
    let v = parseInt(el.value, 10);
    if (isNaN(v) || v < 3) v = 3;
    if (v > 50) v = 50;
    el.value = String(v);
    paramSteps = v;
  });

  if (scanResults.length > 0) {
    renderGraph(scanResults, scanDone);
  }
  if (scanDone) {
    renderResultSummary();
    setupScanSlider();
    updateScanSliderVis();
  } else {
    // Show molecule at midpoint parameter
    updateMolVis([{ param: (paramMin + paramMax) / 2 }]);
  }
}

// ── Broken-symmetry initial guess ────────────────────────────────────

/** Mix HOMO and LUMO from RHF coefficients to create alpha≠beta densities */
function buildBrokenSymmetryGuess(uhf: UHF, C: Matrix, n: number, nocc: number): void {
  if (nocc <= 0 || nocc >= n) return;

  const theta = Math.PI / 4;
  const cosT = Math.cos(theta);
  const sinT = Math.sin(theta);
  const homo = nocc - 1;
  const lumo = nocc;

  // C_alpha: HOMO → cos θ · HOMO + sin θ · LUMO
  const Ca = C.clone();
  for (let r = 0; r < n; r++) {
    Ca.set(r, homo, cosT * C.get(r, homo) + sinT * C.get(r, lumo));
  }
  // C_beta: HOMO → cos θ · HOMO − sin θ · LUMO
  const Cb = C.clone();
  for (let r = 0; r < n; r++) {
    Cb.set(r, homo, cosT * C.get(r, homo) - sinT * C.get(r, lumo));
  }

  // P_alpha = Ca_occ · Ca_occ^T,  P_beta = Cb_occ · Cb_occ^T
  const Pa = new Matrix(n, n);
  const Pb = new Matrix(n, n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let sa = 0, sb = 0;
      for (let a = 0; a < nocc; a++) {
        sa += Ca.get(i, a) * Ca.get(j, a);
        sb += Cb.get(i, a) * Cb.get(j, a);
      }
      Pa.set(i, j, sa);
      Pb.set(i, j, sb);
    }
  }
  uhf.setInitialDensityGuessAlphaBeta(Pa, Pb);
}

// ── Run scan ─────────────────────────────────────────────────────────

async function runScan(): Promise<void> {
  if (running) return;
  // Guard against expensive combinations: DFT + UHF + many points = 5+ minutes.
  // Warn before launching so the user can opt out.
  if (theoryChoice !== 'HF' && enableUHF && enableRHF && paramSteps >= 15) {
    const est = Math.round(paramSteps * 2 * 8 / 60);  // ~8 s/SCF × 2 passes × steps
    if (!confirm(`This scan combines DFT × RHF+UHF × ${paramSteps} points → estimated ~${est} min.\nProceed?`)) return;
  }
  running = true;
  stopRequested = false;
  scanResults = [];
  scanDone = false;
  render();

  const steps = Math.max(3, Math.min(50, paramSteps));
  const values: number[] = [];
  for (let i = 0; i <= steps; i++) {
    values.push(paramMin + (paramMax - paramMin) * i / steps);
  }

  const betaToAlpha = Math.floor((multiplicity - 1) / 2);
  const startTime = performance.now();

  const totalPoints = (enableRHF ? values.length : 0) + (enableUHF ? values.length : 0);
  let completedPoints = 0;

  try {
    const basis = await loadBasis(selectedBasis);
    let prevRHFDensity: Matrix | null = null;

    // ── Pass 1: RHF (forward: small R → large R) ──
    if (enableRHF) {
      for (let i = 0; i < values.length; i++) {
        if (stopRequested) break;

        const param = values[i];
        const atoms = parseXYZ(selectedScenario.generateXYZ(param));
        const mol = new Molecular(atoms, basis, charge, betaToAlpha);
        const rhf = await buildHFOrDFT(mol, basis, theoryChoice, 'RHF') as RHF;
        if (prevRHFDensity) rhf.setInitialDensityGuess(prevRHFDensity);
        const rhfEnergy = await rhf.solve({ eriBackend: 'js' });
        prevRHFDensity = rhf.density.clone();

        scanResults.push({ param, rhf: rhfEnergy });
        completedPoints++;
        updateProgress(completedPoints, totalPoints);
        renderGraph(scanResults, false);
        {
          const entries: MolVisEntry[] = [{ param, label: theoryMethodLabel(theoryChoice, true), color: getSeriesColor('rhf') }];
          if (enableUHF) entries.push({ param: (paramMin + paramMax) / 2, label: theoryMethodLabel(theoryChoice, false), color: getSeriesColor('uhf') });
          updateMolVis(entries);
        }
        await new Promise<void>(r => setTimeout(r, 0));
      }
    } else {
      // Prepare empty results for UHF-only mode
      for (const v of values) scanResults.push({ param: v, rhf: NaN });
    }

    // ── Pass 2: UHF (reverse: large R → small R) ──
    if (enableUHF && !stopRequested) {
      let prevUHFAlpha: Matrix | null = null;
      let prevUHFBeta: Matrix | null = null;

      // Broken-symmetry initial guess from RHF at largest R (even electrons only)
      let rhfCoeffsAtMax: Matrix | null = null;
      const totalE = selectedScenario.defaultCharge !== undefined
        ? (values.length > 0 ? parseXYZ(selectedScenario.generateXYZ(values[0])).reduce((s, a) => s + a.atomicNumber, 0) : 0) - charge
        : 0;
      const isEvenElectrons = totalE % 2 === 0;

      if (isEvenElectrons) {
        const lastAtoms = parseXYZ(selectedScenario.generateXYZ(values[values.length - 1]));
        const lastMol = new Molecular(lastAtoms, basis, charge, betaToAlpha);
        const lastRHF = await buildHFOrDFT(lastMol, basis, theoryChoice, 'RHF') as RHF;
        if (prevRHFDensity) lastRHF.setInitialDensityGuess(prevRHFDensity);
        await lastRHF.solve({ eriBackend: 'js' });
        rhfCoeffsAtMax = lastRHF.coefficients;
      }

      for (let i = values.length - 1; i >= 0; i--) {
        if (stopRequested) break;

        const param = values[i];
        const atoms = parseXYZ(selectedScenario.generateXYZ(param));

        try {
          const mol = new Molecular(atoms, basis, charge, betaToAlpha);
          const uhf = await buildHFOrDFT(mol, basis, theoryChoice, 'UHF') as UHF;

          if (prevUHFAlpha && prevUHFBeta) {
            uhf.setInitialDensityGuessAlphaBeta(prevUHFAlpha, prevUHFBeta);
          } else if (rhfCoeffsAtMax) {
            buildBrokenSymmetryGuess(uhf, rhfCoeffsAtMax, mol.numBasis, mol.numAlphaSpins);
          }

          let lastDeltaE = Infinity;
          const uhfEnergy = await uhf.solve({
            eriBackend: 'js',
            onIteration: (_iter, _e, dE) => { lastDeltaE = Math.abs(dE); },
          });

          if (lastDeltaE < 1e-6) {
            prevUHFAlpha = uhf.densityAlphaMatrix.clone();
            prevUHFBeta = uhf.densityBetaMatrix.clone();
            scanResults[i].uhf = uhfEnergy;
          } else {
            console.warn(`UHF not converged at param=${param}, deltaE=${lastDeltaE}`);
          }
        } catch (err) {
          console.warn('UHF error at param', param, err);
        }

        completedPoints++;
        updateProgress(completedPoints, totalPoints);
        renderGraph(scanResults, false);
        {
          const entries: MolVisEntry[] = [];
          if (enableRHF) {
            const rhfValid = scanResults.filter(p => isFinite(p.rhf));
            const rhfBest = rhfValid.length > 0 ? rhfValid.reduce((a, b) => a.rhf < b.rhf ? a : b) : null;
            entries.push({ param: rhfBest?.param ?? (paramMin + paramMax) / 2, label: theoryMethodLabel(theoryChoice, true), color: getSeriesColor('rhf') });
          }
          entries.push({ param, label: theoryMethodLabel(theoryChoice, false), color: getSeriesColor('uhf') });
          updateMolVis(entries);
        }
        await new Promise<void>(r => setTimeout(r, 0));
      }
    }
  } catch (err) {
    console.error('Scan error:', err);
  }

  scanElapsedMs = performance.now() - startTime;
  running = false;
  scanDone = true;

  // Default slider to minimum energy point
  const primaryKey = enableRHF ? 'rhf' as SeriesKey : 'uhf' as SeriesKey;
  let bestIdx = 0;
  for (let i = 1; i < scanResults.length; i++) {
    const v = scanResults[i][primaryKey];
    const bv = scanResults[bestIdx][primaryKey];
    if (v != null && isFinite(v) && (bv == null || !isFinite(bv) || v < bv)) bestIdx = i;
  }
  selectedScanIdx = bestIdx;

  render();
}

// ── Progress ─────────────────────────────────────────────────────────

function updateProgress(n: number, total: number): void {
  const area = root.querySelector('#progress-area');
  if (!area) return;
  const pct = (n / total * 100).toFixed(0);
  const text = t('opt.progress').replace('{n}', String(n)).replace('{total}', String(total));
  area.innerHTML = `
    <div class="opt-progress">
      <div class="opt-progress-bar" style="width:${pct}%"></div>
    </div>
    <p class="opt-progress-text">${text}</p>`;
}

// ── Result summary ───────────────────────────────────────────────────

function renderResultSummary(): void {
  const area = root.querySelector('#result-summary');
  if (!area || scanResults.length === 0) return;

  const unit = (selectedScenario.paramType === 'angle' || selectedScenario.paramType === 'dihedral') ? '°' : 'Å';
  const elapsed = (scanElapsedMs / 1000).toFixed(1);
  const series = getActiveSeries();

  let rows = '';
  for (const s of series) {
    const vals = scanResults.filter(p => p[s.key] != null && isFinite(p[s.key]!)).map(p => ({ param: p.param, e: p[s.key]! }));
    if (vals.length === 0) continue;
    let best = vals[0];
    for (const v of vals) { if (v.e < best.e) best = v; }
    rows += `<tr>
      <td><span class="opt-dot" style="background:${s.color}"></span>${s.label}</td>
      <td>${best.param.toFixed(4)} ${unit}</td>
      <td><strong>${best.e.toFixed(8)} Eh</strong></td>
    </tr>`;
  }

  area.innerHTML = `
    <div class="opt-summary">
      <h3>${t('opt.done')} (${elapsed} s)</h3>
      <table>
        <tr><th></th><th>${t('opt.resultParam')}</th><th>${t('opt.resultEnergy')}</th></tr>
        ${rows}
      </table>
    </div>`;
}

// ── Scan slider ─────────────────────────────────────────────────────

function setupScanSlider(): void {
  const area = root.querySelector('#scan-slider-area');
  if (!area || scanResults.length === 0) return;
  if (selectedScanIdx >= scanResults.length) selectedScanIdx = 0;

  const isAngle = selectedScenario.paramType === 'angle' || selectedScenario.paramType === 'dihedral';
  const unit = isAngle ? '\u00b0' : ' \u00c5';
  const fmt = isAngle ? scanResults[selectedScanIdx].param.toFixed(1) : scanResults[selectedScanIdx].param.toFixed(3);

  area.innerHTML = `
    <div class="walsh-slider">
      <input type="range" id="scan-slider" min="0" max="${scanResults.length - 1}" value="${selectedScanIdx}" />
      <div class="walsh-slider-label">${fmt}${unit}</div>
    </div>`;

  root.querySelector('#scan-slider')?.addEventListener('input', (e) => {
    selectedScanIdx = parseInt((e.target as HTMLInputElement).value);
    const p = scanResults[selectedScanIdx].param;
    const label = root.querySelector('.walsh-slider-label');
    if (label) label.textContent = `${isAngle ? p.toFixed(1) : p.toFixed(3)}${unit}`;
    renderGraph(scanResults, true);
    updateScanSliderVis();
  });
}

function updateScanSliderVis(): void {
  if (scanResults.length === 0) return;
  if (selectedScanIdx >= scanResults.length) selectedScanIdx = 0;
  const pt = scanResults[selectedScanIdx];
  const series = getActiveSeries();
  const entries: MolVisEntry[] = [];
  for (const s of series) {
    entries.push({ param: pt.param, label: s.label, color: s.color });
  }
  if (entries.length > 0) updateMolVis(entries);
}

// ── Molecule Visualization ───────────────────────────────────────────

const ELEMENT_COLORS: Record<string, string> = {
  H: '#999', He: '#0CC', Li: '#C2C', Be: '#6C0',
  B: '#F90', C: '#555', N: '#35F', O: '#F22', F: '#9E5',
};
const COVALENT_R: Record<string, number> = {
  H: 0.31, He: 0.28, Li: 1.28, Be: 0.96,
  B: 0.84, C: 0.76, N: 0.71, O: 0.66, F: 0.57,
};

type VisAtom = { sym: string; x: number; y: number; z: number };

type ProjAxes = { pU: 'x' | 'y' | 'z'; pV: 'x' | 'y' | 'z'; pW: 'x' | 'y' | 'z' };

/** Determine best 2D projection axes from a set of atoms */
function chooseProjAxes(allAtoms: VisAtom[]): ProjAxes {
  if (selectedScenario.paramType === 'height') {
    return { pU: 'x', pV: 'z', pW: 'y' }; // special: rotated in projectAtoms
  }
  const sp = (['x', 'y', 'z'] as const).map(a => ({
    a, r: Math.max(...allAtoms.map(at => at[a])) - Math.min(...allAtoms.map(at => at[a])),
  })).sort((a, b) => b.r - a.r);
  return { pU: sp[0].a, pV: sp[1].a, pW: sp[2].a };
}

function projectAtoms(atomList: VisAtom[], axes?: ProjAxes): { u: number; v: number }[] {
  const ax = axes ?? chooseProjAxes(atomList);
  if (selectedScenario.paramType === 'height') {
    const ca = Math.cos(Math.PI / 6), sa = Math.sin(Math.PI / 6);
    return atomList.map(a => ({ u: a.x * ca + a.y * sa, v: a.z }));
  }
  return atomList.map(a => ({ u: a[ax.pU], v: a[ax.pV] }));
}

function parseVisAtoms(xyzText: string): VisAtom[] {
  const lines = xyzText.split('\n');
  const n = parseInt(lines[0].trim());
  const atoms: VisAtom[] = [];
  for (let i = 2; i < 2 + n; i++) {
    const p = lines[i].trim().split(/\s+/);
    atoms.push({ sym: p[0], x: +p[1], y: +p[2], z: +p[3] });
  }
  return atoms;
}

/** Compute fixed scale and center from paramMin/paramMax extents */
function getFixedScaleInfo(W: number, H: number, margin: number) {
  // Collect all atoms across param range to determine consistent projection axes
  let allAtoms: VisAtom[] = [];
  for (const rp of [paramMin, paramMax]) {
    allAtoms.push(...parseVisAtoms(selectedScenario.generateXYZ(rp)));
  }
  const axes = chooseProjAxes(allAtoms);

  let allU: number[] = [], allV: number[] = [];
  for (const rp of [paramMin, paramMax]) {
    const pts = projectAtoms(parseVisAtoms(selectedScenario.generateXYZ(rp)), axes);
    allU.push(...pts.map(p => p.u));
    allV.push(...pts.map(p => p.v));
  }
  const uMin = Math.min(...allU), uMax = Math.max(...allU);
  const vMin = Math.min(...allV), vMax = Math.max(...allV);
  const uRange = uMax - uMin || 0.01, vRange = vMax - vMin || 0.01;
  // Extra bottom margin for annotation (arrow + text below atoms)
  const marginTop = margin;
  const marginBottom = margin + 30;
  const scale = Math.min((W - 2 * margin) / uRange, (H - marginTop - marginBottom) / vRange, 120);
  const vMid = (vMin + vMax) / 2;
  // Shift center upward to account for asymmetric vertical margins
  const cySvgOffset = (marginTop - marginBottom) / 2;
  return { scale, uMid: (uMin + uMax) / 2, vMid, cySvgOffset, axes };
}

/** Render a single molecule SVG */
function renderMolSVG(param: number, W: number, H: number, label?: string, labelColor?: string): string {
  const atoms = parseVisAtoms(selectedScenario.generateXYZ(param));
  if (atoms.length === 0) return '';

  const tc = getThemeColors();
  const margin = 35;
  const { scale, uMid, vMid, cySvgOffset, axes } = getFixedScaleInfo(W, H, margin);
  const cxSvg = W / 2, cySvg = H / 2 + cySvgOffset;
  const sx = (u: number) => cxSvg + (u - uMid) * scale;
  const sy = (v: number) => cySvg - (v - vMid) * scale;
  const ar = 18;

  const pts = projectAtoms(atoms, axes);
  let svg = '';

  // Bonds
  for (let i = 0; i < atoms.length; i++) {
    for (let j = i + 1; j < atoms.length; j++) {
      const dx = atoms[i].x - atoms[j].x, dy = atoms[i].y - atoms[j].y, dz = atoms[i].z - atoms[j].z;
      const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
      const thr = 1.4 * ((COVALENT_R[atoms[i].sym] ?? 0.7) + (COVALENT_R[atoms[j].sym] ?? 0.7));
      if (d < thr) {
        svg += `<line x1="${sx(pts[i].u)}" y1="${sy(pts[i].v)}" x2="${sx(pts[j].u)}" y2="${sy(pts[j].v)}" stroke="${tc.grid}" stroke-width="4" stroke-linecap="round"/>`;
      }
    }
  }

  // Atoms (Z-sorted: draw far atoms first so near atoms appear on top)
  const depthOrder = atoms.map((a, i) => i).sort((a, b) => atoms[a][axes.pW] - atoms[b][axes.pW]);
  for (const i of depthOrder) {
    const ax = sx(pts[i].u), ay = sy(pts[i].v);
    const col = ELEMENT_COLORS[atoms[i].sym] ?? '#888';
    svg += `<circle cx="${ax}" cy="${ay}" r="${ar}" fill="${col}" stroke="${tc.axis}" stroke-width="1.2"/>`;
    svg += `<text x="${ax}" y="${ay}" text-anchor="middle" dy="0.38em" font-size="11" font-weight="bold" fill="#fff" stroke="#0003" stroke-width="0.3">${atoms[i].sym}</text>`;
  }

  // Parameter annotation
  const isAngleType = selectedScenario.paramType === 'angle' || selectedScenario.paramType === 'dihedral';
  const paramUnit = isAngleType ? '°' : 'Å';
  const paramFmt = isAngleType ? param.toFixed(1) : param.toFixed(3);
  const paramStr = `${paramFmt} ${paramUnit}`;
  const ac = labelColor ?? tc.accent;

  if (selectedScenario.paramType === 'bond') {
    const [pi, pj] = selectedScenario.paramAtomPair ?? [0, 1];
    const x0 = sx(pts[pi].u), x1 = sx(pts[pj].u);
    const yBase = Math.max(sy(pts[pi].v), sy(pts[pj].v)) + ar + 10;
    svg += `<line x1="${x0}" y1="${yBase}" x2="${x1}" y2="${yBase}" stroke="${ac}" stroke-width="1.2"/>`;
    svg += `<line x1="${x0}" y1="${yBase - 4}" x2="${x0}" y2="${yBase + 4}" stroke="${ac}" stroke-width="1.2"/>`;
    svg += `<line x1="${x1}" y1="${yBase - 4}" x2="${x1}" y2="${yBase + 4}" stroke="${ac}" stroke-width="1.2"/>`;
    const dir = x1 > x0 ? 1 : -1;
    svg += `<polygon points="${x0},${yBase} ${x0 + dir * 6},${yBase - 3} ${x0 + dir * 6},${yBase + 3}" fill="${ac}"/>`;
    svg += `<polygon points="${x1},${yBase} ${x1 - dir * 6},${yBase - 3} ${x1 - dir * 6},${yBase + 3}" fill="${ac}"/>`;
    svg += `<text x="${(x0 + x1) / 2}" y="${yBase + 14}" text-anchor="middle" font-size="11" font-weight="600" fill="${ac}">${paramStr}</text>`;

  } else if (selectedScenario.paramType === 'angle') {
    const vx = sx(pts[0].u), vy = sy(pts[0].v);
    const a1x = sx(pts[1].u), a1y = sy(pts[1].v);
    const a2x = sx(pts[2].u), a2y = sy(pts[2].v);
    const ang1 = Math.atan2(a1y - vy, a1x - vx);
    const ang2 = Math.atan2(a2y - vy, a2x - vx);
    const arcR = 22;
    const sX = vx + arcR * Math.cos(ang1), sY = vy + arcR * Math.sin(ang1);
    const eX = vx + arcR * Math.cos(ang2), eY = vy + arcR * Math.sin(ang2);
    const cross = (sX - vx) * (eY - vy) - (sY - vy) * (eX - vx);
    const sweep = cross > 0 ? 1 : 0;
    svg += `<path d="M ${sX} ${sY} A ${arcR} ${arcR} 0 0 ${sweep} ${eX} ${eY}" fill="none" stroke="${ac}" stroke-width="1.5"/>`;
    let ad = ang2 - ang1; if (ad > Math.PI) ad -= 2 * Math.PI; if (ad < -Math.PI) ad += 2 * Math.PI;
    const midA = ang1 + ad / 2;
    const lR = arcR + 14;
    svg += `<text x="${vx + lR * Math.cos(midA)}" y="${vy + lR * Math.sin(midA)}" text-anchor="middle" dy="0.35em" font-size="11" font-weight="600" fill="${ac}">${paramStr}</text>`;
    const dx01 = atoms[1].x - atoms[0].x, dy01 = atoms[1].y - atoms[0].y, dz01 = atoms[1].z - atoms[0].z;
    const bondLen = Math.sqrt(dx01 * dx01 + dy01 * dy01 + dz01 * dz01);
    const mx = (vx + a1x) / 2, my = (vy + a1y) / 2;
    const bAng = Math.atan2(a1y - vy, a1x - vx);
    const offX = mx + 12 * Math.cos(bAng + Math.PI / 2);
    const offY = my + 12 * Math.sin(bAng + Math.PI / 2);
    svg += `<text x="${offX}" y="${offY}" text-anchor="middle" dy="0.35em" font-size="9" fill="${tc.dim}">${bondLen.toFixed(2)} Å</text>`;

  } else if (selectedScenario.paramType === 'height') {
    const nx = sx(pts[0].u), ny = sy(pts[0].v);
    const py = sy(0);
    const planeLeft = Math.min(...pts.slice(1).map(p => sx(p.u))) - 20;
    const planeRight = Math.max(...pts.slice(1).map(p => sx(p.u))) + 20;
    svg += `<line x1="${planeLeft}" y1="${py}" x2="${planeRight}" y2="${py}" stroke="${tc.dim}" stroke-width="1" stroke-dasharray="5,3"/>`;
    svg += `<text x="${planeRight + 4}" y="${py}" dy="0.35em" font-size="9" fill="${tc.dim}">H\u2083 plane</text>`;
    if (Math.abs(ny - py) > 8) {
      const dx = nx > cxSvg ? -24 : 24;
      svg += `<line x1="${nx + dx}" y1="${ny}" x2="${nx + dx}" y2="${py}" stroke="${ac}" stroke-width="1.5"/>`;
      svg += `<line x1="${nx + dx - 5}" y1="${ny}" x2="${nx + dx + 5}" y2="${ny}" stroke="${ac}" stroke-width="1.5"/>`;
      svg += `<line x1="${nx + dx - 5}" y1="${py}" x2="${nx + dx + 5}" y2="${py}" stroke="${ac}" stroke-width="1.5"/>`;
      const dir = ny < py ? 1 : -1;
      svg += `<polygon points="${nx + dx},${ny} ${nx + dx - 3},${ny + dir * 6} ${nx + dx + 3},${ny + dir * 6}" fill="${ac}"/>`;
      svg += `<polygon points="${nx + dx},${py} ${nx + dx - 3},${py - dir * 6} ${nx + dx + 3},${py - dir * 6}" fill="${ac}"/>`;
      const labelX = nx + dx + (dx > 0 ? 10 : -10);
      const anchor = dx > 0 ? 'start' : 'end';
      svg += `<text x="${labelX}" y="${(ny + py) / 2}" text-anchor="${anchor}" dy="0.35em" font-size="11" font-weight="600" fill="${ac}">h = ${paramStr}</text>`;
    } else {
      svg += `<text x="${nx}" y="${py + ar + 14}" text-anchor="middle" font-size="11" font-weight="600" fill="${ac}">h = ${paramStr}</text>`;
    }

  } else if (selectedScenario.paramType === 'dihedral') {
    // Show dihedral angle label at bottom center
    const yBase = Math.max(...pts.map(p => sy(p.v))) + ar + 12;
    svg += `<text x="${cxSvg}" y="${yBase}" text-anchor="middle" font-size="11" font-weight="600" fill="${ac}">\u03c6 = ${paramStr}</text>`;
  }

  // Label (method name)
  if (label) {
    svg += `<text x="${W / 2}" y="14" text-anchor="middle" font-size="11" font-weight="700" fill="${ac}">${label}</text>`;
  }

  return svg;
}

interface MolVisEntry { param: number; label?: string; color?: string; }

function updateMolVis(entries: MolVisEntry[]): void {
  const el = root.querySelector('#mol-vis');
  if (!el) return;
  if (entries.length === 0) return;

  const totalW = 520;
  const H = 200;

  if (entries.length === 1) {
    const e = entries[0];
    const svg = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${totalW} ${H}" style="width:100%;max-width:${totalW}px;display:block;margin:0 auto;">${renderMolSVG(e.param, totalW, H, e.label, e.color)}</svg>`;
    el.innerHTML = svg;
  } else {
    // Side by side
    const halfW = Math.floor(totalW / entries.length);
    let html = `<div style="display:flex;gap:4px;justify-content:center;">`;
    for (const e of entries) {
      html += `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${halfW} ${H}" style="flex:1;max-width:${halfW}px;">${renderMolSVG(e.param, halfW, H, e.label, e.color)}</svg>`;
    }
    html += '</div>';
    el.innerHTML = html;
  }
}

// ── SVG Graph ────────────────────────────────────────────────────────

function renderGraph(data: ScanPoint[], done: boolean): void {
  const container = root.querySelector('#graph-container');
  if (!container || data.length < 1) return;

  const tc = getThemeColors();
  const xLabelMap: Record<ParamType, string> = {
    bond: t('opt.xBond'), angle: t('opt.xAngle'), height: t('opt.xHeight'), dihedral: t('opt.xDihedral'),
  };
  const xLabel = xLabelMap[selectedScenario.paramType];
  const series = getActiveSeries();

  const width = 520;
  const height = 360;
  const ml = 72;
  const mr = 24;
  const mt = 36;
  const mb = 44;
  const pw = width - ml - mr;
  const ph = height - mt - mb;

  // Compute Y range across all active series
  const allYs: number[] = [];
  for (const s of series) {
    for (const d of data) {
      const v = d[s.key];
      if (v != null && isFinite(v)) allYs.push(v);
    }
  }
  if (allYs.length === 0) allYs.push(0);
  const yMin = Math.min(...allYs);
  const yMax = Math.max(...allYs);
  const yPad = (yMax - yMin) * 0.1 || 0.01;
  const yLow = yMin - yPad;
  const yHigh = yMax + yPad;

  const toX = (x: number) => ml + ((x - paramMin) / (paramMax - paramMin || 1)) * pw;
  const toY = (y: number) => mt + ph - ((y - yLow) / (yHigh - yLow || 1)) * ph;

  let svg = `<svg width="${width}" height="${height}" xmlns="http://www.w3.org/2000/svg" style="display:block;width:100%;height:auto;max-width:${width}px;" viewBox="0 0 ${width} ${height}">`;

  // Background
  svg += `<rect x="${ml}" y="${mt}" width="${pw}" height="${ph}" fill="${tc.surface}" rx="2"/>`;

  // Grid Y
  for (let i = 0; i <= 5; i++) {
    const yVal = yLow + (yHigh - yLow) * i / 5;
    const yy = toY(yVal);
    svg += `<line x1="${ml}" y1="${yy}" x2="${ml + pw}" y2="${yy}" stroke="${tc.grid}" stroke-width="0.5"/>`;
    svg += `<text x="${ml - 6}" y="${yy + 3}" text-anchor="end" font-size="9" font-family="monospace" fill="${tc.dim}">${yVal.toFixed(4)}</text>`;
  }
  // Grid X
  for (let i = 0; i <= 5; i++) {
    const xVal = paramMin + (paramMax - paramMin) * i / 5;
    const xx = toX(xVal);
    svg += `<line x1="${xx}" y1="${mt}" x2="${xx}" y2="${mt + ph}" stroke="${tc.grid}" stroke-width="0.5"/>`;
    svg += `<text x="${xx}" y="${mt + ph + 14}" text-anchor="middle" font-size="9" fill="${tc.dim}">${xVal.toFixed(2)}</text>`;
  }

  // Axes
  svg += `<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt + ph}" stroke="${tc.axis}" stroke-width="1"/>`;
  svg += `<line x1="${ml}" y1="${mt + ph}" x2="${ml + pw}" y2="${mt + ph}" stroke="${tc.axis}" stroke-width="1"/>`;

  // Draw each series
  for (const s of series) {
    const pts = data.filter(d => d[s.key] != null && isFinite(d[s.key]!)).map(d => ({ x: d.param, y: d[s.key]! }));
    if (pts.length < 2) {
      for (const p of pts) {
        svg += `<circle cx="${toX(p.x).toFixed(1)}" cy="${toY(p.y).toFixed(1)}" r="3" fill="${s.color}"/>`;
      }
      continue;
    }

    // Line
    let path = '';
    for (let i = 0; i < pts.length; i++) {
      const px = toX(pts[i].x);
      const py = toY(pts[i].y);
      path += i === 0 ? `M${px.toFixed(1)},${py.toFixed(1)}` : ` L${px.toFixed(1)},${py.toFixed(1)}`;
    }
    svg += `<path d="${path}" fill="none" stroke="${s.color}" stroke-width="2" stroke-linejoin="round"/>`;

    // Points
    for (const p of pts) {
      svg += `<circle cx="${toX(p.x).toFixed(1)}" cy="${toY(p.y).toFixed(1)}" r="3" fill="${s.color}" stroke="${tc.surface}" stroke-width="0.8"/>`;
    }

    // Min marker
    if (pts.length >= 3 && done) {
      let best = pts[0];
      for (const p of pts) { if (p.y < best.y) best = p; }
      const bx = toX(best.x);
      const by = toY(best.y);
      svg += `<line x1="${bx}" y1="${by}" x2="${bx}" y2="${mt + ph}" stroke="${s.color}" stroke-width="1" stroke-dasharray="3,3" opacity="0.5"/>`;
      svg += renderStar(bx, by, 6, s.color);
    }
  }

  // Highlight selected scan point
  if (done && selectedScanIdx >= 0 && selectedScanIdx < data.length) {
    const selPt = data[selectedScanIdx];
    const selX = toX(selPt.param);
    const hlColor = isDark() ? '#ffd700' : '#cc8800';
    svg += `<line x1="${selX}" y1="${mt}" x2="${selX}" y2="${mt + ph}" stroke="${hlColor}" stroke-width="1" stroke-dasharray="3,2" opacity="0.6"/>`;
    for (const s of series) {
      const v = selPt[s.key];
      if (v != null && isFinite(v)) {
        svg += `<circle cx="${selX}" cy="${toY(v)}" r="5" fill="${s.color}" stroke="${hlColor}" stroke-width="2"/>`;
      }
    }
  }

  // Legend
  const legendX = ml + 8;
  let legendY = mt + 14;
  for (const s of series) {
    svg += `<rect x="${legendX}" y="${legendY - 7}" width="10" height="3" rx="1" fill="${s.color}"/>`;
    svg += `<text x="${legendX + 14}" y="${legendY - 3}" font-size="9" font-weight="600" fill="${tc.dim}">${s.label}</text>`;
    legendY += 14;
  }

  // Title & axis labels
  svg += `<text x="${ml + pw / 2}" y="20" text-anchor="middle" font-size="12" font-weight="600" fill="${tc.titleSvg}">${t('opt.graphTitle')}</text>`;
  svg += `<text x="${ml + pw / 2}" y="${height - 4}" text-anchor="middle" font-size="10" fill="${tc.dim}">${xLabel}</text>`;
  svg += `<text x="14" y="${mt + ph / 2}" text-anchor="middle" font-size="10" fill="${tc.dim}" transform="rotate(-90,14,${mt + ph / 2})">${t('opt.yEnergy')}</text>`;

  svg += '</svg>';
  container.innerHTML = svg;
}

function renderStar(cx: number, cy: number, r: number, fill: string): string {
  const pts: string[] = [];
  for (let i = 0; i < 10; i++) {
    const angle = Math.PI / 2 + i * Math.PI / 5;
    const rad = i % 2 === 0 ? r : r * 0.4;
    pts.push(`${(cx + rad * Math.cos(angle)).toFixed(1)},${(cy - rad * Math.sin(angle)).toFixed(1)}`);
  }
  return `<polygon points="${pts.join(' ')}" fill="${fill}" stroke="${fill}" stroke-width="0.5"/>`;
}

// ── Inline styles ────────────────────────────────────────────────────

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
    .opt-category-header {
      font-size: 0.7rem; font-weight: 600; color: var(--color-text-dim);
      text-transform: uppercase; letter-spacing: 0.03em;
      margin-top: 4px;
    }
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
    .opt-scenario-card span { display: none; }

    .opt-basis-row { display: flex; gap: 6px; }
    .opt-basis-btn {
      flex: 1; padding: 6px 0; border: 1px solid var(--color-border); border-radius: 6px;
      background: none; cursor: pointer; font-size: 0.78rem; color: var(--color-text);
      transition: all 0.15s;
    }
    .opt-basis-btn:hover { background: var(--color-surface-alt); }
    .opt-basis-btn.selected {
      border-color: var(--color-accent); color: var(--color-accent);
      font-weight: 600; box-shadow: 0 0 0 1px var(--color-accent);
    }

    .opt-method-row { display: flex; gap: 6px; }
    .opt-method-btn {
      flex: 1; padding: 7px 0; border: 2px solid var(--method-color); border-radius: 6px;
      background: none; cursor: pointer; font-size: 0.82rem; font-weight: 600;
      color: var(--method-color); transition: all 0.15s;
    }
    .opt-method-btn:hover { background: color-mix(in srgb, var(--method-color) 12%, transparent); }
    .opt-method-btn.active {
      background: var(--method-color); color: #fff;
    }

    .opt-charge-row {
      display: flex; gap: 10px; margin-top: 10px;
    }
    .opt-charge-row label {
      flex: 1; display: flex; align-items: center; gap: 6px; font-size: 0.72rem;
      color: var(--color-text-secondary);
    }
    .opt-charge-row input {
      width: 50px; padding: 4px 6px; border: 1px solid var(--color-border-input); border-radius: 5px;
      background: var(--color-input); color: var(--color-text); font-size: 0.8rem;
      outline: none; font-family: monospace; text-align: center;
    }

    .opt-param-form { display: flex; flex-direction: column; gap: 6px; }
    .opt-param-form label {
      display: flex; align-items: center; gap: 6px; font-size: 0.78rem;
    }
    .opt-param-form label span:first-child {
      flex: 0 0 50px; color: var(--color-text-secondary); font-size: 0.72rem;
    }
    .opt-param-form input {
      flex: 1; padding: 5px 8px; border: 1px solid var(--color-border-input); border-radius: 5px;
      background: var(--color-input); color: var(--color-text); font-size: 0.8rem;
      outline: none; font-family: monospace;
    }
    .opt-unit { font-size: 0.72rem; color: var(--color-text-dim); flex: 0 0 16px; }

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
    .opt-summary strong { color: var(--color-text); }
    .opt-dot {
      display: inline-block; width: 8px; height: 8px; border-radius: 50%;
      margin-right: 5px; vertical-align: middle;
    }

    #mol-vis { width: 100%; margin-bottom: 8px; }
    #graph-container { width: 100%; text-align: center; }

    .walsh-slider {
      width: 100%; max-width: 520px; margin: 8px auto 0;
      display: flex; align-items: center; gap: 10px;
    }
    .walsh-slider input[type=range] { flex: 1; cursor: pointer; }
    .walsh-slider-label {
      font-size: 0.85rem; font-weight: 600; min-width: 56px; text-align: center;
      color: var(--color-text);
    }
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
