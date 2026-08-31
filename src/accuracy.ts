/**
 * GANSU Lite — Accuracy comparison page.
 * Unifies Post-HF (RHF→MP2→MP3→CCSD) and DFT (SVWN→BLYP→PBE→B3LYP)
 * into a single page for side-by-side energy comparison.
 */

import type { FunctionalName } from './core/xcFunctional';
import type { AccWorkerRequest, AccWorkerResponse } from './core/accuracyWorker';
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
    id: 'beh2', labelKey: 'ladder.scenBeH2', descKey: 'ladder.descBeH2',
    charge: 0, mult: 1,
    xyz: '3\nBeH2\nBe  0.0  0.0  0.0\nH  0.0  0.0  1.334000\nH  0.0  0.0  -1.334000',
  },
  {
    id: 'lih', labelKey: 'ladder.scenLiH', descKey: 'ladder.descLiH',
    charge: 0, mult: 1,
    xyz: '2\nLiH\nLi  0.0  0.0  0.0\nH   0.0  0.0  1.596000',
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

// ── Basis set options ───────────────────────────────────────────────

const BASIS_OPTIONS = [
  'STO-3G', 'STO-6G', '3-21G', '6-31G', '6-31G(d,p)', 'cc-pVDZ', 'aug-cc-pVDZ', 'def2-SVP', 'cc-pVTZ', 'def2-TZVP',
];

// ── Method definitions ──────────────────────────────────────────────

interface MethodDef {
  id: string;
  label: string;
  category: string;
  type: 'hf' | 'postHF' | 'dft';
  functional?: FunctionalName;
}

const ALL_METHODS: MethodDef[] = [
  { id: 'HF',    label: 'HF',    category: 'Post-HF', type: 'hf' },
  { id: 'MP2',   label: 'MP2',   category: 'Post-HF', type: 'postHF' },
  { id: 'MP3',   label: 'MP3',   category: 'Post-HF', type: 'postHF' },
  { id: 'CCSD',  label: 'CCSD',  category: 'Post-HF', type: 'postHF' },
  { id: 'SVWN',  label: 'SVWN',  category: 'LDA',     type: 'dft', functional: 'SVWN' },
  { id: 'BLYP',  label: 'BLYP',  category: 'GGA',     type: 'dft', functional: 'BLYP' },
  { id: 'PBE',   label: 'PBE',   category: 'GGA',     type: 'dft', functional: 'PBE' },
  { id: 'B3LYP', label: 'B3LYP', category: 'Hybrid',  type: 'dft', functional: 'B3LYP' },
];

interface AccResult {
  method: MethodDef;
  totalEnergy: number;
  diffFromHF: number;   // energy - eHF (mEh)
  timeMs: number;
}

// ── State ────────────────────────────────────────────────────────────

let selectedMol = MOLECULES[0];
let selectedBasis = 'STO-3G';
let selectedMethods = new Set(ALL_METHODS.map(m => m.id)); // all enabled by default
let results: AccResult[] = [];
let running = false;
let stopRequested = false;
let done = false;

// Progress modal state
let progressMethodStates: ('pending' | 'running' | 'done' | 'skipped')[] = [];
let progressMethodTimes: number[] = [];
let progressMethodInfo: string[] = [];
let progressElapsed = 0;
let progressTimer: ReturnType<typeof setInterval> | null = null;
let progressStartTime = 0;
let currentBackend: 'wasm-simd' | 'wasm' | 'js' | 'none' = 'none';

// ── Basis loading (fetch GBS text for worker) ──────────────────────

const gbsCache = new Map<string, string>();

async function loadBasisGBS(name: string): Promise<string> {
  const cached = gbsCache.get(name);
  if (cached) return cached;
  const url = `${import.meta.env.BASE_URL}basis/${name.toLowerCase()}.gbs`;
  const resp = await fetch(url);
  if (!resp.ok) throw new Error(`Failed to load basis set: ${name}`);
  const text = await resp.text();
  gbsCache.set(name, text);
  return text;
}

// ── Method colors ────────────────────────────────────────────────────

const COLORS_LIGHT: Record<string, string> = {
  HF:    '#1f77b4',
  MP2:   '#2ca02c',
  MP3:   '#ff7f0e',
  CCSD:  '#d62728',
  SVWN:  '#9467bd',
  BLYP:  '#8c564b',
  PBE:   '#e377c2',
  B3LYP: '#17becf',
};
const COLORS_DARK: Record<string, string> = {
  HF:    '#4dabf7',
  MP2:   '#51cf66',
  MP3:   '#ffa94d',
  CCSD:  '#ff6b6b',
  SVWN:  '#b197fc',
  BLYP:  '#c4956a',
  PBE:   '#f59fda',
  B3LYP: '#66d9e8',
};

function methodColor(id: string): string {
  return isDark() ? (COLORS_DARK[id] ?? '#ccc') : (COLORS_LIGHT[id] ?? '#666');
}

// ── Progress Modal ───────────────────────────────────────────────────

function getActiveMethods(): MethodDef[] {
  return ALL_METHODS.filter(m => selectedMethods.has(m.id));
}

function showProgressModal(): void {
  removeProgressModal();

  progressMethodStates = ALL_METHODS.map(m => selectedMethods.has(m.id) ? 'pending' : 'skipped');
  progressMethodTimes = ALL_METHODS.map(() => 0);
  progressMethodInfo = ALL_METHODS.map(m => selectedMethods.has(m.id) ? '' : 'skipped');
  progressStartTime = performance.now();
  progressElapsed = 0;
  currentBackend = 'none';

  const overlay = document.createElement('div');
  overlay.id = 'acc-progress-overlay';
  overlay.className = 'acc-prog-overlay';

  const modal = document.createElement('div');
  modal.className = 'acc-prog-modal';
  modal.id = 'acc-progress-modal';

  overlay.appendChild(modal);
  document.body.appendChild(overlay);

  // Start elapsed timer — only update the elapsed text, not full innerHTML
  progressTimer = setInterval(() => {
    progressElapsed = performance.now() - progressStartTime;
    const el = document.getElementById('acc-prog-elapsed');
    if (el) el.textContent = (progressElapsed / 1000).toFixed(1) + 's';
  }, 200);

  renderProgressModal();
}

function renderProgressModal(): void {
  const modal = document.getElementById('acc-progress-modal');
  if (!modal) return;

  const doneCount = progressMethodStates.filter(s => s === 'done' || s === 'skipped').length;
  const pct = (doneCount / ALL_METHODS.length * 100).toFixed(0);
  const elapsedStr = (progressElapsed / 1000).toFixed(1);

  const simdCls = currentBackend === 'wasm-simd' ? 'active' : 'inactive';
  const wasmCls = currentBackend === 'wasm' ? 'active' : 'inactive';
  const jsCls = currentBackend === 'js' ? 'active' : 'inactive';

  let html = `
    <div class="acc-prog-header">
      <span class="acc-prog-title">${t('acc.running')}</span>
      <div class="pt-badges">
        <span class="pt-badge simd ${simdCls}">SIMD</span>
        <span class="pt-badge ${wasmCls}">WASM</span>
        <span class="pt-badge ${jsCls}">JS</span>
      </div>
      <span class="acc-prog-elapsed" id="acc-prog-elapsed">${elapsedStr}s</span>
    </div>
    <div class="acc-prog-bar-wrap">
      <div class="acc-prog-bar" style="width:${pct}%"></div>
    </div>
    <div class="acc-prog-list">`;

  for (let i = 0; i < ALL_METHODS.length; i++) {
    const m = ALL_METHODS[i];
    const st = progressMethodStates[i];
    const color = methodColor(m.id);

    let icon: string;
    let statusText: string;
    let rowClass: string;

    if (st === 'done') {
      const tStr = progressMethodTimes[i] < 1000
        ? `${progressMethodTimes[i].toFixed(0)}ms`
        : `${(progressMethodTimes[i] / 1000).toFixed(1)}s`;
      icon = '<span class="acc-prog-check">&#10003;</span>';
      statusText = tStr;
      rowClass = 'acc-prog-done';
    } else if (st === 'running') {
      icon = '<span class="acc-prog-spinner"></span>';
      statusText = progressMethodInfo[i] || '...';
      rowClass = 'acc-prog-active';
    } else if (st === 'skipped') {
      icon = '<span class="acc-prog-skip">&#8212;</span>';
      statusText = 'stopped';
      rowClass = 'acc-prog-skipped';
    } else {
      icon = '<span class="acc-prog-pending">&#9675;</span>';
      statusText = '';
      rowClass = 'acc-prog-pending-row';
    }

    // Group separator line
    if (i === 4) {
      html += `<div class="acc-prog-sep"></div>`;
    }

    html += `
      <div class="acc-prog-row ${rowClass}">
        <span class="acc-prog-icon">${icon}</span>
        <span class="acc-prog-label" style="color:${color}">${m.label}</span>
        <span class="acc-prog-cat">${m.category}</span>
        <span class="acc-prog-status">${statusText}</span>
      </div>`;
  }

  html += `</div>
    <button id="acc-prog-cancel" class="acc-prog-cancel-btn">${t('acc.stop')}</button>`;

  modal.innerHTML = html;

  document.getElementById('acc-prog-cancel')?.addEventListener('click', () => {
    stopRequested = true;
    if (activeWorker) {
      activeWorker.terminate();
      activeWorker = null;
    }
    markRemaining(0);
    finish();
  });
}

function setMethodState(idx: number, state: 'pending' | 'running' | 'done' | 'skipped', info?: string): void {
  progressMethodStates[idx] = state;
  if (info !== undefined) progressMethodInfo[idx] = info;
  renderProgressModal();
}

function setMethodTime(idx: number, ms: number): void {
  progressMethodTimes[idx] = ms;
}

function setMethodInfo(idx: number, msg: string): void {
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
  else if (msg.includes('SCF iteration')) short = msg.replace(/^.*?(SCF iteration)/, '$1');
  progressMethodInfo[idx] = short;
  renderProgressModal();
}

function removeProgressModal(): void {
  if (progressTimer) {
    clearInterval(progressTimer);
    progressTimer = null;
  }
  document.getElementById('acc-progress-overlay')?.remove();
}

// ── Main render ──────────────────────────────────────────────────────

const root = document.getElementById('app')!;

function render(): void {
  const canUseTZ = selectedMol.id === 'h2';
  if (!canUseTZ && selectedBasis === 'cc-pVTZ') selectedBasis = 'cc-pVDZ';
  const basisOptions = BASIS_OPTIONS.map(b => {
    const dis = (b === 'cc-pVTZ' && !canUseTZ) ? 'disabled' : '';
    return `<option value="${b}" ${b === selectedBasis ? 'selected' : ''} ${dis}>${b}${dis ? ' (H\u2082 only)' : ''}</option>`;
  }).join('');

  root.innerHTML = `
    <div class="opt-page">
      ${renderHeader('accuracy')}

      <div class="opt-content">
        <div class="opt-panel opt-controls">
          <h2>${t('acc.molecule')}</h2>
          <div class="opt-scenario-grid" id="scen-grid"></div>

          <h2>${t('acc.basis')}</h2>
          <select id="basis-select" class="acc-basis-select" ${running ? 'disabled' : ''}>
            ${basisOptions}
          </select>

          <h2>${t('acc.colMethod')}</h2>
          <div class="acc-method-toggles" id="method-toggles"></div>

          <button id="run-btn" class="opt-run-btn" ${running ? 'disabled' : ''}>
            ${running ? t('acc.running') : t('acc.run')}
          </button>

          <div id="summary-area"></div>
        </div>

        <div class="opt-panel opt-graph-panel">
          <div id="graph-container">
            ${!done && !running
              ? `<p class="opt-hint">${t('acc.waiting')}</p>`
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

  // Basis select
  root.querySelector('#basis-select')!.addEventListener('change', (e) => {
    selectedBasis = (e.target as HTMLSelectElement).value;
    render();
  });

  // Method toggle buttons
  const toggles = root.querySelector('#method-toggles')!;
  const postHFGroup = document.createElement('div');
  postHFGroup.className = 'acc-method-group';
  postHFGroup.innerHTML = '<span class="acc-method-group-label">Post-HF</span>';
  const postHFBtns = document.createElement('div');
  postHFBtns.className = 'acc-method-group-btns';
  postHFGroup.appendChild(postHFBtns);
  const dftGroup = document.createElement('div');
  dftGroup.className = 'acc-method-group';
  dftGroup.innerHTML = '<span class="acc-method-group-label">DFT</span>';
  const dftBtns = document.createElement('div');
  dftBtns.className = 'acc-method-group-btns';
  dftGroup.appendChild(dftBtns);

  for (const m of ALL_METHODS) {
    const btn = document.createElement('button');
    const isHF = m.type === 'hf';
    const isOn = selectedMethods.has(m.id);
    btn.className = 'acc-method-btn' + (isOn ? ' active' : '') + (isHF ? ' always-on' : '');
    btn.textContent = m.label;
    btn.disabled = running || isHF;
    btn.style.color = isOn ? methodColor(m.id) : '';
    btn.title = isHF ? 'HF is always required' : '';
    btn.addEventListener('click', () => {
      if (running || isHF) return;
      if (selectedMethods.has(m.id)) {
        selectedMethods.delete(m.id);
      } else {
        selectedMethods.add(m.id);
      }
      render();
    });
    if (m.type === 'hf' || m.type === 'postHF') {
      postHFBtns.appendChild(btn);
    } else {
      dftBtns.appendChild(btn);
    }
  }
  toggles.appendChild(postHFGroup);
  toggles.appendChild(dftGroup);

  // Event listeners
  root.querySelector('#nav-theme')!.addEventListener('click', () => { toggleTheme(); render(); });
  root.querySelector('#nav-lang')!.addEventListener('click', () => { toggleLang(); render(); });
  root.querySelector('#run-btn')!.addEventListener('click', () => { if (!running) runAccuracy(); });

  // Render charts if data exists
  if (results.length > 0) {
    renderGraph();
    renderSummary();
  }
}

// ── Computation (Web Worker) ─────────────────────────────────────────

let activeWorker: Worker | null = null;

async function runAccuracy(): Promise<void> {
  running = true;
  stopRequested = false;
  results = [];
  done = false;
  render();
  showProgressModal();

  const mol = selectedMol;

  try {
    const basisGBS = await loadBasisGBS(selectedBasis);

    const worker = new Worker(
      new URL('./core/accuracyWorker.ts', import.meta.url),
      { type: 'module' },
    );
    activeWorker = worker;

    worker.postMessage({
      type: 'run',
      xyzText: mol.xyz,
      basisGBS,
      charge: mol.charge,
      multiplicity: mol.mult,
      selectedMethods: Array.from(selectedMethods),
      baseUrl: import.meta.env.BASE_URL,
    } satisfies AccWorkerRequest);

    await new Promise<void>((resolve, reject) => {
      worker.onmessage = (e: MessageEvent<AccWorkerResponse>) => {
        const msg = e.data;
        switch (msg.type) {
          case 'method-state':
            setMethodState(msg.index, msg.state, msg.info);
            break;
          case 'method-info':
            setMethodInfo(msg.index, msg.msg);
            break;
          case 'method-time':
            setMethodTime(msg.index, msg.ms);
            break;
          case 'backend':
            currentBackend = msg.backend;
            renderProgressModal();
            break;
          case 'method-result':
            results.push({
              method: ALL_METHODS[msg.methodIndex],
              totalEnergy: msg.totalEnergy,
              diffFromHF: msg.diffFromHF,
              timeMs: msg.timeMs,
            });
            renderGraph();
            renderSummary();
            break;
          case 'done':
            worker.terminate();
            activeWorker = null;
            resolve();
            break;
          case 'error':
            worker.terminate();
            activeWorker = null;
            reject(new Error(msg.message));
            break;
        }
      };
      worker.onerror = (e) => {
        worker.terminate();
        activeWorker = null;
        reject(new Error(e.message || 'Worker failed'));
      };
    });
  } catch (err) {
    if (err instanceof Error && err.message === 'Cancelled') {
      // user cancelled
    } else {
      console.error('Accuracy computation error:', err);
    }
  }

  finish();
}

function markRemaining(fromIdx: number): void {
  for (let i = fromIdx; i < ALL_METHODS.length; i++) {
    if (progressMethodStates[i] === 'pending' || progressMethodStates[i] === 'running') {
      setMethodState(i, 'skipped');
    }
  }
}

function finish(): void {
  running = false;
  done = true;
  // Keep modal visible briefly to show final state, then remove
  setTimeout(() => {
    removeProgressModal();
    render();
  }, 600);
}

// ── Summary ──────────────────────────────────────────────────────────

function renderSummary(): void {
  const el = root.querySelector('#summary-area');
  if (!el || results.length === 0) return;

  const totalTime = results.reduce((s, r) => s + r.timeMs, 0);

  let html = `<div class="opt-summary">
    <h3>${done ? t('acc.done') : t('acc.running')}</h3>
    <table>
      <tr>
        <th>${t('acc.colMethod')}</th>
        <th>${t('acc.colCategory')}</th>
        <th>${t('acc.colEnergy')}</th>
        <th>${t('acc.colDiff')}</th>
        <th>${t('acc.colTime')}</th>
      </tr>`;

  for (const r of results) {
    const eStr = r.totalEnergy.toFixed(6);
    const dStr = r.method.type === 'hf' ? '\u2014' : `${r.diffFromHF.toFixed(2)} mEh`;
    const tStr = r.timeMs < 1000 ? `${r.timeMs.toFixed(0)}ms` : `${(r.timeMs / 1000).toFixed(1)}s`;
    const color = methodColor(r.method.id);
    html += `<tr>
      <td><span style="color:${color};font-weight:600">${r.method.label}</span></td>
      <td style="font-size:0.65rem;color:var(--color-text-dim)">${r.method.category}</td>
      <td style="font-family:monospace">${eStr}</td>
      <td style="font-family:monospace">${dStr}</td>
      <td style="text-align:right">${tStr}</td>
    </tr>`;
  }

  html += `</table>
    ${done ? `<div style="margin-top:6px;font-size:0.72rem;color:var(--color-text-dim)">${t('acc.totalTime')}: ${(totalTime / 1000).toFixed(1)}s</div>` : ''}
  </div>`;
  el.innerHTML = html;
}

// ── SVG Charts ───────────────────────────────────────────────────────

function renderGraph(): void {
  const container = root.querySelector('#graph-container');
  if (!container || results.length === 0) return;

  container.innerHTML = renderEnergyChart() + renderDiffChart();
}

/** Chart 1: Energy level comparison — all 8 methods with group separator */
function renderEnergyChart(): string {
  const tc = getThemeColors();

  const width = 680, height = 300;
  const ml = 82, mr = 20, mt = 36, mb = 52;
  const pw = width - ml - mr, ph = height - mt - mb;

  const energies = results.map(r => r.totalEnergy);
  const eMin = Math.min(...energies);
  const eMax = Math.max(...energies);
  const ePad = Math.max((eMax - eMin) * 0.15, 0.002);
  const eLow = eMin - ePad;
  const eHigh = eMax + ePad;

  const toY = (e: number) => mt + ph - ((e - eLow) / (eHigh - eLow)) * ph;

  const n = results.length;
  const barW = pw / ALL_METHODS.length;

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

  // Group separator: vertical dashed line between CCSD and SVWN (between index 3 and 4)
  const sepX = ml + barW * 4;
  svg += `<line x1="${sepX}" y1="${mt}" x2="${sepX}" y2="${mt + ph}" stroke="${tc.grid}" stroke-width="1" stroke-dasharray="4,3"/>`;

  // Group labels at the bottom
  const postHFCenter = ml + barW * 2;
  const dftCenter = ml + barW * 6;
  svg += `<text x="${postHFCenter}" y="${height - 4}" text-anchor="middle" font-size="9" fill="${tc.dim}" font-weight="600">Post-HF</text>`;
  svg += `<text x="${dftCenter}" y="${height - 4}" text-anchor="middle" font-size="9" fill="${tc.dim}" font-weight="600">DFT</text>`;

  // Axes
  svg += `<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt + ph}" stroke="${tc.axis}" stroke-width="1"/>`;
  svg += `<line x1="${ml}" y1="${mt + ph}" x2="${ml + pw}" y2="${mt + ph}" stroke="${tc.axis}" stroke-width="1"/>`;

  // Draw bars and energy levels for each method
  const hfEnergy = results[0].totalEnergy;
  for (let i = 0; i < n; i++) {
    const r = results[i];
    const color = methodColor(r.method.id);
    const cx = ml + barW * i + barW / 2;
    const y = toY(r.totalEnergy);

    // Filled bar from HF energy to this method's energy (only for non-HF)
    if (r.method.type !== 'hf') {
      const yHF = toY(hfEnergy);
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
    svg += `<text x="${cx}" y="${mt + ph + 14}" text-anchor="middle" font-size="9" fill="${tc.dim}" font-weight="600">${r.method.label}</text>`;

    // Category sub-label for DFT methods
    if (r.method.type === 'dft') {
      svg += `<text x="${cx}" y="${mt + ph + 25}" text-anchor="middle" font-size="7" fill="${tc.dim}" opacity="0.7">${r.method.category}</text>`;
    }

    // Connecting dashed line within same group
    if (i > 0 && i !== 4) {
      const prevY = toY(results[i - 1].totalEnergy);
      const prevCx = ml + barW * (i - 1) + barW / 2;
      svg += `<line x1="${prevCx}" y1="${prevY}" x2="${cx}" y2="${y}" stroke="${tc.grid}" stroke-width="1" stroke-dasharray="4,3"/>`;
    }
  }

  // Title & Y-axis label
  svg += `<text x="${ml + pw / 2}" y="20" text-anchor="middle" font-size="12" font-weight="600" fill="${tc.titleSvg}">${t('acc.graphTitle')}</text>`;
  svg += `<text x="14" y="${mt + ph / 2}" text-anchor="middle" font-size="10" fill="${tc.dim}" transform="rotate(-90,14,${mt + ph / 2})">${t('acc.yEnergy')}</text>`;

  svg += '</svg>';
  return svg;
}

/** Chart 2: ΔE from HF bar chart — all non-HF methods */
function renderDiffChart(): string {
  const tc = getThemeColors();

  const nonHF = results.filter(r => r.method.type !== 'hf');
  if (nonHF.length < 1) return '';

  const width = 680, height = 200;
  const ml = 82, mr = 20, mt = 28, mb = 48;
  const pw = width - ml - mr, ph = height - mt - mb;

  const n = nonHF.length;
  const barW = pw / n;

  // ΔE values in mEh
  const diffs = nonHF.map(r => r.diffFromHF);
  const maxAbs = Math.max(Math.abs(Math.min(...diffs)), Math.abs(Math.max(...diffs)), 1);
  const scale = maxAbs * 1.2;

  const toY = (mEh: number) => mt + ph / 2 - (mEh / scale) * (ph / 2);

  let svg = `<svg width="${width}" height="${height}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" style="max-width:100%;margin-top:8px;">`;

  // Background
  svg += `<rect x="${ml}" y="${mt}" width="${pw}" height="${ph}" fill="${tc.surface}" rx="2"/>`;

  // Group separator between Post-HF and DFT (after CCSD = index 2 in nonHF)
  const postHFCount = nonHF.filter(r => r.method.type === 'postHF').length;
  if (postHFCount > 0 && postHFCount < n) {
    const gSepX = ml + barW * postHFCount;
    svg += `<line x1="${gSepX}" y1="${mt}" x2="${gSepX}" y2="${mt + ph}" stroke="${tc.grid}" stroke-width="1" stroke-dasharray="4,3"/>`;
  }

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
    const r = nonHF[i];
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
    svg += `<text x="${cx}" y="${mt + ph + 14}" text-anchor="middle" font-size="9" fill="${tc.dim}" font-weight="600">${r.method.label}</text>`;

    // Category sub-label for DFT
    if (r.method.type === 'dft') {
      svg += `<text x="${cx}" y="${mt + ph + 25}" text-anchor="middle" font-size="7" fill="${tc.dim}" opacity="0.7">${r.method.category}</text>`;
    }
  }

  // Group labels at the bottom
  if (postHFCount > 0 && postHFCount < n) {
    const phCenter = ml + barW * postHFCount / 2;
    const dCenter = ml + barW * (postHFCount + (n - postHFCount) / 2);
    svg += `<text x="${phCenter}" y="${height - 4}" text-anchor="middle" font-size="9" fill="${tc.dim}" font-weight="600">Post-HF</text>`;
    svg += `<text x="${dCenter}" y="${height - 4}" text-anchor="middle" font-size="9" fill="${tc.dim}" font-weight="600">DFT</text>`;
  }

  // Title & Y-axis
  svg += `<text x="${ml + pw / 2}" y="16" text-anchor="middle" font-size="11" font-weight="600" fill="${tc.titleSvg}">${t('acc.diffTitle')}</text>`;
  svg += `<text x="14" y="${mt + ph / 2}" text-anchor="middle" font-size="9" fill="${tc.dim}" transform="rotate(-90,14,${mt + ph / 2})">${t('acc.yDiff')}</text>`;

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
    .opt-page { max-width: 1040px; margin: 0 auto; padding: 16px 20px; }

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

    .acc-basis-select {
      width: 100%; padding: 6px 8px; border: 1px solid var(--color-border);
      border-radius: 6px; background: var(--color-surface);
      color: var(--color-text); font-size: 0.82rem;
    }

    .opt-run-btn {
      width: 100%; margin-top: 14px; padding: 10px; border: none; border-radius: 8px;
      font-size: 0.85rem; font-weight: 600; cursor: pointer;
      background: var(--color-accent); color: var(--color-accent-on);
      transition: background 0.15s;
    }
    .opt-run-btn:hover:not([disabled]) { background: var(--color-accent-hover); }
    .opt-run-btn[disabled] { opacity: 0.6; cursor: not-allowed; }

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

    /* ── Progress modal ── */
    .acc-prog-overlay {
      position: fixed; inset: 0; z-index: 9000;
      background: rgba(0,0,0,0.45);
      display: flex; align-items: center; justify-content: center;
    }
    .acc-prog-modal {
      background: var(--color-surface, #fff); border: 1px solid var(--color-border, #ddd);
      border-radius: 12px; padding: 20px 24px; width: 380px; max-width: 92vw;
      box-shadow: 0 8px 32px rgba(0,0,0,0.25);
    }
    .acc-prog-header {
      display: flex; justify-content: space-between; align-items: center;
      margin-bottom: 10px;
    }
    .acc-prog-title { font-size: 0.9rem; font-weight: 600; color: var(--color-text); }
    .acc-prog-elapsed { font-size: 0.78rem; color: var(--color-text-dim); font-family: monospace; }

    .acc-prog-bar-wrap {
      height: 6px; background: var(--color-progress-bg, #e0e4ea);
      border-radius: 3px; margin-bottom: 14px; overflow: hidden;
    }
    .acc-prog-bar {
      height: 100%; background: var(--color-accent, #4a90d9);
      border-radius: 3px; transition: width 0.3s ease;
    }

    .acc-prog-list { display: flex; flex-direction: column; gap: 3px; }

    .acc-prog-sep {
      height: 1px; background: var(--color-border, #ddd);
      margin: 4px 0;
    }

    .acc-prog-row {
      display: flex; align-items: center; gap: 8px;
      padding: 4px 6px; border-radius: 5px; font-size: 0.78rem;
    }
    .acc-prog-active {
      background: var(--color-surface-alt, #f5f5f5);
    }
    .acc-prog-icon { width: 18px; text-align: center; flex-shrink: 0; }
    .acc-prog-label { font-weight: 600; width: 50px; }
    .acc-prog-cat { font-size: 0.65rem; color: var(--color-text-dim); width: 52px; }
    .acc-prog-status {
      flex: 1; text-align: right; font-size: 0.7rem;
      color: var(--color-text-dim); font-family: monospace;
      overflow: hidden; text-overflow: ellipsis; white-space: nowrap;
    }

    .acc-prog-check { color: var(--color-converged, #2ecc71); font-weight: 700; }
    .acc-prog-skip { color: var(--color-text-dim, #999); }
    .acc-prog-pending { color: var(--color-text-dim, #ccc); font-size: 0.6rem; }
    .acc-prog-pending-row { opacity: 0.5; }
    .acc-prog-skipped { opacity: 0.4; }

    .acc-prog-spinner {
      display: inline-block; width: 14px; height: 14px;
      border: 2px solid var(--color-border, #ddd);
      border-top-color: var(--color-accent, #4a90d9);
      border-radius: 50%;
      animation: acc-spin 0.8s linear infinite;
    }
    @keyframes acc-spin { to { transform: rotate(360deg); } }

    .acc-prog-cancel-btn {
      display: block; width: 100%; margin-top: 14px;
      padding: 8px; border: 1px solid var(--color-error, #e74c3c);
      border-radius: 8px; font-size: 0.8rem; font-weight: 600;
      cursor: pointer; background: none;
      color: var(--color-error, #e74c3c);
    }
    .acc-prog-cancel-btn:hover { background: rgba(231,76,60,0.08); }

    /* ── Method toggle buttons ── */
    .acc-method-toggles { display: flex; flex-direction: column; gap: 5px; }
    .acc-method-group { display: flex; align-items: center; gap: 4px; }
    .acc-method-group-label {
      font-size: 0.6rem; color: var(--color-text-dim); font-weight: 500;
      letter-spacing: 0.03em; width: 40px; flex-shrink: 0;
    }
    .acc-method-group-btns { display: flex; gap: 3px; flex: 1; }
    .acc-method-btn {
      padding: 3px 0; border: 1px solid var(--color-border);
      border-radius: 4px; font-size: 0.7rem; font-weight: 600;
      cursor: pointer; background: var(--color-surface);
      color: var(--color-text-dim); transition: all 0.15s;
      flex: 1; text-align: center; min-width: 0;
    }
    .acc-method-btn:hover:not([disabled]) { background: var(--color-surface-alt); }
    .acc-method-btn.active {
      border-color: currentColor; background: var(--color-surface-alt);
      box-shadow: 0 0 0 1px currentColor;
    }
    .acc-method-btn.always-on {
      opacity: 0.7; cursor: default;
    }
    .acc-method-btn[disabled]:not(.always-on) { opacity: 0.5; cursor: not-allowed; }

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
