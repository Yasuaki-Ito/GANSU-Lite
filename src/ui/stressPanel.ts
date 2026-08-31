/**
 * Cross-device stress test: push a monotonically growing molecule series until
 * the device gives up, and record *how* it gave up.
 *
 * Each ladder point runs in a throw-away Web Worker, so a failure is observed
 * rather than inferred: an allocation failure comes back as an error message,
 * a job that outlives its budget is terminated as a timeout, and a worker the
 * engine kills is reported as such. Results are checkpointed to localStorage
 * after every point, so even a whole-tab crash leaves a usable record.
 */

import { gatherDeviceInfo, type DeviceInfo } from '../core/deviceInfo';
import { STRESS_SERIES, type StressSeries } from '../core/stressGeometries';
import type { StressRequest, StressResponse, MemoryProbeRequest, ProbeVerdict } from '../core/stressWorker';

const STORAGE_KEY = 'gansu-stress-v1';
const PROBE_KEY = 'gansu-stress-probe-v1';
const DEFAULT_BUDGET_S = 180;
const BASIS_NAME = '6-31g(d,p)';

export type StressStatus = 'ok' | 'timeout' | 'error' | 'worker-died' | 'tab-crash';

export interface StressPoint {
  n: number;
  name: string;
  natoms: number;
  nbasis: number | null;
  status: StressStatus;
  totalMs: number | null;
  scfMs: number | null;
  iterations: number | null;
  converged: boolean | null;
  energy: number | null;
  eriBytes: number | null;
  heapUsedBytes: number | null;
  heapLimitBytes: number | null;
  detail?: string;
}

export interface StressSummary {
  largestCompleted: { name: string; natoms: number; nbasis: number; seconds: number } | null;
  firstFailure: { name: string; natoms: number; nbasis: number | null; mode: string } | null;
  /**
   * Set when a timeout failure sits close enough to the budget that the reported
   * limit is an artefact of the budget rather than a property of the device.
   */
  borderline: { estimatedSeconds: number; budgetSeconds: number; exponent: number } | null;
}

export interface StressOutput {
  device: DeviceInfo;
  backend: string;
  seriesId: string;
  seriesLabel: string;
  method: string;
  basis: string;
  budgetSeconds: number;
  points: StressPoint[];
  summary: StressSummary;
  memoryProbe: ProbeRow[];
  toolUrl: string;
}

interface Checkpoint {
  seriesId: string;
  method: string;
  budgetSeconds: number;
  points: StressPoint[];
  /** Non-null while a point is being computed — a leftover means the tab died. */
  inFlight: { n: number; name: string; natoms: number; startedAt: number } | null;
}

// ── Persistence ──

function loadCheckpoint(): Checkpoint | null {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    return raw ? JSON.parse(raw) as Checkpoint : null;
  } catch { return null; }
}

function saveCheckpoint(cp: Checkpoint): void {
  try { localStorage.setItem(STORAGE_KEY, JSON.stringify(cp)); } catch { /* quota / private mode */ }
}

function clearCheckpoint(): void {
  try { localStorage.removeItem(STORAGE_KEY); } catch { /* ignore */ }
}

// ── Formatting ──

function humanBytes(b: number | null): string {
  if (b == null) return '—';
  if (b < 1024 ** 2) return `${(b / 1024).toFixed(0)} KB`;
  if (b < 1024 ** 3) return `${(b / 1024 ** 2).toFixed(0)} MB`;
  return `${(b / 1024 ** 3).toFixed(2)} GB`;
}

function secs(ms: number | null): string {
  return ms == null ? '—' : (ms / 1000).toFixed(2);
}

const STATUS_LABEL: Record<StressStatus, string> = {
  'ok': 'ok',
  'timeout': 'timeout (exceeded budget)',
  'error': 'error / OOM',
  'worker-died': 'worker killed (OOM)',
  'tab-crash': 'tab crash',
};

// ── One ladder point ──

interface PointOutcome {
  status: StressStatus;
  nbasis: number | null;
  totalMs: number | null;
  scfMs: number | null;
  iterations: number | null;
  converged: boolean | null;
  energy: number | null;
  eriBytes: number | null;
  heapUsedBytes: number | null;
  heapLimitBytes: number | null;
  detail?: string;
}

function runPoint(
  xyzText: string,
  basisGBS: string,
  budgetMs: number,
  dft: boolean,
  onProgress: (msg: string, elapsedMs: number) => void,
  registerWorker: (w: Worker | null) => void,
): Promise<PointOutcome> {
  return new Promise((resolve) => {
    const worker = new Worker(
      new URL('../core/stressWorker.ts', import.meta.url),
      { type: 'module' },
    );
    registerWorker(worker);

    let settled = false;
    const t0 = performance.now();
    // A worker the engine kills for OOM often just goes quiet instead of firing
    // an error event. Tracking the last sign of life lets the timeout branch say
    // whether the job was still working or had already stopped talking.
    let lastMessageAt = t0;
    let lastMessage = 'no progress reported';
    // The worker announces the size before the expensive work starts, so a point
    // that later times out or dies still knows how big it was.
    let observedNbasis: number | null = null;

    const finish = (o: PointOutcome) => {
      if (settled) return;
      settled = true;
      clearTimeout(timer);
      worker.terminate();
      registerWorker(null);
      resolve(o);
    };

    const eriBytesFor = (n: number) => {
      const npair = n * (n + 1) / 2;
      return (npair * (npair + 1) / 2) * 8;
    };
    const partial = () => ({
      nbasis: observedNbasis,
      totalMs: null, scfMs: null, iterations: null, converged: null, energy: null,
      eriBytes: observedNbasis == null ? null : eriBytesFor(observedNbasis),
      heapUsedBytes: null, heapLimitBytes: null,
    });

    const timer = setTimeout(() => {
      const silentS = (performance.now() - lastMessageAt) / 1000;
      finish({
        ...partial(),
        status: 'timeout',
        totalMs: performance.now() - t0,
        detail:
          `exceeded ${(budgetMs / 1000).toFixed(0)} s budget; ` +
          `silent for ${silentS.toFixed(0)} s after "${lastMessage}"` +
          (silentS > 30 ? ' — a long silence usually means the worker was killed for memory' : ''),
      });
    }, budgetMs);

    worker.onmessage = (e: MessageEvent<StressResponse>) => {
      const m = e.data;
      lastMessageAt = performance.now();
      if (m.type === 'progress') {
        const hit = /(\d+) basis functions/.exec(m.message);
        if (hit) observedNbasis = Number(hit[1]);
        lastMessage = m.message;
        onProgress(m.message, m.elapsedMs);
        return;
      }
      if (m.type === 'done') {
        finish({
          status: 'ok',
          nbasis: m.nbasis, totalMs: m.totalMs, scfMs: m.scfMs,
          iterations: m.iterations, converged: m.converged, energy: m.energy,
          eriBytes: m.eriBytes, heapUsedBytes: m.heapUsedBytes, heapLimitBytes: m.heapLimitBytes,
          ...(m.converged ? {} : { detail: 'SCF did not converge' }),
        });
        return;
      }
      if (m.type === 'error') {
        finish({
          ...partial(),
          status: 'error',
          nbasis: m.nbasis ?? observedNbasis,
          totalMs: m.elapsedMs,
          detail: `${m.name}: ${m.message} (phase: ${m.phase})`,
        });
      }
      // probe-* messages belong to the memory probe, not to a ladder point
    };

    // Fires when the worker itself is torn down by the engine (typically OOM).
    worker.onerror = (ev) => {
      finish({
        ...partial(),
        status: 'worker-died',
        totalMs: performance.now() - t0,
        detail: ev.message || 'worker terminated by the browser',
      });
    };

    const req: StressRequest = {
      type: 'stress-run',
      xyzText,
      basisGBS,
      charge: 0,
      eriBackend: 'auto',
      baseUrl: import.meta.env.BASE_URL,
      ...(dft ? { dftConfig: { functional: 'B3LYP', gridLevel: 'medium' as const } } : {}),
    };
    worker.postMessage(req);
  });
}

// ── Memory-ceiling probe ──

export interface ProbeRow { nbasis: number; bytes: number; ok: boolean; ms: number; gbPerSec: number; verdict: ProbeVerdict; detail?: string; }

/**
 * Basis-function counts to probe: the water ladder plus a short tail. The worker
 * stops on its own at the first swapping or failing size, and refuses anything
 * past its own cap, so the tail only needs to reach past a plausible ceiling.
 */
function probeLadder(): number[] {
  const xs = new Set<number>();
  for (const n of STRESS_SERIES[0].ladder) xs.add(n * 25);
  // Stop at 300 basis functions (a 7.6 GB array). Past that the memory answer is
  // moot: the SCF itself would run for well over an hour at N⁴, so there is no
  // configuration in which holding the array is the thing standing in the way.
  return [...xs].filter(nb => nb <= 300).sort((a, b) => a - b);
}

function runMemoryProbe(onRow: (r: ProbeRow) => void): Promise<ProbeRow[]> {
  return new Promise((resolve) => {
    const worker = new Worker(new URL('../core/stressWorker.ts', import.meta.url), { type: 'module' });
    const rows: ProbeRow[] = [];
    const done = () => { worker.terminate(); resolve(rows); };
    worker.onmessage = (e: MessageEvent<StressResponse>) => {
      const m = e.data;
      if (m.type === 'probe-result') {
        const row: ProbeRow = {
          nbasis: m.nbasis, bytes: m.bytes, ok: m.ok, ms: m.ms,
          gbPerSec: m.gbPerSec, verdict: m.verdict,
          ...(m.detail ? { detail: m.detail } : {}),
        };
        rows.push(row);
        onRow(row);
        return;
      }
      if (m.type === 'probe-done') done();
    };
    // If the engine kills the worker mid-probe, the rows gathered so far still stand.
    worker.onerror = () => done();
    const req: MemoryProbeRequest = { type: 'memory-probe', nbasisLadder: probeLadder() };
    worker.postMessage(req);
  });
}

function appendProbeRow(body: HTMLTableSectionElement, r: ProbeRow): void {
  const good = r.verdict === 'resident';
  const label = {
    resident: 'resident in RAM',
    swapping: 'SWAPPING — past usable memory',
    failed: 'ALLOCATION FAILED',
    skipped: 'skipped (past cap)',
  }[r.verdict];
  const tr = document.createElement('tr');
  tr.innerHTML =
    `<td class="num">${r.nbasis}</td><td class="num">${humanBytes(r.bytes)}</td>` +
    `<td class="num">${r.ms.toFixed(0)}</td>` +
    `<td class="num">${r.gbPerSec > 0 ? r.gbPerSec.toFixed(2) : '—'}</td>` +
    `<td${good ? '' : ' class="cb-err"'}>${label}` +
    `${r.detail ? `<br><span class="cb-err">${escapeHtml(r.detail)}</span>` : ''}</td>`;
  body.appendChild(tr);
}

let probeRows: ProbeRow[] = [];

function saveProbeRows(rows: ProbeRow[]): void {
  try { localStorage.setItem(PROBE_KEY, JSON.stringify(rows)); } catch { /* quota / private mode */ }
}

function loadProbeRows(): ProbeRow[] {
  try {
    const raw = localStorage.getItem(PROBE_KEY);
    return raw ? JSON.parse(raw) as ProbeRow[] : [];
  } catch { return []; }
}

async function doProbe(): Promise<void> {
  const btn = document.getElementById('cb-stress-probe') as HTMLButtonElement | null;
  const body = document.getElementById('cb-probe-body') as HTMLTableSectionElement | null;
  const box = document.getElementById('cb-probe-section');
  if (btn) btn.disabled = true;
  if (box) box.style.display = '';
  if (body) body.innerHTML = '';
  probeRows = [];
  setStatus('probing memory ceiling…');

  await runMemoryProbe((r) => {
    probeRows.push(r);
    saveProbeRows(probeRows);
    if (body) appendProbeRow(body, r);
  });

  const resident = probeRows.filter(r => r.verdict === 'resident');
  const stopper = probeRows.find(r => r.verdict !== 'resident');
  const summary = document.getElementById('cb-probe-summary');
  if (summary) {
    const last = resident.length ? resident[resident.length - 1] : null;
    if (!last) {
      summary.textContent = 'Could not hold even the smallest probe array — something is wrong.';
    } else {
      const why = stopper
        ? { swapping: `at ${stopper.nbasis} the array no longer fits in real memory`,
            failed: `allocation failed outright at ${stopper.nbasis}`,
            skipped: `the ladder ran past the probe cap at ${stopper.nbasis}`,
            resident: '' }[stopper.verdict]
        : 'the probe ladder topped out at 300 basis functions without hitting a ceiling — ' +
          'past that the N⁴ SCF cost, not memory, is what stops you';
      summary.textContent =
        `Memory ceiling: ${stopper ? '' : '≥ '}${last.nbasis} basis functions — a ${humanBytes(last.bytes)} ERI array ` +
        `still resident in RAM (${last.gbPerSec.toFixed(1)} GB/s); ${why}.`;
    }
  }
  setStatus('memory probe done');
  if (btn) btn.disabled = false;
}

// ── Summary ──

/**
 * How long the first failing point would have needed, extrapolated from the last
 * two completed points. SCF cost here scales as roughly N^5 once memory traffic
 * is included, so the exponent is fitted from the data rather than assumed.
 */
function estimateFailingTime(ok: StressPoint[], bad: StressPoint): { seconds: number; exponent: number } | null {
  if (ok.length < 2 || bad.nbasis == null) return null;
  const a = ok[ok.length - 2], b = ok[ok.length - 1];
  if (a.nbasis == null || b.nbasis == null || a.totalMs == null || b.totalMs == null) return null;
  if (a.nbasis >= b.nbasis || b.nbasis >= bad.nbasis || a.totalMs <= 0) return null;
  const exponent = Math.log(b.totalMs / a.totalMs) / Math.log(b.nbasis / a.nbasis);
  if (!Number.isFinite(exponent) || exponent <= 0) return null;
  const ms = b.totalMs * Math.pow(bad.nbasis / b.nbasis, exponent);
  return { seconds: ms / 1000, exponent };
}

function summarize(points: StressPoint[], budgetSeconds: number): StressSummary {
  const ok = points.filter(p => p.status === 'ok');
  const last = ok.length ? ok[ok.length - 1] : null;
  const bad = points.find(p => p.status !== 'ok') ?? null;

  let borderline: StressSummary['borderline'] = null;
  if (bad && bad.status === 'timeout') {
    const est = estimateFailingTime(ok, bad);
    // Within 2x of the budget means a modestly larger budget would very likely
    // have let this point through — the rung is not really beyond the device.
    if (est && est.seconds < 2 * budgetSeconds) {
      borderline = { estimatedSeconds: est.seconds, budgetSeconds, exponent: est.exponent };
    }
  }

  return {
    borderline,
    largestCompleted: last && last.nbasis != null && last.totalMs != null
      ? { name: last.name, natoms: last.natoms, nbasis: last.nbasis, seconds: +(last.totalMs / 1000).toFixed(1) }
      : null,
    firstFailure: bad
      ? { name: bad.name, natoms: bad.natoms, nbasis: bad.nbasis, mode: STATUS_LABEL[bad.status] }
      : null,
  };
}

// ── Panel ──

let state: {
  series: StressSeries;
  points: StressPoint[];
  running: boolean;
  stopRequested: boolean;
  activeWorker: Worker | null;
  output: StressOutput | null;
} = {
  series: STRESS_SERIES[0],
  points: [],
  running: false,
  stopRequested: false,
  activeWorker: null,
  output: null,
};

export function stressPanelHTML(): string {
  return `
    <section class="cb-panel">
      <h2>Stress test — how large a system does this device survive?</h2>
      <p class="cb-note">
        Runs a monotonically growing series at RHF/6-31G(d,p) until the device fails, and records
        the failure mode. Each point runs in a fresh Web Worker with a wall-clock budget, so an
        out-of-memory abort or a run that is simply too slow is captured instead of taking the page
        down. Progress is checkpointed to this browser's local storage after every point — if the
        tab does die, reload and the lost point is reported as a tab crash.
      </p>
      <div id="cb-stress-recovered" class="cb-note" style="display:none;color:var(--color-error,#e05050)"></div>
      <div class="cb-stress-controls">
        <label>Series
          <select id="cb-stress-series">
            ${STRESS_SERIES.map(s => `<option value="${s.id}">${s.label}</option>`).join('')}
          </select>
        </label>
        <label>Theory
          <select id="cb-stress-theory">
            <option value="rhf">RHF/6-31G(d,p)</option>
            <option value="b3lyp">B3LYP/6-31G(d,p)</option>
          </select>
        </label>
        <label>Budget per point (s)
          <input id="cb-stress-budget" type="number" min="10" max="3600" step="10" value="${DEFAULT_BUDGET_S}">
        </label>
        <label>Ladder (n)
          <input id="cb-stress-ladder" type="text" value="${STRESS_SERIES[0].ladder.join(', ')}">
        </label>
      </div>
      <div style="margin:10px 0 8px">
        <button id="cb-stress-run">Run stress test</button>
        <button id="cb-stress-stop" disabled>Stop</button>
        <button id="cb-stress-reset">Clear saved progress</button>
        <span id="cb-stress-status" class="cb-wl-status" style="margin-left:10px"></span>
      </div>
      <div style="overflow-x:auto">
        <table class="cb-results-table" style="font-size:0.8rem">
          <thead><tr>
            <th>System</th><th class="num">Atoms</th><th class="num">Basis fns</th>
            <th class="num">ERI array</th><th class="num">Iters</th>
            <th class="num">Time (s)</th><th>Status</th>
          </tr></thead>
          <tbody id="cb-stress-body"></tbody>
        </table>
      </div>
      <h3 style="font-size:0.95rem;margin:16px 0 4px">Hard memory ceiling</h3>
      <p class="cb-note">
        The SCF's dominant allocation is the unique-ERI array (N⁴/8 doubles). Probing that
        allocation directly finds the memory ceiling in seconds, instead of waiting out an SCF that
        would take tens of minutes at the same size. Allocation success alone is not the test — a
        desktop OS will back a multi-gigabyte array with the page file and report success while
        thrashing — so the probe writes one value per page and stops when commit throughput
        collapses against the small-array baseline. Only sizes marked <em>resident in RAM</em> count.
        The probe stops at 300 basis functions and stays inside half of the reported device memory,
        so it can never push the machine into a paging storm.
      </p>
      <button id="cb-stress-probe">Probe memory ceiling</button>
      <div id="cb-probe-section" style="display:none;margin-top:10px">
        <table class="cb-results-table" style="font-size:0.8rem;max-width:520px">
          <thead><tr><th class="num">Basis fns</th><th class="num">ERI array</th><th class="num">ms</th><th class="num">GB/s</th><th>Result</th></tr></thead>
          <tbody id="cb-probe-body"></tbody>
        </table>
        <p id="cb-probe-summary" class="cb-note"></p>
      </div>

      <div id="cb-stress-summary" style="display:none">
        <h3 style="font-size:0.95rem;margin:14px 0 6px">Row for the paper table</h3>
        <pre id="cb-stress-row" class="cb-stress-row"></pre>
        <button id="cb-stress-copy-row">Copy row</button>
        <button id="cb-stress-copy-json">Copy JSON</button>
        <button id="cb-stress-download-json">Download JSON</button>
        <button id="cb-stress-download-csv">Download CSV</button>
      </div>
    </section>
  `;
}

export function wireStressPanel(getBackendLabel: () => string): void {
  const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T | null;

  const seriesSel = $<HTMLSelectElement>('cb-stress-series');
  const ladderInput = $<HTMLInputElement>('cb-stress-ladder');

  seriesSel?.addEventListener('change', () => {
    const s = STRESS_SERIES.find(x => x.id === seriesSel.value);
    if (s && ladderInput) { state.series = s; ladderInput.value = s.ladder.join(', '); }
  });

  $('cb-stress-run')?.addEventListener('click', () => { void runStress(getBackendLabel); });
  $('cb-stress-stop')?.addEventListener('click', () => {
    state.stopRequested = true;
    state.activeWorker?.terminate();
    setStatus('stopping…');
  });
  $('cb-stress-reset')?.addEventListener('click', () => {
    clearCheckpoint();
    try { localStorage.removeItem(PROBE_KEY); } catch { /* ignore */ }
    probeRows = [];
    state.points = [];
    state.output = null;
    renderRows();
    const rec = $('cb-stress-recovered');
    if (rec) rec.style.display = 'none';
    const sum = $('cb-stress-summary');
    if (sum) sum.style.display = 'none';
    setStatus('saved progress cleared');
  });

  $('cb-stress-probe')?.addEventListener('click', () => { void doProbe(); });
  $('cb-stress-copy-row')?.addEventListener('click', () => copy(paperRow()));
  $('cb-stress-copy-json')?.addEventListener('click', () => copy(JSON.stringify(state.output, null, 2)));
  $('cb-stress-download-json')?.addEventListener('click', () =>
    download('stress-test.json', JSON.stringify(state.output, null, 2), 'application/json'));
  $('cb-stress-download-csv')?.addEventListener('click', () =>
    download('stress-test.csv', toCSV(state.output), 'text/csv'));

  recoverCrashedRun();

  const saved = loadProbeRows();
  if (saved.length) {
    probeRows = saved;
    const box = document.getElementById('cb-probe-section');
    const body = document.getElementById('cb-probe-body') as HTMLTableSectionElement | null;
    if (box && body) { box.style.display = ''; for (const r of saved) appendProbeRow(body, r); }
  }
}

/** A leftover in-flight marker means the tab died mid-point last time. */
function recoverCrashedRun(): void {
  const cp = loadCheckpoint();
  if (!cp) return;
  state.points = cp.points ?? [];
  const s = STRESS_SERIES.find(x => x.id === cp.seriesId);
  if (s) state.series = s;

  if (cp.inFlight) {
    const f = cp.inFlight;
    if (!state.points.some(p => p.n === f.n)) {
      state.points.push({
        n: f.n, name: f.name, natoms: f.natoms, nbasis: null,
        status: 'tab-crash', totalMs: null, scfMs: null, iterations: null,
        converged: null, energy: null, eriBytes: null,
        heapUsedBytes: null, heapLimitBytes: null,
        detail: 'tab died while this point was running',
      });
    }
    const rec = document.getElementById('cb-stress-recovered');
    if (rec) {
      rec.style.display = '';
      rec.textContent =
        `Recovered from a previous session: the tab died while running ${f.name}. ` +
        `Results up to that point are restored below and ${f.name} is recorded as a tab crash.`;
    }
    saveCheckpoint({ ...cp, points: state.points, inFlight: null });
  }
  if (state.points.length) renderRows();
}

function setStatus(s: string): void {
  const el = document.getElementById('cb-stress-status');
  if (el) el.textContent = s;
}

async function fetchBasisGBS(): Promise<string> {
  const r = await fetch(`${import.meta.env.BASE_URL}basis/${BASIS_NAME}.gbs`);
  if (!r.ok) throw new Error(`Cannot load basis ${BASIS_NAME}: HTTP ${r.status}`);
  return r.text();
}

function parseLadder(raw: string, fallback: number[]): number[] {
  const xs = raw.split(/[,\s]+/).map(s => parseInt(s, 10)).filter(n => Number.isFinite(n) && n > 0);
  return xs.length ? xs : fallback;
}

async function runStress(getBackendLabel: () => string): Promise<void> {
  if (state.running) return;
  const runBtn = document.getElementById('cb-stress-run') as HTMLButtonElement | null;
  const stopBtn = document.getElementById('cb-stress-stop') as HTMLButtonElement | null;
  const budgetS = Math.max(10, Number((document.getElementById('cb-stress-budget') as HTMLInputElement)?.value) || DEFAULT_BUDGET_S);
  const dft = (document.getElementById('cb-stress-theory') as HTMLSelectElement)?.value === 'b3lyp';
  const ladder = parseLadder((document.getElementById('cb-stress-ladder') as HTMLInputElement)?.value ?? '', state.series.ladder);

  state.running = true;
  state.stopRequested = false;
  state.points = [];
  state.output = null;
  if (runBtn) runBtn.disabled = true;
  if (stopBtn) stopBtn.disabled = false;
  renderRows();

  let basisGBS: string;
  try {
    setStatus('loading basis set…');
    basisGBS = await fetchBasisGBS();
  } catch (e) {
    setStatus(`FAILED: ${e instanceof Error ? e.message : String(e)}`);
    state.running = false;
    if (runBtn) runBtn.disabled = false;
    if (stopBtn) stopBtn.disabled = true;
    return;
  }

  const series = state.series;
  const method = dft ? 'B3LYP' : 'RHF';
  const budgetMs = budgetS * 1000;

  // Warm-up: the first ladder point is run once and discarded so that WASM
  // instantiation and JIT warm-up do not land in a measured number.
  try {
    setStatus('warm-up…');
    await runPoint(series.xyz(ladder[0]), basisGBS, budgetMs, dft, () => {}, w => { state.activeWorker = w; });
  } catch { /* a failing warm-up is reported again by the real first point */ }

  for (const n of ladder) {
    if (state.stopRequested) { setStatus('stopped'); break; }
    const name = series.name(n);
    const natoms = series.natoms(n);

    saveCheckpoint({
      seriesId: series.id, method, budgetSeconds: budgetS,
      points: state.points,
      inFlight: { n, name, natoms, startedAt: Date.now() },
    });

    setStatus(`${name} — starting…`);
    const outcome = await runPoint(
      series.xyz(n), basisGBS, budgetMs, dft,
      (msg, ms) => setStatus(`${name} — ${msg} (${(ms / 1000).toFixed(1)} s)`),
      w => { state.activeWorker = w; },
    );

    const point: StressPoint = { n, name, natoms, ...outcome };
    state.points.push(point);
    renderRows();
    saveCheckpoint({
      seriesId: series.id, method, budgetSeconds: budgetS,
      points: state.points, inFlight: null,
    });

    if (outcome.status !== 'ok') {
      setStatus(`${name} failed: ${STATUS_LABEL[outcome.status]} — stopping ladder`);
      break;
    }
    setStatus(`${name} ok — ${secs(outcome.totalMs)} s`);
    // Give the engine a moment to release the previous point's ERI array.
    await new Promise(r => setTimeout(r, 200));
  }

  state.output = {
    device: gatherDeviceInfo(),
    backend: getBackendLabel(),
    seriesId: series.id,
    seriesLabel: series.label,
    method,
    basis: BASIS_NAME,
    budgetSeconds: budgetS,
    points: state.points,
    summary: summarize(state.points, budgetS),
    memoryProbe: probeRows,
    toolUrl: location.href,
  };
  renderSummary();
  state.running = false;
  state.activeWorker = null;
  if (runBtn) runBtn.disabled = false;
  if (stopBtn) stopBtn.disabled = true;
  if (!state.stopRequested) setStatus('done');
}

function renderRows(): void {
  const body = document.getElementById('cb-stress-body') as HTMLTableSectionElement | null;
  if (!body) return;
  body.innerHTML = '';
  for (const p of state.points) {
    const tr = document.createElement('tr');
    const bad = p.status !== 'ok';
    tr.innerHTML = `
      <td>${p.name}</td>
      <td class="num">${p.natoms}</td>
      <td class="num">${p.nbasis ?? '—'}</td>
      <td class="num">${humanBytes(p.eriBytes)}</td>
      <td class="num">${p.iterations ?? '—'}</td>
      <td class="num">${secs(p.totalMs)}</td>
      <td${bad ? ' class="cb-err"' : ''}>${STATUS_LABEL[p.status]}${p.detail ? `<br><span class="cb-err">${escapeHtml(p.detail)}</span>` : ''}</td>
    `;
    body.appendChild(tr);
  }
}

function escapeHtml(s: string): string {
  return s.replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]!));
}

function paperRow(): string {
  const o = state.output;
  if (!o) return '';
  const lc = o.summary.largestCompleted;
  const ff = o.summary.firstFailure;
  const cells = [
    deviceShortName(o.device),
    lc ? lc.name : '—',
    lc ? String(lc.natoms) : '—',
    lc ? String(lc.nbasis) : '—',
    lc ? lc.seconds.toFixed(1) : '—',
    ff ? ff.name : '(no failure in ladder)',
    ff ? ff.mode : '—',
  ];
  return `| ${cells.join(' | ')} |`;
}

function deviceShortName(d: DeviceInfo): string {
  const ua = d.ua;
  if (/iPhone/.test(ua)) return 'iPhone';
  if (/iPad/.test(ua)) return 'iPad';
  if (/Mac OS X|Macintosh/.test(ua)) return 'Apple-silicon laptop';
  if (/Android/.test(ua)) return 'Android';
  if (/Windows/.test(ua)) return 'Windows';
  return d.platform;
}

function renderSummary(): void {
  const box = document.getElementById('cb-stress-summary');
  const pre = document.getElementById('cb-stress-row');
  if (!box || !pre || !state.output) return;
  box.style.display = '';
  const s = state.output.summary;
  const header = '| Device | Largest completed | Atoms | Basis fns | Time (s) | First failing size | Failure mode |';
  const sep = '|---|---|---|---|---|---|---|';
  const notes: string[] = [];
  if (!s.firstFailure) notes.push('# The ladder finished without a failure — extend it to find the limit.');
  if (s.borderline) {
    const b = s.borderline;
    notes.push(
      `# BORDERLINE: ${s.firstFailure?.name} is extrapolated to need ~${b.estimatedSeconds.toFixed(0)} s ` +
      `(fitted N^${b.exponent.toFixed(1)}) against a ${b.budgetSeconds} s budget.`,
      `# The reported limit reflects the budget, not the device. Re-run with the budget ` +
      `raised to at least ${Math.ceil(b.estimatedSeconds * 1.5 / 60) * 60} s before using this row.`,
    );
  }
  const nonConverged = state.points.filter(p => p.status === 'ok' && p.converged === false);
  if (nonConverged.length) notes.push(`# SCF did not converge for: ${nonConverged.map(p => p.name).join(', ')}`);
  pre.textContent = [header, sep, paperRow(), ...notes].join('\n');
}

function copy(s: string): void {
  navigator.clipboard.writeText(s)
    .then(() => alert('Copied to clipboard'))
    .catch(() => alert('Copy failed — open browser console'));
}

function download(name: string, content: string, mime: string): void {
  const blob = new Blob([content], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function toCSV(out: StressOutput | null): string {
  if (!out) return '';
  const head = ['n', 'system', 'atoms', 'nbasis', 'eri_bytes', 'iterations', 'converged',
                'setup_plus_scf_ms', 'scf_ms', 'energy_hartree', 'status', 'detail'];
  const rows = out.points.map(p => [
    p.n, p.name, p.natoms, p.nbasis ?? '', p.eriBytes ?? '', p.iterations ?? '',
    p.converged == null ? '' : String(p.converged),
    p.totalMs?.toFixed(2) ?? '', p.scfMs?.toFixed(2) ?? '',
    p.energy?.toFixed(8) ?? '', p.status, p.detail ?? '',
  ].map(v => `"${String(v).replace(/"/g, '""')}"`).join(','));
  const meta = [
    `# device.ua,"${out.device.ua.replace(/"/g, '""')}"`,
    `# device.platform,"${out.device.platform}"`,
    `# device.deviceMemoryGB,"${out.device.deviceMemoryGB ?? ''}"`,
    `# device.hardwareConcurrency,"${out.device.hardwareConcurrency}"`,
    `# backend,"${out.backend}"`,
    `# series,"${out.seriesLabel}"`,
    `# method,"${out.method}/${out.basis}"`,
    `# budget_seconds,"${out.budgetSeconds}"`,
    `# timestamp,"${out.device.timestamp}"`,
  ];
  return meta.join('\n') + '\n' + head.join(',') + '\n' + rows.join('\n') + '\n';
}
