/** Progress tracker modal — shows calculation steps as an overlay popup */

import { t } from './i18n';

export interface StepDef {
  id: string;
  label: string;
}

type StepStatus = 'pending' | 'active' | 'done' | 'error';

interface StepState {
  el: HTMLElement;
  iconEl: HTMLElement;
  labelEl: HTMLElement;
  detailEl: HTMLElement;
  timeEl: HTMLElement;
  barEl: HTMLElement;
  status: StepStatus;
  startTime: number;
  timerHandle: number;
}

export class ProgressTracker {
  private overlay: HTMLElement;
  private card: HTMLElement;
  private stepsContainer: HTMLElement;
  private totalTimeEl: HTMLElement;
  private titleEl: HTMLElement;
  private badgeSimd: HTMLElement;
  private badgeWasm: HTMLElement;
  private badgeJS: HTMLElement;
  private cancelBtn: HTMLButtonElement;
  private t0: number;
  private totalTimer: number;
  private steps: Map<string, StepState> = new Map();
  private onCancel?: () => void;

  constructor(stepDefs: StepDef[], onCancel?: () => void) {
    this.onCancel = onCancel;
    this.t0 = performance.now();

    // Overlay
    this.overlay = document.createElement('div');
    this.overlay.className = 'pt-overlay';

    // Card
    this.card = document.createElement('div');
    this.card.className = 'pt-card';

    // Title row
    const titleRow = document.createElement('div');
    titleRow.className = 'pt-title-row';

    this.titleEl = document.createElement('div');
    this.titleEl.className = 'pt-title';
    this.titleEl.innerHTML = `<span class="pt-title-icon"></span> ${t('progress.title')}`;
    titleRow.appendChild(this.titleEl);

    // Backend badges
    const badges = document.createElement('div');
    badges.className = 'pt-badges';
    this.badgeSimd = document.createElement('span');
    this.badgeSimd.className = 'pt-badge simd inactive';
    this.badgeSimd.textContent = 'SIMD';
    this.badgeWasm = document.createElement('span');
    this.badgeWasm.className = 'pt-badge inactive';
    this.badgeWasm.textContent = 'WASM';
    this.badgeJS = document.createElement('span');
    this.badgeJS.className = 'pt-badge inactive';
    this.badgeJS.textContent = 'JS';
    badges.appendChild(this.badgeSimd);
    badges.appendChild(this.badgeWasm);
    badges.appendChild(this.badgeJS);
    titleRow.appendChild(badges);

    this.card.appendChild(titleRow);

    // Steps
    this.stepsContainer = document.createElement('div');
    this.stepsContainer.className = 'pt-steps';

    for (const def of stepDefs) {
      const el = document.createElement('div');
      el.className = 'pt-step pending';

      const iconEl = document.createElement('span');
      iconEl.className = 'pt-icon';

      const labelEl = document.createElement('span');
      labelEl.className = 'pt-label';
      labelEl.textContent = def.label;

      const right = document.createElement('span');
      right.className = 'pt-right';

      const detailEl = document.createElement('span');
      detailEl.className = 'pt-detail';

      const timeEl = document.createElement('span');
      timeEl.className = 'pt-time';

      right.appendChild(detailEl);
      right.appendChild(timeEl);

      // Progress bar (for steps like integrals)
      const barEl = document.createElement('div');
      barEl.className = 'pt-bar';

      el.appendChild(iconEl);
      el.appendChild(labelEl);
      el.appendChild(right);
      el.appendChild(barEl);
      this.stepsContainer.appendChild(el);

      this.steps.set(def.id, {
        el, iconEl, labelEl, detailEl, timeEl, barEl,
        status: 'pending',
        startTime: 0,
        timerHandle: 0,
      });
    }

    this.card.appendChild(this.stepsContainer);

    // Footer: cancel button + total time
    const footer = document.createElement('div');
    footer.className = 'pt-footer';
    this.cancelBtn = document.createElement('button');
    this.cancelBtn.className = 'pt-cancel';
    this.cancelBtn.textContent = t('progress.cancel');
    this.cancelBtn.addEventListener('click', () => this.onCancel?.());
    footer.appendChild(this.cancelBtn);
    this.totalTimeEl = document.createElement('span');
    this.totalTimeEl.className = 'pt-total';
    footer.appendChild(this.totalTimeEl);
    this.card.appendChild(footer);

    this.overlay.appendChild(this.card);
    document.body.appendChild(this.overlay);

    // Animate in
    requestAnimationFrame(() => {
      this.overlay.classList.add('visible');
    });

    // Total elapsed timer
    this.totalTimer = window.setInterval(() => {
      const t = (performance.now() - this.t0) / 1000;
      this.totalTimeEl.textContent = formatTime(t);
    }, 100);
  }

  /** Set which backend is currently active */
  setBackend(backend: 'wasm-simd' | 'wasm' | 'js' | 'none'): void {
    this.badgeSimd.className = `pt-badge simd ${backend === 'wasm-simd' ? 'active' : 'inactive'}`;
    this.badgeWasm.className = `pt-badge ${backend === 'wasm' ? 'active' : 'inactive'}`;
    this.badgeJS.className = `pt-badge ${backend === 'js' ? 'active' : 'inactive'}`;
  }

  startStep(id: string, detail?: string): void {
    const s = this.steps.get(id);
    if (!s) return;
    s.status = 'active';
    s.el.className = 'pt-step active';
    s.startTime = performance.now();
    if (detail) s.detailEl.textContent = detail;

    s.timerHandle = window.setInterval(() => {
      const elapsed = (performance.now() - s.startTime) / 1000;
      s.timeEl.textContent = formatTime(elapsed);
    }, 100);
  }

  updateStep(id: string, detail: string): void {
    const s = this.steps.get(id);
    if (!s) return;
    s.detailEl.textContent = detail;
  }

  completeStep(id: string, detail?: string): void {
    const s = this.steps.get(id);
    if (!s) return;
    clearInterval(s.timerHandle);
    s.status = 'done';
    s.el.className = 'pt-step done';
    if (detail) s.detailEl.textContent = detail;
    const elapsed = (performance.now() - s.startTime) / 1000;
    s.timeEl.textContent = formatTime(elapsed);
  }

  failStep(id: string, detail?: string): void {
    const s = this.steps.get(id);
    if (!s) return;
    clearInterval(s.timerHandle);
    s.status = 'error';
    s.el.className = 'pt-step error';
    if (detail) s.detailEl.textContent = detail;
    if (s.startTime > 0) {
      const elapsed = (performance.now() - s.startTime) / 1000;
      s.timeEl.textContent = formatTime(elapsed);
    }
  }

  failRemaining(): void {
    for (const s of this.steps.values()) {
      if (s.status === 'active') {
        clearInterval(s.timerHandle);
        s.status = 'error';
        s.el.className = 'pt-step error';
      }
    }
  }

  /** Mark as cancelled and close after delay */
  cancel(): void {
    this.titleEl.innerHTML = `<span class="pt-title-icon error"></span> ${t('progress.cancelled')}`;
    this.cancelBtn.disabled = true;
    for (const s of this.steps.values()) {
      if (s.status === 'active') {
        clearInterval(s.timerHandle);
        s.status = 'error';
        s.el.className = 'pt-step error';
        if (s.startTime > 0) {
          const elapsed = (performance.now() - s.startTime) / 1000;
          s.timeEl.textContent = formatTime(elapsed);
        }
      }
    }
    setTimeout(() => this.close(), 1200);
  }

  /** Close modal with fade-out */
  close(): void {
    clearInterval(this.totalTimer);
    for (const s of this.steps.values()) clearInterval(s.timerHandle);
    this.overlay.classList.remove('visible');
    this.overlay.addEventListener('transitionend', () => this.overlay.remove(), { once: true });
    // Fallback removal
    setTimeout(() => { if (this.overlay.parentNode) this.overlay.remove(); }, 500);
  }
}

function formatTime(seconds: number): string {
  if (seconds < 0.01) return '<0.01s';
  if (seconds < 10) return seconds.toFixed(2) + 's';
  if (seconds < 60) return seconds.toFixed(1) + 's';
  const m = Math.floor(seconds / 60);
  const s = seconds % 60;
  return `${m}m ${s.toFixed(0)}s`;
}
