/** SVG UV-Vis absorption spectrum: Gaussian-broadened oscillator strengths vs energy */

import { getThemeColors } from './theme';
import { t } from './i18n';
import type { CISExcitedState } from '../core/cis';

const HA_TO_EV = 27.211386245988;

/** Render UV-Vis absorption spectrum as SVG.
 *  X-axis: energy (eV), Y-axis: molar absorptivity ε (arb. units).
 *  Each transition is broadened with a Gaussian peak. */
export function renderSpectrumChart(
  container: HTMLElement,
  states: CISExcitedState[],
  isTriplet: boolean,
): void {
  if (!states || states.length === 0) {
    container.innerHTML = `<p style="color:var(--color-text-dim);font-size:0.85rem;">${t('spectrum.noData')}</p>`;
    return;
  }

  const tc = getThemeColors();

  const width = 440;
  const height = 240;
  const ml = 48;
  const mr = 16;
  const mt = 28;
  const mb = 38;
  const pw = width - ml - mr;
  const ph = height - mt - mb;

  // Energy range (eV)
  const energies = states.map(s => s.energyEV);
  const eMin = Math.max(0, Math.min(...energies) - 2);
  const eMax = Math.max(...energies) + 2;

  // Gaussian broadening: σ = 0.4 eV (FWHM ≈ 0.94 eV)
  const sigma = 0.4;
  const nPts = 200;
  const de = (eMax - eMin) / nPts;

  // Build spectrum curve
  const spectrum = new Float64Array(nPts);
  for (let k = 0; k < nPts; k++) {
    const E = eMin + k * de;
    let val = 0.0;
    for (const s of states) {
      const f = s.oscillatorStrength;
      if (f > 0) {
        const diff = E - s.energyEV;
        val += f * Math.exp(-diff * diff / (2 * sigma * sigma));
      }
    }
    spectrum[k] = val;
  }

  const yMax = Math.max(...spectrum) * 1.15 || 1.0;

  const toX = (e: number) => ml + ((e - eMin) / (eMax - eMin || 1)) * pw;
  const toY = (y: number) => mt + ph - (y / yMax) * ph;

  let svg = `<svg width="${width}" height="${height}" xmlns="http://www.w3.org/2000/svg" style="display:block;width:100%;height:auto;" viewBox="0 0 ${width} ${height}">`;

  // Grid lines
  const eStep = Math.ceil((eMax - eMin) / 6);
  for (let ge = Math.ceil(eMin); ge <= eMax; ge += Math.max(1, eStep)) {
    const xx = toX(ge);
    if (xx < ml || xx > ml + pw) continue;
    svg += `<line x1="${xx}" y1="${mt}" x2="${xx}" y2="${mt + ph}" stroke="${tc.grid}" stroke-width="0.5"/>`;
    svg += `<text x="${xx}" y="${mt + ph + 14}" text-anchor="middle" font-size="9" fill="${tc.dim}">${ge}</text>`;
  }

  // Axes
  svg += `<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt + ph}" stroke="${tc.axis}" stroke-width="1"/>`;
  svg += `<line x1="${ml}" y1="${mt + ph}" x2="${ml + pw}" y2="${mt + ph}" stroke="${tc.axis}" stroke-width="1"/>`;

  // Spectrum curve (filled area)
  if (!isTriplet) {
    let pathArea = `M${ml},${mt + ph}`;
    for (let k = 0; k < nPts; k++) {
      pathArea += ` L${toX(eMin + k * de).toFixed(1)},${toY(spectrum[k]).toFixed(1)}`;
    }
    pathArea += ` L${toX(eMax).toFixed(1)},${mt + ph} Z`;
    svg += `<path d="${pathArea}" fill="${tc.accent}" fill-opacity="0.15" stroke="none"/>`;

    // Spectrum line
    let pathLine = '';
    for (let k = 0; k < nPts; k++) {
      const px = toX(eMin + k * de);
      const py = toY(spectrum[k]);
      pathLine += k === 0 ? `M${px.toFixed(1)},${py.toFixed(1)}` : ` L${px.toFixed(1)},${py.toFixed(1)}`;
    }
    svg += `<path d="${pathLine}" fill="none" stroke="${tc.accent}" stroke-width="1.5"/>`;
  }

  // Stick spectrum (vertical lines for each transition)
  for (const s of states) {
    const xx = toX(s.energyEV);
    if (xx < ml || xx > ml + pw) continue;
    const stickHeight = isTriplet ? ph * 0.5 : (s.oscillatorStrength / yMax) * ph;
    const yy = mt + ph - stickHeight;
    const color = isTriplet ? tc.error : tc.accent;
    svg += `<line x1="${xx.toFixed(1)}" y1="${(mt + ph).toFixed(1)}" x2="${xx.toFixed(1)}" y2="${yy.toFixed(1)}" stroke="${color}" stroke-width="1.5" stroke-opacity="0.7"/>`;
    svg += `<circle cx="${xx.toFixed(1)}" cy="${yy.toFixed(1)}" r="2" fill="${color}"/>`;
  }

  // Title
  const titleText = isTriplet ? t('spectrum.titleTriplet') : t('spectrum.title');
  svg += `<text x="${width / 2}" y="16" text-anchor="middle" font-size="11" fill="${tc.titleSvg}">${titleText}</text>`;

  // Axis labels
  svg += `<text x="${ml + pw / 2}" y="${height - 4}" text-anchor="middle" font-size="9" fill="${tc.dim}">${t('spectrum.xAxis')}</text>`;
  if (!isTriplet) {
    svg += `<text x="12" y="${mt + ph / 2}" text-anchor="middle" font-size="9" fill="${tc.dim}" transform="rotate(-90,12,${mt + ph / 2})">${t('spectrum.yAxis')}</text>`;
  }

  svg += '</svg>';

  // Excited states table
  let table = `<table class="result-table"><tr><th>${t('spectrum.col.state')}</th><th>${t('spectrum.col.energy')}</th><th>${t('spectrum.col.ev')}</th><th>${t('spectrum.col.f')}</th><th>${t('spectrum.col.transition')}</th></tr>`;

  for (let n = 0; n < states.length; n++) {
    const s = states[n];
    const transStr = s.dominantTransitions
      .map(tr => `${tr.i + 1}→${tr.a + 1} (${tr.coeff.toFixed(2)})`)
      .join(', ');
    table += `<tr>
      <td>${n + 1}</td>
      <td>${s.energy.toFixed(6)}</td>
      <td>${s.energyEV.toFixed(3)}</td>
      <td>${s.oscillatorStrength.toFixed(4)}</td>
      <td style="font-size:0.8rem">${transStr}</td>
    </tr>`;
  }
  table += '</table>';

  container.innerHTML = `<h2>${t('spectrum.heading')}</h2>${svg}${table}`;
}
