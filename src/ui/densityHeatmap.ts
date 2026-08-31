/** SVG-based density matrix heatmap with blue-white-red divergent colormap */

import { getThemeColors } from './theme';
import { t } from './i18n';

export function renderDensityHeatmap(
  container: HTMLElement,
  densityData: number[],
  numBasis: number,
): void {
  if (numBasis === 0 || densityData.length < numBasis * numBasis) {
    container.innerHTML = `<p style="color:var(--color-text-dim);font-size:0.85rem;">${t('density.noData')}</p>`;
    return;
  }

  const tc = getThemeColors();

  // Find max absolute value for symmetric scaling
  let maxAbs = 0;
  for (let i = 0; i < numBasis * numBasis; i++) {
    const a = Math.abs(densityData[i]);
    if (a > maxAbs) maxAbs = a;
  }
  if (maxAbs === 0) maxAbs = 1;

  const size = Math.min(280, Math.max(140, numBasis * 12));
  const ml = 28; // margin left for labels
  const mt = 28; // margin top for labels
  const svgW = size + ml + 20; // extra room for color bar
  const svgH = size + mt + 8;
  const cellW = size / numBasis;
  const cellH = size / numBasis;

  let svg = `<svg width="${svgW}" height="${svgH}" xmlns="http://www.w3.org/2000/svg" style="display:block;margin:8px auto;">`;

  // Background
  svg += `<rect width="${svgW}" height="${svgH}" fill="${tc.surface}"/>`;

  // Heatmap cells
  for (let i = 0; i < numBasis; i++) {
    for (let j = 0; j < numBasis; j++) {
      const raw = densityData[i * numBasis + j];
      const v = raw / maxAbs; // [-1, 1]
      const x = ml + j * cellW;
      const y = mt + i * cellH;
      svg += `<rect x="${x}" y="${y}" width="${Math.ceil(cellW)}" height="${Math.ceil(cellH)}" fill="${divergentColor(v)}">`;
      svg += `<title>P[${i + 1},${j + 1}] = ${raw.toFixed(6)}</title>`;
      svg += `</rect>`;
    }
  }

  // Axis labels
  const fontSize = Math.min(9, Math.max(6, cellH * 0.8));
  const step = numBasis <= 10 ? 1 : numBasis <= 30 ? 5 : 10;
  for (let i = 0; i < numBasis; i += step) {
    // Top labels (columns)
    svg += `<text x="${ml + (i + 0.5) * cellW}" y="${mt - 4}" text-anchor="middle" font-size="${fontSize}" fill="${tc.dim}">${i + 1}</text>`;
    // Left labels (rows)
    svg += `<text x="${ml - 4}" y="${mt + (i + 0.5) * cellH + fontSize / 3}" text-anchor="end" font-size="${fontSize}" fill="${tc.dim}">${i + 1}</text>`;
  }

  // Color bar
  const barX = ml + size + 6;
  const barW = 6;
  const barSteps = 40;
  const barStepH = size / barSteps;
  for (let k = 0; k < barSteps; k++) {
    const v = 1 - 2 * k / barSteps; // +1 at top, -1 at bottom
    svg += `<rect x="${barX}" y="${mt + k * barStepH}" width="${barW}" height="${Math.ceil(barStepH)}" fill="${divergentColor(v)}"/>`;
  }

  svg += '</svg>';
  container.innerHTML = `<h2 title="${t('density.titleTip')}">${t('density.title')}</h2>${svg}`;
}

/** Blue-White-Red divergent colormap. v in [-1, 1]. */
function divergentColor(v: number): string {
  v = Math.max(-1, Math.min(1, v));

  let r: number, g: number, b: number;
  if (v >= 0) {
    // White → Red
    r = 255;
    g = Math.round(255 * (1 - v));
    b = Math.round(255 * (1 - v));
  } else {
    // Blue → White
    const t = -v;
    r = Math.round(255 * (1 - t));
    g = Math.round(255 * (1 - t));
    b = 255;
  }
  return `rgb(${r},${g},${b})`;
}
