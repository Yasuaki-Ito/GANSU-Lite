/** SVG convergence graph: log10(|ΔE|) vs iteration */

import { getThemeColors } from './theme';
import { t } from './i18n';

export function renderConvergenceGraph(
  container: HTMLElement,
  iterations: Array<{ iter: number; deltaE: number }>,
  threshold: number,
): void {
  const title = t('graph.title');
  if (iterations.length < 2) {
    container.innerHTML = `<p style="color:var(--color-text-dim);font-size:0.85rem;">${t('graph.noIter')}</p>`;
    return;
  }

  // Filter out iter=0 (deltaE = energy - 0, not a real change) and zero deltaE
  const pts = iterations
    .filter(d => d.iter >= 1 && d.deltaE !== 0)
    .map(d => ({ x: d.iter, y: Math.log10(Math.abs(d.deltaE)) }));

  if (pts.length < 2) {
    container.innerHTML = `<p style="color:var(--color-text-dim);font-size:0.85rem;">${t('graph.noData')}</p>`;
    return;
  }

  const tc = getThemeColors();

  const width = 320;
  const height = 200;
  const ml = 52; // margin left
  const mr = 12;
  const mt = 24;
  const mb = 32;
  const pw = width - ml - mr;   // plot width
  const ph = height - mt - mb;  // plot height

  const xMin = pts[0].x;
  const xMax = pts[pts.length - 1].x;
  const yMin = Math.min(...pts.map(p => p.y), Math.log10(threshold)) - 1;
  const yMax = Math.max(pts[0].y, 0) + 0.5;

  const toX = (x: number) => ml + ((x - xMin) / (xMax - xMin || 1)) * pw;
  const toY = (y: number) => mt + ph - ((y - yMin) / (yMax - yMin || 1)) * ph;

  let svg = `<svg width="${width}" height="${height}" xmlns="http://www.w3.org/2000/svg" style="display:block;width:100%;height:auto;" viewBox="0 0 ${width} ${height}">`;

  // Background grid lines (Y axis)
  const yStep = 2; // every 2 orders of magnitude
  const yGridStart = Math.ceil(yMin / yStep) * yStep;
  for (let gy = yGridStart; gy <= yMax; gy += yStep) {
    const yy = toY(gy);
    if (yy < mt || yy > mt + ph) continue;
    svg += `<line x1="${ml}" y1="${yy}" x2="${ml + pw}" y2="${yy}" stroke="${tc.grid}" stroke-width="0.5"/>`;
    svg += `<text x="${ml - 4}" y="${yy + 3}" text-anchor="end" font-size="9" fill="${tc.dim}">1e${gy}</text>`;
  }

  // Axes
  svg += `<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt + ph}" stroke="${tc.axis}" stroke-width="1"/>`;
  svg += `<line x1="${ml}" y1="${mt + ph}" x2="${ml + pw}" y2="${mt + ph}" stroke="${tc.axis}" stroke-width="1"/>`;

  // X axis labels
  const xStep = Math.max(1, Math.round((xMax - xMin) / 5));
  for (let gx = xMin; gx <= xMax; gx += xStep) {
    const xx = toX(gx);
    svg += `<text x="${xx}" y="${mt + ph + 14}" text-anchor="middle" font-size="9" fill="${tc.dim}">${gx}</text>`;
  }

  // Convergence threshold line
  const threshY = toY(Math.log10(threshold));
  if (threshY >= mt && threshY <= mt + ph) {
    svg += `<line x1="${ml}" y1="${threshY}" x2="${ml + pw}" y2="${threshY}" stroke="${tc.error}" stroke-width="1" stroke-dasharray="4,3"/>`;
    svg += `<text x="${ml + pw + 2}" y="${threshY + 3}" font-size="8" fill="${tc.error}">${threshold.toExponential(0)}</text>`;
  }

  // Data line
  let path = '';
  for (let i = 0; i < pts.length; i++) {
    const px = toX(pts[i].x);
    const py = toY(pts[i].y);
    path += i === 0 ? `M${px},${py}` : ` L${px},${py}`;
  }
  svg += `<path d="${path}" fill="none" stroke="${tc.accent}" stroke-width="1.5"/>`;

  // Data points
  for (const p of pts) {
    svg += `<circle cx="${toX(p.x)}" cy="${toY(p.y)}" r="2.5" fill="${tc.accent}"/>`;
  }

  // Title
  svg += `<text x="${width / 2}" y="14" text-anchor="middle" font-size="11" fill="${tc.titleSvg}">${title}</text>`;

  // Axis labels
  svg += `<text x="${ml + pw / 2}" y="${height - 2}" text-anchor="middle" font-size="9" fill="${tc.dim}">${t('graph.xAxis')}</text>`;
  svg += `<text x="12" y="${mt + ph / 2}" text-anchor="middle" font-size="9" fill="${tc.dim}" transform="rotate(-90,12,${mt + ph / 2})">${t('graph.yAxis')}</text>`;

  svg += '</svg>';
  container.innerHTML = `<h2 title="${t('graph.titleTip')}">${title}</h2>${svg}`;
}
