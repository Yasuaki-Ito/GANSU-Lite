/** SVG orbital energy diagram — zoomable, with label de-overlap and degeneracy grouping */

import type { FloatArray } from '../linalg/matrix';
import { getThemeColors } from './theme';
import { t } from './i18n';

const MIN_LABEL_GAP = 12; // minimum px between label centers
const DEGEN_THRESHOLD = 1e-4; // Hartree — orbitals closer than this are grouped as degenerate
const ORB_GAP = 6; // px gap between degenerate orbital lines

interface OrbGroup {
  indices: number[];
  energy: number; // representative energy (average)
}

/** Group consecutive orbitals with nearly identical energies */
function groupByDegeneracy(energies: FloatArray, n: number): OrbGroup[] {
  const groups: OrbGroup[] = [];
  let i = 0;
  while (i < n) {
    const indices = [i];
    const e0 = energies[i];
    while (i + 1 < n && Math.abs(energies[i + 1] - e0) < DEGEN_THRESHOLD) {
      i++;
      indices.push(i);
    }
    const avg = indices.reduce((s, idx) => s + energies[idx], 0) / indices.length;
    groups.push({ indices, energy: avg });
    i++;
  }
  return groups;
}

/** Compute the right edge x of a group drawn at xBase with full-width lines */
function groupRightEdge(xBase: number, k: number, lineLen: number): number {
  return xBase + k * lineLen + (k - 1) * ORB_GAP;
}

/** Push label Y positions apart so no two are closer than minGap */
function deOverlap(positions: number[], minGap: number): number[] {
  const out = positions.slice();
  const idx = out.map((_, i) => i);
  idx.sort((a, b) => out[a] - out[b]);

  for (let k = 1; k < idx.length; k++) {
    const prev = idx[k - 1];
    const curr = idx[k];
    if (out[curr] - out[prev] < minGap) {
      out[curr] = out[prev] + minGap;
    }
  }
  return out;
}

/** Create a zoomable/scrollable wrapper with Ctrl+wheel zoom and +/- buttons */
function makeZoomableContainer(
  container: HTMLElement,
  buildSvg: (zoom: number) => string,
  baseHeight: number,
): void {
  let zoom = 1;
  const maxZoom = 8;
  const minZoom = 0.5;
  const scrollH = Math.min(baseHeight, 520);

  // Outer: holds toolbar + scrollable area
  const outer = document.createElement('div');
  outer.style.position = 'relative';

  // Zoom toolbar
  const toolbar = document.createElement('div');
  toolbar.style.cssText = 'display:flex;align-items:center;justify-content:center;gap:6px;padding:2px 0 4px;font-size:11px;color:var(--color-text-dim,#888)';
  const btnCss = 'width:22px;height:22px;border:1px solid var(--color-border,#ccc);border-radius:4px;background:var(--color-surface,#fff);color:var(--color-text,#333);cursor:pointer;font-size:14px;line-height:1;display:flex;align-items:center;justify-content:center';
  const zoomOutBtn = document.createElement('button');
  zoomOutBtn.style.cssText = btnCss;
  zoomOutBtn.textContent = '\u2212'; // minus sign
  zoomOutBtn.title = t('orb.zoomOut');
  const zoomInBtn = document.createElement('button');
  zoomInBtn.style.cssText = btnCss;
  zoomInBtn.textContent = '+';
  zoomInBtn.title = t('orb.zoomIn');
  const zoomLabel = document.createElement('span');
  zoomLabel.style.cssText = 'min-width:40px;text-align:center;font-size:10px';

  toolbar.appendChild(zoomOutBtn);
  toolbar.appendChild(zoomLabel);
  toolbar.appendChild(zoomInBtn);
  outer.appendChild(toolbar);

  const wrapper = document.createElement('div');
  wrapper.style.overflowY = 'auto';
  wrapper.style.maxHeight = `${scrollH}px`;
  wrapper.style.position = 'relative';
  outer.appendChild(wrapper);

  function applyZoom(newZoom: number) {
    newZoom = Math.max(minZoom, Math.min(maxZoom, newZoom));
    if (newZoom === zoom) return;
    const scrollRatio = wrapper.scrollHeight > scrollH
      ? wrapper.scrollTop / (wrapper.scrollHeight - scrollH)
      : 0;
    zoom = newZoom;
    render();
    const newScrollMax = wrapper.scrollHeight - scrollH;
    if (newScrollMax > 0) {
      wrapper.scrollTop = scrollRatio * newScrollMax;
    }
  }

  function render() {
    wrapper.innerHTML = buildSvg(zoom);
    zoomLabel.textContent = `${Math.round(zoom * 100)}%`;
    zoomOutBtn.disabled = zoom <= minZoom;
    zoomInBtn.disabled = zoom >= maxZoom;
  }

  // Ctrl+wheel zoom (regular scroll passes through for page scrolling)
  wrapper.addEventListener('wheel', (e: WheelEvent) => {
    if (!e.ctrlKey && !e.metaKey) return; // let normal scroll through
    e.preventDefault();
    const factor = e.deltaY < 0 ? 1.2 : 1 / 1.2;
    applyZoom(zoom * factor);
  }, { passive: false });

  // +/- buttons
  zoomInBtn.addEventListener('click', () => applyZoom(zoom * 1.3));
  zoomOutBtn.addEventListener('click', () => applyZoom(zoom / 1.3));

  render();
  container.appendChild(outer);
}

/** Draw up+down arrow pair (doubly occupied) centered at cx, above y */
function drawDoubleArrows(cx: number, y: number, color: string): string {
  const hw = 3;
  const off = hw + 0.5;
  let svg = `<polygon points="${cx - off - hw},${y - 2} ${cx - off},${y - 10} ${cx - off + hw},${y - 2}" fill="${color}"/>`;
  svg += `<polygon points="${cx + off - hw},${y - 10} ${cx + off},${y - 2} ${cx + off + hw},${y - 10}" fill="${color}"/>`;
  return svg;
}

/** Draw single up arrow centered at cx, above y */
function drawUpArrow(cx: number, y: number, color: string): string {
  return `<polygon points="${cx - 3},${y - 2} ${cx},${y - 10} ${cx + 3},${y - 2}" fill="${color}"/>`;
}

/** Draw single down arrow centered at cx, above y */
function drawDownArrow(cx: number, y: number, color: string): string {
  return `<polygon points="${cx - 3},${y - 10} ${cx},${y - 2} ${cx + 3},${y - 10}" fill="${color}"/>`;
}

/** Draw HOMO-LUMO gap annotation at a given X position */
function drawGapAnnotation(
  homoE: number, lumoE: number, toY: (e: number) => number,
  gapX: number, tc: ReturnType<typeof getThemeColors>,
): string {
  const gapEv = (lumoE - homoE) * 27.2114;
  if (Math.abs(gapEv) < 0.01) return '';
  const yHOMO = toY(homoE);
  const yLUMO = toY(lumoE);

  let svg = `<line x1="${gapX}" y1="${yHOMO}" x2="${gapX}" y2="${yLUMO}" stroke="${tc.gap}" stroke-width="1.5" stroke-dasharray="3,2"/>`;
  svg += `<polygon points="${gapX - 3},${yHOMO - 4} ${gapX},${yHOMO} ${gapX + 3},${yHOMO - 4}" fill="${tc.gap}"/>`;
  svg += `<polygon points="${gapX - 3},${yLUMO + 4} ${gapX},${yLUMO} ${gapX + 3},${yLUMO + 4}" fill="${tc.gap}"/>`;
  const gapMidY = (yHOMO + yLUMO) / 2;
  svg += `<text x="${gapX + 6}" y="${gapMidY + 4}" font-size="10" fill="${tc.gap}">${gapEv.toFixed(2)} eV</text>`;
  return svg;
}

/** RHF orbital diagram (single column, up+down arrows for occupied, degeneracy grouping) */
export function renderOrbitalDiagram(
  container: HTMLElement,
  orbitalEnergies: FloatArray,
  numOccupied: number,
): void {
  const n = orbitalEnergies.length;
  if (n === 0) return;

  const groups = groupByDegeneracy(orbitalEnergies, n);

  const marginTop = 28;
  const marginBottom = 24;
  const marginLeft = 70;
  const lineLen = 40;
  const maxK = Math.max(...groups.map(g => g.indices.length));
  const maxRight = groupRightEdge(marginLeft, maxK, lineLen);
  const width = Math.max(300, maxRight + 100);

  const eMin = orbitalEnergies[0];
  const eMax = orbitalEnergies[n - 1];
  const range = eMax - eMin || 1;
  const padded = range * 0.1;
  const yMin = eMin - padded;
  const yMax = eMax + padded;
  const basePlotH = Math.max(180, Math.min(480, groups.length * 40));

  function buildSvg(zoom: number): string {
    const tc = getThemeColors();
    const plotH = basePlotH * zoom;
    const height = plotH + marginTop + marginBottom;
    const toY = (e: number) => marginTop + plotH - ((e - yMin) / (yMax - yMin)) * plotH;

    const rawGroupY = groups.map(g => toY(g.energy));
    const groupLabelY = deOverlap(rawGroupY, MIN_LABEL_GAP);

    let svg = `<svg width="${width}" height="${height}" xmlns="http://www.w3.org/2000/svg" style="display:block;margin:0 auto;">`;

    // Axis line
    svg += `<line x1="${marginLeft - 8}" y1="${marginTop}" x2="${marginLeft - 8}" y2="${marginTop + plotH}" stroke="${tc.axis}" stroke-width="1"/>`;

    for (let gi = 0; gi < groups.length; gi++) {
      const group = groups[gi];
      const y = toY(group.energy);
      const ly = groupLabelY[gi];
      const k = group.indices.length;
      const gRight = groupRightEdge(marginLeft, k, lineLen);

      const groupHasHOMO = numOccupied > 0 && group.indices.includes(numOccupied - 1);
      const groupHasLUMO = numOccupied < n && group.indices.includes(numOccupied);

      for (let g = 0; g < k; g++) {
        const idx = group.indices[g];
        const x0 = marginLeft + g * (lineLen + ORB_GAP);
        const isOcc = idx < numOccupied;
        const isSpecial = (idx === numOccupied - 1) || (idx === numOccupied);
        const color = isOcc ? tc.occupied : tc.virtual;
        const strokeW = isSpecial ? 2.5 : 1.5;

        svg += `<line x1="${x0}" y1="${y}" x2="${x0 + lineLen}" y2="${y}" stroke="${color}" stroke-width="${strokeW}"/>`;

        if (isOcc) {
          svg += drawDoubleArrows(x0 + lineLen / 2, y, color);
        }
      }

      // Leader line
      if (Math.abs(ly - y) > 3) {
        svg += `<line x1="${marginLeft - 10}" y1="${y}" x2="${marginLeft - 14}" y2="${ly}" stroke="${tc.leader}" stroke-width="0.5"/>`;
      }

      // Energy label
      const degLabel = k > 1 ? ` (\u00d7${k})` : '';
      svg += `<text x="${marginLeft - 16}" y="${ly + 4}" text-anchor="end" font-size="10" fill="${tc.label}">${group.energy.toFixed(3)}${degLabel}</text>`;

      // HOMO / LUMO labels (after group's right edge)
      if (groupHasHOMO) {
        svg += `<text x="${gRight + 6}" y="${y + 4}" font-size="11" fill="${tc.occupied}" font-weight="bold">HOMO</text>`;
      }
      if (groupHasLUMO && !groupHasHOMO) {
        svg += `<text x="${gRight + 6}" y="${y + 4}" font-size="11" fill="${tc.dim}" font-weight="bold">LUMO</text>`;
      }
    }

    // HOMO-LUMO gap annotation
    if (numOccupied > 0 && numOccupied < n) {
      svg += drawGapAnnotation(
        orbitalEnergies[numOccupied - 1], orbitalEnergies[numOccupied],
        toY, maxRight + 50, tc,
      );
    }

    svg += `<text x="${width / 2}" y="14" text-anchor="middle" font-size="10" fill="${tc.hint}">${t('orb.zoomHint')}</text>`;

    svg += '</svg>';
    return svg;
  }

  container.innerHTML = '';
  const baseH = basePlotH + marginTop + marginBottom;
  makeZoomableContainer(container, buildSvg, baseH);
}

/** ROHF orbital diagram — single column with doubly/singly occupied distinction, degeneracy grouping */
export function renderROHFOrbitalDiagram(
  container: HTMLElement,
  orbitalEnergies: FloatArray,
  numOccAlpha: number,
  numOccBeta: number,
): void {
  const n = orbitalEnergies.length;
  if (n === 0) return;

  const groups = groupByDegeneracy(orbitalEnergies, n);

  const marginTop = 28;
  const marginBottom = 24;
  const marginLeft = 70;
  const lineLen = 40;
  const maxK = Math.max(...groups.map(g => g.indices.length));
  const maxRight = groupRightEdge(marginLeft, maxK, lineLen);
  const width = Math.max(300, maxRight + 100);

  const eMin = orbitalEnergies[0];
  const eMax = orbitalEnergies[n - 1];
  const range = eMax - eMin || 1;
  const padded = range * 0.1;
  const yMin = eMin - padded;
  const yMax = eMax + padded;
  const basePlotH = Math.max(180, Math.min(480, groups.length * 40));

  function buildSvg(zoom: number): string {
    const tc = getThemeColors();
    const plotH = basePlotH * zoom;
    const height = plotH + marginTop + marginBottom;
    const toY = (e: number) => marginTop + plotH - ((e - yMin) / (yMax - yMin)) * plotH;

    const rawGroupY = groups.map(g => toY(g.energy));
    const groupLabelY = deOverlap(rawGroupY, MIN_LABEL_GAP);

    let svg = `<svg width="${width}" height="${height}" xmlns="http://www.w3.org/2000/svg" style="display:block;margin:0 auto;">`;

    // Axis line
    svg += `<line x1="${marginLeft - 8}" y1="${marginTop}" x2="${marginLeft - 8}" y2="${marginTop + plotH}" stroke="${tc.axis}" stroke-width="1"/>`;

    for (let gi = 0; gi < groups.length; gi++) {
      const group = groups[gi];
      const y = toY(group.energy);
      const ly = groupLabelY[gi];
      const k = group.indices.length;
      const gRight = groupRightEdge(marginLeft, k, lineLen);

      const groupHasSOMO = group.indices.some(idx => idx >= numOccBeta && idx < numOccAlpha);
      const groupHasHOMO = numOccAlpha > 0 && group.indices.includes(numOccAlpha - 1);
      const groupHasLUMO = numOccAlpha < n && group.indices.includes(numOccAlpha);

      for (let g = 0; g < k; g++) {
        const idx = group.indices[g];
        const x0 = marginLeft + g * (lineLen + ORB_GAP);
        const isDoublyOcc = idx < numOccBeta;
        const isSinglyOcc = idx >= numOccBeta && idx < numOccAlpha;
        const isOcc = idx < numOccAlpha;
        const isSpecial = (idx === numOccAlpha - 1) || (idx === numOccAlpha) || isSinglyOcc;
        const color = isOcc ? tc.occupied : tc.virtual;
        const strokeW = isSpecial ? 2.5 : 1.5;

        svg += `<line x1="${x0}" y1="${y}" x2="${x0 + lineLen}" y2="${y}" stroke="${color}" stroke-width="${strokeW}"/>`;

        if (isDoublyOcc) {
          svg += drawDoubleArrows(x0 + lineLen / 2, y, color);
        } else if (isSinglyOcc) {
          svg += drawUpArrow(x0 + lineLen / 2, y, color);
        }
      }

      // Leader line
      if (Math.abs(ly - y) > 3) {
        svg += `<line x1="${marginLeft - 10}" y1="${y}" x2="${marginLeft - 14}" y2="${ly}" stroke="${tc.leader}" stroke-width="0.5"/>`;
      }

      // Energy label
      const degLabel = k > 1 ? ` (\u00d7${k})` : '';
      svg += `<text x="${marginLeft - 16}" y="${ly + 4}" text-anchor="end" font-size="10" fill="${tc.label}">${group.energy.toFixed(3)}${degLabel}</text>`;

      // Orbital labels
      if (groupHasSOMO && !groupHasHOMO) {
        svg += `<text x="${gRight + 6}" y="${y + 4}" font-size="10" fill="${tc.alpha}" font-weight="bold">SOMO</text>`;
      }
      if (groupHasHOMO) {
        const lbl = groupHasSOMO ? 'SOMO' : 'HOMO';
        svg += `<text x="${gRight + 6}" y="${y + 4}" font-size="11" fill="${tc.occupied}" font-weight="bold">${lbl}</text>`;
      }
      if (groupHasLUMO && !groupHasHOMO) {
        svg += `<text x="${gRight + 6}" y="${y + 4}" font-size="11" fill="${tc.dim}" font-weight="bold">LUMO</text>`;
      }
    }

    // HOMO-LUMO gap annotation
    if (numOccAlpha > 0 && numOccAlpha < n) {
      svg += drawGapAnnotation(
        orbitalEnergies[numOccAlpha - 1], orbitalEnergies[numOccAlpha],
        toY, maxRight + 50, tc,
      );
    }

    svg += `<text x="${width / 2}" y="14" text-anchor="middle" font-size="10" fill="${tc.hint}">${t('orb.zoomHint')}</text>`;

    svg += '</svg>';
    return svg;
  }

  container.innerHTML = '';
  const baseH = basePlotH + marginTop + marginBottom;
  makeZoomableContainer(container, buildSvg, baseH);
}

/** UHF orbital diagram — alpha (left) and beta (right) side by side, with degeneracy grouping */
export function renderUHFOrbitalDiagram(
  container: HTMLElement,
  alphaEnergies: FloatArray,
  numOccAlpha: number,
  betaEnergies: FloatArray,
  numOccBeta: number,
): void {
  const nA = alphaEnergies.length;
  const nB = betaEnergies.length;
  if (nA === 0 && nB === 0) return;

  const groupsA = groupByDegeneracy(alphaEnergies, nA);
  const groupsB = groupByDegeneracy(betaEnergies, nB);

  const marginTop = 28;
  const marginBottom = 24;
  const marginLeft = 60;
  const lineLen = 36;
  const colSep = 50; // space between alpha right edge and beta left

  const maxKa = Math.max(1, ...groupsA.map(g => g.indices.length));
  const maxKb = Math.max(1, ...groupsB.map(g => g.indices.length));
  const alphaMaxRight = groupRightEdge(marginLeft, maxKa, lineLen);
  const betaX = alphaMaxRight + colSep;
  const betaMaxRight = groupRightEdge(betaX, maxKb, lineLen);
  const width = Math.max(420, betaMaxRight + 80);

  // Unified energy scale across both spins
  const allE: number[] = [];
  for (let i = 0; i < nA; i++) allE.push(alphaEnergies[i]);
  for (let i = 0; i < nB; i++) allE.push(betaEnergies[i]);
  allE.sort((a, b) => a - b);
  const eMin = allE[0];
  const eMax = allE[allE.length - 1];
  const range = eMax - eMin || 1;
  const padded = range * 0.1;
  const yMin = eMin - padded;
  const yMax = eMax + padded;
  const nMaxGroups = Math.max(groupsA.length, groupsB.length);
  const basePlotH = Math.max(180, Math.min(480, nMaxGroups * 40));

  function buildSvg(zoom: number): string {
    const tc = getThemeColors();
    const plotH = basePlotH * zoom;
    const height = plotH + marginTop + marginBottom;
    const toY = (e: number) => marginTop + plotH - ((e - yMin) / (yMax - yMin)) * plotH;

    let svg = `<svg width="${width}" height="${height}" xmlns="http://www.w3.org/2000/svg" style="display:block;margin:0 auto;">`;

    // Title
    svg += `<text x="${width / 2}" y="14" text-anchor="middle" font-size="10" fill="${tc.hint}">${t('orb.zoomHint')}</text>`;

    function drawColumn(
      groups: OrbGroup[], numOrb: number, numOcc: number,
      xBase: number, maxGroupRight: number,
      color: string, arrowUp: boolean, label: string,
    ) {
      // Column header (centered over widest group area)
      const headerX = (xBase + maxGroupRight) / 2;
      svg += `<text x="${headerX}" y="${marginTop - 4}" text-anchor="middle" font-size="10" fill="${color}" font-weight="bold">${label}</text>`;

      // Axis line
      svg += `<line x1="${xBase - 8}" y1="${marginTop}" x2="${xBase - 8}" y2="${marginTop + plotH}" stroke="${tc.axis}" stroke-width="1"/>`;

      const rawGroupY = groups.map(g => toY(g.energy));
      const groupLabelY = deOverlap(rawGroupY, MIN_LABEL_GAP);

      for (let gi = 0; gi < groups.length; gi++) {
        const group = groups[gi];
        const y = toY(group.energy);
        const ly = groupLabelY[gi];
        const k = group.indices.length;
        const gRight = groupRightEdge(xBase, k, lineLen);

        const groupHasHOMO = numOcc > 0 && group.indices.includes(numOcc - 1);
        const groupHasLUMO = numOcc < numOrb && group.indices.includes(numOcc);

        for (let g = 0; g < k; g++) {
          const idx = group.indices[g];
          const x0 = xBase + g * (lineLen + ORB_GAP);
          const isOcc = idx < numOcc;
          const isSpecial = (idx === numOcc - 1) || (idx === numOcc);
          const lineColor = isOcc ? color : tc.virtual;
          const strokeW = isSpecial ? 2.5 : 1.5;

          svg += `<line x1="${x0}" y1="${y}" x2="${x0 + lineLen}" y2="${y}" stroke="${lineColor}" stroke-width="${strokeW}"/>`;

          if (isOcc) {
            const cx = x0 + lineLen / 2;
            if (arrowUp) {
              svg += drawUpArrow(cx, y, lineColor);
            } else {
              svg += drawDownArrow(cx, y, lineColor);
            }
          }
        }

        // Energy label with leader line
        if (arrowUp) {
          // Alpha: labels on left
          if (Math.abs(ly - y) > 3) {
            svg += `<line x1="${xBase - 10}" y1="${y}" x2="${xBase - 14}" y2="${ly}" stroke="${tc.leader}" stroke-width="0.5"/>`;
          }
          const degLabel = k > 1 ? ` (\u00d7${k})` : '';
          svg += `<text x="${xBase - 16}" y="${ly + 4}" text-anchor="end" font-size="9" fill="${tc.label}">${group.energy.toFixed(3)}${degLabel}</text>`;
        } else {
          // Beta: labels on right of group
          if (Math.abs(ly - y) > 3) {
            svg += `<line x1="${gRight + 2}" y1="${y}" x2="${gRight + 6}" y2="${ly}" stroke="${tc.leader}" stroke-width="0.5"/>`;
          }
          const degLabel = k > 1 ? ` (\u00d7${k})` : '';
          svg += `<text x="${gRight + 8}" y="${ly + 4}" text-anchor="start" font-size="9" fill="${tc.label}">${group.energy.toFixed(3)}${degLabel}</text>`;
        }

        // HOMO / LUMO label
        if (groupHasHOMO) {
          const lx = arrowUp ? gRight + 4 : xBase - 12;
          const anchor = arrowUp ? 'start' : 'end';
          svg += `<text x="${lx}" y="${y + 4}" text-anchor="${anchor}" font-size="9" fill="${color}" font-weight="bold">HOMO</text>`;
        }
        if (groupHasLUMO && !groupHasHOMO) {
          const lx = arrowUp ? gRight + 4 : xBase - 12;
          const anchor = arrowUp ? 'start' : 'end';
          svg += `<text x="${lx}" y="${y + 4}" text-anchor="${anchor}" font-size="9" fill="${tc.dim}" font-weight="bold">LUMO</text>`;
        }
      }
    }

    drawColumn(groupsA, nA, numOccAlpha, marginLeft, alphaMaxRight, tc.alpha, true, '\u03B1');
    drawColumn(groupsB, nB, numOccBeta, betaX, betaMaxRight, tc.beta, false, '\u03B2');

    svg += '</svg>';
    return svg;
  }

  container.innerHTML = '';
  const baseH = basePlotH + marginTop + marginBottom;
  makeZoomableContainer(container, buildSvg, baseH);
}
