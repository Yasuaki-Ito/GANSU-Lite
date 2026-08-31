/** Interactive 3D molecule viewer — Three.js WebGL */

import * as THREE from 'three';
import type { Atom } from '../core/types';
import { BOHR_TO_ANGSTROM, ELEMENT_NAME_TO_ATOMIC_NUMBER } from '../core/constants';
import { isDark } from './theme';
import { t } from './i18n';

// CPK element colors
const CPK_COLORS: Record<number, number> = {
  1: 0xFFFFFF, 2: 0xD9FFFF, 3: 0xCC80FF, 4: 0xC2FF00, 5: 0xFFB5B5,
  6: 0x909090, 7: 0x3050F8, 8: 0xFF0D0D, 9: 0x90E050, 10: 0xB3E3F5,
  11: 0xAB5CF2, 12: 0x8AFF00, 13: 0xBFA6A6, 14: 0xF0C8A0, 15: 0xFF8000,
  16: 0xFFFF30, 17: 0x1FF01F, 18: 0x80D1E3, 19: 0x8F40D4, 20: 0x3DFF00,
  26: 0xE06633, 29: 0xC88033, 30: 0x7D80B0, 35: 0xA62929, 53: 0x940094,
};

// Covalent radii in Angstrom
const COVALENT_RADII: Record<number, number> = {
  1: 0.31, 2: 0.28, 3: 1.28, 4: 0.96, 5: 0.84, 6: 0.76, 7: 0.71, 8: 0.66,
  9: 0.57, 10: 0.58, 11: 1.66, 12: 1.41, 13: 1.21, 14: 1.11, 15: 1.07, 16: 1.05,
  17: 1.02, 18: 1.06, 19: 2.03, 20: 1.76, 26: 1.32, 29: 1.32, 30: 1.22,
  35: 1.20, 53: 1.39,
};

// Display radii for atoms (Angstrom)
const ATOM_RADII: Record<number, number> = {
  1: 0.25, 6: 0.40, 7: 0.38, 8: 0.36, 9: 0.33, 15: 0.45, 16: 0.45, 17: 0.43,
};

function getColor(z: number): number { return CPK_COLORS[z] ?? 0xFF69B4; }
function getCovalentRadius(z: number): number { return COVALENT_RADII[z] ?? 1.5; }
function getAtomRadius(z: number): number { return ATOM_RADII[z] ?? 0.40; }

interface Vec3 { x: number; y: number; z: number; }
interface Bond { i: number; j: number; }

/** Options for customising the 3D viewer (e.g. charge colouring). */
export interface ViewerOptions {
  /** Per-atom 0xRRGGBB colour override (replaces CPK). */
  atomColors?: number[];
  /** Per-atom text labels shown below each atom (e.g. charge values). */
  atomLabels?: string[];
  /** Dipole moment arrow (atomic units, centred at origin). */
  dipoleArrow?: { x: number; y: number; z: number; label: string };
  /** Show a colour-bar legend overlay. */
  colorBar?: { topLabel: string; bottomLabel: string; topColor: string; bottomColor: string };
  /** Disable auto-rotation on load (default: true = auto-rotate). */
  autoRotate?: boolean;
  /** Per-atom force/gradient arrows (e.g. for geometry optimisation). */
  forceArrows?: { x: number; y: number; z: number }[];
  /** Force arrow magnitude threshold — arrows below this are hidden. */
  forceThreshold?: number;
}

function detectBonds(positions: Vec3[], atomicNumbers: number[]): Bond[] {
  const bonds: Bond[] = [];
  const tolerance = 0.4;
  for (let i = 0; i < positions.length; i++) {
    for (let j = i + 1; j < positions.length; j++) {
      const dx = positions[i].x - positions[j].x;
      const dy = positions[i].y - positions[j].y;
      const dz = positions[i].z - positions[j].z;
      const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
      const maxDist = getCovalentRadius(atomicNumbers[i]) + getCovalentRadius(atomicNumbers[j]) + tolerance;
      if (dist < maxDist && dist > 0.1) bonds.push({ i, j });
    }
  }
  return bonds;
}

// Shared geometries (reused across instances)
const sphereGeo = new THREE.SphereGeometry(1, 24, 16);
const cylinderGeo = new THREE.CylinderGeometry(1, 1, 1, 12);

/** Create a text sprite using Canvas2D → Texture. */
function makeTextSprite(text: string, color = '#ffffff'): THREE.Sprite {
  const canvas = document.createElement('canvas');
  canvas.width = 256; canvas.height = 64;
  const ctx = canvas.getContext('2d')!;
  ctx.font = 'bold 36px -apple-system, BlinkMacSystemFont, sans-serif';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'middle';
  // Shadow for readability
  ctx.shadowColor = 'rgba(0,0,0,0.8)';
  ctx.shadowBlur = 4;
  ctx.fillStyle = color;
  ctx.fillText(text, 128, 32);
  const texture = new THREE.CanvasTexture(canvas);
  texture.minFilter = THREE.LinearFilter;
  const mat = new THREE.SpriteMaterial({ map: texture, transparent: true, depthTest: false });
  const sprite = new THREE.Sprite(mat);
  sprite.scale.set(0.8, 0.2, 1);
  return sprite;
}

function createViewer(
  container: HTMLElement,
  positions: Vec3[],
  atomicNumbers: number[],
  options?: ViewerOptions,
) {
  if (positions.length === 0) {
    container.innerHTML = `<p>${t('viewer.noAtoms')}</p>`;
    return;
  }

  // Center molecule
  let cx = 0, cy = 0, cz = 0;
  for (const p of positions) { cx += p.x; cy += p.y; cz += p.z; }
  cx /= positions.length; cy /= positions.length; cz /= positions.length;
  const centered = positions.map(p => ({ x: p.x - cx, y: p.y - cy, z: p.z - cz }));
  const bonds = detectBonds(centered, atomicNumbers);

  // Extent for camera distance
  let maxExtent = 0;
  for (const p of centered) {
    const r = Math.sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
    if (r > maxExtent) maxExtent = r;
  }
  const camDist = Math.max(maxExtent * 2.5, options?.atomLabels ? 2.0 : 3);

  // State
  let atomScale = 1.0;
  let bondScale = 1.0;

  // Build UI
  container.innerHTML = '';

  // Three.js setup
  const size = 400;
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(45, 1, 0.1, 1000);
  camera.position.set(0, 0, camDist);

  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setSize(size, size);
  renderer.setPixelRatio(window.devicePixelRatio);
  const canvas = renderer.domElement;
  canvas.style.width = '100%';
  canvas.style.height = '100%';
  canvas.style.cursor = 'grab';
  canvas.style.touchAction = 'none';
  canvas.style.display = 'block';

  // Wrapper: canvas + overlay controls
  const wrapper = document.createElement('div');
  wrapper.style.cssText = 'position:relative;width:100%;aspect-ratio:1;border-radius:8px;overflow:hidden';
  wrapper.appendChild(canvas);

  // Overlay controls (bottom of canvas)
  const controls = document.createElement('div');
  controls.style.cssText = 'position:absolute;bottom:0;left:0;right:0;display:flex;align-items:center;justify-content:center;gap:6px;padding:6px 10px;background:rgba(0,0,0,0.35);backdrop-filter:blur(4px);font-size:11px;color:#eee;pointer-events:auto';
  const sliderStyle = 'width:60px;accent-color:#7ec8e3;height:3px';
  const btnStyle = 'padding:2px 8px;font-size:11px;cursor:pointer;border:1px solid rgba(255,255,255,0.3);border-radius:4px;background:rgba(255,255,255,0.1);color:#eee';
  controls.innerHTML =
    `<span style="opacity:0.7">${t('viewer.atom')}</span><input type="range" min="30" max="200" value="100" style="${sliderStyle}" id="mol-atom-sl" title="${t('viewer.atomTip')}">` +
    `<span style="opacity:0.7">${t('viewer.bond')}</span><input type="range" min="30" max="200" value="100" style="${sliderStyle}" id="mol-bond-sl" title="${t('viewer.bondTip')}">` +
    `<button style="${btnStyle}" id="mol-rot-l" title="${t('viewer.rotLeft')}">\u25C0</button>` +
    `<button style="${btnStyle}" id="mol-rot-r" title="${t('viewer.rotRight')}">\u25B6</button>` +
    `<button style="${btnStyle}" id="mol-rot-u" title="${t('viewer.rotUp')}">\u25B2</button>` +
    `<button style="${btnStyle}" id="mol-rot-d" title="${t('viewer.rotDown')}">\u25BC</button>`;
  wrapper.appendChild(controls);

  // Colour-bar legend overlay (charge mode)
  if (options?.colorBar) {
    const cb = options.colorBar;
    const bar = document.createElement('div');
    bar.style.cssText = `position:absolute;top:12px;right:12px;display:flex;flex-direction:column;align-items:center;gap:2px;pointer-events:none`;
    bar.innerHTML =
      `<span style="font-size:9px;color:#eee;text-shadow:0 0 3px rgba(0,0,0,0.8)">${cb.topLabel}</span>` +
      `<div style="width:14px;height:80px;border-radius:3px;border:1px solid rgba(255,255,255,0.3);background:linear-gradient(to bottom,${cb.topColor},#ffffff,${cb.bottomColor})"></div>` +
      `<span style="font-size:9px;color:#eee;text-shadow:0 0 3px rgba(0,0,0,0.8)">${cb.bottomLabel}</span>`;
    wrapper.appendChild(bar);
  }

  container.appendChild(wrapper);

  const atomSlider = controls.querySelector<HTMLInputElement>('#mol-atom-sl')!;
  const bondSlider = controls.querySelector<HTMLInputElement>('#mol-bond-sl')!;
  const rotLBtn = controls.querySelector<HTMLButtonElement>('#mol-rot-l')!;
  const rotRBtn = controls.querySelector<HTMLButtonElement>('#mol-rot-r')!;
  const rotUBtn = controls.querySelector<HTMLButtonElement>('#mol-rot-u')!;
  const rotDBtn = controls.querySelector<HTMLButtonElement>('#mol-rot-d')!;

  // Resize renderer to match wrapper
  const ro = new ResizeObserver(() => {
    const w = wrapper.clientWidth;
    if (w > 0) {
      renderer.setSize(w, w);
      renderer.render(scene, camera);
    }
  });
  ro.observe(wrapper);

  // Lighting — bright enough to see colors clearly
  scene.add(new THREE.AmbientLight(0xffffff, 1.2));
  const dirLight = new THREE.DirectionalLight(0xffffff, 1.0);
  dirLight.position.set(2, 3, 4);
  scene.add(dirLight);
  const dirLight2 = new THREE.DirectionalLight(0xffffff, 0.6);
  dirLight2.position.set(-2, -1, -2);
  scene.add(dirLight2);

  // Molecule group (rotate this)
  const molGroup = new THREE.Group();
  scene.add(molGroup);

  // Atom meshes
  const atomMeshes: THREE.Mesh[] = [];
  for (let i = 0; i < centered.length; i++) {
    const p = centered[i];
    const z = atomicNumbers[i];
    const r = getAtomRadius(z);
    const col = options?.atomColors ? options.atomColors[i] : getColor(z);
    const mat = new THREE.MeshPhongMaterial({ color: col, shininess: 60 });
    const mesh = new THREE.Mesh(sphereGeo, mat);
    mesh.position.set(p.x, p.y, p.z);
    mesh.scale.setScalar(r);
    molGroup.add(mesh);
    atomMeshes.push(mesh);
  }

  // Atom text labels (sprites)
  if (options?.atomLabels) {
    const labelScale = Math.max(0.35, Math.min(0.8, maxExtent * 0.4));
    for (let i = 0; i < centered.length; i++) {
      const label = options.atomLabels[i];
      if (!label) continue;
      const sprite = makeTextSprite(label);
      sprite.scale.set(labelScale, labelScale * 0.25, 1);
      const p = centered[i];
      const r = getAtomRadius(atomicNumbers[i]);
      sprite.position.set(p.x, p.y - r - labelScale * 0.15, p.z);
      molGroup.add(sprite);
    }
  }

  // Dipole arrow (always rendered in front)
  if (options?.dipoleArrow) {
    const da = options.dipoleArrow;
    const mag = Math.sqrt(da.x * da.x + da.y * da.y + da.z * da.z);
    if (mag > 1e-6) {
      const dir = new THREE.Vector3(da.x, da.y, da.z).normalize();
      const arrowLen = Math.min(maxExtent * 1.2, mag * 1.5);
      const arrow = new THREE.ArrowHelper(dir, new THREE.Vector3(0, 0, 0), arrowLen, 0xFFD700, arrowLen * 0.18, arrowLen * 0.10);
      // Render arrow on top of everything
      arrow.renderOrder = 999;
      arrow.traverse(child => {
        if (child instanceof THREE.Mesh || child instanceof THREE.Line) {
          child.renderOrder = 999;
          (child.material as THREE.Material).depthTest = false;
        }
      });
      molGroup.add(arrow);
      // Dipole label sprite
      const tipPos = dir.clone().multiplyScalar(arrowLen + 0.2);
      const dSprite = makeTextSprite(da.label, '#FFD700');
      dSprite.position.copy(tipPos);
      molGroup.add(dSprite);
    }
  }

  // Force/gradient arrows per atom
  if (options?.forceArrows && options.forceArrows.length === positions.length) {
    const threshold = options.forceThreshold ?? 1e-4;
    // Compute max magnitude for scaling
    let maxMag = 0;
    const mags: number[] = [];
    for (const f of options.forceArrows) {
      const m = Math.sqrt(f.x * f.x + f.y * f.y + f.z * f.z);
      mags.push(m);
      if (m > maxMag) maxMag = m;
    }
    if (maxMag > 1e-8) {
      const arrowMaxLen = maxExtent * 0.6;
      const arrowMinLen = 0.15;
      const logTh = Math.log10(Math.max(threshold, 1e-10));
      const logRange = 3; // 3 orders of magnitude
      for (let i = 0; i < options.forceArrows.length; i++) {
        if (mags[i] < threshold * 0.1) continue;
        const f = options.forceArrows[i];
        const dir = new THREE.Vector3(f.x, f.y, f.z).normalize();
        // Log-scaled length
        const logF = Math.log10(Math.max(mags[i], 1e-10));
        const frac = Math.max(0, Math.min(1, (logF - logTh) / logRange));
        const arrowLen = arrowMinLen + frac * (arrowMaxLen - arrowMinLen);
        const origin = new THREE.Vector3(centered[i].x, centered[i].y, centered[i].z);
        const arrow = new THREE.ArrowHelper(dir, origin, arrowLen, 0x22c55e, arrowLen * 0.22, arrowLen * 0.12);
        arrow.renderOrder = 998;
        arrow.traverse(child => {
          if (child instanceof THREE.Mesh || child instanceof THREE.Line) {
            child.renderOrder = 998;
            (child.material as THREE.Material).depthTest = false;
          }
        });
        molGroup.add(arrow);
      }
    }
  }

  // Bond meshes (half-bonds: 2 cylinders per bond, each colored by its atom)
  const bondRadius = 0.08;
  const bondMeshes: THREE.Mesh[] = [];
  for (const bond of bonds) {
    const pi = centered[bond.i], pj = centered[bond.j];
    const a = new THREE.Vector3(pi.x, pi.y, pi.z);
    const b = new THREE.Vector3(pj.x, pj.y, pj.z);
    const mid = a.clone().add(b).multiplyScalar(0.5);

    // Half i→mid
    const ci = options?.atomColors ? options.atomColors[bond.i] : getColor(atomicNumbers[bond.i]);
    addCylinder(a, mid, ci);
    // Half mid→j
    const cj = options?.atomColors ? options.atomColors[bond.j] : getColor(atomicNumbers[bond.j]);
    addCylinder(mid, b, cj);
  }

  function addCylinder(start: THREE.Vector3, end: THREE.Vector3, color: number) {
    const dir = end.clone().sub(start);
    const len = dir.length();
    if (len < 0.001) return;
    const mat = new THREE.MeshPhongMaterial({ color, shininess: 40 });
    const mesh = new THREE.Mesh(cylinderGeo, mat);
    mesh.position.copy(start.clone().add(end).multiplyScalar(0.5));
    mesh.scale.set(bondRadius, len, bondRadius);
    // Align cylinder (default Y-axis) to bond direction
    const axis = new THREE.Vector3(0, 1, 0);
    const quat = new THREE.Quaternion().setFromUnitVectors(axis, dir.normalize());
    mesh.quaternion.copy(quat);
    molGroup.add(mesh);
    bondMeshes.push(mesh);
  }

  // Initial tilt (quaternion — fixed world axes)
  const initQ = new THREE.Quaternion();
  initQ.setFromEuler(new THREE.Euler(-0.35, 0.52, 0, 'YXZ'));
  molGroup.quaternion.copy(initQ);

  // Theme
  function updateBackground() {
    scene.background = new THREE.Color(isDark() ? 0x1e1e2e : 0xeef1f5);
  }
  updateBackground();

  // Render loop
  let autoRotH = 0; // horizontal (left/right around Y)
  let autoRotV = 0; // vertical (up/down around X)
  let animId = 0;

  const _worldY = new THREE.Vector3(0, 1, 0);
  const _worldX = new THREE.Vector3(1, 0, 0);
  const _tmpQ = new THREE.Quaternion();

  function isAnimating() { return autoRotH !== 0 || autoRotV !== 0; }

  function render() {
    if (autoRotH !== 0) {
      _tmpQ.setFromAxisAngle(_worldY, 0.012 * autoRotH);
      molGroup.quaternion.premultiply(_tmpQ);
    }
    if (autoRotV !== 0) {
      _tmpQ.setFromAxisAngle(_worldX, 0.012 * autoRotV);
      molGroup.quaternion.premultiply(_tmpQ);
    }
    renderer.render(scene, camera);
    if (isAnimating()) {
      animId = requestAnimationFrame(render);
    }
  }

  function requestRender() {
    if (!isAnimating()) renderer.render(scene, camera);
  }

  // Trackball drag
  let dragging = false;
  let lastX = 0, lastY = 0;

  canvas.addEventListener('pointerdown', (e) => {
    dragging = true;
    lastX = e.clientX;
    lastY = e.clientY;
    canvas.setPointerCapture(e.pointerId);
    canvas.style.cursor = 'grabbing';
  });

  canvas.addEventListener('pointermove', (e) => {
    if (!dragging) return;
    const dx = e.clientX - lastX;
    const dy = e.clientY - lastY;
    lastX = e.clientX;
    lastY = e.clientY;
    // Rotate around fixed world axes (not molecule-local axes)
    _tmpQ.setFromAxisAngle(_worldY, dx * 0.01);
    molGroup.quaternion.premultiply(_tmpQ);
    _tmpQ.setFromAxisAngle(_worldX, dy * 0.01);
    molGroup.quaternion.premultiply(_tmpQ);
    requestRender();
  });

  canvas.addEventListener('pointerup', () => { dragging = false; canvas.style.cursor = 'grab'; });
  canvas.addEventListener('pointercancel', () => { dragging = false; canvas.style.cursor = 'grab'; });

  // Zoom
  canvas.addEventListener('wheel', (e) => {
    e.preventDefault();
    const factor = e.deltaY > 0 ? 1.1 : 0.9;
    camera.position.z = Math.max(1, Math.min(camDist * 4, camera.position.z * factor));
    requestRender();
  }, { passive: false });

  // Scale sliders
  function rebuildScales() {
    for (let i = 0; i < atomMeshes.length; i++) {
      const r = getAtomRadius(atomicNumbers[i]) * atomScale;
      atomMeshes[i].scale.setScalar(r);
    }
    for (const m of bondMeshes) {
      const sy = m.scale.y; // keep length
      m.scale.set(bondRadius * bondScale, sy, bondRadius * bondScale);
    }
    requestRender();
  }

  atomSlider.addEventListener('input', () => {
    atomScale = parseInt(atomSlider.value, 10) / 100;
    rebuildScales();
  });
  bondSlider.addEventListener('input', () => {
    bondScale = parseInt(bondSlider.value, 10) / 100;
    rebuildScales();
  });

  // Auto-rotate buttons
  const activeRotStyle = 'background:rgba(126,200,227,0.5);color:#fff;border-color:rgba(126,200,227,0.6)';
  const inactiveRotStyle = 'background:rgba(255,255,255,0.1);color:#eee;border-color:rgba(255,255,255,0.3)';

  function updateRotBtns() {
    rotLBtn.style.cssText = btnStyle + ';' + (autoRotH === -1 ? activeRotStyle : inactiveRotStyle);
    rotRBtn.style.cssText = btnStyle + ';' + (autoRotH === 1 ? activeRotStyle : inactiveRotStyle);
    rotUBtn.style.cssText = btnStyle + ';' + (autoRotV === -1 ? activeRotStyle : inactiveRotStyle);
    rotDBtn.style.cssText = btnStyle + ';' + (autoRotV === 1 ? activeRotStyle : inactiveRotStyle);
  }

  function setAutoRotH(dir: number) {
    autoRotH = autoRotH === dir ? 0 : dir;
    updateRotBtns();
    cancelAnimationFrame(animId);
    if (isAnimating()) animId = requestAnimationFrame(render);
    else requestRender();
  }

  function setAutoRotV(dir: number) {
    autoRotV = autoRotV === dir ? 0 : dir;
    updateRotBtns();
    cancelAnimationFrame(animId);
    if (isAnimating()) animId = requestAnimationFrame(render);
    else requestRender();
  }

  rotLBtn.addEventListener('click', () => setAutoRotH(-1));
  rotRBtn.addEventListener('click', () => setAutoRotH(1));
  rotUBtn.addEventListener('click', () => setAutoRotV(-1));
  rotDBtn.addEventListener('click', () => setAutoRotV(1));
  if (options?.autoRotate === false) {
    setAutoRotH(0);          // start stopped
  } else {
    setAutoRotH(1);          // default: auto-rotate right
  }

  // Theme change
  const observer = new MutationObserver(() => { updateBackground(); requestRender(); });
  observer.observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
}

/** Render viewer from Atom[] (coordinates in Bohr) */
export function renderMoleculeViewer(container: HTMLElement, atoms: Atom[]): void {
  const positions = atoms.map(a => ({
    x: a.coordinate.x * BOHR_TO_ANGSTROM,
    y: a.coordinate.y * BOHR_TO_ANGSTROM,
    z: a.coordinate.z * BOHR_TO_ANGSTROM,
  }));
  createViewer(container, positions, atoms.map(a => a.atomicNumber));
}

/** Render charge-coloured viewer from ChargeResult-like data. */
export function renderChargeViewer(
  container: HTMLElement,
  data: {
    atomSymbols: string[];
    atomCoords: { x: number; y: number; z: number }[];
    charges: number[];
    dipole: { x: number; y: number; z: number; debye: number };
  },
): void {
  const positions = data.atomCoords;
  const atomicNumbers = data.atomSymbols.map(
    sym => ELEMENT_NAME_TO_ATOMIC_NUMBER[sym] ?? 1,
  );

  // Build charge colours (blue + → white 0 → red −)
  const atomColors = data.charges.map(q => chargeToHex(q));

  // Build labels: "El\n±0.123"
  const atomLabels = data.atomSymbols.map(
    (sym, i) => `${sym} ${data.charges[i] >= 0 ? '+' : ''}${data.charges[i].toFixed(3)}`,
  );

  const opts: ViewerOptions = {
    atomColors,
    atomLabels,
    colorBar: { topLabel: '+', bottomLabel: '\u2212', topColor: '#4466ff', bottomColor: '#ff4444' },
    autoRotate: false,
  };

  if (data.dipole.debye > 0.1) {
    opts.dipoleArrow = { x: data.dipole.x, y: data.dipole.y, z: data.dipole.z, label: `${data.dipole.debye.toFixed(2)} D` };
  }

  createViewer(container, positions, atomicNumbers, opts);
}

/** Convert Mulliken charge to 0xRRGGBB. Blue(+) → White(0) → Red(−). */
function chargeToHex(q: number): number {
  const clamp = Math.max(-0.6, Math.min(0.6, q));
  const t = clamp / 0.6; // -1..+1
  let r: number, g: number, b: number;
  if (t > 0) {
    // Positive → blue
    r = Math.round(255 * (1 - t * 0.6));
    g = Math.round(255 * (1 - t * 0.6));
    b = 255;
  } else {
    const s = -t;
    r = 255;
    g = Math.round(255 * (1 - s * 0.6));
    b = Math.round(255 * (1 - s * 0.6));
  }
  return (r << 16) | (g << 8) | b;
}

/** Render preview from raw XYZ text (coordinates in Angstrom) */
export function renderMoleculePreview(container: HTMLElement, xyzText: string): void {
  const lines = xyzText.split(/\r?\n/);
  if (lines.length < 3) { container.innerHTML = ''; return; }
  const atomCount = parseInt(lines[0].trim(), 10);
  if (isNaN(atomCount) || atomCount < 1) { container.innerHTML = ''; return; }

  const positions: Vec3[] = [];
  const atomicNumbers: number[] = [];
  for (let i = 2; i < 2 + atomCount && i < lines.length; i++) {
    const parts = lines[i]?.trim().split(/\s+/);
    if (!parts || parts.length < 4) continue;
    const z = ELEMENT_NAME_TO_ATOMIC_NUMBER[parts[0]];
    if (z === undefined) continue;
    const x = parseFloat(parts[1]);
    const y = parseFloat(parts[2]);
    const zz = parseFloat(parts[3]);
    if (isNaN(x) || isNaN(y) || isNaN(zz)) continue;
    positions.push({ x, y, z: zz });
    atomicNumbers.push(z);
  }

  createViewer(container, positions, atomicNumbers, { autoRotate: false });
}

/** Persistent geometry optimisation viewer — camera survives step changes. */
export interface GeomOptViewer {
  update(coords: number[], gradient: number[], forceThreshold: number): void;
}

export function createGeomOptViewer(
  container: HTMLElement,
  atomicNumbers: number[],
  initialCoords: number[],
  initialGradient: number[],
  forceThreshold: number,
): GeomOptViewer {
  const n = atomicNumbers.length;

  // ── Build scene infrastructure (once) ─────────────────
  container.innerHTML = '';

  const size = 400;
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(45, 1, 0.1, 1000);

  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setSize(size, size);
  renderer.setPixelRatio(window.devicePixelRatio);
  const canvas = renderer.domElement;
  canvas.style.cssText = 'width:100%;height:100%;cursor:grab;touch-action:none;display:block';

  const wrapper = document.createElement('div');
  wrapper.style.cssText = 'position:relative;width:100%;aspect-ratio:1;border-radius:8px;overflow:hidden';
  wrapper.appendChild(canvas);
  container.appendChild(wrapper);

  const ro = new ResizeObserver(() => {
    const w = wrapper.clientWidth;
    if (w > 0) { renderer.setSize(w, w); renderer.render(scene, camera); }
  });
  ro.observe(wrapper);

  // Lighting
  scene.add(new THREE.AmbientLight(0xffffff, 1.2));
  const dl1 = new THREE.DirectionalLight(0xffffff, 1.0); dl1.position.set(2, 3, 4); scene.add(dl1);
  const dl2 = new THREE.DirectionalLight(0xffffff, 0.6); dl2.position.set(-2, -1, -2); scene.add(dl2);

  // Background
  function updateBg() { scene.background = new THREE.Color(isDark() ? 0x1e1e2e : 0xeef1f5); }
  updateBg();
  const themeObs = new MutationObserver(() => { updateBg(); renderer.render(scene, camera); });
  themeObs.observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });

  // Molecule group — replaced on each update
  let molGroup = new THREE.Group();
  scene.add(molGroup);

  // Camera initialised once: compute initial orientation from first geometry
  let cameraInitialised = false;

  function initCamera(positions: Vec3[]) {
    // Center
    let cx = 0, cy = 0, cz = 0;
    for (const p of positions) { cx += p.x; cy += p.y; cz += p.z; }
    cx /= n; cy /= n; cz /= n;
    const centered = positions.map(p => ({ x: p.x - cx, y: p.y - cy, z: p.z - cz }));

    let maxExtent = 0;
    for (const p of centered) { const r = Math.sqrt(p.x * p.x + p.y * p.y + p.z * p.z); if (r > maxExtent) maxExtent = r; }
    const camDist = Math.max(maxExtent * 2.5, 3);
    camera.position.set(0, 0, camDist);

    // For diatomics: rotate so bond axis is horizontal (X)
    if (n === 2) {
      const bondDir = new THREE.Vector3(
        centered[1].x - centered[0].x,
        centered[1].y - centered[0].y,
        centered[1].z - centered[0].z,
      ).normalize();
      const xAxis = new THREE.Vector3(1, 0, 0);
      const q = new THREE.Quaternion().setFromUnitVectors(bondDir, xAxis);
      molGroup.quaternion.copy(q);
    } else {
      molGroup.quaternion.setFromEuler(new THREE.Euler(-0.35, 0.52, 0, 'YXZ'));
    }
    cameraInitialised = true;
  }

  // Trackball drag
  let dragging = false;
  let lastX = 0, lastY = 0;
  const _worldY = new THREE.Vector3(0, 1, 0);
  const _worldX = new THREE.Vector3(1, 0, 0);
  const _tmpQ = new THREE.Quaternion();

  canvas.addEventListener('pointerdown', (e) => {
    dragging = true; lastX = e.clientX; lastY = e.clientY;
    canvas.setPointerCapture(e.pointerId); canvas.style.cursor = 'grabbing';
  });
  canvas.addEventListener('pointermove', (e) => {
    if (!dragging) return;
    const dx = e.clientX - lastX, dy = e.clientY - lastY;
    lastX = e.clientX; lastY = e.clientY;
    _tmpQ.setFromAxisAngle(_worldY, dx * 0.01); molGroup.quaternion.premultiply(_tmpQ);
    _tmpQ.setFromAxisAngle(_worldX, dy * 0.01); molGroup.quaternion.premultiply(_tmpQ);
    renderer.render(scene, camera);
  });
  canvas.addEventListener('pointerup', () => { dragging = false; canvas.style.cursor = 'grab'; });
  canvas.addEventListener('pointercancel', () => { dragging = false; canvas.style.cursor = 'grab'; });

  // Zoom
  canvas.addEventListener('wheel', (e) => {
    e.preventDefault();
    const factor = e.deltaY > 0 ? 1.1 : 0.9;
    camera.position.z = Math.max(1, Math.min(camera.position.z * 4, camera.position.z * factor));
    renderer.render(scene, camera);
  }, { passive: false });

  // ── Rebuild molecule content ──────────────────────────
  function buildMolecule(coords: number[], gradient: number[], threshold: number) {
    // Remove old molecule content, keep quaternion
    const savedQ = molGroup.quaternion.clone();
    scene.remove(molGroup);
    molGroup = new THREE.Group();
    molGroup.quaternion.copy(savedQ);
    scene.add(molGroup);

    // Convert coords (Bohr → Angstrom)
    const positions: Vec3[] = [];
    for (let i = 0; i < n; i++) {
      positions.push({
        x: coords[i * 3] * BOHR_TO_ANGSTROM,
        y: coords[i * 3 + 1] * BOHR_TO_ANGSTROM,
        z: coords[i * 3 + 2] * BOHR_TO_ANGSTROM,
      });
    }

    // Init camera only on first call
    if (!cameraInitialised) initCamera(positions);

    // Center molecule at origin
    let cx = 0, cy = 0, cz = 0;
    for (const p of positions) { cx += p.x; cy += p.y; cz += p.z; }
    cx /= n; cy /= n; cz /= n;
    const centered = positions.map(p => ({ x: p.x - cx, y: p.y - cy, z: p.z - cz }));

    // Atom spheres
    for (let i = 0; i < n; i++) {
      const z = atomicNumbers[i];
      const r = getAtomRadius(z);
      const mat = new THREE.MeshPhongMaterial({ color: getColor(z), shininess: 60 });
      const mesh = new THREE.Mesh(sphereGeo, mat);
      mesh.position.set(centered[i].x, centered[i].y, centered[i].z);
      mesh.scale.setScalar(r);
      molGroup.add(mesh);
    }

    // Bonds
    const bonds = detectBonds(centered, atomicNumbers);
    const bondRadius = 0.08;
    for (const bond of bonds) {
      const pi = centered[bond.i], pj = centered[bond.j];
      const a = new THREE.Vector3(pi.x, pi.y, pi.z);
      const b = new THREE.Vector3(pj.x, pj.y, pj.z);
      const mid = a.clone().add(b).multiplyScalar(0.5);
      addCyl(a, mid, getColor(atomicNumbers[bond.i]));
      addCyl(mid, b, getColor(atomicNumbers[bond.j]));
    }

    function addCyl(start: THREE.Vector3, end: THREE.Vector3, color: number) {
      const dir = end.clone().sub(start);
      const len = dir.length();
      if (len < 0.001) return;
      const mat = new THREE.MeshPhongMaterial({ color, shininess: 40 });
      const mesh = new THREE.Mesh(cylinderGeo, mat);
      mesh.position.copy(start.clone().add(end).multiplyScalar(0.5));
      mesh.scale.set(bondRadius, len, bondRadius);
      mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir.normalize());
      molGroup.add(mesh);
    }

    // Force arrows
    let maxExtent = 0;
    for (const p of centered) { const r = Math.sqrt(p.x * p.x + p.y * p.y + p.z * p.z); if (r > maxExtent) maxExtent = r; }

    const mags: number[] = [];
    for (let i = 0; i < n; i++) {
      const gx = gradient[i * 3], gy = gradient[i * 3 + 1], gz = gradient[i * 3 + 2];
      mags.push(Math.sqrt(gx * gx + gy * gy + gz * gz));
    }
    const maxMag = Math.max(...mags);
    if (maxMag > 1e-8) {
      const arrowMaxLen = Math.max(maxExtent * 0.6, 0.5);
      const arrowMinLen = 0.15;
      const logTh = Math.log10(Math.max(threshold, 1e-10));
      const logRange = 3;
      for (let i = 0; i < n; i++) {
        if (mags[i] < threshold * 0.1) continue;
        // Force = -gradient
        const dir = new THREE.Vector3(-gradient[i * 3], -gradient[i * 3 + 1], -gradient[i * 3 + 2]).normalize();
        const logF = Math.log10(Math.max(mags[i], 1e-10));
        const frac = Math.max(0, Math.min(1, (logF - logTh) / logRange));
        const arrowLen = arrowMinLen + frac * (arrowMaxLen - arrowMinLen);
        const origin = new THREE.Vector3(centered[i].x, centered[i].y, centered[i].z);
        const arrow = new THREE.ArrowHelper(dir, origin, arrowLen, 0x22c55e, arrowLen * 0.22, arrowLen * 0.12);
        arrow.renderOrder = 998;
        arrow.traverse(child => {
          if (child instanceof THREE.Mesh || child instanceof THREE.Line) {
            child.renderOrder = 998;
            (child.material as THREE.Material).depthTest = false;
          }
        });
        molGroup.add(arrow);
      }
    }

    renderer.render(scene, camera);
  }

  // Build initial state
  buildMolecule(initialCoords, initialGradient, forceThreshold);

  return {
    update(coords: number[], gradient: number[], ft: number) {
      buildMolecule(coords, gradient, ft);
    },
  };
}
