/** Periodic table modal — pick an element from a visual grid. */

import { t } from './i18n';
import { enhanceRadioGroup } from './radioGroup';

/* ── Covalent radii (Å) for scan-range estimation ── */

export const COVALENT_RADII: Record<number, number> = {
  1: 0.31, 2: 0.28, 3: 1.28, 4: 0.96, 5: 0.84, 6: 0.76, 7: 0.71, 8: 0.66,
  9: 0.57, 10: 0.58, 11: 1.66, 12: 1.41, 13: 1.21, 14: 1.11, 15: 1.07, 16: 1.05,
  17: 1.02, 18: 1.06, 19: 2.03, 20: 1.76, 21: 1.70, 22: 1.60, 23: 1.53, 24: 1.39,
  25: 1.39, 26: 1.32, 27: 1.26, 28: 1.24, 29: 1.32, 30: 1.22, 31: 1.22, 32: 1.20,
  33: 1.19, 34: 1.20, 35: 1.20, 36: 1.16,
};

/* ── Periodic table layout (period 1-4, 18 columns) ── */
// Each entry: [row, col, Z, symbol]   (0-indexed row/col)

const LAYOUT: [number, number, number, string][] = [
  // Period 1
  [0, 0, 1, 'H'], [0, 17, 2, 'He'],
  // Period 2
  [1, 0, 3, 'Li'], [1, 1, 4, 'Be'],
  [1, 12, 5, 'B'], [1, 13, 6, 'C'], [1, 14, 7, 'N'], [1, 15, 8, 'O'], [1, 16, 9, 'F'], [1, 17, 10, 'Ne'],
  // Period 3
  [2, 0, 11, 'Na'], [2, 1, 12, 'Mg'],
  [2, 12, 13, 'Al'], [2, 13, 14, 'Si'], [2, 14, 15, 'P'], [2, 15, 16, 'S'], [2, 16, 17, 'Cl'], [2, 17, 18, 'Ar'],
  // Period 4 (displayed but disabled beyond maxZ)
  [3, 0, 19, 'K'], [3, 1, 20, 'Ca'],
  [3, 2, 21, 'Sc'], [3, 3, 22, 'Ti'], [3, 4, 23, 'V'], [3, 5, 24, 'Cr'], [3, 6, 25, 'Mn'],
  [3, 7, 26, 'Fe'], [3, 8, 27, 'Co'], [3, 9, 28, 'Ni'], [3, 10, 29, 'Cu'], [3, 11, 30, 'Zn'],
  [3, 12, 31, 'Ga'], [3, 13, 32, 'Ge'], [3, 14, 33, 'As'], [3, 15, 34, 'Se'], [3, 16, 35, 'Br'], [3, 17, 36, 'Kr'],
];

/** Build a single periodic table grid element.
 *  When a cell is clicked, `onPick(sym, z)` is called. */
function buildGrid(
  onPick: (sym: string, z: number) => void,
  label: string,
  key: string,
  selectedZ?: number,
  maxZ = 18,
): HTMLDivElement {
  const grid = document.createElement('div');
  grid.className = 'pte-grid';

  for (const [row, col, z, sym] of LAYOUT) {
    const cell = document.createElement('div');
    cell.className = 'pte-cell';
    cell.style.gridRow = String(row + 1);
    cell.style.gridColumn = String(col + 1);

    const disabled = z > maxZ;
    if (disabled) cell.classList.add('pte-disabled');

    // Colour code by block
    const blockClass = z <= 2 ? 'pte-s'
      : (z <= 4 || (z >= 11 && z <= 12) || (z >= 19 && z <= 20)) ? 'pte-s'
      : (z >= 21 && z <= 30) ? 'pte-d'
      : 'pte-p';
    cell.classList.add(blockClass);

    if (z === selectedZ) cell.classList.add('pte-selected');

    cell.innerHTML = `<span class="pte-z">${z}</span><span class="pte-sym">${sym}</span>`;
    if (!disabled) cell.addEventListener('click', () => onPick(sym, z));

    grid.appendChild(cell);
  }

  enhanceRadioGroup(grid, {
    label,
    key,
    itemSelector: '.pte-cell',
    selectedClass: 'pte-selected',
    disabledClass: 'pte-disabled',
    // The layout is a real periodic table, so up/down should change period.
    coords: (el) => ({
      row: Number(el.style.gridRow) || 0,
      col: Number(el.style.gridColumn) || 0,
    }),
  });

  return grid;
}

/**
 * Show a dual periodic table modal (side-by-side) for picking two atoms.
 * Calls `onSelect(symA, zA, symB, zB)` once both are chosen.
 */
export function showDualPeriodicTable(
  onSelect: (symA: string, zA: number, symB: string, zB: number) => void,
  initialA?: { symbol: string; z: number } | null,
  initialB?: { symbol: string; z: number } | null,
): void {
  let pickA: { symbol: string; z: number } | null = initialA ?? null;
  let pickB: { symbol: string; z: number } | null = initialB ?? null;

  const overlay = document.createElement('div');
  overlay.className = 'pte-overlay';

  const modal = document.createElement('div');
  modal.className = 'pte-modal pte-dual';
  modal.setAttribute('role', 'dialog');
  modal.setAttribute('aria-modal', 'true');
  modal.setAttribute('aria-label', t('opt.customSelect'));

  function rebuild() {
    modal.innerHTML = '';

    // Tables row
    const tablesRow = document.createElement('div');
    tablesRow.className = 'pte-tables';

    // Left: Atom A
    const colA = document.createElement('div');
    colA.className = 'pte-col';
    const labelA = document.createElement('div');
    labelA.className = 'pte-title';
    labelA.textContent = `${t('opt.customAtomA')}${pickA ? ': ' + pickA.symbol : ''}`;
    colA.appendChild(labelA);
    colA.appendChild(buildGrid((sym, z) => {
      pickA = { symbol: sym, z };
      rebuild();
    }, t('opt.customAtomA'), 'pte-atom-a', pickA?.z, 20));  // Period 1-3 + K,Ca
    tablesRow.appendChild(colA);

    // Right: Atom B
    const colB = document.createElement('div');
    colB.className = 'pte-col';
    const labelB = document.createElement('div');
    labelB.className = 'pte-title';
    labelB.textContent = `${t('opt.customAtomB')}${pickB ? ': ' + pickB.symbol : ''}`;
    colB.appendChild(labelB);
    colB.appendChild(buildGrid((sym, z) => {
      pickB = { symbol: sym, z };
      rebuild();
    }, t('opt.customAtomB'), 'pte-atom-b', pickB?.z, 12));  // Period 1-2 + Na,Mg
    tablesRow.appendChild(colB);

    modal.appendChild(tablesRow);

    // OK / Cancel buttons
    const btnRow = document.createElement('div');
    btnRow.className = 'pte-btn-row';

    const cancelBtn = document.createElement('button');
    cancelBtn.className = 'pte-btn pte-btn-cancel';
    cancelBtn.textContent = 'Cancel';
    cancelBtn.addEventListener('click', () => close());
    btnRow.appendChild(cancelBtn);

    const okBtn = document.createElement('button');
    okBtn.className = 'pte-btn pte-btn-ok';
    okBtn.textContent = 'OK';
    okBtn.disabled = !pickA || !pickB;
    okBtn.addEventListener('click', () => {
      if (pickA && pickB) {
        close();
        onSelect(pickA.symbol, pickA.z, pickB.symbol, pickB.z);
      }
    });
    btnRow.appendChild(okBtn);

    modal.appendChild(btnRow);
  }

  // Without this the dialog is a dead end for anyone not using a mouse: there is
  // nothing focusable to reach, no way to pick an element and no way back out.
  const previouslyFocused = document.activeElement as HTMLElement | null;

  function close(): void {
    overlay.remove();
    document.removeEventListener('keydown', onKeydown, true);
    previouslyFocused?.focus?.();
  }

  /** Everything the dialog offers the keyboard: the two grids' tab stops and the buttons. */
  function focusables(): HTMLElement[] {
    return Array.from(
      modal.querySelectorAll<HTMLElement>('[tabindex="0"], button:not([disabled])'),
    );
  }

  function onKeydown(e: KeyboardEvent): void {
    if (e.key === 'Escape') {
      e.preventDefault();
      close();
      return;
    }
    if (e.key !== 'Tab') return;
    // Keep Tab inside the dialog; otherwise focus wanders into the page behind
    // it, which is still there and still clickable.
    const stops = focusables();
    if (stops.length === 0) return;
    const first = stops[0], last = stops[stops.length - 1];
    const active = document.activeElement as HTMLElement | null;
    if (!modal.contains(active)) {
      e.preventDefault();
      (e.shiftKey ? last : first).focus();
    } else if (e.shiftKey && active === first) {
      e.preventDefault();
      last.focus();
    } else if (!e.shiftKey && active === last) {
      e.preventDefault();
      first.focus();
    }
  }

  rebuild();
  overlay.appendChild(modal);

  overlay.addEventListener('click', (e) => {
    if (e.target === overlay) close();
  });

  document.addEventListener('keydown', onKeydown, true);
  document.body.appendChild(overlay);
  focusables()[0]?.focus();
}
