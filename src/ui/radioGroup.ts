/**
 * Makes the option grids keyboard-operable.
 *
 * The grids are built independently on nine pages — a `<div class="opt-scenario-card">`
 * per option with a click listener — so there is no component to fix. Rather than
 * rewrite ten call sites, this upgrades a grid *after* it has been built: the
 * existing click handlers stay untouched and keyboard activation is routed
 * through them via `click()`.
 *
 * Implements the ARIA radiogroup pattern: arrow keys move and select, Home/End
 * jump to the ends, Enter/Space selects, and only the checked option is in the
 * tab order (roving tabindex) so Tab moves past the whole group rather than
 * through every molecule.
 */

/**
 * Selecting an option re-renders the page, which destroys the focused element.
 * Remembering which group was being driven from the keyboard lets the next
 * enhancement of that group restore focus, so arrow-key navigation survives the
 * re-render instead of dumping focus on `<body>`.
 */
let pendingFocusGroup: string | null = null;

export interface RadioGroupOptions {
  /** Accessible name for the group as a whole. */
  label: string;
  /** Which descendants are the options. Defaults to the shared card class. */
  itemSelector?: string;
  /**
   * Identity of the group across re-renders, used to restore focus. Defaults to
   * `label`; pass it explicitly when one page has two grids with the same label.
   */
  key?: string;
  /** Class marking the checked option. Defaults to the shared `selected`. */
  selectedClass?: string;
  /** Class marking an unselectable option. Defaults to the shared `disabled`. */
  disabledClass?: string;
  /**
   * Row/column of an option, when the group is genuinely two-dimensional (the
   * periodic table). Given this, up/down move between rows instead of stepping
   * through DOM order, which in a periodic table would just repeat left/right.
   */
  coords?: (el: HTMLElement) => { row: number; col: number };
}

function makeIsDisabled(cls: string) {
  return (el: HTMLElement): boolean =>
    el.classList.contains(cls) || el.getAttribute('aria-disabled') === 'true';
}

/**
 * Upgrade an already-populated option grid in place. Safe to call on every
 * render, and a no-op for an empty grid.
 */
export function enhanceRadioGroup(container: HTMLElement, options: RadioGroupOptions): void {
  injectStyles();

  const items = Array.from(
    container.querySelectorAll<HTMLElement>(options.itemSelector ?? '.opt-scenario-card'),
  );
  if (items.length === 0) return;

  const key = options.key ?? options.label;
  const selectedClass = options.selectedClass ?? 'selected';
  const isDisabled = makeIsDisabled(options.disabledClass ?? 'disabled');
  container.setAttribute('role', 'radiogroup');
  container.setAttribute('aria-label', options.label);

  // Several grids put layout rows and category captions inside the container.
  // A radiogroup is only supposed to contain radios, so strip the semantics of
  // everything else rather than leaving stray generics in the accessibility
  // tree; `presentation` keeps their contents but drops the element itself.
  for (const child of Array.from(container.children)) {
    if (!child.hasAttribute('role')) child.setAttribute('role', 'presentation');
  }

  const checkedIndex = items.findIndex(el => el.classList.contains(selectedClass));
  // With nothing selected yet, the first enabled option is the one Tab reaches.
  const tabStop = checkedIndex >= 0 ? checkedIndex : items.findIndex(el => !isDisabled(el));

  items.forEach((el, i) => {
    el.setAttribute('role', 'radio');
    el.setAttribute('aria-checked', String(el.classList.contains(selectedClass)));
    if (isDisabled(el)) el.setAttribute('aria-disabled', 'true');
    el.tabIndex = i === tabStop ? 0 : -1;

    el.addEventListener('keydown', (ev: KeyboardEvent) => {
      if (ev.altKey || ev.ctrlKey || ev.metaKey) return;

      if (ev.key === 'Enter' || ev.key === ' ' || ev.key === 'Spacebar') {
        ev.preventDefault();
        if (isDisabled(el)) return;
        pendingFocusGroup = key;
        el.click();
        return;
      }

      const step =
        ev.key === 'ArrowRight' || ev.key === 'ArrowDown' ? 1 :
        ev.key === 'ArrowLeft' || ev.key === 'ArrowUp' ? -1 : 0;

      const vertical = ev.key === 'ArrowUp' || ev.key === 'ArrowDown';

      let target: HTMLElement | undefined;
      if (step !== 0 && vertical && options.coords) {
        // Nearest enabled cell in the next occupied row, matching column as
        // closely as possible — a periodic table has ragged rows, so "the cell
        // above" is not a fixed offset in DOM order.
        const here = options.coords(el);
        let bestRow = Infinity;
        for (const other of items) {
          if (isDisabled(other)) continue;
          const c = options.coords(other);
          const ahead = step > 0 ? c.row > here.row : c.row < here.row;
          if (!ahead) continue;
          if (Math.abs(c.row - here.row) < Math.abs(bestRow - here.row)) bestRow = c.row;
        }
        if (bestRow !== Infinity) {
          let bestCol = Infinity;
          for (const other of items) {
            if (isDisabled(other)) continue;
            const c = options.coords(other);
            if (c.row !== bestRow) continue;
            if (Math.abs(c.col - here.col) < Math.abs(bestCol - here.col)) {
              bestCol = c.col;
              target = other;
            }
          }
        }
      } else if (step !== 0) {
        // Skip disabled options, and stop at the ends rather than wrapping —
        // these grids are visually two-dimensional, so wrapping would jump the
        // focus ring across the panel and read as a glitch.
        for (let j = i + step; j >= 0 && j < items.length; j += step) {
          if (!isDisabled(items[j])) { target = items[j]; break; }
        }
      } else if (ev.key === 'Home') {
        target = items.find(x => !isDisabled(x));
      } else if (ev.key === 'End') {
        target = [...items].reverse().find(x => !isDisabled(x));
      } else {
        return;
      }

      ev.preventDefault();
      if (!target || target === el) return;
      // The radiogroup pattern selects on arrow, and selecting re-renders, so
      // move focus first for the case where the click is a no-op.
      target.tabIndex = 0;
      el.tabIndex = -1;
      target.focus();
      pendingFocusGroup = key;
      target.click();
    });
  });

  if (pendingFocusGroup === key) {
    pendingFocusGroup = null;
    const focusIndex = checkedIndex >= 0 ? checkedIndex : tabStop;
    const toFocus = items[focusIndex];
    // Callers may enhance a grid before inserting it (the periodic table builds
    // its tables then appends them), and focus() on a detached node is a no-op.
    // Deferring until the current task finishes means the node is in place.
    if (toFocus) {
      if (toFocus.isConnected) toFocus.focus();
      else requestAnimationFrame(() => toFocus.focus());
    }
  }
}

let stylesInjected = false;

function injectStyles(): void {
  if (stylesInjected) return;
  stylesInjected = true;
  if (document.getElementById('radio-group-styles')) return;
  const style = document.createElement('style');
  style.id = 'radio-group-styles';
  // :focus-visible so pointer users never see the ring; the 2px accent outline
  // matches how selection is already drawn elsewhere in the app.
  style.textContent = `
    [role="radio"]:focus-visible {
      outline: 2px solid var(--color-accent);
      outline-offset: 2px;
    }
    [role="radio"][aria-disabled="true"] {
      cursor: not-allowed;
    }
  `;
  document.head.appendChild(style);
}
