/**
 * Citation metadata and the "How to cite" dialog.
 *
 * Mirrors MOrbVis (web_viewer/src/core/citation.ts) so the two tools present
 * their references the same way: one record here drives the formatted
 * reference, the BibTeX entry and the dialog, and the README / CITATION.cff in
 * the public repository carry the same fields.
 */

import { t } from './i18n';

interface CitationInfo {
  bibtexKey: string;
  authors: string[];
  title: string;
  journal: string;
  volume: string;
  pages: string;
  year: number;
  issn: string;
  doi: string;
  url: string;
}

export const CITATION: CitationInfo = {
  bibtexKey: 'ITO2026103046',
  authors: ['Yasuaki Ito', 'Haruto Fujii', 'Satoki Tsuji', 'Koji Nakano', 'Akihiko Kasagi'],
  title: 'GANSU lite: A zero-install, browser-based quantum chemistry platform',
  journal: 'SoftwareX',
  volume: '36',
  pages: '103046',
  year: 2026,
  issn: '2352-7110',
  doi: '10.1016/j.softx.2026.103046',
  url: 'https://www.sciencedirect.com/science/article/pii/S2352711026005376',
};

/** Single-line reference, e.g. "... SoftwareX 2026, 36, 103046. DOI: ..." */
export function formattedCitation(): string {
  const c = CITATION;
  return `${c.authors.join(', ')}. "${c.title}". ` +
    `${c.journal} ${c.year}, ${c.volume}, ${c.pages}. DOI: ${c.doi}`;
}

export function bibtexEntry(): string {
  const c = CITATION;
  // Double braces keep "GANSU" capitalised under styles that down-case titles.
  return [
    `@article{${c.bibtexKey},`,
    `  author  = {${c.authors.join(' and ')}},`,
    `  title   = {{${c.title}}},`,
    `  journal = {${c.journal}},`,
    `  volume  = {${c.volume}},`,
    `  pages   = {${c.pages}},`,
    `  year    = {${c.year}},`,
    `  issn    = {${c.issn}},`,
    `  doi     = {${c.doi}},`,
    `  url     = {${c.url}},`,
    '}',
  ].join('\n');
}

// ── Dialog ──

/**
 * Copy text, falling back to selecting it when the Clipboard API is refused
 * (insecure context, or permission denied) so the user can still press Ctrl+C.
 */
async function copyText(text: string, fallbackTarget: HTMLElement): Promise<boolean> {
  try {
    await navigator.clipboard.writeText(text);
    return true;
  } catch {
    const range = document.createRange();
    range.selectNodeContents(fallbackTarget);
    const sel = window.getSelection();
    sel?.removeAllRanges();
    sel?.addRange(range);
    return false;
  }
}

function escapeHtml(s: string): string {
  return s.replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]!));
}

export function showCiteDialog(): void {
  if (document.querySelector('.cite-overlay')) return;
  const previouslyFocused = document.activeElement as HTMLElement | null;

  const overlay = document.createElement('div');
  overlay.className = 'cite-overlay';
  overlay.innerHTML = `
    <div class="cite-modal" role="dialog" aria-modal="true" aria-labelledby="cite-title">
      <div class="cite-head">
        <h2 id="cite-title">${t('cite.title')}</h2>
        <button type="button" class="cite-close" aria-label="${t('cite.close')}">✕</button>
      </div>
      <p class="cite-intro">${t('cite.intro')}</p>
      <pre class="cite-ref" id="cite-ref">${escapeHtml(formattedCitation())}</pre>
      <p class="cite-intro">BibTeX</p>
      <pre class="cite-bib" id="cite-bib">${escapeHtml(bibtexEntry())}</pre>
      <div class="cite-actions">
        <button type="button" data-act="copy-ref">${t('cite.copyText')}</button>
        <button type="button" data-act="copy-bib">${t('cite.copyBibtex')}</button>
        <button type="button" data-act="download" class="cite-primary">${t('cite.downloadBibtex')}</button>
        <a href="https://doi.org/${CITATION.doi}" target="_blank" rel="noopener noreferrer">${t('cite.openDoi')}</a>
      </div>
      <p class="cite-status" role="status" aria-live="polite"></p>
    </div>`;

  const modal = overlay.querySelector<HTMLElement>('.cite-modal')!;
  const status = overlay.querySelector<HTMLElement>('.cite-status')!;
  const refPre = overlay.querySelector<HTMLElement>('#cite-ref')!;
  const bibPre = overlay.querySelector<HTMLElement>('#cite-bib')!;

  function close(): void {
    overlay.remove();
    document.removeEventListener('keydown', onKeydown, true);
    previouslyFocused?.focus?.();
  }

  function focusables(): HTMLElement[] {
    return Array.from(modal.querySelectorAll<HTMLElement>('button, a[href]'));
  }

  function onKeydown(e: KeyboardEvent): void {
    if (e.key === 'Escape') { e.preventDefault(); close(); return; }
    if (e.key !== 'Tab') return;
    const stops = focusables();
    const first = stops[0], last = stops[stops.length - 1];
    const active = document.activeElement as HTMLElement | null;
    if (!modal.contains(active)) { e.preventDefault(); (e.shiftKey ? last : first).focus(); }
    else if (e.shiftKey && active === first) { e.preventDefault(); last.focus(); }
    else if (!e.shiftKey && active === last) { e.preventDefault(); first.focus(); }
  }

  async function copy(text: string, target: HTMLElement): Promise<void> {
    const ok = await copyText(text, target);
    status.textContent = ok ? t('cite.copied') : t('cite.copyFallback');
  }

  overlay.addEventListener('click', (e) => { if (e.target === overlay) close(); });
  overlay.querySelector('.cite-close')!.addEventListener('click', close);
  overlay.querySelector('[data-act="copy-ref"]')!.addEventListener('click', () => { void copy(formattedCitation(), refPre); });
  overlay.querySelector('[data-act="copy-bib"]')!.addEventListener('click', () => { void copy(bibtexEntry(), bibPre); });
  overlay.querySelector('[data-act="download"]')!.addEventListener('click', () => {
    const blob = new Blob([bibtexEntry() + '\n'], { type: 'application/x-bibtex' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${CITATION.bibtexKey}.bib`;
    a.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  });

  document.addEventListener('keydown', onKeydown, true);
  document.body.appendChild(overlay);
  overlay.querySelector<HTMLElement>('.cite-close')!.focus();
}
