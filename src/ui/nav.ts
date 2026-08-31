/** Shared page navigation — tab bar header + footer for all GANSU Lite pages. */

import { t, getLang } from './i18n';
import { isDark } from './theme';

const FOOTER_ID = 'gansu-footer';
const REPO_URL = 'https://github.com/Yasuaki-Ito/GANSU-Lite';
const BOOK_URL_JA = 'https://yasuaki-ito.github.io/book/qcbook/';
const BOOK_URL_EN = 'https://yasuaki-ito.github.io/book/en/qcbook/';

/** GitHub Octocat mark, inline SVG (currentColor for theme support). */
const GITHUB_ICON = `<svg viewBox="0 0 16 16" width="14" height="14" fill="currentColor" aria-hidden="true" style="vertical-align:-2px;margin-right:3px;">
  <path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"/>
</svg>`;

/** Octicon "law" (gavel + scales) for license link. */
const LAW_ICON = `<svg viewBox="0 0 16 16" width="14" height="14" fill="currentColor" aria-hidden="true" style="vertical-align:-2px;margin-right:3px;">
  <path d="M8.75.75V2h.985c.304 0 .603.08.867.231l1.29.736c.038.022.08.033.124.033h2.234a.75.75 0 0 1 0 1.5h-.427l2.111 4.692a.75.75 0 0 1-.154.838l-.53-.53.529.531-.001.002-.002.002-.006.006-.006.005-.01.01-.045.04c-.21.176-.441.327-.686.45C14.556 10.78 13.88 11 13 11a4.498 4.498 0 0 1-2.023-.454 3.544 3.544 0 0 1-.686-.45l-.045-.04-.016-.015-.006-.006-.004-.004v-.001a.75.75 0 0 1-.154-.838L12.178 4.5h-.162c-.305 0-.604-.079-.868-.231l-1.29-.736a.245.245 0 0 0-.124-.033H8.75V13h2.5a.75.75 0 0 1 0 1.5h-6.5a.75.75 0 0 1 0-1.5h2.5V3.5h-.984a.245.245 0 0 0-.124.033l-1.289.737c-.265.15-.564.23-.869.23h-.162l2.112 4.692a.75.75 0 0 1-.154.838l-.53-.53.529.531-.001.002-.002.002-.006.006-.016.015-.045.04c-.21.176-.441.327-.686.45C4.556 10.78 3.88 11 3 11a4.498 4.498 0 0 1-2.023-.454 3.544 3.544 0 0 1-.686-.45l-.045-.04-.016-.015-.006-.006-.004-.004v-.001a.75.75 0 0 1-.154-.838L2.178 4.5H1.75a.75.75 0 0 1 0-1.5h2.234a.249.249 0 0 0 .125-.033l1.288-.737c.265-.15.564-.23.869-.23h.984V.75a.75.75 0 0 1 1.5 0Zm2.945 8.477c.285.135.718.273 1.305.273s1.02-.138 1.305-.273L13 6.327Zm-10 0c.285.135.718.273 1.305.273s1.02-.138 1.305-.273L3 6.327Z"/>
</svg>`;

/** Build footer HTML using the current UI language (from i18n.getLang). */
function footerHTML(): string {
  const isJa = getLang() === 'ja';
  const disclaimer = isJa
    ? '本ソフトウェアは教育・研究目的で提供されており、結果の正確性は保証しません。'
    : 'For educational and research use. No warranty of correctness.';
  const licenseLabel = isJa ? 'BSD 3-Clause ライセンス' : 'BSD 3-Clause License';
  const bookLabel = isJa ? '解説書' : 'Companion textbook';
  const bookHref = isJa ? BOOK_URL_JA : BOOK_URL_EN;
  return `
    <span>© 2026 Yasuaki Ito</span>
    <span class="sep">·</span>
    <a href="${bookHref}" target="_blank" rel="noopener" class="book-link">📖 ${bookLabel}</a>
    <span class="sep">·</span>
    <a href="${REPO_URL}/blob/main/LICENSE" target="_blank" rel="noopener">${LAW_ICON}${licenseLabel}</a>
    <span class="sep">·</span>
    <a href="${REPO_URL}" target="_blank" rel="noopener">${GITHUB_ICON}GitHub</a>
    <span class="sep">·</span>
    <span class="disclaimer">${disclaimer}</span>
  `;
}

/** Auto-mounted on DOMContentLoaded — appends a small footer to <body>.
 *  Re-renders text on language change without rebuilding the element. */
function mountFooter(): void {
  let footer = document.getElementById(FOOTER_ID);
  if (!footer) {
    footer = document.createElement('footer');
    footer.id = FOOTER_ID;
    footer.className = 'gansu-footer';
    document.body.appendChild(footer);
  }
  footer.innerHTML = footerHTML();
}

if (typeof document !== 'undefined') {
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', mountFooter);
  } else {
    mountFooter();
  }
  // Update on language toggle
  document.addEventListener('gansu-lang-change', mountFooter);
}

export type PageId = 'calc' | 'optimize' | 'geomopt' | 'freqanalysis' | 'walsh' | 'accuracy' | 'charges' | 'convergence';

// Tab order follows the pedagogical "basic -> advanced" sequence of the
// companion paper's worked examples; geometry optimization and atomic
// charges (not in the worked examples) come last.
const PAGES: { id: PageId; href: string; labelKey: string }[] = [
  { id: 'calc',        href: './',               labelKey: 'nav.calc' },
  { id: 'convergence', href: './convergence.html', labelKey: 'nav.convergence' },
  { id: 'walsh',       href: './walsh.html',      labelKey: 'nav.walsh' },
  { id: 'freqanalysis', href: './freqanalysis.html', labelKey: 'nav.freq' },
  { id: 'optimize',    href: './optimize.html',   labelKey: 'nav.optimize' },
  { id: 'accuracy',    href: './accuracy.html',   labelKey: 'nav.accuracy' },
  { id: 'geomopt',     href: './geomopt.html',    labelKey: 'nav.geomopt' },
  { id: 'charges',     href: './charges.html',    labelKey: 'nav.charges' },
];

/**
 * Returns HTML string for the shared header + tab bar.
 * Each page calls this in its render() and inserts into the top of the page.
 */
export function renderHeader(currentPage: PageId | string): string {
  const langLabel = getLang() === 'en' ? 'JA' : 'EN';
  const themeIcon = isDark() ? '\u2600' : '\u263E';

  const tabs = PAGES.map(p => {
    const cls = p.id === currentPage ? 'active' : '';
    return `<a href="${p.href}" class="${cls}">${t(p.labelKey)}</a>`;
  }).join('');

  return `
    <header class="gansu-header">
      <div class="gansu-header-top">
        <span class="gansu-header-title">GANSU Lite</span>
        <div class="gansu-header-right">
          <button id="nav-theme" title="${t('header.themeToggle')}">${themeIcon}</button>
          <button id="nav-lang" title="${t('header.langToggle')}">${langLabel}</button>
        </div>
      </div>
      <nav class="gansu-nav">${tabs}</nav>
    </header>`;
}
