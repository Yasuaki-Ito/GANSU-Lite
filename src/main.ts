import { initTheme } from './ui/theme';
import { initApp } from './ui/app';

initTheme();

const root = document.getElementById('app');
if (root) {
  initApp(root);
}
