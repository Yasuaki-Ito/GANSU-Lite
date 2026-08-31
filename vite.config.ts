import { defineConfig } from 'vite'

// The built site is committed to docs/ and served by GitHub Pages
// (Settings -> Pages -> Deploy from a branch -> main / docs).
// Public URL: https://yasuaki-ito.github.io/GANSU-Lite/
export default defineConfig(({ command, isPreview }) => ({
  base: process.env.VITE_BASE || (command === 'build' || isPreview ? '/GANSU-Lite/' : '/'),
  build: {
    outDir: 'docs',
    rollupOptions: {
      input: {
        main: 'index.html',
        optimize: 'optimize.html',
        geomopt: 'geomopt.html',
        convergence: 'convergence.html',
        walsh: 'walsh.html',
        accuracy: 'accuracy.html',
        charges: 'charges.html',
        freqanalysis: 'freqanalysis.html',
        crossbench: 'crossbench.html',
      },
    },
  },
  server: {
    strictPort: true,
  },
}))
