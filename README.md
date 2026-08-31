# GANSU Lite — Quantum Chemistry in Your Browser

A client-side, zero-install quantum-chemistry playground that runs Hartree–Fock, DFT, post-HF, excited-state, gradient and Hessian calculations directly in your browser. Built with TypeScript + Vite, accelerated by WebAssembly (optionally SIMD).

**Live demo:** https://yasuaki-ito.github.io/GANSU-Lite/

No server, no upload, no install — all computations run locally in JavaScript / WebAssembly.

## Features

- **SCF:** RHF / UHF / ROHF, RKS / UKS / ROKS
- **DFT functionals:**
  - LDA: SVWN
  - GGA: BLYP, PBE
  - meta-GGA: TPSS (full PKZB correlation)
  - Hybrid: B3LYP, PBE0
  - Range-separated hybrid: CAM-B3LYP, ωB97X-D (RSH-lite, SR-LDA + SR-B88-lite)
- **Dispersion:** Grimme D2, D3(BJ)-lite (Becke-Johnson rational damping with C8 term)
- **Post-HF:** MP2, MP3, CCSD (R / U / RO variants)
- **Excited states:** CIS, ADC(2), TDDFT-TDA, Full Casida (B-coupled RPA)
- **Derivatives:**
  - Analytic nuclear gradient — RHF / RKS (DFT)
  - Numerical Hessian via analytic gradient + TR-mode projection
  - Vibrational frequencies, IR intensities, thermochemistry
- **Geometry optimisation:** SD, CG (FR/PR/HS/DY), BFGS, DFP, SR1, GDIIS
- **Properties:** Mulliken / Löwdin charges, Wiberg bond orders, dipole moment, ⟨S²⟩, energy decomposition, Molden export
- **Basis sets:** STO-3G, 3-21G, 6-31G, cc-pVDZ, aug-cc-pVDZ, def2-SVP, def2-TZVP
- **RI-J:** automatic for pure DFT (auto-generated or optimised cc-pVxZ-RIFIT auxiliary basis)
- **Performance:** WebAssembly + SIMD acceleration, Web Worker for non-blocking UI

## Pages

Every page provides a **Theory selector** to switch between HF and the DFT functionals listed above.

| Page | Description |
|------|-------------|
| **Calculator** (`index.html`) | Full HF/DFT/post-HF analysis on any XYZ molecule |
| **PES Scan** (`optimize.html`) | 1-D potential-energy surface scans (stretch, bend, dihedral); RHF/UHF + DFT |
| **Walsh** (`walsh.html`) | Orbital energies vs. bend angle (e.g. why H₂O bends to 104.5°) |
| **Accuracy** (`accuracy.html`) | Compare HF + MP2 + MP3 + CCSD + DFT side-by-side |
| **Charges** (`charges.html`) | Mulliken/Löwdin atomic charges and dipole on the 3D structure |
| **Basis Set** (`convergence.html`) | Energy convergence vs. basis-set size (STO-3G → def2-TZVP) |
| **Geometry Opt.** (`geomopt.html`) | Interactive optimisation with 3D force arrows; HF + DFT |
| **Vibrations** (`freqanalysis.html`) | Frequencies, IR spectrum, thermochemistry; HF + DFT |

### Heavy-combination guards
- Iterative DFT pages (Geom Opt, Freq, PES Scan, Walsh, Basis Set) hide TPSS — its FD-based V_xc is too slow for many SCF runs.
- Basis Set page auto-skips def2-TZVP for DFT (~1 min per basis × molecule).
- Freq Analysis disables ≥4-atom scenarios for DFT (6N grad evals × DFT × big molecule = several minutes).
- PES Scan warns before starting DFT × UHF × ≥15 points.

## Validation

PySCF cross-check (grid level 5, 17 cases):

| Functional family | ΔE (vs PySCF) | ΔHOMO (vs PySCF) |
|-------------------|--------------:|-----------------:|
| LDA / GGA / standard hybrid | < 0.2 mH | < 0.05 eV |
| TPSS (full PKZB) | 17–60 mH | < 0.2 eV |
| CAM-B3LYP, ωB97X-D (RSH-lite) | 80–150 mH | 0.6–1.0 eV (expected from simplification) |

## Performance

Hot paths are compiled to **WebAssembly (Rust)** with optional **SIMD** (`f64x2`). The runtime auto-detects SIMD support and loads the appropriate binary.

| Kernel | JS | WASM | WASM+SIMD |
|--------|:--:|:----:|:---------:|
| ERI (4-index 2-electron integrals) | ✔ | ✔ | ✔ |
| Fock matrix (RHF/UHF) | ✔ | ✔ | ✔ |
| MP2 / MP3 / CCSD (R/U/RO) | ✔ | ✔ | ✔ |
| RI-J / RI-K / RI-MP2 | ✔ | ✔ | ✔ |
| Nuclear gradient (2-electron + V_xc) | ✔ | ✔ | ✔ |
| Analytical Hessian (2-electron, RHF) | ✔ | ✔ | ✔ |
| CPHF (MO-ERI transform + CG solver) | ✔ | ✔ | ✔ |
| XC numerical integration (Becke grid) | ✔ | ✔ | ✔ |

The SCF pipeline (including post-HF) runs inside a **Web Worker** so the UI thread stays responsive.

## Privacy

GANSU Lite is fully client-side. **No data is sent anywhere** — no analytics, no server-side computation, no tracking. Open the browser DevTools network tab to verify: every calculation runs locally.

## Repository layout

This repository contains both the full source and the deployed build.

```
.
├── index.html, optimize.html, walsh.html, ...   Vite entry points (one per page)
├── src/                    TypeScript sources
│   ├── core/               SCF driver, integrals, DFT (grid + functionals),
│   │                       post-HF (MP2/MP3/CCSD), gradients, Hessian, workers
│   ├── linalg/             eigensolver / linear algebra
│   ├── data/               element data, basis-set tables
│   └── ui/                 page controllers, 3D viewer, charts, i18n, styles
├── wasm-eri/               Rust crate compiled to WebAssembly
│   ├── Cargo.toml / Cargo.lock
│   ├── build-wasm.sh       builds both the plain and the SIMD binary
│   └── src/                ERI (MD/OS/Rys), Fock, MP2/MP3/CCSD, RI,
│                           gradient, Hessian, CPHF, XC integration, Boys
├── public/                 static assets copied verbatim into the build
│   ├── wasm/               pre-built .wasm binaries (committed)
│   ├── basis/              basis-set files (.gbs)
│   ├── xyz/                sample geometries
│   ├── shaders/            WebGPU compute shaders
│   └── tests/              reference values for validation
├── docs/                   >>> BUILD OUTPUT — served by GitHub Pages <<<
├── package.json / package-lock.json
├── tsconfig.json
└── vite.config.ts
```

`docs/` is generated by `npm run build` and is committed only so that GitHub Pages
can serve it (Settings → Pages → *Deploy from a branch* → `main` / `/docs`).
**Do not edit `docs/` by hand** — edit `src/` and rebuild.

## Building from source

Requirements: **Node.js 20+** (tested on 24) and npm.
A Rust toolchain is needed *only* if you want to rebuild the WebAssembly kernels —
the pre-built `.wasm` binaries are committed under `public/wasm/`.

### 1. Web application

```bash
npm ci            # install exact dependency versions from package-lock.json
npm run dev       # dev server (http://localhost:5173/)
npm run build     # type-check (tsc -b) + bundle (vite build) -> docs/
npm run preview   # serve docs/ locally at http://localhost:4173/GANSU-Lite/
```

`vite.config.ts` sets `base` to `/GANSU-Lite/` for production builds (matching the
GitHub Pages URL) and `/` for the dev server. Override with the `VITE_BASE`
environment variable if you deploy under a different path:

```bash
VITE_BASE=/ npm run build
```

### 2. WebAssembly kernels (optional)

```bash
rustup target add wasm32-unknown-unknown
cargo install wasm-pack

cd wasm-eri
./build-wasm.sh   # -> ../public/wasm/wasm_eri_bg.wasm  (baseline)
                  # -> ../public/wasm/wasm_eri_simd_bg.wasm  (-C target-feature=+simd128)
```

The runtime feature-detects SIMD support and loads the appropriate binary.
After rebuilding the wasm, re-run `npm run build` so the new binaries are copied
into `docs/`.

### 3. Deployment

`docs/` is the published site. After `npm run build`, commit the regenerated
`docs/` together with your source changes and push to `main`; GitHub Pages picks
it up automatically.

## License

BSD 3-Clause License. See [LICENSE](LICENSE).

## Acknowledgements

- Reference values cross-checked against [PySCF](https://pyscf.org/).
- Auxiliary basis sets (cc-pVxZ-RIFIT) derived from PySCF's bundled Weigend / Hattig basis tables.
- Boys function tables derived from standard quantum-chemistry literature.

## Companion textbook

GANSU Lite is used in the exercises of a companion textbook: [日本語](https://yasuaki-ito.github.io/book/qcbook/) · [English](https://yasuaki-ito.github.io/book/en/qcbook/)
