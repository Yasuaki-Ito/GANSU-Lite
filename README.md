# GANSU Lite — Quantum Chemistry in Your Browser

A client-side, zero-install quantum-chemistry playground that runs Hartree–Fock, DFT, post-HF, excited-state, gradient and Hessian calculations directly in your browser. Built with TypeScript + Vite, accelerated by WebAssembly (optionally SIMD).

**Live demo:** https://yasuaki-ito.github.io/GANSU-Lite/

**Companion textbook:** [日本語](https://yasuaki-ito.github.io/book/qcbook/) · [English](https://yasuaki-ito.github.io/book/en/qcbook/)

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
- **Basis sets:** STO-3G, STO-6G, 3-21G, 6-31G, 6-31G(d,p), cc-pVDZ, cc-pVTZ, aug-cc-pVDZ, def2-SVP, def2-TZVP
- **RI-J:** automatic for pure DFT (auto-generated or optimised cc-pVxZ-RIFIT auxiliary basis)
- **Performance:** WebAssembly + SIMD acceleration, Web Worker for non-blocking UI

## Pages

Every page (except DFT functional ladder) provides a **Theory selector** to switch between HF and the DFT functionals listed above.

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
| **DFT** (`dft.html`) | Functional ladder demo (LDA → GGA → meta-GGA → hybrid → RSH) |

### Heavy-combination guards
- Iterative DFT pages (Geom Opt, Freq, PES Scan, Walsh, Basis Set) hide TPSS — its FD-based V_xc is too slow for many SCF runs.
- Basis Set page auto-skips cc-pVTZ / def2-TZVP for DFT (~1 min per basis × molecule).
- Freq Analysis disables ≥4-atom scenarios for DFT (6N grad evals × DFT × big molecule = several minutes).
- PES Scan warns before starting DFT × UHF × ≥15 points.

## Validation

PySCF cross-check (grid level 5, 17 cases):

| Functional family | |E_GANSU − E_PySCF| | |HOMO_GANSU − HOMO_PySCF| |
|-------------------|-----:|-----:|
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

## License

BSD 3-Clause License. See [LICENSE](LICENSE).

## Acknowledgements

- Reference values cross-checked against [PySCF](https://pyscf.org/).
- Auxiliary basis sets (cc-pVxZ-RIFIT) derived from PySCF's bundled Weigend / Hattig basis tables.
- Boys function tables derived from standard quantum-chemistry literature.
