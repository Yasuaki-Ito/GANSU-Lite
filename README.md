# GANSU Lite — Quantum Chemistry in Your Browser

A client-side, zero-install quantum-chemistry playground that runs Hartree–Fock, DFT, post-HF, and property calculations directly in your browser. Built with TypeScript + Vite, accelerated by WebAssembly (optionally SIMD).

**Live demo:** https://yasuaki-ito.github.io/GANSU-Lite/

No server, no upload, no install — all computations run locally in JavaScript/WASM.

## Features

- **SCF:** RHF / UHF / ROHF / RKS / UKS
- **DFT:** SVWN (LDA), BLYP, PBE (GGA), B3LYP (hybrid)
- **Post-HF:** MP2, MP3, CCSD (R/U/RO variants)
- **Excited states:** CIS (singlet/triplet), ADC(2)
- **Derivatives:** Analytic gradient, numerical/analytical Hessian, vibrational frequencies + IR spectrum
- **Geometry optimisation:** SD, CG, BFGS, DFP, SR1, GDIIS
- **Properties:** Mulliken/Löwdin charges, Wiberg bond orders, dipole moment, Molden export
- **Basis sets:** STO-3G, STO-6G, 3-21G, 6-31G, 6-31G(d,p), cc-pVDZ, cc-pVTZ
- **Performance:** WebAssembly + SIMD acceleration, Web Worker for non-blocking UI

## Pages

| Page | Description |
|------|-------------|
| **Calculator** | Free-form HF/DFT/post-HF calculations with 3D molecule viewer |
| **PES Scan** | Potential energy surface scans (stretch, bend, dihedral) |
| **Walsh** | Walsh diagrams — orbital energies vs. geometry |
| **Accuracy** | Side-by-side comparison of 8 methods (HF/MP2/MP3/CCSD/SVWN/BLYP/PBE/B3LYP) |
| **Charges** | Atomic charges and dipole visualisation |
| **Basis Set** | Basis-set convergence demo |
| **Geometry Opt.** | Interactive geometry optimisation with 3D force arrows |
| **Vibrations** | Vibrational frequencies, IR spectrum, thermochemistry |

## License

BSD 3-Clause License. See [LICENSE](LICENSE).
