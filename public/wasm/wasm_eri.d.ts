/* tslint:disable */
/* eslint-disable */

/**
 * Compute 2-electron gradient contribution.
 * shells_flat: Float64Array with 8 fields per shell (ext format with atomIndex)
 * density: n×n density matrix (row-major)
 * norms: CGTO normalization factors
 * num_basis: number of basis functions
 * num_atoms: number of atoms
 * Returns: Float64Array of gradient (3*num_atoms)
 */
export function compute_2e_gradient_wasm(shells_flat: Float64Array, density: Float64Array, norms: Float64Array, num_basis: number, num_atoms: number): Float64Array;

/**
 * Compute 2-electron skeleton Hessian.
 * shells_flat: Float64Array with 8 fields per shell (ext format with atomIndex)
 * density: n×n density matrix (row-major)
 * norms: CGTO normalization factors
 * Returns: Hessian (3*num_atoms)² row-major
 */
export function compute_2e_hessian_wasm(shells_flat: Float64Array, density: Float64Array, norms: Float64Array, num_basis: number, num_atoms: number): Float64Array;

/**
 * Compute CCSD correlation energy (RHF spin-orbital formulation).
 * Returns f64 scalar (correlation energy).
 */
export function compute_ccsd_wasm(eri: Float64Array, coeff: Float64Array, epsilon: Float64Array, nocc: number, num_basis: number): number;

/**
 * Compute all ERIs.
 *
 * `shells_flat`: Float64Array with 7 fields per shell (see FIELDS_PER_SHELL)
 * `norm_factors`: Float64Array of CGTO normalization factors
 * `num_basis`: number of basis functions
 * `schwarz_threshold`: screening threshold (e.g. 1e-10)
 *
 * Returns Float64Array of ERI values in 8-fold symmetric storage.
 */
export function compute_eris_wasm(shells_flat: Float64Array, norm_factors: Float64Array, num_basis: number, schwarz_threshold: number): Float64Array;

/**
 * Compute RHF Fock matrix.
 * `eri`: 8-fold symmetric ERI data
 * `density`: P matrix (N*N, row-major)
 * `core_h`: H matrix (N*N, row-major)
 * `num_basis`: number of basis functions
 * Returns F matrix (N*N, row-major)
 */
export function compute_fock_rhf(eri: Float64Array, density: Float64Array, core_h: Float64Array, num_basis: number): Float64Array;

/**
 * Compute UHF Fock matrices (alpha and beta).
 * Returns concatenated [Fa..., Fb...] (2*N*N elements)
 */
export function compute_fock_uhf(eri: Float64Array, density_alpha: Float64Array, density_beta: Float64Array, density_total: Float64Array, core_h: Float64Array, num_basis: number): Float64Array;

/**
 * Compute MP2 correlation energy.
 * `eri`: 8-fold symmetric ERI data
 * `coeff`: MO coefficient matrix (N*N, row-major)
 * `epsilon`: orbital energies (N elements)
 * `nocc`: number of occupied spatial orbitals
 * `num_basis`: number of basis functions
 * Returns MP2 correlation energy (f64 scalar)
 */
export function compute_mp2_wasm(eri: Float64Array, coeff: Float64Array, epsilon: Float64Array, nocc: number, num_basis: number): number;

/**
 * Compute MP3 correlation energy (RHF spatial-orbital formulas).
 * Returns Vec<f64> with 2 elements: [mp2_corr, mp3_corr]
 */
export function compute_mp3_wasm(eri: Float64Array, coeff: Float64Array, epsilon: Float64Array, nocc: number, num_basis: number): Float64Array;

/**
 * Compute RI Coulomb matrix J.
 * `b_matrix`: B[P*N*N + mu*N + nu], length = naux * N * N
 * `density`: P matrix (N*N, row-major)
 * Returns J matrix (N*N, row-major)
 */
export function compute_ri_coulomb_wasm(b_matrix: Float64Array, density: Float64Array, naux: number, num_basis: number): Float64Array;

/**
 * Compute RI Exchange matrix K.
 * `coefficients`: MO coefficient matrix (N*N, row-major)
 * `nocc`: number of occupied orbitals
 * Returns K matrix (N*N, row-major)
 */
export function compute_ri_exchange_wasm(b_matrix: Float64Array, coefficients: Float64Array, naux: number, num_basis: number, nocc: number): Float64Array;

/**
 * Compute RI-RHF Fock matrix: F = H + J - K
 */
export function compute_ri_fock_rhf_wasm(b_matrix: Float64Array, density: Float64Array, core_h: Float64Array, coefficients: Float64Array, naux: number, num_basis: number, nocc: number): Float64Array;

/**
 * Compute RI-UHF Fock matrices: Fα = H + J(Ptot) - Kα, Fβ = H + J(Ptot) - Kβ
 * Returns concatenated [Fα..., Fβ...] (2*N*N elements)
 */
export function compute_ri_fock_uhf_wasm(b_matrix: Float64Array, density_total: Float64Array, core_h: Float64Array, coeff_alpha: Float64Array, nocc_alpha: number, coeff_beta: Float64Array, nocc_beta: number, naux: number, num_basis: number): Float64Array;

/**
 * Compute RI-MP2 correlation energy (RHF).
 */
export function compute_ri_mp2_wasm(b_matrix: Float64Array, coefficients: Float64Array, epsilon: Float64Array, naux: number, nocc: number, num_basis: number): number;

/**
 * Compute RI Setup: 2c/3c integrals, Cholesky, B matrix.
 * Returns Vec<f64> with layout [naux, B_matrix...]
 */
export function compute_ri_setup_wasm(primary_shells_flat: Float64Array, primary_norm: Float64Array, nbasis: number, aux_shells_flat: Float64Array, aux_norm: Float64Array, naux: number): Float64Array;

/**
 * Compute RI-UMP2 correlation energy.
 */
export function compute_ri_ump2_wasm(b_matrix: Float64Array, ca: Float64Array, cb: Float64Array, epsilon_a: Float64Array, epsilon_b: Float64Array, naux: number, nocc_a: number, nocc_b: number, num_basis: number): number;

/**
 * Compute UCCSD correlation energy from converged UHF results.
 * Returns f64 scalar (correlation energy).
 */
export function compute_uccsd_wasm(eri: Float64Array, ca: Float64Array, cb: Float64Array, epsilon_a: Float64Array, epsilon_b: Float64Array, nocc_a: number, nocc_b: number, num_basis: number): number;

/**
 * Compute UMP2 correlation energy from converged UHF results.
 */
export function compute_ump2_wasm(eri: Float64Array, ca: Float64Array, cb: Float64Array, epsilon_a: Float64Array, epsilon_b: Float64Array, nocc_a: number, nocc_b: number, num_basis: number): number;

/**
 * Compute UMP3 correlation energy from converged UHF results.
 * Returns Vec<f64> with 2 elements: [mp2_corr, mp3_corr]
 */
export function compute_ump3_wasm(eri: Float64Array, ca: Float64Array, cb: Float64Array, epsilon_a: Float64Array, epsilon_b: Float64Array, nocc_a: number, nocc_b: number, num_basis: number): Float64Array;

/**
 * Compute XC energy and V_xc matrix on a numerical grid.
 *
 * Returns Vec<f64> with layout: [exc, num_electrons, vxc_a..., vxc_b...]
 * `grid_flat`: [x0, y0, z0, w0, x1, y1, z1, w1, ...] (4 values per point)
 * `density_b`: empty for RKS (restricted); N*N for UKS (unrestricted)
 * `func_id`: 0=SVWN, 1=BLYP, 2=PBE, 3=B3LYP
 * `need_grad`: true for GGA/hybrid functionals
 */
export function compute_xc_wasm(shells_flat: Float64Array, norm_factors: Float64Array, density_a: Float64Array, density_b: Float64Array, grid_flat: Float64Array, num_basis: number, func_id: number, need_grad: boolean): Float64Array;

/**
 * Solve CPHF equations via preconditioned CG.
 * mo_eri: MO-basis ERIs (n⁴)
 * eps: orbital energies (nmo)
 * rhs: right-hand side [n_pert × dim], dim = nocc*nvir
 * Returns: solution U^x [n_pert × dim]
 */
export function solve_cphf_wasm(mo_eri: Float64Array, eps: Float64Array, rhs: Float64Array, nocc: number, nvir: number, nmo: number, n_pert: number, tol: number, max_iter: number): Float64Array;

/**
 * AO → MO integral transform (4-quarter, O(N⁵)).
 * c_flat: MO coefficient matrix (n×n, row-major)
 * eri_sym: AO ERIs in 8-fold symmetric packed storage
 * num_basis: basis set size
 * Returns: full (pq|rs) as n⁴ flat array
 */
export function transform_eri_to_mo_wasm(c_flat: Float64Array, eri_sym: Float64Array, num_basis: number): Float64Array;

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly compute_2e_gradient_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => [number, number];
    readonly compute_2e_hessian_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => [number, number];
    readonly compute_ccsd_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => number;
    readonly compute_eris_wasm: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number];
    readonly compute_fock_rhf: (a: number, b: number, c: number, d: number, e: number, f: number, g: number) => [number, number];
    readonly compute_fock_uhf: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number, k: number) => [number, number];
    readonly compute_mp2_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => number;
    readonly compute_mp3_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => [number, number];
    readonly compute_ri_coulomb_wasm: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number];
    readonly compute_ri_exchange_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number) => [number, number];
    readonly compute_ri_fock_rhf_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number, k: number) => [number, number];
    readonly compute_ri_fock_uhf_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number, k: number, l: number, m: number, n: number) => [number, number];
    readonly compute_ri_mp2_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number) => number;
    readonly compute_ri_setup_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number) => [number, number];
    readonly compute_ri_ump2_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number, k: number, l: number, m: number, n: number) => number;
    readonly compute_uccsd_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number, k: number, l: number, m: number) => number;
    readonly compute_ump2_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number, k: number, l: number, m: number) => number;
    readonly compute_ump3_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number, k: number, l: number, m: number) => [number, number];
    readonly compute_xc_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number, k: number, l: number, m: number) => [number, number];
    readonly solve_cphf_wasm: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number, k: number, l: number) => [number, number];
    readonly transform_eri_to_mo_wasm: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
