mod boys;
mod boys_grid;
mod fock;
mod cphf;
mod grad;
mod hessian;
mod integrals;
mod integrals_os;
mod integrals_rys;
mod ccsd;
mod mp2;
mod mp3;
mod ri;
mod simd_utils;
mod xc;

use wasm_bindgen::prelude::*;
use integrals::PrimShell;

/// Flat-packed primitive shell data from JS.
/// Layout per shell: [exponent, coefficient, x, y, z, shellType, basisIndex]
/// Total length = numShells * 7
const FIELDS_PER_SHELL: usize = 7;

/// Extended layout with atomIndex: [exponent, coefficient, x, y, z, shellType, basisIndex, atomIndex]
const FIELDS_PER_SHELL_EXT: usize = 8;

fn unpack_shells(flat: &[f64]) -> Vec<PrimShell> {
    let n = flat.len() / FIELDS_PER_SHELL;
    let mut shells = Vec::with_capacity(n);
    for i in 0..n {
        let off = i * FIELDS_PER_SHELL;
        shells.push(PrimShell {
            exponent: flat[off],
            coefficient: flat[off + 1],
            x: flat[off + 2],
            y: flat[off + 3],
            z: flat[off + 4],
            shell_type: flat[off + 5] as u32,
            basis_index: flat[off + 6] as u32,
            atom_index: 0,
        });
    }
    shells
}

fn unpack_shells_ext(flat: &[f64]) -> Vec<PrimShell> {
    let n = flat.len() / FIELDS_PER_SHELL_EXT;
    let mut shells = Vec::with_capacity(n);
    for i in 0..n {
        let off = i * FIELDS_PER_SHELL_EXT;
        shells.push(PrimShell {
            exponent: flat[off],
            coefficient: flat[off + 1],
            x: flat[off + 2],
            y: flat[off + 3],
            z: flat[off + 4],
            shell_type: flat[off + 5] as u32,
            basis_index: flat[off + 6] as u32,
            atom_index: flat[off + 7] as u32,
        });
    }
    shells
}

/// Compute all ERIs.
///
/// `shells_flat`: Float64Array with 7 fields per shell (see FIELDS_PER_SHELL)
/// `norm_factors`: Float64Array of CGTO normalization factors
/// `num_basis`: number of basis functions
/// `schwarz_threshold`: screening threshold (e.g. 1e-10)
///
/// Returns Float64Array of ERI values in 8-fold symmetric storage.
#[wasm_bindgen]
pub fn compute_eris_wasm(
    shells_flat: &[f64],
    norm_factors: &[f64],
    num_basis: usize,
    schwarz_threshold: f64,
) -> Vec<f64> {
    let shells = unpack_shells(shells_flat);
    integrals::compute_eris(&shells, norm_factors, num_basis, schwarz_threshold)
}

/// Compute RHF Fock matrix.
/// `eri`: 8-fold symmetric ERI data
/// `density`: P matrix (N*N, row-major)
/// `core_h`: H matrix (N*N, row-major)
/// `num_basis`: number of basis functions
/// Returns F matrix (N*N, row-major)
#[wasm_bindgen]
pub fn compute_fock_rhf(
    eri: &[f64],
    density: &[f64],
    core_h: &[f64],
    num_basis: usize,
) -> Vec<f64> {
    fock::compute_fock_rhf(eri, density, core_h, num_basis)
}

/// Compute UHF Fock matrices (alpha and beta).
/// Returns concatenated [Fa..., Fb...] (2*N*N elements)
#[wasm_bindgen]
pub fn compute_fock_uhf(
    eri: &[f64],
    density_alpha: &[f64],
    density_beta: &[f64],
    density_total: &[f64],
    core_h: &[f64],
    num_basis: usize,
) -> Vec<f64> {
    fock::compute_fock_uhf(eri, density_alpha, density_beta, density_total, core_h, num_basis)
}

/// Compute MP2 correlation energy.
/// `eri`: 8-fold symmetric ERI data
/// `coeff`: MO coefficient matrix (N*N, row-major)
/// `epsilon`: orbital energies (N elements)
/// `nocc`: number of occupied spatial orbitals
/// `num_basis`: number of basis functions
/// Returns MP2 correlation energy (f64 scalar)
#[wasm_bindgen]
pub fn compute_mp2_wasm(
    eri: &[f64],
    coeff: &[f64],
    epsilon: &[f64],
    nocc: usize,
    num_basis: usize,
) -> f64 {
    mp2::compute_mp2(eri, coeff, epsilon, nocc, num_basis)
}

/// Compute MP3 correlation energy (RHF spatial-orbital formulas).
/// Returns Vec<f64> with 2 elements: [mp2_corr, mp3_corr]
#[wasm_bindgen]
pub fn compute_mp3_wasm(
    eri: &[f64],
    coeff: &[f64],
    epsilon: &[f64],
    nocc: usize,
    num_basis: usize,
) -> Vec<f64> {
    mp3::compute_mp3(eri, coeff, epsilon, nocc, num_basis)
}

/// Compute CCSD correlation energy (RHF spin-orbital formulation).
/// Returns f64 scalar (correlation energy).
#[wasm_bindgen]
pub fn compute_ccsd_wasm(
    eri: &[f64],
    coeff: &[f64],
    epsilon: &[f64],
    nocc: usize,
    num_basis: usize,
) -> f64 {
    ccsd::compute_ccsd(eri, coeff, epsilon, nocc, num_basis)
}

// ---------------------------------------------------------------------------
// UHF post-HF entry points
// ---------------------------------------------------------------------------

/// Compute UMP2 correlation energy from converged UHF results.
#[wasm_bindgen]
pub fn compute_ump2_wasm(
    eri: &[f64],
    ca: &[f64],
    cb: &[f64],
    epsilon_a: &[f64],
    epsilon_b: &[f64],
    nocc_a: usize,
    nocc_b: usize,
    num_basis: usize,
) -> f64 {
    mp2::compute_ump2(eri, ca, cb, epsilon_a, epsilon_b, nocc_a, nocc_b, num_basis)
}

/// Compute UMP3 correlation energy from converged UHF results.
/// Returns Vec<f64> with 2 elements: [mp2_corr, mp3_corr]
#[wasm_bindgen]
pub fn compute_ump3_wasm(
    eri: &[f64],
    ca: &[f64],
    cb: &[f64],
    epsilon_a: &[f64],
    epsilon_b: &[f64],
    nocc_a: usize,
    nocc_b: usize,
    num_basis: usize,
) -> Vec<f64> {
    mp3::compute_ump3(eri, ca, cb, epsilon_a, epsilon_b, nocc_a, nocc_b, num_basis)
}

/// Compute UCCSD correlation energy from converged UHF results.
/// Returns f64 scalar (correlation energy).
#[wasm_bindgen]
pub fn compute_uccsd_wasm(
    eri: &[f64],
    ca: &[f64],
    cb: &[f64],
    epsilon_a: &[f64],
    epsilon_b: &[f64],
    nocc_a: usize,
    nocc_b: usize,
    num_basis: usize,
) -> f64 {
    ccsd::compute_uccsd(eri, ca, cb, epsilon_a, epsilon_b, nocc_a, nocc_b, num_basis)
}

// ---------------------------------------------------------------------------
// RI (Resolution of Identity) entry points
// ---------------------------------------------------------------------------

/// Compute RI Setup: 2c/3c integrals, Cholesky, B matrix.
/// Returns Vec<f64> with layout [naux, B_matrix...]
#[wasm_bindgen]
pub fn compute_ri_setup_wasm(
    primary_shells_flat: &[f64],
    primary_norm: &[f64],
    nbasis: usize,
    aux_shells_flat: &[f64],
    aux_norm: &[f64],
    naux: usize,
) -> Vec<f64> {
    let primary_shells = unpack_shells(primary_shells_flat);
    let aux_shells = unpack_shells(aux_shells_flat);
    ri::compute_ri_setup(&primary_shells, primary_norm, nbasis, &aux_shells, aux_norm, naux)
}

/// Compute RI Coulomb matrix J.
/// `b_matrix`: B[P*N*N + mu*N + nu], length = naux * N * N
/// `density`: P matrix (N*N, row-major)
/// Returns J matrix (N*N, row-major)
#[wasm_bindgen]
pub fn compute_ri_coulomb_wasm(
    b_matrix: &[f64],
    density: &[f64],
    naux: usize,
    num_basis: usize,
) -> Vec<f64> {
    ri::compute_ri_coulomb(b_matrix, density, naux, num_basis)
}

/// Compute RI Exchange matrix K.
/// `coefficients`: MO coefficient matrix (N*N, row-major)
/// `nocc`: number of occupied orbitals
/// Returns K matrix (N*N, row-major)
#[wasm_bindgen]
pub fn compute_ri_exchange_wasm(
    b_matrix: &[f64],
    coefficients: &[f64],
    naux: usize,
    num_basis: usize,
    nocc: usize,
) -> Vec<f64> {
    ri::compute_ri_exchange(b_matrix, coefficients, naux, num_basis, nocc)
}

/// Compute RI-RHF Fock matrix: F = H + J - K
#[wasm_bindgen]
pub fn compute_ri_fock_rhf_wasm(
    b_matrix: &[f64],
    density: &[f64],
    core_h: &[f64],
    coefficients: &[f64],
    naux: usize,
    num_basis: usize,
    nocc: usize,
) -> Vec<f64> {
    ri::compute_ri_fock_rhf(b_matrix, density, core_h, coefficients, naux, num_basis, nocc)
}

/// Compute RI-UHF Fock matrices: Fα = H + J(Ptot) - Kα, Fβ = H + J(Ptot) - Kβ
/// Returns concatenated [Fα..., Fβ...] (2*N*N elements)
#[wasm_bindgen]
pub fn compute_ri_fock_uhf_wasm(
    b_matrix: &[f64],
    density_total: &[f64],
    core_h: &[f64],
    coeff_alpha: &[f64],
    nocc_alpha: usize,
    coeff_beta: &[f64],
    nocc_beta: usize,
    naux: usize,
    num_basis: usize,
) -> Vec<f64> {
    ri::compute_ri_fock_uhf(b_matrix, density_total, core_h, coeff_alpha, nocc_alpha, coeff_beta, nocc_beta, naux, num_basis)
}

/// Compute RI-MP2 correlation energy (RHF).
#[wasm_bindgen]
pub fn compute_ri_mp2_wasm(
    b_matrix: &[f64],
    coefficients: &[f64],
    epsilon: &[f64],
    naux: usize,
    nocc: usize,
    num_basis: usize,
) -> f64 {
    ri::compute_ri_mp2(b_matrix, coefficients, epsilon, naux, nocc, num_basis)
}

/// Compute RI-UMP2 correlation energy.
#[wasm_bindgen]
pub fn compute_ri_ump2_wasm(
    b_matrix: &[f64],
    ca: &[f64],
    cb: &[f64],
    epsilon_a: &[f64],
    epsilon_b: &[f64],
    naux: usize,
    nocc_a: usize,
    nocc_b: usize,
    num_basis: usize,
) -> f64 {
    ri::compute_ri_ump2(b_matrix, ca, cb, epsilon_a, epsilon_b, naux, nocc_a, nocc_b, num_basis)
}

// ---------------------------------------------------------------------------
// XC integration
// ---------------------------------------------------------------------------

/// Compute XC energy and V_xc matrix on a numerical grid.
///
/// Returns Vec<f64> with layout: [exc, num_electrons, vxc_a..., vxc_b...]
/// `grid_flat`: [x0, y0, z0, w0, x1, y1, z1, w1, ...] (4 values per point)
/// `density_b`: empty for RKS (restricted); N*N for UKS (unrestricted)
/// `func_id`: 0=SVWN, 1=BLYP, 2=PBE, 3=B3LYP
/// `need_grad`: true for GGA/hybrid functionals
#[wasm_bindgen]
pub fn compute_xc_wasm(
    shells_flat: &[f64],
    norm_factors: &[f64],
    density_a: &[f64],
    density_b: &[f64],
    grid_flat: &[f64],
    num_basis: usize,
    func_id: u32,
    need_grad: bool,
) -> Vec<f64> {
    let shells = unpack_shells(shells_flat);
    xc::compute_xc(&shells, norm_factors, density_a, density_b, grid_flat, num_basis, func_id, need_grad)
}

// ---------------------------------------------------------------------------
// Nuclear gradient (2-electron part)
// ---------------------------------------------------------------------------

/// Compute 2-electron gradient contribution.
/// shells_flat: Float64Array with 8 fields per shell (ext format with atomIndex)
/// density: n×n density matrix (row-major)
/// norms: CGTO normalization factors
/// num_basis: number of basis functions
/// num_atoms: number of atoms
/// Returns: Float64Array of gradient (3*num_atoms)
#[wasm_bindgen]
pub fn compute_2e_gradient_wasm(
    shells_flat: &[f64],
    density: &[f64],
    norms: &[f64],
    num_basis: usize,
    num_atoms: usize,
) -> Vec<f64> {
    let shells = unpack_shells_ext(shells_flat);
    grad::compute_2e_gradient(&shells, density, norms, num_basis, num_atoms)
}

// ---------------------------------------------------------------------------
// 2-electron skeleton Hessian
// ---------------------------------------------------------------------------

/// Compute 2-electron skeleton Hessian.
/// shells_flat: Float64Array with 8 fields per shell (ext format with atomIndex)
/// density: n×n density matrix (row-major)
/// norms: CGTO normalization factors
/// Returns: Hessian (3*num_atoms)² row-major
#[wasm_bindgen]
pub fn compute_2e_hessian_wasm(
    shells_flat: &[f64],
    density: &[f64],
    norms: &[f64],
    num_basis: usize,
    num_atoms: usize,
) -> Vec<f64> {
    let shells = unpack_shells_ext(shells_flat);
    hessian::compute_2e_hessian(&shells, density, norms, num_basis, num_atoms)
}

// ---------------------------------------------------------------------------
// CPHF: MO-ERI transform + CG solver
// ---------------------------------------------------------------------------

/// AO → MO integral transform (4-quarter, O(N⁵)).
/// c_flat: MO coefficient matrix (n×n, row-major)
/// eri_sym: AO ERIs in 8-fold symmetric packed storage
/// num_basis: basis set size
/// Returns: full (pq|rs) as n⁴ flat array
#[wasm_bindgen]
pub fn transform_eri_to_mo_wasm(
    c_flat: &[f64],
    eri_sym: &[f64],
    num_basis: usize,
) -> Vec<f64> {
    cphf::transform_eri_to_mo(c_flat, eri_sym, num_basis)
}

/// Solve CPHF equations via preconditioned CG.
/// mo_eri: MO-basis ERIs (n⁴)
/// eps: orbital energies (nmo)
/// rhs: right-hand side [n_pert × dim], dim = nocc*nvir
/// Returns: solution U^x [n_pert × dim]
#[wasm_bindgen]
pub fn solve_cphf_wasm(
    mo_eri: &[f64],
    eps: &[f64],
    rhs: &[f64],
    nocc: usize,
    nvir: usize,
    nmo: usize,
    n_pert: usize,
    tol: f64,
    max_iter: usize,
) -> Vec<f64> {
    cphf::solve_cphf(mo_eri, eps, rhs, nocc, nvir, nmo, n_pert, tol, max_iter)
}
