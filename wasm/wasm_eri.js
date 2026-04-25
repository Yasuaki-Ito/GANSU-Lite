/* @ts-self-types="./wasm_eri.d.ts" */

/**
 * Compute 2-electron gradient contribution.
 * shells_flat: Float64Array with 8 fields per shell (ext format with atomIndex)
 * density: n×n density matrix (row-major)
 * norms: CGTO normalization factors
 * num_basis: number of basis functions
 * num_atoms: number of atoms
 * Returns: Float64Array of gradient (3*num_atoms)
 * @param {Float64Array} shells_flat
 * @param {Float64Array} density
 * @param {Float64Array} norms
 * @param {number} num_basis
 * @param {number} num_atoms
 * @returns {Float64Array}
 */
export function compute_2e_gradient_wasm(shells_flat, density, norms, num_basis, num_atoms) {
    const ptr0 = passArrayF64ToWasm0(shells_flat, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(density, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArrayF64ToWasm0(norms, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ret = wasm.compute_2e_gradient_wasm(ptr0, len0, ptr1, len1, ptr2, len2, num_basis, num_atoms);
    var v4 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
    return v4;
}

/**
 * Compute 2-electron skeleton Hessian.
 * shells_flat: Float64Array with 8 fields per shell (ext format with atomIndex)
 * density: n×n density matrix (row-major)
 * norms: CGTO normalization factors
 * Returns: Hessian (3*num_atoms)² row-major
 * @param {Float64Array} shells_flat
 * @param {Float64Array} density
 * @param {Float64Array} norms
 * @param {number} num_basis
 * @param {number} num_atoms
 * @returns {Float64Array}
 */
export function compute_2e_hessian_wasm(shells_flat, density, norms, num_basis, num_atoms) {
    const ptr0 = passArrayF64ToWasm0(shells_flat, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(density, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArrayF64ToWasm0(norms, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ret = wasm.compute_2e_hessian_wasm(ptr0, len0, ptr1, len1, ptr2, len2, num_basis, num_atoms);
    var v4 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
    return v4;
}

/**
 * Compute CCSD correlation energy (RHF spin-orbital formulation).
 * Returns f64 scalar (correlation energy).
 * @param {Float64Array} eri
 * @param {Float64Array} coeff
 * @param {Float64Array} epsilon
 * @param {number} nocc
 * @param {number} num_basis
 * @returns {number}
 */
export function compute_ccsd_wasm(eri, coeff, epsilon, nocc, num_basis) {
    const ptr0 = passArrayF64ToWasm0(eri, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(coeff, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArrayF64ToWasm0(epsilon, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ret = wasm.compute_ccsd_wasm(ptr0, len0, ptr1, len1, ptr2, len2, nocc, num_basis);
    return ret;
}

/**
 * Compute all ERIs.
 *
 * `shells_flat`: Float64Array with 7 fields per shell (see FIELDS_PER_SHELL)
 * `norm_factors`: Float64Array of CGTO normalization factors
 * `num_basis`: number of basis functions
 * `schwarz_threshold`: screening threshold (e.g. 1e-10)
 *
 * Returns Float64Array of ERI values in 8-fold symmetric storage.
 * @param {Float64Array} shells_flat
 * @param {Float64Array} norm_factors
 * @param {number} num_basis
 * @param {number} schwarz_threshold
 * @returns {Float64Array}
 */
export function compute_eris_wasm(shells_flat, norm_factors, num_basis, schwarz_threshold) {
    const ptr0 = passArrayF64ToWasm0(shells_flat, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(norm_factors, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ret = wasm.compute_eris_wasm(ptr0, len0, ptr1, len1, num_basis, schwarz_threshold);
    var v3 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
    return v3;
}

/**
 * Compute RHF Fock matrix.
 * `eri`: 8-fold symmetric ERI data
 * `density`: P matrix (N*N, row-major)
 * `core_h`: H matrix (N*N, row-major)
 * `num_basis`: number of basis functions
 * Returns F matrix (N*N, row-major)
 * @param {Float64Array} eri
 * @param {Float64Array} density
 * @param {Float64Array} core_h
 * @param {number} num_basis
 * @returns {Float64Array}
 */
export function compute_fock_rhf(eri, density, core_h, num_basis) {
    const ptr0 = passArrayF64ToWasm0(eri, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(density, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArrayF64ToWasm0(core_h, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ret = wasm.compute_fock_rhf(ptr0, len0, ptr1, len1, ptr2, len2, num_basis);
    var v4 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
    return v4;
}

/**
 * Compute UHF Fock matrices (alpha and beta).
 * Returns concatenated [Fa..., Fb...] (2*N*N elements)
 * @param {Float64Array} eri
 * @param {Float64Array} density_alpha
 * @param {Float64Array} density_beta
 * @param {Float64Array} density_total
 * @param {Float64Array} core_h
 * @param {number} num_basis
 * @returns {Float64Array}
 */
export function compute_fock_uhf(eri, density_alpha, density_beta, density_total, core_h, num_basis) {
    const ptr0 = passArrayF64ToWasm0(eri, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(density_alpha, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArrayF64ToWasm0(density_beta, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ptr3 = passArrayF64ToWasm0(density_total, wasm.__wbindgen_malloc);
    const len3 = WASM_VECTOR_LEN;
    const ptr4 = passArrayF64ToWasm0(core_h, wasm.__wbindgen_malloc);
    const len4 = WASM_VECTOR_LEN;
    const ret = wasm.compute_fock_uhf(ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3, ptr4, len4, num_basis);
    var v6 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
    return v6;
}

/**
 * Compute MP2 correlation energy.
 * `eri`: 8-fold symmetric ERI data
 * `coeff`: MO coefficient matrix (N*N, row-major)
 * `epsilon`: orbital energies (N elements)
 * `nocc`: number of occupied spatial orbitals
 * `num_basis`: number of basis functions
 * Returns MP2 correlation energy (f64 scalar)
 * @param {Float64Array} eri
 * @param {Float64Array} coeff
 * @param {Float64Array} epsilon
 * @param {number} nocc
 * @param {number} num_basis
 * @returns {number}
 */
export function compute_mp2_wasm(eri, coeff, epsilon, nocc, num_basis) {
    const ptr0 = passArrayF64ToWasm0(eri, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(coeff, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArrayF64ToWasm0(epsilon, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ret = wasm.compute_mp2_wasm(ptr0, len0, ptr1, len1, ptr2, len2, nocc, num_basis);
    return ret;
}

/**
 * Compute MP3 correlation energy (RHF spatial-orbital formulas).
 * Returns Vec<f64> with 2 elements: [mp2_corr, mp3_corr]
 * @param {Float64Array} eri
 * @param {Float64Array} coeff
 * @param {Float64Array} epsilon
 * @param {number} nocc
 * @param {number} num_basis
 * @returns {Float64Array}
 */
export function compute_mp3_wasm(eri, coeff, epsilon, nocc, num_basis) {
    const ptr0 = passArrayF64ToWasm0(eri, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(coeff, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArrayF64ToWasm0(epsilon, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ret = wasm.compute_mp3_wasm(ptr0, len0, ptr1, len1, ptr2, len2, nocc, num_basis);
    var v4 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
    return v4;
}

/**
 * Compute RI Coulomb matrix J.
 * `b_matrix`: B[P*N*N + mu*N + nu], length = naux * N * N
 * `density`: P matrix (N*N, row-major)
 * Returns J matrix (N*N, row-major)
 * @param {Float64Array} b_matrix
 * @param {Float64Array} density
 * @param {number} naux
 * @param {number} num_basis
 * @returns {Float64Array}
 */
export function compute_ri_coulomb_wasm(b_matrix, density, naux, num_basis) {
    const ptr0 = passArrayF64ToWasm0(b_matrix, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(density, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ret = wasm.compute_ri_coulomb_wasm(ptr0, len0, ptr1, len1, naux, num_basis);
    var v3 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
    return v3;
}

/**
 * Compute RI Exchange matrix K.
 * `coefficients`: MO coefficient matrix (N*N, row-major)
 * `nocc`: number of occupied orbitals
 * Returns K matrix (N*N, row-major)
 * @param {Float64Array} b_matrix
 * @param {Float64Array} coefficients
 * @param {number} naux
 * @param {number} num_basis
 * @param {number} nocc
 * @returns {Float64Array}
 */
export function compute_ri_exchange_wasm(b_matrix, coefficients, naux, num_basis, nocc) {
    const ptr0 = passArrayF64ToWasm0(b_matrix, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(coefficients, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ret = wasm.compute_ri_exchange_wasm(ptr0, len0, ptr1, len1, naux, num_basis, nocc);
    var v3 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
    return v3;
}

/**
 * Compute RI-RHF Fock matrix: F = H + J - K
 * @param {Float64Array} b_matrix
 * @param {Float64Array} density
 * @param {Float64Array} core_h
 * @param {Float64Array} coefficients
 * @param {number} naux
 * @param {number} num_basis
 * @param {number} nocc
 * @returns {Float64Array}
 */
export function compute_ri_fock_rhf_wasm(b_matrix, density, core_h, coefficients, naux, num_basis, nocc) {
    const ptr0 = passArrayF64ToWasm0(b_matrix, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(density, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArrayF64ToWasm0(core_h, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ptr3 = passArrayF64ToWasm0(coefficients, wasm.__wbindgen_malloc);
    const len3 = WASM_VECTOR_LEN;
    const ret = wasm.compute_ri_fock_rhf_wasm(ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3, naux, num_basis, nocc);
    var v5 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
    return v5;
}

/**
 * Compute RI-UHF Fock matrices: Fα = H + J(Ptot) - Kα, Fβ = H + J(Ptot) - Kβ
 * Returns concatenated [Fα..., Fβ...] (2*N*N elements)
 * @param {Float64Array} b_matrix
 * @param {Float64Array} density_total
 * @param {Float64Array} core_h
 * @param {Float64Array} coeff_alpha
 * @param {number} nocc_alpha
 * @param {Float64Array} coeff_beta
 * @param {number} nocc_beta
 * @param {number} naux
 * @param {number} num_basis
 * @returns {Float64Array}
 */
export function compute_ri_fock_uhf_wasm(b_matrix, density_total, core_h, coeff_alpha, nocc_alpha, coeff_beta, nocc_beta, naux, num_basis) {
    const ptr0 = passArrayF64ToWasm0(b_matrix, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(density_total, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArrayF64ToWasm0(core_h, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ptr3 = passArrayF64ToWasm0(coeff_alpha, wasm.__wbindgen_malloc);
    const len3 = WASM_VECTOR_LEN;
    const ptr4 = passArrayF64ToWasm0(coeff_beta, wasm.__wbindgen_malloc);
    const len4 = WASM_VECTOR_LEN;
    const ret = wasm.compute_ri_fock_uhf_wasm(ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3, nocc_alpha, ptr4, len4, nocc_beta, naux, num_basis);
    var v6 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
    return v6;
}

/**
 * Compute RI-MP2 correlation energy (RHF).
 * @param {Float64Array} b_matrix
 * @param {Float64Array} coefficients
 * @param {Float64Array} epsilon
 * @param {number} naux
 * @param {number} nocc
 * @param {number} num_basis
 * @returns {number}
 */
export function compute_ri_mp2_wasm(b_matrix, coefficients, epsilon, naux, nocc, num_basis) {
    const ptr0 = passArrayF64ToWasm0(b_matrix, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(coefficients, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArrayF64ToWasm0(epsilon, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ret = wasm.compute_ri_mp2_wasm(ptr0, len0, ptr1, len1, ptr2, len2, naux, nocc, num_basis);
    return ret;
}

/**
 * Compute RI Setup: 2c/3c integrals, Cholesky, B matrix.
 * Returns Vec<f64> with layout [naux, B_matrix...]
 * @param {Float64Array} primary_shells_flat
 * @param {Float64Array} primary_norm
 * @param {number} nbasis
 * @param {Float64Array} aux_shells_flat
 * @param {Float64Array} aux_norm
 * @param {number} naux
 * @returns {Float64Array}
 */
export function compute_ri_setup_wasm(primary_shells_flat, primary_norm, nbasis, aux_shells_flat, aux_norm, naux) {
    const ptr0 = passArrayF64ToWasm0(primary_shells_flat, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(primary_norm, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArrayF64ToWasm0(aux_shells_flat, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ptr3 = passArrayF64ToWasm0(aux_norm, wasm.__wbindgen_malloc);
    const len3 = WASM_VECTOR_LEN;
    const ret = wasm.compute_ri_setup_wasm(ptr0, len0, ptr1, len1, nbasis, ptr2, len2, ptr3, len3, naux);
    var v5 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
    return v5;
}

/**
 * Compute RI-UMP2 correlation energy.
 * @param {Float64Array} b_matrix
 * @param {Float64Array} ca
 * @param {Float64Array} cb
 * @param {Float64Array} epsilon_a
 * @param {Float64Array} epsilon_b
 * @param {number} naux
 * @param {number} nocc_a
 * @param {number} nocc_b
 * @param {number} num_basis
 * @returns {number}
 */
export function compute_ri_ump2_wasm(b_matrix, ca, cb, epsilon_a, epsilon_b, naux, nocc_a, nocc_b, num_basis) {
    const ptr0 = passArrayF64ToWasm0(b_matrix, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(ca, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArrayF64ToWasm0(cb, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ptr3 = passArrayF64ToWasm0(epsilon_a, wasm.__wbindgen_malloc);
    const len3 = WASM_VECTOR_LEN;
    const ptr4 = passArrayF64ToWasm0(epsilon_b, wasm.__wbindgen_malloc);
    const len4 = WASM_VECTOR_LEN;
    const ret = wasm.compute_ri_ump2_wasm(ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3, ptr4, len4, naux, nocc_a, nocc_b, num_basis);
    return ret;
}

/**
 * Compute UCCSD correlation energy from converged UHF results.
 * Returns f64 scalar (correlation energy).
 * @param {Float64Array} eri
 * @param {Float64Array} ca
 * @param {Float64Array} cb
 * @param {Float64Array} epsilon_a
 * @param {Float64Array} epsilon_b
 * @param {number} nocc_a
 * @param {number} nocc_b
 * @param {number} num_basis
 * @returns {number}
 */
export function compute_uccsd_wasm(eri, ca, cb, epsilon_a, epsilon_b, nocc_a, nocc_b, num_basis) {
    const ptr0 = passArrayF64ToWasm0(eri, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(ca, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArrayF64ToWasm0(cb, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ptr3 = passArrayF64ToWasm0(epsilon_a, wasm.__wbindgen_malloc);
    const len3 = WASM_VECTOR_LEN;
    const ptr4 = passArrayF64ToWasm0(epsilon_b, wasm.__wbindgen_malloc);
    const len4 = WASM_VECTOR_LEN;
    const ret = wasm.compute_uccsd_wasm(ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3, ptr4, len4, nocc_a, nocc_b, num_basis);
    return ret;
}

/**
 * Compute UMP2 correlation energy from converged UHF results.
 * @param {Float64Array} eri
 * @param {Float64Array} ca
 * @param {Float64Array} cb
 * @param {Float64Array} epsilon_a
 * @param {Float64Array} epsilon_b
 * @param {number} nocc_a
 * @param {number} nocc_b
 * @param {number} num_basis
 * @returns {number}
 */
export function compute_ump2_wasm(eri, ca, cb, epsilon_a, epsilon_b, nocc_a, nocc_b, num_basis) {
    const ptr0 = passArrayF64ToWasm0(eri, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(ca, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArrayF64ToWasm0(cb, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ptr3 = passArrayF64ToWasm0(epsilon_a, wasm.__wbindgen_malloc);
    const len3 = WASM_VECTOR_LEN;
    const ptr4 = passArrayF64ToWasm0(epsilon_b, wasm.__wbindgen_malloc);
    const len4 = WASM_VECTOR_LEN;
    const ret = wasm.compute_ump2_wasm(ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3, ptr4, len4, nocc_a, nocc_b, num_basis);
    return ret;
}

/**
 * Compute UMP3 correlation energy from converged UHF results.
 * Returns Vec<f64> with 2 elements: [mp2_corr, mp3_corr]
 * @param {Float64Array} eri
 * @param {Float64Array} ca
 * @param {Float64Array} cb
 * @param {Float64Array} epsilon_a
 * @param {Float64Array} epsilon_b
 * @param {number} nocc_a
 * @param {number} nocc_b
 * @param {number} num_basis
 * @returns {Float64Array}
 */
export function compute_ump3_wasm(eri, ca, cb, epsilon_a, epsilon_b, nocc_a, nocc_b, num_basis) {
    const ptr0 = passArrayF64ToWasm0(eri, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(ca, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArrayF64ToWasm0(cb, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ptr3 = passArrayF64ToWasm0(epsilon_a, wasm.__wbindgen_malloc);
    const len3 = WASM_VECTOR_LEN;
    const ptr4 = passArrayF64ToWasm0(epsilon_b, wasm.__wbindgen_malloc);
    const len4 = WASM_VECTOR_LEN;
    const ret = wasm.compute_ump3_wasm(ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3, ptr4, len4, nocc_a, nocc_b, num_basis);
    var v6 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
    return v6;
}

/**
 * Compute XC energy and V_xc matrix on a numerical grid.
 *
 * Returns Vec<f64> with layout: [exc, num_electrons, vxc_a..., vxc_b...]
 * `grid_flat`: [x0, y0, z0, w0, x1, y1, z1, w1, ...] (4 values per point)
 * `density_b`: empty for RKS (restricted); N*N for UKS (unrestricted)
 * `func_id`: 0=SVWN, 1=BLYP, 2=PBE, 3=B3LYP
 * `need_grad`: true for GGA/hybrid functionals
 * @param {Float64Array} shells_flat
 * @param {Float64Array} norm_factors
 * @param {Float64Array} density_a
 * @param {Float64Array} density_b
 * @param {Float64Array} grid_flat
 * @param {number} num_basis
 * @param {number} func_id
 * @param {boolean} need_grad
 * @returns {Float64Array}
 */
export function compute_xc_wasm(shells_flat, norm_factors, density_a, density_b, grid_flat, num_basis, func_id, need_grad) {
    const ptr0 = passArrayF64ToWasm0(shells_flat, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(norm_factors, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArrayF64ToWasm0(density_a, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ptr3 = passArrayF64ToWasm0(density_b, wasm.__wbindgen_malloc);
    const len3 = WASM_VECTOR_LEN;
    const ptr4 = passArrayF64ToWasm0(grid_flat, wasm.__wbindgen_malloc);
    const len4 = WASM_VECTOR_LEN;
    const ret = wasm.compute_xc_wasm(ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3, ptr4, len4, num_basis, func_id, need_grad);
    var v6 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
    return v6;
}

/**
 * Solve CPHF equations via preconditioned CG.
 * mo_eri: MO-basis ERIs (n⁴)
 * eps: orbital energies (nmo)
 * rhs: right-hand side [n_pert × dim], dim = nocc*nvir
 * Returns: solution U^x [n_pert × dim]
 * @param {Float64Array} mo_eri
 * @param {Float64Array} eps
 * @param {Float64Array} rhs
 * @param {number} nocc
 * @param {number} nvir
 * @param {number} nmo
 * @param {number} n_pert
 * @param {number} tol
 * @param {number} max_iter
 * @returns {Float64Array}
 */
export function solve_cphf_wasm(mo_eri, eps, rhs, nocc, nvir, nmo, n_pert, tol, max_iter) {
    const ptr0 = passArrayF64ToWasm0(mo_eri, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(eps, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArrayF64ToWasm0(rhs, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ret = wasm.solve_cphf_wasm(ptr0, len0, ptr1, len1, ptr2, len2, nocc, nvir, nmo, n_pert, tol, max_iter);
    var v4 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
    return v4;
}

/**
 * AO → MO integral transform (4-quarter, O(N⁵)).
 * c_flat: MO coefficient matrix (n×n, row-major)
 * eri_sym: AO ERIs in 8-fold symmetric packed storage
 * num_basis: basis set size
 * Returns: full (pq|rs) as n⁴ flat array
 * @param {Float64Array} c_flat
 * @param {Float64Array} eri_sym
 * @param {number} num_basis
 * @returns {Float64Array}
 */
export function transform_eri_to_mo_wasm(c_flat, eri_sym, num_basis) {
    const ptr0 = passArrayF64ToWasm0(c_flat, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayF64ToWasm0(eri_sym, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ret = wasm.transform_eri_to_mo_wasm(ptr0, len0, ptr1, len1, num_basis);
    var v3 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
    return v3;
}

function __wbg_get_imports() {
    const import0 = {
        __proto__: null,
        __wbindgen_init_externref_table: function() {
            const table = wasm.__wbindgen_externrefs;
            const offset = table.grow(4);
            table.set(0, undefined);
            table.set(offset + 0, undefined);
            table.set(offset + 1, null);
            table.set(offset + 2, true);
            table.set(offset + 3, false);
        },
    };
    return {
        __proto__: null,
        "./wasm_eri_bg.js": import0,
    };
}

function getArrayF64FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getFloat64ArrayMemory0().subarray(ptr / 8, ptr / 8 + len);
}

let cachedFloat64ArrayMemory0 = null;
function getFloat64ArrayMemory0() {
    if (cachedFloat64ArrayMemory0 === null || cachedFloat64ArrayMemory0.byteLength === 0) {
        cachedFloat64ArrayMemory0 = new Float64Array(wasm.memory.buffer);
    }
    return cachedFloat64ArrayMemory0;
}

function passArrayF64ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 8, 8) >>> 0;
    getFloat64ArrayMemory0().set(arg, ptr / 8);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
}

let WASM_VECTOR_LEN = 0;

let wasmModule, wasm;
function __wbg_finalize_init(instance, module) {
    wasm = instance.exports;
    wasmModule = module;
    cachedFloat64ArrayMemory0 = null;
    wasm.__wbindgen_start();
    return wasm;
}

async function __wbg_load(module, imports) {
    if (typeof Response === 'function' && module instanceof Response) {
        if (typeof WebAssembly.instantiateStreaming === 'function') {
            try {
                return await WebAssembly.instantiateStreaming(module, imports);
            } catch (e) {
                const validResponse = module.ok && expectedResponseType(module.type);

                if (validResponse && module.headers.get('Content-Type') !== 'application/wasm') {
                    console.warn("`WebAssembly.instantiateStreaming` failed because your server does not serve Wasm with `application/wasm` MIME type. Falling back to `WebAssembly.instantiate` which is slower. Original error:\n", e);

                } else { throw e; }
            }
        }

        const bytes = await module.arrayBuffer();
        return await WebAssembly.instantiate(bytes, imports);
    } else {
        const instance = await WebAssembly.instantiate(module, imports);

        if (instance instanceof WebAssembly.Instance) {
            return { instance, module };
        } else {
            return instance;
        }
    }

    function expectedResponseType(type) {
        switch (type) {
            case 'basic': case 'cors': case 'default': return true;
        }
        return false;
    }
}

function initSync(module) {
    if (wasm !== undefined) return wasm;


    if (module !== undefined) {
        if (Object.getPrototypeOf(module) === Object.prototype) {
            ({module} = module)
        } else {
            console.warn('using deprecated parameters for `initSync()`; pass a single object instead')
        }
    }

    const imports = __wbg_get_imports();
    if (!(module instanceof WebAssembly.Module)) {
        module = new WebAssembly.Module(module);
    }
    const instance = new WebAssembly.Instance(module, imports);
    return __wbg_finalize_init(instance, module);
}

async function __wbg_init(module_or_path) {
    if (wasm !== undefined) return wasm;


    if (module_or_path !== undefined) {
        if (Object.getPrototypeOf(module_or_path) === Object.prototype) {
            ({module_or_path} = module_or_path)
        } else {
            console.warn('using deprecated parameters for the initialization function; pass a single object instead')
        }
    }

    if (module_or_path === undefined) {
        module_or_path = new URL('wasm_eri_bg.wasm', import.meta.url);
    }
    const imports = __wbg_get_imports();

    if (typeof module_or_path === 'string' || (typeof Request === 'function' && module_or_path instanceof Request) || (typeof URL === 'function' && module_or_path instanceof URL)) {
        module_or_path = fetch(module_or_path);
    }

    const { instance, module } = await __wbg_load(await module_or_path, imports);

    return __wbg_finalize_init(instance, module);
}

export { initSync, __wbg_init as default };
