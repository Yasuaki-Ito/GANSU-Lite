/// WASM-accelerated Fock matrix construction

use crate::integrals::eri_1d_index;
use crate::simd_utils;

/// RHF Fock matrix: F[i][j] = H[i][j] + Σ_{kl} P[k][l] * ((ij|kl) - 0.5*(ik|jl))
///
/// `eri`: 8-fold symmetric ERI data (flat f64)
/// `density`: P matrix (N*N, row-major)
/// `core_h`: H matrix (N*N, row-major)
/// `n`: number of basis functions
///
/// Returns F matrix (N*N, row-major)
pub fn compute_fock_rhf(eri: &[f64], density: &[f64], core_h: &[f64], n: usize) -> Vec<f64> {
    let mut f = vec![0.0_f64; n * n];
    let nn = n * n;
    let mut g_buf = vec![0.0_f64; nn];

    for i in 0..n {
        for j in i..n {
            // Pre-build combined integral array: G[k*n+l] = (ij|kl) - 0.5*(ik|jl)
            for k in 0..n {
                let k_off = k * n;
                for l in 0..n {
                    g_buf[k_off + l] = eri[eri_1d_index(i, j, k, l)]
                        - 0.5 * eri[eri_1d_index(i, k, j, l)];
                }
            }
            // SIMD dot product: g = Σ_{kl} P[k,l] * G[k,l]
            let g = simd_utils::dot(&density[..nn], &g_buf[..nn]);
            let val = core_h[i * n + j] + g;
            f[i * n + j] = val;
            f[j * n + i] = val;
        }
    }

    f
}

/// UHF Fock matrices:
///   Fa[i,j] = H[i,j] + Σ_{kl} Ptot[k,l]*(ij|kl) - Pa[k,l]*(ik|jl)
///   Fb[i,j] = H[i,j] + Σ_{kl} Ptot[k,l]*(ij|kl) - Pb[k,l]*(ik|jl)
///
/// Returns concatenated [Fa..., Fb...] (2*N*N elements)
pub fn compute_fock_uhf(
    eri: &[f64],
    density_alpha: &[f64],
    density_beta: &[f64],
    density_total: &[f64],
    core_h: &[f64],
    n: usize,
) -> Vec<f64> {
    let nn = n * n;
    let mut result = vec![0.0_f64; 2 * nn];
    let (fa, fb) = result.split_at_mut(nn);

    let mut j_buf = vec![0.0_f64; nn]; // Coulomb integrals
    let mut k_buf = vec![0.0_f64; nn]; // Exchange integrals

    for i in 0..n {
        for j in i..n {
            // Pre-build J and K integral arrays
            for k in 0..n {
                let k_off = k * n;
                for l in 0..n {
                    j_buf[k_off + l] = eri[eri_1d_index(i, j, k, l)];
                    k_buf[k_off + l] = eri[eri_1d_index(i, k, j, l)];
                }
            }
            // Fa: Σ Ptot*(ij|kl) - Pa*(ik|jl)
            let ga = simd_utils::dot(&density_total[..nn], &j_buf[..nn])
                   - simd_utils::dot(&density_alpha[..nn], &k_buf[..nn]);
            // Fb: Σ Ptot*(ij|kl) - Pb*(ik|jl)
            let gb = simd_utils::dot(&density_total[..nn], &j_buf[..nn])
                   - simd_utils::dot(&density_beta[..nn], &k_buf[..nn]);

            let ha = core_h[i * n + j] + ga;
            let hb = core_h[i * n + j] + gb;
            fa[i * n + j] = ha;
            fa[j * n + i] = ha;
            fb[i * n + j] = hb;
            fb[j * n + i] = hb;
        }
    }

    result
}
