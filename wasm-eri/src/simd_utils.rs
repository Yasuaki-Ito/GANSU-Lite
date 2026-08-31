/// SIMD-accelerated dot-product and DAXPY utilities.
///
/// When compiled with `-C target-feature=+simd128`, uses wasm f64x2 intrinsics.
/// Otherwise falls back to scalar loops (auto-vectorization may still apply).

// ---------------------------------------------------------------------------
// SIMD path (wasm32 + simd128)
// ---------------------------------------------------------------------------

#[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
pub mod inner {
    use std::arch::wasm32::*;

    /// Dot product: Σ a[i] * b[i]
    #[inline]
    pub fn dot(a: &[f64], b: &[f64]) -> f64 {
        let n = a.len().min(b.len());
        let chunks = n / 4; // process 4 elements per iteration (2 x f64x2)
        let mut acc0 = f64x2_splat(0.0);
        let mut acc1 = f64x2_splat(0.0);

        let a_ptr = a.as_ptr();
        let b_ptr = b.as_ptr();

        for i in 0..chunks {
            let offset = i * 4;
            unsafe {
                let va0 = v128_load(a_ptr.add(offset) as *const v128);
                let vb0 = v128_load(b_ptr.add(offset) as *const v128);
                let va1 = v128_load(a_ptr.add(offset + 2) as *const v128);
                let vb1 = v128_load(b_ptr.add(offset + 2) as *const v128);
                acc0 = f64x2_add(acc0, f64x2_mul(va0, vb0));
                acc1 = f64x2_add(acc1, f64x2_mul(va1, vb1));
            }
        }

        let combined = f64x2_add(acc0, acc1);
        let mut sum = f64x2_extract_lane::<0>(combined) + f64x2_extract_lane::<1>(combined);

        // Remainder
        for i in (chunks * 4)..n {
            sum += a[i] * b[i];
        }
        sum
    }

    /// DAXPY: y[i] += alpha * x[i]
    #[inline]
    pub fn daxpy(alpha: f64, x: &[f64], y: &mut [f64]) {
        let n = x.len().min(y.len());
        let chunks = n / 2;
        let va = f64x2_splat(alpha);

        let x_ptr = x.as_ptr();
        let y_ptr = y.as_mut_ptr();

        for i in 0..chunks {
            let offset = i * 2;
            unsafe {
                let vx = v128_load(x_ptr.add(offset) as *const v128);
                let vy = v128_load(y_ptr.add(offset) as *const v128);
                let result = f64x2_add(vy, f64x2_mul(va, vx));
                v128_store(y_ptr.add(offset) as *mut v128, result);
            }
        }

        // Remainder
        for i in (chunks * 2)..n {
            y[i] += alpha * x[i];
        }
    }

}

// ---------------------------------------------------------------------------
// Scalar fallback
// ---------------------------------------------------------------------------

#[cfg(not(all(target_arch = "wasm32", target_feature = "simd128")))]
pub mod inner {
    /// Dot product: Σ a[i] * b[i]
    #[inline]
    pub fn dot(a: &[f64], b: &[f64]) -> f64 {
        let n = a.len().min(b.len());
        let mut sum = 0.0;
        for i in 0..n {
            sum += a[i] * b[i];
        }
        sum
    }

    /// DAXPY: y[i] += alpha * x[i]
    #[inline]
    pub fn daxpy(alpha: f64, x: &[f64], y: &mut [f64]) {
        let n = x.len().min(y.len());
        for i in 0..n {
            y[i] += alpha * x[i];
        }
    }
}

// Re-export for convenience
pub use inner::{daxpy, dot};
