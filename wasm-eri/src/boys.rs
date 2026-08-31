/// Boys function — port of GANSU boys.hpp / boys.ts

use crate::boys_grid::BOYS_GRID;

const A_TR: f64 = 0.352905920120321;
const B_TR: f64 = 0.015532762923351;
const A_RS: f64 = 0.064048916778075;
const B_RS: f64 = 28.487431543672;
const LUT_XI_INTERVAL: f64 = 0.03125;
const LUT_NUM_XI: usize = 1024;
const LUT_K_MAX: usize = 5;
const LUT_N_RANGE: usize = 24;

/// Approximate erf(x) via Abramowitz-Stegun 7.1.26
fn erf_approx(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let t = 1.0 / (1.0 + 0.3275911 * x.abs());
    let y = 1.0
        - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t
            + 0.254829592)
            * t
            * (-x * x).exp();
    sign * y
}

/// Compute Boys functions F_0(x) through F_n(x) into `boys` slice.
pub fn boys_all(n: usize, x: f64, boys: &mut [f64]) {
    if x == 0.0 {
        for j in 0..=n {
            boys[j] = 1.0 / (2 * j + 1) as f64;
        }
        return;
    }

    if x < A_TR * n as f64 + B_TR {
        taylor_interpolation_all(n, x, boys);
    } else {
        recurrence_semi_infinite_all(n, x, boys);
    }
}

fn taylor_interpolation_all(n: usize, x: f64, boys: &mut [f64]) {
    let x_idx = ((x / LUT_XI_INTERVAL).round() as usize).min(LUT_NUM_XI - 1);
    let delta_x = x - LUT_XI_INTERVAL * x_idx as f64;

    for j in 0..=n {
        let mut numerator = 1.0_f64;
        let mut factorial = 1_u64;
        let mut fx = BOYS_GRID[LUT_NUM_XI * j + x_idx];

        let k_max = LUT_K_MAX.min(LUT_N_RANGE - j);
        for k in 1..=k_max {
            numerator *= -delta_x;
            factorial *= k as u64;
            fx += BOYS_GRID[LUT_NUM_XI * (j + k) + x_idx] * numerator / factorial as f64;
        }
        boys[j] = fx;
    }
}

fn recurrence_semi_infinite_all(n: usize, x: f64, boys: &mut [f64]) {
    let mut exp_neg_x = 0.0_f64;
    let recip2x = 1.0 / (2.0 * x);
    let mut fx = 0.5 * (std::f64::consts::PI / x).sqrt();

    if x < A_RS * n as f64 + B_RS {
        exp_neg_x = (-x).exp();
        fx *= erf_approx(x.sqrt());
    }

    boys[0] = fx;
    for j in 1..=n {
        fx = ((2 * j - 1) as f64 * fx - exp_neg_x) * recip2x;
        boys[j] = fx;
    }
}
