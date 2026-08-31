/** Boys function — port of GANSU boys.hpp */

import { BOYS_GRID } from '../data/boysGrid';
import { A_TR, B_TR, LUT_XI_INTERVAL, LUT_NUM_XI, LUT_K_MAX } from './constants';

const A_RS = 0.064048916778075;
const B_RS = 28.487431543672;
// Maclaurin expansion parameters (reserved for future use)
// const AA = 0.03768724;
// const BB = 0.60549623;
// const CC = 6.32743473;
// const DD = 10.350421;

/** Approximate erf(x) via Abramowitz-Stegun 7.1.26 (max error ~1.5e-7) */
function erf(x: number): number {
  const sign = x < 0 ? -1 : 1;
  const t = 1.0 / (1.0 + 0.3275911 * Math.abs(x));
  const y = 1.0 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
  return sign * y;
}

/** Single Boys function F_n(x) */
export function boysFunction(n: number, x: number): number {
  if (x === 0.0) {
    return 1.0 / (2 * n + 1);
  }
  if (x < A_TR * n + B_TR) {
    return taylorInterpolation(n, x);
  }
  return recurrenceSemiInfinite(n, x);
}

/** Range-separated Boys function for erf(ωr)/r (long-range) Coulomb operator.
 *  F_n^{erf}(x, η, ω) = ξ^(n+1/2) · F_n(ξ x),  where ξ = ω²/(η + ω²).
 *  Reduces to standard F_n when ω → ∞, and to 0 when ω → 0.
 *  @param n  highest order requested (returns F_0 ... F_n)
 *  @param x  standard Boys argument (η · |R_PQ|²)
 *  @param eta η = pq/(p+q) for the shell quartet
 *  @param omega range-separation parameter (Bohr⁻¹) */
export function boysAllOmega(n: number, x: number, eta: number, omega: number): Float64Array {
  if (omega <= 0) return new Float64Array(n + 1); // erf with ω=0 is zero
  const omega2 = omega * omega;
  const xi = omega2 / (eta + omega2);
  const xPrime = xi * x;
  const fStandard = boysAll(n, xPrime);
  const result = new Float64Array(n + 1);
  let factor = Math.sqrt(xi); // ξ^(0+1/2)
  for (let j = 0; j <= n; j++) {
    result[j] = factor * fStandard[j];
    factor *= xi; // ξ^(j+1+1/2) for next iteration
  }
  return result;
}

/** Compute Boys functions F_0(x), F_1(x), ..., F_n(x) */
export function boysAll(n: number, x: number): Float64Array {
  const boys = new Float64Array(n + 1);

  if (x === 0.0) {
    for (let j = 0; j <= n; j++) boys[j] = 1.0 / (2 * j + 1);
    return boys;
  }

  if (x < A_TR * n + B_TR) {
    // Taylor interpolation for all orders
    taylorInterpolationAll(n, x, boys);
  } else {
    // Recurrence from semi-infinite
    recurrenceSemiInfiniteAll(n, x, boys);
  }

  return boys;
}

function taylorInterpolation(n: number, x: number): number {
  const xIdx = Math.min(Math.round(x / LUT_XI_INTERVAL), LUT_NUM_XI - 1);
  const deltaX = x - LUT_XI_INTERVAL * xIdx;
  let numerator = 1.0;
  let factorial = 1;
  let Fx = BOYS_GRID[LUT_NUM_XI * n + xIdx];

  for (let k = 1; k <= LUT_K_MAX; k++) {
    numerator *= -deltaX;
    factorial *= k;
    Fx += (BOYS_GRID[LUT_NUM_XI * (n + k) + xIdx] * numerator) / factorial;
  }
  return Fx;
}

function taylorInterpolationAll(n: number, x: number, boys: Float64Array) {
  const xIdx = Math.min(Math.round(x / LUT_XI_INTERVAL), LUT_NUM_XI - 1);
  const deltaX = x - LUT_XI_INTERVAL * xIdx;

  for (let j = 0; j <= n; j++) {
    let numerator = 1.0;
    let factorial = 1;
    let Fx = BOYS_GRID[LUT_NUM_XI * j + xIdx];

    for (let k = 1; k <= LUT_K_MAX; k++) {
      numerator *= -deltaX;
      factorial *= k;
      Fx += (BOYS_GRID[LUT_NUM_XI * (j + k) + xIdx] * numerator) / factorial;
    }
    boys[j] = Fx;
  }
}

function recurrenceSemiInfinite(n: number, x: number): number {
  let expNegX = 0.0;
  const recip2x = 1.0 / (2.0 * x);
  let Fx = 0.5 * Math.sqrt(Math.PI / x);

  if (x < A_RS * n + B_RS) {
    expNegX = Math.exp(-x);
    Fx *= erf(Math.sqrt(x));
  }

  for (let j = 1; j <= n; j++) {
    Fx = ((2 * j - 1) * Fx - expNegX) * recip2x;
  }
  return Fx;
}

function recurrenceSemiInfiniteAll(n: number, x: number, boys: Float64Array) {
  let expNegX = 0.0;
  const recip2x = 1.0 / (2.0 * x);
  let Fx = 0.5 * Math.sqrt(Math.PI / x);

  if (x < A_RS * n + B_RS) {
    expNegX = Math.exp(-x);
    Fx *= erf(Math.sqrt(x));
  }

  boys[0] = Fx;
  for (let j = 1; j <= n; j++) {
    Fx = ((2 * j - 1) * Fx - expNegX) * recip2x;
    boys[j] = Fx;
  }
}
