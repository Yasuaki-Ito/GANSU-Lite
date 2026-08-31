/** Small linear system solver (Gauss elimination with partial pivoting) for DIIS */

/** Solve Ax = b in-place. A is n x n, b is length n.
 *  Returns x as a new Float64Array. */
export function solveLinearSystem(A: Float64Array, b: Float64Array, n: number): Float64Array {
  // Make copies to avoid modifying inputs
  const a = new Float64Array(A);
  const x = new Float64Array(b);

  // Forward elimination with partial pivoting
  for (let col = 0; col < n; col++) {
    // Find pivot
    let maxVal = Math.abs(a[col * n + col]);
    let maxRow = col;
    for (let row = col + 1; row < n; row++) {
      const v = Math.abs(a[row * n + col]);
      if (v > maxVal) {
        maxVal = v;
        maxRow = row;
      }
    }

    // Swap rows
    if (maxRow !== col) {
      for (let j = 0; j < n; j++) {
        const tmp = a[col * n + j];
        a[col * n + j] = a[maxRow * n + j];
        a[maxRow * n + j] = tmp;
      }
      const tmp = x[col];
      x[col] = x[maxRow];
      x[maxRow] = tmp;
    }

    const pivot = a[col * n + col];
    if (Math.abs(pivot) < 1e-30) continue;

    // Eliminate below
    for (let row = col + 1; row < n; row++) {
      const factor = a[row * n + col] / pivot;
      for (let j = col; j < n; j++) {
        a[row * n + j] -= factor * a[col * n + j];
      }
      x[row] -= factor * x[col];
    }
  }

  // Back substitution
  for (let row = n - 1; row >= 0; row--) {
    for (let j = row + 1; j < n; j++) {
      x[row] -= a[row * n + j] * x[j];
    }
    const diag = a[row * n + row];
    if (Math.abs(diag) > 1e-30) {
      x[row] /= diag;
    } else {
      x[row] = 0;
    }
  }

  return x;
}
