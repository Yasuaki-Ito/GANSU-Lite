/** Matrix multiplication — naive O(n^3) for small matrices */

import { Matrix } from './matrix';

/** C = A * B */
export function matmul(A: Matrix, B: Matrix): Matrix {
  const m = A.rows, n = B.cols, k = A.cols;
  const C = new Matrix(m, n);
  for (let i = 0; i < m; i++) {
    for (let p = 0; p < k; p++) {
      const aip = A.data[i * k + p];
      if (aip === 0) continue;
      for (let j = 0; j < n; j++) {
        C.data[i * n + j] += aip * B.data[p * n + j];
      }
    }
  }
  return C;
}

/** C = A^T * B */
export function matmulAtB(A: Matrix, B: Matrix): Matrix {
  const m = A.cols, n = B.cols, k = A.rows;
  const C = new Matrix(m, n);
  for (let p = 0; p < k; p++) {
    for (let i = 0; i < m; i++) {
      const api = A.data[p * A.cols + i];
      if (api === 0) continue;
      for (let j = 0; j < n; j++) {
        C.data[i * n + j] += api * B.data[p * n + j];
      }
    }
  }
  return C;
}

/** C = A * B^T */
export function matmulABt(A: Matrix, B: Matrix): Matrix {
  const m = A.rows, n = B.rows, k = A.cols;
  const C = new Matrix(m, n);
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      let sum = 0;
      for (let p = 0; p < k; p++) {
        sum += A.data[i * k + p] * B.data[j * k + p];
      }
      C.data[i * n + j] = sum;
    }
  }
  return C;
}

/** Tr(A * B) without allocating the full product */
export function traceProduct(A: Matrix, B: Matrix): number {
  const n = A.rows;
  let tr = 0;
  for (let i = 0; i < n; i++) {
    for (let k = 0; k < n; k++) {
      tr += A.data[i * n + k] * B.data[k * n + i];
    }
  }
  return tr;
}
