/** Dense matrix class — Float64 only (row-major) */

export type FloatArray = Float64Array;

export function createFloatArray(n: number): Float64Array {
  return new Float64Array(n);
}

export function copyToFloatArray(src: ArrayLike<number>): Float64Array {
  return new Float64Array(src);
}

export class Matrix {
  readonly rows: number;
  readonly cols: number;
  readonly data: Float64Array;

  constructor(rows: number, cols: number, data?: Float64Array) {
    this.rows = rows;
    this.cols = cols;
    this.data = data ?? new Float64Array(rows * cols);
  }

  static zeros(n: number): Matrix {
    return new Matrix(n, n);
  }

  static identity(n: number): Matrix {
    const m = new Matrix(n, n);
    for (let i = 0; i < n; i++) m.data[i * n + i] = 1.0;
    return m;
  }

  get(i: number, j: number): number {
    return this.data[i * this.cols + j];
  }

  set(i: number, j: number, v: number) {
    this.data[i * this.cols + j] = v;
  }

  add(i: number, j: number, v: number) {
    this.data[i * this.cols + j] += v;
  }

  clone(): Matrix {
    return new Matrix(this.rows, this.cols, new Float64Array(this.data));
  }

  trace(): number {
    const n = Math.min(this.rows, this.cols);
    let s = 0;
    for (let i = 0; i < n; i++) s += this.data[i * this.cols + i];
    return s;
  }

  /** Transpose (returns new matrix) */
  transpose(): Matrix {
    const r = new Matrix(this.cols, this.rows);
    for (let i = 0; i < this.rows; i++) {
      for (let j = 0; j < this.cols; j++) {
        r.data[j * this.rows + i] = this.data[i * this.cols + j];
      }
    }
    return r;
  }

  /** Scale all elements */
  scale(s: number): Matrix {
    const r = new Matrix(this.rows, this.cols);
    for (let k = 0; k < this.data.length; k++) r.data[k] = this.data[k] * s;
    return r;
  }

  /** Element-wise addition */
  addMatrix(other: Matrix): Matrix {
    const r = new Matrix(this.rows, this.cols);
    for (let k = 0; k < this.data.length; k++) r.data[k] = this.data[k] + other.data[k];
    return r;
  }

  /** Element-wise subtraction */
  subMatrix(other: Matrix): Matrix {
    const r = new Matrix(this.rows, this.cols);
    for (let k = 0; k < this.data.length; k++) r.data[k] = this.data[k] - other.data[k];
    return r;
  }
}
