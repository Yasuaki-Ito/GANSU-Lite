/** Jacobi eigendecomposition for real symmetric matrices */

import { Matrix, createFloatArray, type FloatArray } from './matrix';

/** Returns { eigenvalues, eigenvectors } where eigenvectors columns are the eigenvectors.
 *  eigenvalues are sorted ascending. */
export function jacobiEigen(A: Matrix): { eigenvalues: FloatArray; eigenvectors: Matrix } {
  const n = A.rows;
  const S = A.clone(); // working copy
  const V = Matrix.identity(n);
  const maxIter = 100 * n * n;
  const tol = 1e-12;

  for (let iter = 0; iter < maxIter; iter++) {
    // Find the largest off-diagonal element
    let maxVal = 0;
    let p = 0, q = 1;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const v = Math.abs(S.get(i, j));
        if (v > maxVal) {
          maxVal = v;
          p = i;
          q = j;
        }
      }
    }

    if (maxVal < tol) break;

    // Compute rotation
    const app = S.get(p, p);
    const aqq = S.get(q, q);
    const apq = S.get(p, q);

    let theta: number;
    if (Math.abs(app - aqq) < 1e-30) {
      theta = Math.PI / 4;
    } else {
      theta = 0.5 * Math.atan2(2 * apq, app - aqq);
    }

    const c = Math.cos(theta);
    const s = Math.sin(theta);

    // Apply rotation to S: S' = G^T * S * G
    // Update rows/cols p and q
    for (let i = 0; i < n; i++) {
      if (i === p || i === q) continue;
      const sip = S.get(i, p);
      const siq = S.get(i, q);
      S.set(i, p, c * sip + s * siq);
      S.set(p, i, c * sip + s * siq);
      S.set(i, q, -s * sip + c * siq);
      S.set(q, i, -s * sip + c * siq);
    }

    const newApp = c * c * app + 2 * s * c * apq + s * s * aqq;
    const newAqq = s * s * app - 2 * s * c * apq + c * c * aqq;
    S.set(p, p, newApp);
    S.set(q, q, newAqq);
    S.set(p, q, 0);
    S.set(q, p, 0);

    // Update eigenvector matrix V
    for (let i = 0; i < n; i++) {
      const vip = V.get(i, p);
      const viq = V.get(i, q);
      V.set(i, p, c * vip + s * viq);
      V.set(i, q, -s * vip + c * viq);
    }
  }

  // Extract eigenvalues
  const eigenvalues = createFloatArray(n);
  for (let i = 0; i < n; i++) eigenvalues[i] = S.get(i, i);

  // Sort by eigenvalue (ascending) and reorder eigenvectors
  const indices = Array.from({ length: n }, (_, i) => i);
  indices.sort((a, b) => eigenvalues[a] - eigenvalues[b]);

  const sortedEvals = createFloatArray(n);
  const sortedEvecs = new Matrix(n, n);
  for (let j = 0; j < n; j++) {
    sortedEvals[j] = eigenvalues[indices[j]];
    for (let i = 0; i < n; i++) {
      sortedEvecs.set(i, j, V.get(i, indices[j]));
    }
  }

  return { eigenvalues: sortedEvals, eigenvectors: sortedEvecs };
}

/** Async variant of jacobiEigen that yields to the event loop periodically.
 *  Use for large matrices on the main thread to avoid browser timeout. */
export async function jacobiEigenAsync(
  A: Matrix,
  yieldInterval = 50,
): Promise<{ eigenvalues: FloatArray; eigenvectors: Matrix }> {
  const n = A.rows;
  const S = A.clone();
  const V = Matrix.identity(n);
  const maxIter = 100 * n * n;
  const tol = 1e-12;

  for (let iter = 0; iter < maxIter; iter++) {
    let maxVal = 0;
    let p = 0, q = 1;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const v = Math.abs(S.get(i, j));
        if (v > maxVal) {
          maxVal = v;
          p = i;
          q = j;
        }
      }
    }

    if (maxVal < tol) break;

    const app = S.get(p, p);
    const aqq = S.get(q, q);
    const apq = S.get(p, q);

    let theta: number;
    if (Math.abs(app - aqq) < 1e-30) {
      theta = Math.PI / 4;
    } else {
      theta = 0.5 * Math.atan2(2 * apq, app - aqq);
    }

    const c = Math.cos(theta);
    const s = Math.sin(theta);

    for (let i = 0; i < n; i++) {
      if (i === p || i === q) continue;
      const sip = S.get(i, p);
      const siq = S.get(i, q);
      S.set(i, p, c * sip + s * siq);
      S.set(p, i, c * sip + s * siq);
      S.set(i, q, -s * sip + c * siq);
      S.set(q, i, -s * sip + c * siq);
    }

    const newApp = c * c * app + 2 * s * c * apq + s * s * aqq;
    const newAqq = s * s * app - 2 * s * c * apq + c * c * aqq;
    S.set(p, p, newApp);
    S.set(q, q, newAqq);
    S.set(p, q, 0);
    S.set(q, p, 0);

    for (let i = 0; i < n; i++) {
      const vip = V.get(i, p);
      const viq = V.get(i, q);
      V.set(i, p, c * vip + s * viq);
      V.set(i, q, -s * vip + c * viq);
    }

    if (iter % yieldInterval === yieldInterval - 1) {
      await new Promise<void>(r => setTimeout(r, 0));
    }
  }

  const eigenvalues = createFloatArray(n);
  for (let i = 0; i < n; i++) eigenvalues[i] = S.get(i, i);

  const indices = Array.from({ length: n }, (_, i) => i);
  indices.sort((a, b) => eigenvalues[a] - eigenvalues[b]);

  const sortedEvals = createFloatArray(n);
  const sortedEvecs = new Matrix(n, n);
  for (let j = 0; j < n; j++) {
    sortedEvals[j] = eigenvalues[indices[j]];
    for (let i = 0; i < n; i++) {
      sortedEvecs.set(i, j, V.get(i, indices[j]));
    }
  }

  return { eigenvalues: sortedEvals, eigenvectors: sortedEvecs };
}
