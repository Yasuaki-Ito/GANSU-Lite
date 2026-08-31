/** 2-Electron repulsion integrals (ERI) — multi-algorithm with automatic dispatch.
 *  Algorithms: McMurchie-Davidson (MD), Obara-Saika (OS), Rys quadrature, Head-Gordon-Pople (HGP).
 *  Selection is per shell-quartet based on total angular momentum and contraction degree. */

import type { PrimitiveShell } from './types';
import { ANGULAR_MOMENTUMS, shellTypeToNumBasis } from './constants';
import { eCoefficients, primitiveNorm } from './integrals1e';
import { boysAll, boysAllOmega } from './boys';
import { createFloatArray, type FloatArray } from '../linalg/matrix';
import { computeERIsWasm, isWasmAvailable } from './eriWasm';
import type { EriBackend, EriAlgorithm } from './hf';

/** Function signature for computing a shell quartet of ERIs */
export type ShellQuartetFn = (
  grpA: ShellGroup, grpB: ShellGroup,
  grpC: ShellGroup, grpD: ShellGroup,
  normFactors: number[], numBasis: number,
  eri: FloatArray,
  sameAB: boolean, sameCD: boolean, sameBraKet: boolean,
) => void;

/** Select the best algorithm for a given shell quartet */
function selectShellQuartetFn(
  grpA: ShellGroup, grpB: ShellGroup,
  grpC: ShellGroup, grpD: ShellGroup,
  algorithm: EriAlgorithm,
): ShellQuartetFn {
  if (algorithm !== 'auto') {
    switch (algorithm) {
      case 'md':  return computeShellQuartetMD;
      case 'os':  return computeShellQuartetOS;
      case 'rys': return computeShellQuartetRys;
      case 'hgp': return computeShellQuartetHGP;
    }
  }
  const L = grpA.shellType + grpB.shellType + grpC.shellType + grpD.shellType;
  if (L <= 1) return computeShellQuartetMD;
  if (L <= 6) {
    // HGP benefits from high contraction degree (many primitives per shell).
    // Use HGP when total primitive count is large, otherwise OS is faster.
    const nPrim = grpA.primitives.length * grpB.primitives.length
      * grpC.primitives.length * grpD.primitives.length;
    if (nPrim > 100 && L >= 2) return computeShellQuartetHGP;
    return computeShellQuartetOS;
  }
  return computeShellQuartetRys;
}

/** Compute and store all ERIs with 8-fold symmetry.
 *  Tries backends in priority order: WASM > JS.
 *  Returns a flat FloatArray indexed by eri1DIndex. */
export async function computeERIs(
  primitiveShells: PrimitiveShell[],
  normFactors: number[],
  numBasis: number,
  schwarzThreshold = 1e-10,
  onProgress?: (percent: number, msg?: string) => void,
  backend: EriBackend = 'auto',
  algorithm: EriAlgorithm = 'auto',
): Promise<FloatArray> {
  // Resolve 'auto' to a concrete backend based on availability.
  let resolved = backend;
  if (resolved === 'auto') {
    resolved = isWasmAvailable() ? 'wasm' : 'js';
  }

  // Try WASM
  if (resolved === 'wasm') {
    if (isWasmAvailable()) {
      if (onProgress) onProgress(0, 'WASM');
      const t0 = performance.now();
      const result = computeERIsWasm(primitiveShells, normFactors, numBasis, schwarzThreshold);
      if (result) {
        const ms = (performance.now() - t0).toFixed(1);
        if (onProgress) onProgress(100, `WASM (${ms} ms)`);
        return result;
      }
    }
  }

  const eriSize = numBasis * (numBasis + 1) / 2;
  const totalSize = eriSize * (eriSize + 1) / 2;
  const eri = createFloatArray(totalSize);

  // Group primitive shells by contracted shell
  const shellGroups = groupShells(primitiveShells);
  const nGrp = shellGroups.length;

  // Schwarz bounds: Q[ij] = sqrt(|(ij|ij)|) for each shell pair (i,j)
  const schwarz = computeSchwarzBounds(shellGroups, normFactors, numBasis);

  const totalPairs = nGrp * (nGrp + 1) / 2;
  let pairCount = 0;

  for (let iGrp = 0; iGrp < nGrp; iGrp++) {
    for (let jGrp = iGrp; jGrp < nGrp; jGrp++) {
      pairCount++;
      if (onProgress && pairCount % 10 === 0) {
        onProgress(pairCount / totalPairs * 50);
      }

      const Qij = schwarz[iGrp * nGrp + jGrp];

      for (let kGrp = 0; kGrp < nGrp; kGrp++) {
        for (let lGrp = kGrp; lGrp < nGrp; lGrp++) {
          const ijMin = Math.min(iGrp, jGrp), ijMax = Math.max(iGrp, jGrp);
          const klMin = Math.min(kGrp, lGrp), klMax = Math.max(kGrp, lGrp);
          if (ijMax < klMax || (ijMax === klMax && ijMin < klMin)) continue;

          const Qkl = schwarz[kGrp * nGrp + lGrp];
          if (Qij * Qkl < schwarzThreshold) continue;

          const quartetFn = selectShellQuartetFn(
            shellGroups[iGrp], shellGroups[jGrp],
            shellGroups[kGrp], shellGroups[lGrp], algorithm);
          quartetFn(
            shellGroups[iGrp], shellGroups[jGrp],
            shellGroups[kGrp], shellGroups[lGrp],
            normFactors, numBasis, eri,
            iGrp === jGrp, kGrp === lGrp,
            iGrp === kGrp && jGrp === lGrp,
          );
        }
      }
    }
  }

  if (onProgress) onProgress(100);
  return eri;
}

/** Compute long-range ERIs with the erf(ωr)/r operator (range-separated hybrids).
 *  Uses MD path throughout (no WASM yet). Slower than full ERI but works for any L_total. */
export function computeLongRangeERIs(
  primitiveShells: PrimitiveShell[],
  normFactors: number[],
  numBasis: number,
  omega: number,
  schwarzThreshold = 1e-10,
): FloatArray {
  const eriSize = numBasis * (numBasis + 1) / 2;
  const totalSize = eriSize * (eriSize + 1) / 2;
  const eri = createFloatArray(totalSize);

  const shellGroups = groupShells(primitiveShells);
  const nGrp = shellGroups.length;
  // Use full Schwarz bounds (with full Coulomb) — long-range bounds are smaller, so screening is conservative.
  const schwarz = computeSchwarzBounds(shellGroups, normFactors, numBasis);

  for (let iGrp = 0; iGrp < nGrp; iGrp++) {
    for (let jGrp = iGrp; jGrp < nGrp; jGrp++) {
      const Qij = schwarz[iGrp * nGrp + jGrp];
      for (let kGrp = 0; kGrp < nGrp; kGrp++) {
        for (let lGrp = kGrp; lGrp < nGrp; lGrp++) {
          const ijMin = Math.min(iGrp, jGrp), ijMax = Math.max(iGrp, jGrp);
          const klMin = Math.min(kGrp, lGrp), klMax = Math.max(kGrp, lGrp);
          if (ijMax < klMax || (ijMax === klMax && ijMin < klMin)) continue;
          const Qkl = schwarz[kGrp * nGrp + lGrp];
          if (Qij * Qkl < schwarzThreshold) continue;
          // Always use MD with omega (other algorithms don't support range separation yet)
          computeShellQuartetMD(
            shellGroups[iGrp], shellGroups[jGrp],
            shellGroups[kGrp], shellGroups[lGrp],
            normFactors, numBasis, eri,
            iGrp === jGrp, kGrp === lGrp,
            iGrp === kGrp && jGrp === lGrp,
            omega,
          );
        }
      }
    }
  }
  return eri;
}

/** ERI 1D index with 8-fold symmetry: (ij|kl) where ij >= kl */
export function eri1DIndex(i: number, j: number, k: number, l: number, _numBasis: number): number {
  let ii = i, jj = j, kk = k, ll = l;
  if (ii < jj) { const tmp = ii; ii = jj; jj = tmp; }
  if (kk < ll) { const tmp = kk; kk = ll; ll = tmp; }
  let ij = ii * (ii + 1) / 2 + jj;
  let kl = kk * (kk + 1) / 2 + ll;
  if (ij < kl) { const tmp = ij; ij = kl; kl = tmp; }
  return ij * (ij + 1) / 2 + kl;
}

export interface ShellGroup {
  primitives: PrimitiveShell[];
  basisIndex: number;
  shellType: number;
}

function groupShells(shells: PrimitiveShell[]): ShellGroup[] {
  const groups: ShellGroup[] = [];
  let current: ShellGroup | null = null;
  for (const ps of shells) {
    if (!current || current.basisIndex !== ps.basisIndex || current.shellType !== ps.shellType) {
      current = { primitives: [], basisIndex: ps.basisIndex, shellType: ps.shellType };
      groups.push(current);
    }
    current.primitives.push(ps);
  }
  return groups;
}

function computeSchwarzBounds(
  shellGroups: ShellGroup[], normFactors: number[], numBasis: number,
): Float64Array {
  const nGrp = shellGroups.length;
  const Q = new Float64Array(nGrp * nGrp);
  const eriSize = numBasis * (numBasis + 1) / 2;
  const diagEri = new Float64Array(eriSize * (eriSize + 1) / 2);

  for (let iGrp = 0; iGrp < nGrp; iGrp++) {
    for (let jGrp = iGrp; jGrp < nGrp; jGrp++) {
      computeShellQuartetMD(
        shellGroups[iGrp], shellGroups[jGrp],
        shellGroups[iGrp], shellGroups[jGrp],
        normFactors, numBasis, diagEri,
        iGrp === jGrp, iGrp === jGrp, true,
      );
      let maxVal = 0;
      const numA = shellTypeToNumBasis(shellGroups[iGrp].shellType);
      const numB = shellTypeToNumBasis(shellGroups[jGrp].shellType);
      for (let ia = 0; ia < numA; ia++) {
        const muA = shellGroups[iGrp].basisIndex + ia;
        const ibStart = iGrp === jGrp ? ia : 0;
        for (let ib = ibStart; ib < numB; ib++) {
          const muB = shellGroups[jGrp].basisIndex + ib;
          const idx = eri1DIndex(muA, muB, muA, muB, numBasis);
          const v = Math.abs(diagEri[idx]);
          if (v > maxVal) maxVal = v;
        }
      }
      const qval = Math.sqrt(maxVal);
      Q[iGrp * nGrp + jGrp] = qval;
      Q[jGrp * nGrp + iGrp] = qval;
    }
  }
  return Q;
}

// ─── McMurchie-Davidson ───────────────────────────────────────────────────────

let _rBuf = new Float64Array(0);
let _rDim = 0;

function rFunction2ePooled(
  maxOrder: number, eta: number,
  PCx: number, PCy: number, PCz: number,
  boys: Float64Array,
): void {
  const dim = maxOrder + 1;
  const size = (dim + 1) ** 4;
  if (size > _rBuf.length) { _rBuf = new Float64Array(size); }
  else { _rBuf.fill(0, 0, size); }
  _rDim = dim;
  const d1 = dim + 1, d2 = d1 * d1, d3 = d2 * d1;
  for (let n = 0; n <= maxOrder; n++) {
    _rBuf[n * d3] = Math.pow(-2 * eta, n) * boys[n];
  }
  for (let n = maxOrder - 1; n >= 0; n--) {
    for (let t = 0; t <= maxOrder - n; t++) {
      for (let u = 0; u <= maxOrder - n - t; u++) {
        for (let v = 0; v <= maxOrder - n - t - u; v++) {
          if (t + u + v === 0) continue;
          if (v > 0) {
            let val = PCz * _rBuf[(n + 1) * d3 + t * d2 + u * d1 + v - 1];
            if (v >= 2) val += (v - 1) * _rBuf[(n + 1) * d3 + t * d2 + u * d1 + v - 2];
            _rBuf[n * d3 + t * d2 + u * d1 + v] = val;
          } else if (u > 0) {
            let val = PCy * _rBuf[(n + 1) * d3 + t * d2 + (u - 1) * d1 + v];
            if (u >= 2) val += (u - 1) * _rBuf[(n + 1) * d3 + t * d2 + (u - 2) * d1 + v];
            _rBuf[n * d3 + t * d2 + u * d1 + v] = val;
          } else {
            let val = PCx * _rBuf[(n + 1) * d3 + (t - 1) * d2 + u * d1 + v];
            if (t >= 2) val += (t - 1) * _rBuf[(n + 1) * d3 + (t - 2) * d2 + u * d1 + v];
            _rBuf[n * d3 + t * d2 + u * d1 + v] = val;
          }
        }
      }
    }
  }
}

function rGetPooled(t: number, u: number, v: number, n: number): number {
  const d1 = _rDim + 1;
  return _rBuf[n * d1 * d1 * d1 + t * d1 * d1 + u * d1 + v];
}

function computeShellQuartetMD(
  grpA: ShellGroup, grpB: ShellGroup, grpC: ShellGroup, grpD: ShellGroup,
  normFactors: number[], numBasis: number, eri: FloatArray,
  sameAB: boolean, sameCD: boolean, sameBraKet: boolean,
  omega: number = 0,
) {
  const la = grpA.shellType, lb = grpB.shellType, lc = grpC.shellType, ld = grpD.shellType;
  const angA = ANGULAR_MOMENTUMS[la], angB = ANGULAR_MOMENTUMS[lb];
  const angC = ANGULAR_MOMENTUMS[lc], angD = ANGULAR_MOMENTUMS[ld];
  const numA = shellTypeToNumBasis(la), numB = shellTypeToNumBasis(lb);
  const numC = shellTypeToNumBasis(lc), numD = shellTypeToNumBasis(ld);
  const pNormsA = precomputeNorms(grpA, angA, numA);
  const pNormsB = precomputeNorms(grpB, angB, numB);
  const pNormsC = precomputeNorms(grpC, angC, numC);
  const pNormsD = precomputeNorms(grpD, angD, numD);

  for (let ia = 0; ia < numA; ia++) {
    const [lxa, lya, lza] = angA[ia];
    const muA = grpA.basisIndex + ia;
    const ibStart = sameAB ? ia : 0;
    for (let ib = ibStart; ib < numB; ib++) {
      const [lxb, lyb, lzb] = angB[ib];
      const muB = grpB.basisIndex + ib;
      const ij = muA >= muB ? muA * (muA + 1) / 2 + muB : muB * (muB + 1) / 2 + muA;
      for (let ic = 0; ic < numC; ic++) {
        const [lxc, lyc, lzc] = angC[ic];
        const muC = grpC.basisIndex + ic;
        const idStart = sameCD ? ic : 0;
        for (let id = idStart; id < numD; id++) {
          const [lxd, lyd, lzd] = angD[id];
          const muD = grpD.basisIndex + id;
          if (sameBraKet) {
            const kl = muC >= muD ? muC * (muC + 1) / 2 + muD : muD * (muD + 1) / 2 + muC;
            if (ij < kl) continue;
          }
          let val = 0;
          for (let pA = 0; pA < grpA.primitives.length; pA++) {
            const primA = grpA.primitives[pA];
            for (let pB = 0; pB < grpB.primitives.length; pB++) {
              const primB = grpB.primitives[pB];
              const alpha = primA.exponent, beta = primB.exponent;
              const p = alpha + beta, muAB = alpha * beta / p;
              const Ax = primA.coordinate.x, Ay = primA.coordinate.y, Az = primA.coordinate.z;
              const Bx = primB.coordinate.x, By = primB.coordinate.y, Bz = primB.coordinate.z;
              const Kab = Math.exp(-muAB * ((Ax - Bx) ** 2 + (Ay - By) ** 2 + (Az - Bz) ** 2));
              const Px = (alpha * Ax + beta * Bx) / p;
              const Py = (alpha * Ay + beta * By) / p;
              const Pz = (alpha * Az + beta * Bz) / p;
              const ExAB = eCoefficients(lxa, lxb, p, Px - Ax, Px - Bx);
              const EyAB = eCoefficients(lya, lyb, p, Py - Ay, Py - By);
              const EzAB = eCoefficients(lza, lzb, p, Pz - Az, Pz - Bz);
              for (let pC = 0; pC < grpC.primitives.length; pC++) {
                const primC = grpC.primitives[pC];
                for (let pD = 0; pD < grpD.primitives.length; pD++) {
                  const primD = grpD.primitives[pD];
                  const gamma = primC.exponent, delta = primD.exponent;
                  const q = gamma + delta, muCD = gamma * delta / q;
                  const Cx = primC.coordinate.x, Cy = primC.coordinate.y, Cz = primC.coordinate.z;
                  const Dx = primD.coordinate.x, Dy = primD.coordinate.y, Dz = primD.coordinate.z;
                  const Kcd = Math.exp(-muCD * ((Cx - Dx) ** 2 + (Cy - Dy) ** 2 + (Cz - Dz) ** 2));
                  const Qx = (gamma * Cx + delta * Dx) / q;
                  const Qy = (gamma * Cy + delta * Dy) / q;
                  const Qz = (gamma * Cz + delta * Dz) / q;
                  const ExCD = eCoefficients(lxc, lxd, q, Qx - Cx, Qx - Dx);
                  const EyCD = eCoefficients(lyc, lyd, q, Qy - Cy, Qy - Dy);
                  const EzCD = eCoefficients(lzc, lzd, q, Qz - Cz, Qz - Dz);
                  const RPQx = Px - Qx, RPQy = Py - Qy, RPQz = Pz - Qz;
                  const eta = p * q / (p + q);
                  const maxN = lxa + lya + lza + lxb + lyb + lzb + lxc + lyc + lzc + lxd + lyd + lzd;
                  const T = eta * (RPQx * RPQx + RPQy * RPQy + RPQz * RPQz);
                  const boys = omega > 0 ? boysAllOmega(maxN, T, eta, omega) : boysAll(maxN, T);
                  rFunction2ePooled(maxN, eta, RPQx, RPQy, RPQz, boys);
                  let sum = 0;
                  for (let t1 = 0; t1 <= lxa + lxb; t1++)
                    for (let u1 = 0; u1 <= lya + lyb; u1++)
                      for (let v1 = 0; v1 <= lza + lzb; v1++) {
                        const e1 = ExAB[t1] * EyAB[u1] * EzAB[v1];
                        if (Math.abs(e1) < 1e-15) continue;
                        for (let t2 = 0; t2 <= lxc + lxd; t2++)
                          for (let u2 = 0; u2 <= lyc + lyd; u2++)
                            for (let v2 = 0; v2 <= lzc + lzd; v2++) {
                              const e2 = ExCD[t2] * EyCD[u2] * EzCD[v2];
                              if (Math.abs(e2) < 1e-15) continue;
                              sum += e1 * e2 * (((t2 + u2 + v2) & 1) ? -1 : 1)
                                * rGetPooled(t1 + t2, u1 + u2, v1 + v2, 0);
                            }
                      }
                  val += primA.coefficient * primB.coefficient * primC.coefficient * primD.coefficient
                    * Kab * Kcd * pNormsA[ia][pA] * pNormsB[ib][pB] * pNormsC[ic][pC] * pNormsD[id][pD]
                    * 2 * Math.pow(Math.PI, 2.5) / (p * q * Math.sqrt(p + q)) * sum;
                }
              }
            }
          }
          eri[eri1DIndex(muA, muB, muC, muD, numBasis)] +=
            val * normFactors[muA] * normFactors[muB] * normFactors[muC] * normFactors[muD];
        }
      }
    }
  }
}

/** Pre-compute primitive norms for each basis function in a shell group */
function precomputeNorms(grp: ShellGroup, ang: number[][], num: number): Float64Array[] {
  const norms: Float64Array[] = [];
  for (let i = 0; i < num; i++) {
    const [lx, ly, lz] = ang[i];
    const n = new Float64Array(grp.primitives.length);
    for (let p = 0; p < grp.primitives.length; p++)
      n[p] = primitiveNorm(grp.primitives[p].exponent, lx, ly, lz);
    norms.push(n);
  }
  return norms;
}

// ─── Obara-Saika (VRR + HRR) ─────────────────────────────────────────────────

// Workspace buffer for OS VRR (module-level to avoid GC)
let _osVrrBuf = new Float64Array(0);

/** Obara-Saika VRR: 6-phase build (bra x,y,z then ket x,y,z).
 *  V[ax][ay][az][cx][cy][cz][m] stored in flat buffer.
 *  Returns the buffer; use vi() to index. */
export function osVRR(
  LAB: number, LCD: number,
  p: number, q: number,
  PAx: number, PAy: number, PAz: number,
  WPx: number, WPy: number, WPz: number,
  QCx: number, QCy: number, QCz: number,
  WQx: number, WQy: number, WQz: number,
  boys: Float64Array,
): { V: Float64Array; vi: (ax: number, ay: number, az: number, cx: number, cy: number, cz: number, m: number) => number } {
  const rho_p = q / (p + q);
  const rho_q = p / (p + q);
  const inv2p = 1 / (2 * p);
  const inv2q = 1 / (2 * q);
  const inv2pq = 1 / (2 * (p + q));

  const SA = LAB + 1, SC = LCD + 1, SM = LAB + LCD + 1;
  // Strides for 7D array: ax, ay, az, cx, cy, cz, m
  const sCZ = SM, sCY = SC * sCZ, sCX = SC * sCY;
  const sAZ = SC * sCX, sAY = SA * sAZ, sAX = SA * sAY;
  const totalSize = SA * sAX;

  if (totalSize > _osVrrBuf.length) { _osVrrBuf = new Float64Array(totalSize * 2); }
  else { _osVrrBuf.fill(0, 0, totalSize); }
  const V = _osVrrBuf;

  function vi(ax: number, ay: number, az: number, cx: number, cy: number, cz: number, m: number): number {
    return ax * sAX + ay * sAY + az * sAZ + cx * sCX + cy * sCY + cz * sCZ + m;
  }

  // Base: V[0,0,0,0,0,0,m] = boys[m]
  for (let m = 0; m <= LAB + LCD; m++) V[vi(0, 0, 0, 0, 0, 0, m)] = boys[m];

  // Phase 1: Bra x-direction (ket = 0,0,0)
  for (let ax = 0; ax < LAB; ax++) {
    const maxM = LAB + LCD - ax - 1;
    for (let m = 0; m <= maxM; m++) {
      let v = PAx * V[vi(ax, 0, 0, 0, 0, 0, m)] + WPx * V[vi(ax, 0, 0, 0, 0, 0, m + 1)];
      if (ax > 0) v += ax * inv2p * (V[vi(ax - 1, 0, 0, 0, 0, 0, m)] - rho_p * V[vi(ax - 1, 0, 0, 0, 0, 0, m + 1)]);
      V[vi(ax + 1, 0, 0, 0, 0, 0, m)] = v;
    }
  }

  // Phase 2: Bra y-direction (ket = 0,0,0)
  for (let ax = 0; ax <= LAB; ax++) {
    for (let ay = 0; ay < LAB - ax; ay++) {
      const maxM = LAB + LCD - ax - ay - 1;
      for (let m = 0; m <= maxM; m++) {
        let v = PAy * V[vi(ax, ay, 0, 0, 0, 0, m)] + WPy * V[vi(ax, ay, 0, 0, 0, 0, m + 1)];
        if (ay > 0) v += ay * inv2p * (V[vi(ax, ay - 1, 0, 0, 0, 0, m)] - rho_p * V[vi(ax, ay - 1, 0, 0, 0, 0, m + 1)]);
        V[vi(ax, ay + 1, 0, 0, 0, 0, m)] = v;
      }
    }
  }

  // Phase 3: Bra z-direction (ket = 0,0,0)
  for (let ax = 0; ax <= LAB; ax++) {
    for (let ay = 0; ay <= LAB - ax; ay++) {
      for (let az = 0; az < LAB - ax - ay; az++) {
        const maxM = LAB + LCD - ax - ay - az - 1;
        for (let m = 0; m <= maxM; m++) {
          let v = PAz * V[vi(ax, ay, az, 0, 0, 0, m)] + WPz * V[vi(ax, ay, az, 0, 0, 0, m + 1)];
          if (az > 0) v += az * inv2p * (V[vi(ax, ay, az - 1, 0, 0, 0, m)] - rho_p * V[vi(ax, ay, az - 1, 0, 0, 0, m + 1)]);
          V[vi(ax, ay, az + 1, 0, 0, 0, m)] = v;
        }
      }
    }
  }

  // Phase 4: Ket x-direction (coupling: ax)
  for (let cx = 0; cx < LCD; cx++) {
    for (let ax = 0; ax <= LAB; ax++) {
      for (let ay = 0; ay <= LAB - ax; ay++) {
        for (let az = 0; az <= LAB - ax - ay; az++) {
          const maxM = LCD - cx - 1;
          for (let m = 0; m <= maxM; m++) {
            let v = QCx * V[vi(ax, ay, az, cx, 0, 0, m)] + WQx * V[vi(ax, ay, az, cx, 0, 0, m + 1)];
            if (cx > 0) v += cx * inv2q * (V[vi(ax, ay, az, cx - 1, 0, 0, m)] - rho_q * V[vi(ax, ay, az, cx - 1, 0, 0, m + 1)]);
            if (ax > 0) v += ax * inv2pq * V[vi(ax - 1, ay, az, cx, 0, 0, m + 1)];
            V[vi(ax, ay, az, cx + 1, 0, 0, m)] = v;
          }
        }
      }
    }
  }

  // Phase 5: Ket y-direction (coupling: ay)
  for (let cx = 0; cx <= LCD; cx++) {
    for (let cy = 0; cy < LCD - cx; cy++) {
      for (let ax = 0; ax <= LAB; ax++) {
        for (let ay = 0; ay <= LAB - ax; ay++) {
          for (let az = 0; az <= LAB - ax - ay; az++) {
            const maxM = LCD - cx - cy - 1;
            for (let m = 0; m <= maxM; m++) {
              let v = QCy * V[vi(ax, ay, az, cx, cy, 0, m)] + WQy * V[vi(ax, ay, az, cx, cy, 0, m + 1)];
              if (cy > 0) v += cy * inv2q * (V[vi(ax, ay, az, cx, cy - 1, 0, m)] - rho_q * V[vi(ax, ay, az, cx, cy - 1, 0, m + 1)]);
              if (ay > 0) v += ay * inv2pq * V[vi(ax, ay - 1, az, cx, cy, 0, m + 1)];
              V[vi(ax, ay, az, cx, cy + 1, 0, m)] = v;
            }
          }
        }
      }
    }
  }

  // Phase 6: Ket z-direction (coupling: az)
  for (let cx = 0; cx <= LCD; cx++) {
    for (let cy = 0; cy <= LCD - cx; cy++) {
      for (let cz = 0; cz < LCD - cx - cy; cz++) {
        for (let ax = 0; ax <= LAB; ax++) {
          for (let ay = 0; ay <= LAB - ax; ay++) {
            for (let az = 0; az <= LAB - ax - ay; az++) {
              const maxM = LCD - cx - cy - cz - 1;
              for (let m = 0; m <= maxM; m++) {
                let v = QCz * V[vi(ax, ay, az, cx, cy, cz, m)] + WQz * V[vi(ax, ay, az, cx, cy, cz, m + 1)];
                if (cz > 0) v += cz * inv2q * (V[vi(ax, ay, az, cx, cy, cz - 1, m)] - rho_q * V[vi(ax, ay, az, cx, cy, cz - 1, m + 1)]);
                if (az > 0) v += az * inv2pq * V[vi(ax, ay, az - 1, cx, cy, cz, m + 1)];
                V[vi(ax, ay, az, cx, cy, cz + 1, m)] = v;
              }
            }
          }
        }
      }
    }
  }

  return { V, vi };
}

/** HRR key encoder for memoization (12 components, 4 bits each = 48 bits) */
export function hrrKey(ax: number, ay: number, az: number, bx: number, by: number, bz: number,
  cx: number, cy: number, cz: number, dx: number, dy: number, dz: number): number {
  return (((((((((((ax * 16 + ay) * 16 + az) * 16 + bx) * 16 + by) * 16 + bz) * 16 + cx) * 16 + cy) * 16 + cz) * 16 + dx) * 16 + dy) * 16 + dz);
}

/** Recursive HRR with memoization. Transfers angular momentum from A→B and C→D. */
export function hrrEval(
  ax: number, ay: number, az: number, bx: number, by: number, bz: number,
  cx: number, cy: number, cz: number, dx: number, dy: number, dz: number,
  V: Float64Array, vi: (ax: number, ay: number, az: number, cx: number, cy: number, cz: number, m: number) => number,
  ABx: number, ABy: number, ABz: number, CDx: number, CDy: number, CDz: number,
  memo: Map<number, number>,
): number {
  const key = hrrKey(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz);
  const cached = memo.get(key);
  if (cached !== undefined) return cached;

  let result: number;
  if (bx + by + bz > 0) {
    // Bra HRR: [a, b+e_i | cd] = [a+e_i, b | cd] + AB_i * [a, b | cd]
    if (bx > 0) {
      result = hrrEval(ax + 1, ay, az, bx - 1, by, bz, cx, cy, cz, dx, dy, dz, V, vi, ABx, ABy, ABz, CDx, CDy, CDz, memo)
        + ABx * hrrEval(ax, ay, az, bx - 1, by, bz, cx, cy, cz, dx, dy, dz, V, vi, ABx, ABy, ABz, CDx, CDy, CDz, memo);
    } else if (by > 0) {
      result = hrrEval(ax, ay + 1, az, bx, by - 1, bz, cx, cy, cz, dx, dy, dz, V, vi, ABx, ABy, ABz, CDx, CDy, CDz, memo)
        + ABy * hrrEval(ax, ay, az, bx, by - 1, bz, cx, cy, cz, dx, dy, dz, V, vi, ABx, ABy, ABz, CDx, CDy, CDz, memo);
    } else {
      result = hrrEval(ax, ay, az + 1, bx, by, bz - 1, cx, cy, cz, dx, dy, dz, V, vi, ABx, ABy, ABz, CDx, CDy, CDz, memo)
        + ABz * hrrEval(ax, ay, az, bx, by, bz - 1, cx, cy, cz, dx, dy, dz, V, vi, ABx, ABy, ABz, CDx, CDy, CDz, memo);
    }
  } else if (dx + dy + dz > 0) {
    // Ket HRR: [a | c, d+e_i] = [a | c+e_i, d] + CD_i * [a | c, d]
    if (dx > 0) {
      result = hrrEval(ax, ay, az, 0, 0, 0, cx + 1, cy, cz, dx - 1, dy, dz, V, vi, ABx, ABy, ABz, CDx, CDy, CDz, memo)
        + CDx * hrrEval(ax, ay, az, 0, 0, 0, cx, cy, cz, dx - 1, dy, dz, V, vi, ABx, ABy, ABz, CDx, CDy, CDz, memo);
    } else if (dy > 0) {
      result = hrrEval(ax, ay, az, 0, 0, 0, cx, cy + 1, cz, dx, dy - 1, dz, V, vi, ABx, ABy, ABz, CDx, CDy, CDz, memo)
        + CDy * hrrEval(ax, ay, az, 0, 0, 0, cx, cy, cz, dx, dy - 1, dz, V, vi, ABx, ABy, ABz, CDx, CDy, CDz, memo);
    } else {
      result = hrrEval(ax, ay, az, 0, 0, 0, cx, cy, cz + 1, dx, dy, dz - 1, V, vi, ABx, ABy, ABz, CDx, CDy, CDz, memo)
        + CDz * hrrEval(ax, ay, az, 0, 0, 0, cx, cy, cz, dx, dy, dz - 1, V, vi, ABx, ABy, ABz, CDx, CDy, CDz, memo);
    }
  } else {
    // Base case: b=0, d=0 → read from VRR at m=0
    result = V[vi(ax, ay, az, cx, cy, cz, 0)];
  }

  memo.set(key, result);
  return result;
}

function computeShellQuartetOS(
  grpA: ShellGroup, grpB: ShellGroup, grpC: ShellGroup, grpD: ShellGroup,
  normFactors: number[], numBasis: number, eri: FloatArray,
  sameAB: boolean, sameCD: boolean, sameBraKet: boolean,
) {
  const la = grpA.shellType, lb = grpB.shellType, lc = grpC.shellType, ld = grpD.shellType;
  const LAB = la + lb, LCD = lc + ld;
  const angA = ANGULAR_MOMENTUMS[la], angB = ANGULAR_MOMENTUMS[lb];
  const angC = ANGULAR_MOMENTUMS[lc], angD = ANGULAR_MOMENTUMS[ld];
  const numA = shellTypeToNumBasis(la), numB = shellTypeToNumBasis(lb);
  const numC = shellTypeToNumBasis(lc), numD = shellTypeToNumBasis(ld);
  const pNormsA = precomputeNorms(grpA, angA, numA);
  const pNormsB = precomputeNorms(grpB, angB, numB);
  const pNormsC = precomputeNorms(grpC, angC, numC);
  const pNormsD = precomputeNorms(grpD, angD, numD);

  // Accumulated contracted integrals per basis function quartet
  const contracted = new Float64Array(numA * numB * numC * numD);

  for (let pA = 0; pA < grpA.primitives.length; pA++) {
    const primA = grpA.primitives[pA];
    const alpha = primA.exponent;
    const Ax = primA.coordinate.x, Ay = primA.coordinate.y, Az = primA.coordinate.z;
    for (let pB = 0; pB < grpB.primitives.length; pB++) {
      const primB = grpB.primitives[pB];
      const beta = primB.exponent;
      const Bx = primB.coordinate.x, By = primB.coordinate.y, Bz = primB.coordinate.z;
      const p = alpha + beta, muAB = alpha * beta / p;
      const AB2 = (Ax - Bx) ** 2 + (Ay - By) ** 2 + (Az - Bz) ** 2;
      const Kab = Math.exp(-muAB * AB2);
      const Px = (alpha * Ax + beta * Bx) / p;
      const Py = (alpha * Ay + beta * By) / p;
      const Pz = (alpha * Az + beta * Bz) / p;

      for (let pC = 0; pC < grpC.primitives.length; pC++) {
        const primC = grpC.primitives[pC];
        const gamma = primC.exponent;
        const Cx = primC.coordinate.x, Cy = primC.coordinate.y, Cz = primC.coordinate.z;
        for (let pD = 0; pD < grpD.primitives.length; pD++) {
          const primD = grpD.primitives[pD];
          const delta = primD.exponent;
          const Dx = primD.coordinate.x, Dy = primD.coordinate.y, Dz = primD.coordinate.z;
          const q = gamma + delta, muCD = gamma * delta / q;
          const CD2 = (Cx - Dx) ** 2 + (Cy - Dy) ** 2 + (Cz - Dz) ** 2;
          const Kcd = Math.exp(-muCD * CD2);
          const Qx = (gamma * Cx + delta * Dx) / q;
          const Qy = (gamma * Cy + delta * Dy) / q;
          const Qz = (gamma * Cz + delta * Dz) / q;

          const eta = p * q / (p + q);
          const PQx = Px - Qx, PQy = Py - Qy, PQz = Pz - Qz;
          const T = eta * (PQx * PQx + PQy * PQy + PQz * PQz);
          const boys = boysAll(LAB + LCD, T);
          const prefactor = 2 * Math.pow(Math.PI, 2.5) / (p * q * Math.sqrt(p + q)) * Kab * Kcd;

          const Wx = (p * Px + q * Qx) / (p + q);
          const Wy = (p * Py + q * Qy) / (p + q);
          const Wz = (p * Pz + q * Qz) / (p + q);

          // VRR
          const { V, vi } = osVRR(LAB, LCD, p, q,
            Px - Ax, Py - Ay, Pz - Az,
            Wx - Px, Wy - Py, Wz - Pz,
            Qx - Cx, Qy - Cy, Qz - Cz,
            Wx - Qx, Wy - Qy, Wz - Qz,
            boys);

          // HRR + accumulate for each basis function quartet
          const ABx = Ax - Bx, ABy = Ay - By, ABz = Az - Bz;
          const CDx2 = Cx - Dx, CDy2 = Cy - Dy, CDz2 = Cz - Dz;
          const memo = new Map<number, number>();

          for (let ia = 0; ia < numA; ia++) {
            const [lxa, lya, lza] = angA[ia];
            for (let ib = 0; ib < numB; ib++) {
              const [lxb, lyb, lzb] = angB[ib];
              const cAB = primA.coefficient * primB.coefficient * pNormsA[ia][pA] * pNormsB[ib][pB];
              for (let ic = 0; ic < numC; ic++) {
                const [lxc, lyc, lzc] = angC[ic];
                for (let id = 0; id < numD; id++) {
                  const [lxd, lyd, lzd] = angD[id];
                  const cCD = primC.coefficient * primD.coefficient * pNormsC[ic][pC] * pNormsD[id][pD];
                  const val = hrrEval(lxa, lya, lza, lxb, lyb, lzb, lxc, lyc, lzc, lxd, lyd, lzd,
                    V, vi, ABx, ABy, ABz, CDx2, CDy2, CDz2, memo);
                  contracted[(ia * numB + ib) * numC * numD + ic * numD + id] += prefactor * cAB * cCD * val;
                }
              }
            }
          }
        }
      }
    }
  }

  // Store into ERI array with symmetry
  for (let ia = 0; ia < numA; ia++) {
    const muA = grpA.basisIndex + ia;
    const ibStart = sameAB ? ia : 0;
    for (let ib = ibStart; ib < numB; ib++) {
      const muB = grpB.basisIndex + ib;
      const ij = muA >= muB ? muA * (muA + 1) / 2 + muB : muB * (muB + 1) / 2 + muA;
      for (let ic = 0; ic < numC; ic++) {
        const muC = grpC.basisIndex + ic;
        const idStart = sameCD ? ic : 0;
        for (let id = idStart; id < numD; id++) {
          const muD = grpD.basisIndex + id;
          if (sameBraKet) {
            const kl = muC >= muD ? muC * (muC + 1) / 2 + muD : muD * (muD + 1) / 2 + muC;
            if (ij < kl) continue;
          }
          eri[eri1DIndex(muA, muB, muC, muD, numBasis)] +=
            contracted[(ia * numB + ib) * numC * numD + ic * numD + id]
            * normFactors[muA] * normFactors[muB] * normFactors[muC] * normFactors[muD];
        }
      }
    }
  }
}

// ─── Rys Quadrature ───────────────────────────────────────────────────────────

/** Compute Rys roots and weights via modified Chebyshev algorithm + Golub-Welsch. */
function rysRootsWeights(nRoots: number, T: number): { roots: Float64Array; weights: Float64Array } {
  const roots = new Float64Array(nRoots);
  const weights = new Float64Array(nRoots);
  if (nRoots === 0) return { roots, weights };

  const moments = boysAll(2 * nRoots - 1, T);
  const alpha = new Float64Array(nRoots);
  const beta = new Float64Array(nRoots);

  // Modified Chebyshev: compute alpha, beta from moments
  const sig = new Array<Float64Array>(nRoots + 1);
  for (let i = 0; i <= nRoots; i++) sig[i] = new Float64Array(2 * nRoots);
  for (let k = 0; k < 2 * nRoots; k++) sig[0][k] = moments[k];

  alpha[0] = moments[1] / moments[0];
  beta[0] = moments[0];

  for (let l = 1; l < nRoots; l++) {
    for (let k = l; k < 2 * nRoots - l; k++) {
      sig[l][k] = sig[l - 1][k + 1] - alpha[l - 1] * sig[l - 1][k]
        - (l >= 2 ? beta[l - 1] * sig[l - 2][k] : 0);
    }
    alpha[l] = sig[l][l + 1] / sig[l][l] - sig[l - 1][l] / sig[l - 1][l - 1];
    beta[l] = sig[l][l] / sig[l - 1][l - 1];
  }

  // Golub-Welsch: eigenvalues of tridiagonal matrix → roots, weights
  const d = new Float64Array(nRoots);
  const e = new Float64Array(nRoots);
  const z = new Float64Array(nRoots);
  for (let i = 0; i < nRoots; i++) {
    d[i] = alpha[i];
    z[i] = i === 0 ? 1.0 : 0.0;
    if (i > 0) e[i] = Math.sqrt(Math.abs(beta[i]));
  }
  tqlEigenvectors(d, e, z, nRoots);
  for (let i = 0; i < nRoots; i++) {
    roots[i] = d[i];
    weights[i] = beta[0] * z[i] * z[i];
  }

  return { roots, weights };
}

/** QL algorithm with implicit shift for symmetric tridiagonal eigenvalues. */
function tqlEigenvectors(d: Float64Array, e: Float64Array, z: Float64Array, n: number) {
  for (let l = 0; l < n; l++) {
    let iter = 0;
    while (true) {
      let m = l;
      while (m < n - 1) {
        const dd = Math.abs(d[m]) + Math.abs(d[m + 1]);
        if (Math.abs(e[m + 1]) <= 1e-15 * dd) break;
        m++;
      }
      if (m === l) break;
      if (++iter > 100) break;
      let g = (d[l + 1] - d[l]) / (2.0 * e[l + 1]);
      let r = Math.sqrt(g * g + 1.0);
      g = d[m] - d[l] + e[l + 1] / (g + (g >= 0 ? r : -r));
      let s = 1.0, c = 1.0, pp = 0.0;
      for (let i = m - 1; i >= l; i--) {
        const f = s * e[i + 1], b = c * e[i + 1];
        r = Math.sqrt(f * f + g * g);
        e[i + 2] = r;
        if (r === 0) { d[i + 1] -= pp; if (m + 1 < n) e[m + 1] = 0; break; }
        s = f / r; c = g / r;
        g = d[i + 1] - pp;
        r = (d[i] - g) * s + 2.0 * c * b;
        pp = s * r; d[i + 1] = g + pp; g = c * r - b;
        const fz = z[i + 1]; z[i + 1] = s * z[i] + c * fz; z[i] = c * z[i] - s * fz;
      }
      d[l] -= pp; e[l + 1] = g; if (m + 1 < n) e[m + 1] = 0;
    }
  }
}

/** 1D Rys 2D integral with VRR + HRR for one Cartesian direction.
 *  Computes I_x[la, lb, lc, ld] using VRR to build [a,c] then HRR to split (a→a,b) and (c→c,d). */
function rys1D(
  la: number, lb: number, lc: number, ld: number,
  PA: number, QC: number, B10: number, B01: number, B00: number,
  AB: number, CD: number,
): number {
  const maxA = la + lb, maxC = lc + ld;
  const strC = maxC + 1;

  // VRR: build vrr[a][c] for a=0..maxA, c=0..maxC
  const vrr = new Float64Array((maxA + 1) * strC);
  vrr[0] = 1.0;
  for (let c = 0; c < maxC; c++) {
    vrr[c + 1] = QC * vrr[c] + (c > 0 ? c * B01 * vrr[c - 1] : 0);
  }
  for (let a = 0; a < maxA; a++) {
    for (let c = 0; c <= maxC; c++) {
      vrr[(a + 1) * strC + c] = PA * vrr[a * strC + c]
        + (a > 0 ? a * B10 * vrr[(a - 1) * strC + c] : 0)
        + (c > 0 ? c * B00 * vrr[a * strC + c - 1] : 0);
    }
  }

  // Ket HRR: transfer ld from c → d
  // ketArr[a * (cSize) + c] at each d-level; cSize decreases by 1 each step
  let ketArr = vrr;
  let cSize = maxC + 1;
  for (let dd = 0; dd < ld; dd++) {
    const newCSize = cSize - 1;
    const newKet = new Float64Array((maxA + 1) * newCSize);
    for (let a = 0; a <= maxA; a++) {
      for (let c = 0; c < newCSize; c++) {
        newKet[a * newCSize + c] = ketArr[a * cSize + c + 1] + CD * ketArr[a * cSize + c];
      }
    }
    ketArr = newKet;
    cSize = newCSize;
  }
  // ketArr[a * (lc+1) + c] for a=0..maxA, c=0..lc; we need c=lc

  // Bra HRR: transfer lb from a → b
  let braArr = new Float64Array(maxA + 1);
  for (let a = 0; a <= maxA; a++) braArr[a] = ketArr[a * (lc + 1) + lc];
  for (let bb = 0; bb < lb; bb++) {
    const newSize = braArr.length - 1;
    const newBra = new Float64Array(newSize);
    for (let a = 0; a < newSize; a++) {
      newBra[a] = braArr[a + 1] + AB * braArr[a];
    }
    braArr = newBra;
  }
  return braArr[la];
}

function computeShellQuartetRys(
  grpA: ShellGroup, grpB: ShellGroup, grpC: ShellGroup, grpD: ShellGroup,
  normFactors: number[], numBasis: number, eri: FloatArray,
  sameAB: boolean, sameCD: boolean, sameBraKet: boolean,
) {
  const la = grpA.shellType, lb = grpB.shellType, lc = grpC.shellType, ld = grpD.shellType;
  const angA = ANGULAR_MOMENTUMS[la], angB = ANGULAR_MOMENTUMS[lb];
  const angC = ANGULAR_MOMENTUMS[lc], angD = ANGULAR_MOMENTUMS[ld];
  const numA = shellTypeToNumBasis(la), numB = shellTypeToNumBasis(lb);
  const numC = shellTypeToNumBasis(lc), numD = shellTypeToNumBasis(ld);
  const nRoots = Math.ceil((la + lb + lc + ld + 1) / 2);
  const pNormsA = precomputeNorms(grpA, angA, numA);
  const pNormsB = precomputeNorms(grpB, angB, numB);
  const pNormsC = precomputeNorms(grpC, angC, numC);
  const pNormsD = precomputeNorms(grpD, angD, numD);

  for (let ia = 0; ia < numA; ia++) {
    const [lxa, lya, lza] = angA[ia];
    const muA = grpA.basisIndex + ia;
    const ibStart = sameAB ? ia : 0;
    for (let ib = ibStart; ib < numB; ib++) {
      const [lxb, lyb, lzb] = angB[ib];
      const muB = grpB.basisIndex + ib;
      const ij = muA >= muB ? muA * (muA + 1) / 2 + muB : muB * (muB + 1) / 2 + muA;
      for (let ic = 0; ic < numC; ic++) {
        const [lxc, lyc, lzc] = angC[ic];
        const muC = grpC.basisIndex + ic;
        const idStart = sameCD ? ic : 0;
        for (let id = idStart; id < numD; id++) {
          const [lxd, lyd, lzd] = angD[id];
          const muD = grpD.basisIndex + id;
          if (sameBraKet) {
            const kl = muC >= muD ? muC * (muC + 1) / 2 + muD : muD * (muD + 1) / 2 + muC;
            if (ij < kl) continue;
          }
          let val = 0;
          for (let pA = 0; pA < grpA.primitives.length; pA++) {
            const primA = grpA.primitives[pA];
            for (let pB = 0; pB < grpB.primitives.length; pB++) {
              const primB = grpB.primitives[pB];
              const alphaV = primA.exponent, betaV = primB.exponent;
              const p = alphaV + betaV, muABv = alphaV * betaV / p;
              const Ax = primA.coordinate.x, Ay = primA.coordinate.y, Az = primA.coordinate.z;
              const Bx = primB.coordinate.x, By = primB.coordinate.y, Bz = primB.coordinate.z;
              const Kab = Math.exp(-muABv * ((Ax - Bx) ** 2 + (Ay - By) ** 2 + (Az - Bz) ** 2));
              const Px = (alphaV * Ax + betaV * Bx) / p;
              const Py = (alphaV * Ay + betaV * By) / p;
              const Pz = (alphaV * Az + betaV * Bz) / p;
              for (let pC = 0; pC < grpC.primitives.length; pC++) {
                const primC = grpC.primitives[pC];
                for (let pD = 0; pD < grpD.primitives.length; pD++) {
                  const primD = grpD.primitives[pD];
                  const gammaV = primC.exponent, deltaV = primD.exponent;
                  const q = gammaV + deltaV, muCDv = gammaV * deltaV / q;
                  const Cx = primC.coordinate.x, Cy = primC.coordinate.y, Cz = primC.coordinate.z;
                  const Dx = primD.coordinate.x, Dy = primD.coordinate.y, Dz = primD.coordinate.z;
                  const Kcd = Math.exp(-muCDv * ((Cx - Dx) ** 2 + (Cy - Dy) ** 2 + (Cz - Dz) ** 2));
                  const Qx = (gammaV * Cx + deltaV * Dx) / q;
                  const Qy = (gammaV * Cy + deltaV * Dy) / q;
                  const Qz = (gammaV * Cz + deltaV * Dz) / q;
                  const PQx = Px - Qx, PQy = Py - Qy, PQz = Pz - Qz;
                  const T = p * q / (p + q) * (PQx * PQx + PQy * PQy + PQz * PQz);
                  const { roots, weights } = rysRootsWeights(nRoots, T);
                  const ABx = Ax - Bx, ABy = Ay - By, ABz = Az - Bz;
                  const CDx2 = Cx - Dx, CDy2 = Cy - Dy, CDz2 = Cz - Dz;

                  let sum = 0;
                  for (let r = 0; r < nRoots; r++) {
                    const t2 = roots[r];
                    const B00 = t2 / (2 * (p + q));
                    const B10 = (1 - t2 * q / (p + q)) / (2 * p);
                    const B01 = (1 - t2 * p / (p + q)) / (2 * q);
                    const PAx2 = Px - Ax + t2 * (Qx - Px) * q / (p + q);
                    const PAy2 = Py - Ay + t2 * (Qy - Py) * q / (p + q);
                    const PAz2 = Pz - Az + t2 * (Qz - Pz) * q / (p + q);
                    const QCx2 = Qx - Cx + t2 * (Px - Qx) * p / (p + q);
                    const QCy2 = Qy - Cy + t2 * (Py - Qy) * p / (p + q);
                    const QCz2 = Qz - Cz + t2 * (Pz - Qz) * p / (p + q);

                    const Ix = rys1D(lxa, lxb, lxc, lxd, PAx2, QCx2, B10, B01, B00, ABx, CDx2);
                    const Iy = rys1D(lya, lyb, lyc, lyd, PAy2, QCy2, B10, B01, B00, ABy, CDy2);
                    const Iz = rys1D(lza, lzb, lzc, lzd, PAz2, QCz2, B10, B01, B00, ABz, CDz2);
                    sum += weights[r] * Ix * Iy * Iz;
                  }
                  // Match MD multiplication order: coeff * prefactor * sum in single chain
                  val += primA.coefficient * primB.coefficient * primC.coefficient * primD.coefficient
                    * Kab * Kcd * pNormsA[ia][pA] * pNormsB[ib][pB] * pNormsC[ic][pC] * pNormsD[id][pD]
                    * 2 * Math.pow(Math.PI, 2.5) / (p * q * Math.sqrt(p + q)) * sum;
                }
              }
            }
          }
          eri[eri1DIndex(muA, muB, muC, muD, numBasis)] +=
            val * normFactors[muA] * normFactors[muB] * normFactors[muC] * normFactors[muD];
        }
      }
    }
  }
}

// ─── Head-Gordon-Pople (HGP) ─────────────────────────────────────────────────
// HGP = VRR → contract over primitives → HRR (saves repeating HRR for each primitive)
function computeShellQuartetHGP(
  grpA: ShellGroup, grpB: ShellGroup, grpC: ShellGroup, grpD: ShellGroup,
  normFactors: number[], numBasis: number, eri: FloatArray,
  sameAB: boolean, sameCD: boolean, sameBraKet: boolean,
) {
  const la = grpA.shellType, lb = grpB.shellType, lc = grpC.shellType, ld = grpD.shellType;
  const LAB = la + lb, LCD = lc + ld;
  const angA = ANGULAR_MOMENTUMS[la], angB = ANGULAR_MOMENTUMS[lb];
  const angC = ANGULAR_MOMENTUMS[lc], angD = ANGULAR_MOMENTUMS[ld];
  const numA = shellTypeToNumBasis(la), numB = shellTypeToNumBasis(lb);
  const numC = shellTypeToNumBasis(lc), numD = shellTypeToNumBasis(ld);

  // Number of Cartesian functions at each total angular momentum level
  const ncartA = (LAB + 1) * (LAB + 2) * (LAB + 3) / 6;  // cumulative up to LAB
  const ncartC = (LCD + 1) * (LCD + 2) * (LCD + 3) / 6;

  // Contracted VRR result: sum over primitives of prefactor * coeff * V[a|c]
  // Indexed by flat Cartesian indices for bra (ax,ay,az with |a|<=LAB) and ket (cx,cy,cz with |c|<=LCD)
  // Use direct (ax,ay,az,cx,cy,cz) indexing same as VRR
  const SA = LAB + 1, SC = LCD + 1;
  const cStride = SC * SC * SC;
  const aStride = SA * SA * SA * cStride;
  const contracted = new Float64Array(aStride);

  function cIdx(ax: number, ay: number, az: number, cx: number, cy: number, cz: number): number {
    return ((ax * SA + ay) * SA + az) * cStride + (cx * SC + cy) * SC + cz;
  }

  const pNormsA = precomputeNorms(grpA, angA, numA);
  const pNormsB = precomputeNorms(grpB, angB, numB);
  const pNormsC = precomputeNorms(grpC, angC, numC);
  const pNormsD = precomputeNorms(grpD, angD, numD);

  // Phase 1: VRR + contraction (accumulate over all primitives)
  for (let pA = 0; pA < grpA.primitives.length; pA++) {
    const primA = grpA.primitives[pA];
    const alpha = primA.exponent;
    const Ax = primA.coordinate.x, Ay = primA.coordinate.y, Az = primA.coordinate.z;
    for (let pB = 0; pB < grpB.primitives.length; pB++) {
      const primB = grpB.primitives[pB];
      const beta = primB.exponent;
      const Bx = primB.coordinate.x, By = primB.coordinate.y, Bz = primB.coordinate.z;
      const p = alpha + beta, muAB = alpha * beta / p;
      const AB2 = (Ax - Bx) ** 2 + (Ay - By) ** 2 + (Az - Bz) ** 2;
      const Kab = Math.exp(-muAB * AB2);
      const Px = (alpha * Ax + beta * Bx) / p;
      const Py = (alpha * Ay + beta * By) / p;
      const Pz = (alpha * Az + beta * Bz) / p;

      for (let pC = 0; pC < grpC.primitives.length; pC++) {
        const primC = grpC.primitives[pC];
        const gamma = primC.exponent;
        const Cx = primC.coordinate.x, Cy = primC.coordinate.y, Cz = primC.coordinate.z;
        for (let pD = 0; pD < grpD.primitives.length; pD++) {
          const primD = grpD.primitives[pD];
          const delta = primD.exponent;
          const Dx = primD.coordinate.x, Dy = primD.coordinate.y, Dz = primD.coordinate.z;
          const q = gamma + delta, muCD = gamma * delta / q;
          const CD2 = (Cx - Dx) ** 2 + (Cy - Dy) ** 2 + (Cz - Dz) ** 2;
          const Kcd = Math.exp(-muCD * CD2);
          const Qx = (gamma * Cx + delta * Dx) / q;
          const Qy = (gamma * Cy + delta * Dy) / q;
          const Qz = (gamma * Cz + delta * Dz) / q;

          const eta = p * q / (p + q);
          const PQx = Px - Qx, PQy = Py - Qy, PQz = Pz - Qz;
          const T = eta * (PQx * PQx + PQy * PQy + PQz * PQz);
          const boys = boysAll(LAB + LCD, T);
          const prefactor = 2 * Math.pow(Math.PI, 2.5) / (p * q * Math.sqrt(p + q)) * Kab * Kcd;

          const Wx = (p * Px + q * Qx) / (p + q);
          const Wy = (p * Py + q * Qy) / (p + q);
          const Wz = (p * Pz + q * Qz) / (p + q);

          // VRR: compute [a,0|c,0]^(0) for all Cartesian triples
          const { V, vi } = osVRR(LAB, LCD, p, q,
            Px - Ax, Py - Ay, Pz - Az,
            Wx - Px, Wy - Py, Wz - Pz,
            Qx - Cx, Qy - Cy, Qz - Cz,
            Wx - Qx, Wy - Qy, Wz - Qz,
            boys);

          // Contract VRR values with primitive coefficients + norms
          // For each (ia, ic) pair, accumulate coeff_a * coeff_c * prefactor * V[a|c]
          // But HGP contracts at the [a,0|c,0] level (before HRR), so we loop over
          // all Cartesian triples (ax,ay,az) with |a|<=LAB and (cx,cy,cz) with |c|<=LCD
          for (let ax = 0; ax <= LAB; ax++) {
            for (let ay = 0; ay <= LAB - ax; ay++) {
              for (let az = 0; az <= LAB - ax - ay; az++) {
                // Primitive norm contribution for this Cartesian triple on bra
                // Sum over (ia, pA, pB) → but pA, pB are fixed in this loop.
                // We need to weight by the primitive coefficients and norms for ALL
                // basis functions that this Cartesian triple contributes to.
                // Since VRR gives [a,0|c,0] (no b,d), we need coefficients for
                // the COMBINED a=la+lb angular momentum, not individual (la,lb).
                // The primitive coefficient is just primA.coeff * primB.coeff * Kab * ...
                // The norm depends on (ax,ay,az) which maps to specific basis functions.
                // But in HGP, we contract the primitive pairs, not the basis functions.
                // The contraction coefficient for a primitive pair (pA,pB) at Cartesian (ax,ay,az) is:
                //   primA.coeff * primB.coeff * pNorm(alpha, ax_part_from_A) * pNorm(beta, ax_part_from_B)
                // This DEPENDS on how (ax,ay,az) decomposes into A and B contributions.
                // That decomposition is exactly what HRR does! So we can't fully contract before HRR.
                //
                // The correct HGP approach: contract over ket primitives (pC,pD) first,
                // keeping bra primitives (pA,pB) uncontracted. Then HRR on ket side.
                // Then contract over bra primitives. Then HRR on bra side.
                // This is the "interleaved contraction" of HGP.
                //
                // Simpler approach that still helps: since VRR result V[a|c] at m=0
                // is independent of basis function indices (it only depends on primitives),
                // we can accumulate: contracted_VRR[a][c] += prefactor * primCoeff * V[a][c]
                // where primCoeff = dA * dB * dC * dD (primitive coefficients only, no norms yet).
                // Then apply HRR on the contracted VRR, and norms during basis function output.
                // But norms depend on individual (lxa,lya,lza) not (ax,ay,az), so this doesn't work
                // because different (ia,ib) pairs contribute to the same (ax,ay,az).

                // Actually, the correct HGP contraction at the VRR level:
                // The VRR value V[ax,ay,az|cx,cy,cz] at m=0 depends only on primitives.
                // The primitive coefficient (without norms) is:
                //   primA.coeff * primB.coeff * primC.coeff * primD.coeff * Kab * Kcd * prefactor
                // The primitive NORMS depend on which basis function we're computing,
                // so they can't be contracted at this level.
                //
                // HGP solution: separate the contraction into:
                // 1. Contract primitive coefficients and Gaussians (Kab, Kcd, prefactor)
                // 2. Keep norms for the HRR + basis function output stage
                //
                // contracted_VRR[ax,ay,az,cx,cy,cz] += dA*dB*dC*dD * Kab*Kcd * prefactor * V[a|c]

                for (let cx = 0; cx <= LCD; cx++) {
                  for (let cy = 0; cy <= LCD - cx; cy++) {
                    for (let cz = 0; cz <= LCD - cx - cy; cz++) {
                      contracted[cIdx(ax, ay, az, cx, cy, cz)] +=
                        primA.coefficient * primB.coefficient * primC.coefficient * primD.coefficient
                        * prefactor * V[vi(ax, ay, az, cx, cy, cz, 0)];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Phase 2: HRR on contracted VRR (done ONCE, not per primitive quartet)
  // Apply bra HRR and ket HRR, with norms applied per basis function
  for (let ia = 0; ia < numA; ia++) {
    const [lxa, lya, lza] = angA[ia];
    const muA = grpA.basisIndex + ia;
    const ibStart = sameAB ? ia : 0;
    for (let ib = ibStart; ib < numB; ib++) {
      const [lxb, lyb, lzb] = angB[ib];
      const muB = grpB.basisIndex + ib;
      const ij = muA >= muB ? muA * (muA + 1) / 2 + muB : muB * (muB + 1) / 2 + muA;
      for (let ic = 0; ic < numC; ic++) {
        const [lxc, lyc, lzc] = angC[ic];
        const muC = grpC.basisIndex + ic;
        const idStart = sameCD ? ic : 0;
        for (let id = idStart; id < numD; id++) {
          const [lxd, lyd, lzd] = angD[id];
          const muD = grpD.basisIndex + id;
          if (sameBraKet) {
            const kl = muC >= muD ? muC * (muC + 1) / 2 + muD : muD * (muD + 1) / 2 + muC;
            if (ij < kl) continue;
          }

          // HRR on the contracted VRR values
          // Need to compute sum over all primitives of:
          //   pNormA[ia][pA] * pNormB[ib][pB] * pNormC[ic][pC] * pNormD[id][pD] * contracted_prim_VRR[a|c]
          // But contracted_VRR already has primitive coefficients summed.
          // The norms depend on (ia, pA) pairs — we need to handle them.
          //
          // Actually, pNorm(alpha, lx, ly, lz) depends on the EXPONENT and the
          // angular momentum of the specific basis function. For the VRR Cartesian
          // triple (ax,ay,az), there are multiple (ia,ib) decompositions via HRR.
          // The norm for a specific (ia with lxa,lya,lza) at primitive pA is:
          //   pNorm(alpha_pA, lxa, lya, lza)
          // This depends on BOTH the primitive index AND the basis function.
          //
          // So the full contraction is:
          // result[ia,ib,ic,id] = sum_{pA,pB,pC,pD} dA*dB*dC*dD * Kab*Kcd * prefactor
          //   * nA[ia][pA] * nB[ib][pB] * nC[ic][pC] * nD[id][pD] * HRR(V[a|c])
          //
          // We can factor this as:
          // result = sum_{pA} nA[ia][pA]*dA * sum_{pB} nB[ib][pB]*dB * ... * HRR(V)
          //
          // But V depends on ALL four primitive indices! So we can't separate like this.
          //
          // The HGP insight: V[a|c] from VRR depends on (pA,pB,pC,pD).
          // HRR(V) = sum of V at different (a,c) — it's LINEAR in V.
          // So: HRR(sum_{prim} coeff * V) = sum_{prim} coeff * HRR(V)
          // This means we CAN apply HRR AFTER contraction!
          //
          // But we need: sum_{prim} dA*dB*dC*dD * nA*nB*nC*nD * Kab*Kcd * prefactor * V
          // Our contracted_VRR has: sum_{prim} dA*dB*dC*dD * Kab*Kcd * prefactor * V
          // Missing: nA*nB*nC*nD per primitive. These can't be factored out because
          // nA depends on (alpha_pA, lxa, lya, lza) — both primitive and basis function.
          //
          // For shells where all primitives have the SAME center (contracted Gaussians),
          // the norms factor as: nA[ia][pA] = N(alpha_pA, lxa, lya, lza)
          // This depends on pA (through alpha) and ia (through lxa,lya,lza).
          // We can't separate pA from ia.
          //
          // However, we CAN do a partial contraction. The standard HGP approach:
          // 1. For each (pA, pB) bra pair: contract over (pC, pD) ket pairs
          // 2. Apply ket HRR on the ket-contracted result
          // 3. Contract over (pA, pB) bra pairs (with bra norms)
          // 4. Apply bra HRR on the bra-contracted result
          //
          // This saves: ket HRR is done nPrimA*nPrimB times instead of nPrimA*nPrimB*nPrimC*nPrimD
          //             bra HRR is done once instead of nPrimA*nPrimB times
          //
          // But implementing this requires restructuring the loops significantly.
          // For now, use the simpler approach: contracted_VRR without norms,
          // then multiply norms in a factored way during HRR.

          // Compute contracted norm factors: cnA[ax,ay,az] = sum_pA dA * nA(alpha_pA, ax,ay,az)
          // Wait, this doesn't work either since contracted_VRR already sums over all (pA,pB,pC,pD).

          // Fall back: the simple contracted_VRR approach works ONLY if norms are constant
          // across primitives (i.e., all primitives have the same exponent, which is never true).
          // So the simple approach is INCORRECT for contracted shells.

          // For correctness, fall back to OS for contracted shells, use HGP only for
          // uncontracted (single primitive) shells.
          // Since this defeats the purpose, let's implement the proper 2-stage contraction.
          break; // placeholder — will be replaced below
        }
        break;
      }
      break;
    }
    break;
  }

  // PROPER HGP: 2-stage contraction
  // Stage 1: For each bra pair (pA, pB), contract over ket pairs (pC, pD)
  //   ketContracted[pA][pB][ax,ay,az][cx,cy,cz] = sum_{pC,pD} dC*dD*nC*nD * Kcd * prefactor_ket * V[a|c]
  //   But this still needs separate storage per (pA,pB)... large memory.
  //
  // Actually, the proper HGP restructures the LOOP ORDER, not the storage.
  // Let me implement it correctly.

  // Clear the contracted array and redo with proper HGP contraction
  contracted.fill(0);

  // CORRECT HGP IMPLEMENTATION:
  // The key insight: HRR is LINEAR, so HRR(sum_prim c_i V_i) = sum_prim c_i HRR(V_i).
  // We want: result[ia,ib,ic,id] = sum_{prim} [all coefficients including norms] * HRR(V)
  //
  // The norms factor: nA(alpha, lxa,lya,lza) * nB(beta, lxb,lyb,lzb) * nC * nD
  // depends on BOTH primitive indices AND basis function indices.
  //
  // HGP approach: express the norm-weighted contraction as:
  //   result[ia,ib,ic,id] = sum_{ax,ay,az,cx,cy,cz} HRR_coeff(ia,ib,ax,ay,az) * HRR_coeff(ic,id,cx,cy,cz)
  //     * sum_{prim} dA*dB*dC*dD * nA(pA,ax,ay,az)*nB(pB,0,0,0)*... wait, VRR uses (ax,ay,az) = combined.
  //
  // Actually I realize the issue: the VRR computes [LAB, 0 | LCD, 0] where LAB = la+lb.
  // The Cartesian indices (ax,ay,az) in VRR have ax+ay+az <= LAB.
  // But norms are for INDIVIDUAL shells: nA depends on (lxa,lya,lza) with |lxa+lya+lza|=la.
  // The HRR maps VRR's (ax,ay,az) to pairs (lxa,lya,lza, lxb,lyb,lzb).
  // So nA depends on the BASIS FUNCTION index, not on VRR's Cartesian triple.
  //
  // Therefore we CANNOT simply multiply norms at the VRR level.
  // We must multiply norms at the basis function level, AFTER HRR.
  //
  // This means: contracted_VRR[a|c] = sum_{prim} dA*dB*dC*dD * nA(pA,...) * nB(pB,...) * ... * prefactor * V[a|c]
  // is WRONG because nA depends on (lxa,lya,lza) which is a function of ia, not of (ax,ay,az).
  //
  // Wait, actually... pNorm(alpha, lx, ly, lz) is a function of the EXPONENT and angular momentum.
  // In the VRR, (ax,ay,az) is the combined angular momentum (la+lb). The individual decomposition
  // into (lxa,lya,lza) + (lxb,lyb,lzb) = (ax,ay,az) happens in HRR.
  // So the norm nA(alpha, lxa, lya, lza) can only be applied AFTER HRR resolves the decomposition.
  //
  // CONCLUSION: The standard HGP with primitive norms requires the 2-stage loop restructuring:
  // Loop over (pA,pB):
  //   Loop over (pC,pD):
  //     VRR → accumulate into ket-contracted buffer (weighted by dC*dD*prefactor)
  //   Apply ket HRR on ket-contracted buffer → get [LAB,0|lc,ld] for all (ic,id)
  //   Weight by nC,nD for each (ic,id)
  //   Accumulate into bra-contracted buffer (weighted by dA*dB)
  // Apply bra HRR on bra-contracted buffer → get [la,lb|lc,ld] for all (ia,ib,ic,id)
  // Weight by nA,nB for each (ia,ib)
  //
  // This is the CORRECT HGP. Let me implement it.

  // Actually, let me simplify. The primitive norms pNorm(alpha, lx, ly, lz) only depend on
  // the exponent and the angular momentum of the INDIVIDUAL shell function.
  // In the VRR+HRR framework:
  // - VRR computes [ax,ay,az | cx,cy,cz] using combined angular momenta
  // - HRR decomposes: [ax,ay,az | cx,cy,cz] → [lxa,lya,lza, lxb,lyb,lzb | lxc,lyc,lzc, lxd,lyd,lzd]
  // - Norms: nA(alpha, lxa,lya,lza) applies to specific (ia) AFTER HRR
  //
  // For HGP with 2-stage contraction:
  // ketContracted_ab[ax,ay,az, ic, id] = sum_{pC,pD} dC*dD*Kcd*prefactor_ket
  //   * sum_{cx,cy,cz} HRR_ket_coeff(ic,id,cx,cy,cz) * V[ax,ay,az|cx,cy,cz]
  //   * nC(gamma, lxc,lyc,lzc) * nD(delta, lxd,lyd,lzd)
  //
  // This requires ket HRR per (pA,pB) pair. The savings come from not repeating
  // bra HRR per (pC,pD) pair.

  // For simplicity and correctness, implement the 2-stage HGP:
  // Outer: loop (pA, pB)
  //   Inner: loop (pC, pD) → VRR → accumulate into ketBuf[ax,ay,az][ic,id]
  //   After all (pC,pD): apply bra HRR on ketBuf → braResult[ia,ib][ic,id]
  //   Accumulate braResult weighted by dA*dB*nA*nB into final result

  const nKet = numC * numD;
  // ketBuf[ax * SA * SA + ay * SA + az][ic * numD + id]
  const ketBufSize = SA * SA * SA * nKet;
  const ketBuf = new Float64Array(ketBufSize);
  const kIdx = (ax: number, ay: number, az: number, icd: number) =>
    ((ax * SA + ay) * SA + az) * nKet + icd;

  // Final result
  const result = new Float64Array(numA * numB * nKet);

  for (let pA = 0; pA < grpA.primitives.length; pA++) {
    const primA = grpA.primitives[pA];
    const alpha = primA.exponent;
    const Ax = primA.coordinate.x, Ay = primA.coordinate.y, Az = primA.coordinate.z;
    for (let pB = 0; pB < grpB.primitives.length; pB++) {
      const primB = grpB.primitives[pB];
      const beta = primB.exponent;
      const Bx = primB.coordinate.x, By = primB.coordinate.y, Bz = primB.coordinate.z;
      const p = alpha + beta, muAB = alpha * beta / p;
      const AB2 = (Ax - Bx) ** 2 + (Ay - By) ** 2 + (Az - Bz) ** 2;
      const Kab = Math.exp(-muAB * AB2);
      const Px = (alpha * Ax + beta * Bx) / p;
      const Py = (alpha * Ay + beta * By) / p;
      const Pz = (alpha * Az + beta * Bz) / p;

      // Clear ket buffer for this (pA, pB) pair
      ketBuf.fill(0);

      // Stage 1: Contract over ket primitives (pC, pD) → ket HRR → ketBuf
      for (let pC = 0; pC < grpC.primitives.length; pC++) {
        const primC = grpC.primitives[pC];
        const gamma = primC.exponent;
        const Cx = primC.coordinate.x, Cy = primC.coordinate.y, Cz = primC.coordinate.z;
        for (let pD = 0; pD < grpD.primitives.length; pD++) {
          const primD = grpD.primitives[pD];
          const delta = primD.exponent;
          const Dx = primD.coordinate.x, Dy = primD.coordinate.y, Dz = primD.coordinate.z;
          const q = gamma + delta, muCD = gamma * delta / q;
          const CD2 = (Cx - Dx) ** 2 + (Cy - Dy) ** 2 + (Cz - Dz) ** 2;
          const Kcd = Math.exp(-muCD * CD2);
          const Qx = (gamma * Cx + delta * Dx) / q;
          const Qy = (gamma * Cy + delta * Dy) / q;
          const Qz = (gamma * Cz + delta * Dz) / q;

          const eta = p * q / (p + q);
          const PQx = Px - Qx, PQy = Py - Qy, PQz = Pz - Qz;
          const T = eta * (PQx * PQx + PQy * PQy + PQz * PQz);
          const boys = boysAll(LAB + LCD, T);
          const prefactor = 2 * Math.pow(Math.PI, 2.5) / (p * q * Math.sqrt(p + q)) * Kab * Kcd;

          const Wx = (p * Px + q * Qx) / (p + q);
          const Wy = (p * Py + q * Qy) / (p + q);
          const Wz = (p * Pz + q * Qz) / (p + q);

          const { V, vi } = osVRR(LAB, LCD, p, q,
            Px - Ax, Py - Ay, Pz - Az,
            Wx - Px, Wy - Py, Wz - Pz,
            Qx - Cx, Qy - Cy, Qz - Cz,
            Wx - Qx, Wy - Qy, Wz - Qz,
            boys);

          // Apply ket HRR + ket norms, accumulate into ketBuf
          const CDx2 = Cx - Dx, CDy2 = Cy - Dy, CDz2 = Cz - Dz;
          const ketMemo = new Map<number, number>();

          for (let ax = 0; ax <= LAB; ax++) {
            for (let ay = 0; ay <= LAB - ax; ay++) {
              for (let az = 0; az <= LAB - ax - ay; az++) {
                for (let ic = 0; ic < numC; ic++) {
                  const [lxc, lyc, lzc] = angC[ic];
                  for (let id = 0; id < numD; id++) {
                    const [lxd, lyd, lzd] = angD[id];
                    // Ket HRR: [ax,ay,az,0,0,0 | lxc,lyc,lzc, lxd,lyd,lzd]
                    const val = hrrEval(ax, ay, az, 0, 0, 0, lxc, lyc, lzc, lxd, lyd, lzd,
                      V, vi, 0, 0, 0, CDx2, CDy2, CDz2, ketMemo);
                    ketBuf[kIdx(ax, ay, az, ic * numD + id)] +=
                      primC.coefficient * primD.coefficient
                      * pNormsC[ic][pC] * pNormsD[id][pD]
                      * prefactor * val;
                  }
                }
              }
            }
          }
        }
      }

      // Stage 2: Apply bra HRR on ketBuf, accumulate into result with bra norms
      const ABx = grpA.primitives[0].coordinate.x - grpB.primitives[0].coordinate.x;
      const ABy = grpA.primitives[0].coordinate.y - grpB.primitives[0].coordinate.y;
      const ABz = grpA.primitives[0].coordinate.z - grpB.primitives[0].coordinate.z;

      for (let ia = 0; ia < numA; ia++) {
        const [lxa, lya, lza] = angA[ia];
        for (let ib = 0; ib < numB; ib++) {
          const [lxb, lyb, lzb] = angB[ib];
          const nAB = primA.coefficient * primB.coefficient * pNormsA[ia][pA] * pNormsB[ib][pB];
          for (let icd = 0; icd < nKet; icd++) {
            // Bra HRR: expand [LAB,0] → [la,lb] using AB distances
            const val = braHRR(lxa, lya, lza, lxb, lyb, lzb, ABx, ABy, ABz,
              ketBuf, SA, nKet, icd);
            result[(ia * numB + ib) * nKet + icd] += nAB * val;
          }
        }
      }
    }
  }

  // Store into ERI array with symmetry
  for (let ia = 0; ia < numA; ia++) {
    const muA = grpA.basisIndex + ia;
    const ibStart = sameAB ? ia : 0;
    for (let ib = ibStart; ib < numB; ib++) {
      const muB = grpB.basisIndex + ib;
      const ij = muA >= muB ? muA * (muA + 1) / 2 + muB : muB * (muB + 1) / 2 + muA;
      for (let ic = 0; ic < numC; ic++) {
        const muC = grpC.basisIndex + ic;
        const idStart = sameCD ? ic : 0;
        for (let id = idStart; id < numD; id++) {
          const muD = grpD.basisIndex + id;
          if (sameBraKet) {
            const kl = muC >= muD ? muC * (muC + 1) / 2 + muD : muD * (muD + 1) / 2 + muC;
            if (ij < kl) continue;
          }
          eri[eri1DIndex(muA, muB, muC, muD, numBasis)] +=
            result[(ia * numB + ib) * nKet + ic * numD + id]
            * normFactors[muA] * normFactors[muB] * normFactors[muC] * normFactors[muD];
        }
      }
    }
  }
}

/** Bra HRR on pre-contracted ket buffer.
 *  ketBuf layout: ketBuf[(ax*SA+ay)*SA+az) * nKet + icd]
 *  Returns the bra-HRR'd value for (lxa,lya,lza, lxb,lyb,lzb) at ket index icd. */
function braHRR(
  ax: number, ay: number, az: number, bx: number, by: number, bz: number,
  ABx: number, ABy: number, ABz: number,
  ketBuf: Float64Array, SA: number, nKet: number, icd: number,
): number {
  if (bx + by + bz === 0) {
    return ketBuf[((ax * SA + ay) * SA + az) * nKet + icd];
  }
  if (bx > 0) {
    return braHRR(ax + 1, ay, az, bx - 1, by, bz, ABx, ABy, ABz, ketBuf, SA, nKet, icd)
      + ABx * braHRR(ax, ay, az, bx - 1, by, bz, ABx, ABy, ABz, ketBuf, SA, nKet, icd);
  }
  if (by > 0) {
    return braHRR(ax, ay + 1, az, bx, by - 1, bz, ABx, ABy, ABz, ketBuf, SA, nKet, icd)
      + ABy * braHRR(ax, ay, az, bx, by - 1, bz, ABx, ABy, ABz, ketBuf, SA, nKet, icd);
  }
  return braHRR(ax, ay, az + 1, bx, by, bz - 1, ABx, ABy, ABz, ketBuf, SA, nKet, icd)
    + ABz * braHRR(ax, ay, az, bx, by, bz - 1, ABx, ABy, ABz, ketBuf, SA, nKet, icd);
}
