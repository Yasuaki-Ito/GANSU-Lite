// McMurchie-Davidson ERI compute shader for GANSU Lite
// Each thread computes one basis function quartet within a shell quartet.

// === Constants ===
const PI: f32 = 3.14159265358979;
const LUT_NUM_XI: u32 = 1024u;
const LUT_XI_INTERVAL: f32 = 0.03125;
const LUT_K_MAX: u32 = 5u;
const LUT_N_RANGE: u32 = 24u;
const A_TR: f32 = 0.352905920120321;
const B_TR: f32 = 0.015532762923351;
const A_RS: f32 = 0.064048916778075;
const B_RS: f32 = 28.487431543672;

// === Bindings ===
struct Params {
  num_basis: u32,
  num_quartets: u32,
  schwarz_threshold: f32,
  _pad: u32,
}

// Quartet metadata: [iGrp, jGrp, kGrp, lGrp, numA, numB, numC, numD,
//                     basisA, basisB, basisC, basisD, shellTypeA, shellTypeB, shellTypeC, shellTypeD]
struct Quartet {
  i_grp: u32, j_grp: u32, k_grp: u32, l_grp: u32,
  num_a: u32, num_b: u32, num_c: u32, num_d: u32,
  basis_a: u32, basis_b: u32, basis_c: u32, basis_d: u32,
  shell_type_a: u32, shell_type_b: u32, shell_type_c: u32, shell_type_d: u32,
  same_ab: u32, same_cd: u32, same_bra_ket: u32, _pad: u32,
}

// Primitive shell: [exponent, coefficient, x, y, z, shellType, basisIndex]
struct PrimShell {
  exponent: f32, coefficient: f32,
  x: f32, y: f32, z: f32,
  shell_type: u32, basis_index: u32,
}

// Shell group index: [primStart, primCount]
struct ShellGroupInfo {
  prim_start: u32,
  prim_count: u32,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> quartets: array<Quartet>;
@group(0) @binding(2) var<storage, read> prims: array<PrimShell>;
@group(0) @binding(3) var<storage, read> shell_groups: array<ShellGroupInfo>;
@group(0) @binding(4) var<storage, read> norm_factors: array<f32>;
@group(0) @binding(5) var<storage, read> boys_lut: array<f32>;
@group(0) @binding(6) var<storage, read_write> eri_out: array<atomic<u32>>;

// === Angular momentum tables (Cartesian, up to f-type) ===
// s: (0,0,0)
// p: (1,0,0),(0,1,0),(0,0,1)
// d: (2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)
// f: (3,0,0),(0,3,0),(0,0,3),(1,2,0),(2,1,0),(2,0,1),(1,0,2),(0,1,2),(0,2,1),(1,1,1)
// Packed as 3 arrays of length 20
const ANG_LX: array<i32, 20> = array<i32, 20>(
  0,  1,0,0,  2,0,0,1,1,0,  3,0,0,1,2,2,1,0,0,1);
const ANG_LY: array<i32, 20> = array<i32, 20>(
  0,  0,1,0,  0,2,0,1,0,1,  0,3,0,2,1,0,0,1,2,1);
const ANG_LZ: array<i32, 20> = array<i32, 20>(
  0,  0,0,1,  0,0,2,0,1,1,  0,0,3,0,0,1,2,2,1,1);
// Offsets for each shell type in the angular momentum table
const ANG_OFFSET: array<u32, 5> = array<u32, 5>(0u, 1u, 4u, 10u, 20u);

fn ang_lx(shell_type: u32, idx: u32) -> i32 { return ANG_LX[ANG_OFFSET[shell_type] + idx]; }
fn ang_ly(shell_type: u32, idx: u32) -> i32 { return ANG_LY[ANG_OFFSET[shell_type] + idx]; }
fn ang_lz(shell_type: u32, idx: u32) -> i32 { return ANG_LZ[ANG_OFFSET[shell_type] + idx]; }

// === Helper functions ===

fn double_factorial(n: i32) -> f32 {
  if (n <= 1) { return 1.0; }
  var result = 1.0;
  var i = n;
  while (i >= 2) {
    result *= f32(i);
    i -= 2;
  }
  return result;
}

fn primitive_norm(alpha: f32, lx: i32, ly: i32, lz: i32) -> f32 {
  let l = lx + ly + lz;
  let dbl = double_factorial(2*lx - 1) * double_factorial(2*ly - 1) * double_factorial(2*lz - 1);
  return pow(2.0, f32(l)) / sqrt(dbl)
       * pow(2.0 / PI, 0.75)
       * pow(alpha, f32(2*l + 3) / 4.0);
}

fn eri_1d_index(i_in: u32, j_in: u32, k_in: u32, l_in: u32) -> u32 {
  var ii = i_in; var jj = j_in;
  if (ii < jj) { let t = ii; ii = jj; jj = t; }
  var kk = k_in; var ll = l_in;
  if (kk < ll) { let t = kk; kk = ll; ll = t; }
  var ij = ii * (ii + 1u) / 2u + jj;
  var kl = kk * (kk + 1u) / 2u + ll;
  if (ij < kl) { let t = ij; ij = kl; kl = t; }
  return ij * (ij + 1u) / 2u + kl;
}

// Atomic f32 add via compare-and-swap
fn atomic_add_f32(idx: u32, value: f32) {
  var old_val = atomicLoad(&eri_out[idx]);
  loop {
    let new_val = bitcast<f32>(old_val) + value;
    let result = atomicCompareExchangeWeak(&eri_out[idx], old_val, bitcast<u32>(new_val));
    if (result.exchanged) { break; }
    old_val = result.old_value;
  }
}

// === E-coefficients ===
// Returns E^{i,j}_t in a flat array. Max size: (i+j+1) entries.
// We use var<private> storage for the 3D work array.
var<private> e_work: array<f32, 128>; // (i+1)*(j+1)*(i+j+2) max: f+f = (4)(4)(8) = 128

fn e_idx(a: i32, b: i32, t: i32, nj: i32, nt: i32) -> i32 {
  return a * nj * nt + b * nt + t;
}

fn compute_e_coefficients(i: i32, j: i32, p: f32, xpa: f32, xpb: f32, result: ptr<private, array<f32, 8>>) {
  let ni = i + 1;
  let nj = j + 1;
  let max_t = i + j;
  let nt = max_t + 2;
  let size = ni * nj * nt;

  // Zero work array
  for (var k = 0; k < size; k++) { e_work[k] = 0.0; }

  e_work[e_idx(0, 0, 0, nj, nt)] = 1.0;

  // Build in a direction
  for (var a = 0; a < i; a++) {
    for (var b = 0; b <= j; b++) {
      for (var t = 0; t <= a + b + 1; t++) {
        var val = xpa * select(0.0, e_work[e_idx(a, b, t, nj, nt)], t <= a + b);
        if (t > 0) { val += (1.0 / (2.0 * p)) * e_work[e_idx(a, b, t - 1, nj, nt)]; }
        if (t <= a + b) { val += f32(t + 1) * select(0.0, e_work[e_idx(a, b, t + 1, nj, nt)], t + 1 <= a + b); }
        e_work[e_idx(a + 1, b, t, nj, nt)] += val;
      }
    }
  }

  // Build in b direction
  for (var b = 0; b < j; b++) {
    for (var t = 0; t <= i + b + 1; t++) {
      var val = xpb * select(0.0, e_work[e_idx(i, b, t, nj, nt)], t <= i + b);
      if (t > 0) { val += (1.0 / (2.0 * p)) * e_work[e_idx(i, b, t - 1, nj, nt)]; }
      if (t <= i + b) { val += f32(t + 1) * select(0.0, e_work[e_idx(i, b, t + 1, nj, nt)], t + 1 <= i + b); }
      e_work[e_idx(i, b + 1, t, nj, nt)] += val;
    }
  }

  // Extract result
  for (var t = 0; t <= max_t; t++) {
    (*result)[t] = e_work[e_idx(i, j, t, nj, nt)];
  }
}

// === Boys function ===
fn erf_approx(x: f32) -> f32 {
  let s = select(1.0, -1.0, x < 0.0);
  let t = 1.0 / (1.0 + 0.3275911 * abs(x));
  let y = 1.0 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * exp(-x * x);
  return s * y;
}

var<private> boys_vals: array<f32, 16>; // max order ~12 for s+p+d+f

// E-coefficient output buffers (private so they can be passed to compute_e_coefficients)
var<private> ex_ab: array<f32, 8>;
var<private> ey_ab: array<f32, 8>;
var<private> ez_ab: array<f32, 8>;
var<private> ex_cd: array<f32, 8>;
var<private> ey_cd: array<f32, 8>;
var<private> ez_cd: array<f32, 8>;

fn compute_boys(n: u32, x: f32) {
  if (x == 0.0) {
    for (var j = 0u; j <= n; j++) {
      boys_vals[j] = 1.0 / f32(2u * j + 1u);
    }
    return;
  }

  if (x < A_TR * f32(n) + B_TR) {
    // Taylor interpolation from LUT
    let x_idx = min(u32(round(x / LUT_XI_INTERVAL)), LUT_NUM_XI - 1u);
    let delta_x = x - LUT_XI_INTERVAL * f32(x_idx);

    for (var j = 0u; j <= n; j++) {
      var numerator = 1.0;
      var factorial = 1.0;
      var fx = boys_lut[LUT_NUM_XI * j + x_idx];

      let k_max = min(LUT_K_MAX, LUT_N_RANGE - j);
      for (var k = 1u; k <= k_max; k++) {
        numerator *= -delta_x;
        factorial *= f32(k);
        fx += boys_lut[LUT_NUM_XI * (j + k) + x_idx] * numerator / factorial;
      }
      boys_vals[j] = fx;
    }
  } else {
    // Recurrence from semi-infinite
    var exp_neg_x = 0.0f;
    let recip2x = 1.0 / (2.0 * x);
    var fx = 0.5 * sqrt(PI / x);

    if (x < A_RS * f32(n) + B_RS) {
      exp_neg_x = exp(-x);
      fx *= erf_approx(sqrt(x));
    }

    boys_vals[0] = fx;
    for (var j = 1u; j <= n; j++) {
      fx = (f32(2u * j - 1u) * fx - exp_neg_x) * recip2x;
      boys_vals[j] = fx;
    }
  }
}

// === R-function (Hermite Coulomb integrals) ===
// Max order for shells: ss=0, pp=4, dd=8, ff=12
// dim = max_order+1, d1 = dim+1, size = d1^4
// dd|dd: max_order=8, d1=10, size=10000 (f-shells use WASM path)
var<private> r_buf: array<f32, 10000>;
var<private> r_dim: u32;

fn r_flat_idx(n: u32, t: u32, u: u32, v: u32) -> u32 {
  let d1 = r_dim + 1u;
  return n * d1 * d1 * d1 + t * d1 * d1 + u * d1 + v;
}

fn compute_r_function(max_order: u32, eta: f32, pcx: f32, pcy: f32, pcz: f32) {
  let dim = max_order + 1u;
  r_dim = dim;
  let d1 = dim + 1u;
  let size = d1 * d1 * d1 * d1;

  for (var k = 0u; k < size; k++) { r_buf[k] = 0.0; }

  for (var n = 0u; n <= max_order; n++) {
    // pow() with negative base is NaN on GPU; compute (-2*eta)^n manually
    let sign_n = select(1.0, -1.0, (n & 1u) == 1u);
    r_buf[r_flat_idx(n, 0u, 0u, 0u)] = sign_n * pow(2.0 * eta, f32(n)) * boys_vals[n];
  }

  // Must iterate n from max_order-1 down to 0, but WGSL for-loop doesn't do negative steps easily
  if (max_order > 0u) {
    for (var ni = 0u; ni < max_order; ni++) {
      let n = max_order - 1u - ni; // n goes from max_order-1 down to 0
      for (var t = 0u; t <= max_order - n; t++) {
        for (var u2 = 0u; u2 <= max_order - n - t; u2++) {
          for (var v = 0u; v <= max_order - n - t - u2; v++) {
            if (t + u2 + v == 0u) { continue; }
            if (v > 0u) {
              var val = pcz * r_buf[r_flat_idx(n + 1u, t, u2, v - 1u)];
              if (v >= 2u) { val += f32(v - 1u) * r_buf[r_flat_idx(n + 1u, t, u2, v - 2u)]; }
              r_buf[r_flat_idx(n, t, u2, v)] = val;
            } else if (u2 > 0u) {
              var val = pcy * r_buf[r_flat_idx(n + 1u, t, u2 - 1u, v)];
              if (u2 >= 2u) { val += f32(u2 - 1u) * r_buf[r_flat_idx(n + 1u, t, u2 - 2u, v)]; }
              r_buf[r_flat_idx(n, t, u2, v)] = val;
            } else {
              var val = pcx * r_buf[r_flat_idx(n + 1u, t - 1u, u2, v)];
              if (t >= 2u) { val += f32(t - 1u) * r_buf[r_flat_idx(n + 1u, t - 2u, u2, v)]; }
              r_buf[r_flat_idx(n, t, u2, v)] = val;
            }
          }
        }
      }
    }
  }
}

// === Main compute kernel ===
// One workgroup per shell quartet. Each thread loops with stride to cover all basis sub-quartets.
// For (pp|pp) = 81 quartets > workgroup_size 64, threads iterate multiple times.
// 2D dispatch: quartet_idx = workgroup_id.x + workgroup_id.y * 65535

const MAX_X_WORKGROUPS: u32 = 65535u;
const WORKGROUP_SIZE: u32 = 64u;

@compute @workgroup_size(64)
fn main(
  @builtin(workgroup_id) wid: vec3u,
  @builtin(local_invocation_id) lid: vec3u,
) {
  let quartet_idx = wid.x + wid.y * MAX_X_WORKGROUPS;
  let local_idx = lid.x;

  if (quartet_idx >= params.num_quartets) { return; }

  let q = quartets[quartet_idx];
  let total_basis = q.num_a * q.num_b * q.num_c * q.num_d;

  // Get shell group info (shared across all basis sub-quartets)
  let grp_a = shell_groups[q.i_grp];
  let grp_b = shell_groups[q.j_grp];
  let grp_c = shell_groups[q.k_grp];
  let grp_d = shell_groups[q.l_grp];

  // Strided loop: each thread handles work_idx, work_idx + 64, work_idx + 128, ...
  for (var work_idx = local_idx; work_idx < total_basis; work_idx += WORKGROUP_SIZE) {

    // Decode (ia, ib, ic, id) from work_idx
    let id_val = work_idx % q.num_d;
    let ic_val = (work_idx / q.num_d) % q.num_c;
    let ib_val = (work_idx / (q.num_d * q.num_c)) % q.num_b;
    let ia_val = work_idx / (q.num_d * q.num_c * q.num_b);

    // Angular momenta
    let lxa = ang_lx(q.shell_type_a, ia_val);
    let lya = ang_ly(q.shell_type_a, ia_val);
    let lza = ang_lz(q.shell_type_a, ia_val);
    let lxb = ang_lx(q.shell_type_b, ib_val);
    let lyb = ang_ly(q.shell_type_b, ib_val);
    let lzb = ang_lz(q.shell_type_b, ib_val);
    let lxc = ang_lx(q.shell_type_c, ic_val);
    let lyc = ang_ly(q.shell_type_c, ic_val);
    let lzc = ang_lz(q.shell_type_c, ic_val);
    let lxd = ang_lx(q.shell_type_d, id_val);
    let lyd = ang_ly(q.shell_type_d, id_val);
    let lzd = ang_lz(q.shell_type_d, id_val);

    let mu_a = q.basis_a + ia_val;
    let mu_b = q.basis_b + ib_val;
    let mu_c = q.basis_c + ic_val;
    let mu_d = q.basis_d + id_val;

    // Symmetry checks — skip this sub-quartet if redundant
    let ib_start = select(0u, ia_val, q.same_ab != 0u);
    if (ib_val < ib_start) { continue; }

    let id_start = select(0u, ic_val, q.same_cd != 0u);
    if (id_val < id_start) { continue; }

    if (q.same_bra_ket != 0u) {
      let ij = select(mu_b * (mu_b + 1u) / 2u + mu_a, mu_a * (mu_a + 1u) / 2u + mu_b, mu_a >= mu_b);
      let kl = select(mu_d * (mu_d + 1u) / 2u + mu_c, mu_c * (mu_c + 1u) / 2u + mu_d, mu_c >= mu_d);
      if (ij < kl) { continue; }
    }

    var val = 0.0f;

    // Loop over primitive quartets
    for (var pa = 0u; pa < grp_a.prim_count; pa++) {
      let prim_a = prims[grp_a.prim_start + pa];
      let pnorm_a = primitive_norm(prim_a.exponent, lxa, lya, lza);

      for (var pb = 0u; pb < grp_b.prim_count; pb++) {
        let prim_b = prims[grp_b.prim_start + pb];
        let pnorm_b = primitive_norm(prim_b.exponent, lxb, lyb, lzb);

        let alpha = prim_a.exponent;
        let beta = prim_b.exponent;
        let p = alpha + beta;
        let mu_ab_val = alpha * beta / p;
        let ab2 = (prim_a.x - prim_b.x) * (prim_a.x - prim_b.x)
                + (prim_a.y - prim_b.y) * (prim_a.y - prim_b.y)
                + (prim_a.z - prim_b.z) * (prim_a.z - prim_b.z);
        let kab = exp(-mu_ab_val * ab2);
        let px = (alpha * prim_a.x + beta * prim_b.x) / p;
        let py = (alpha * prim_a.y + beta * prim_b.y) / p;
        let pz = (alpha * prim_a.z + beta * prim_b.z) / p;

        compute_e_coefficients(lxa, lxb, p, px - prim_a.x, px - prim_b.x, &ex_ab);
        compute_e_coefficients(lya, lyb, p, py - prim_a.y, py - prim_b.y, &ey_ab);
        compute_e_coefficients(lza, lzb, p, pz - prim_a.z, pz - prim_b.z, &ez_ab);

        for (var pc = 0u; pc < grp_c.prim_count; pc++) {
          let prim_c = prims[grp_c.prim_start + pc];
          let pnorm_c = primitive_norm(prim_c.exponent, lxc, lyc, lzc);

          for (var pd = 0u; pd < grp_d.prim_count; pd++) {
            let prim_d = prims[grp_d.prim_start + pd];
            let pnorm_d = primitive_norm(prim_d.exponent, lxd, lyd, lzd);

            let gamma_val = prim_c.exponent;
            let delta_val = prim_d.exponent;
            let q_val = gamma_val + delta_val;
            let mu_cd_val = gamma_val * delta_val / q_val;
            let cd2 = (prim_c.x - prim_d.x) * (prim_c.x - prim_d.x)
                     + (prim_c.y - prim_d.y) * (prim_c.y - prim_d.y)
                     + (prim_c.z - prim_d.z) * (prim_c.z - prim_d.z);
            let kcd = exp(-mu_cd_val * cd2);
            let qx = (gamma_val * prim_c.x + delta_val * prim_d.x) / q_val;
            let qy = (gamma_val * prim_c.y + delta_val * prim_d.y) / q_val;
            let qz = (gamma_val * prim_c.z + delta_val * prim_d.z) / q_val;

            compute_e_coefficients(lxc, lxd, q_val, qx - prim_c.x, qx - prim_d.x, &ex_cd);
            compute_e_coefficients(lyc, lyd, q_val, qy - prim_c.y, qy - prim_d.y, &ey_cd);
            compute_e_coefficients(lzc, lzd, q_val, qz - prim_c.z, qz - prim_d.z, &ez_cd);

            let rpqx = px - qx;
            let rpqy = py - qy;
            let rpqz = pz - qz;
            let rpq2 = rpqx * rpqx + rpqy * rpqy + rpqz * rpqz;
            let eta = p * q_val / (p + q_val);

            let max_n = u32(lxa + lya + lza + lxb + lyb + lzb + lxc + lyc + lzc + lxd + lyd + lzd);
            compute_boys(max_n, eta * rpq2);
            compute_r_function(max_n, eta, rpqx, rpqy, rpqz);

            var sum = 0.0f;
            for (var t1 = 0; t1 <= lxa + lxb; t1++) {
              for (var u1 = 0; u1 <= lya + lyb; u1++) {
                for (var v1 = 0; v1 <= lza + lzb; v1++) {
                  let e1 = ex_ab[t1] * ey_ab[u1] * ez_ab[v1];
                  if (abs(e1) < 1e-10) { continue; }
                  for (var t2 = 0; t2 <= lxc + lxd; t2++) {
                    for (var u22 = 0; u22 <= lyc + lyd; u22++) {
                      for (var v2 = 0; v2 <= lzc + lzd; v2++) {
                        let e2 = ex_cd[t2] * ey_cd[u22] * ez_cd[v2];
                        if (abs(e2) < 1e-10) { continue; }
                        let sign = select(1.0, -1.0, ((t2 + u22 + v2) & 1) == 1);
                        sum += e1 * e2 * sign * r_buf[r_flat_idx(0u, u32(t1 + t2), u32(u1 + u22), u32(v1 + v2))];
                      }
                    }
                  }
                }
              }
            }

            let coeff = prim_a.coefficient * prim_b.coefficient
                      * prim_c.coefficient * prim_d.coefficient
                      * kab * kcd
                      * pnorm_a * pnorm_b * pnorm_c * pnorm_d;
            let prefactor = 2.0 * pow(PI, 2.5) / (p * q_val * sqrt(p + q_val));

            val += coeff * prefactor * sum;
          }
        }
      }
    }

    // Apply CGTO normalization
    val *= norm_factors[mu_a] * norm_factors[mu_b] * norm_factors[mu_c] * norm_factors[mu_d];

    // Write to ERI output with atomic add
    if (abs(val) > 1e-15) {
      let idx = eri_1d_index(mu_a, mu_b, mu_c, mu_d);
      atomic_add_f32(idx, val);
    }

  } // end strided loop
}
