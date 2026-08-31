// Fock matrix construction on GPU
// Each workgroup computes one (i,j) upper-triangle element
// Threads within workgroup split the (k,l) summation + reduction

struct Params {
  num_basis: u32,
  mode: u32,      // 0 = RHF, 1 = UHF
  _pad2: u32,
  _pad3: u32,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> eri: array<f32>;
@group(0) @binding(2) var<storage, read> density_a: array<f32>;  // RHF: P, UHF: Pa
@group(0) @binding(3) var<storage, read> density_b: array<f32>;  // UHF: Pb (RHF: dummy)
@group(0) @binding(4) var<storage, read> density_t: array<f32>;  // UHF: Ptot (RHF: dummy)
@group(0) @binding(5) var<storage, read> core_h: array<f32>;
@group(0) @binding(6) var<storage, read_write> fock_a: array<f32>;  // RHF: F, UHF: Fa
@group(0) @binding(7) var<storage, read_write> fock_b: array<f32>;  // UHF: Fb (RHF: dummy)

var<workgroup> shared_ga: array<f32, 64>;
var<workgroup> shared_gb: array<f32, 64>;

fn eri_idx(i: u32, j: u32, k: u32, l: u32) -> u32 {
  var ii = i; var jj = j;
  if (ii < jj) { let t = ii; ii = jj; jj = t; }
  var kk = k; var ll = l;
  if (kk < ll) { let t = kk; kk = ll; ll = t; }
  var ij = ii * (ii + 1u) / 2u + jj;
  var kl = kk * (kk + 1u) / 2u + ll;
  if (ij < kl) { let t = ij; ij = kl; kl = t; }
  return ij * (ij + 1u) / 2u + kl;
}

@compute @workgroup_size(64)
fn main(
  @builtin(workgroup_id) wid: vec3u,
  @builtin(local_invocation_index) tid: u32,
) {
  let n = params.num_basis;
  let pair_idx = wid.x + wid.y * 65535u;
  let total_pairs = n * (n + 1u) / 2u;
  if (pair_idx >= total_pairs) {
    return;
  }

  // Decode (i,j) from upper-triangle pair index
  // Pairs: (0,0),(0,1),...,(0,n-1),(1,1),(1,2),...
  var i = 0u;
  var j = 0u;
  var start = 0u;
  for (var ii = 0u; ii < n; ii++) {
    let row_size = n - ii;
    if (start + row_size > pair_idx) {
      i = ii;
      j = ii + (pair_idx - start);
      break;
    }
    start += row_size;
  }

  let total_kl = n * n;
  var partial_ga = 0.0f;
  var partial_gb = 0.0f;

  if (params.mode == 0u) {
    // RHF: g += P[kl] * ((ij|kl) - 0.5*(ik|jl))
    for (var kl = tid; kl < total_kl; kl += 64u) {
      let k = kl / n;
      let l = kl % n;
      let pkl = density_a[k * n + l];
      if (abs(pkl) < 1e-10f) { continue; }
      let j_int = eri[eri_idx(i, j, k, l)];
      let k_int = eri[eri_idx(i, k, j, l)];
      partial_ga += pkl * (j_int - 0.5f * k_int);
    }
  } else {
    // UHF: ga += Ptot*J - Pa*K,  gb += Ptot*J - Pb*K
    for (var kl = tid; kl < total_kl; kl += 64u) {
      let k = kl / n;
      let l = kl % n;
      let ptot = density_t[k * n + l];
      let pa = density_a[k * n + l];
      let pb = density_b[k * n + l];
      if (abs(ptot) < 1e-10f && abs(pa) < 1e-10f) { continue; }
      let j_int = eri[eri_idx(i, j, k, l)];
      let k_int = eri[eri_idx(i, k, j, l)];
      partial_ga += ptot * j_int - pa * k_int;
      partial_gb += ptot * j_int - pb * k_int;
    }
  }

  // Workgroup reduction for ga
  shared_ga[tid] = partial_ga;
  shared_gb[tid] = partial_gb;
  workgroupBarrier();

  for (var s = 32u; s > 0u; s >>= 1u) {
    if (tid < s) {
      shared_ga[tid] += shared_ga[tid + s];
      shared_gb[tid] += shared_gb[tid + s];
    }
    workgroupBarrier();
  }

  if (tid == 0u) {
    let h_ij = core_h[i * n + j];
    let va = h_ij + shared_ga[0u];
    fock_a[i * n + j] = va;
    if (i != j) {
      fock_a[j * n + i] = va;
    }

    if (params.mode == 1u) {
      let vb = h_ij + shared_gb[0u];
      fock_b[i * n + j] = vb;
      if (i != j) {
        fock_b[j * n + i] = vb;
      }
    }
  }
}
