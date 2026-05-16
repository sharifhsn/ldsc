//! Wasmtime microbenchmark for `wasm_simd::simd128::gemm_tn_f32_slice_raw`.
//!
//! Goal: iterate on the f32 TN GEMM kernel locally without round-
//! tripping through the browser deploy. Three benches:
//!
//! - `bb_dot`: m=200, n=200, k=2490 (within-chunk B^T·B at 1000G N=2,490)
//! - `ab_dot`: m=8000, n=200, k=2490 (cross-window A^T·B at 1000G EUR)
//! - `bb_biobank`: m=200, n=200, k=50000 (within-chunk at biobank N=50K)
//!
//! Each bench reports wall-time per call + an effective GFLOPS number
//! (2*m*n*k flops). Run multiple iterations to amortise JIT warmup.
//!
//! The single-output dot kernel is structurally bound by memory
//! bandwidth on ab_dot (A is 8000 columns × 2490 rows × 4 = 80 MB;
//! way past L2). The relevant figure of merit is therefore wall-time,
//! not "theoretical SIMD throughput".

// Stub `crate::bed` module so the path-included `wasm_simd.rs` (which
// references `crate::bed::IidPos` in its `scatter` submodule) compiles
// inside this self-contained bench crate. Layout must match the
// canonical definition in `src/bed.rs::IidPos` (same field types and
// order) since `scatter_one_column_*` pointer-casts to this type.
mod bed {
    #[derive(Clone, Copy)]
    pub(crate) struct IidPos {
        pub byte_idx: usize,
        pub shift: u8,
    }
}

#[path = "../../src/wasm_simd.rs"]
mod wasm_simd;

use faer::Mat;
use std::time::Instant;

/// Splitmix64 deterministic PRNG — same as the unit tests in
/// wasm_simd.rs so the bench shapes have the same numeric character.
fn random_mat(seed: u64, rows: usize, cols: usize) -> Mat<f32> {
    let mut state = seed;
    let mut m = Mat::<f32>::zeros(rows, cols);
    for j in 0..cols {
        for i in 0..rows {
            state = state.wrapping_add(0x9E3779B97F4A7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
            z ^= z >> 31;
            let u = (z >> 40) as u32; // 24 bits
            let f = (u as f32) / ((1u32 << 24) as f32);
            m[(i, j)] = f * 2.0 - 1.0;
        }
    }
    m
}

/// Kernel variant under test.
#[derive(Clone, Copy)]
enum Kernel {
    Baseline,
    Tiled4,
    Tiled2x4,
    Tiled4x4,
    Tiled2x4Fma,
    Tiled4x4Fma,
    Tiled4x4FmaBlock,
}

impl Kernel {
    fn label(&self) -> &'static str {
        match self {
            Kernel::Baseline => "base",
            Kernel::Tiled4 => "tile4",
            Kernel::Tiled2x4 => "tile2x4",
            Kernel::Tiled4x4 => "tile4x4",
            Kernel::Tiled2x4Fma => "2x4_fma",
            Kernel::Tiled4x4Fma => "4x4_fma",
            Kernel::Tiled4x4FmaBlock => "4x4_fma_blk",
        }
    }
}

/// One bench iteration: compute `dst = lhs^T · rhs` via the kernel
/// under test. Returns the wall-time in seconds.
fn bench_one(dst: &mut Mat<f32>, lhs: &Mat<f32>, rhs: &Mat<f32>, kernel: Kernel) -> f64 {
    let m = dst.nrows();
    let n = dst.ncols();
    let k = lhs.nrows();
    debug_assert_eq!(lhs.ncols(), m);
    debug_assert_eq!(rhs.nrows(), k);
    debug_assert_eq!(rhs.ncols(), n);

    let dst_ptr = dst.as_ptr_mut();
    let dst_col = dst.col_stride();
    let a_ptr = lhs.as_ptr();
    let a_col = lhs.col_stride();
    let b_ptr = rhs.as_ptr();
    let b_col = rhs.col_stride();

    let t0 = Instant::now();
    // SAFETY: dst/lhs/rhs are all f32 col-major, row_stride==1 (faer
    // default), shapes match. Output slice is the full [0, n) range.
    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
    unsafe {
        match kernel {
            Kernel::Baseline => wasm_simd::simd128::gemm_tn_f32_slice_raw(
                dst_ptr, dst_col, m, a_ptr, a_col, k, b_ptr, b_col, 0, n,
            ),
            Kernel::Tiled4 => wasm_simd::simd128::gemm_tn_f32_slice_raw_tiled(
                dst_ptr, dst_col, m, a_ptr, a_col, k, b_ptr, b_col, 0, n,
            ),
            Kernel::Tiled2x4 => wasm_simd::simd128::gemm_tn_f32_slice_raw_tiled_2x4(
                dst_ptr, dst_col, m, a_ptr, a_col, k, b_ptr, b_col, 0, n,
            ),
            Kernel::Tiled4x4 => wasm_simd::simd128::gemm_tn_f32_slice_raw_tiled_4x4(
                dst_ptr, dst_col, m, a_ptr, a_col, k, b_ptr, b_col, 0, n,
            ),
            #[cfg(target_feature = "relaxed-simd")]
            Kernel::Tiled2x4Fma => wasm_simd::simd128::gemm_tn_f32_slice_raw_tiled_2x4_fma(
                dst_ptr, dst_col, m, a_ptr, a_col, k, b_ptr, b_col, 0, n,
            ),
            #[cfg(target_feature = "relaxed-simd")]
            Kernel::Tiled4x4Fma => wasm_simd::simd128::gemm_tn_f32_slice_raw_tiled_4x4_fma(
                dst_ptr, dst_col, m, a_ptr, a_col, k, b_ptr, b_col, 0, n,
            ),
            #[cfg(target_feature = "relaxed-simd")]
            Kernel::Tiled4x4FmaBlock => {
                wasm_simd::simd128::gemm_tn_f32_slice_raw_tiled_4x4_fma_block(
                    dst_ptr, dst_col, m, a_ptr, a_col, k, b_ptr, b_col, 0, n,
                )
            }
            #[cfg(not(target_feature = "relaxed-simd"))]
            Kernel::Tiled2x4Fma | Kernel::Tiled4x4Fma | Kernel::Tiled4x4FmaBlock => {
                panic!("FMA variants need +relaxed-simd")
            }
        }
    }
    #[cfg(not(all(target_arch = "wasm32", target_feature = "simd128")))]
    {
        // Native fallback: use the scalar reference so we can sanity-
        // check the bench harness locally before wasmtime.
        let _ = (dst_ptr, dst_col, a_ptr, a_col, b_ptr, b_col, kernel);
        wasm_simd::gemm_tn_f32_scalar(dst.as_mut(), lhs.as_ref(), rhs.as_ref());
    }
    t0.elapsed().as_secs_f64()
}

/// Run a bench shape for a given number of iterations + warmup, print
/// min / median / mean wall-times and effective GFLOPS.
fn bench_shape(
    label: &str,
    kernel: Kernel,
    m: usize,
    n: usize,
    k: usize,
    iters: usize,
    warmup: usize,
) {
    let lhs = random_mat(1, k, m);
    let rhs = random_mat(2, k, n);
    let mut dst = Mat::<f32>::zeros(m, n);

    // Warmup.
    for _ in 0..warmup {
        let _ = bench_one(&mut dst, &lhs, &rhs, kernel);
    }

    let mut times = Vec::with_capacity(iters);
    for _ in 0..iters {
        let t = bench_one(&mut dst, &lhs, &rhs, kernel);
        times.push(t);
    }
    times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let min = times[0];
    let median = times[times.len() / 2];

    // Defeat dead-code elim: read one element of dst.
    let checksum = dst[(0, 0)];

    let flops = 2.0 * (m as f64) * (n as f64) * (k as f64);
    let gflops_at = |t: f64| flops / t / 1e9;
    println!(
        "{:14} [{}] m={m:>5} n={n:>4} k={k:>5}  iters={iters:>2}  \
         min={:>8.4}s ({:>6.2} GF/s)  median={:>8.4}s ({:>6.2} GF/s)  \
         checksum={:>12.4}",
        label,
        kernel.label(),
        min,
        gflops_at(min),
        median,
        gflops_at(median),
        checksum,
    );
}

fn run_all_shapes(kernel: Kernel) {
    bench_shape("bb_dot/1000G", kernel, 200, 200, 2490, 50, 3);
    bench_shape("ab_dot/1000G", kernel, 8000, 200, 2490, 5, 1);
    bench_shape("ab_dot/1000G-mid", kernel, 2000, 200, 2490, 10, 2);
    bench_shape("bb_dot/biobank", kernel, 200, 200, 50000, 5, 1);
}

/// Compare a kernel variant's output to the baseline (single-j K-16
/// kernel) on a battery of shapes. Reports max abs error and
/// relative error; should be ≤ ~1e-4 because summation order
/// differs across variants but the operation is the same dot product
/// with f32 precision (~sqrt(k) * 1e-7 relative error).
fn verify_against_baseline(kernel: Kernel, m: usize, n: usize, k: usize) {
    let lhs = random_mat(11, k, m);
    let rhs = random_mat(22, k, n);

    let mut expected = Mat::<f32>::zeros(m, n);
    let mut got = Mat::<f32>::zeros(m, n);

    // Reference: the baseline (single-j K-16) SIMD kernel. Both
    // variants are SIMD so we compare wasm-to-wasm — apples to
    // apples on accumulation precision.
    let _ = bench_one(&mut expected, &lhs, &rhs, Kernel::Baseline);
    let _ = bench_one(&mut got, &lhs, &rhs, kernel);

    let mut max_abs = 0.0f32;
    let mut max_rel = 0.0f32;
    let mut max_loc = (0usize, 0usize);
    for j in 0..n {
        for i in 0..m {
            let e = expected[(i, j)];
            let g = got[(i, j)];
            let abs = (g - e).abs();
            let rel = abs / e.abs().max(1e-6);
            if abs > max_abs {
                max_abs = abs;
                max_loc = (i, j);
            }
            if rel > max_rel {
                max_rel = rel;
            }
        }
    }
    let pass = max_rel < 1e-3 || max_abs < 1e-3;
    println!(
        "verify [{}] m={m:>5} n={n:>4} k={k:>5}  max_abs={max_abs:>12.2e} max_rel={max_rel:>12.2e} at ({},{}) -> {}",
        kernel.label(),
        max_loc.0,
        max_loc.1,
        if pass { "PASS" } else { "FAIL" },
    );
}

// ── Scatter parallelism microbench ───────────────────────────────────
//
// Measures the per-chunk wall of `wasm_simd::scatter::scatter_one_column_f32`
// — the fused CountSketch (BED-decode → normalize → scatter-add →
// renorm) hot path that Workstream H parallelised in the browser.
// Single-threaded under wasmtime gives us the SERIAL per-chunk cost
// so we can extrapolate the post-H per-outer-Worker wall as roughly
// `serial_per_chunk × chunks_per_chr / inner_workers`.
//
// Splits c columns evenly across `n_workers` *sequential* slices to
// also measure dispatch overhead — when n_workers=1 this is just a
// straight serial run. Above 1, this approximates what one outer
// Worker's pool dispatches will look like if all inner Workers
// ran on the SAME core in sequence (so it's a *lower-bound* model
// of wall reduction; actual cores will overlap).

/// Deterministic synthetic BED packed bytes for `c` SNPs × `n_indiv`
/// individuals. Same statistical character as `make_bed_fixture` in
/// the wasm_simd tests (~6 % missing, ~31 % each real code).
fn make_synthetic_bed(seed: u64, c: usize, n_indiv: usize) -> Vec<u8> {
    let bytes_per_snp = n_indiv.div_ceil(4);
    let mut state = seed;
    let mut raw = vec![0u8; c * bytes_per_snp];
    for j in 0..c {
        for b in 0..bytes_per_snp {
            let mut byte = 0u8;
            for slot in 0..4 {
                state = state.wrapping_mul(0x9E3779B97F4A7C15).wrapping_add(1);
                let r = (state >> 32) as u32;
                let code: u8 = match r % 100 {
                    0..=5 => 0b01,
                    6..=37 => 0b00,
                    38..=69 => 0b10,
                    _ => 0b11,
                };
                byte |= code << (2 * slot);
            }
            raw[j * bytes_per_snp + b] = byte;
        }
    }
    raw
}

/// CountSketch projection state: bucket[i] ∈ [0, d), sign[i] ∈ {+1, -1}.
struct CsProj {
    bucket: Vec<u32>,
    sign: Vec<f32>,
}

fn make_cs_proj(seed: u64, n_indiv: usize, d: usize) -> CsProj {
    let mut state = seed;
    let bucket: Vec<u32> = (0..n_indiv)
        .map(|_| {
            state = state.wrapping_mul(0x9E3779B97F4A7C15).wrapping_add(1);
            ((state >> 32) as u32) % d as u32
        })
        .collect();
    let sign: Vec<f32> = (0..n_indiv)
        .map(|_| {
            state = state.wrapping_mul(0x9E3779B97F4A7C15).wrapping_add(1);
            if (state >> 63) & 1 == 0 { 1.0 } else { -1.0 }
        })
        .collect();
    CsProj { bucket, sign }
}

/// Per-SNP mean / inv_std computed the same way Pass 1 of
/// `countsketch_fused_project_f32` does. Returns parallel arrays
/// so the bench can pass raw pointers like the pool dispatch.
fn compute_stats(raw: &[u8], c: usize, n_indiv: usize) -> (Vec<f32>, Vec<f32>) {
    let bytes_per_snp = n_indiv.div_ceil(4);
    let mut mean_buf = vec![0f32; c];
    let mut inv_std_buf = vec![0f32; c];
    for j in 0..c {
        let mut sum = 0.0f64;
        let mut sum_sq = 0.0f64;
        let mut count = 0u64;
        for i in 0..n_indiv {
            let byte = raw[j * bytes_per_snp + i / 4];
            let code = (byte >> ((i % 4) * 2)) & 0b11;
            let g = match code {
                0b00 => Some(2.0_f64),
                0b01 => None,
                0b10 => Some(1.0_f64),
                _ => Some(0.0_f64),
            };
            if let Some(g) = g {
                sum += g;
                sum_sq += g * g;
                count += 1;
            }
        }
        let mean = if count > 0 { sum / count as f64 } else { 0.0 };
        let centered = sum_sq - count as f64 * mean * mean;
        let var = centered / n_indiv as f64;
        mean_buf[j] = mean as f32;
        inv_std_buf[j] = if var > 0.0 { (1.0 / var.sqrt()) as f32 } else { 0.0 };
    }
    (mean_buf, inv_std_buf)
}

/// Run `scatter_one_column_f32` over a single contiguous slice
/// `[j_start, j_end)` of output columns. Equivalent to one worker's
/// share when n_workers > 1.
#[allow(clippy::too_many_arguments)]
fn run_scatter_slice(
    j_start: usize,
    j_end: usize,
    raw: &[u8],
    bytes_per_snp: usize,
    n_indiv: usize,
    cs: &CsProj,
    n_norm: f64,
    out_ptr: *mut f32,
    out_col_stride: isize,
    d: usize,
    mean_buf: &[f32],
    inv_std_buf: &[f32],
) {
    // SAFETY: each call covers a disjoint column slice; out_ptr +
    // j*col_stride for j ∈ [j_start, j_end) is exclusively owned by
    // this slice.
    for j in j_start..j_end {
        unsafe {
            wasm_simd::scatter::scatter_one_column_f32(
                j,
                raw.as_ptr(),
                bytes_per_snp,
                n_indiv,
                cs.bucket.as_ptr(),
                cs.sign.as_ptr(),
                core::ptr::null::<bed::IidPos>(),
                0,
                true, // all_iids
                mean_buf[j],
                inv_std_buf[j],
                n_norm,
                out_ptr,
                out_col_stride,
                d,
            );
        }
    }
}

/// Time one full chunk-scatter: c output columns, n_workers serial
/// slices (when n_workers=1, one whole serial pass).
fn bench_scatter_one(
    raw: &[u8],
    n_indiv: usize,
    c: usize,
    d: usize,
    cs: &CsProj,
    mean_buf: &[f32],
    inv_std_buf: &[f32],
    n_workers: usize,
) -> f64 {
    let bytes_per_snp = n_indiv.div_ceil(4);
    let n_norm = n_indiv as f64;
    let mut out = Mat::<f32>::zeros(d, c);
    let out_ptr = out.as_mut().as_ptr_mut();
    let out_col_stride = out.col_stride();

    let chunk = c.div_ceil(n_workers);

    let t0 = Instant::now();
    for w in 0..n_workers {
        let j_start = (w * chunk).min(c);
        let j_end = ((w + 1) * chunk).min(c);
        if j_start >= j_end {
            continue;
        }
        run_scatter_slice(
            j_start,
            j_end,
            raw,
            bytes_per_snp,
            n_indiv,
            cs,
            n_norm,
            out_ptr,
            out_col_stride,
            d,
            mean_buf,
            inv_std_buf,
        );
    }
    let t = t0.elapsed().as_secs_f64();
    // Defeat dead-code elim: read one element.
    std::hint::black_box(out[(0, 0)]);
    t
}

/// Sweep `n_workers` ∈ {1, 4} for a single (n_indiv, c, d) shape.
/// `n_workers=1` is the pre-H per-chunk wall (one outer doing the
/// whole scatter); `n_workers=4` is the post-H ceiling on a 4-inner
/// pool when all slices run on the same physical core in sequence.
///
/// Notes:
/// - Wasmtime is single-threaded, so this can't measure actual
///   parallel speedup — only the SERIAL cost of c columns. The
///   reported per-slice time at `n_workers=4` is wall/4 (the time
///   each inner Worker would spend if perfectly load-balanced and
///   the kernel scales linearly with cores).
/// - Real-world speedup will land *below* this ceiling because of
///   shared L2/L3 bandwidth contention. The Workstream H expected
///   per-outer speedup is ~3-3.5× on 4 cores, so use that as the
///   realistic target instead.
fn bench_scatter_shape(label: &str, n_indiv: usize, c: usize, d: usize, iters: usize, warmup: usize) {
    let raw = make_synthetic_bed(42, c, n_indiv);
    let cs = make_cs_proj(7, n_indiv, d);
    let (mean_buf, inv_std_buf) = compute_stats(&raw, c, n_indiv);

    // Warmup at n_workers=1.
    for _ in 0..warmup {
        let _ = bench_scatter_one(&raw, n_indiv, c, d, &cs, &mean_buf, &inv_std_buf, 1);
    }

    for &n_workers in &[1usize, 4] {
        let mut times = Vec::with_capacity(iters);
        for _ in 0..iters {
            let t = bench_scatter_one(&raw, n_indiv, c, d, &cs, &mean_buf, &inv_std_buf, n_workers);
            times.push(t);
        }
        times.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let min = times[0];
        let median = times[times.len() / 2];
        let per_col = min / c as f64;
        // Per-slice ceiling: equivalent wall if all `n_workers`
        // slices ran in parallel on separate cores (linear scaling).
        let per_slice_ceiling = min / n_workers as f64;
        println!(
            "{:14} [n_workers={n_workers}] N={n_indiv:>6} c={c:>3} d={d:>4}  iters={iters:>2}  \
             min={:>8.5}s  median={:>8.5}s  per_col={:>9.3}µs  \
             ⇒ ideal_parallel_wall≈{:>8.5}s",
            label,
            min,
            median,
            per_col * 1e6,
            per_slice_ceiling,
        );
    }
}

fn run_scatter_benches() {
    println!("\n== scatter parallelism microbench (Workstream H) ==");
    println!("(wasmtime single-threaded; n_workers=4 row is the per-slice");
    println!(" ceiling = serial_wall / 4, an upper bound on what the SAB");
    println!(" pool will achieve with 4 inner Workers.)");

    // 1000G N=2,490: per-chunk scatter at the default c=200, d=200
    // sketch. Browser per-chunk is ~5-15 ms range.
    bench_scatter_shape("scatter/1000G", 2_490, 200, 200, 30, 3);
    // Biobank N=50K, d=200: the case Workstream H is built for.
    // Pre-H browser per-outer-wall on biobank is what dominates.
    bench_scatter_shape("scatter/biobank", 50_000, 200, 200, 10, 2);
    // Biobank N=50K with chunk-size 500 (CLI default for very long
    // chrs): sometimes worth checking if larger chunks shift the
    // dispatch overhead.
    bench_scatter_shape("scatter/biobank-c500", 50_000, 500, 200, 5, 1);
}

fn main() {
    println!("== ldsc wasm GEMM microbench ==");
    println!(
        "target_arch={}  simd128={}",
        std::env::consts::ARCH,
        cfg!(target_feature = "simd128"),
    );

    println!("\n-- baseline kernel (single-j, K-unroll-16, 4 accs) --");
    run_all_shapes(Kernel::Baseline);

    println!("\n-- tiled kernel (4 j's per A col load, 4 indep accs) --");
    run_all_shapes(Kernel::Tiled4);

    println!("\n-- tiled 2x4 kernel (2 i's × 4 j's per K iter, 8 accs) --");
    run_all_shapes(Kernel::Tiled2x4);

    println!("\n-- tiled 4x4 kernel (4 i's × 4 j's per K iter, 16 accs) --");
    run_all_shapes(Kernel::Tiled4x4);

    if cfg!(target_feature = "relaxed-simd") {
        println!("\n-- 2x4 + relaxed-SIMD FMA (1 op vs mul+add pair) --");
        run_all_shapes(Kernel::Tiled2x4Fma);

        println!("\n-- 4x4 + relaxed-SIMD FMA (16 accs, fused) --");
        run_all_shapes(Kernel::Tiled4x4Fma);

        println!("\n-- 4x4 FMA + ip-outer cache-block (better for big m) --");
        run_all_shapes(Kernel::Tiled4x4FmaBlock);
    }

    println!("\n-- correctness sweep (tile2x4 vs baseline) --");
    verify_against_baseline(Kernel::Tiled2x4, 200, 200, 2490);
    verify_against_baseline(Kernel::Tiled2x4, 8000, 200, 2490);
    verify_against_baseline(Kernel::Tiled2x4, 200, 200, 50000);
    // Stress shapes: odd j-tail and odd i-tail.
    verify_against_baseline(Kernel::Tiled2x4, 199, 201, 2490);
    verify_against_baseline(Kernel::Tiled2x4, 7999, 199, 2490);

    if cfg!(target_feature = "relaxed-simd") {
        println!("\n-- correctness sweep (2x4 FMA vs baseline) --");
        verify_against_baseline(Kernel::Tiled2x4Fma, 200, 200, 2490);
        verify_against_baseline(Kernel::Tiled2x4Fma, 8000, 200, 2490);
        verify_against_baseline(Kernel::Tiled2x4Fma, 200, 200, 50000);
        verify_against_baseline(Kernel::Tiled4x4Fma, 8000, 200, 2490);
        verify_against_baseline(Kernel::Tiled4x4Fma, 200, 200, 50000);
        verify_against_baseline(Kernel::Tiled4x4FmaBlock, 200, 200, 2490);
        verify_against_baseline(Kernel::Tiled4x4FmaBlock, 8000, 200, 2490);
        verify_against_baseline(Kernel::Tiled4x4FmaBlock, 200, 200, 50000);
        verify_against_baseline(Kernel::Tiled4x4FmaBlock, 199, 201, 2490);
        verify_against_baseline(Kernel::Tiled4x4FmaBlock, 7999, 199, 2490);
    }

    run_scatter_benches();

    println!("\n== done ==");
}
