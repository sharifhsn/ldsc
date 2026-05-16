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
}

impl Kernel {
    fn label(&self) -> &'static str {
        match self {
            Kernel::Baseline => "base",
            Kernel::Tiled4 => "tile4",
            Kernel::Tiled2x4 => "tile2x4",
            Kernel::Tiled4x4 => "tile4x4",
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

    println!("\n-- correctness sweep (tile2x4 vs baseline) --");
    verify_against_baseline(Kernel::Tiled2x4, 200, 200, 2490);
    verify_against_baseline(Kernel::Tiled2x4, 8000, 200, 2490);
    verify_against_baseline(Kernel::Tiled2x4, 200, 200, 50000);
    // Stress shapes: odd j-tail and odd i-tail.
    verify_against_baseline(Kernel::Tiled2x4, 199, 201, 2490);
    verify_against_baseline(Kernel::Tiled2x4, 7999, 199, 2490);

    println!("\n== done ==");
}
