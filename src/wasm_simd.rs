//! Hand-rolled GEMM kernels for `wasm32-unknown-unknown`.
//!
//! Why this module exists: faer 0.24's SIMD backend (`pulp`) does not
//! support wasm32 — on wasm faer falls back to a fully scalar GEMM, and
//! the resulting compute speed is the entire ~30-50× wall-clock gap
//! between the in-browser ldsc-web demo and the native CLI on full
//! 1000G EUR (see /Users/sharif/.claude/plans/yes-absolutely-make-a-zazzy-hellman.md
//! workstream F for the full plan and rationale).
//!
//! This module hosts a hand-rolled f32 TN GEMM microkernel
//! (`c = a^T · b`, all column-major) for the two hot-path call sites
//! in `l2/compute.rs::compute_ldscore_global`. Dispatch happens in
//! [`crate::la::matmul_tn_to_f32`] and is gated on
//! `cfg(all(target_arch = "wasm32", target_feature = "simd128"))`, so
//! native builds are byte-identical to before — they continue to call
//! faer's `matmul` directly.
//!
//! # Status (workstream F.2)
//!
//! Two implementations live here:
//!
//! - [`gemm_tn_f32_scalar`] — portable, safe, runs on every target.
//!   Used as the test oracle (compared against faer's `matmul` in the
//!   property tests below). NOT used on the wasm hot path anymore.
//!
//! - [`simd128::gemm_tn_f32`] — wasm32-simd128 kernel using v128
//!   intrinsics with 4 parallel `f32x4` accumulators and a K-unroll
//!   of 16. Active dispatch target on wasm32+simd128 builds.
//!
//! The two implementations are **IEEE-754 bit-identical by
//! construction**: same per-lane accumulation, same `(l0+l1)+(l2+l3)`
//! horizontal reduction, same `(h0+h1)+(h2+h3)` group combine, same
//! left-to-right scalar tail. Verified at the L2-score level by the
//! browser smoke against the CLI baseline (the scalar wasm fallback
//! already matched CLI bit-identically, and the SIMD kernel preserves
//! that match).
//!
//! # API contract
//!
//! `gemm_tn_f32_scalar(dst, lhs, rhs)` computes
//! `dst = lhs^T · rhs` with `alpha = 1.0` and replace-not-add
//! accumulation. Both call sites in compute.rs use exactly these
//! parameters; other shapes (rare on wasm — the f64 hot path and the
//! cold annotation matmuls) stay on faer.

use faer::{MatMut, MatRef};

/// Portable scalar reference for `dst = lhs^T · rhs`, all column-major,
/// f32. `alpha` is implicitly 1.0; accumulation mode is implicitly
/// `Replace` (dst is overwritten, not added to).
///
/// The reduction order is deliberately structured to match the
/// upcoming wasm32-simd128 kernel: 4 lane-groups (mimicking 4 parallel
/// `f32x4` accumulators) × an inner K-unroll of 16, with a scalar tail
/// for the trailing `k % 16` elements. Lane-internal sums proceed
/// left-to-right (matching SIMD lane order); the 4-way reduction at
/// the end uses `(a0 + a1) + (a2 + a3)` (matching the typical
/// `add_pairs(add_pairs(a0,a1), add_pairs(a2,a3))` SIMD reduction
/// pattern). When the F.2 SIMD kernel lands, scalar and SIMD outputs
/// should be bit-identical.
///
/// # Panics
///
/// Debug-asserts that `lhs.nrows() == rhs.nrows()` (the K dimension
/// matches), `lhs.ncols() == dst.nrows()` (M dimension matches), and
/// `rhs.ncols() == dst.ncols()` (N dimension matches). Release builds
/// rely on faer's bounds-checked indexing.
pub fn gemm_tn_f32_scalar(dst: MatMut<'_, f32>, lhs: MatRef<'_, f32>, rhs: MatRef<'_, f32>) {
    let m = dst.nrows();
    let n = dst.ncols();
    let k = lhs.nrows();
    debug_assert_eq!(lhs.ncols(), m, "lhs.ncols() must equal dst.nrows()");
    debug_assert_eq!(rhs.nrows(), k, "rhs.nrows() must equal lhs.nrows()");
    debug_assert_eq!(rhs.ncols(), n, "rhs.ncols() must equal dst.ncols()");

    let mut dst = dst;
    for j in 0..n {
        for i in 0..m {
            // Mimic the future SIMD kernel: 4 lane-groups, each
            // simulating an f32x4 lane-vector accumulator. Within a
            // lane-group we store 4 separate accumulators (one per
            // lane); lane-internal sums proceed left-to-right.
            let mut a0 = [0.0f32; 4];
            let mut a1 = [0.0f32; 4];
            let mut a2 = [0.0f32; 4];
            let mut a3 = [0.0f32; 4];

            let mut p = 0usize;
            while p + 16 <= k {
                for q in 0..4 {
                    a0[q] += lhs[(p + q, i)] * rhs[(p + q, j)];
                    a1[q] += lhs[(p + 4 + q, i)] * rhs[(p + 4 + q, j)];
                    a2[q] += lhs[(p + 8 + q, i)] * rhs[(p + 8 + q, j)];
                    a3[q] += lhs[(p + 12 + q, i)] * rhs[(p + 12 + q, j)];
                }
                p += 16;
            }

            // Horizontal-reduce each lane-group, then combine groups.
            // ((l0+l1)+(l2+l3)) mirrors the typical wasm v128 pairwise
            // reduction (no native horizontal-sum).
            let h0 = (a0[0] + a0[1]) + (a0[2] + a0[3]);
            let h1 = (a1[0] + a1[1]) + (a1[2] + a1[3]);
            let h2 = (a2[0] + a2[1]) + (a2[2] + a2[3]);
            let h3 = (a3[0] + a3[1]) + (a3[2] + a3[3]);
            let mut acc = (h0 + h1) + (h2 + h3);

            // Scalar tail for `k % 16` trailing elements. wasm has no
            // masked loads, so a scalar tail beats any blend trick at
            // these sizes (max 15 elements).
            while p < k {
                acc += lhs[(p, i)] * rhs[(p, j)];
                p += 1;
            }

            dst[(i, j)] = acc;
        }
    }
}

/// wasm32-simd128 kernel module. Only compiled when both the target
/// arch is `wasm32` AND the `simd128` target feature is enabled (gated
/// by `ldsc-web/.cargo/config.toml`'s `+simd128` rustflag).
///
/// Native callers see this module as empty — there are no symbols to
/// link against. Wasm callers go through
/// [`crate::la::matmul_tn_to_f32`], which dispatches into
/// [`simd128::gemm_tn_f32`] when alpha=1.0 and accum=Replace (the
/// only shape used by `l2/compute.rs::compute_ldscore_global`).
#[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
pub mod simd128 {
    use core::arch::wasm32::{
        f32x4_add, f32x4_extract_lane, f32x4_mul, f32x4_splat, v128, v128_load,
    };
    use faer::{MatMut, MatRef};

    /// In-place TN GEMM for f32 column-major matrices: `dst = lhs^T · rhs`.
    /// `alpha` is implicit 1.0; accumulation mode is implicit `Replace`.
    ///
    /// # Algorithm
    ///
    /// Single-output dot-product layout (dot kernel, not panel-tiled):
    /// for each output cell `(i, j)`, walk K with 4 parallel `f32x4`
    /// accumulators unrolled by 16. The 4 accumulators eliminate the
    /// 2-cycle dep chain on `f32x4_mul` + `f32x4_add` (wasm v128 has
    /// no FMA in baseline simd128), saturating the issue queue.
    ///
    /// Why dot-kernel and not BLIS-style register-tiled: K=2,490 keeps
    /// the inner loop compute-bound, and `B[:, j]` stays hot in L1d
    /// across all `i` iterations (10 KB at 1000G N=2,490; 200 KB at
    /// biobank N=50K, comfortably in L2). Published wasm SIMD GEMM
    /// benchmarks show only 1.3-1.6× upside from register tiling at
    /// these K sizes — not worth the ~250 LOC and a packing impl.
    ///
    /// Why scalar tail (not masked): wasm has no masked loads. A
    /// blend-via-mask trick costs more than the worst-case 15-element
    /// scalar tail.
    ///
    /// # Numerics
    ///
    /// Bit-identical to [`super::gemm_tn_f32_scalar`] by construction:
    /// same per-lane multiplications (lane q in iteration p maps to
    /// `lhs[p+q*step, i] * rhs[p+q*step, j]` in the same order), same
    /// `(l0+l1)+(l2+l3)` per-accumulator horizontal reduction, same
    /// `(h0+h1)+(h2+h3)` cross-accumulator combine, same left-to-right
    /// scalar tail.
    ///
    /// # Safety
    ///
    /// Caller must ensure:
    /// - `lhs.row_stride() == 1` (column-major, K-contiguous)
    /// - `rhs.row_stride() == 1`
    /// - `dst.row_stride() == 1`
    /// - `lhs.ncols() == dst.nrows()` (M dimension)
    /// - `rhs.nrows() == lhs.nrows()` (K dimension)
    /// - `rhs.ncols() == dst.ncols()` (N dimension)
    /// - column strides are non-negative (always true for `Mat<f32>`
    ///   and submatrices thereof; this is a faer invariant)
    ///
    /// All preconditions are debug-asserted; release builds elide them.
    /// Violation is undefined behavior — the v128 loads do raw pointer
    /// arithmetic with no bounds checks.
    ///
    /// The dispatch in `la.rs::matmul_tn_to_f32` upholds all of these:
    /// both compute.rs call sites pass column-major contiguous views
    /// produced by faer's `submatrix` / `Mat::zeros`.
    pub unsafe fn gemm_tn_f32(dst: MatMut<'_, f32>, lhs: MatRef<'_, f32>, rhs: MatRef<'_, f32>) {
        let m = dst.nrows();
        let n = dst.ncols();
        let k = lhs.nrows();
        debug_assert_eq!(lhs.ncols(), m);
        debug_assert_eq!(rhs.nrows(), k);
        debug_assert_eq!(rhs.ncols(), n);
        debug_assert_eq!(
            lhs.row_stride(),
            1,
            "lhs must be col-major (row_stride == 1)"
        );
        debug_assert_eq!(
            rhs.row_stride(),
            1,
            "rhs must be col-major (row_stride == 1)"
        );
        debug_assert_eq!(
            dst.row_stride(),
            1,
            "dst must be col-major (row_stride == 1)"
        );

        // SAFETY for the entire body: all preconditions documented on
        // this function (col-major, shapes consistent, non-negative
        // strides) are debug-asserted above and upheld by the
        // dispatch in `crate::la::matmul_tn_to_f32`. Each unsafe op
        // below is in-bounds by construction; specific reasoning is
        // inlined as `// SAFETY:` comments at the load/store sites.
        unsafe {
            let a_ptr = lhs.as_ptr();
            let a_col = lhs.col_stride();
            let b_ptr = rhs.as_ptr();
            let b_col = rhs.col_stride();

            // Extract the write pointer + leading dimension. faer's
            // `as_ptr_mut` takes `&self` in 0.24 (the `mut` in the
            // name refers to the returned pointer's mutability, not
            // the receiver), so no `mut dst` rebind is needed. The
            // MatMut is moved into this block; lifetime ensures the
            // backing storage stays valid for the kernel.
            let dst_ptr = dst.as_ptr_mut();
            let dst_col = dst.col_stride();

            let zero: v128 = f32x4_splat(0.0);

            for j in 0..n {
                // SAFETY: column j is in-bounds; b_col is the leading
                // dimension of `rhs` so j*b_col stays within rhs's
                // allocation.
                let b_col_ptr = b_ptr.offset(j as isize * b_col);

                for i in 0..m {
                    // SAFETY: column i is in-bounds; a_col is the
                    // leading dimension of `lhs`.
                    let a_col_ptr = a_ptr.offset(i as isize * a_col);

                    // 4 parallel accumulators. Each is an f32x4 SIMD
                    // vector accumulating one of 4 K-blocks per unroll
                    // iteration.
                    let mut acc0 = zero;
                    let mut acc1 = zero;
                    let mut acc2 = zero;
                    let mut acc3 = zero;

                    let mut p: usize = 0;
                    while p + 16 <= k {
                        // SAFETY: row_stride==1 (debug-asserted) so
                        // a_col_ptr.add(p+δ) for δ ∈ {0,4,8,12} maps
                        // to 16 contiguous f32 elements. p+16 ≤ k
                        // (loop guard) keeps the highest load within
                        // the column. v128_load handles unaligned
                        // addresses (faer's allocator typically aligns
                        // to ≥16, but the spec doesn't require it).
                        let a0 = v128_load(a_col_ptr.add(p) as *const v128);
                        let a1 = v128_load(a_col_ptr.add(p + 4) as *const v128);
                        let a2 = v128_load(a_col_ptr.add(p + 8) as *const v128);
                        let a3 = v128_load(a_col_ptr.add(p + 12) as *const v128);
                        let b0 = v128_load(b_col_ptr.add(p) as *const v128);
                        let b1 = v128_load(b_col_ptr.add(p + 4) as *const v128);
                        let b2 = v128_load(b_col_ptr.add(p + 8) as *const v128);
                        let b3 = v128_load(b_col_ptr.add(p + 12) as *const v128);

                        acc0 = f32x4_add(acc0, f32x4_mul(a0, b0));
                        acc1 = f32x4_add(acc1, f32x4_mul(a1, b1));
                        acc2 = f32x4_add(acc2, f32x4_mul(a2, b2));
                        acc3 = f32x4_add(acc3, f32x4_mul(a3, b3));

                        p += 16;
                    }

                    // Horizontal-reduce each accumulator using
                    // (l0+l1)+(l2+l3). wasm v128 has no native
                    // horizontal sum instruction — pairwise extract +
                    // scalar add is the standard pattern. Matches the
                    // scalar reference's reduction order exactly.
                    let h0 = (f32x4_extract_lane::<0>(acc0) + f32x4_extract_lane::<1>(acc0))
                        + (f32x4_extract_lane::<2>(acc0) + f32x4_extract_lane::<3>(acc0));
                    let h1 = (f32x4_extract_lane::<0>(acc1) + f32x4_extract_lane::<1>(acc1))
                        + (f32x4_extract_lane::<2>(acc1) + f32x4_extract_lane::<3>(acc1));
                    let h2 = (f32x4_extract_lane::<0>(acc2) + f32x4_extract_lane::<1>(acc2))
                        + (f32x4_extract_lane::<2>(acc2) + f32x4_extract_lane::<3>(acc2));
                    let h3 = (f32x4_extract_lane::<0>(acc3) + f32x4_extract_lane::<1>(acc3))
                        + (f32x4_extract_lane::<2>(acc3) + f32x4_extract_lane::<3>(acc3));
                    let mut acc = (h0 + h1) + (h2 + h3);

                    // Scalar tail (≤15 elements). Same left-to-right
                    // order as the scalar reference.
                    while p < k {
                        // SAFETY: p < k ≤ column length;
                        // row_stride==1.
                        acc += *a_col_ptr.add(p) * *b_col_ptr.add(p);
                        p += 1;
                    }

                    // SAFETY: (i, j) is in-bounds for dst;
                    // row_stride==1 so the row offset is just `i`.
                    *dst_ptr.offset(j as isize * dst_col + i as isize) = acc;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::linalg::matmul::matmul;
    use faer::{Accum, Mat, Par};

    /// Splitmix64 — small, deterministic, no dep. Good enough to
    /// generate test data with predictable seeds.
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
                // Map to [-1, 1).
                let u = (z >> 40) as u32; // 24 bits
                let f = (u as f32) / ((1u32 << 24) as f32);
                m[(i, j)] = f * 2.0 - 1.0;
            }
        }
        m
    }

    /// Compare scalar reference against faer's matmul for a given
    /// (m, n, k). Tolerance is loose enough to absorb the small
    /// summation-order differences between the scalar reference's
    /// 4-lane reduction and faer's blocked accumulator.
    fn assert_matches_faer(m: usize, n: usize, k: usize) {
        let lhs = random_mat(1, k, m);
        let rhs = random_mat(2, k, n);

        let mut got = Mat::<f32>::zeros(m, n);
        gemm_tn_f32_scalar(got.as_mut(), lhs.as_ref(), rhs.as_ref());

        let mut expected = Mat::<f32>::zeros(m, n);
        matmul(
            expected.as_mut(),
            Accum::Replace,
            lhs.as_ref().transpose(),
            rhs.as_ref(),
            1.0f32,
            Par::Seq,
        );

        for j in 0..n {
            for i in 0..m {
                let g = got[(i, j)];
                let e = expected[(i, j)];
                let abs_err = (g - e).abs();
                // Relative error in the same order as the dot-product
                // magnitude: ~sqrt(k) * 1e-7 for IEEE-754 f32.
                let denom = e.abs().max(1e-6);
                let rel_err = abs_err / denom;
                assert!(
                    abs_err < 1e-3 || rel_err < 1e-4,
                    "mismatch at ({i},{j}) for shape ({m},{n},{k}): got {g}, expected {e}, abs_err {abs_err}, rel_err {rel_err}",
                );
            }
        }
    }

    #[test]
    fn scalar_matches_faer_small_square() {
        assert_matches_faer(50, 50, 100);
    }

    #[test]
    fn scalar_matches_faer_realistic_chunk() {
        // (m=200, n=200, k=2490) is the dominant within-chunk B^T·B
        // shape on 1000G N=2,490 with the default --chunk-size 200.
        assert_matches_faer(200, 200, 2490);
    }

    #[test]
    fn scalar_matches_faer_cross_window() {
        // (m=8000, n=200, k=2490) is a typical cross-window A^T·B at
        // --ld-wind-kb 1000 on 1000G EUR (window can carry ~8K SNPs).
        // Capped smaller here to keep CI fast.
        assert_matches_faer(2000, 200, 2490);
    }

    #[test]
    fn scalar_matches_faer_k_tail_15() {
        // K = 16*q + 15: maximum scalar tail. Verifies the tail loop.
        assert_matches_faer(64, 64, 2495);
    }

    #[test]
    fn scalar_matches_faer_k_under_16() {
        // K shorter than the unroll factor: only the tail loop runs.
        assert_matches_faer(32, 32, 7);
    }

    #[test]
    fn scalar_matches_faer_k_zero() {
        // Degenerate K=0: every output cell should be 0.
        let m = 8;
        let n = 8;
        let k = 0;
        let lhs = Mat::<f32>::zeros(k, m);
        let rhs = Mat::<f32>::zeros(k, n);
        let mut got = Mat::<f32>::zeros(m, n);
        gemm_tn_f32_scalar(got.as_mut(), lhs.as_ref(), rhs.as_ref());
        for j in 0..n {
            for i in 0..m {
                assert_eq!(got[(i, j)], 0.0);
            }
        }
    }
}
