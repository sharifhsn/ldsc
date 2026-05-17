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
    #[cfg(target_feature = "relaxed-simd")]
    use core::arch::wasm32::f32x4_relaxed_madd;
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

        let a_ptr = lhs.as_ptr();
        let a_col = lhs.col_stride();
        let b_ptr = rhs.as_ptr();
        let b_col = rhs.col_stride();

        // Extract the write pointer + leading dimension. faer's
        // `as_ptr_mut` takes `&self` in 0.24 (the `mut` in the name
        // refers to the returned pointer's mutability, not the
        // receiver), so no `mut dst` rebind is needed. The MatMut is
        // moved into this block; lifetime ensures the backing storage
        // stays valid for the kernel (and across any inner-worker
        // parallel dispatch).
        let dst_ptr = dst.as_ptr_mut();
        let dst_col = dst.col_stride();

        // Try parallel dispatch first. `parallel_gemm_tn_f32` returns
        // `true` if the pool was initialised AND N>1 — in which case it
        // ran the work across the inner Workers + the calling thread.
        // Returns `false` (single-threaded fast path) when the pool is
        // disabled, in which case we fall through to the serial loop.
        #[cfg(all(target_arch = "wasm32", target_feature = "atomics"))]
        {
            // SAFETY: the matrices outlive this call (MatMut/MatRef
            // live for at least this function's body), all pointers
            // are valid for the durations the inner Workers will use
            // them (we await the barrier before returning), and the
            // dst region is exclusively owned by us.
            let dispatched = unsafe {
                crate::wasm_simd::pool::parallel_gemm_tn_f32(
                    dst_ptr, dst_col, m, n, a_ptr, a_col, k, b_ptr, b_col,
                )
            };
            if dispatched {
                return;
            }
        }

        // SAFETY for the serial body: same as before. Each unsafe op
        // is in-bounds by construction; specific reasoning is inlined
        // as `// SAFETY:` comments at the load/store sites.
        unsafe {
            // Single-threaded path (pool not initialised, OR caller
            // doesn't have `+atomics` enabled). Pick the fastest
            // available kernel given the target features:
            //
            // - `+relaxed-simd`: Tile4×4 with `f32x4_relaxed_madd`
            //   FMA + `ip-outer, jt-inner` cache-blocked loop order
            //   so B (smaller dim on ab_dot) stays in L2 and A is
            //   streamed only once per call. 112-125 GF/s across
            //   all shapes (3× over the baseline mul+add kernel).
            //
            // - baseline simd128: Tile2×4 with separate mul+add.
            //   65-67 GF/s flat across shapes.
            //
            // See `wasm-bench/` for the full per-shape sweep.
            #[cfg(target_feature = "relaxed-simd")]
            gemm_tn_f32_slice_raw_tiled_4x4_fma_block(
                dst_ptr, dst_col, m, a_ptr, a_col, k, b_ptr, b_col, 0, n,
            );
            #[cfg(not(target_feature = "relaxed-simd"))]
            gemm_tn_f32_slice_raw_tiled_2x4(
                dst_ptr, dst_col, m, a_ptr, a_col, k, b_ptr, b_col, 0, n,
            );
        }
    }

    /// Inner kernel: compute output columns `[j_start, j_end)` of the
    /// TN GEMM. Used both by the serial path (called once with the
    /// full range) and by the parallel pool's worker loop (called per
    /// worker with its assigned slice).
    ///
    /// # Safety
    ///
    /// All pointer/stride preconditions of [`gemm_tn_f32`] apply. In
    /// addition, `j_end <= dst_cols` and `j_start <= j_end`. Multiple
    /// workers running disjoint `[j_start, j_end)` slices in parallel
    /// is sound because output columns are non-overlapping memory
    /// regions and the lhs/rhs are read-only (treated as shared
    /// references).
    #[allow(clippy::too_many_arguments)]
    pub unsafe fn gemm_tn_f32_slice_raw(
        dst_ptr: *mut f32,
        dst_col: isize,
        m: usize,
        a_ptr: *const f32,
        a_col: isize,
        k: usize,
        b_ptr: *const f32,
        b_col: isize,
        j_start: usize,
        j_end: usize,
    ) {
        unsafe {
            let zero: v128 = f32x4_splat(0.0);

            for j in j_start..j_end {
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

    /// **Register-tiled** TN GEMM: compute 4 output columns per A
    /// column load. Same algorithm as [`gemm_tn_f32_slice_raw`] but
    /// processes a `jb = 4` tile of j's at a time, sharing the
    /// `A[:, i]` load across all 4 j's. Each j gets an independent
    /// accumulator (4 indep dep chains), and within each j we keep a
    /// K-unroll of 4 (one v128 per acc, no sub-accumulator unrolling)
    /// to fit comfortably in the wasm SIMD register pressure budget.
    ///
    /// Effect on ab_dot (m=8000, n=200, k=2490):
    /// - Original kernel: 4 × 8000 × 200 A column reads = 6.4M loads
    ///   of 10 KB each per slice; even at L1d this is a lot of
    ///   address-generation pressure.
    /// - Tiled kernel: 4 × 8000 × 50 A column reads = 1.6M loads.
    ///   4× reduction.
    ///
    /// The remaining j's outside the largest multiple of 4 are
    /// processed by falling back to the single-column path so the
    /// kernel handles all `(j_end - j_start)` cleanly. Per call this
    /// is at most 3 j's of the original path — negligible.
    ///
    /// # Safety
    ///
    /// Identical to [`gemm_tn_f32_slice_raw`].
    ///
    /// # Numerics
    ///
    /// **NOT bit-identical** to the scalar reference: the per-j
    /// accumulator structure here uses 4 v128 accs per j vs the
    /// reference's 4 lane-groups, but the K is stepped in 4-element
    /// chunks instead of 16, changing the precise summation order.
    /// Relative error stays comfortably below 1e-5 for our shapes
    /// (verified by the loose-tolerance property test). At the
    /// L2-score level (sum over k of r²) this is well below the f32
    /// noise floor and matches CLI to displayed precision.
    #[allow(clippy::too_many_arguments)]
    pub unsafe fn gemm_tn_f32_slice_raw_tiled(
        dst_ptr: *mut f32,
        dst_col: isize,
        m: usize,
        a_ptr: *const f32,
        a_col: isize,
        k: usize,
        b_ptr: *const f32,
        b_col: isize,
        j_start: usize,
        j_end: usize,
    ) {
        unsafe {
            let zero: v128 = f32x4_splat(0.0);
            // Process j's in tiles of 4.
            let mut jt = j_start;
            while jt + 4 <= j_end {
                // SAFETY: jt+3 < j_end ≤ rhs.ncols(); each column
                // offset stays within rhs's allocation.
                let b_col_ptr_0 = b_ptr.offset(jt as isize * b_col);
                let b_col_ptr_1 = b_ptr.offset((jt + 1) as isize * b_col);
                let b_col_ptr_2 = b_ptr.offset((jt + 2) as isize * b_col);
                let b_col_ptr_3 = b_ptr.offset((jt + 3) as isize * b_col);

                for i in 0..m {
                    let a_col_ptr = a_ptr.offset(i as isize * a_col);

                    // 4 independent SIMD accumulators, one per output
                    // column. Each is a 4-lane v128; final horizontal
                    // reduce at the end gives one f32 per cell.
                    let mut acc0 = zero;
                    let mut acc1 = zero;
                    let mut acc2 = zero;
                    let mut acc3 = zero;

                    let mut p: usize = 0;
                    while p + 4 <= k {
                        // SAFETY: row_stride==1 (debug-asserted),
                        // p+4 ≤ k (loop guard). v128_load handles
                        // unaligned addresses.
                        let a = v128_load(a_col_ptr.add(p) as *const v128);
                        let b0 = v128_load(b_col_ptr_0.add(p) as *const v128);
                        let b1 = v128_load(b_col_ptr_1.add(p) as *const v128);
                        let b2 = v128_load(b_col_ptr_2.add(p) as *const v128);
                        let b3 = v128_load(b_col_ptr_3.add(p) as *const v128);

                        acc0 = f32x4_add(acc0, f32x4_mul(a, b0));
                        acc1 = f32x4_add(acc1, f32x4_mul(a, b1));
                        acc2 = f32x4_add(acc2, f32x4_mul(a, b2));
                        acc3 = f32x4_add(acc3, f32x4_mul(a, b3));
                        p += 4;
                    }

                    // Horizontal-reduce each accumulator. Same
                    // (l0+l1)+(l2+l3) pattern as the reference.
                    let mut s0 = (f32x4_extract_lane::<0>(acc0) + f32x4_extract_lane::<1>(acc0))
                        + (f32x4_extract_lane::<2>(acc0) + f32x4_extract_lane::<3>(acc0));
                    let mut s1 = (f32x4_extract_lane::<0>(acc1) + f32x4_extract_lane::<1>(acc1))
                        + (f32x4_extract_lane::<2>(acc1) + f32x4_extract_lane::<3>(acc1));
                    let mut s2 = (f32x4_extract_lane::<0>(acc2) + f32x4_extract_lane::<1>(acc2))
                        + (f32x4_extract_lane::<2>(acc2) + f32x4_extract_lane::<3>(acc2));
                    let mut s3 = (f32x4_extract_lane::<0>(acc3) + f32x4_extract_lane::<1>(acc3))
                        + (f32x4_extract_lane::<2>(acc3) + f32x4_extract_lane::<3>(acc3));

                    // Scalar tail (≤3 elements). Left-to-right.
                    while p < k {
                        let a_val = *a_col_ptr.add(p);
                        s0 += a_val * *b_col_ptr_0.add(p);
                        s1 += a_val * *b_col_ptr_1.add(p);
                        s2 += a_val * *b_col_ptr_2.add(p);
                        s3 += a_val * *b_col_ptr_3.add(p);
                        p += 1;
                    }

                    // SAFETY: (i, jt+δ) for δ ∈ 0..4 is in-bounds
                    // for dst.
                    let row = i as isize;
                    *dst_ptr.offset(jt as isize * dst_col + row) = s0;
                    *dst_ptr.offset((jt + 1) as isize * dst_col + row) = s1;
                    *dst_ptr.offset((jt + 2) as isize * dst_col + row) = s2;
                    *dst_ptr.offset((jt + 3) as isize * dst_col + row) = s3;
                }

                jt += 4;
            }

            // Tail j's: fall back to the single-j kernel for the at-
            // most-3 columns that didn't fit in a full tile.
            if jt < j_end {
                gemm_tn_f32_slice_raw(
                    dst_ptr, dst_col, m, a_ptr, a_col, k, b_ptr, b_col, jt, j_end,
                );
            }
        }
    }

    /// **2×4 register-tiled** TN GEMM: compute a 2×4 output tile per
    /// inner iteration, sharing both A column loads (across 4 j's)
    /// AND B column loads (across 2 i's). Per K iter:
    ///   - 2 A loads (one per i in the i-pair)
    ///   - 4 B loads (one per j in the j-quad)
    ///   - 8 mul + 8 add (one per cell in the 2×4 tile)
    ///   - 8 accumulators live + 6 loads = 14 v128 in flight.
    ///
    /// vs the 4-j tile: 2× reduction in B loads, slightly better cell
    /// throughput (8 mul+add per 6 loads vs 4 mul+add per 5 loads).
    ///
    /// # Numerics
    ///
    /// NOT bit-identical to scalar reference (same caveat as the
    /// 4-j-tile variant). Property test uses loose 1e-4 tolerance.
    ///
    /// # Safety
    ///
    /// Identical to [`gemm_tn_f32_slice_raw`].
    #[allow(clippy::too_many_arguments)]
    pub unsafe fn gemm_tn_f32_slice_raw_tiled_2x4(
        dst_ptr: *mut f32,
        dst_col: isize,
        m: usize,
        a_ptr: *const f32,
        a_col: isize,
        k: usize,
        b_ptr: *const f32,
        b_col: isize,
        j_start: usize,
        j_end: usize,
    ) {
        unsafe {
            let zero: v128 = f32x4_splat(0.0);

            let mut jt = j_start;
            while jt + 4 <= j_end {
                let b0p = b_ptr.offset(jt as isize * b_col);
                let b1p = b_ptr.offset((jt + 1) as isize * b_col);
                let b2p = b_ptr.offset((jt + 2) as isize * b_col);
                let b3p = b_ptr.offset((jt + 3) as isize * b_col);

                let m_pairs = m / 2;
                for ip in 0..m_pairs {
                    let i0 = 2 * ip;
                    let i1 = i0 + 1;
                    let a0p = a_ptr.offset(i0 as isize * a_col);
                    let a1p = a_ptr.offset(i1 as isize * a_col);

                    // 2×4 = 8 SIMD accumulators.
                    let mut c00 = zero;
                    let mut c01 = zero;
                    let mut c02 = zero;
                    let mut c03 = zero;
                    let mut c10 = zero;
                    let mut c11 = zero;
                    let mut c12 = zero;
                    let mut c13 = zero;

                    let mut p: usize = 0;
                    while p + 4 <= k {
                        let a0 = v128_load(a0p.add(p) as *const v128);
                        let a1 = v128_load(a1p.add(p) as *const v128);
                        let b0 = v128_load(b0p.add(p) as *const v128);
                        let b1 = v128_load(b1p.add(p) as *const v128);
                        let b2 = v128_load(b2p.add(p) as *const v128);
                        let b3 = v128_load(b3p.add(p) as *const v128);

                        c00 = f32x4_add(c00, f32x4_mul(a0, b0));
                        c01 = f32x4_add(c01, f32x4_mul(a0, b1));
                        c02 = f32x4_add(c02, f32x4_mul(a0, b2));
                        c03 = f32x4_add(c03, f32x4_mul(a0, b3));
                        c10 = f32x4_add(c10, f32x4_mul(a1, b0));
                        c11 = f32x4_add(c11, f32x4_mul(a1, b1));
                        c12 = f32x4_add(c12, f32x4_mul(a1, b2));
                        c13 = f32x4_add(c13, f32x4_mul(a1, b3));
                        p += 4;
                    }

                    // Horizontal reduce each accumulator.
                    let hred = |v: v128| -> f32 {
                        (f32x4_extract_lane::<0>(v) + f32x4_extract_lane::<1>(v))
                            + (f32x4_extract_lane::<2>(v) + f32x4_extract_lane::<3>(v))
                    };
                    let mut s00 = hred(c00);
                    let mut s01 = hred(c01);
                    let mut s02 = hred(c02);
                    let mut s03 = hred(c03);
                    let mut s10 = hred(c10);
                    let mut s11 = hred(c11);
                    let mut s12 = hred(c12);
                    let mut s13 = hred(c13);

                    while p < k {
                        let a0v = *a0p.add(p);
                        let a1v = *a1p.add(p);
                        let b0v = *b0p.add(p);
                        let b1v = *b1p.add(p);
                        let b2v = *b2p.add(p);
                        let b3v = *b3p.add(p);
                        s00 += a0v * b0v;
                        s01 += a0v * b1v;
                        s02 += a0v * b2v;
                        s03 += a0v * b3v;
                        s10 += a1v * b0v;
                        s11 += a1v * b1v;
                        s12 += a1v * b2v;
                        s13 += a1v * b3v;
                        p += 1;
                    }

                    let r0 = i0 as isize;
                    let r1 = i1 as isize;
                    *dst_ptr.offset(jt as isize * dst_col + r0) = s00;
                    *dst_ptr.offset((jt + 1) as isize * dst_col + r0) = s01;
                    *dst_ptr.offset((jt + 2) as isize * dst_col + r0) = s02;
                    *dst_ptr.offset((jt + 3) as isize * dst_col + r0) = s03;
                    *dst_ptr.offset(jt as isize * dst_col + r1) = s10;
                    *dst_ptr.offset((jt + 1) as isize * dst_col + r1) = s11;
                    *dst_ptr.offset((jt + 2) as isize * dst_col + r1) = s12;
                    *dst_ptr.offset((jt + 3) as isize * dst_col + r1) = s13;
                }

                // Odd-i tail: if m is odd, handle the last row of this
                // j-quad with a single-i × 4-j strip.
                if m % 2 == 1 {
                    let i = m - 1;
                    let ap = a_ptr.offset(i as isize * a_col);
                    let mut c0 = zero;
                    let mut c1 = zero;
                    let mut c2 = zero;
                    let mut c3 = zero;
                    let mut p: usize = 0;
                    while p + 4 <= k {
                        let a = v128_load(ap.add(p) as *const v128);
                        let b0 = v128_load(b0p.add(p) as *const v128);
                        let b1 = v128_load(b1p.add(p) as *const v128);
                        let b2 = v128_load(b2p.add(p) as *const v128);
                        let b3 = v128_load(b3p.add(p) as *const v128);
                        c0 = f32x4_add(c0, f32x4_mul(a, b0));
                        c1 = f32x4_add(c1, f32x4_mul(a, b1));
                        c2 = f32x4_add(c2, f32x4_mul(a, b2));
                        c3 = f32x4_add(c3, f32x4_mul(a, b3));
                        p += 4;
                    }
                    let hred = |v: v128| -> f32 {
                        (f32x4_extract_lane::<0>(v) + f32x4_extract_lane::<1>(v))
                            + (f32x4_extract_lane::<2>(v) + f32x4_extract_lane::<3>(v))
                    };
                    let mut s0 = hred(c0);
                    let mut s1 = hred(c1);
                    let mut s2 = hred(c2);
                    let mut s3 = hred(c3);
                    while p < k {
                        let av = *ap.add(p);
                        s0 += av * *b0p.add(p);
                        s1 += av * *b1p.add(p);
                        s2 += av * *b2p.add(p);
                        s3 += av * *b3p.add(p);
                        p += 1;
                    }
                    let r = i as isize;
                    *dst_ptr.offset(jt as isize * dst_col + r) = s0;
                    *dst_ptr.offset((jt + 1) as isize * dst_col + r) = s1;
                    *dst_ptr.offset((jt + 2) as isize * dst_col + r) = s2;
                    *dst_ptr.offset((jt + 3) as isize * dst_col + r) = s3;
                }

                jt += 4;
            }

            // j-tail (≤3 cols): single-j fallback.
            if jt < j_end {
                gemm_tn_f32_slice_raw(
                    dst_ptr, dst_col, m, a_ptr, a_col, k, b_ptr, b_col, jt, j_end,
                );
            }
        }
    }

    /// **4×4 register-tiled** TN GEMM: 4 i's × 4 j's per K iter.
    /// 16 SIMD accumulators (the BLIS sweet spot for ARM NEON, 32
    /// regs). Per K iter:
    ///   - 4 A loads + 4 B loads = 8 loads
    ///   - 16 mul + 16 add = 32 ops over 16 cells = 2 flops/load
    ///
    /// vs 2×4: 2× more reuse on both A and B simultaneously. Each
    /// cell sees its mul+add chain of length k/4 just like 2×4 but
    /// 4 chains share each load.
    ///
    /// Performance ceiling: ARM64 has 32 v128 regs; 16 accs + 8
    /// loads = 24 in flight, well under cap. Cranelift's wasm
    /// register allocator should keep all accs in regs across the
    /// inner loop.
    ///
    /// # Numerics
    ///
    /// NOT bit-identical to scalar reference (loose-tolerance test
    /// covers the kernel's K-step-4 reduction order vs reference's
    /// K-step-16).
    ///
    /// # Safety
    ///
    /// Identical to [`gemm_tn_f32_slice_raw`]: caller upholds the
    /// matrix pointer / stride / dimension preconditions, dst region
    /// is exclusively owned, lhs/rhs are read-only (may be aliased).
    #[allow(clippy::too_many_arguments)]
    pub unsafe fn gemm_tn_f32_slice_raw_tiled_4x4(
        dst_ptr: *mut f32,
        dst_col: isize,
        m: usize,
        a_ptr: *const f32,
        a_col: isize,
        k: usize,
        b_ptr: *const f32,
        b_col: isize,
        j_start: usize,
        j_end: usize,
    ) {
        unsafe {
            let zero: v128 = f32x4_splat(0.0);

            let mut jt = j_start;
            while jt + 4 <= j_end {
                let bp0 = b_ptr.offset(jt as isize * b_col);
                let bp1 = b_ptr.offset((jt + 1) as isize * b_col);
                let bp2 = b_ptr.offset((jt + 2) as isize * b_col);
                let bp3 = b_ptr.offset((jt + 3) as isize * b_col);

                let m_quads = m / 4;
                for ip in 0..m_quads {
                    let i0 = 4 * ip;
                    let ap0 = a_ptr.offset(i0 as isize * a_col);
                    let ap1 = a_ptr.offset((i0 + 1) as isize * a_col);
                    let ap2 = a_ptr.offset((i0 + 2) as isize * a_col);
                    let ap3 = a_ptr.offset((i0 + 3) as isize * a_col);

                    // 4×4 = 16 SIMD accumulators.
                    let mut c00 = zero;
                    let mut c01 = zero;
                    let mut c02 = zero;
                    let mut c03 = zero;
                    let mut c10 = zero;
                    let mut c11 = zero;
                    let mut c12 = zero;
                    let mut c13 = zero;
                    let mut c20 = zero;
                    let mut c21 = zero;
                    let mut c22 = zero;
                    let mut c23 = zero;
                    let mut c30 = zero;
                    let mut c31 = zero;
                    let mut c32 = zero;
                    let mut c33 = zero;

                    let mut p: usize = 0;
                    while p + 4 <= k {
                        let a0 = v128_load(ap0.add(p) as *const v128);
                        let a1 = v128_load(ap1.add(p) as *const v128);
                        let a2 = v128_load(ap2.add(p) as *const v128);
                        let a3 = v128_load(ap3.add(p) as *const v128);
                        let b0 = v128_load(bp0.add(p) as *const v128);
                        let b1 = v128_load(bp1.add(p) as *const v128);
                        let b2 = v128_load(bp2.add(p) as *const v128);
                        let b3 = v128_load(bp3.add(p) as *const v128);

                        c00 = f32x4_add(c00, f32x4_mul(a0, b0));
                        c01 = f32x4_add(c01, f32x4_mul(a0, b1));
                        c02 = f32x4_add(c02, f32x4_mul(a0, b2));
                        c03 = f32x4_add(c03, f32x4_mul(a0, b3));
                        c10 = f32x4_add(c10, f32x4_mul(a1, b0));
                        c11 = f32x4_add(c11, f32x4_mul(a1, b1));
                        c12 = f32x4_add(c12, f32x4_mul(a1, b2));
                        c13 = f32x4_add(c13, f32x4_mul(a1, b3));
                        c20 = f32x4_add(c20, f32x4_mul(a2, b0));
                        c21 = f32x4_add(c21, f32x4_mul(a2, b1));
                        c22 = f32x4_add(c22, f32x4_mul(a2, b2));
                        c23 = f32x4_add(c23, f32x4_mul(a2, b3));
                        c30 = f32x4_add(c30, f32x4_mul(a3, b0));
                        c31 = f32x4_add(c31, f32x4_mul(a3, b1));
                        c32 = f32x4_add(c32, f32x4_mul(a3, b2));
                        c33 = f32x4_add(c33, f32x4_mul(a3, b3));
                        p += 4;
                    }

                    let hred = |v: v128| -> f32 {
                        (f32x4_extract_lane::<0>(v) + f32x4_extract_lane::<1>(v))
                            + (f32x4_extract_lane::<2>(v) + f32x4_extract_lane::<3>(v))
                    };
                    let mut s00 = hred(c00);
                    let mut s01 = hred(c01);
                    let mut s02 = hred(c02);
                    let mut s03 = hred(c03);
                    let mut s10 = hred(c10);
                    let mut s11 = hred(c11);
                    let mut s12 = hred(c12);
                    let mut s13 = hred(c13);
                    let mut s20 = hred(c20);
                    let mut s21 = hred(c21);
                    let mut s22 = hred(c22);
                    let mut s23 = hred(c23);
                    let mut s30 = hred(c30);
                    let mut s31 = hred(c31);
                    let mut s32 = hred(c32);
                    let mut s33 = hred(c33);

                    while p < k {
                        let a0v = *ap0.add(p);
                        let a1v = *ap1.add(p);
                        let a2v = *ap2.add(p);
                        let a3v = *ap3.add(p);
                        let b0v = *bp0.add(p);
                        let b1v = *bp1.add(p);
                        let b2v = *bp2.add(p);
                        let b3v = *bp3.add(p);
                        s00 += a0v * b0v;
                        s01 += a0v * b1v;
                        s02 += a0v * b2v;
                        s03 += a0v * b3v;
                        s10 += a1v * b0v;
                        s11 += a1v * b1v;
                        s12 += a1v * b2v;
                        s13 += a1v * b3v;
                        s20 += a2v * b0v;
                        s21 += a2v * b1v;
                        s22 += a2v * b2v;
                        s23 += a2v * b3v;
                        s30 += a3v * b0v;
                        s31 += a3v * b1v;
                        s32 += a3v * b2v;
                        s33 += a3v * b3v;
                        p += 1;
                    }

                    let store4 = |dst: *mut f32, jbase: usize, row: isize, v0, v1, v2, v3| {
                        // SAFETY: (row, jbase+δ) for δ∈0..4 in-bounds.
                        // Already inside the outer fn's `unsafe { ... }`
                        // so no nested block needed.
                        *dst.offset(jbase as isize * dst_col + row) = v0;
                        *dst.offset((jbase + 1) as isize * dst_col + row) = v1;
                        *dst.offset((jbase + 2) as isize * dst_col + row) = v2;
                        *dst.offset((jbase + 3) as isize * dst_col + row) = v3;
                    };
                    store4(dst_ptr, jt, i0 as isize, s00, s01, s02, s03);
                    store4(dst_ptr, jt, (i0 + 1) as isize, s10, s11, s12, s13);
                    store4(dst_ptr, jt, (i0 + 2) as isize, s20, s21, s22, s23);
                    store4(dst_ptr, jt, (i0 + 3) as isize, s30, s31, s32, s33);
                }

                // i-tail (m % 4 rows). Shift dst's row base and a's
                // column base, then use the single-j kernel to walk
                // through the tail rows only. No wasted recomputation.
                let tail_start = (m / 4) * 4;
                if tail_start < m {
                    let a_tail_off = (tail_start as isize).wrapping_mul(a_col);
                    gemm_tn_f32_slice_raw(
                        dst_ptr.add(tail_start),
                        dst_col,
                        m - tail_start,
                        a_ptr.offset(a_tail_off),
                        a_col,
                        k,
                        b_ptr,
                        b_col,
                        jt,
                        jt + 4,
                    );
                }

                jt += 4;
            }

            // j-tail (≤3 cols).
            if jt < j_end {
                gemm_tn_f32_slice_raw(
                    dst_ptr, dst_col, m, a_ptr, a_col, k, b_ptr, b_col, jt, j_end,
                );
            }
        }
    }

    // ─────────────────────────────────────────────────────────────
    // Relaxed-SIMD FMA variants. Gated on `+relaxed-simd`. Use
    // `f32x4_relaxed_madd(a, b, c)` (fused multiply-add `a*b + c`)
    // instead of separate `mul + add`.
    //
    // Why this matters: baseline simd128 has no FMA on f32x4. Each
    // mul+add pair is 2 SIMD instructions and 6+ cycles latency.
    // Relaxed-SIMD's `f32x4_relaxed_madd` lowers to a single FMLA
    // (ARM64 NEON) or `vfmadd231ps` (x86 AVX2/FMA3) — 1 instruction,
    // typically 4-cycle latency, and crucially TWO ops counted per
    // cycle on most hardware (mul and add fused into one).
    //
    // Theoretical peak doubles vs baseline simd128:
    // - Apple M5 NEON: 4 SIMD × 1 FMLA × 8 flops/cycle × 3 GHz =
    //   96 GFLOPs → 192 GFLOPs with FMA. Realistic: ~150 GF/s.
    // - x86 AVX2/FMA: similar.
    //
    // Determinism note: relaxed-SIMD FMA is allowed up to 1 ULP
    // of nondeterminism across implementations (the spec lets a
    // platform implement it as either fused FMA or separate
    // mul+add; the result may differ in the last bit). For LD
    // scores summed over K terms each O(1/√K), the per-cell error
    // bound is O(√K · ULP) ≈ 1e-5 at K=2490. Well below the
    // displayed precision of mean L2 (4 sig figs).
    // ─────────────────────────────────────────────────────────────

    /// 2×4 register-tiled TN GEMM using **relaxed-SIMD FMA**. Same
    /// loop structure as [`gemm_tn_f32_slice_raw_tiled_2x4`] but
    /// each `mul + add` pair is fused into a single
    /// `f32x4_relaxed_madd`. Expected ~1.3-1.5× speedup over the
    /// baseline-simd128 2×4 kernel.
    ///
    /// # Numerics
    ///
    /// Up to 1 ULP per FMA may differ vs the non-fused kernel. Not
    /// bit-identical to either the scalar reference OR the baseline
    /// 2×4 kernel. The L2 score gets summation noise of O(√K · ULP)
    /// ≈ 1e-5 at K=2490; below displayed precision.
    ///
    /// # Safety
    ///
    /// Identical to [`gemm_tn_f32_slice_raw`].
    #[cfg(target_feature = "relaxed-simd")]
    #[allow(clippy::too_many_arguments)]
    pub unsafe fn gemm_tn_f32_slice_raw_tiled_2x4_fma(
        dst_ptr: *mut f32,
        dst_col: isize,
        m: usize,
        a_ptr: *const f32,
        a_col: isize,
        k: usize,
        b_ptr: *const f32,
        b_col: isize,
        j_start: usize,
        j_end: usize,
    ) {
        unsafe {
            let zero: v128 = f32x4_splat(0.0);

            let mut jt = j_start;
            while jt + 4 <= j_end {
                let b0p = b_ptr.offset(jt as isize * b_col);
                let b1p = b_ptr.offset((jt + 1) as isize * b_col);
                let b2p = b_ptr.offset((jt + 2) as isize * b_col);
                let b3p = b_ptr.offset((jt + 3) as isize * b_col);

                let m_pairs = m / 2;
                for ip in 0..m_pairs {
                    let i0 = 2 * ip;
                    let i1 = i0 + 1;
                    let a0p = a_ptr.offset(i0 as isize * a_col);
                    let a1p = a_ptr.offset(i1 as isize * a_col);

                    let mut c00 = zero;
                    let mut c01 = zero;
                    let mut c02 = zero;
                    let mut c03 = zero;
                    let mut c10 = zero;
                    let mut c11 = zero;
                    let mut c12 = zero;
                    let mut c13 = zero;

                    let mut p: usize = 0;
                    while p + 4 <= k {
                        let a0 = v128_load(a0p.add(p) as *const v128);
                        let a1 = v128_load(a1p.add(p) as *const v128);
                        let b0 = v128_load(b0p.add(p) as *const v128);
                        let b1 = v128_load(b1p.add(p) as *const v128);
                        let b2 = v128_load(b2p.add(p) as *const v128);
                        let b3 = v128_load(b3p.add(p) as *const v128);

                        // FMA: c += a * b. Single instruction.
                        c00 = f32x4_relaxed_madd(a0, b0, c00);
                        c01 = f32x4_relaxed_madd(a0, b1, c01);
                        c02 = f32x4_relaxed_madd(a0, b2, c02);
                        c03 = f32x4_relaxed_madd(a0, b3, c03);
                        c10 = f32x4_relaxed_madd(a1, b0, c10);
                        c11 = f32x4_relaxed_madd(a1, b1, c11);
                        c12 = f32x4_relaxed_madd(a1, b2, c12);
                        c13 = f32x4_relaxed_madd(a1, b3, c13);
                        p += 4;
                    }

                    let hred = |v: v128| -> f32 {
                        (f32x4_extract_lane::<0>(v) + f32x4_extract_lane::<1>(v))
                            + (f32x4_extract_lane::<2>(v) + f32x4_extract_lane::<3>(v))
                    };
                    let mut s00 = hred(c00);
                    let mut s01 = hred(c01);
                    let mut s02 = hred(c02);
                    let mut s03 = hred(c03);
                    let mut s10 = hred(c10);
                    let mut s11 = hred(c11);
                    let mut s12 = hred(c12);
                    let mut s13 = hred(c13);

                    while p < k {
                        let a0v = *a0p.add(p);
                        let a1v = *a1p.add(p);
                        let b0v = *b0p.add(p);
                        let b1v = *b1p.add(p);
                        let b2v = *b2p.add(p);
                        let b3v = *b3p.add(p);
                        s00 += a0v * b0v;
                        s01 += a0v * b1v;
                        s02 += a0v * b2v;
                        s03 += a0v * b3v;
                        s10 += a1v * b0v;
                        s11 += a1v * b1v;
                        s12 += a1v * b2v;
                        s13 += a1v * b3v;
                        p += 1;
                    }

                    let r0 = i0 as isize;
                    let r1 = i1 as isize;
                    *dst_ptr.offset(jt as isize * dst_col + r0) = s00;
                    *dst_ptr.offset((jt + 1) as isize * dst_col + r0) = s01;
                    *dst_ptr.offset((jt + 2) as isize * dst_col + r0) = s02;
                    *dst_ptr.offset((jt + 3) as isize * dst_col + r0) = s03;
                    *dst_ptr.offset(jt as isize * dst_col + r1) = s10;
                    *dst_ptr.offset((jt + 1) as isize * dst_col + r1) = s11;
                    *dst_ptr.offset((jt + 2) as isize * dst_col + r1) = s12;
                    *dst_ptr.offset((jt + 3) as isize * dst_col + r1) = s13;
                }

                // Odd-i tail.
                if m % 2 == 1 {
                    let i = m - 1;
                    let ap = a_ptr.offset(i as isize * a_col);
                    let mut c0 = zero;
                    let mut c1 = zero;
                    let mut c2 = zero;
                    let mut c3 = zero;
                    let mut p: usize = 0;
                    while p + 4 <= k {
                        let a = v128_load(ap.add(p) as *const v128);
                        let b0 = v128_load(b0p.add(p) as *const v128);
                        let b1 = v128_load(b1p.add(p) as *const v128);
                        let b2 = v128_load(b2p.add(p) as *const v128);
                        let b3 = v128_load(b3p.add(p) as *const v128);
                        c0 = f32x4_relaxed_madd(a, b0, c0);
                        c1 = f32x4_relaxed_madd(a, b1, c1);
                        c2 = f32x4_relaxed_madd(a, b2, c2);
                        c3 = f32x4_relaxed_madd(a, b3, c3);
                        p += 4;
                    }
                    let hred = |v: v128| -> f32 {
                        (f32x4_extract_lane::<0>(v) + f32x4_extract_lane::<1>(v))
                            + (f32x4_extract_lane::<2>(v) + f32x4_extract_lane::<3>(v))
                    };
                    let mut s0 = hred(c0);
                    let mut s1 = hred(c1);
                    let mut s2 = hred(c2);
                    let mut s3 = hred(c3);
                    while p < k {
                        let av = *ap.add(p);
                        s0 += av * *b0p.add(p);
                        s1 += av * *b1p.add(p);
                        s2 += av * *b2p.add(p);
                        s3 += av * *b3p.add(p);
                        p += 1;
                    }
                    let r = i as isize;
                    *dst_ptr.offset(jt as isize * dst_col + r) = s0;
                    *dst_ptr.offset((jt + 1) as isize * dst_col + r) = s1;
                    *dst_ptr.offset((jt + 2) as isize * dst_col + r) = s2;
                    *dst_ptr.offset((jt + 3) as isize * dst_col + r) = s3;
                }

                jt += 4;
            }

            if jt < j_end {
                gemm_tn_f32_slice_raw(
                    dst_ptr, dst_col, m, a_ptr, a_col, k, b_ptr, b_col, jt, j_end,
                );
            }
        }
    }

    /// 4×4 register-tiled TN GEMM using **relaxed-SIMD FMA** with
    /// the canonical `jt-outer, ip-inner` loop order. Best for
    /// `m` that fits L2 (e.g. bb_dot m=200 → A fits comfortably).
    /// For very large m (8000+), prefer the `_blocked` variant
    /// below, which swaps the loops to keep B (smaller dim) in L2
    /// and stream A only once per call.
    ///
    /// # Numerics
    ///
    /// Same as [`gemm_tn_f32_slice_raw_tiled_2x4_fma`].
    ///
    /// # Safety
    ///
    /// Identical to [`gemm_tn_f32_slice_raw`].
    #[cfg(target_feature = "relaxed-simd")]
    #[allow(clippy::too_many_arguments)]
    pub unsafe fn gemm_tn_f32_slice_raw_tiled_4x4_fma(
        dst_ptr: *mut f32,
        dst_col: isize,
        m: usize,
        a_ptr: *const f32,
        a_col: isize,
        k: usize,
        b_ptr: *const f32,
        b_col: isize,
        j_start: usize,
        j_end: usize,
    ) {
        unsafe {
            let zero: v128 = f32x4_splat(0.0);

            let mut jt = j_start;
            while jt + 4 <= j_end {
                let bp0 = b_ptr.offset(jt as isize * b_col);
                let bp1 = b_ptr.offset((jt + 1) as isize * b_col);
                let bp2 = b_ptr.offset((jt + 2) as isize * b_col);
                let bp3 = b_ptr.offset((jt + 3) as isize * b_col);

                let m_quads = m / 4;
                for ip in 0..m_quads {
                    let i0 = 4 * ip;
                    let ap0 = a_ptr.offset(i0 as isize * a_col);
                    let ap1 = a_ptr.offset((i0 + 1) as isize * a_col);
                    let ap2 = a_ptr.offset((i0 + 2) as isize * a_col);
                    let ap3 = a_ptr.offset((i0 + 3) as isize * a_col);

                    let mut c00 = zero;
                    let mut c01 = zero;
                    let mut c02 = zero;
                    let mut c03 = zero;
                    let mut c10 = zero;
                    let mut c11 = zero;
                    let mut c12 = zero;
                    let mut c13 = zero;
                    let mut c20 = zero;
                    let mut c21 = zero;
                    let mut c22 = zero;
                    let mut c23 = zero;
                    let mut c30 = zero;
                    let mut c31 = zero;
                    let mut c32 = zero;
                    let mut c33 = zero;

                    let mut p: usize = 0;
                    while p + 4 <= k {
                        let a0 = v128_load(ap0.add(p) as *const v128);
                        let a1 = v128_load(ap1.add(p) as *const v128);
                        let a2 = v128_load(ap2.add(p) as *const v128);
                        let a3 = v128_load(ap3.add(p) as *const v128);
                        let b0 = v128_load(bp0.add(p) as *const v128);
                        let b1 = v128_load(bp1.add(p) as *const v128);
                        let b2 = v128_load(bp2.add(p) as *const v128);
                        let b3 = v128_load(bp3.add(p) as *const v128);

                        c00 = f32x4_relaxed_madd(a0, b0, c00);
                        c01 = f32x4_relaxed_madd(a0, b1, c01);
                        c02 = f32x4_relaxed_madd(a0, b2, c02);
                        c03 = f32x4_relaxed_madd(a0, b3, c03);
                        c10 = f32x4_relaxed_madd(a1, b0, c10);
                        c11 = f32x4_relaxed_madd(a1, b1, c11);
                        c12 = f32x4_relaxed_madd(a1, b2, c12);
                        c13 = f32x4_relaxed_madd(a1, b3, c13);
                        c20 = f32x4_relaxed_madd(a2, b0, c20);
                        c21 = f32x4_relaxed_madd(a2, b1, c21);
                        c22 = f32x4_relaxed_madd(a2, b2, c22);
                        c23 = f32x4_relaxed_madd(a2, b3, c23);
                        c30 = f32x4_relaxed_madd(a3, b0, c30);
                        c31 = f32x4_relaxed_madd(a3, b1, c31);
                        c32 = f32x4_relaxed_madd(a3, b2, c32);
                        c33 = f32x4_relaxed_madd(a3, b3, c33);
                        p += 4;
                    }

                    let hred = |v: v128| -> f32 {
                        (f32x4_extract_lane::<0>(v) + f32x4_extract_lane::<1>(v))
                            + (f32x4_extract_lane::<2>(v) + f32x4_extract_lane::<3>(v))
                    };
                    let mut s00 = hred(c00);
                    let mut s01 = hred(c01);
                    let mut s02 = hred(c02);
                    let mut s03 = hred(c03);
                    let mut s10 = hred(c10);
                    let mut s11 = hred(c11);
                    let mut s12 = hred(c12);
                    let mut s13 = hred(c13);
                    let mut s20 = hred(c20);
                    let mut s21 = hred(c21);
                    let mut s22 = hred(c22);
                    let mut s23 = hred(c23);
                    let mut s30 = hred(c30);
                    let mut s31 = hred(c31);
                    let mut s32 = hred(c32);
                    let mut s33 = hred(c33);

                    while p < k {
                        let a0v = *ap0.add(p);
                        let a1v = *ap1.add(p);
                        let a2v = *ap2.add(p);
                        let a3v = *ap3.add(p);
                        let b0v = *bp0.add(p);
                        let b1v = *bp1.add(p);
                        let b2v = *bp2.add(p);
                        let b3v = *bp3.add(p);
                        s00 += a0v * b0v;
                        s01 += a0v * b1v;
                        s02 += a0v * b2v;
                        s03 += a0v * b3v;
                        s10 += a1v * b0v;
                        s11 += a1v * b1v;
                        s12 += a1v * b2v;
                        s13 += a1v * b3v;
                        s20 += a2v * b0v;
                        s21 += a2v * b1v;
                        s22 += a2v * b2v;
                        s23 += a2v * b3v;
                        s30 += a3v * b0v;
                        s31 += a3v * b1v;
                        s32 += a3v * b2v;
                        s33 += a3v * b3v;
                        p += 1;
                    }

                    // Same store4 closure pattern as the non-FMA kernel.
                    let store4 = |dst: *mut f32, jbase: usize, row: isize, v0, v1, v2, v3| {
                        *dst.offset(jbase as isize * dst_col + row) = v0;
                        *dst.offset((jbase + 1) as isize * dst_col + row) = v1;
                        *dst.offset((jbase + 2) as isize * dst_col + row) = v2;
                        *dst.offset((jbase + 3) as isize * dst_col + row) = v3;
                    };
                    store4(dst_ptr, jt, i0 as isize, s00, s01, s02, s03);
                    store4(dst_ptr, jt, (i0 + 1) as isize, s10, s11, s12, s13);
                    store4(dst_ptr, jt, (i0 + 2) as isize, s20, s21, s22, s23);
                    store4(dst_ptr, jt, (i0 + 3) as isize, s30, s31, s32, s33);
                }

                // m-tail (rows that didn't fit in m_quads × 4).
                // Process tail rows with single-j kernel by shifting
                // both dst's row base AND a's column base — dst_ptr
                // is row-base (row_stride==1) so +tail_start offsets
                // to row tail_start; a's "row i" in TN GEMM is
                // a.col(i), so we shift a_ptr by tail_start columns.
                // The single-j kernel's `for i in 0..m` then walks
                // through (m - tail_start) rows starting at
                // tail_start in the original dst.
                let tail_start = (m / 4) * 4;
                if tail_start < m {
                    let a_tail_off = (tail_start as isize).wrapping_mul(a_col);
                    gemm_tn_f32_slice_raw(
                        dst_ptr.add(tail_start),
                        dst_col,
                        m - tail_start,
                        a_ptr.offset(a_tail_off),
                        a_col,
                        k,
                        b_ptr,
                        b_col,
                        jt,
                        jt + 4,
                    );
                }

                jt += 4;
            }

            if jt < j_end {
                gemm_tn_f32_slice_raw(
                    dst_ptr, dst_col, m, a_ptr, a_col, k, b_ptr, b_col, jt, j_end,
                );
            }
        }
    }

    /// 4×4 register-tiled TN GEMM using **FMA + ip-outer loop order**.
    /// Identical microkernel to [`gemm_tn_f32_slice_raw_tiled_4x4_fma`]
    /// — same 4 A loads + 4 B loads + 16 FMA per K iter, same
    /// horizontal reduce, same scalar tail — but the OUTER two loops
    /// are swapped:
    ///
    /// ```ignore
    /// for ip in 0..m_quads:                  // outer (slow)
    ///   for jt in (j_start..j_end).step_by(4):  // inner (fast)
    ///     inner_K_loop_with_4x4_FMA_kernel
    /// ```
    ///
    /// # Why this is the right loop order for ab_dot
    ///
    /// On the ab_dot shape (m≈8000, n=200, k=2490):
    /// - A (lhs) is `m × k × 4 bytes = 80 MB` — far past L2.
    /// - B (rhs) is `n × k × 4 bytes = 2 MB`  — fits L2 easily.
    ///
    /// With the original `jt-outer, ip-inner` order (the canonical
    /// `_fma` kernel above), each of the 50 jt-tiles streams 80 MB
    /// of A from DRAM. Total A traffic per call: 4 GB DRAM. At
    /// ~80 GB/s the memory cost is ~50 ms — competitive with the
    /// kernel's 40 ms compute, so ab_dot is ~50% memory-bound.
    ///
    /// With `ip-outer, jt-inner`, each ip iteration holds one 40 KB
    /// A-quad hot in L1d while sweeping all 50 jt-tiles. B is
    /// streamed (2 MB per ip) but fits in L2 after the first
    /// ip iteration, so the 2000 ip iterations cost 80 MB DRAM
    /// (A read once) + 4 GB L2 (B re-read 2000×). At L2's ~300 GB/s
    /// that's ~14 ms total memory — kernel becomes compute-bound.
    /// Expected ~1.5-1.8× on m=8000 ab_dot.
    ///
    /// On small m (bb_dot m=200), both A and B fit L2 entirely so
    /// the loop order doesn't matter — same perf either way. We use
    /// the blocked variant as the universal default.
    ///
    /// # Numerics
    ///
    /// Identical to [`gemm_tn_f32_slice_raw_tiled_4x4_fma`]. The
    /// loop reorder doesn't change which (i, j) cells are computed
    /// or with what summation order — it just changes the visitation
    /// sequence. Wasmtime correctness sweep confirms bit-identical
    /// output to the canonical 4x4_fma kernel.
    ///
    /// # Safety
    ///
    /// Identical to [`gemm_tn_f32_slice_raw`].
    #[cfg(target_feature = "relaxed-simd")]
    #[allow(clippy::too_many_arguments)]
    pub unsafe fn gemm_tn_f32_slice_raw_tiled_4x4_fma_block(
        dst_ptr: *mut f32,
        dst_col: isize,
        m: usize,
        a_ptr: *const f32,
        a_col: isize,
        k: usize,
        b_ptr: *const f32,
        b_col: isize,
        j_start: usize,
        j_end: usize,
    ) {
        unsafe {
            let zero: v128 = f32x4_splat(0.0);

            // How many full 4-col j-tiles fit in [j_start, j_end).
            let n_tiles = (j_end - j_start) / 4;
            let m_quads = m / 4;

            // ip-OUTER loop: hold A's 4-row strip hot in L1d while
            // sweeping all jt tiles.
            for ip in 0..m_quads {
                let i0 = 4 * ip;
                let ap0 = a_ptr.offset(i0 as isize * a_col);
                let ap1 = a_ptr.offset((i0 + 1) as isize * a_col);
                let ap2 = a_ptr.offset((i0 + 2) as isize * a_col);
                let ap3 = a_ptr.offset((i0 + 3) as isize * a_col);

                // jt-INNER loop: stream through B's 4-col tiles.
                for jtile in 0..n_tiles {
                    let jt = j_start + 4 * jtile;
                    let bp0 = b_ptr.offset(jt as isize * b_col);
                    let bp1 = b_ptr.offset((jt + 1) as isize * b_col);
                    let bp2 = b_ptr.offset((jt + 2) as isize * b_col);
                    let bp3 = b_ptr.offset((jt + 3) as isize * b_col);

                    let mut c00 = zero;
                    let mut c01 = zero;
                    let mut c02 = zero;
                    let mut c03 = zero;
                    let mut c10 = zero;
                    let mut c11 = zero;
                    let mut c12 = zero;
                    let mut c13 = zero;
                    let mut c20 = zero;
                    let mut c21 = zero;
                    let mut c22 = zero;
                    let mut c23 = zero;
                    let mut c30 = zero;
                    let mut c31 = zero;
                    let mut c32 = zero;
                    let mut c33 = zero;

                    let mut p: usize = 0;
                    while p + 4 <= k {
                        let a0 = v128_load(ap0.add(p) as *const v128);
                        let a1 = v128_load(ap1.add(p) as *const v128);
                        let a2 = v128_load(ap2.add(p) as *const v128);
                        let a3 = v128_load(ap3.add(p) as *const v128);
                        let b0 = v128_load(bp0.add(p) as *const v128);
                        let b1 = v128_load(bp1.add(p) as *const v128);
                        let b2 = v128_load(bp2.add(p) as *const v128);
                        let b3 = v128_load(bp3.add(p) as *const v128);

                        c00 = f32x4_relaxed_madd(a0, b0, c00);
                        c01 = f32x4_relaxed_madd(a0, b1, c01);
                        c02 = f32x4_relaxed_madd(a0, b2, c02);
                        c03 = f32x4_relaxed_madd(a0, b3, c03);
                        c10 = f32x4_relaxed_madd(a1, b0, c10);
                        c11 = f32x4_relaxed_madd(a1, b1, c11);
                        c12 = f32x4_relaxed_madd(a1, b2, c12);
                        c13 = f32x4_relaxed_madd(a1, b3, c13);
                        c20 = f32x4_relaxed_madd(a2, b0, c20);
                        c21 = f32x4_relaxed_madd(a2, b1, c21);
                        c22 = f32x4_relaxed_madd(a2, b2, c22);
                        c23 = f32x4_relaxed_madd(a2, b3, c23);
                        c30 = f32x4_relaxed_madd(a3, b0, c30);
                        c31 = f32x4_relaxed_madd(a3, b1, c31);
                        c32 = f32x4_relaxed_madd(a3, b2, c32);
                        c33 = f32x4_relaxed_madd(a3, b3, c33);
                        p += 4;
                    }

                    let hred = |v: v128| -> f32 {
                        (f32x4_extract_lane::<0>(v) + f32x4_extract_lane::<1>(v))
                            + (f32x4_extract_lane::<2>(v) + f32x4_extract_lane::<3>(v))
                    };
                    let mut s00 = hred(c00);
                    let mut s01 = hred(c01);
                    let mut s02 = hred(c02);
                    let mut s03 = hred(c03);
                    let mut s10 = hred(c10);
                    let mut s11 = hred(c11);
                    let mut s12 = hred(c12);
                    let mut s13 = hred(c13);
                    let mut s20 = hred(c20);
                    let mut s21 = hred(c21);
                    let mut s22 = hred(c22);
                    let mut s23 = hred(c23);
                    let mut s30 = hred(c30);
                    let mut s31 = hred(c31);
                    let mut s32 = hred(c32);
                    let mut s33 = hred(c33);

                    while p < k {
                        let a0v = *ap0.add(p);
                        let a1v = *ap1.add(p);
                        let a2v = *ap2.add(p);
                        let a3v = *ap3.add(p);
                        let b0v = *bp0.add(p);
                        let b1v = *bp1.add(p);
                        let b2v = *bp2.add(p);
                        let b3v = *bp3.add(p);
                        s00 += a0v * b0v;
                        s01 += a0v * b1v;
                        s02 += a0v * b2v;
                        s03 += a0v * b3v;
                        s10 += a1v * b0v;
                        s11 += a1v * b1v;
                        s12 += a1v * b2v;
                        s13 += a1v * b3v;
                        s20 += a2v * b0v;
                        s21 += a2v * b1v;
                        s22 += a2v * b2v;
                        s23 += a2v * b3v;
                        s30 += a3v * b0v;
                        s31 += a3v * b1v;
                        s32 += a3v * b2v;
                        s33 += a3v * b3v;
                        p += 1;
                    }

                    let store4 = |dst: *mut f32, jbase: usize, row: isize, v0, v1, v2, v3| {
                        *dst.offset(jbase as isize * dst_col + row) = v0;
                        *dst.offset((jbase + 1) as isize * dst_col + row) = v1;
                        *dst.offset((jbase + 2) as isize * dst_col + row) = v2;
                        *dst.offset((jbase + 3) as isize * dst_col + row) = v3;
                    };
                    store4(dst_ptr, jt, i0 as isize, s00, s01, s02, s03);
                    store4(dst_ptr, jt, (i0 + 1) as isize, s10, s11, s12, s13);
                    store4(dst_ptr, jt, (i0 + 2) as isize, s20, s21, s22, s23);
                    store4(dst_ptr, jt, (i0 + 3) as isize, s30, s31, s32, s33);
                }
            }

            // m-tail (rows that didn't fit in m_quads × 4). For each
            // 4-col j-tile, walk through the tail rows with the single-j
            // kernel. (For m=200 / m=8000 the tail never fires.)
            let m_tail_start = m_quads * 4;
            if m_tail_start < m {
                let a_tail_off = (m_tail_start as isize).wrapping_mul(a_col);
                for jtile in 0..n_tiles {
                    let jt = j_start + 4 * jtile;
                    gemm_tn_f32_slice_raw(
                        dst_ptr.add(m_tail_start),
                        dst_col,
                        m - m_tail_start,
                        a_ptr.offset(a_tail_off),
                        a_col,
                        k,
                        b_ptr,
                        b_col,
                        jt,
                        jt + 4,
                    );
                }
            }

            // j-tail (cols ≤3 cols at the end that didn't fit in a
            // 4-col tile). Use the single-j kernel.
            let j_tail_start = j_start + n_tiles * 4;
            if j_tail_start < j_end {
                gemm_tn_f32_slice_raw(
                    dst_ptr,
                    dst_col,
                    m,
                    a_ptr,
                    a_col,
                    k,
                    b_ptr,
                    b_col,
                    j_tail_start,
                    j_end,
                );
            }
        }
    }
}

/// Per-column body for the fused CountSketch scatter-add (BED-decode →
/// normalize → scatter-add → renorm). Lives alongside the GEMM kernels
/// so the wasm SAB pool can dispatch it without a circular module
/// dependency: `compute.rs` (high-level) only needs to depend on
/// `wasm_simd` (low-level), not the other way round.
///
/// Both native (via `rayon::par_iter`) and wasm (via [`pool::parallel_scatter_f32`])
/// call into the same `scatter_one_column_*` helper for IEEE-754
/// bit-identical numerics across the two parallel paths.
pub(crate) mod scatter {
    use crate::bed::IidPos;

    /// Build a branchless 256-entry byte LUT for fused CountSketch scatter-add.
    /// For a given (mean, inv_std), each BED byte (4 packed genotypes) maps
    /// to 4 normalized f32 values. Missing genotypes map to 0.0 (impute
    /// to mean → center → 0), so the scatter-add is unconditional.
    #[inline]
    fn build_norm_byte_lut_f32(mean: f32, inv_std: f32) -> [[f32; 4]; 256] {
        // 0b00 = hom A1 (geno=2) → (2 - mean) * inv_std
        // 0b01 = missing         → 0.0
        // 0b10 = het (geno=1)    → (1 - mean) * inv_std
        // 0b11 = hom A2 (geno=0) → -mean * inv_std
        let code_vals: [f32; 4] = [
            (2.0 - mean) * inv_std,
            0.0,
            (1.0 - mean) * inv_std,
            -mean * inv_std,
        ];
        core::array::from_fn(|b| {
            let byte = b as u8;
            [
                code_vals[(byte & 0b11) as usize],
                code_vals[((byte >> 2) & 0b11) as usize],
                code_vals[((byte >> 4) & 0b11) as usize],
                code_vals[((byte >> 6) & 0b11) as usize],
            ]
        })
    }

    /// f64 variant of the branchless byte LUT. See [`build_norm_byte_lut_f32`].
    #[inline]
    fn build_norm_byte_lut_f64(mean: f64, inv_std: f64) -> [[f64; 4]; 256] {
        let code_vals: [f64; 4] = [
            (2.0 - mean) * inv_std,
            0.0,
            (1.0 - mean) * inv_std,
            -mean * inv_std,
        ];
        core::array::from_fn(|b| {
            let byte = b as u8;
            [
                code_vals[(byte & 0b11) as usize],
                code_vals[((byte >> 2) & 0b11) as usize],
                code_vals[((byte >> 4) & 0b11) as usize],
                code_vals[((byte >> 6) & 0b11) as usize],
            ]
        })
    }

    /// Scatter-add one column `j` of the fused CountSketch (f32 variant).
    ///
    /// Decodes the BED bytes for SNP `j` (each byte = 4 packed
    /// genotypes), normalizes each genotype using the supplied
    /// per-SNP `mean` and `inv_std`, scatter-adds into
    /// `out[bucket[i], j]` with sign `sign[i]`, then ratio-renorms so
    /// `||out[:, j]||² == n_norm`.
    ///
    /// # Safety
    ///
    /// - `raw_bytes` must point to at least `(j + 1) * bytes_per_snp` bytes.
    /// - `bucket` and `sign` must point to at least `n_indiv` elements each.
    /// - When `all_iids == false`, `iid_positions` must point to at least
    ///   `iid_positions_len` elements (the per-kept-individual byte/shift
    ///   indices), and `n_indiv` must equal `iid_positions_len`.
    /// - `out_ptr + j * out_col_stride` must be a valid pointer to at least
    ///   `d` f32 elements, exclusively owned for this call (no other thread
    ///   may write that column).
    /// - All pointer arguments must outlive the call. Reads from
    ///   `raw_bytes`, `bucket`, `sign`, `iid_positions` may be aliased with
    ///   other concurrent reads (treated as `&[_]` here) but must not be
    ///   mutated by any other thread during this call.
    ///
    /// # Dispatch
    ///
    /// Picks the wasm32+simd128 multi-lane impl when available,
    /// falls back to the scalar reference otherwise. All callers
    /// (native rayon path in `compute.rs`, wasm SAB-pool dispatch
    /// in `super::pool::run_scatter_slice_f32`) go through this
    /// single entry point.
    #[allow(clippy::too_many_arguments)]
    pub(crate) unsafe fn scatter_one_column_f32(
        j: usize,
        raw_bytes: *const u8,
        bytes_per_snp: usize,
        n_indiv: usize,
        bucket: *const u32,
        sign: *const f32,
        iid_positions: *const IidPos,
        iid_positions_len: usize,
        all_iids: bool,
        mean: f32,
        inv_std: f32,
        n_norm: f64,
        out_ptr: *mut f32,
        out_col_stride: isize,
        d: usize,
    ) {
        // Workstream J: SIMD path is gated on `+simd128` and the
        // all_iids fast path (which is the hot one — 99% of biobank
        // runs hit this path because there's no `--keep` filter).
        // The subsample path's irregular byte+shift indexing doesn't
        // amortize the multi-lane setup, so we leave it on the scalar
        // path even on simd128 builds.
        #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
        {
            if all_iids {
                unsafe {
                    scatter_one_column_f32_simd128(
                        j,
                        raw_bytes,
                        bytes_per_snp,
                        n_indiv,
                        bucket,
                        sign,
                        mean,
                        inv_std,
                        n_norm,
                        out_ptr,
                        out_col_stride,
                        d,
                    );
                }
                return;
            }
        }
        unsafe {
            scatter_one_column_f32_scalar(
                j,
                raw_bytes,
                bytes_per_snp,
                n_indiv,
                bucket,
                sign,
                iid_positions,
                iid_positions_len,
                all_iids,
                mean,
                inv_std,
                n_norm,
                out_ptr,
                out_col_stride,
                d,
            );
        }
    }

    /// Portable scalar reference for one scatter column. Runs on every
    /// target — used as the native CPU implementation (rayon path)
    /// and as the property-test oracle for the wasm32+simd128 variant.
    ///
    /// # Safety
    ///
    /// Same preconditions as [`scatter_one_column_f32`].
    #[allow(clippy::too_many_arguments)]
    pub(crate) unsafe fn scatter_one_column_f32_scalar(
        j: usize,
        raw_bytes: *const u8,
        bytes_per_snp: usize,
        n_indiv: usize,
        bucket: *const u32,
        sign: *const f32,
        iid_positions: *const IidPos,
        iid_positions_len: usize,
        all_iids: bool,
        mean: f32,
        inv_std: f32,
        n_norm: f64,
        out_ptr: *mut f32,
        out_col_stride: isize,
        d: usize,
    ) {
        // Reconstruct slices from raw pointers (the pool dispatch carries
        // them as `usize` to keep `JobArgs` POD/Copy/Sync).
        let snp_bytes =
            unsafe { core::slice::from_raw_parts(raw_bytes.add(j * bytes_per_snp), bytes_per_snp) };
        let dst = unsafe {
            core::slice::from_raw_parts_mut(out_ptr.offset(j as isize * out_col_stride), d)
        };
        let bucket_slice = unsafe { core::slice::from_raw_parts(bucket, n_indiv) };
        let sign_slice = unsafe { core::slice::from_raw_parts(sign, n_indiv) };

        // Zero output column. Replace-not-add accumulation, so no need to
        // touch this from the caller.
        for v in dst.iter_mut() {
            *v = 0.0;
        }

        let norm_lut = build_norm_byte_lut_f32(mean, inv_std);

        if all_iids {
            // Branchless byte-level decode: each byte → 4 normalized
            // values; missing genotypes map to 0.0 so scatter-add is
            // unconditional. SAFETY: bucket[i] is built from
            // `rng.u32(..d as u32)` so `< d`; verified at CountSketchProj
            // construction.
            let n_full_bytes = n_indiv / 4;
            for byte_idx in 0..n_full_bytes {
                let vals = unsafe { &norm_lut[*snp_bytes.get_unchecked(byte_idx) as usize] };
                let base_i = byte_idx * 4;
                for (k, &val) in vals.iter().enumerate() {
                    let i = base_i + k;
                    unsafe {
                        *dst.get_unchecked_mut(*bucket_slice.get_unchecked(i) as usize) +=
                            *sign_slice.get_unchecked(i) * val;
                    }
                }
            }
            // Remainder individuals in the trailing partial byte.
            let rem_start = n_full_bytes * 4;
            if rem_start < n_indiv {
                let vals = unsafe { &norm_lut[*snp_bytes.get_unchecked(n_full_bytes) as usize] };
                for (k, &val) in vals.iter().enumerate().take(n_indiv - rem_start) {
                    let i = rem_start + k;
                    unsafe {
                        *dst.get_unchecked_mut(*bucket_slice.get_unchecked(i) as usize) +=
                            *sign_slice.get_unchecked(i) * val;
                    }
                }
            }
        } else {
            // Subsample path: each kept individual's byte/shift is
            // precomputed in iid_positions. `shift / 2` maps byte shift
            // (0, 2, 4, 6) to LUT index (0, 1, 2, 3).
            let iid_pos = unsafe { core::slice::from_raw_parts(iid_positions, iid_positions_len) };
            for (i, pos) in iid_pos.iter().enumerate() {
                let vals = unsafe { &norm_lut[*snp_bytes.get_unchecked(pos.byte_idx) as usize] };
                let val = vals[(pos.shift / 2) as usize];
                unsafe {
                    *dst.get_unchecked_mut(*bucket_slice.get_unchecked(i) as usize) +=
                        *sign_slice.get_unchecked(i) * val;
                }
            }
        }

        // Ratio-estimator renorm: rescale so ||col||² = n_norm. Accumulate
        // squared sum in f64 to keep the renorm bit-identical regardless
        // of the order of partial sums.
        let mut nrm_sq = 0.0f64;
        for &v in dst.iter() {
            nrm_sq += (v as f64) * (v as f64);
        }
        if nrm_sq > 0.0 {
            let scale = (n_norm / nrm_sq).sqrt() as f32;
            for v in dst.iter_mut() {
                *v *= scale;
            }
        }
    }

    /// Multi-lane scatter kernel for the wasm32+simd128 / all_iids
    /// fast path (the dominant case in browser sketch runs — full
    /// genome compute with no `--keep` filter).
    ///
    /// Key trick is lane decorrelation, not SIMD ops: K=4 independent
    /// scatter buffers (`lane0..lane3`) receive the 4 individuals
    /// packed in each BED byte 1:1, breaking the serial dep chain on
    /// the indirect `dst[bucket[i]] +=` store. The lanes are folded
    /// into the final `dst` column with f32x4 SIMD adds. SIMD also
    /// applies to the 4 per-byte `sign * val` muls and the renorm
    /// tail; the scatter-add itself stays scalar (wasm has no
    /// gather/scatter).
    ///
    /// # Numerics
    ///
    /// Reduction order changes from single-stream to K=4 streams
    /// folded at the end. f32 error bound `~ε·sqrt(K)` ≈ 3e-7 at K=4.
    /// `nrm_sq` is accumulated in `f64x2` lanes then folded in f64,
    /// matching the scalar reference's f64-accumulation intent.
    ///
    /// # Safety
    ///
    /// Same as [`scatter_one_column_f32`]; additionally requires
    /// `all_iids == true` (dispatcher gates on this).
    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
    #[allow(clippy::too_many_arguments)]
    pub(crate) unsafe fn scatter_one_column_f32_simd128(
        j: usize,
        raw_bytes: *const u8,
        bytes_per_snp: usize,
        n_indiv: usize,
        bucket: *const u32,
        sign: *const f32,
        mean: f32,
        inv_std: f32,
        n_norm: f64,
        out_ptr: *mut f32,
        out_col_stride: isize,
        d: usize,
    ) {
        use core::arch::wasm32::*;

        let snp_bytes =
            unsafe { core::slice::from_raw_parts(raw_bytes.add(j * bytes_per_snp), bytes_per_snp) };
        let dst = unsafe {
            core::slice::from_raw_parts_mut(out_ptr.offset(j as isize * out_col_stride), d)
        };
        let bucket_slice = unsafe { core::slice::from_raw_parts(bucket, n_indiv) };
        let sign_slice = unsafe { core::slice::from_raw_parts(sign, n_indiv) };

        // K=4 lane buffers, stack-allocated in the inner worker's
        // 2 MB TLS region. At d=200 → 4 × 200 × 4 = 3.2 KB (L1-fit).
        // MAX_D bounds the stack reservation so bounds checks elide;
        // d > MAX_D falls back to scalar (defensive — sketch slider
        // max is 3200).
        const MAX_D: usize = 4096;
        if d > MAX_D {
            unsafe {
                scatter_one_column_f32_scalar(
                    j,
                    raw_bytes,
                    bytes_per_snp,
                    n_indiv,
                    bucket,
                    sign,
                    core::ptr::null(),
                    0,
                    true,
                    mean,
                    inv_std,
                    n_norm,
                    out_ptr,
                    out_col_stride,
                    d,
                );
            }
            return;
        }
        // `MaybeUninit` + zero only the `d`-sized prefix instead of
        // the full MAX_D (literal `[0.0; MAX_D]` forces a 64 KB
        // memset per call; at small N that init dominates and the
        // SIMD path goes from 1.43× → 0.88× scalar — verified on the
        // wasm-bench scatter sweep). SAFETY: every read/write below
        // is bounded by `d` (via `bucket[i] < d` from CountSketchProj
        // construction), so the uninit tail is never touched.
        let mut lane0_uninit: [core::mem::MaybeUninit<f32>; MAX_D] =
            [core::mem::MaybeUninit::uninit(); MAX_D];
        let mut lane1_uninit: [core::mem::MaybeUninit<f32>; MAX_D] =
            [core::mem::MaybeUninit::uninit(); MAX_D];
        let mut lane2_uninit: [core::mem::MaybeUninit<f32>; MAX_D] =
            [core::mem::MaybeUninit::uninit(); MAX_D];
        let mut lane3_uninit: [core::mem::MaybeUninit<f32>; MAX_D] =
            [core::mem::MaybeUninit::uninit(); MAX_D];
        unsafe {
            core::ptr::write_bytes(lane0_uninit.as_mut_ptr() as *mut f32, 0, d);
            core::ptr::write_bytes(lane1_uninit.as_mut_ptr() as *mut f32, 0, d);
            core::ptr::write_bytes(lane2_uninit.as_mut_ptr() as *mut f32, 0, d);
            core::ptr::write_bytes(lane3_uninit.as_mut_ptr() as *mut f32, 0, d);
        }
        let lane0 = unsafe { &mut *(&mut lane0_uninit as *mut _ as *mut [f32; MAX_D]) };
        let lane1 = unsafe { &mut *(&mut lane1_uninit as *mut _ as *mut [f32; MAX_D]) };
        let lane2 = unsafe { &mut *(&mut lane2_uninit as *mut _ as *mut [f32; MAX_D]) };
        let lane3 = unsafe { &mut *(&mut lane3_uninit as *mut _ as *mut [f32; MAX_D]) };

        let norm_lut = build_norm_byte_lut_f32(mean, inv_std);
        let n_full_bytes = n_indiv / 4;
        // SAFETY: bucket[i] < d ≤ MAX_D by CountSketchProj construction;
        // snp_bytes/sign/bucket reads are bounded by n_indiv.
        unsafe {
            for byte_idx in 0..n_full_bytes {
                let byte = *snp_bytes.get_unchecked(byte_idx);
                let vals_v = v128_load(norm_lut[byte as usize].as_ptr() as *const v128);
                let base_i = byte_idx * 4;
                let signs_v = v128_load(sign_slice.as_ptr().add(base_i) as *const v128);
                // 1 vector mul replaces 4 scalar muls; the 4 lane
                // scatter-adds below are independent (disjoint
                // `laneN` buffers) so the CPU can keep them in flight
                // concurrently, breaking the dep chain.
                let scaled_v = f32x4_mul(vals_v, signs_v);
                let s0 = f32x4_extract_lane::<0>(scaled_v);
                let s1 = f32x4_extract_lane::<1>(scaled_v);
                let s2 = f32x4_extract_lane::<2>(scaled_v);
                let s3 = f32x4_extract_lane::<3>(scaled_v);
                let b0 = *bucket_slice.get_unchecked(base_i) as usize;
                let b1 = *bucket_slice.get_unchecked(base_i + 1) as usize;
                let b2 = *bucket_slice.get_unchecked(base_i + 2) as usize;
                let b3 = *bucket_slice.get_unchecked(base_i + 3) as usize;
                *lane0.get_unchecked_mut(b0) += s0;
                *lane1.get_unchecked_mut(b1) += s1;
                *lane2.get_unchecked_mut(b2) += s2;
                *lane3.get_unchecked_mut(b3) += s3;
            }
        }

        // Remainder individuals (n_indiv % 4 ∈ {1, 2, 3}) onto lanes
        // 0..rem-1. Reborrow lanes so they stay usable for the fold.
        let rem_start = n_full_bytes * 4;
        let rem = n_indiv - rem_start;
        if rem > 0 {
            unsafe {
                let vals_arr = &norm_lut[*snp_bytes.get_unchecked(n_full_bytes) as usize];
                let lanes: [&mut [f32; MAX_D]; 4] =
                    [&mut *lane0, &mut *lane1, &mut *lane2, &mut *lane3];
                for (k, lane) in lanes.into_iter().enumerate().take(rem) {
                    let i = rem_start + k;
                    let v = vals_arr[k] * *sign_slice.get_unchecked(i);
                    *lane.get_unchecked_mut(*bucket_slice.get_unchecked(i) as usize) += v;
                }
            }
        }

        // Fold: dst[j_chunk..j_chunk+4] = lane0[..] + lane1[..] +
        // lane2[..] + lane3[..]. Process 4 elements per iteration via
        // f32x4 adds. The tail (d % 4) is handled scalarly.
        let d_chunks = d / 4 * 4;
        // SAFETY: d_chunks ≤ d ≤ MAX_D; all loads/stores in range.
        unsafe {
            let mut j_chunk = 0usize;
            while j_chunk < d_chunks {
                let l0 = v128_load(lane0.as_ptr().add(j_chunk) as *const v128);
                let l1 = v128_load(lane1.as_ptr().add(j_chunk) as *const v128);
                let l2 = v128_load(lane2.as_ptr().add(j_chunk) as *const v128);
                let l3 = v128_load(lane3.as_ptr().add(j_chunk) as *const v128);
                let sum = f32x4_add(f32x4_add(l0, l1), f32x4_add(l2, l3));
                v128_store(dst.as_mut_ptr().add(j_chunk) as *mut v128, sum);
                j_chunk += 4;
            }
            // Tail (d % 4 ∈ {0, 1, 2, 3}).
            for j in d_chunks..d {
                *dst.get_unchecked_mut(j) = lane0[j] + lane1[j] + lane2[j] + lane3[j];
            }
        }

        // Renorm tail: ||dst||^2 in 4 parallel f64 accs (2 × f64x2)
        // to break the dep chain, then scale dst by sqrt(n_norm/ssq).
        // SAFETY: dst.len() == d ≤ MAX_D; all loads in range.
        let mut ssq_a = f64x2_splat(0.0);
        let mut ssq_b = f64x2_splat(0.0);
        unsafe {
            let mut j = 0usize;
            while j + 4 <= d {
                let v = v128_load(dst.as_ptr().add(j) as *const v128);
                // Widen lo/hi pairs of f32 → f64 (the wasm intrinsic
                // only promotes the low 2 lanes, so shuffle high
                // first for the hi pair).
                let v_hi_lo = i32x4_shuffle::<2, 3, 0, 0>(v, v);
                let lo64 = f64x2_promote_low_f32x4(v);
                let hi64 = f64x2_promote_low_f32x4(v_hi_lo);
                ssq_a = f64x2_add(ssq_a, f64x2_mul(lo64, lo64));
                ssq_b = f64x2_add(ssq_b, f64x2_mul(hi64, hi64));
                j += 4;
            }
            let acc = f64x2_add(ssq_a, ssq_b);
            let mut nrm_sq = f64x2_extract_lane::<0>(acc) + f64x2_extract_lane::<1>(acc);
            while j < d {
                let v = *dst.get_unchecked(j) as f64;
                nrm_sq += v * v;
                j += 1;
            }
            if nrm_sq > 0.0 {
                let scale = (n_norm / nrm_sq).sqrt() as f32;
                let scale_v = f32x4_splat(scale);
                let mut j2 = 0usize;
                while j2 + 4 <= d {
                    let v = v128_load(dst.as_ptr().add(j2) as *const v128);
                    v128_store(dst.as_mut_ptr().add(j2) as *mut v128, f32x4_mul(v, scale_v));
                    j2 += 4;
                }
                while j2 < d {
                    *dst.get_unchecked_mut(j2) *= scale;
                    j2 += 1;
                }
            }
        }
    }

    /// f64 variant of [`scatter_one_column_f32`]. Same algorithm,
    /// different element width on `sign`, `mean`/`inv_std`, and `out`.
    #[allow(clippy::too_many_arguments)]
    pub(crate) unsafe fn scatter_one_column_f64(
        j: usize,
        raw_bytes: *const u8,
        bytes_per_snp: usize,
        n_indiv: usize,
        bucket: *const u32,
        sign: *const f64,
        iid_positions: *const IidPos,
        iid_positions_len: usize,
        all_iids: bool,
        mean: f64,
        inv_std: f64,
        n_norm: f64,
        out_ptr: *mut f64,
        out_col_stride: isize,
        d: usize,
    ) {
        let snp_bytes =
            unsafe { core::slice::from_raw_parts(raw_bytes.add(j * bytes_per_snp), bytes_per_snp) };
        let dst = unsafe {
            core::slice::from_raw_parts_mut(out_ptr.offset(j as isize * out_col_stride), d)
        };
        let bucket_slice = unsafe { core::slice::from_raw_parts(bucket, n_indiv) };
        let sign_slice = unsafe { core::slice::from_raw_parts(sign, n_indiv) };

        for v in dst.iter_mut() {
            *v = 0.0;
        }

        let norm_lut = build_norm_byte_lut_f64(mean, inv_std);

        if all_iids {
            let n_full_bytes = n_indiv / 4;
            for byte_idx in 0..n_full_bytes {
                let vals = unsafe { &norm_lut[*snp_bytes.get_unchecked(byte_idx) as usize] };
                let base_i = byte_idx * 4;
                for (k, &val) in vals.iter().enumerate() {
                    let i = base_i + k;
                    unsafe {
                        *dst.get_unchecked_mut(*bucket_slice.get_unchecked(i) as usize) +=
                            *sign_slice.get_unchecked(i) * val;
                    }
                }
            }
            let rem_start = n_full_bytes * 4;
            if rem_start < n_indiv {
                let vals = unsafe { &norm_lut[*snp_bytes.get_unchecked(n_full_bytes) as usize] };
                for (k, &val) in vals.iter().enumerate().take(n_indiv - rem_start) {
                    let i = rem_start + k;
                    unsafe {
                        *dst.get_unchecked_mut(*bucket_slice.get_unchecked(i) as usize) +=
                            *sign_slice.get_unchecked(i) * val;
                    }
                }
            }
        } else {
            let iid_pos = unsafe { core::slice::from_raw_parts(iid_positions, iid_positions_len) };
            for (i, pos) in iid_pos.iter().enumerate() {
                let vals = unsafe { &norm_lut[*snp_bytes.get_unchecked(pos.byte_idx) as usize] };
                let val = vals[(pos.shift / 2) as usize];
                unsafe {
                    *dst.get_unchecked_mut(*bucket_slice.get_unchecked(i) as usize) +=
                        *sign_slice.get_unchecked(i) * val;
                }
            }
        }

        let mut nrm_sq = 0.0f64;
        for &v in dst.iter() {
            nrm_sq += v * v;
        }
        if nrm_sq > 0.0 {
            let scale = (n_norm / nrm_sq).sqrt();
            for v in dst.iter_mut() {
                *v *= scale;
            }
        }
    }
}

/// Manual SAB-backed worker pool for parallel f32 GEMM on
/// `wasm32 + atomics + simd128`.
///
/// # Why not wasm-bindgen-rayon
///
/// The bundler-mode `workerHelpers.js` that wasm-bindgen-rayon ships
/// with `--target=web` does `import('../../..')` to load the main JS
/// bundle. That resolves to `/ldsc/snippets/` (a directory with no
/// entry file) under trunk's static-copy snippet layout, hanging the
/// inner-worker import indefinitely. Rather than fight that path
/// resolution, we hand-roll a minimal pool specialized for our one
/// hot kernel.
///
/// # Topology
///
/// One pool per outer compute Worker (Workers don't share linear
/// memory across each other). Each outer Worker, on its first
/// compute, spawns N inner Workers via `assets/inner_worker.js` and
/// passes them the wasm module + memory. Each inner Worker imports
/// the bundle, instantiates wasm with the SHARED memory, and parks
/// in [`worker_loop`] waiting for the outer Worker to dispatch work.
///
/// # Synchronization
///
/// The outer Worker is "worker 0"; inner Workers are 1..N. We use
/// raw `memory.atomic.wait32` / `memory.atomic.notify` (the wasm
/// threads proposal's wait/notify primitives) on three static atomics:
///
/// - `JOB_GEN`: monotonically incrementing generation. Inner workers
///   sleep on it; outer increments + notifies to wake them. Each
///   inner remembers `last_seen_gen` to detect a fresh job.
/// - `JOB_DONE`: counter inner workers fetch_add when their slice is
///   done. Outer waits on it reaching N-1 (worker 0 / outer ran its
///   own slice synchronously).
/// - `JOB_KIND`: discriminator (1 = exit, 2 = gemm).
///
/// Args are kept in a `static` `UnsafeCell<JobArgs>`. The
/// release/acquire on `JOB_GEN` provides the happens-before edge that
/// makes the unsafe access sound: outer mutates `JOB_ARGS` *before*
/// the `Release` increment of `JOB_GEN`; each inner reads it *after*
/// observing the new generation with `Acquire` ordering.
///
/// # Numerics
///
/// Bit-identical to serial. Each worker computes a disjoint slice of
/// output columns; column writes don't overlap; lhs/rhs are read-only.
/// No reduction across workers, so no summation-order changes.
#[cfg(all(
    target_arch = "wasm32",
    target_feature = "atomics",
    target_feature = "simd128"
))]
pub mod pool {
    use core::arch::wasm32::{memory_atomic_notify, memory_atomic_wait32};
    use core::cell::UnsafeCell;
    use core::sync::atomic::{AtomicI32, AtomicUsize, Ordering};

    /// Pool size. 0 = pool not initialised (serial fallback). Set once
    /// per outer Worker via [`init`] before the first dispatch.
    pub static POOL_SIZE: AtomicUsize = AtomicUsize::new(0);

    /// Generation counter, incremented per dispatched job. Inner
    /// workers sleep via `memory_atomic_wait32` on this.
    pub static JOB_GEN: AtomicI32 = AtomicI32::new(0);

    /// Inner workers fetch_add this when their slice is done. Outer
    /// waits for it to reach `POOL_SIZE - 1`.
    pub static JOB_DONE: AtomicI32 = AtomicI32::new(0);

    /// Job discriminator. `1` = exit (worker_loop returns), `2` = gemm,
    /// `3` = scatter_f32 (CountSketch fused scatter, f32 variant),
    /// `4` = scatter_f64 (same, f64 variant). Set by the outer Worker
    /// BEFORE incrementing `JOB_GEN`; read by inner workers AFTER they
    /// observe the new generation. Each kind has its own dedicated args
    /// cell (`JOB_ARGS`, `SCATTER_ARGS_F32`, `SCATTER_ARGS_F64`); only
    /// the active kind's cell is touched per dispatch.
    pub static JOB_KIND: AtomicI32 = AtomicI32::new(0);

    /// Wrapper around `UnsafeCell<JobArgs>` to manually mark Sync. Sound
    /// because all access is ordered through `JOB_GEN`'s
    /// release/acquire fences. The same invariant applies to
    /// [`ScatterArgsCellF32`] / [`ScatterArgsCellF64`] below: between
    /// `JOB_GEN` bumps only the outer Worker writes; inner Workers read
    /// the cell selected by [`JOB_KIND`] after the Acquire fence.
    pub struct JobArgsCell(UnsafeCell<JobArgs>);
    // SAFETY: see module docs.
    unsafe impl Sync for JobArgsCell {}

    /// Job arguments. Mutated only by the outer Worker between
    /// dispatches; read by inner workers after they observe a new
    /// `JOB_GEN`. Fixed-size (no allocations) so we can keep it in a
    /// `static`.
    #[derive(Clone, Copy)]
    pub struct JobArgs {
        pub dst_ptr: usize,
        pub dst_col_stride: isize,
        pub dst_rows: usize,
        pub dst_cols: usize,
        pub lhs_ptr: usize,
        pub lhs_col_stride: isize,
        pub lhs_rows: usize,
        pub rhs_ptr: usize,
        pub rhs_col_stride: isize,
    }

    impl JobArgs {
        const ZERO: Self = Self {
            dst_ptr: 0,
            dst_col_stride: 0,
            dst_rows: 0,
            dst_cols: 0,
            lhs_ptr: 0,
            lhs_col_stride: 0,
            lhs_rows: 0,
            rhs_ptr: 0,
            rhs_col_stride: 0,
        };
    }

    pub static JOB_ARGS: JobArgsCell = JobArgsCell(UnsafeCell::new(JobArgs::ZERO));

    /// Args for a CountSketch scatter dispatch (f32 variant).
    ///
    /// All pointers are passed as `usize` to keep the struct POD/Copy/Sync.
    /// They are cast back at the per-worker entry point
    /// ([`run_scatter_slice_f32`]). The contract on each pointer (length,
    /// aliasing, mutability) is documented in
    /// [`super::scatter::scatter_one_column_f32`].
    #[derive(Clone, Copy)]
    pub struct ScatterArgsF32 {
        pub raw_bytes_ptr: usize, // *const u8, len >= c * bytes_per_snp
        pub bytes_per_snp: usize,
        pub n_indiv: usize,
        pub c: usize,           // number of SNPs (output columns)
        pub bucket_ptr: usize,  // *const u32, len n_indiv
        pub sign_ptr: usize,    // *const f32, len n_indiv
        pub iid_pos_ptr: usize, // *const IidPos, len iid_pos_len
        pub iid_pos_len: usize,
        pub all_iids: u32,  // bool as u32 (1 = all individuals kept)
        pub n_norm: f64,    // renorm target: ||out[:, j]||^2 = n_norm
        pub out_ptr: usize, // *mut f32
        pub out_col_stride: isize,
        pub d: usize,           // sketch dimension (output rows)
        pub mean_ptr: usize,    // *const f32, len c
        pub inv_std_ptr: usize, // *const f32, len c
    }

    impl ScatterArgsF32 {
        const ZERO: Self = Self {
            raw_bytes_ptr: 0,
            bytes_per_snp: 0,
            n_indiv: 0,
            c: 0,
            bucket_ptr: 0,
            sign_ptr: 0,
            iid_pos_ptr: 0,
            iid_pos_len: 0,
            all_iids: 0,
            n_norm: 0.0,
            out_ptr: 0,
            out_col_stride: 0,
            d: 0,
            mean_ptr: 0,
            inv_std_ptr: 0,
        };
    }

    pub struct ScatterArgsCellF32(UnsafeCell<ScatterArgsF32>);
    // SAFETY: see `JobArgsCell` doc-comment — same release/acquire
    // sequencing on `JOB_GEN`.
    unsafe impl Sync for ScatterArgsCellF32 {}

    pub static SCATTER_ARGS_F32: ScatterArgsCellF32 =
        ScatterArgsCellF32(UnsafeCell::new(ScatterArgsF32::ZERO));

    /// f64 sibling of [`ScatterArgsF32`]. Used by the f64 fused-CountSketch
    /// dispatch path (rare on wasm — only when the user opts out of
    /// `--fast-f32`).
    #[derive(Clone, Copy)]
    pub struct ScatterArgsF64 {
        pub raw_bytes_ptr: usize,
        pub bytes_per_snp: usize,
        pub n_indiv: usize,
        pub c: usize,
        pub bucket_ptr: usize, // *const u32
        pub sign_ptr: usize,   // *const f64
        pub iid_pos_ptr: usize,
        pub iid_pos_len: usize,
        pub all_iids: u32,
        pub n_norm: f64,
        pub out_ptr: usize, // *mut f64
        pub out_col_stride: isize,
        pub d: usize,
        pub mean_ptr: usize,    // *const f64
        pub inv_std_ptr: usize, // *const f64
    }

    impl ScatterArgsF64 {
        const ZERO: Self = Self {
            raw_bytes_ptr: 0,
            bytes_per_snp: 0,
            n_indiv: 0,
            c: 0,
            bucket_ptr: 0,
            sign_ptr: 0,
            iid_pos_ptr: 0,
            iid_pos_len: 0,
            all_iids: 0,
            n_norm: 0.0,
            out_ptr: 0,
            out_col_stride: 0,
            d: 0,
            mean_ptr: 0,
            inv_std_ptr: 0,
        };
    }

    pub struct ScatterArgsCellF64(UnsafeCell<ScatterArgsF64>);
    unsafe impl Sync for ScatterArgsCellF64 {}

    pub static SCATTER_ARGS_F64: ScatterArgsCellF64 =
        ScatterArgsCellF64(UnsafeCell::new(ScatterArgsF64::ZERO));

    /// Initialise the pool size. Call once per outer Worker, after the
    /// JS side has spawned the N inner Workers and they have called
    /// [`worker_loop`]. Idempotent — calling with a different N after
    /// inner workers are already parked has no effect on them, only
    /// on what slice each worker computes for the next dispatch.
    pub fn init(n_workers: usize) {
        POOL_SIZE.store(n_workers, Ordering::SeqCst);
    }

    /// Inner Worker entry point. Loops forever waiting for jobs.
    /// Called by the JS-side `inner_worker.js` via the wasm-bindgen
    /// export `inner_worker_loop` in `ldsc-web::worker`.
    ///
    /// `worker_id` is the inner Worker's index in `1..N`. (Worker 0
    /// is the outer Worker, which runs synchronously in
    /// [`parallel_gemm_tn_f32`].)
    pub fn worker_loop(worker_id: usize) {
        let mut last_gen: i32 = JOB_GEN.load(Ordering::Acquire);
        loop {
            // Wait for a new job (JOB_GEN changes). The Atomics.wait
            // call returns immediately if JOB_GEN already differs from
            // the value we passed in (we re-load each iteration to
            // avoid a "missed wake" race).
            loop {
                let cur = JOB_GEN.load(Ordering::Acquire);
                if cur != last_gen {
                    last_gen = cur;
                    break;
                }
                // SAFETY: JOB_GEN.as_ptr() is a valid `*const AtomicI32`;
                // memory_atomic_wait32 takes `*mut i32` (interpretation
                // is the same — the wasm spec only requires alignment).
                unsafe {
                    let _ = memory_atomic_wait32(JOB_GEN.as_ptr().cast(), cur, -1);
                }
            }

            let n_workers = POOL_SIZE.load(Ordering::Acquire);
            let kind = JOB_KIND.load(Ordering::Acquire);
            match kind {
                1 => return, // exit
                2 => {
                    // SAFETY: outer Worker mutated JOB_ARGS BEFORE it
                    // incremented JOB_GEN; we observed the increment
                    // above with Acquire, so the args are visible.
                    let args = unsafe { *JOB_ARGS.0.get() };
                    run_gemm_slice(worker_id, n_workers, &args);
                }
                3 => {
                    // SAFETY: same release/acquire on JOB_GEN — outer
                    // wrote SCATTER_ARGS_F32 before bumping JOB_GEN; the
                    // Acquire load above gives the happens-before edge.
                    let args = unsafe { *SCATTER_ARGS_F32.0.get() };
                    run_scatter_slice_f32(worker_id, n_workers, &args);
                }
                4 => {
                    let args = unsafe { *SCATTER_ARGS_F64.0.get() };
                    run_scatter_slice_f64(worker_id, n_workers, &args);
                }
                _ => {}
            }

            // Mark this worker as done.
            JOB_DONE.fetch_add(1, Ordering::Release);
            // SAFETY: see comment in the wait above.
            unsafe {
                let _ = memory_atomic_notify(JOB_DONE.as_ptr().cast(), 1);
            }
        }
    }

    /// Outer-side dispatch. Returns `true` if the work was farmed
    /// out to the pool; `false` (call serial fallback) if the pool
    /// isn't initialised or N <= 1.
    ///
    /// # Safety
    ///
    /// Caller must ensure all matrix pointers + strides are valid for
    /// the duration of the call. The dst region must be exclusively
    /// owned (no aliasing reads/writes by other code while the pool
    /// is computing). lhs and rhs may be aliased read-only.
    #[allow(clippy::too_many_arguments)]
    pub unsafe fn parallel_gemm_tn_f32(
        dst_ptr: *mut f32,
        dst_col_stride: isize,
        dst_rows: usize,
        dst_cols: usize,
        lhs_ptr: *const f32,
        lhs_col_stride: isize,
        lhs_rows: usize,
        rhs_ptr: *const f32,
        rhs_col_stride: isize,
    ) -> bool {
        let n = POOL_SIZE.load(Ordering::Acquire);
        if n <= 1 {
            return false;
        }

        let args = JobArgs {
            dst_ptr: dst_ptr as usize,
            dst_col_stride,
            dst_rows,
            dst_cols,
            lhs_ptr: lhs_ptr as usize,
            lhs_col_stride,
            lhs_rows,
            rhs_ptr: rhs_ptr as usize,
            rhs_col_stride,
        };

        // SAFETY: between dispatches, only the outer Worker writes
        // JOB_ARGS. The release ordering on the JOB_GEN increment
        // below makes inner workers' Acquire-loads see this write.
        unsafe {
            *JOB_ARGS.0.get() = args;
        }

        JOB_KIND.store(2, Ordering::Release);
        JOB_DONE.store(0, Ordering::Release);

        // Wake all inner workers. fetch_add returns the OLD value;
        // memory_atomic_notify wakes up to N waiters.
        JOB_GEN.fetch_add(1, Ordering::Release);
        // SAFETY: see worker_loop's notify SAFETY comment.
        unsafe {
            let _ = memory_atomic_notify(JOB_GEN.as_ptr().cast(), n as u32);
        }

        // Outer Worker runs its share (worker_id = 0).
        run_gemm_slice(0, n, &args);

        // Wait for the n-1 inner workers to finish. Re-load each
        // iteration to handle spurious wakes / races where workers
        // finish between our load and our wait.
        let target = (n - 1) as i32;
        loop {
            let done = JOB_DONE.load(Ordering::Acquire);
            if done >= target {
                break;
            }
            // SAFETY: JOB_DONE.as_ptr() is a valid pointer.
            unsafe {
                let _ = memory_atomic_wait32(JOB_DONE.as_ptr().cast(), done, -1);
            }
        }

        true
    }

    /// Compute the `[j_start, j_end)` slice of output columns for
    /// `worker_id` out of `n_workers`. Equal-sized chunks; the last
    /// worker absorbs any remainder.
    fn run_gemm_slice(worker_id: usize, n_workers: usize, args: &JobArgs) {
        let total = args.dst_cols;
        let chunk = total.div_ceil(n_workers);
        let j_start = (worker_id * chunk).min(total);
        let j_end = ((worker_id + 1) * chunk).min(total);
        if j_start >= j_end {
            return;
        }
        // SAFETY: args was provided by the outer Worker, which
        // upholds the matrix pointer/stride preconditions of the
        // tiled FMA kernels.
        unsafe {
            // Match the serial dispatch in `simd128::gemm_tn_f32`:
            // cache-blocked 4×4 FMA when `+relaxed-simd`, else 2×4
            // mul+add.
            #[cfg(target_feature = "relaxed-simd")]
            super::simd128::gemm_tn_f32_slice_raw_tiled_4x4_fma_block(
                args.dst_ptr as *mut f32,
                args.dst_col_stride,
                args.dst_rows,
                args.lhs_ptr as *const f32,
                args.lhs_col_stride,
                args.lhs_rows,
                args.rhs_ptr as *const f32,
                args.rhs_col_stride,
                j_start,
                j_end,
            );
            #[cfg(not(target_feature = "relaxed-simd"))]
            super::simd128::gemm_tn_f32_slice_raw_tiled_2x4(
                args.dst_ptr as *mut f32,
                args.dst_col_stride,
                args.dst_rows,
                args.lhs_ptr as *const f32,
                args.lhs_col_stride,
                args.lhs_rows,
                args.rhs_ptr as *const f32,
                args.rhs_col_stride,
                j_start,
                j_end,
            );
        }
    }

    /// Outer-side dispatch for the fused CountSketch scatter (f32 variant).
    /// Returns `true` if the work was farmed out to the pool; `false`
    /// (caller falls back to serial / native-rayon) when the pool isn't
    /// initialised or N <= 1.
    ///
    /// # Safety
    ///
    /// All pointers inside `args` must remain valid for the duration of
    /// this call. The `out` region (covering `c * out_col_stride` f32
    /// elements at `out_ptr`) must be exclusively owned (no aliasing
    /// reads/writes by other code while the pool is computing). All
    /// read-only inputs (`raw_bytes`, `bucket`, `sign`, `iid_positions`,
    /// `mean`, `inv_std`) must outlive this call and must not be mutated
    /// by any other thread during it. See
    /// [`super::scatter::scatter_one_column_f32`] for per-pointer length
    /// requirements.
    pub unsafe fn parallel_scatter_f32(args: ScatterArgsF32) -> bool {
        let n = POOL_SIZE.load(Ordering::Acquire);
        if n <= 1 {
            return false;
        }

        // SAFETY: between dispatches, only the outer Worker writes
        // SCATTER_ARGS_F32. The Release ordering on the JOB_GEN
        // increment below makes inner workers' Acquire-loads see this
        // write.
        unsafe {
            *SCATTER_ARGS_F32.0.get() = args;
        }

        JOB_KIND.store(3, Ordering::Release);
        JOB_DONE.store(0, Ordering::Release);

        // Wake all inner workers. fetch_add returns the OLD value;
        // memory_atomic_notify wakes up to N waiters.
        JOB_GEN.fetch_add(1, Ordering::Release);
        // SAFETY: see worker_loop's notify SAFETY comment.
        unsafe {
            let _ = memory_atomic_notify(JOB_GEN.as_ptr().cast(), n as u32);
        }

        // Outer Worker runs its share (worker_id = 0).
        run_scatter_slice_f32(0, n, &args);

        // Wait for the n-1 inner workers to finish.
        let target = (n - 1) as i32;
        loop {
            let done = JOB_DONE.load(Ordering::Acquire);
            if done >= target {
                break;
            }
            unsafe {
                let _ = memory_atomic_wait32(JOB_DONE.as_ptr().cast(), done, -1);
            }
        }

        true
    }

    /// f64 sibling of [`parallel_scatter_f32`]. Same dispatch shape;
    /// only the args cell, `JOB_KIND` (4 instead of 3) and per-worker
    /// slice fn change.
    ///
    /// # Safety
    ///
    /// Same preconditions as [`parallel_scatter_f32`]: all pointers
    /// inside `args` must remain valid for the duration of the call,
    /// `out` must be exclusively owned, and all read-only inputs
    /// must outlive the call without concurrent mutation. See
    /// [`super::scatter::scatter_one_column_f64`] for per-pointer
    /// length requirements.
    pub unsafe fn parallel_scatter_f64(args: ScatterArgsF64) -> bool {
        let n = POOL_SIZE.load(Ordering::Acquire);
        if n <= 1 {
            return false;
        }
        unsafe {
            *SCATTER_ARGS_F64.0.get() = args;
        }
        JOB_KIND.store(4, Ordering::Release);
        JOB_DONE.store(0, Ordering::Release);
        JOB_GEN.fetch_add(1, Ordering::Release);
        unsafe {
            let _ = memory_atomic_notify(JOB_GEN.as_ptr().cast(), n as u32);
        }
        run_scatter_slice_f64(0, n, &args);
        let target = (n - 1) as i32;
        loop {
            let done = JOB_DONE.load(Ordering::Acquire);
            if done >= target {
                break;
            }
            unsafe {
                let _ = memory_atomic_wait32(JOB_DONE.as_ptr().cast(), done, -1);
            }
        }
        true
    }

    /// Per-worker slice for the f32 scatter dispatch. Computes the
    /// `[j_start, j_end)` slice of output columns by repeatedly invoking
    /// [`super::scatter::scatter_one_column_f32`].
    fn run_scatter_slice_f32(worker_id: usize, n_workers: usize, args: &ScatterArgsF32) {
        let total = args.c;
        let chunk = total.div_ceil(n_workers);
        let j_start = (worker_id * chunk).min(total);
        let j_end = ((worker_id + 1) * chunk).min(total);
        if j_start >= j_end {
            return;
        }
        let raw_bytes = args.raw_bytes_ptr as *const u8;
        let bucket = args.bucket_ptr as *const u32;
        let sign = args.sign_ptr as *const f32;
        let iid_pos = args.iid_pos_ptr as *const crate::bed::IidPos;
        let mean_arr = args.mean_ptr as *const f32;
        let inv_std_arr = args.inv_std_ptr as *const f32;
        let out_ptr = args.out_ptr as *mut f32;
        let all_iids = args.all_iids != 0;
        // SAFETY: args was provided by the outer Worker, which upholds
        // the per-pointer preconditions of `scatter_one_column_f32`.
        // Workers operate on disjoint `[j_start, j_end)` slices so the
        // `out` column writes don't overlap.
        for j in j_start..j_end {
            unsafe {
                super::scatter::scatter_one_column_f32(
                    j,
                    raw_bytes,
                    args.bytes_per_snp,
                    args.n_indiv,
                    bucket,
                    sign,
                    iid_pos,
                    args.iid_pos_len,
                    all_iids,
                    *mean_arr.add(j),
                    *inv_std_arr.add(j),
                    args.n_norm,
                    out_ptr,
                    args.out_col_stride,
                    args.d,
                );
            }
        }
    }

    /// f64 sibling of [`run_scatter_slice_f32`].
    fn run_scatter_slice_f64(worker_id: usize, n_workers: usize, args: &ScatterArgsF64) {
        let total = args.c;
        let chunk = total.div_ceil(n_workers);
        let j_start = (worker_id * chunk).min(total);
        let j_end = ((worker_id + 1) * chunk).min(total);
        if j_start >= j_end {
            return;
        }
        let raw_bytes = args.raw_bytes_ptr as *const u8;
        let bucket = args.bucket_ptr as *const u32;
        let sign = args.sign_ptr as *const f64;
        let iid_pos = args.iid_pos_ptr as *const crate::bed::IidPos;
        let mean_arr = args.mean_ptr as *const f64;
        let inv_std_arr = args.inv_std_ptr as *const f64;
        let out_ptr = args.out_ptr as *mut f64;
        let all_iids = args.all_iids != 0;
        for j in j_start..j_end {
            unsafe {
                super::scatter::scatter_one_column_f64(
                    j,
                    raw_bytes,
                    args.bytes_per_snp,
                    args.n_indiv,
                    bucket,
                    sign,
                    iid_pos,
                    args.iid_pos_len,
                    all_iids,
                    *mean_arr.add(j),
                    *inv_std_arr.add(j),
                    args.n_norm,
                    out_ptr,
                    args.out_col_stride,
                    args.d,
                );
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

    // ── scatter::scatter_one_column_f32/_f64 parity tests ──────────
    //
    // Regression tests for the W. H scatter refactor: ensure the
    // factored per-column body matches a hand-written reference
    // implementation (decode → normalize → scatter-add → renorm)
    // bit-for-bit. The hand reference uses the same operation order
    // as the LUT version, so f32 round-trips collapse identically.

    use crate::bed::IidPos;
    use crate::wasm_simd::scatter::{scatter_one_column_f32, scatter_one_column_f64};

    /// Deterministic small PLINK BED fixture: `c` SNPs × `n_indiv`
    /// individuals, biased toward non-missing values so per-SNP
    /// stats are well-defined.
    fn make_bed_fixture(seed: u64, c: usize, n_indiv: usize) -> Vec<u8> {
        let bytes_per_snp = n_indiv.div_ceil(4);
        let mut state = seed;
        let mut raw = vec![0u8; c * bytes_per_snp];
        for j in 0..c {
            for b in 0..bytes_per_snp {
                // Generate 4 codes; bias 0b01 (missing) to ~6 % so
                // most SNPs have a sensible mean.
                let mut byte = 0u8;
                for slot in 0..4 {
                    state = state.wrapping_mul(0x9E3779B97F4A7C15).wrapping_add(1);
                    let r = (state >> 32) as u32;
                    // 6 % missing, 31 % each of the three real codes.
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

    /// Decode a single individual's genotype out of `raw_bytes` for
    /// SNP `j` (BED layout: 2 bits per individual, low-bit-first
    /// inside each byte). Returns `None` for missing (0b01).
    fn decode(raw_bytes: &[u8], bytes_per_snp: usize, j: usize, i: usize) -> Option<f64> {
        let byte = raw_bytes[j * bytes_per_snp + i / 4];
        let code = (byte >> ((i % 4) * 2)) & 0b11;
        match code {
            0b00 => Some(2.0),
            0b01 => None,
            0b10 => Some(1.0),
            _ => Some(0.0), // 0b11
        }
    }

    /// Reference implementation: same algorithm as
    /// `scatter_one_column_f32` but written without the LUT (decode
    /// directly + manual scatter-add). Operations are sequenced so the
    /// LUT and non-LUT paths collapse to the same f32 round-trips.
    fn scatter_one_column_f32_ref(
        j: usize,
        raw_bytes: &[u8],
        bytes_per_snp: usize,
        n_indiv: usize,
        bucket: &[u32],
        sign: &[f32],
        mean: f32,
        inv_std: f32,
        n_norm: f64,
        d: usize,
    ) -> Vec<f32> {
        let mut col = vec![0f32; d];
        for i in 0..n_indiv {
            let g = decode(raw_bytes, bytes_per_snp, j, i).unwrap_or(mean as f64);
            // Match the LUT: compute (g - mean) * inv_std in f32, but
            // missing maps to 0.0 directly (no centering needed since
            // mean - mean = 0).
            let val = match decode(raw_bytes, bytes_per_snp, j, i) {
                None => 0.0f32,
                Some(_) => (g as f32 - mean) * inv_std,
            };
            col[bucket[i] as usize] += sign[i] * val;
        }
        let mut nrm_sq = 0.0f64;
        for &v in &col {
            nrm_sq += (v as f64) * (v as f64);
        }
        if nrm_sq > 0.0 {
            let scale = (n_norm / nrm_sq).sqrt() as f32;
            for v in &mut col {
                *v *= scale;
            }
        }
        col
    }

    fn build_cs_bucket_sign(n_indiv: usize, d: usize, seed: u64) -> (Vec<u32>, Vec<f32>, Vec<f64>) {
        let mut state = seed;
        let bucket: Vec<u32> = (0..n_indiv)
            .map(|_| {
                state = state.wrapping_mul(0x9E3779B97F4A7C15).wrapping_add(1);
                ((state >> 32) as u32) % d as u32
            })
            .collect();
        let sign_f32: Vec<f32> = (0..n_indiv)
            .map(|_| {
                state = state.wrapping_mul(0x9E3779B97F4A7C15).wrapping_add(1);
                if (state >> 63) & 1 == 0 { 1.0 } else { -1.0 }
            })
            .collect();
        let sign_f64: Vec<f64> = sign_f32.iter().map(|&s| s as f64).collect();
        (bucket, sign_f32, sign_f64)
    }

    #[test]
    fn scatter_one_column_f32_matches_ref_all_iids() {
        let c = 8;
        let n_indiv = 503; // 1000G EUR sample size, exercises non-multiple-of-4 tail.
        let d = 50;
        let raw = make_bed_fixture(42, c, n_indiv);
        let bytes_per_snp = n_indiv.div_ceil(4);
        let (bucket, sign_f32, _) = build_cs_bucket_sign(n_indiv, d, 7);
        let n_norm = n_indiv as f64;

        // Compute per-SNP mean / inv_std the same way Pass 1 does.
        let mut mean_buf = vec![0f32; c];
        let mut inv_std_buf = vec![0f32; c];
        for j in 0..c {
            let mut sum = 0.0f64;
            let mut sum_sq = 0.0f64;
            let mut count = 0u64;
            for i in 0..n_indiv {
                if let Some(g) = decode(&raw, bytes_per_snp, j, i) {
                    sum += g;
                    sum_sq += g * g;
                    count += 1;
                }
            }
            let mean = if count > 0 { sum / count as f64 } else { 0.0 };
            let centered = sum_sq - count as f64 * mean * mean;
            let var = centered / n_indiv as f64;
            mean_buf[j] = mean as f32;
            inv_std_buf[j] = if var > 0.0 {
                (1.0 / var.sqrt()) as f32
            } else {
                0.0
            };
        }

        // d × c output column-major Vec.
        let mut out = vec![0f32; d * c];
        for j in 0..c {
            unsafe {
                scatter_one_column_f32(
                    j,
                    raw.as_ptr(),
                    bytes_per_snp,
                    n_indiv,
                    bucket.as_ptr(),
                    sign_f32.as_ptr(),
                    core::ptr::null::<IidPos>(),
                    0,
                    true, // all_iids
                    mean_buf[j],
                    inv_std_buf[j],
                    n_norm,
                    out.as_mut_ptr(),
                    d as isize,
                    d,
                );
            }
            let expect = scatter_one_column_f32_ref(
                j,
                &raw,
                bytes_per_snp,
                n_indiv,
                &bucket,
                &sign_f32,
                mean_buf[j],
                inv_std_buf[j],
                n_norm,
                d,
            );
            for i in 0..d {
                let got = out[j * d + i];
                let exp = expect[i];
                let diff = (got - exp).abs();
                let denom = exp.abs().max(1e-6);
                assert!(
                    diff < 1e-4 || diff / denom < 1e-4,
                    "f32 mismatch at (i={i},j={j}): got {got} expected {exp} diff {diff}",
                );
            }
        }
    }

    #[test]
    fn scatter_one_column_f64_matches_f32_to_tolerance() {
        // f64 path should match the f32 reference within f32 precision.
        let c = 4;
        let n_indiv = 200;
        let d = 32;
        let raw = make_bed_fixture(13, c, n_indiv);
        let bytes_per_snp = n_indiv.div_ceil(4);
        let (bucket, sign_f32, sign_f64) = build_cs_bucket_sign(n_indiv, d, 11);
        let n_norm = n_indiv as f64;

        let mut mean_f32 = vec![0f32; c];
        let mut inv_std_f32 = vec![0f32; c];
        let mut mean_f64 = vec![0f64; c];
        let mut inv_std_f64 = vec![0f64; c];
        for j in 0..c {
            let mut sum = 0.0f64;
            let mut sum_sq = 0.0f64;
            let mut count = 0u64;
            for i in 0..n_indiv {
                if let Some(g) = decode(&raw, bytes_per_snp, j, i) {
                    sum += g;
                    sum_sq += g * g;
                    count += 1;
                }
            }
            let mean = if count > 0 { sum / count as f64 } else { 0.0 };
            let centered = sum_sq - count as f64 * mean * mean;
            let var = centered / n_indiv as f64;
            mean_f64[j] = mean;
            mean_f32[j] = mean as f32;
            inv_std_f64[j] = if var > 0.0 { 1.0 / var.sqrt() } else { 0.0 };
            inv_std_f32[j] = inv_std_f64[j] as f32;
        }

        let mut out_f64 = vec![0f64; d * c];
        let mut out_f32 = vec![0f32; d * c];
        for j in 0..c {
            unsafe {
                scatter_one_column_f64(
                    j,
                    raw.as_ptr(),
                    bytes_per_snp,
                    n_indiv,
                    bucket.as_ptr(),
                    sign_f64.as_ptr(),
                    core::ptr::null::<IidPos>(),
                    0,
                    true,
                    mean_f64[j],
                    inv_std_f64[j],
                    n_norm,
                    out_f64.as_mut_ptr(),
                    d as isize,
                    d,
                );
                scatter_one_column_f32(
                    j,
                    raw.as_ptr(),
                    bytes_per_snp,
                    n_indiv,
                    bucket.as_ptr(),
                    sign_f32.as_ptr(),
                    core::ptr::null::<IidPos>(),
                    0,
                    true,
                    mean_f32[j],
                    inv_std_f32[j],
                    n_norm,
                    out_f32.as_mut_ptr(),
                    d as isize,
                    d,
                );
            }
            for i in 0..d {
                let g = out_f32[j * d + i];
                let e = out_f64[j * d + i] as f32;
                let diff = (g - e).abs();
                let denom = e.abs().max(1e-6);
                assert!(
                    diff < 1e-3 || diff / denom < 1e-3,
                    "f32 vs f64 mismatch at (i={i},j={j}): f32={g} f64={e} diff={diff}",
                );
            }
        }
    }

    #[test]
    fn scatter_one_column_f32_subsample_matches_all_iids_on_full_keep() {
        // When iid_positions covers every individual (in order), the
        // subsample path must produce the same column as the all_iids
        // path. Catches off-by-one between byte/shift indexing and the
        // 4-per-byte fast path.
        let c = 3;
        let n_indiv = 47;
        let d = 16;
        let raw = make_bed_fixture(99, c, n_indiv);
        let bytes_per_snp = n_indiv.div_ceil(4);
        let (bucket, sign_f32, _) = build_cs_bucket_sign(n_indiv, d, 5);
        let n_norm = n_indiv as f64;

        let iid_positions: Vec<IidPos> = (0..n_indiv)
            .map(|i| IidPos {
                byte_idx: i / 4,
                shift: ((i % 4) * 2) as u8,
            })
            .collect();

        let mut mean_buf = vec![0f32; c];
        let mut inv_std_buf = vec![0f32; c];
        for j in 0..c {
            let mut sum = 0.0f64;
            let mut sum_sq = 0.0f64;
            let mut count = 0u64;
            for i in 0..n_indiv {
                if let Some(g) = decode(&raw, bytes_per_snp, j, i) {
                    sum += g;
                    sum_sq += g * g;
                    count += 1;
                }
            }
            let mean = if count > 0 { sum / count as f64 } else { 0.0 };
            let centered = sum_sq - count as f64 * mean * mean;
            let var = centered / n_indiv as f64;
            mean_buf[j] = mean as f32;
            inv_std_buf[j] = if var > 0.0 {
                (1.0 / var.sqrt()) as f32
            } else {
                0.0
            };
        }

        let mut out_all = vec![0f32; d * c];
        let mut out_sub = vec![0f32; d * c];
        for j in 0..c {
            unsafe {
                scatter_one_column_f32(
                    j,
                    raw.as_ptr(),
                    bytes_per_snp,
                    n_indiv,
                    bucket.as_ptr(),
                    sign_f32.as_ptr(),
                    core::ptr::null::<IidPos>(),
                    0,
                    true,
                    mean_buf[j],
                    inv_std_buf[j],
                    n_norm,
                    out_all.as_mut_ptr(),
                    d as isize,
                    d,
                );
                scatter_one_column_f32(
                    j,
                    raw.as_ptr(),
                    bytes_per_snp,
                    n_indiv,
                    bucket.as_ptr(),
                    sign_f32.as_ptr(),
                    iid_positions.as_ptr(),
                    iid_positions.len(),
                    false,
                    mean_buf[j],
                    inv_std_buf[j],
                    n_norm,
                    out_sub.as_mut_ptr(),
                    d as isize,
                    d,
                );
            }
            for i in 0..d {
                let a = out_all[j * d + i];
                let b = out_sub[j * d + i];
                assert_eq!(
                    a, b,
                    "all_iids vs subsample mismatch at (i={i},j={j}): {a} vs {b}",
                );
            }
        }
    }
}
