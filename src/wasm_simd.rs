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
    #[cfg(target_feature = "relaxed-simd")]
    use core::arch::wasm32::f32x4_relaxed_madd;
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
                    dst_ptr,
                    dst_col,
                    m,
                    n,
                    a_ptr,
                    a_col,
                    k,
                    b_ptr,
                    b_col,
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
                    let mut s0 = (f32x4_extract_lane::<0>(acc0)
                        + f32x4_extract_lane::<1>(acc0))
                        + (f32x4_extract_lane::<2>(acc0) + f32x4_extract_lane::<3>(acc0));
                    let mut s1 = (f32x4_extract_lane::<0>(acc1)
                        + f32x4_extract_lane::<1>(acc1))
                        + (f32x4_extract_lane::<2>(acc1) + f32x4_extract_lane::<3>(acc1));
                    let mut s2 = (f32x4_extract_lane::<0>(acc2)
                        + f32x4_extract_lane::<1>(acc2))
                        + (f32x4_extract_lane::<2>(acc2) + f32x4_extract_lane::<3>(acc2));
                    let mut s3 = (f32x4_extract_lane::<0>(acc3)
                        + f32x4_extract_lane::<1>(acc3))
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

    /// Job discriminator. `1` = exit (worker_loop returns), `2` = gemm.
    pub static JOB_KIND: AtomicI32 = AtomicI32::new(0);

    /// Wrapper around `UnsafeCell<JobArgs>` to manually mark Sync. Sound
    /// because all access is ordered through `JOB_GEN`'s
    /// release/acquire fences.
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
