use faer::linalg::matmul::matmul;
use faer::{Accum, Mat, MatMut, MatRef, Par};

pub type MatF = Mat<f64>;
pub type MatF32 = Mat<f32>;
pub type ColF = Mat<f64>; // column vector as n x 1 matrix

/// Default `Par` setting for compute-heavy faer kernels.
///
/// On native targets this returns `Par::rayon(0)` (use the global rayon
/// pool with all cores). On `wasm32-unknown-unknown` this returns
/// `Par::Seq` — the WASM build disables faer's `rayon` feature
/// (`spindle` / `atomic-wait` don't compile to wasm), so any
/// `Par::rayon(_)` value would panic at runtime. The caller-facing API
/// is identical so the GEMM sites don't need their own cfg gates.
#[inline]
pub fn par_default() -> Par {
    #[cfg(not(target_arch = "wasm32"))]
    {
        Par::rayon(0)
    }
    #[cfg(target_arch = "wasm32")]
    {
        Par::Seq
    }
}

#[inline]
pub fn mat_zeros(nrows: usize, ncols: usize) -> MatF {
    Mat::zeros(nrows, ncols)
}

#[inline]
pub fn col_zeros(nrows: usize) -> ColF {
    Mat::zeros(nrows, 1)
}

#[inline]
pub fn col_from_vec(values: Vec<f64>) -> ColF {
    let n = values.len();
    let mut out = Mat::zeros(n, 1);
    for (i, v) in values.into_iter().enumerate() {
        out[(i, 0)] = v;
    }
    out
}

#[inline]
pub fn col_len(col: &ColF) -> usize {
    col.nrows()
}

#[inline]
pub fn col_sum(col: &ColF) -> f64 {
    let n = col.nrows();
    let mut sum = 0.0;
    for i in 0..n {
        sum += col[(i, 0)];
    }
    sum
}

#[inline]
pub fn col_mean(col: &ColF) -> f64 {
    let n = col.nrows();
    if n == 0 {
        return f64::NAN;
    }
    col_sum(col) / n as f64
}

#[inline]
pub fn mat_add_in_place(mut dst: MatMut<'_, f64>, src: MatRef<'_, f64>) {
    let nrows = dst.nrows();
    let ncols = dst.ncols();
    for j in 0..ncols {
        for i in 0..nrows {
            dst[(i, j)] += src[(i, j)];
        }
    }
}

#[inline]
pub fn mat_copy_from(mut dst: MatMut<'_, f64>, src: MatRef<'_, f64>) {
    dst.copy_from(src);
}

#[inline]
pub fn matmul_to(
    dst: MatMut<'_, f64>,
    lhs: MatRef<'_, f64>,
    rhs: MatRef<'_, f64>,
    alpha: f64,
    accum: Accum,
    par: Par,
) {
    matmul(dst, accum, lhs, rhs, alpha, par);
}

#[inline]
pub fn matmul_tn_to(
    dst: MatMut<'_, f64>,
    lhs: MatRef<'_, f64>,
    rhs: MatRef<'_, f64>,
    alpha: f64,
    accum: Accum,
    par: Par,
) {
    matmul(dst, accum, lhs.transpose(), rhs, alpha, par);
}

#[inline]
pub fn matmul_tn_to_f32(
    dst: MatMut<'_, f32>,
    lhs: MatRef<'_, f32>,
    rhs: MatRef<'_, f32>,
    alpha: f32,
    accum: Accum,
    par: Par,
) {
    // On wasm32 with simd128, route the f32 hot-path TN matmul into
    // our hand-rolled kernel — faer's `pulp` SIMD backend doesn't
    // support wasm32 and otherwise lowers to a fully scalar GEMM
    // (the entire ~30-50× wall-time gap vs native CLI on full 1000G).
    //
    // Only the alpha=1.0 / accum=Replace fast path is custom; both
    // hot-path call sites in `l2/compute.rs::compute_ldscore_global`
    // (within-chunk B^T·B at compute.rs:1718, cross-window A^T·B at
    // compute.rs:2027) use exactly these parameters. Other shapes
    // (rare on wasm) fall through to faer.
    //
    // F.2 dispatches to the v128 SIMD microkernel
    // (`simd128::gemm_tn_f32`): 4 parallel f32x4 accumulators,
    // K-unroll-16, scalar tail. Bit-identical to the scalar reference
    // (and to the previous faer scalar fallback at L2-score precision)
    // but ~10× faster — closing most of the remaining gap to native.
    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
    {
        // Comparison to the exact f32 literal `1.0` is intentional:
        // both compute.rs call sites pass `1.0f32` literally, so the
        // bit pattern matches. Anything else falls to faer.
        #[allow(clippy::float_cmp)]
        let is_unit_alpha = alpha == 1.0;
        // Also require column-major (row_stride==1) on all three
        // operands. The kernel debug-asserts this; the explicit
        // pre-check here means non-contiguous shapes (rare in the
        // wild — we don't construct them in compute.rs) silently fall
        // through to faer instead of UB-ing in release builds.
        let all_col_major = dst.row_stride() == 1 && lhs.row_stride() == 1 && rhs.row_stride() == 1;
        if is_unit_alpha && matches!(accum, Accum::Replace) && all_col_major {
            // SAFETY: the four preconditions of `simd128::gemm_tn_f32`
            // are all upheld here:
            //   - col-major (just checked)
            //   - shapes consistent (faer matmul has the same
            //     requirement; the kernel debug-asserts identically)
            //   - col strides non-negative (faer Mat / submatrix
            //     invariant)
            unsafe {
                crate::wasm_simd::simd128::gemm_tn_f32(dst, lhs, rhs);
            }
            // `par` is intentionally unused on this path: our kernel
            // is single-threaded (faer's wasm fallback is too — the
            // `rayon` feature is disabled on wasm because spindle and
            // atomic-wait don't compile to wasm32-unknown-unknown).
            let _ = par;
            return;
        }
    }
    matmul(dst, accum, lhs.transpose(), rhs, alpha, par);
}

#[inline]
pub fn matmul_nt_to(
    dst: MatMut<'_, f64>,
    lhs: MatRef<'_, f64>,
    rhs: MatRef<'_, f64>,
    alpha: f64,
    accum: Accum,
    par: Par,
) {
    matmul(dst, accum, lhs, rhs.transpose(), alpha, par);
}

#[inline]
pub fn mat_slice<'a>(
    mat: MatRef<'a, f64>,
    rows: std::ops::Range<usize>,
    cols: std::ops::Range<usize>,
) -> MatRef<'a, f64> {
    let nrows = rows.end.saturating_sub(rows.start);
    let ncols = cols.end.saturating_sub(cols.start);
    mat.submatrix(rows.start, cols.start, nrows, ncols)
}

#[inline]
pub fn mat_slice_f32<'a>(
    mat: MatRef<'a, f32>,
    rows: std::ops::Range<usize>,
    cols: std::ops::Range<usize>,
) -> MatRef<'a, f32> {
    let nrows = rows.end.saturating_sub(rows.start);
    let ncols = cols.end.saturating_sub(cols.start);
    mat.submatrix(rows.start, cols.start, nrows, ncols)
}

#[inline]
pub fn mat_slice_mut<'a>(
    mat: MatMut<'a, f64>,
    rows: std::ops::Range<usize>,
    cols: std::ops::Range<usize>,
) -> MatMut<'a, f64> {
    let nrows = rows.end.saturating_sub(rows.start);
    let ncols = cols.end.saturating_sub(cols.start);
    mat.submatrix_mut(rows.start, cols.start, nrows, ncols)
}

#[inline]
pub fn mat_slice_mut_f32<'a>(
    mat: MatMut<'a, f32>,
    rows: std::ops::Range<usize>,
    cols: std::ops::Range<usize>,
) -> MatMut<'a, f32> {
    let nrows = rows.end.saturating_sub(rows.start);
    let ncols = cols.end.saturating_sub(cols.start);
    mat.submatrix_mut(rows.start, cols.start, nrows, ncols)
}
