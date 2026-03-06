use faer::linalg::matmul::matmul;
use faer::{Accum, Mat, MatMut, MatRef, Par};

pub type MatF = Mat<f64>;
pub type MatF32 = Mat<f32>;
pub type ColF = Mat<f64>; // column vector as n x 1 matrix

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
    for i in 0..nrows {
        for j in 0..ncols {
            dst[(i, j)] += src[(i, j)];
        }
    }
}

#[inline]
pub fn mat_copy_from(mut dst: MatMut<'_, f64>, src: MatRef<'_, f64>) {
    let nrows = dst.nrows();
    let ncols = dst.ncols();
    for i in 0..nrows {
        for j in 0..ncols {
            dst[(i, j)] = src[(i, j)];
        }
    }
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
