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
#[allow(dead_code)]
pub fn mat_zeros_f32(nrows: usize, ncols: usize) -> MatF32 {
    Mat::zeros(nrows, ncols)
}

#[inline]
#[allow(dead_code)]
pub fn mat_from_fn(nrows: usize, ncols: usize, f: impl FnMut(usize, usize) -> f64) -> MatF {
    Mat::from_fn(nrows, ncols, f)
}

#[inline]
pub fn col_zeros(nrows: usize) -> ColF {
    Mat::zeros(nrows, 1)
}

#[inline]
#[allow(dead_code)]
pub fn col_ones(nrows: usize) -> ColF {
    let mut out = Mat::zeros(nrows, 1);
    for i in 0..nrows {
        out[(i, 0)] = 1.0;
    }
    out
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
#[allow(dead_code)]
pub fn col_get(col: &ColF, i: usize) -> f64 {
    col[(i, 0)]
}

#[inline]
#[allow(dead_code)]
pub fn col_set(col: &mut ColF, i: usize, v: f64) {
    col[(i, 0)] = v;
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
#[allow(dead_code)]
pub fn col_map_in_place(col: &mut ColF, f: impl Fn(f64) -> f64) {
    let n = col.nrows();
    for i in 0..n {
        col[(i, 0)] = f(col[(i, 0)]);
    }
}

#[inline]
#[allow(dead_code)]
pub fn col_mul_in_place(col: &mut ColF, other: &ColF) {
    let n = col.nrows();
    for i in 0..n {
        col[(i, 0)] *= other[(i, 0)];
    }
}

#[inline]
#[allow(dead_code)]
pub fn col_add_in_place(col: &mut ColF, other: &ColF) {
    let n = col.nrows();
    for i in 0..n {
        col[(i, 0)] += other[(i, 0)];
    }
}

#[inline]
#[allow(dead_code)]
pub fn col_sub_in_place(col: &mut ColF, other: &ColF) {
    let n = col.nrows();
    for i in 0..n {
        col[(i, 0)] -= other[(i, 0)];
    }
}

#[inline]
#[allow(dead_code)]
pub fn mat_row_scale_in_place(mat: &mut MatF, scales: &ColF) {
    let nrows = mat.nrows();
    let ncols = mat.ncols();
    for i in 0..nrows {
        let s = scales[(i, 0)];
        for j in 0..ncols {
            mat[(i, j)] *= s;
        }
    }
}

#[inline]
#[allow(dead_code)]
pub fn mat_fill(mut mat: MatMut<'_, f64>, value: f64) {
    let nrows = mat.nrows();
    let ncols = mat.ncols();
    for i in 0..nrows {
        for j in 0..ncols {
            mat[(i, j)] = value;
        }
    }
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

#[cfg(feature = "fast-f32")]
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

#[cfg(feature = "fast-f32")]
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

#[cfg(feature = "fast-f32")]
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
