/// Iteratively Re-Weighted Least Squares.
///
/// Weight matrices are pre-allocated and updated in-place across IRWLS iterations.
use anyhow::Result;
use ndarray::{Array1, Array2};
use ndarray_linalg::LeastSquaresSvd;

/// Result of an IRWLS fit, mirroring the Python IrwlsResult namedtuple.
#[derive(Debug)]
pub struct IrwlsResult {
    /// Point estimates (one per predictor column)
    pub est: Array1<f64>,
    /// Jackknife SE (filled in by jackknife.rs after regression)
    pub jknife_se: Option<Array1<f64>>,
    /// Jackknife variance (SE² × n_blocks)
    #[allow(dead_code)]
    pub jknife_var: Option<Array1<f64>>,
    /// Jackknife covariance matrix
    pub jknife_cov: Option<Array2<f64>>,
    /// Per-block delete values (n_blocks × n_params)
    pub delete_values: Option<Array2<f64>>,
}

/// Fit IRWLS regression of `y` on `x` using initial `weights`.
///
/// # Arguments
/// * `x`        — design matrix (n_obs × n_params), MUST include intercept column
/// * `y`        — response vector (n_obs,)
/// * `weights`  — initial weights (n_obs,); updated in-place each iteration
/// * `n_iter`   — number of reweighting iterations (typically 2)
pub fn irwls(
    x: &Array2<f64>,
    y: &Array1<f64>,
    weights: &mut Array1<f64>,
    n_iter: usize,
) -> Result<IrwlsResult> {
    let n = x.nrows();
    let p = x.ncols();

    // Pre-allocate weighted matrices — reused across iterations (zero-alloc loop).
    let mut xw = Array2::<f64>::zeros((n, p));
    let mut yw = Array1::<f64>::zeros(n);
    let mut est = Array1::<f64>::zeros(p);

    for _ in 0..n_iter {
        // Apply weights: xw[i,:] = x[i,:] * sqrt(w[i]), yw[i] = y[i] * sqrt(w[i])
        for i in 0..n {
            let sw = weights[i].sqrt();
            xw.row_mut(i).assign(&(&x.row(i) * sw));
            yw[i] = y[i] * sw;
        }

        // Solve weighted OLS via SVD.
        let result = xw.least_squares(&yw)?;
        est.assign(&result.solution);

        // Update weights: w[i] = 1 / fitted[i]².
        let fitted = x.dot(&est);
        for i in 0..n {
            let f = fitted[i].max(1e-9); // guard against divide-by-zero
            weights[i] = 1.0 / (f * f);
        }
    }

    Ok(IrwlsResult {
        est,
        jknife_se: None,
        jknife_var: None,
        jknife_cov: None,
        delete_values: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    /// Smoke test: IRWLS on a trivial linear system should recover slope/intercept.
    #[test]
    fn test_irwls_trivial() {
        // y = 2x + 1  (slope=2, intercept=1)
        let x = array![[1.0, 0.0], [1.0, 1.0], [1.0, 2.0], [1.0, 3.0]]; // [intercept, x]
        let y = array![1.0, 3.0, 5.0, 7.0];
        let mut w = Array1::ones(4);

        let res = irwls(&x, &y, &mut w, 2).unwrap();
        let intercept = res.est[0];
        let slope = res.est[1];

        assert!((intercept - 1.0).abs() < 1e-6, "intercept={}", intercept);
        assert!((slope - 2.0).abs() < 1e-6, "slope={}", slope);
    }

    /// 1D intercept-only system: x = 1, y = 1 → coef = 1.
    #[test]
    fn test_irwls_1d_constant() {
        let x = array![[1.0], [1.0], [1.0], [1.0]];
        let y = array![1.0, 1.0, 1.0, 1.0];
        let mut w = Array1::ones(4);
        let res = irwls(&x, &y, &mut w, 2).unwrap();
        assert!((res.est[0] - 1.0).abs() < 1e-9, "coef={}", res.est[0]);
    }

    /// 2D system where y = col1 + col2 → coefficients should be [1, 1].
    #[test]
    fn test_irwls_2d_sum() {
        // x = [[1,1],[1,4],[1,3],[1,2]], y = row sums = [2,5,4,3]
        let x = array![[1.0, 1.0], [1.0, 4.0], [1.0, 3.0], [1.0, 2.0]];
        let y = array![2.0, 5.0, 4.0, 3.0];
        let mut w = Array1::ones(4);
        let res = irwls(&x, &y, &mut w, 2).unwrap();
        assert!((res.est[0] - 1.0).abs() < 1e-6, "coef0={}", res.est[0]);
        assert!((res.est[1] - 1.0).abs() < 1e-6, "coef1={}", res.est[1]);
    }

    /// IRWLS with non-uniform initial weights: zero-weight rows contribute sqrt(0)=0
    /// to xw/yw and are effectively excluded from the fit.
    #[test]
    fn test_irwls_nonuniform_weights() {
        let x = array![[1.0], [1.0], [1.0], [1.0]];
        let y = array![1.0, 1.0, 1.0, 1.0];
        let mut w = array![1.0, 0.0, 0.0, 1.0];
        let res = irwls(&x, &y, &mut w, 2).unwrap();
        assert!((res.est[0] - 1.0).abs() < 1e-6, "coef={}", res.est[0]);
    }
}
