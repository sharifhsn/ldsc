/// Block jackknife variance estimation using parallel IRWLS folds.
use anyhow::Result;
use faer::{Accum, Par};
use rayon::prelude::*;

use crate::irwls::{IrwlsResult, irwls};
use crate::la::{ColF, MatF, col_len, col_zeros, matmul_tn_to};

/// Run block jackknife over `n_blocks` contiguous genomic windows.
///
/// For each block `k`, the estimator is refit on all data *except* block `k`.
/// Pseudo-values are then combined to produce a jackknife SE.
///
/// # Arguments
/// * `x`        — design matrix (n_obs × n_params)
/// * `y`        — response vector (n_obs,)
/// * `weights`  — regression weights (n_obs,)
/// * `n_blocks` — number of jackknife blocks (typically ~200 for LDSC)
/// * `n_iter`   — IRWLS iterations per fold
#[allow(dead_code)]
pub fn jackknife(
    x: &MatF,
    y: &ColF,
    weights: &ColF,
    n_blocks: usize,
    n_iter: usize,
) -> Result<IrwlsResult> {
    let n = x.nrows();
    let p = x.ncols();
    let block_size = n / n_blocks;

    // Fit the full-data estimator first.
    let mut w_full = weights.clone();
    let full_fit = irwls(x, y, &mut w_full, n_iter)?;
    let n_params = col_len(&full_fit.est);

    // Parallel jackknife: each block independently removed, IRWLS refit.
    let delete_values: Vec<ColF> = (0..n_blocks)
        .into_par_iter()
        .map(|k| {
            let lo = k * block_size;
            let hi = ((k + 1) * block_size).min(n);
            let keep = n - (hi - lo);
            let mut x_jack = MatF::zeros(keep, p);
            let mut y_jack = col_zeros(keep);
            let mut w_jack = col_zeros(keep);
            let mut row = 0usize;
            for i in 0..n {
                if i < lo || i >= hi {
                    for j in 0..p {
                        x_jack[(row, j)] = x[(i, j)];
                    }
                    y_jack[(row, 0)] = y[(i, 0)];
                    w_jack[(row, 0)] = weights[(i, 0)];
                    row += 1;
                }
            }

            irwls(&x_jack, &y_jack, &mut w_jack, n_iter)
                .expect("irwls on jackknife fold")
                .est
        })
        .collect();

    // Stack delete values into (n_blocks × n_params).
    let n_b = delete_values.len();
    let mut delete_mat = MatF::zeros(n_b, n_params);
    for (k, dv) in delete_values.iter().enumerate() {
        for j in 0..n_params {
            delete_mat[(k, j)] = dv[(j, 0)];
        }
    }

    // Pseudo-values: theta_full + (n_blocks - 1) * (theta_full - theta_k)
    // Jackknife SE = std(pseudo-values) / sqrt(n_blocks)
    let theta = &full_fit.est;
    let mut pseudo = MatF::zeros(n_b, n_params);
    for k in 0..n_b {
        for j in 0..n_params {
            let pv = theta[(j, 0)] + (theta[(j, 0)] - delete_mat[(k, j)]) * (n_b as f64 - 1.0);
            pseudo[(k, j)] = pv;
        }
    }

    let mut mean_pv = vec![0.0f64; n_params];
    for j in 0..n_params {
        let mut sum = 0.0;
        for k in 0..n_b {
            sum += pseudo[(k, j)];
        }
        mean_pv[j] = sum / n_b as f64;
    }

    // Per-column sample variance: var_pv[j] = Σ_k (pseudo[k,j] − mean_pv[j])² / (n_b − 1)
    let mut var_pv = col_zeros(n_params);
    for j in 0..n_params {
        let m = mean_pv[j];
        let mut sum = 0.0;
        for k in 0..n_b {
            sum += (pseudo[(k, j)] - m).powi(2);
        }
        var_pv[(j, 0)] = sum / (n_b as f64 - 1.0);
    }

    // Jackknife SE = sqrt(sample_variance / n_blocks)  [Bulik-Sullivan 2015 eq. S4]
    let mut se = col_zeros(n_params);
    for j in 0..n_params {
        se[(j, 0)] = (var_pv[(j, 0)] / n_b as f64).sqrt();
    }

    // Covariance matrix of the jackknife estimator via BLAS DGEMM.
    // cov = centered.T @ centered / ((n_b - 1) * n_b)
    // where centered[k, j] = pseudo[k, j] - mean_pv[j]
    let mut centered = MatF::zeros(n_b, n_params);
    for k in 0..n_b {
        for j in 0..n_params {
            centered[(k, j)] = pseudo[(k, j)] - mean_pv[j];
        }
    }
    let mut jknife_cov = MatF::zeros(n_params, n_params);
    matmul_tn_to(
        jknife_cov.as_mut(),
        centered.as_ref(),
        centered.as_ref(),
        1.0,
        Accum::Replace,
        Par::rayon(0),
    );
    let denom = (n_b as f64 - 1.0) * n_b as f64;
    for i in 0..n_params {
        for j in 0..n_params {
            jknife_cov[(i, j)] /= denom;
        }
    }

    Ok(IrwlsResult {
        est: full_fit.est,
        jknife_se: Some(se),
        jknife_var: Some(var_pv),
        jknife_cov: Some(jknife_cov),
        delete_values: Some(delete_mat),
    })
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::la::col_ones;

    /// Jackknife on a perfectly linear system should produce near-zero SE
    /// (every block gives the same estimate, so variance of pseudo-values ≈ 0).
    #[test]
    fn test_jackknife_trivial_se() {
        // y = 2x + 1  (no noise) across 100 observations
        let n = 100usize;
        let mut x = MatF::zeros(n, 2);
        for i in 0..n {
            x[(i, 0)] = i as f64;
            x[(i, 1)] = 1.0;
        }
        let mut y = col_zeros(n);
        for i in 0..n {
            y[(i, 0)] = 2.0 * i as f64 + 1.0;
        }
        let w = col_ones(n);

        let res = jackknife(&x, &y, &w, 10, 2).unwrap();

        // Point estimates should recover slope=2, intercept=1.
        assert!(
            (res.est[(0, 0)] - 2.0).abs() < 1e-6,
            "slope={}",
            res.est[(0, 0)]
        );
        assert!(
            (res.est[(1, 0)] - 1.0).abs() < 1e-6,
            "intercept={}",
            res.est[(1, 0)]
        );

        // SEs should be tiny for a noiseless system.
        let se = res.jknife_se.unwrap();
        assert!(se[(0, 0)] < 1e-6, "SE(slope) too large: {}", se[(0, 0)]);
        assert!(se[(1, 0)] < 1e-6, "SE(intercept) too large: {}", se[(1, 0)]);
    }

    /// Jackknife SE for each parameter should be computed independently
    /// (tests the per-column variance fix).
    #[test]
    fn test_jackknife_se_per_column() {
        // Two-parameter system where the two params have very different magnitudes.
        // slope ≈ 1000, intercept ≈ 0.001 — ensures the column-mixing bug would
        // produce obviously wrong SEs for one of them.
        let n = 80usize;
        let mut x = MatF::zeros(n, 2);
        for i in 0..n {
            x[(i, 0)] = i as f64 * 1000.0;
            x[(i, 1)] = 1.0;
        }
        let mut y = col_zeros(n);
        for i in 0..n {
            y[(i, 0)] = i as f64 * 1000.0 + 0.001;
        }
        let w = col_ones(n);

        let res = jackknife(&x, &y, &w, 8, 2).unwrap();
        let se = res.jknife_se.unwrap();

        // Both SEs should be very small for a noiseless system,
        // regardless of the parameter magnitudes.
        assert!(se[(0, 0)] < 1e-6, "SE(slope)={} too large", se[(0, 0)]);
        assert!(se[(1, 0)] < 1e-2, "SE(intercept)={} too large", se[(1, 0)]);
    }

    /// delete_values matrix should have shape (n_blocks, n_params).
    #[test]
    fn test_jackknife_delete_vals_shape() {
        let n = 20usize;
        let n_blocks = 5;
        let mut x = MatF::zeros(n, 2);
        for i in 0..n {
            x[(i, 0)] = i as f64;
            x[(i, 1)] = 1.0;
        }
        let mut y = col_zeros(n);
        for i in 0..n {
            y[(i, 0)] = 2.0 * i as f64 + 1.0;
        }
        let w = col_ones(n);

        let res = jackknife(&x, &y, &w, n_blocks, 2).unwrap();
        let dv = res.delete_values.unwrap();
        assert_eq!(dv.nrows(), n_blocks);
        assert_eq!(dv.ncols(), 2);
    }

    /// Covariance matrix should be (n_params × n_params) and symmetric with
    /// non-negative diagonal entries.
    #[test]
    fn test_jackknife_cov_symmetry() {
        let n = 20usize;
        let mut x = MatF::zeros(n, 2);
        for i in 0..n {
            x[(i, 0)] = i as f64;
            x[(i, 1)] = 1.0;
        }
        let mut y = col_zeros(n);
        for i in 0..n {
            y[(i, 0)] = 2.0 * i as f64 + 1.0;
        }
        let w = col_ones(n);

        let res = jackknife(&x, &y, &w, 5, 2).unwrap();
        let cov = res.jknife_cov.unwrap();
        assert_eq!(cov.nrows(), 2);
        assert_eq!(cov.ncols(), 2);
        assert!(cov[(0, 0)] >= 0.0, "cov[0,0] should be ≥ 0");
        assert!(cov[(1, 1)] >= 0.0, "cov[1,1] should be ≥ 0");
        assert!(
            (cov[(0, 1)] - cov[(1, 0)]).abs() < 1e-12,
            "cov not symmetric"
        );
    }

    /// Known SE: 10 obs, 1 param (intercept), y = 0..9, n_blocks = 10.
    ///
    /// Leave-one-out estimates: theta_{-k} = (45 - k) / 9.
    /// Pseudo-values: pv_k = theta + 9*(theta - theta_{-k}) = k.
    /// Sample variance of [0..9] = 82.5/9; SE = sqrt(82.5/9/10) ≈ 0.9574.
    #[test]
    fn test_jackknife_known_se() {
        let n = 10usize;
        let mut x = MatF::zeros(n, 1);
        for i in 0..n {
            x[(i, 0)] = 1.0;
        }
        let mut y = col_zeros(n);
        for i in 0..n {
            y[(i, 0)] = i as f64;
        }
        let w = col_ones(n);

        let res = jackknife(&x, &y, &w, n, 2).unwrap();

        assert!(
            (res.est[(0, 0)] - 4.5).abs() < 1e-6,
            "est={}",
            res.est[(0, 0)]
        );

        let se = res.jknife_se.unwrap();
        assert!((se[(0, 0)] - 0.95742).abs() < 0.001, "se={}", se[(0, 0)]);

        // Python jknife_var = SE² ≈ 0.9167; in Rust jknife_var is sample variance.
        // Verify SE² = jknife_var / n_blocks ≈ 0.9167.
        let var = res.jknife_var.unwrap();
        let se2 = var[(0, 0)] / n as f64;
        assert!((se2 - 0.91667).abs() < 0.001, "SE²={}", se2);
    }

    /// Verifies the leave-one-block-out delete_values are computed with the correct formula.
    /// With n=10, n_blocks=10 (one obs per block), x=ones, y=[0..9]:
    ///   delete_values[k] = mean(y excluding obs k) = (45 - k) / 9.
    #[test]
    fn test_jackknife_delete_values_correctness() {
        let n = 10usize;
        let mut x = MatF::zeros(n, 1);
        for i in 0..n {
            x[(i, 0)] = 1.0;
        }
        let mut y = col_zeros(n);
        for i in 0..n {
            y[(i, 0)] = i as f64;
        }
        let w = col_ones(n);

        let res = jackknife(&x, &y, &w, n, 2).unwrap();
        let dv = res.delete_values.unwrap();

        assert_eq!(dv.nrows(), n);
        assert_eq!(dv.ncols(), 1);
        for k in 0..n {
            let expected = (45.0 - k as f64) / 9.0;
            assert!(
                (dv[(k, 0)] - expected).abs() < 1e-9,
                "delete_values[{k}]={:.9} expected {expected:.9}",
                dv[(k, 0)]
            );
        }
    }

    /// Verifies the pseudo-value formula: pv_k = theta + (n_b-1)*(theta - delete_k).
    /// For x=ones, y=[0..9], n_blocks=10: pseudo[k] = k exactly.
    #[test]
    fn test_jackknife_pseudovalue_formula() {
        let n = 10usize;
        let mut x = MatF::zeros(n, 1);
        for i in 0..n {
            x[(i, 0)] = 1.0;
        }
        let mut y = col_zeros(n);
        for i in 0..n {
            y[(i, 0)] = i as f64;
        }
        let w = col_ones(n);

        let res = jackknife(&x, &y, &w, n, 2).unwrap();
        let theta = res.est[(0, 0)]; // 4.5
        let dv = res.delete_values.as_ref().unwrap();

        // Reconstruct pseudo-values and verify they equal [0, 1, ..., 9]
        for k in 0..n {
            let pv = theta + (n as f64 - 1.0) * (theta - dv[(k, 0)]);
            assert!(
                (pv - k as f64).abs() < 1e-9,
                "pseudo[{k}]={pv:.9} expected {k}"
            );
        }
    }
}
