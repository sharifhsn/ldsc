use faer::{Accum, Par};

use crate::la::{ColF, MatF, col_zeros, matmul_tn_to};

/// Jackknife SE and covariance from pseudovalue expansion.
///
/// Given the full-sample estimate `est` (p×1) and a `delete_values` matrix
/// (n_blocks × p) of leave-one-block-out estimates, compute the jackknife
/// standard errors and covariance matrix.
pub fn jackknife_se_cov(est: &ColF, delete_values: &MatF, n_blocks: usize) -> (ColF, MatF) {
    let p = est.nrows();
    let n_b = n_blocks;

    let mut pseudo = MatF::zeros(n_b, p);
    for k in 0..n_b {
        for j in 0..p {
            pseudo[(k, j)] =
                est[(j, 0)] * n_b as f64 - delete_values[(k, j)] * (n_b as f64 - 1.0);
        }
    }

    let mut mean_pv = vec![0.0f64; p];
    for j in 0..p {
        let mut sum = 0.0;
        for k in 0..n_b {
            sum += pseudo[(k, j)];
        }
        mean_pv[j] = sum / n_b as f64;
    }

    let mut centered = MatF::zeros(n_b, p);
    for k in 0..n_b {
        for j in 0..p {
            centered[(k, j)] = pseudo[(k, j)] - mean_pv[j];
        }
    }

    let mut cov = MatF::zeros(p, p);
    matmul_tn_to(
        cov.as_mut(),
        centered.as_ref(),
        centered.as_ref(),
        1.0,
        Accum::Replace,
        Par::rayon(0),
    );
    let denom = (n_b as f64 - 1.0) * n_b as f64;
    for i in 0..p {
        for j in 0..p {
            cov[(i, j)] /= denom;
        }
    }

    let mut se = col_zeros(p);
    for j in 0..p {
        se[(j, 0)] = cov[(j, j)].sqrt();
    }

    (se, cov)
}


