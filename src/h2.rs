//! LDSC h2 regression logic aligned to the Python implementation.
use crate::jackknife::jackknife_se_cov;
use crate::la::{
    ColF, MatF, col_from_vec, col_len, col_mean, col_sum, col_zeros, mat_add_in_place, mat_slice,
    matmul_tn_to,
};
use anyhow::{Result, bail};
use faer::linalg::solvers::{PartialPivLu, Solve, SolveLstsq, Svd};
use faer::{Accum, Par};

#[derive(Debug, Clone)]
pub struct JackknifeResult {
    pub est: ColF,
    pub jknife_se: ColF,
    pub jknife_cov: MatF,
    pub delete_values: MatF,
    pub separators: Vec<usize>,
}

#[derive(Debug, Clone)]
pub struct H2Result {
    pub h2: f64,
    pub h2_se: f64,
    pub intercept: f64,
    pub intercept_se: f64,
    pub mean_chi2: f64,
    pub lambda_gc: f64,
    pub ratio: Option<(f64, f64)>,
}

fn get_separators(n: usize, n_blocks: usize) -> Vec<usize> {
    let mut seps = Vec::with_capacity(n_blocks + 1);
    for i in 0..=n_blocks {
        let v = ((i as f64) * (n as f64) / (n_blocks as f64)).floor() as usize;
        seps.push(v);
    }
    if let Some(last) = seps.last_mut() {
        *last = n;
    }
    seps
}

pub(crate) fn update_separators(seps: &[usize], mask: &[bool]) -> Vec<usize> {
    let maplist: Vec<usize> = mask
        .iter()
        .enumerate()
        .filter_map(|(i, keep)| if *keep { Some(i) } else { None })
        .collect();
    let mut out = Vec::with_capacity(seps.len());
    out.push(0);
    if seps.len() > 2 {
        for &s in &seps[1..seps.len() - 1] {
            out.push(maplist[s]);
        }
    }
    out.push(mask.len());
    out
}

pub(crate) fn weight_xy(x: &MatF, y: &ColF, w: &ColF) -> Result<(MatF, ColF)> {
    let n = x.nrows();
    let p = x.ncols();
    if col_len(y) != n || col_len(w) != n {
        bail!("weight_xy: shape mismatch");
    }
    for i in 0..n {
        if w[(i, 0)] <= 0.0 {
            bail!("weights must be > 0");
        }
    }
    let sum_w = col_sum(w);
    let mut xw = MatF::zeros(n, p);
    let mut yw = col_zeros(n);
    for i in 0..n {
        let wi = w[(i, 0)] / sum_w;
        for j in 0..p {
            xw[(i, j)] = x[(i, j)] * wi;
        }
        yw[(i, 0)] = y[(i, 0)] * wi;
    }
    Ok((xw, yw))
}

fn wls(x: &MatF, y: &ColF, w: &ColF) -> Result<ColF> {
    let (xw, yw) = weight_xy(x, y, w)?;
    let svd = Svd::new(xw.as_ref()).map_err(|err| anyhow::anyhow!("svd for wls: {:?}", err))?;
    let mut rhs = yw.clone();
    svd.solve_lstsq_in_place(rhs.as_mut());
    Ok(rhs)
}

fn solve_square(xtx: &MatF, xty: &ColF) -> Result<ColF> {
    let lu = PartialPivLu::new(xtx.as_ref());
    let mut rhs = xty.clone();
    lu.solve_in_place(rhs.as_mut());
    Ok(rhs)
}

pub(crate) fn jackknife_fast(
    x: &MatF,
    y: &ColF,
    n_blocks: usize,
    separators: Option<&[usize]>,
) -> Result<JackknifeResult> {
    let n = x.nrows();
    let p = x.ncols();
    if col_len(y) != n {
        bail!("jackknife_fast: shape mismatch");
    }
    let seps = if let Some(s) = separators {
        if s.len() < 2 || *s.last().unwrap() != n {
            bail!("invalid separators");
        }
        s.to_vec()
    } else {
        get_separators(n, n_blocks)
    };
    let n_blocks = seps.len() - 1;
    if n_blocks > n {
        bail!("More blocks than data points.");
    }

    let mut xty_blocks = MatF::zeros(n_blocks, p);
    let mut xtx_blocks: Vec<MatF> = Vec::with_capacity(n_blocks);
    for i in 0..n_blocks {
        let lo = seps[i];
        let hi = seps[i + 1];
        let xb = mat_slice(x.as_ref(), lo..hi, 0..p);
        let yb = mat_slice(y.as_ref(), lo..hi, 0..1);

        let mut xty = MatF::zeros(p, 1);
        matmul_tn_to(xty.as_mut(), xb, yb, 1.0, Accum::Replace, Par::rayon(0));
        for j in 0..p {
            xty_blocks[(i, j)] = xty[(j, 0)];
        }

        let mut xtx = MatF::zeros(p, p);
        matmul_tn_to(xtx.as_mut(), xb, xb, 1.0, Accum::Replace, Par::rayon(0));
        xtx_blocks.push(xtx);
    }

    let mut xty_tot = col_zeros(p);
    for i in 0..n_blocks {
        for j in 0..p {
            xty_tot[(j, 0)] += xty_blocks[(i, j)];
        }
    }

    let mut xtx_tot = MatF::zeros(p, p);
    for xtx in &xtx_blocks {
        mat_add_in_place(xtx_tot.as_mut(), xtx.as_ref());
    }
    let est = solve_square(&xtx_tot, &xty_tot)?;

    let mut delete_values = MatF::zeros(n_blocks, p);
    for (i, xtx_block) in xtx_blocks.iter().enumerate() {
        let mut xty_del = col_zeros(p);
        for j in 0..p {
            xty_del[(j, 0)] = xty_tot[(j, 0)] - xty_blocks[(i, j)];
        }
        let mut xtx_del = MatF::zeros(p, p);
        for r in 0..p {
            for c in 0..p {
                xtx_del[(r, c)] = xtx_tot[(r, c)] - xtx_block[(r, c)];
            }
        }
        let coef = solve_square(&xtx_del, &xty_del)?;
        for j in 0..p {
            delete_values[(i, j)] = coef[(j, 0)];
        }
    }

    let (jknife_se, jknife_cov) = jackknife_se_cov(&est, &delete_values, n_blocks);

    Ok(JackknifeResult {
        est,
        jknife_se,
        jknife_cov,
        delete_values,
        separators: seps,
    })
}

pub(crate) fn ldsc_weights(
    ld: &ColF,
    w_ld: &ColF,
    n_vec: &ColF,
    m: f64,
    hsq: f64,
    intercept: f64,
) -> ColF {
    let mut hsq = hsq;
    hsq = hsq.clamp(0.0, 1.0);
    let n = col_len(ld);
    let mut out = col_zeros(n);
    for i in 0..n {
        let ld_i = ld[(i, 0)].max(1.0);
        let w_i = w_ld[(i, 0)].max(1.0);
        let c = hsq * n_vec[(i, 0)] / m;
        let het_w = 1.0 / (2.0 * (intercept + c * ld_i).powi(2));
        let oc_w = 1.0 / w_i;
        out[(i, 0)] = het_w * oc_w;
    }
    out
}

pub(crate) fn aggregate(y: &ColF, x: &ColF, n_vec: &ColF, m: f64, intercept: f64) -> f64 {
    let num = m * (col_mean(y) - intercept);
    let n = col_len(x);
    let mut denom_vec = col_zeros(n);
    for i in 0..n {
        denom_vec[(i, 0)] = x[(i, 0)] * n_vec[(i, 0)];
    }
    let denom = col_mean(&denom_vec);
    num / denom
}

pub(crate) fn combine_twostep(
    step1: &JackknifeResult,
    step2: &JackknifeResult,
    c: f64,
) -> Result<JackknifeResult> {
    let n_blocks = step1.delete_values.nrows();
    let n_annot = step2.delete_values.ncols();
    let step1_int = step1.est[(n_annot, 0)];

    let mut est = col_zeros(n_annot + 1);
    for j in 0..n_annot {
        est[(j, 0)] = step2.est[(j, 0)];
    }
    est[(n_annot, 0)] = step1_int;

    let mut delete_values = MatF::zeros(n_blocks, n_annot + 1);
    for k in 0..n_blocks {
        delete_values[(k, n_annot)] = step1.delete_values[(k, n_annot)];
    }

    for k in 0..n_blocks {
        for j in 0..n_annot {
            let adj = c * (step1.delete_values[(k, n_annot)] - step1_int);
            delete_values[(k, j)] = step2.delete_values[(k, j)] - adj;
        }
    }

    let (jknife_se, jknife_cov) = jackknife_se_cov(&est, &delete_values, n_blocks);

    Ok(JackknifeResult {
        est,
        jknife_se,
        jknife_cov,
        delete_values,
        separators: step2.separators.clone(),
    })
}

#[allow(clippy::too_many_arguments)]
pub fn run_h2_ldsc(
    chi2: &ColF,
    ref_l2: &ColF,
    w_l2: &ColF,
    n_vec: &ColF,
    m_snps: f64,
    n_blocks: usize,
    two_step: Option<f64>,
    fixed_intercept: Option<f64>,
) -> Result<H2Result> {
    let n = col_len(chi2);
    if col_len(ref_l2) != n || col_len(w_l2) != n || col_len(n_vec) != n {
        bail!("run_h2_ldsc: input length mismatch");
    }
    let nbar = col_mean(n_vec);

    let intercept0 = fixed_intercept.unwrap_or(1.0);
    let tot_agg = aggregate(chi2, ref_l2, n_vec, m_snps, intercept0);
    let initial_w = ldsc_weights(ref_l2, w_l2, n_vec, m_snps, tot_agg, intercept0);

    let mut x = MatF::zeros(n, 1);
    for i in 0..n {
        x[(i, 0)] = n_vec[(i, 0)] * ref_l2[(i, 0)] / nbar;
    }

    let (y_reg, x_reg) = if let Some(intercept) = fixed_intercept {
        let mut y_adj = chi2.clone();
        for i in 0..n {
            y_adj[(i, 0)] -= intercept;
        }
        (y_adj, x)
    } else {
        let mut x_i = MatF::zeros(n, 2);
        for i in 0..n {
            x_i[(i, 0)] = x[(i, 0)];
            x_i[(i, 1)] = 1.0;
        }
        (chi2.clone(), x_i)
    };

    if two_step.is_some() && fixed_intercept.is_some() {
        bail!("twostep is not compatible with constrained intercept");
    }

    let jknife = if let Some(twostep) = two_step {
        let mut mask: Vec<bool> = Vec::with_capacity(n);
        for i in 0..n {
            mask.push(chi2[(i, 0)] < twostep);
        }
        let n1 = mask.iter().filter(|&&b| b).count();

        let mut x1 = MatF::zeros(n1, x_reg.ncols());
        let mut y1 = col_zeros(n1);
        let mut w1 = col_zeros(n1);
        let mut n1v = col_zeros(n1);
        let mut idx = 0usize;
        for i in 0..n {
            if mask[i] {
                for j in 0..x_reg.ncols() {
                    x1[(idx, j)] = x_reg[(i, j)];
                }
                y1[(idx, 0)] = y_reg[(i, 0)];
                w1[(idx, 0)] = w_l2[(i, 0)];
                n1v[(idx, 0)] = n_vec[(i, 0)];
                idx += 1;
            }
        }
        let mut initial_w1_vec = Vec::with_capacity(n1);
        for (i, keep) in mask.iter().enumerate() {
            if *keep {
                initial_w1_vec.push(initial_w[(i, 0)]);
            }
        }
        let initial_w1 = col_from_vec(initial_w1_vec);

        let mut x1_col0 = col_zeros(n1);
        for i in 0..n1 {
            x1_col0[(i, 0)] = x1[(i, 0)];
        }
        let update_func1 = move |coef: &ColF| {
            let hsq = m_snps * coef[(0, 0)] / nbar;
            let intercept = coef[(1, 0)];
            ldsc_weights(&x1_col0, &w1, &n1v, m_snps, hsq, intercept)
        };

        let step1 = irwls_ldsc(&x1, &y1, &initial_w1, update_func1, n_blocks, None)?;
        let step1_int = step1.est[(1, 0)];

        let mut y2 = y_reg.clone();
        for i in 0..n {
            y2[(i, 0)] -= step1_int;
        }
        let mut x2 = MatF::zeros(n, 1);
        for i in 0..n {
            x2[(i, 0)] = x_reg[(i, 0)];
        }

        let update_func2 = |coef: &ColF| {
            let hsq = m_snps * coef[(0, 0)] / nbar;
            ldsc_weights(ref_l2, w_l2, n_vec, m_snps, hsq, step1_int)
        };

        let seps = update_separators(&step1.separators, &mask);
        let step2 = irwls_ldsc(&x2, &y2, &initial_w, update_func2, n_blocks, Some(&seps))?;

        let mut num = 0.0f64;
        let mut denom = 0.0f64;
        for i in 0..n {
            let w = initial_w[(i, 0)];
            let x0 = x2[(i, 0)];
            num += w * x0;
            denom += w * x0 * x0;
        }
        let c = num / denom;
        combine_twostep(&step1, &step2, c)?
    } else {
        let update_func = |coef: &ColF| {
            let hsq = m_snps * coef[(0, 0)] / nbar;
            let intercept = fixed_intercept.unwrap_or(coef[(1, 0)]);
            ldsc_weights(ref_l2, w_l2, n_vec, m_snps, hsq, intercept)
        };
        irwls_ldsc(&x_reg, &y_reg, &initial_w, update_func, n_blocks, None)?
    };

    let coef = jknife.est[(0, 0)] / nbar;
    let coef_cov = jknife.jknife_cov[(0, 0)] / (nbar * nbar);
    let h2 = m_snps * coef;
    let h2_se = (m_snps * m_snps * coef_cov).sqrt();

    let (intercept, intercept_se) = if let Some(fixed) = fixed_intercept {
        (fixed, f64::NAN)
    } else {
        (jknife.est[(1, 0)], jknife.jknife_se[(1, 0)])
    };

    let mean_chi2 = col_mean(chi2);
    let lambda_gc = {
        let mut vals = Vec::with_capacity(n);
        for i in 0..n {
            vals.push(chi2[(i, 0)]);
        }
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let n = vals.len();
        let mid = if n == 0 {
            f64::NAN
        } else if n % 2 == 1 {
            vals[n / 2]
        } else {
            let hi = vals[n / 2];
            let lo = vals[n / 2 - 1];
            (lo + hi) / 2.0
        };
        mid / 0.4549
    };
    let ratio = if mean_chi2 > 1.0 && fixed_intercept.is_none() {
        let ratio = (intercept - 1.0) / (mean_chi2 - 1.0);
        let ratio_se = intercept_se / (mean_chi2 - 1.0);
        Some((ratio, ratio_se))
    } else {
        None
    };

    Ok(H2Result {
        h2,
        h2_se,
        intercept,
        intercept_se,
        mean_chi2,
        lambda_gc,
        ratio,
    })
}

pub(crate) fn irwls_ldsc(
    x: &MatF,
    y: &ColF,
    initial_w: &ColF,
    update_func: impl Fn(&ColF) -> ColF,
    n_blocks: usize,
    separators: Option<&[usize]>,
) -> Result<JackknifeResult> {
    let mut w = initial_w.clone();
    for i in 0..col_len(&w) {
        w[(i, 0)] = w[(i, 0)].sqrt();
    }
    for _ in 0..2 {
        let coef = wls(x, y, &w)?;
        let new_w = update_func(&coef);
        w = new_w;
        for i in 0..col_len(&w) {
            w[(i, 0)] = w[(i, 0)].sqrt();
        }
    }

    let (xw, yw) = weight_xy(x, y, &w)?;
    jackknife_fast(&xw, &yw, n_blocks, separators)
}
