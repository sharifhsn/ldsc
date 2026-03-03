//! LDSC h2 regression logic aligned to the Python implementation.
use anyhow::{Result, bail};
use ndarray::{Array1, Array2, Axis, s};
use ndarray_linalg::{LeastSquaresSvd, Solve};

#[derive(Debug, Clone)]
pub struct JackknifeResult {
    pub est: Array1<f64>,
    pub jknife_se: Array1<f64>,
    pub jknife_cov: Array2<f64>,
    pub delete_values: Array2<f64>,
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

fn mean(arr: &Array1<f64>) -> f64 {
    arr.mean().unwrap_or(f64::NAN)
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

pub(crate) fn weight_xy(
    x: &Array2<f64>,
    y: &Array1<f64>,
    w: &Array1<f64>,
) -> Result<(Array2<f64>, Array1<f64>)> {
    let n = x.nrows();
    if y.len() != n || w.len() != n {
        bail!("weight_xy: shape mismatch");
    }
    if w.iter().any(|&v| v <= 0.0) {
        bail!("weights must be > 0");
    }
    let sum_w = w.sum();
    let mut xw = Array2::<f64>::zeros(x.raw_dim());
    let mut yw = Array1::<f64>::zeros(n);
    for i in 0..n {
        let wi = w[i] / sum_w;
        xw.row_mut(i).assign(&(&x.row(i) * wi));
        yw[i] = y[i] * wi;
    }
    Ok((xw, yw))
}

fn wls(x: &Array2<f64>, y: &Array1<f64>, w: &Array1<f64>) -> Result<Array1<f64>> {
    let (xw, yw) = weight_xy(x, y, w)?;
    let res = xw.least_squares(&yw)?;
    Ok(res.solution)
}

pub(crate) fn jackknife_fast(
    x: &Array2<f64>,
    y: &Array1<f64>,
    n_blocks: usize,
    separators: Option<&[usize]>,
) -> Result<JackknifeResult> {
    let n = x.nrows();
    let p = x.ncols();
    if y.len() != n {
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

    let mut xty_blocks = Array2::<f64>::zeros((n_blocks, p));
    let mut xtx_blocks = Vec::with_capacity(n_blocks);
    for i in 0..n_blocks {
        let lo = seps[i];
        let hi = seps[i + 1];
        let xb = x.slice(s![lo..hi, ..]);
        let yb = y.slice(s![lo..hi]);
        let xty = xb.t().dot(&yb);
        xty_blocks.row_mut(i).assign(&xty);
        let xtx = xb.t().dot(&xb);
        xtx_blocks.push(xtx);
    }

    let xty_tot = xty_blocks.sum_axis(Axis(0));
    let mut xtx_tot = Array2::<f64>::zeros((p, p));
    for xtx in &xtx_blocks {
        xtx_tot += xtx;
    }
    let est = xtx_tot.solve_into(xty_tot.clone())?;

    let mut delete_values = Array2::<f64>::zeros((n_blocks, p));
    for (i, xtx_block) in xtx_blocks.iter().enumerate() {
        let xty_del = &xty_tot - &xty_blocks.row(i);
        let xtx_del = &xtx_tot - xtx_block;
        let coef = xtx_del.solve_into(xty_del)?;
        delete_values.row_mut(i).assign(&coef);
    }

    let mut pseudovalues = Array2::<f64>::zeros((n_blocks, p));
    for i in 0..n_blocks {
        let pv = &est * (n_blocks as f64) - &delete_values.row(i) * ((n_blocks - 1) as f64);
        pseudovalues.row_mut(i).assign(&pv);
    }

    let mean_pv = pseudovalues.mean_axis(Axis(0)).unwrap();
    let centered = Array2::from_shape_fn((n_blocks, p), |(i, j)| pseudovalues[[i, j]] - mean_pv[j]);
    let jknife_cov = centered.t().dot(&centered) / ((n_blocks - 1) as f64 * n_blocks as f64);
    let mut jknife_se = Array1::<f64>::zeros(p);
    for j in 0..p {
        jknife_se[j] = jknife_cov[[j, j]].sqrt();
    }

    Ok(JackknifeResult {
        est,
        jknife_se,
        jknife_cov,
        delete_values,
        separators: seps,
    })
}

pub(crate) fn ldsc_weights(
    ld: &Array1<f64>,
    w_ld: &Array1<f64>,
    n_vec: &Array1<f64>,
    m: f64,
    hsq: f64,
    intercept: f64,
) -> Array1<f64> {
    let mut hsq = hsq;
    hsq = hsq.clamp(0.0, 1.0);
    let mut out = Array1::<f64>::zeros(ld.len());
    for i in 0..ld.len() {
        let ld_i = ld[i].max(1.0);
        let w_i = w_ld[i].max(1.0);
        let c = hsq * n_vec[i] / m;
        let het_w = 1.0 / (2.0 * (intercept + c * ld_i).powi(2));
        let oc_w = 1.0 / w_i;
        out[i] = het_w * oc_w;
    }
    out
}

pub(crate) fn aggregate(
    y: &Array1<f64>,
    x: &Array1<f64>,
    n_vec: &Array1<f64>,
    m: f64,
    intercept: f64,
) -> f64 {
    let num = m * (mean(y) - intercept);
    let denom = mean(&(x * n_vec));
    num / denom
}

pub(crate) fn combine_twostep(
    step1: &JackknifeResult,
    step2: &JackknifeResult,
    c: f64,
) -> Result<JackknifeResult> {
    let n_blocks = step1.delete_values.nrows();
    let n_annot = step2.delete_values.ncols();
    let step1_int = step1.est[n_annot];

    let mut est = Array1::<f64>::zeros(n_annot + 1);
    for j in 0..n_annot {
        est[j] = step2.est[j];
    }
    est[n_annot] = step1_int;

    let mut delete_values = Array2::<f64>::zeros((n_blocks, n_annot + 1));
    delete_values
        .column_mut(n_annot)
        .assign(&step1.delete_values.column(n_annot));

    for k in 0..n_blocks {
        for j in 0..n_annot {
            let adj = c * (step1.delete_values[[k, n_annot]] - step1_int);
            delete_values[[k, j]] = step2.delete_values[[k, j]] - adj;
        }
    }

    let mut pseudovalues = Array2::<f64>::zeros((n_blocks, n_annot + 1));
    for k in 0..n_blocks {
        let pv = &est * (n_blocks as f64) - &delete_values.row(k) * ((n_blocks - 1) as f64);
        pseudovalues.row_mut(k).assign(&pv);
    }
    let mean_pv = pseudovalues.mean_axis(Axis(0)).unwrap();
    let centered = Array2::from_shape_fn((n_blocks, n_annot + 1), |(i, j)| {
        pseudovalues[[i, j]] - mean_pv[j]
    });
    let jknife_cov = centered.t().dot(&centered) / ((n_blocks - 1) as f64 * n_blocks as f64);
    let mut jknife_se = Array1::<f64>::zeros(n_annot + 1);
    for j in 0..n_annot + 1 {
        jknife_se[j] = jknife_cov[[j, j]].sqrt();
    }

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
    chi2: &Array1<f64>,
    ref_l2: &Array1<f64>,
    w_l2: &Array1<f64>,
    n_vec: &Array1<f64>,
    m_snps: f64,
    n_blocks: usize,
    two_step: Option<f64>,
    fixed_intercept: Option<f64>,
) -> Result<H2Result> {
    let n = chi2.len();
    if ref_l2.len() != n || w_l2.len() != n || n_vec.len() != n {
        bail!("run_h2_ldsc: input length mismatch");
    }
    let nbar = mean(n_vec);

    let intercept0 = fixed_intercept.unwrap_or(1.0);
    let tot_agg = aggregate(chi2, ref_l2, n_vec, m_snps, intercept0);
    let initial_w = ldsc_weights(ref_l2, w_l2, n_vec, m_snps, tot_agg, intercept0);

    let mut x = Array2::<f64>::zeros((n, 1));
    for i in 0..n {
        x[[i, 0]] = n_vec[i] * ref_l2[i] / nbar;
    }

    let (y_reg, x_reg) = if let Some(intercept) = fixed_intercept {
        let y_adj = chi2.mapv(|c| c - intercept);
        (y_adj, x)
    } else {
        let mut x_i = Array2::<f64>::zeros((n, 2));
        x_i.column_mut(0).assign(&x.column(0));
        x_i.column_mut(1).fill(1.0);
        (chi2.clone(), x_i)
    };

    if two_step.is_some() && fixed_intercept.is_some() {
        bail!("twostep is not compatible with constrained intercept");
    }

    let jknife = if let Some(twostep) = two_step {
        let mask: Vec<bool> = chi2.iter().map(|&c| c < twostep).collect();
        let n1 = mask.iter().filter(|&&b| b).count();

        let mut x1 = Array2::<f64>::zeros((n1, x_reg.ncols()));
        let mut y1 = Array1::<f64>::zeros(n1);
        let mut w1 = Array1::<f64>::zeros(n1);
        let mut n1v = Array1::<f64>::zeros(n1);
        let mut idx = 0usize;
        for i in 0..n {
            if mask[i] {
                x1.row_mut(idx).assign(&x_reg.row(i));
                y1[idx] = y_reg[i];
                w1[idx] = w_l2[i];
                n1v[idx] = n_vec[i];
                idx += 1;
            }
        }
        let initial_w1_vec: Vec<f64> = initial_w
            .iter()
            .zip(mask.iter())
            .filter_map(|(&v, keep)| if *keep { Some(v) } else { None })
            .collect();
        let initial_w1 = Array1::from_vec(initial_w1_vec);

        let x1_col0 = x1.column(0).to_owned();
        let update_func1 = move |coef: &Array1<f64>| {
            let hsq = m_snps * coef[0] / nbar;
            let intercept = coef[1];
            ldsc_weights(&x1_col0, &w1, &n1v, m_snps, hsq, intercept)
        };

        let step1 = irwls_ldsc(&x1, &y1, &initial_w1, update_func1, n_blocks, None)?;
        let step1_int = step1.est[1];

        let y2 = y_reg.mapv(|c| c - step1_int);
        let x2 = x_reg.slice(s![.., 0..1]).to_owned(); // remove intercept column

        let update_func2 = |coef: &Array1<f64>| {
            let hsq = m_snps * coef[0] / nbar;
            ldsc_weights(ref_l2, w_l2, n_vec, m_snps, hsq, step1_int)
        };

        let seps = update_separators(&step1.separators, &mask);
        let step2 = irwls_ldsc(&x2, &y2, &initial_w, update_func2, n_blocks, Some(&seps))?;

        let num = (&initial_w * &x2.column(0)).sum();
        let denom = (&initial_w * &x2.column(0).mapv(|v| v * v)).sum();
        let c = num / denom;
        combine_twostep(&step1, &step2, c)?
    } else {
        let update_func = |coef: &Array1<f64>| {
            let hsq = m_snps * coef[0] / nbar;
            let intercept = fixed_intercept.unwrap_or(coef[1]);
            ldsc_weights(ref_l2, w_l2, n_vec, m_snps, hsq, intercept)
        };
        irwls_ldsc(&x_reg, &y_reg, &initial_w, update_func, n_blocks, None)?
    };

    let coef = jknife.est[0] / nbar;
    let coef_cov = jknife.jknife_cov[[0, 0]] / (nbar * nbar);
    let h2 = m_snps * coef;
    let h2_se = (m_snps * m_snps * coef_cov).sqrt();

    let (intercept, intercept_se) = if let Some(fixed) = fixed_intercept {
        (fixed, f64::NAN)
    } else {
        (jknife.est[1], jknife.jknife_se[1])
    };

    let mean_chi2 = mean(chi2);
    let lambda_gc = {
        let mut vals = chi2.to_vec();
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
    x: &Array2<f64>,
    y: &Array1<f64>,
    initial_w: &Array1<f64>,
    update_func: impl Fn(&Array1<f64>) -> Array1<f64>,
    n_blocks: usize,
    separators: Option<&[usize]>,
) -> Result<JackknifeResult> {
    let mut w = initial_w.mapv(|v| v.sqrt());
    for _ in 0..2 {
        let coef = wls(x, y, &w)?;
        let new_w = update_func(&coef);
        w = new_w.mapv(|v| v.sqrt());
    }

    let (xw, yw) = weight_xy(x, y, &w)?;
    jackknife_fast(&xw, &yw, n_blocks, separators)
}
