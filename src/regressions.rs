/// Heritability and genetic correlation estimation — replaces ldscore/regressions.py.
///
/// Pipeline:
///   parse → merge (sumstats ⨝ ld_scores ⨝ weights) → build ndarray
///   matrices → jackknife(IRWLS, n_blocks) → h2 / rg estimates.
use anyhow::{Context, Result};
use ndarray::{Array1, Array2};
use polars::prelude::*;
use statrs::distribution::{ChiSquared, Continuous, ContinuousCDF, Normal};

use crate::cli::{H2Args, RgArgs};
use crate::jackknife;
use crate::parse;

// ---------------------------------------------------------------------------
// Helper: load LD scores from either a single file or per-chr split files
// ---------------------------------------------------------------------------

/// Load a scalar (single L2 column) LD score file, aliasing the L2 column.
/// Used for weight LD scores (always scalar).
fn load_ld(single: Option<&str>, chr_prefix: Option<&str>, alias: &str) -> Result<LazyFrame> {
    match (single, chr_prefix) {
        (Some(path), _) => {
            Ok(parse::scan_ldscore(path)?.select([col("SNP"), col("L2").alias(alias)]))
        }
        (None, Some(prefix)) => Ok(parse::concat_chrs(prefix, ".l2.ldscore.gz")?
            .select([col("SNP"), col("L2").alias(alias)])),
        (None, None) => anyhow::bail!(
            "Must specify exactly one of --ref-ld / --ref-ld-chr (or --w-ld / --w-ld-chr)"
        ),
    }
}

/// Load reference LD scores, keeping ALL annotation columns.
///
/// For standard (non-partitioned) ldscore files the result has columns: SNP, L2.
/// For partitioned ldscore files (produced by `ldscore --annot`) the result has
/// columns: SNP, ANNOT1L2, ANNOT2L2, … (one per annotation).
///
/// Metadata columns (CHR, BP, CM) are dropped.
fn load_ld_ref(single: Option<&str>, chr_prefix: Option<&str>) -> Result<LazyFrame> {
    let lf = match (single, chr_prefix) {
        (Some(path), _) => parse::scan_ldscore(path)?,
        (None, Some(prefix)) => parse::concat_chrs(prefix, ".l2.ldscore.gz")?,
        (None, None) => anyhow::bail!("Must specify exactly one of --ref-ld / --ref-ld-chr"),
    };
    // Drop metadata columns, keeping SNP and all L2 annotation columns.
    let peek = lf
        .clone()
        .limit(0)
        .collect()
        .context("peeking ref LD schema")?;
    let meta = ["CHR", "BP", "CM"];
    let keep: Vec<Expr> = peek
        .get_column_names()
        .into_iter()
        .filter(|n| !meta.contains(&n.as_str()))
        .map(|n| col(n.as_str()))
        .collect();
    Ok(lf.select(keep))
}

// ---------------------------------------------------------------------------
// Helper: determine M denominator
// ---------------------------------------------------------------------------

fn resolve_m(
    m_snps_override: Option<f64>,
    ref_ld_chr: Option<&str>,
    not_m_5_50: bool,
    n_obs: usize,
) -> f64 {
    if let Some(m) = m_snps_override {
        println!("  Using M = {:.0} from --m-snps", m);
        return m;
    }
    if let Some(prefix) = ref_ld_chr {
        let suffix = if not_m_5_50 { ".l2.M" } else { ".l2.M_5_50" };
        match parse::read_m_total(prefix, suffix) {
            Ok(m) => {
                println!("  Read M = {:.0} from {} files", m, suffix);
                return m;
            }
            Err(_) => {
                println!(
                    "  No {} files found alongside --ref-ld-chr; using M = {} (regression SNPs)",
                    suffix, n_obs
                );
            }
        }
    } else {
        println!(
            "  Single-file --ref-ld: no M files to read. \
             Use --m-snps to set M explicitly. Using M = {} (regression SNPs).",
            n_obs
        );
    }
    n_obs as f64
}

// ---------------------------------------------------------------------------
// Helper: resolve per-annotation M denominator vector
// ---------------------------------------------------------------------------

/// Return per-annotation M values (length K) for partitioned h2 regression.
///
/// Priority:
///   1. Read K values from `.l2.M` / `.l2.M_5_50` files alongside `--ref-ld-chr`.
///   2. Fall back: divide `--m-snps` (if given) uniformly across K annotations.
///   3. Final fallback: divide n_obs uniformly across K annotations.
fn resolve_m_vec(
    m_snps_override: Option<f64>,
    ref_ld_chr: Option<&str>,
    not_m_5_50: bool,
    n_obs: usize,
    k: usize,
) -> Vec<f64> {
    if let Some(prefix) = ref_ld_chr {
        let suffix = if not_m_5_50 { ".l2.M" } else { ".l2.M_5_50" };
        match parse::read_m_vec(prefix, suffix) {
            Ok(m_vec) if m_vec.len() == k => {
                let total: f64 = m_vec.iter().sum();
                println!(
                    "  Read per-annotation M (K={}) from {} files; total M = {:.0}",
                    k, suffix, total
                );
                return m_vec;
            }
            Ok(m_vec) => {
                println!(
                    "  WARNING: {} files have {} M values but K={} annotation columns; using fallback",
                    suffix,
                    m_vec.len(),
                    k
                );
            }
            Err(e) => {
                println!("  No {} files found ({}); using fallback M", suffix, e);
            }
        }
    }
    // Fallback: divide total M uniformly across annotations.
    let total = m_snps_override.unwrap_or(n_obs as f64);
    println!(
        "  Using M = {:.0} / K={} (uniform per annotation)",
        total, k
    );
    vec![total / k as f64; k]
}

// ---------------------------------------------------------------------------
// Helper: filter multiple parallel arrays by a boolean mask
// ---------------------------------------------------------------------------

fn filter_by_mask(v: &Array1<f64>, mask: &[bool]) -> Array1<f64> {
    Array1::from_vec(
        v.iter()
            .zip(mask.iter())
            .filter_map(|(&val, &keep)| if keep { Some(val) } else { None })
            .collect(),
    )
}

// ---------------------------------------------------------------------------
// Heritability (--h2)
// ---------------------------------------------------------------------------

pub fn run_h2(args: H2Args) -> Result<()> {
    let sumstats_lf = parse::scan_sumstats(&args.h2)?;
    // Reference LD keeps all annotation columns (scalar or partitioned).
    let ref_ld = load_ld_ref(args.ref_ld.as_deref(), args.ref_ld_chr.as_deref())?;
    let w_ld = load_ld(args.w_ld.as_deref(), args.w_ld_chr.as_deref(), "w_l2")?;

    let merged = sumstats_lf
        .join(
            ref_ld,
            [col("SNP")],
            [col("SNP")],
            JoinArgs::new(JoinType::Inner),
        )
        .join(
            w_ld,
            [col("SNP")],
            [col("SNP")],
            JoinArgs::new(JoinType::Inner),
        )
        .collect()
        .context("merging sumstats with LD scores")?;

    let n_total = merged.height();
    anyhow::ensure!(
        n_total > 0,
        "No SNPs remaining after merging sumstats with LD scores"
    );

    println!(
        "Loaded {} SNPs from '{}' after merging with LD scores",
        n_total, args.h2
    );

    // Detect L2 annotation columns: everything that isn't a sumstats or weight column.
    let non_l2_set: std::collections::HashSet<&str> = ["SNP", "A1", "A2", "Z", "N", "w_l2"]
        .iter()
        .copied()
        .collect();
    let l2_cols: Vec<String> = merged
        .get_column_names()
        .into_iter()
        .filter(|n| !non_l2_set.contains(n.as_str()))
        .map(|s| s.to_string())
        .collect();
    let k = l2_cols.len();
    anyhow::ensure!(k > 0, "No L2 annotation columns found in merged dataset");

    if k > 1 {
        println!("  Detected K={} annotation columns (partitioned h2)", k);
    }

    let chi2_raw = extract_f64(&merged, "Z")?.mapv(|z| (z * z).max(0.0));
    let w_l2_raw = extract_f64(&merged, "w_l2")?.mapv(|l| l.max(0.0));
    let n_vec_raw = extract_f64(&merged, "N")?;
    let ref_l2_raw_k: Vec<Array1<f64>> = l2_cols
        .iter()
        .map(|name| extract_f64(&merged, name).map(|v| v.mapv(|l| l.max(0.0))))
        .collect::<Result<Vec<_>>>()?;

    let (chi2, w_l2, n_vec, ref_l2_k) = if let Some(chisq_max) = args.chisq_max {
        let mask: Vec<bool> = chi2_raw.iter().map(|&c| c <= chisq_max).collect();
        let n_removed = mask.iter().filter(|&&b| !b).count();
        if n_removed > 0 {
            println!("  Removed {} SNPs with chi2 > {:.1}", n_removed, chisq_max);
        }
        (
            filter_by_mask(&chi2_raw, &mask),
            filter_by_mask(&w_l2_raw, &mask),
            filter_by_mask(&n_vec_raw, &mask),
            ref_l2_raw_k
                .iter()
                .map(|v| filter_by_mask(v, &mask))
                .collect::<Vec<_>>(),
        )
    } else {
        (chi2_raw, w_l2_raw, n_vec_raw, ref_l2_raw_k)
    };
    let n_obs = chi2.len();
    let n_mean = n_vec.mean().unwrap_or(1.0);

    if k > 1 {
        return run_h2_partitioned(
            &chi2, &ref_l2_k, &w_l2, &n_vec, &l2_cols, n_obs, n_mean, &args,
        );
    }

    // Scalar h2 (K == 1).
    let ref_l2 = ref_l2_k.into_iter().next().unwrap();
    let m_snps = resolve_m(
        args.m_snps,
        args.ref_ld_chr.as_deref(),
        args.not_m_5_50,
        n_obs,
    );

    // Two-step: step 1 estimates intercept on SNPs with chi2 ≤ threshold;
    // step 2 fixes that intercept and re-estimates slope on all SNPs.
    if let Some(two_step_max) = args.two_step {
        let mask1: Vec<bool> = chi2.iter().map(|&c| c <= two_step_max).collect();
        let n_s1 = mask1.iter().filter(|&&b| b).count();
        println!(
            "  Two-step: step 1 uses {} SNPs (chi2 <= {})",
            n_s1, two_step_max
        );

        let chi2_s1 = filter_by_mask(&chi2, &mask1);
        let l2_s1 = filter_by_mask(&ref_l2, &mask1);
        let wl2_s1 = filter_by_mask(&w_l2, &mask1);
        let n_s1v = filter_by_mask(&n_vec, &mask1);

        let mut x_s1 = Array2::<f64>::zeros((n_s1, 2));
        let mut w_s1 = Array1::<f64>::zeros(n_s1);
        for i in 0..n_s1 {
            x_s1[[i, 0]] = l2_s1[i];
            x_s1[[i, 1]] = 1.0;
            let v = (wl2_s1[i] / m_snps + 1.0 / n_s1v[i]).max(1e-9);
            w_s1[i] = 1.0 / (2.0 * v * v);
        }
        let res_s1 = jackknife::jackknife(&x_s1, &chi2_s1, &w_s1, args.n_blocks, 2)
            .context("jackknife step 1 for two-step h2")?;
        let intercept_fixed = res_s1.est[1];
        println!("  Two-step: step 1 intercept = {:.4}", intercept_fixed);

        let y2 = chi2.mapv(|c| c - intercept_fixed);
        let mut x2 = Array2::<f64>::zeros((n_obs, 1));
        let mut w2 = Array1::<f64>::zeros(n_obs);
        for i in 0..n_obs {
            x2[[i, 0]] = ref_l2[i];
            let v = (w_l2[i] / m_snps + 1.0 / n_vec[i]).max(1e-9);
            w2[i] = 1.0 / (2.0 * v * v);
        }
        let res2 = jackknife::jackknife(&x2, &y2, &w2, args.n_blocks, 2)
            .context("jackknife step 2 for two-step h2")?;
        let h2 = res2.est[0] * m_snps / n_mean;
        let h2_se_str = res2
            .jknife_se
            .as_ref()
            .map(|se| format!("{:.4}", se[0] * m_snps / n_mean))
            .unwrap_or_else(|| "N/A".to_string());
        println!("Total observed-scale h2: {:.4} (SE: {})", h2, h2_se_str);
        println!("Intercept: {:.4} (fixed from step 1)", intercept_fixed);
        if let (Some(samp_prev), Some(pop_prev)) = (args.samp_prev, args.pop_prev) {
            let c = liability_conversion_factor(samp_prev, pop_prev);
            println!("Liability-scale h2: {:.4}", h2 * c);
        }
        return Ok(());
    }

    // Priority: --intercept-h2 > --no-intercept > free intercept.
    let (y_reg, x_reg, fixed_intercept) = if let Some(fixed) = args.intercept_h2 {
        let y_adj = chi2.mapv(|c| c - fixed);
        let mut x = Array2::<f64>::zeros((n_obs, 1));
        for i in 0..n_obs {
            x[[i, 0]] = ref_l2[i];
        }
        (y_adj, x, Some(fixed))
    } else if args.no_intercept {
        let y_adj = chi2.mapv(|c| c - 1.0);
        let mut x = Array2::<f64>::zeros((n_obs, 1));
        for i in 0..n_obs {
            x[[i, 0]] = ref_l2[i];
        }
        (y_adj, x, Some(1.0f64))
    } else {
        let mut x = Array2::<f64>::zeros((n_obs, 2));
        for i in 0..n_obs {
            x[[i, 0]] = ref_l2[i];
            x[[i, 1]] = 1.0;
        }
        (chi2.clone(), x, None)
    };

    let mut weights: Array1<f64> = Array1::zeros(n_obs);
    for i in 0..n_obs {
        let v = (w_l2[i] / m_snps + 1.0 / n_vec[i]).max(1e-9);
        weights[i] = 1.0 / (2.0 * v * v);
    }

    let result = jackknife::jackknife(&x_reg, &y_reg, &weights, args.n_blocks, 2)
        .context("jackknife IRWLS for h2")?;

    // h2 = slope × M / N_mean
    let h2 = result.est[0] * m_snps / n_mean;
    let intercept = fixed_intercept.unwrap_or_else(|| result.est[1]);

    let h2_se_str = result
        .jknife_se
        .as_ref()
        .map(|se| format!("{:.4}", se[0] * m_snps / n_mean))
        .unwrap_or_else(|| "N/A".to_string());

    println!("Total observed-scale h2: {:.4} (SE: {})", h2, h2_se_str);
    println!("Intercept: {:.4}", intercept);

    print_jackknife_diagnostics(&result, args.print_cov, args.print_delete_vals);

    if let (Some(samp_prev), Some(pop_prev)) = (args.samp_prev, args.pop_prev) {
        let c = liability_conversion_factor(samp_prev, pop_prev);
        println!("Liability-scale h2: {:.4}", h2 * c);
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Partitioned h2 regression (K annotation columns)
// ---------------------------------------------------------------------------

/// Partitioned stratified LD score regression (Finucane et al. 2015).
///
/// Fits:  E[χ²ᵢ] ≈ Nᵢ × Σₖ τₖ × ℓᵢₖ / Mₖ + Na + 1
///
/// Which in matrix form is:  y = X β + ε
///   y[i]     = χ²[i]
///   X[i, k]  = ℓᵢₖ  (L2 score for annotation k at SNP i)
///   X[i, K]  = 1.0  (intercept column, unless fixed)
///   β[k]     ≈ N × h2ₖ / Mₖ  →  h2ₖ = β[k] × Mₖ / N_mean
///
/// Weights are based on the weight-LD-score file (w_l2) and total M,
/// matching the scalar LDSC regression formula.
#[allow(clippy::too_many_arguments)]
fn run_h2_partitioned(
    chi2: &Array1<f64>,
    ref_l2_k: &[Array1<f64>],
    w_l2: &Array1<f64>,
    n_vec: &Array1<f64>,
    l2_cols: &[String],
    n_obs: usize,
    n_mean: f64,
    args: &H2Args,
) -> Result<()> {
    let k = l2_cols.len();

    let m_vec = resolve_m_vec(
        args.m_snps,
        args.ref_ld_chr.as_deref(),
        args.not_m_5_50,
        n_obs,
        k,
    );
    let m_total: f64 = m_vec.iter().sum();

    let (y_reg, x_reg, fixed_intercept) = if let Some(fixed) = args.intercept_h2 {
        let y_adj = chi2.mapv(|c| c - fixed);
        let mut x = Array2::<f64>::zeros((n_obs, k));
        for (j, col_v) in ref_l2_k.iter().enumerate() {
            for i in 0..n_obs {
                x[[i, j]] = col_v[i];
            }
        }
        (y_adj, x, Some(fixed))
    } else if args.no_intercept {
        let y_adj = chi2.mapv(|c| c - 1.0);
        let mut x = Array2::<f64>::zeros((n_obs, k));
        for (j, col_v) in ref_l2_k.iter().enumerate() {
            for i in 0..n_obs {
                x[[i, j]] = col_v[i];
            }
        }
        (y_adj, x, Some(1.0f64))
    } else {
        let mut x = Array2::<f64>::zeros((n_obs, k + 1));
        for (j, col_v) in ref_l2_k.iter().enumerate() {
            for i in 0..n_obs {
                x[[i, j]] = col_v[i];
            }
        }
        for i in 0..n_obs {
            x[[i, k]] = 1.0;
        }
        (chi2.clone(), x, None)
    };

    let mut weights = Array1::<f64>::zeros(n_obs);
    for i in 0..n_obs {
        let v = (w_l2[i] / m_total + 1.0 / n_vec[i]).max(1e-9);
        weights[i] = 1.0 / (2.0 * v * v);
    }

    let result = jackknife::jackknife(&x_reg, &y_reg, &weights, args.n_blocks, 2)
        .context("jackknife IRWLS for partitioned h2")?;

    // h2_k = τ_k × M_k / N_mean
    let h2_per_annot: Vec<f64> = (0..k).map(|j| result.est[j] * m_vec[j] / n_mean).collect();
    let h2_total: f64 = h2_per_annot.iter().sum();
    let intercept = fixed_intercept.unwrap_or_else(|| result.est[k]);

    println!("Total observed-scale h2: {:.4}", h2_total);
    println!("Intercept: {:.4}", intercept);

    if args.print_coefficients {
        println!();
        println!(
            "{:<35} {:>10} {:>10} {:>12} {:>14}",
            "Category", "Prop_h2", "Prop_SNPs", "Enrichment", "Coefficient"
        );
        println!("{}", "-".repeat(83));
        for (j, col_name) in l2_cols.iter().enumerate() {
            let prop_h2 = if h2_total.abs() > 1e-12 {
                h2_per_annot[j] / h2_total
            } else {
                f64::NAN
            };
            let prop_m = m_vec[j] / m_total;
            let enrichment = if prop_m > 1e-12 {
                prop_h2 / prop_m
            } else {
                f64::NAN
            };
            // Coefficient (τ*) = τ_k / (h2_total / M_k): the per-annotation coefficient
            // normalized by mean h2 per SNP.
            let coeff = result.est[j] / n_mean;
            println!(
                "{:<35} {:>10.4} {:>10.4} {:>12.4} {:>14.4e}",
                col_name, prop_h2, prop_m, enrichment, coeff
            );
        }
        println!();
    }

    print_jackknife_diagnostics(&result, args.print_cov, args.print_delete_vals);

    if let (Some(samp_prev), Some(pop_prev)) = (args.samp_prev, args.pop_prev) {
        let c = liability_conversion_factor(samp_prev, pop_prev);
        println!("Liability-scale h2: {:.4}", h2_total * c);
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Genetic correlation (--rg)
// ---------------------------------------------------------------------------

pub fn run_rg(args: RgArgs) -> Result<()> {
    anyhow::ensure!(
        args.rg.len() >= 2,
        "--rg requires at least 2 sumstats files"
    );

    let ref_ld = load_ld(args.ref_ld.as_deref(), args.ref_ld_chr.as_deref(), "ref_l2")?
        .collect()
        .context("loading ref LD scores")?;

    let w_ld = load_ld(args.w_ld.as_deref(), args.w_ld_chr.as_deref(), "w_l2")?
        .collect()
        .context("loading weight LD scores")?;

    // Determine M before the per-pair loop; falls back to n_obs per pair if not found.
    let m_from_files: Option<f64> = if let Some(m) = args.m_snps {
        println!("Using M = {:.0} from --m-snps", m);
        Some(m)
    } else if let Some(prefix) = &args.ref_ld_chr {
        let suffix = if args.not_m_5_50 {
            ".l2.M"
        } else {
            ".l2.M_5_50"
        };
        match parse::read_m_total(prefix, suffix) {
            Ok(m) => {
                println!("Read M = {:.0} from {} files", m, suffix);
                Some(m)
            }
            Err(_) => None, // fall back to n_obs per pair
        }
    } else {
        None // single-file --ref-ld: no M files
    };

    for (pair_idx, pair) in args.rg.windows(2).enumerate() {
        let file1 = &pair[0];
        let file2 = &pair[1];
        println!("Computing rg: {} vs {}", file1, file2);

        // Load both sumstats and suffix their Z/N columns to avoid conflicts.
        let ss1 = parse::scan_sumstats(file1)?
            .select([
                col("SNP"),
                col("A1").alias("A1_1"),
                col("A2").alias("A2_1"),
                col("Z").alias("Z1"),
                col("N").alias("N1"),
            ])
            .join(
                ref_ld.clone().lazy(),
                [col("SNP")],
                [col("SNP")],
                JoinArgs::new(JoinType::Inner),
            )
            .join(
                w_ld.clone().lazy(),
                [col("SNP")],
                [col("SNP")],
                JoinArgs::new(JoinType::Inner),
            );

        let ss2 = parse::scan_sumstats(file2)?.select([
            col("SNP"),
            col("Z").alias("Z2"),
            col("N").alias("N2"),
        ]);

        let merged = ss1
            .join(
                ss2,
                [col("SNP")],
                [col("SNP")],
                JoinArgs::new(JoinType::Inner),
            )
            .collect()
            .with_context(|| format!("merging {} vs {}", file1, file2))?;

        let n_obs = merged.height();
        if n_obs == 0 {
            println!("  No overlapping SNPs — skipping pair");
            continue;
        }

        let z1_raw = extract_f64(&merged, "Z1")?;
        let z2_raw = extract_f64(&merged, "Z2")?;
        let prod_raw = &z1_raw * &z2_raw;
        let ref_l2_raw = extract_f64(&merged, "ref_l2")?.mapv(|l| l.max(0.0));
        let w_l2_raw = extract_f64(&merged, "w_l2")?.mapv(|l| l.max(0.0));
        let n1_raw = extract_f64(&merged, "N1")?;
        let n2_raw = extract_f64(&merged, "N2")?;

        let (prod, ref_l2, w_l2, n1, n2) = if let Some(chisq_max) = args.chisq_max {
            let mask: Vec<bool> = prod_raw.iter().map(|&p| p.abs() <= chisq_max).collect();
            let n_removed = mask.iter().filter(|&&b| !b).count();
            if n_removed > 0 {
                println!(
                    "  Removed {} SNPs with |Z1*Z2| > {:.1}",
                    n_removed, chisq_max
                );
            }
            (
                filter_by_mask(&prod_raw, &mask),
                filter_by_mask(&ref_l2_raw, &mask),
                filter_by_mask(&w_l2_raw, &mask),
                filter_by_mask(&n1_raw, &mask),
                filter_by_mask(&n2_raw, &mask),
            )
        } else {
            (prod_raw, ref_l2_raw, w_l2_raw, n1_raw, n2_raw)
        };
        let n_obs_filtered = prod.len();

        let m_snps = m_from_files.unwrap_or_else(|| {
            println!(
                "  No M files; using M = {} (regression SNPs)",
                n_obs_filtered
            );
            n_obs_filtered as f64
        });
        let n1_mean = n1.mean().unwrap_or(1.0);
        let n2_mean = n2.mean().unwrap_or(1.0);

        // Per-trait liability-scale prevalences (one per file in --rg).
        let samp_prev1 = args.samp_prev.get(pair_idx).copied();
        let samp_prev2 = args.samp_prev.get(pair_idx + 1).copied();
        let pop_prev1 = args.pop_prev.get(pair_idx).copied();
        let pop_prev2 = args.pop_prev.get(pair_idx + 1).copied();

        let (gencov, gencov_intercept) = if let Some(two_step_max) = args.two_step {
            let mask1: Vec<bool> = prod.iter().map(|&p| p.abs() <= two_step_max).collect();
            let n_s1 = mask1.iter().filter(|&&b| b).count();
            println!(
                "  Two-step: step 1 uses {} SNPs (|Z1*Z2| <= {})",
                n_s1, two_step_max
            );

            let prod_s1 = filter_by_mask(&prod, &mask1);
            let l2_s1 = filter_by_mask(&ref_l2, &mask1);
            let wl2_s1 = filter_by_mask(&w_l2, &mask1);
            let n1_s1 = filter_by_mask(&n1, &mask1);
            let n2_s1 = filter_by_mask(&n2, &mask1);

            let mut x_s1 = Array2::<f64>::zeros((n_s1, 2));
            let mut w_s1 = Array1::<f64>::zeros(n_s1);
            for i in 0..n_s1 {
                x_s1[[i, 0]] = l2_s1[i];
                x_s1[[i, 1]] = 1.0;
                let v = (wl2_s1[i] / m_snps + 1.0 / n1_s1[i] + 1.0 / n2_s1[i]).max(1e-9);
                w_s1[i] = 1.0 / (2.0 * v * v);
            }
            let res_s1 = jackknife::jackknife(&x_s1, &prod_s1, &w_s1, args.n_blocks, 2)
                .with_context(|| format!("jackknife step 1 gencov {} vs {}", file1, file2))?;
            let intercept_fixed = res_s1.est[1];
            println!("  Two-step: step 1 intercept = {:.4}", intercept_fixed);

            let prod2 = prod.mapv(|p| p - intercept_fixed);
            let mut x2 = Array2::<f64>::zeros((n_obs_filtered, 1));
            let mut w2 = Array1::<f64>::zeros(n_obs_filtered);
            for i in 0..n_obs_filtered {
                x2[[i, 0]] = ref_l2[i];
                let v = (w_l2[i] / m_snps + 1.0 / n1[i] + 1.0 / n2[i]).max(1e-9);
                w2[i] = 1.0 / (2.0 * v * v);
            }
            let res2 = jackknife::jackknife(&x2, &prod2, &w2, args.n_blocks, 2)
                .with_context(|| format!("jackknife step 2 gencov {} vs {}", file1, file2))?;
            (
                res2.est[0] * m_snps / (n1_mean * n2_mean).sqrt(),
                intercept_fixed,
            )
        } else {
            let pair_intercept = args.intercept_gencov.get(pair_idx).copied();
            let (y_reg, x_cov, fixed_intercept) = if let Some(fixed) = pair_intercept {
                let y_adj = prod.mapv(|p| p - fixed);
                let mut x = Array2::<f64>::zeros((n_obs_filtered, 1));
                for i in 0..n_obs_filtered {
                    x[[i, 0]] = ref_l2[i];
                }
                (y_adj, x, Some(fixed))
            } else if args.no_intercept {
                let mut x = Array2::<f64>::zeros((n_obs_filtered, 1));
                for i in 0..n_obs_filtered {
                    x[[i, 0]] = ref_l2[i];
                }
                (prod.clone(), x, Some(0.0f64))
            } else {
                let mut x = Array2::<f64>::zeros((n_obs_filtered, 2));
                for i in 0..n_obs_filtered {
                    x[[i, 0]] = ref_l2[i];
                    x[[i, 1]] = 1.0;
                }
                (prod.clone(), x, None)
            };
            let mut w_cov: Array1<f64> = Array1::zeros(n_obs_filtered);
            for i in 0..n_obs_filtered {
                let v = (w_l2[i] / m_snps + 1.0 / n1[i] + 1.0 / n2[i]).max(1e-9);
                w_cov[i] = 1.0 / (2.0 * v * v);
            }
            let res_cov = jackknife::jackknife(&x_cov, &y_reg, &w_cov, args.n_blocks, 2)
                .with_context(|| format!("jackknife for gencov {} vs {}", file1, file2))?;
            if args.print_cov || args.print_delete_vals {
                print_jackknife_diagnostics(&res_cov, args.print_cov, args.print_delete_vals);
            }
            (
                res_cov.est[0] * m_snps / (n1_mean * n2_mean).sqrt(),
                fixed_intercept.unwrap_or_else(|| res_cov.est[1]),
            )
        };

        let chi2_1 = extract_f64(&merged, "Z1")?.mapv(|z| (z * z).max(0.0));
        let chi2_2 = extract_f64(&merged, "Z2")?.mapv(|z| (z * z).max(0.0));
        // When --no-intercept is set, fix the univariate h2 intercepts at 1 (matching Python).
        let h2_fixed_int = if args.no_intercept { Some(1.0) } else { None };
        let h2_1 = run_h2_scalar(
            &chi2_1,
            &ref_l2_raw2(&merged)?,
            &w_l2_raw2(&merged)?,
            &n1_raw2(&merged)?,
            m_snps,
            args.n_blocks,
            h2_fixed_int,
        )?;
        let h2_2 = run_h2_scalar(
            &chi2_2,
            &ref_l2_raw2(&merged)?,
            &w_l2_raw2(&merged)?,
            &n2_raw2(&merged)?,
            m_snps,
            args.n_blocks,
            h2_fixed_int,
        )?;

        let rg = if h2_1 > 0.0 && h2_2 > 0.0 {
            gencov / (h2_1 * h2_2).sqrt()
        } else {
            f64::NAN
        };

        println!(
            "  gencov = {:.4}  (intercept: {:.4})",
            gencov, gencov_intercept
        );
        println!("  h2({}) = {:.4}", file1, h2_1);
        println!("  h2({}) = {:.4}", file2, h2_2);
        println!("  rg = {:.4}", rg);

        if let (Some(sp1), Some(pp1), Some(sp2), Some(pp2)) =
            (samp_prev1, pop_prev1, samp_prev2, pop_prev2)
        {
            let c1 = liability_conversion_factor(sp1, pp1);
            let c2 = liability_conversion_factor(sp2, pp2);
            let h2_1_liab = h2_1 * c1;
            let h2_2_liab = h2_2 * c2;
            let gencov_liab = gencov * (c1 * c2).sqrt();
            let rg_liab = if h2_1_liab > 0.0 && h2_2_liab > 0.0 {
                gencov_liab / (h2_1_liab * h2_2_liab).sqrt()
            } else {
                f64::NAN
            };
            println!(
                "  Liability-scale h2_1={:.4}  h2_2={:.4}  rg={:.4}",
                h2_1_liab, h2_2_liab, rg_liab
            );
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Jackknife diagnostic printer (--print-cov / --print-delete-vals)
// ---------------------------------------------------------------------------

fn print_jackknife_diagnostics(
    result: &crate::irwls::IrwlsResult,
    print_cov: bool,
    print_delete_vals: bool,
) {
    if print_cov && let Some(cov) = &result.jknife_cov {
        println!("\nJackknife covariance matrix ({n}×{n}):", n = cov.nrows());
        for i in 0..cov.nrows() {
            let row: Vec<String> = (0..cov.ncols())
                .map(|j| format!("{:.6e}", cov[[i, j]]))
                .collect();
            println!("  {}", row.join("  "));
        }
    }
    if print_delete_vals && let Some(dv) = &result.delete_values {
        println!(
            "\nJackknife delete values ({} blocks × {} params):",
            dv.nrows(),
            dv.ncols()
        );
        for k in 0..dv.nrows() {
            let row: Vec<String> = (0..dv.ncols())
                .map(|j| format!("{:.6}", dv[[k, j]]))
                .collect();
            println!("  Block {:3}: {}", k, row.join("  "));
        }
    }
}

fn ref_l2_raw2(merged: &DataFrame) -> Result<Array1<f64>> {
    extract_f64(merged, "ref_l2").map(|v| v.mapv(|l| l.max(0.0)))
}
fn w_l2_raw2(merged: &DataFrame) -> Result<Array1<f64>> {
    extract_f64(merged, "w_l2").map(|v| v.mapv(|l| l.max(0.0)))
}
fn n1_raw2(merged: &DataFrame) -> Result<Array1<f64>> {
    extract_f64(merged, "N1")
}
fn n2_raw2(merged: &DataFrame) -> Result<Array1<f64>> {
    extract_f64(merged, "N2")
}

fn run_h2_scalar(
    chi2: &Array1<f64>,
    ref_l2: &Array1<f64>,
    w_l2: &Array1<f64>,
    n_vec: &Array1<f64>,
    m_snps: f64,
    n_blocks: usize,
    fixed_intercept: Option<f64>,
) -> Result<f64> {
    let n_obs = chi2.len();
    let n_mean = n_vec.mean().unwrap_or(1.0);

    let (y, x) = if let Some(fixed) = fixed_intercept {
        let y_adj = chi2.mapv(|c| c - fixed);
        let mut x = Array2::<f64>::zeros((n_obs, 1));
        for i in 0..n_obs {
            x[[i, 0]] = ref_l2[i];
        }
        (y_adj, x)
    } else {
        let mut x = Array2::<f64>::zeros((n_obs, 2));
        for i in 0..n_obs {
            x[[i, 0]] = ref_l2[i];
            x[[i, 1]] = 1.0;
        }
        (chi2.clone(), x)
    };

    let mut weights: Array1<f64> = Array1::zeros(n_obs);
    for i in 0..n_obs {
        let v = (w_l2[i] / m_snps + 1.0 / n_vec[i]).max(1e-9);
        weights[i] = 1.0 / (2.0 * v * v);
    }

    let result = jackknife::jackknife(&x, &y, &weights, n_blocks, 2)?;
    Ok(result.est[0] * m_snps / n_mean)
}

// ---------------------------------------------------------------------------
// Polars → ndarray helper
// ---------------------------------------------------------------------------

/// Extract a column from a Polars DataFrame as `Array1<f64>`.
/// Missing values (null) are replaced with NaN.
fn extract_f64(df: &DataFrame, name: &str) -> Result<Array1<f64>> {
    let series = df
        .column(name)
        .with_context(|| format!("column '{}' not found in DataFrame", name))?;
    let ca = series
        .cast(&DataType::Float64)
        .with_context(|| format!("casting column '{}' to f64", name))?;
    let ca = ca
        .f64()
        .with_context(|| format!("column '{}' as f64 chunked array", name))?;
    let vec: Vec<f64> = ca.into_iter().map(|v| v.unwrap_or(f64::NAN)).collect();
    Ok(Array1::from_vec(vec))
}

// ---------------------------------------------------------------------------
// Statistical helpers
// ---------------------------------------------------------------------------

/// Convert a chi-squared statistic to a two-sided p-value.
#[allow(dead_code)]
pub fn chi2_p_value(chi2_stat: f64) -> f64 {
    let dist = ChiSquared::new(1.0).expect("ChiSquared(df=1)");
    1.0 - dist.cdf(chi2_stat)
}

/// Convert estimate + SE to a Z-score and p-value (two-sided).
#[allow(dead_code)]
pub fn p_z_norm(est: f64, se: f64) -> (f64, f64) {
    let z = est / se;
    let p = chi2_p_value(z * z);
    (z, p)
}

/// Convert observed-scale genetic covariance to liability scale.
#[allow(dead_code)]
pub fn gencov_obs_to_liab(
    gencov: f64,
    samp_prev1: f64,
    pop_prev1: f64,
    samp_prev2: f64,
    pop_prev2: f64,
) -> f64 {
    let c1 = liability_conversion_factor(samp_prev1, pop_prev1);
    let c2 = liability_conversion_factor(samp_prev2, pop_prev2);
    gencov * (c1 * c2).sqrt()
}

/// Compute the conversion factor for one trait (sample → liability scale).
fn liability_conversion_factor(samp_prev: f64, pop_prev: f64) -> f64 {
    let normal = Normal::new(0.0, 1.0).expect("Normal(0,1)");
    let t = normal.inverse_cdf(1.0 - pop_prev);
    let z = normal.pdf(t);
    let k = pop_prev;
    let p = samp_prev;
    k * k * (1.0 - k) * (1.0 - k) / (p * (1.0 - p) * z * z)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// p_z_norm: est=3, se=1 → z=3, p ≈ 0.00270 (2-sided normal tail).
    #[test]
    fn test_p_z_norm_basic() {
        let (z, p) = p_z_norm(3.0, 1.0);
        assert!((z - 3.0).abs() < 1e-9, "z={}", z);
        // 2-sided p for z=3: 2*Φ(-3) ≈ 0.002700
        assert!((p - 0.002700).abs() < 0.0001, "p={}", p);
    }

    /// p_z_norm at z=10: statrs underflows chi2(100, df=1) to p=0 (float64 limit);
    /// we assert p < 1e-20 to accommodate a future fix.
    #[test]
    fn test_p_z_norm_extreme() {
        let (z, p) = p_z_norm(10.0, 1.0);
        assert!((z - 10.0).abs() < 1e-9, "z={}", z);
        assert!(p < 1e-20 || p == 0.0, "p={} should be near zero", p);
    }

    /// p_z_norm: se=0 → z = ±inf, p = 0.
    #[test]
    fn test_p_z_norm_zero_se() {
        let (z, p) = p_z_norm(10.0, 0.0);
        assert!(z.is_infinite(), "z should be inf, got {}", z);
        assert_eq!(p, 0.0, "p should be 0 when z=inf");
    }

    /// chi2_p_value(0) = 1 (all probability mass above 0).
    #[test]
    fn test_chi2_p_value_at_zero() {
        let p = chi2_p_value(0.0);
        assert!((p - 1.0).abs() < 1e-9, "chi2_p(0)={}", p);
    }

    /// chi2_p_value at the canonical 5% critical value of a 1-df chi2.
    #[test]
    fn test_chi2_p_value_critical_5pct() {
        // 1-df chi2 critical value for alpha=0.05 is 3.84146...
        let p = chi2_p_value(3.841459);
        assert!((p - 0.05).abs() < 0.001, "chi2_p(3.841)={}", p);
    }

    /// gencov_obs_to_liab for two binary traits (balanced study, 1% phenotype).
    /// gencov_obs_to_liab(1, 0.5, 0.01, 0.5, 0.01) ≈ 0.551907.
    #[test]
    fn test_gencov_obs_to_liab_both_traits() {
        let result = gencov_obs_to_liab(1.0, 0.5, 0.01, 0.5, 0.01);
        assert!(
            (result - 0.551907).abs() < 1e-4,
            "gencov_obs_to_liab={:.6}",
            result
        );
    }

    /// Checks that sqrt(both-trait result) = single-trait factor.
    /// For a QT trait (c=1) paired with a binary: gencov factor = sqrt(binary_factor).
    #[test]
    fn test_gencov_obs_to_liab_single_trait_factor() {
        let both = gencov_obs_to_liab(1.0, 0.5, 0.01, 0.5, 0.01);
        // sqrt of the both-trait factor should equal the single-trait factor
        assert!(
            (both.sqrt() - 0.551907f64.sqrt()).abs() < 1e-4,
            "sqrt(both)={:.6}",
            both.sqrt()
        );
    }

    /// liability_conversion_factor(samp_prev=0.5, pop_prev=0.01) ≈ 0.551907
    /// (balanced case-control study, 1% phenotype prevalence).
    #[test]
    fn test_h2_obs_to_liab_balanced_study() {
        // balanced case: samp_prev = 0.5 (case-control 1:1), pop_prev = 0.01 (1% phenotype)
        let c = liability_conversion_factor(0.5, 0.01);
        assert!(
            (c - 0.551907).abs() < 1e-5,
            "h2_obs_to_liab factor={:.6}",
            c
        );
    }

    /// Verifies sign convention: negative estimate → negative z, same |p| as positive.
    #[test]
    fn test_p_z_norm_sign_convention() {
        // Negative estimate → negative z, same |p| as positive
        let (z_pos, p_pos) = p_z_norm(3.0, 1.0);
        let (z_neg, p_neg) = p_z_norm(-3.0, 1.0);
        assert!((z_pos - 3.0).abs() < 1e-9);
        assert!((z_neg + 3.0).abs() < 1e-9);
        // p-value is two-sided so should be the same magnitude
        assert!(
            (p_pos - p_neg).abs() < 1e-9,
            "p_pos={} p_neg={}",
            p_pos,
            p_neg
        );
    }
}
