/// Heritability and genetic correlation estimation — replaces ldscore/regressions.py.
///
/// Pipeline:
///   parse → merge (sumstats ⨝ ld_scores ⨝ weights) → build ndarray
///   matrices → jackknife(IRWLS, n_blocks) → h2 / rg estimates.
use anyhow::{Context, Result};
use ndarray::{Array1, Array2, Axis, s};
use polars::prelude::*;
use statrs::distribution::StudentsT;
use statrs::distribution::{ChiSquared, Continuous, ContinuousCDF, Normal};
use std::fs::File;
use std::io::{BufWriter, Write};

use crate::cli::{H2Args, RgArgs};
use crate::irwls::IrwlsResult;
use crate::jackknife;
use crate::parse;

// ---------------------------------------------------------------------------
// Helper: load LD scores from either a single file or per-chr split files
// ---------------------------------------------------------------------------

/// Load a scalar (single L2 column) LD score file, aliasing the L2 column.
/// Used for weight LD scores (always scalar). Comma-separated lists are not allowed.
fn load_ld(single: Option<&str>, chr_prefix: Option<&str>, alias: &str) -> Result<LazyFrame> {
    if single.map(|s| s.contains(',')).unwrap_or(false)
        || chr_prefix.map(|s| s.contains(',')).unwrap_or(false)
    {
        anyhow::bail!("--w-ld/--w-ld-chr must point to a single fileset (no commas allowed)");
    }
    match (single, chr_prefix) {
        (Some(path), _) => {
            Ok(parse::scan_ldscore(path)?.select([col("SNP"), col("L2").alias(alias)]))
        }
        (None, Some(prefix)) => Ok(parse::concat_chrs_any(
            prefix,
            &[".l2.ldscore.gz", ".l2.ldscore.bz2", ".l2.ldscore"],
        )?
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
fn split_paths(raw: &str) -> Vec<String> {
    raw.split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect()
}

fn load_ld_ref_single_path(path: &str) -> Result<LazyFrame> {
    let lf = parse::scan_ldscore(path)?;
    drop_ld_meta(lf)
}

fn load_ld_ref_single_chr(prefix: &str) -> Result<LazyFrame> {
    let lf = parse::concat_chrs_any(
        prefix,
        &[".l2.ldscore.gz", ".l2.ldscore.bz2", ".l2.ldscore"],
    )?;
    drop_ld_meta(lf)
}

fn drop_ld_meta(lf: LazyFrame) -> Result<LazyFrame> {
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

fn load_ld_ref(single: Option<&str>, chr_prefix: Option<&str>) -> Result<LazyFrame> {
    match (single, chr_prefix) {
        (Some(path), None) => {
            let paths = split_paths(path);
            if paths.len() == 1 {
                load_ld_ref_single_path(&paths[0])
            } else {
                load_ld_ref_multi_paths(&paths, false)
            }
        }
        (None, Some(prefix)) => {
            let prefixes = split_paths(prefix);
            if prefixes.len() == 1 {
                load_ld_ref_single_chr(&prefixes[0])
            } else {
                load_ld_ref_multi_paths(&prefixes, true)
            }
        }
        (None, None) => anyhow::bail!("Must specify exactly one of --ref-ld / --ref-ld-chr"),
        _ => anyhow::bail!("Cannot set both --ref-ld and --ref-ld-chr"),
    }
}

fn load_ld_ref_multi_paths(paths: &[String], chr_split: bool) -> Result<LazyFrame> {
    anyhow::ensure!(!paths.is_empty(), "Empty LD score list");
    let mut dfs: Vec<DataFrame> = Vec::new();
    for p in paths {
        let lf = if chr_split {
            load_ld_ref_single_chr(p)?
        } else {
            load_ld_ref_single_path(p)?
        };
        dfs.push(
            lf.collect()
                .with_context(|| format!("loading LD scores '{}'", p))?,
        );
    }
    let mut base = dfs.remove(0);
    // Match Python: suffix columns from the first file with _0.
    let base_cols: Vec<String> = base
        .get_column_names()
        .iter()
        .filter(|c| c.as_str() != "SNP")
        .map(|c| c.to_string())
        .collect();
    for c in base_cols {
        let new_name = format!("{}{}", c, "_0");
        base.rename(c.as_str(), new_name.into())?;
    }
    let base_snps: Vec<Option<String>> = base
        .column("SNP")?
        .as_materialized_series()
        .str()?
        .into_iter()
        .map(|opt| opt.map(|s| s.to_string()))
        .collect();
    for (idx, df) in dfs.into_iter().enumerate() {
        let snps: Vec<Option<String>> = df
            .column("SNP")?
            .as_materialized_series()
            .str()?
            .into_iter()
            .map(|opt| opt.map(|s| s.to_string()))
            .collect();
        anyhow::ensure!(
            snps == base_snps,
            "LD Scores for concatenation must have identical SNP columns"
        );
        let mut df = df.clone();
        df.drop_in_place("SNP")?;
        // Match Python: add _i suffix to each column from the i-th file.
        let suffix = format!("_{}", idx + 1);
        let cols: Vec<String> = df
            .get_column_names()
            .iter()
            .map(|c| c.to_string())
            .collect();
        for c in cols {
            let new_name = format!("{}{}", c, suffix);
            df.rename(c.as_str(), new_name.into())?;
        }
        base.hstack_mut(df.columns())?;
    }
    Ok(base.lazy())
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
        let prefixes = split_paths(prefix);
        let m_vec_res = if prefixes.len() == 1 {
            parse::read_m_vec(prefixes[0].as_str(), suffix)
        } else {
            parse::read_m_vec_list(&prefixes, suffix)
        };
        match m_vec_res {
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

fn ensure_sumstats_have_alleles(path: &str) -> Result<()> {
    let lf = parse::scan_sumstats(path)?;
    let header = lf
        .limit(0)
        .collect()
        .with_context(|| format!("peeking sumstats header '{}'", path))?;
    let cols = header.get_column_names();
    let has = |n: &str| cols.iter().any(|c| *c == n);
    anyhow::ensure!(
        has("A1") && has("A2"),
        "Sumstats file '{}' is missing A1/A2 columns (required unless --no-check-alleles)",
        path
    );
    Ok(())
}

fn is_valid_snp(a1: &str, a2: &str) -> bool {
    if a1.len() != 1 || a2.len() != 1 {
        return false;
    }
    let a1 = a1.chars().next().unwrap();
    let a2 = a2.chars().next().unwrap();
    let valid = |c| matches!(c, 'A' | 'C' | 'G' | 'T');
    if !valid(a1) || !valid(a2) || a1 == a2 {
        return false;
    }
    // strand-ambiguous
    !matches!((a1, a2), ('A', 'T') | ('T', 'A') | ('C', 'G') | ('G', 'C'))
}

fn complement_base(c: char) -> char {
    match c {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => c,
    }
}

fn match_and_flip(a1: char, a2: char, b1: char, b2: char) -> Option<bool> {
    let b1c = complement_base(b1);
    let b2c = complement_base(b2);
    if (a1 == b1 && a2 == b2) || (a1 == b1c && a2 == b2c) {
        return Some(false);
    }
    if (a1 == b2 && a2 == b1) || (a1 == b2c && a2 == b1c) {
        return Some(true);
    }
    None
}

fn align_rg_alleles(merged: &DataFrame, z2: &Array1<f64>) -> Result<(Vec<bool>, Vec<f64>, usize)> {
    let a1_1 = merged.column("A1_1")?.as_materialized_series().str()?;
    let a2_1 = merged.column("A2_1")?.as_materialized_series().str()?;
    let a1_2 = merged.column("A1_2")?.as_materialized_series().str()?;
    let a2_2 = merged.column("A2_2")?.as_materialized_series().str()?;

    let mut mask: Vec<bool> = Vec::with_capacity(merged.height());
    let mut z2_aligned: Vec<f64> = Vec::new();
    let mut removed = 0usize;

    for i in 0..merged.height() {
        let a1 = a1_1.get(i).unwrap_or("").to_ascii_uppercase();
        let a2 = a2_1.get(i).unwrap_or("").to_ascii_uppercase();
        let b1 = a1_2.get(i).unwrap_or("").to_ascii_uppercase();
        let b2 = a2_2.get(i).unwrap_or("").to_ascii_uppercase();

        if !is_valid_snp(&a1, &a2) || !is_valid_snp(&b1, &b2) {
            mask.push(false);
            removed += 1;
            continue;
        }
        let a1c = a1.chars().next().unwrap();
        let a2c = a2.chars().next().unwrap();
        let b1c = b1.chars().next().unwrap();
        let b2c = b2.chars().next().unwrap();

        match match_and_flip(a1c, a2c, b1c, b2c) {
            Some(flip) => {
                mask.push(true);
                let z = z2[i];
                z2_aligned.push(if flip { -z } else { z });
            }
            None => {
                mask.push(false);
                removed += 1;
            }
        }
    }

    Ok((mask, z2_aligned, removed))
}

// ---------------------------------------------------------------------------
// Heritability (--h2)
// ---------------------------------------------------------------------------

pub fn run_h2(args: H2Args) -> Result<()> {
    if args.h2_cts.is_some() {
        return run_h2_cts(&args);
    }

    if args.overlap_annot {
        if args.not_m_5_50 {
            if args.frqfile.is_some() || args.frqfile_chr.is_some() {
                println!("  WARNING: --frqfile is ignored when --not-m-5-50 is set");
            }
        } else {
            let ok = (args.frqfile.is_some() && args.ref_ld.is_some())
                || (args.frqfile_chr.is_some() && args.ref_ld_chr.is_some());
            anyhow::ensure!(
                ok,
                "Must set either --frqfile and --ref-ld or --frqfile-chr and --ref-ld-chr"
            );
        }
    } else if args.frqfile.is_some() || args.frqfile_chr.is_some() {
        println!("  WARNING: --frqfile is ignored without --overlap-annot");
    }

    let h2_path = args
        .h2
        .as_ref()
        .context("--h2 is required unless --h2-cts is set")?;
    let sumstats_lf = parse::scan_sumstats(h2_path)?;
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
        n_total, h2_path
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
        let fit = run_h2_partitioned(
            &chi2, &ref_l2_k, &w_l2, &n_vec, &l2_cols, n_obs, n_mean, &args,
        )?;

        if args.overlap_annot {
            let chr_split = args.ref_ld_chr.is_some();
            let prefixes = if let Some(prefix) = args.ref_ld_chr.as_deref() {
                split_paths(prefix)
            } else {
                split_paths(args.ref_ld.as_deref().unwrap_or_default())
            };
            let (overlap, m_tot, overlap_names) = parse::read_overlap_matrix(
                &prefixes,
                if args.not_m_5_50 {
                    None
                } else {
                    args.frqfile.as_deref()
                },
                if args.not_m_5_50 {
                    None
                } else {
                    args.frqfile_chr.as_deref()
                },
                chr_split,
            )?;
            anyhow::ensure!(
                overlap_names.len() == l2_cols.len(),
                "Overlap annotations (K={}) do not match LD score columns (K={})",
                overlap_names.len(),
                l2_cols.len()
            );
            write_overlap_results(
                &fit,
                &l2_cols,
                &overlap,
                m_tot,
                args.print_coefficients,
                &args.out,
            )?;
        }
        return Ok(());
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

fn read_ldcts(path: &str) -> Result<Vec<(String, String)>> {
    let content =
        std::fs::read_to_string(path).with_context(|| format!("reading ldcts '{}'", path))?;
    let mut out = Vec::new();
    for (i, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        anyhow::ensure!(
            parts.len() == 2,
            "ldcts line {} must have 2 columns (name + prefixes)",
            i + 1
        );
        out.push((parts[0].to_string(), parts[1].to_string()));
    }
    Ok(out)
}

fn run_h2_cts(args: &H2Args) -> Result<()> {
    let sumstats_path = args
        .h2_cts
        .as_ref()
        .context("--h2-cts requires a sumstats file")?;
    let ldcts_path = args
        .ref_ld_chr_cts
        .as_ref()
        .context("--h2-cts requires --ref-ld-chr-cts")?;
    anyhow::ensure!(args.ref_ld_chr.is_some(), "--h2-cts requires --ref-ld-chr");
    anyhow::ensure!(
        args.ref_ld.is_none(),
        "--h2-cts requires --ref-ld-chr (not --ref-ld)"
    );
    anyhow::ensure!(args.w_ld_chr.is_some(), "--h2-cts requires --w-ld-chr");
    anyhow::ensure!(
        args.w_ld.is_none(),
        "--h2-cts requires --w-ld-chr (not --w-ld)"
    );

    let sumstats_lf = parse::scan_sumstats(sumstats_path)?;
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
        n_total, sumstats_path
    );

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
    let k_base = l2_cols.len();
    anyhow::ensure!(
        k_base > 0,
        "No L2 annotation columns found in merged dataset"
    );

    let chi2_raw = extract_f64(&merged, "Z")?.mapv(|z| (z * z).max(0.0));
    let w_l2_raw = extract_f64(&merged, "w_l2")?.mapv(|l| l.max(0.0));
    let n_vec_raw = extract_f64(&merged, "N")?;
    let ref_l2_raw_k: Vec<Array1<f64>> = l2_cols
        .iter()
        .map(|name| extract_f64(&merged, name).map(|v| v.mapv(|l| l.max(0.0))))
        .collect::<Result<Vec<_>>>()?;

    let chisq_max = if let Some(c) = args.chisq_max {
        Some(c)
    } else {
        let max_n = n_vec_raw.iter().cloned().fold(0.0f64, f64::max);
        Some((0.001 * max_n).max(80.0))
    };

    let (chi2, w_l2, n_vec, ref_l2_k, keep_snps) = if let Some(cmax) = chisq_max {
        let mask: Vec<bool> = chi2_raw.iter().map(|&c| c <= cmax).collect();
        let n_removed = mask.iter().filter(|&&b| !b).count();
        if n_removed > 0 {
            println!("  Removed {} SNPs with chi2 > {:.1}", n_removed, cmax);
        }
        let keep_snps_all: Vec<String> = merged
            .column("SNP")?
            .as_materialized_series()
            .str()?
            .into_iter()
            .map(|opt| opt.unwrap_or("").to_string())
            .collect();
        let keep_snps: Vec<String> = keep_snps_all
            .into_iter()
            .zip(mask.iter())
            .filter_map(|(s, keep)| if *keep { Some(s) } else { None })
            .collect();
        (
            filter_by_mask(&chi2_raw, &mask),
            filter_by_mask(&w_l2_raw, &mask),
            filter_by_mask(&n_vec_raw, &mask),
            ref_l2_raw_k
                .iter()
                .map(|v| filter_by_mask(v, &mask))
                .collect::<Vec<_>>(),
            keep_snps,
        )
    } else {
        let keep_snps: Vec<String> = merged
            .column("SNP")?
            .as_materialized_series()
            .str()?
            .into_iter()
            .map(|opt| opt.unwrap_or("").to_string())
            .collect();
        (chi2_raw, w_l2_raw, n_vec_raw, ref_l2_raw_k, keep_snps)
    };

    let n_obs = chi2.len();
    let n_mean = n_vec.mean().unwrap_or(1.0);

    let m_base = resolve_m_vec(
        args.m_snps,
        args.ref_ld_chr.as_deref(),
        args.not_m_5_50,
        n_obs,
        k_base,
    );
    anyhow::ensure!(
        m_base.len() == ref_l2_k.len(),
        "Baseline M length ({}) does not match LD score columns ({})",
        m_base.len(),
        ref_l2_k.len()
    );

    let ldcts = read_ldcts(ldcts_path)?;
    let norm = Normal::new(0.0, 1.0)?;
    let mut results: Vec<(String, f64, f64, f64)> = Vec::new();

    for (name, ct_ld_chr) in ldcts {
        let ct_prefixes = split_paths(&ct_ld_chr);
        let suffix = if args.not_m_5_50 {
            ".l2.M"
        } else {
            ".l2.M_5_50"
        };
        let m_cts = parse::read_m_vec_list(&ct_prefixes, suffix)
            .with_context(|| format!("reading M files for '{}'", ct_ld_chr))?;

        let cts_ld = load_ld_ref(None, Some(ct_ld_chr.as_str()))
            .context("loading CTS LD scores")?
            .collect()
            .context("collecting CTS LD scores")?;
        let cts_cols: Vec<String> = cts_ld
            .get_column_names()
            .iter()
            .filter(|c| c.as_str() != "SNP")
            .map(|c| c.to_string())
            .collect();

        let keep_df = DataFrame::new_infer_height(vec![
            Series::new("SNP".into(), keep_snps.clone()).into(),
            Series::new(
                "idx".into(),
                (0..keep_snps.len()).map(|i| i as u32).collect::<Vec<u32>>(),
            )
            .into(),
        ])?;
        let mut joined = keep_df.join(
            &cts_ld,
            ["SNP"],
            ["SNP"],
            JoinArgs::new(JoinType::Left),
            None,
        )?;
        joined = joined.sort(["idx"], SortMultipleOptions::default())?;

        let mut cts_arrays: Vec<Array1<f64>> = Vec::new();
        for col in &cts_cols {
            let s = joined
                .column(col)
                .with_context(|| format!("CTS column '{}'", col))?
                .cast(&DataType::Float64)
                .with_context(|| format!("casting CTS column '{}'", col))?;
            let ca = s
                .f64()
                .with_context(|| format!("CTS column '{}' as f64", col))?;
            let mut v = Vec::with_capacity(n_obs);
            for val in ca.into_iter() {
                let val = val.ok_or_else(|| {
                    anyhow::anyhow!(
                        "Missing some LD scores from cts files. Are you sure all SNPs in ref-ld-chr are also in ref-ld-chr-cts?"
                    )
                })?;
                v.push(val.max(0.0));
            }
            cts_arrays.push(Array1::from(v));
        }

        anyhow::ensure!(
            m_cts.len() == cts_arrays.len(),
            "M vector length ({}) does not match CTS LD columns ({}) for '{}'",
            m_cts.len(),
            cts_arrays.len(),
            ct_ld_chr
        );

        let mut ref_l2_k_final = Vec::new();
        ref_l2_k_final.extend(cts_arrays);
        ref_l2_k_final.extend(ref_l2_k.clone());

        let mut m_vec = Vec::new();
        m_vec.extend(m_cts);
        m_vec.extend(m_base.clone());

        let fit = fit_h2_partitioned(
            &chi2,
            &ref_l2_k_final,
            &w_l2,
            &n_vec,
            n_obs,
            n_mean,
            args,
            Some(&m_vec),
        )?;

        let coef = fit.est[0] / n_mean;
        let coef_var = fit.jknife_cov[[0, 0]] / (n_mean * n_mean);
        let coef_se = coef_var.max(0.0).sqrt();
        let z = if coef_se == 0.0 {
            f64::NAN
        } else {
            coef / coef_se
        };
        let p = if z.is_nan() { f64::NAN } else { norm.sf(z) };
        results.push((name.clone(), coef, coef_se, p));

        if args.print_all_cts {
            for i in 1..ct_prefixes.len() {
                if i >= fit.est.len() {
                    break;
                }
                let coef = fit.est[i] / n_mean;
                let coef_var = fit.jknife_cov[[i, i]] / (n_mean * n_mean);
                let coef_se = coef_var.max(0.0).sqrt();
                let z = if coef_se == 0.0 {
                    f64::NAN
                } else {
                    coef / coef_se
                };
                let p = if z.is_nan() { f64::NAN } else { norm.sf(z) };
                results.push((format!("{}_{}", name, i), coef, coef_se, p));
            }
        }
    }

    results.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal));

    let out_path = format!("{}.cell_type_results.txt", args.out);
    let mut w = BufWriter::new(
        File::create(&out_path)
            .with_context(|| format!("creating CTS results file '{}'", out_path))?,
    );
    writeln!(
        w,
        "Name\tCoefficient\tCoefficient_std_error\tCoefficient_P_value"
    )?;
    for (name, coef, se, p) in results {
        writeln!(w, "{}\t{:.6e}\t{:.6e}\t{:.6e}", name, coef, se, p)?;
    }
    println!("Results printed to {}", out_path);
    Ok(())
}

// ---------------------------------------------------------------------------
// Partitioned h2 regression (K annotation columns)
// ---------------------------------------------------------------------------

struct PartitionedFit {
    est: Array1<f64>,
    jknife_cov: Array2<f64>,
    delete_values: Array2<f64>,
    m_vec: Vec<f64>,
    n_mean: f64,
    h2_per_annot: Vec<f64>,
    h2_total: f64,
    intercept: f64,
    n_blocks: usize,
}

fn ratio_jackknife(
    est: &Array1<f64>,
    numer_delete: &Array2<f64>,
    denom_delete: &Array2<f64>,
) -> (Array2<f64>, Array1<f64>) {
    let n_blocks = numer_delete.nrows();
    let k = numer_delete.ncols();
    let mut pseudo = Array2::<f64>::zeros((n_blocks, k));
    for b in 0..n_blocks {
        for j in 0..k {
            let denom = denom_delete[[b, j]];
            let numer = numer_delete[[b, j]];
            pseudo[[b, j]] = if denom == 0.0 {
                f64::NAN
            } else {
                n_blocks as f64 * est[j] - (n_blocks as f64 - 1.0) * (numer / denom)
            };
        }
    }
    let mean = pseudo.mean_axis(Axis(0)).unwrap();
    let centered = Array2::from_shape_fn((n_blocks, k), |(b, j)| pseudo[[b, j]] - mean[j]);
    let cov = centered.t().dot(&centered) / ((n_blocks as f64 - 1.0) * n_blocks as f64);
    let mut se = Array1::<f64>::zeros(k);
    for j in 0..k {
        se[j] = cov[[j, j]].max(0.0).sqrt();
    }
    (cov, se)
}

#[allow(clippy::too_many_arguments)]
fn fit_h2_partitioned(
    chi2: &Array1<f64>,
    ref_l2_k: &[Array1<f64>],
    w_l2: &Array1<f64>,
    n_vec: &Array1<f64>,
    n_obs: usize,
    n_mean: f64,
    args: &H2Args,
    m_override: Option<&[f64]>,
) -> Result<PartitionedFit> {
    let k = ref_l2_k.len();
    let m_vec = if let Some(m) = m_override {
        anyhow::ensure!(
            m.len() == k,
            "M vector length {} does not match K={}",
            m.len(),
            k
        );
        m.to_vec()
    } else {
        resolve_m_vec(
            args.m_snps,
            args.ref_ld_chr.as_deref(),
            args.not_m_5_50,
            n_obs,
            k,
        )
    };
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

    let result: IrwlsResult = jackknife::jackknife(&x_reg, &y_reg, &weights, args.n_blocks, 2)
        .context("jackknife IRWLS for partitioned h2")?;

    let h2_per_annot: Vec<f64> = (0..k).map(|j| result.est[j] * m_vec[j] / n_mean).collect();
    let h2_total: f64 = h2_per_annot.iter().sum();
    let intercept = fixed_intercept.unwrap_or_else(|| result.est[k]);

    let jknife_cov = result
        .jknife_cov
        .clone()
        .context("missing jackknife covariance")?;
    let delete_values = result
        .delete_values
        .clone()
        .context("missing jackknife delete values")?;

    Ok(PartitionedFit {
        est: result.est,
        jknife_cov,
        delete_values,
        m_vec,
        n_mean,
        h2_per_annot,
        h2_total,
        intercept,
        n_blocks: args.n_blocks,
    })
}

fn print_partitioned_summary(fit: &PartitionedFit, l2_cols: &[String], args: &H2Args) {
    println!("Total observed-scale h2: {:.4}", fit.h2_total);
    println!("Intercept: {:.4}", fit.intercept);

    if args.print_coefficients && !args.overlap_annot {
        let m_total: f64 = fit.m_vec.iter().sum();
        println!();
        println!(
            "{:<35} {:>10} {:>10} {:>12} {:>14}",
            "Category", "Prop_h2", "Prop_SNPs", "Enrichment", "Coefficient"
        );
        println!("{}", "-".repeat(83));
        for (j, col_name) in l2_cols.iter().enumerate() {
            let prop_h2 = if fit.h2_total.abs() > 1e-12 {
                fit.h2_per_annot[j] / fit.h2_total
            } else {
                f64::NAN
            };
            let prop_m = fit.m_vec[j] / m_total;
            let enrichment = if prop_m > 1e-12 {
                prop_h2 / prop_m
            } else {
                f64::NAN
            };
            let coeff = fit.est[j] / fit.n_mean;
            println!(
                "{:<35} {:>10.4} {:>10.4} {:>12.4} {:>14.4e}",
                col_name, prop_h2, prop_m, enrichment, coeff
            );
        }
        println!();
    }
}

fn write_overlap_results(
    fit: &PartitionedFit,
    category_names: &[String],
    overlap: &Array2<f64>,
    m_tot: usize,
    print_coefficients: bool,
    out_prefix: &str,
) -> Result<()> {
    let k = category_names.len();
    anyhow::ensure!(
        overlap.nrows() == k && overlap.ncols() == k,
        "Overlap matrix shape {}x{} does not match K={}",
        overlap.nrows(),
        overlap.ncols(),
        k
    );

    let m_vec = &fit.m_vec;
    let m_tot_f = m_tot as f64;

    let coef = fit.est.slice(s![..k]).to_owned() / fit.n_mean;
    let coef_cov = fit.jknife_cov.slice(s![..k, ..k]).to_owned() / (fit.n_mean * fit.n_mean);
    let mut coef_se = Array1::<f64>::zeros(k);
    for j in 0..k {
        coef_se[j] = coef_cov[[j, j]].max(0.0).sqrt();
    }

    let prop = Array1::from(
        fit.h2_per_annot
            .iter()
            .map(|&h| {
                if fit.h2_total.abs() > 1e-12 {
                    h / fit.h2_total
                } else {
                    f64::NAN
                }
            })
            .collect::<Vec<_>>(),
    );
    let (prop_cov, _prop_se) = {
        let n_blocks = fit.delete_values.nrows();
        let mut numer = Array2::<f64>::zeros((n_blocks, k));
        for b in 0..n_blocks {
            for j in 0..k {
                numer[[b, j]] = m_vec[j] * fit.delete_values[[b, j]] / fit.n_mean;
            }
        }
        let mut denom = Array2::<f64>::zeros((n_blocks, k));
        for b in 0..n_blocks {
            let sum = numer.row(b).sum();
            for j in 0..k {
                denom[[b, j]] = sum;
            }
        }
        ratio_jackknife(&prop, &numer, &denom)
    };

    let mut overlap_prop = Array2::<f64>::zeros((k, k));
    for i in 0..k {
        for j in 0..k {
            overlap_prop[[i, j]] = overlap[[i, j]] / m_vec[j];
        }
    }

    let prop_h2_overlap = overlap_prop.dot(&prop);
    let prop_h2_overlap_var = overlap_prop.dot(&prop_cov).dot(&overlap_prop.t());
    let mut prop_h2_overlap_se = Array1::<f64>::zeros(k);
    for j in 0..k {
        prop_h2_overlap_se[j] = prop_h2_overlap_var[[j, j]].max(0.0).sqrt();
    }

    let prop_m_overlap = Array1::from(
        m_vec
            .iter()
            .map(|&m| if m_tot_f > 0.0 { m / m_tot_f } else { f64::NAN })
            .collect::<Vec<_>>(),
    );
    let enrichment = &prop_h2_overlap / &prop_m_overlap;
    let enrichment_se = &prop_h2_overlap_se / &prop_m_overlap;

    let mut overlap_diff = Array2::<f64>::zeros((k, k));
    for i in 0..k {
        if (m_tot_f - m_vec[i]).abs() < 1e-12 {
            continue;
        }
        for j in 0..k {
            let term1 = overlap[[i, j]] / m_vec[i];
            let term2 = (m_vec[j] - overlap[[i, j]]) / (m_tot_f - m_vec[i]);
            overlap_diff[[i, j]] = term1 - term2;
        }
    }

    let diff_est = overlap_diff.dot(&coef);
    let diff_cov: Array2<f64> = overlap_diff.dot(&coef_cov).dot(&overlap_diff.t());
    let mut diff_se = Array1::<f64>::zeros(k);
    for j in 0..k {
        diff_se[j] = diff_cov[[j, j]].max(0.0).sqrt();
    }
    let tdist =
        StudentsT::new(0.0, 1.0, fit.n_blocks as f64).context("constructing t distribution")?;
    let diff_p: Vec<String> = (0..k)
        .map(|i| {
            if diff_se[i] == 0.0 {
                "NA".to_string()
            } else {
                let t = (diff_est[i] / diff_se[i]).abs();
                format!("{}", 2.0 * tdist.sf(t))
            }
        })
        .collect();

    let out_path = format!("{}.results", out_prefix);
    let mut w = BufWriter::new(
        File::create(&out_path)
            .with_context(|| format!("creating overlap results file '{}'", out_path))?,
    );

    if print_coefficients {
        writeln!(
            w,
            "Category\tProp._SNPs\tProp._h2\tProp._h2_std_error\tEnrichment\tEnrichment_std_error\tEnrichment_p\tCoefficient\tCoefficient_std_error\tCoefficient_z-score"
        )?;
    } else {
        writeln!(
            w,
            "Category\tProp._SNPs\tProp._h2\tProp._h2_std_error\tEnrichment\tEnrichment_std_error\tEnrichment_p"
        )?;
    }

    for i in 0..k {
        if print_coefficients {
            let z = if coef_se[i] == 0.0 {
                f64::NAN
            } else {
                coef[i] / coef_se[i]
            };
            writeln!(
                w,
                "{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{:.6e}\t{:.6e}\t{:.6}",
                category_names[i],
                prop_m_overlap[i],
                prop_h2_overlap[i],
                prop_h2_overlap_se[i],
                enrichment[i],
                enrichment_se[i],
                diff_p[i],
                coef[i],
                coef_se[i],
                z
            )?;
        } else {
            writeln!(
                w,
                "{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}",
                category_names[i],
                prop_m_overlap[i],
                prop_h2_overlap[i],
                prop_h2_overlap_se[i],
                enrichment[i],
                enrichment_se[i],
                diff_p[i]
            )?;
        }
    }

    println!("Results printed to {}", out_path);
    Ok(())
}

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
) -> Result<PartitionedFit> {
    let fit = fit_h2_partitioned(chi2, ref_l2_k, w_l2, n_vec, n_obs, n_mean, args, None)?;

    print_partitioned_summary(&fit, l2_cols, args);
    print_jackknife_diagnostics(
        &IrwlsResult {
            est: fit.est.clone(),
            jknife_se: None,
            jknife_var: None,
            jknife_cov: Some(fit.jknife_cov.clone()),
            delete_values: Some(fit.delete_values.clone()),
        },
        args.print_cov,
        args.print_delete_vals,
    );

    if let (Some(samp_prev), Some(pop_prev)) = (args.samp_prev, args.pop_prev) {
        let c = liability_conversion_factor(samp_prev, pop_prev);
        println!("Liability-scale h2: {:.4}", fit.h2_total * c);
    }

    Ok(fit)
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
        let prefixes = split_paths(prefix);
        let m_res = if prefixes.len() == 1 {
            parse::read_m_total(prefixes[0].as_str(), suffix)
        } else {
            parse::read_m_total_list(&prefixes, suffix)
        };
        match m_res {
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

        if !args.no_check_alleles {
            ensure_sumstats_have_alleles(file1)?;
            ensure_sumstats_have_alleles(file2)?;
        }

        // Load both sumstats and suffix their Z/N columns to avoid conflicts.
        let ss1 = if args.no_check_alleles {
            parse::scan_sumstats(file1)?.select([
                col("SNP"),
                col("Z").alias("Z1"),
                col("N").alias("N1"),
            ])
        } else {
            parse::scan_sumstats(file1)?.select([
                col("SNP"),
                col("A1").alias("A1_1"),
                col("A2").alias("A2_1"),
                col("Z").alias("Z1"),
                col("N").alias("N1"),
            ])
        }
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

        let ss2 = if args.no_check_alleles {
            parse::scan_sumstats(file2)?.select([
                col("SNP"),
                col("Z").alias("Z2"),
                col("N").alias("N2"),
            ])
        } else {
            parse::scan_sumstats(file2)?.select([
                col("SNP"),
                col("A1").alias("A1_2"),
                col("A2").alias("A2_2"),
                col("Z").alias("Z2"),
                col("N").alias("N2"),
            ])
        };

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
        let ref_l2_raw = extract_f64(&merged, "ref_l2")?.mapv(|l| l.max(0.0));
        let w_l2_raw = extract_f64(&merged, "w_l2")?.mapv(|l| l.max(0.0));
        let n1_raw = extract_f64(&merged, "N1")?;
        let n2_raw = extract_f64(&merged, "N2")?;

        let (z1_raw, z2_raw, ref_l2_raw, w_l2_raw, n1_raw, n2_raw) = if args.no_check_alleles {
            (z1_raw, z2_raw, ref_l2_raw, w_l2_raw, n1_raw, n2_raw)
        } else {
            let (mask, z2_aligned, n_removed) =
                align_rg_alleles(&merged, &z2_raw).context("aligning alleles")?;
            if n_removed > 0 {
                println!("  Removed {} SNPs with incompatible alleles", n_removed);
            }
            (
                filter_by_mask(&z1_raw, &mask),
                Array1::from_vec(z2_aligned),
                filter_by_mask(&ref_l2_raw, &mask),
                filter_by_mask(&w_l2_raw, &mask),
                filter_by_mask(&n1_raw, &mask),
                filter_by_mask(&n2_raw, &mask),
            )
        };

        let prod_raw = &z1_raw * &z2_raw;

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
    use polars::prelude::df;

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

    #[test]
    fn test_load_ld_ref_multi_paths_suffixes() {
        let dir = tempfile::tempdir().unwrap();
        let f1 = dir.path().join("ld1.l2.ldscore");
        let f2 = dir.path().join("ld2.l2.ldscore");
        std::fs::write(&f1, "SNP\tL2\nrs1\t1.0\nrs2\t2.0\n").unwrap();
        std::fs::write(&f2, "SNP\tL2\nrs1\t3.0\nrs2\t4.0\n").unwrap();

        let ref_ld = load_ld_ref(
            Some(&format!(
                "{},{}",
                f1.to_str().unwrap(),
                f2.to_str().unwrap()
            )),
            None,
        )
        .unwrap();
        let df = ref_ld.collect().unwrap();
        let cols: Vec<String> = df
            .get_column_names()
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert!(cols.contains(&"SNP".to_string()));
        assert!(cols.contains(&"L2_0".to_string()));
        assert!(cols.contains(&"L2_1".to_string()));
        assert_eq!(df.height(), 2);
    }

    #[test]
    fn test_align_rg_alleles_flip() {
        let df = df![
            "A1_1" => &["A", "A"],
            "A2_1" => &["G", "G"],
            "A1_2" => &["A", "G"],
            "A2_2" => &["G", "A"],
        ]
        .unwrap();
        let z2 = Array1::from_vec(vec![1.0, 2.0]);
        let (mask, z2_aligned, removed) = align_rg_alleles(&df, &z2).unwrap();
        assert_eq!(removed, 0);
        assert_eq!(mask, vec![true, true]);
        assert_eq!(z2_aligned, vec![1.0, -2.0]);
    }
}
