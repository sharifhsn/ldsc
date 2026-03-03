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
use std::time::Instant;
use tracing::{debug, trace};

use crate::cli::{H2Args, RgArgs};
use crate::h2::{
    JackknifeResult, aggregate, combine_twostep, irwls_ldsc, jackknife_fast, ldsc_weights,
    run_h2_ldsc, update_separators, weight_xy,
};
use crate::irwls::IrwlsResult;
use crate::parse;

const MIN_SNPS_WARN: usize = 200_000;

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

fn warn_small_snp_count(n_snps: usize) {
    if n_snps < MIN_SNPS_WARN {
        println!("WARNING: number of SNPs less than 200k; this is almost always bad.");
    }
}

fn smart_merge_on_snp(mut left: DataFrame, mut right: DataFrame) -> Result<DataFrame> {
    let same_order = if left.height() == right.height() {
        let left_snp = left.column("SNP")?;
        let right_snp = right.column("SNP")?;
        let left_snp = if matches!(left_snp.dtype(), DataType::String) {
            left_snp.clone()
        } else {
            left_snp.cast(&DataType::String)?
        };
        let right_snp = if matches!(right_snp.dtype(), DataType::String) {
            right_snp.clone()
        } else {
            right_snp.cast(&DataType::String)?
        };
        left_snp.equals_missing(&right_snp)
    } else {
        false
    };

    if same_order {
        right.drop_in_place("SNP")?;
        let right_cols = right.columns().to_vec();
        left.hstack_mut(&right_cols)?;
        Ok(left)
    } else {
        let mut args = JoinArgs::new(JoinType::Inner);
        args.maintain_order = MaintainOrderJoin::Left;
        Ok(left.join(&right, ["SNP"], ["SNP"], args, None)?)
    }
}

fn trace_array_stats(label: &str, arr: &Array1<f64>) {
    let mut count = 0usize;
    let mut sum = 0.0f64;
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for &v in arr.iter() {
        if v.is_nan() {
            continue;
        }
        count += 1;
        sum += v;
        if v < min {
            min = v;
        }
        if v > max {
            max = v;
        }
    }
    if count == 0 {
        trace!("{label}: empty or all-NaN");
        return;
    }
    let mean = sum / count as f64;
    trace!("{label}: n={count} mean={mean:.6} min={min:.6} max={max:.6}");
}

fn trace_ref_l2_stats(labels: &[String], cols: &[Array1<f64>]) {
    if labels.is_empty() || cols.is_empty() {
        return;
    }
    let k = labels.len().min(cols.len());
    let max_log = 5usize;
    for i in 0..k.min(max_log) {
        trace_array_stats(&format!("ref_l2[{}]", labels[i]), &cols[i]);
    }
    if k > max_log {
        trace!("ref_l2: {} columns total (showing first {})", k, max_log);
    }
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
        println!("  Using M = {:.0} from --M", m);
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
                    "  No {} files found alongside --ref-ld-chr; using M = {} (regression SNPs), may not be correct because M should be # of SNPs used for estimating LD score",
                    suffix, n_obs
                );
            }
        }
    } else {
        println!(
            "  Single-file --ref-ld: no M files to read. \
             Use --M to set M explicitly. Using M = {} (regression SNPs), may not be correct because M should be # of SNPs used for estimating LD score",
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
    keep_idx: Option<&[usize]>,
    k_full: Option<usize>,
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
            Ok(m_vec) => {
                if m_vec.len() == k {
                    let total: f64 = m_vec.iter().sum();
                    println!(
                        "  Read per-annotation M (K={}) from {} files; total M = {:.0}",
                        k, suffix, total
                    );
                    return m_vec;
                }
                if let (Some(idx), Some(k_full)) = (keep_idx, k_full) && m_vec.len() == k_full {
                    let total: f64 = m_vec.iter().sum();
                    let mut kept = Vec::with_capacity(idx.len());
                    for &i in idx {
                        if i >= m_vec.len() {
                            break;
                        }
                        kept.push(m_vec[i]);
                    }
                    println!(
                        "  Read per-annotation M (K={} original, K={} kept) from {} files; total M = {:.0}",
                        k_full,
                        kept.len(),
                        suffix,
                        total
                    );
                    return kept;
                }
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
// Helper: drop zero-variance LD score columns (partitioned h2 parity).
// ---------------------------------------------------------------------------

fn column_is_zero_variance(values: &Array1<f64>) -> bool {
    let mut count = 0usize;
    let mut first: Option<f64> = None;
    for &v in values.iter() {
        if v.is_nan() {
            continue;
        }
        count += 1;
        if let Some(f) = first {
            if v != f {
                return false;
            }
        } else {
            first = Some(v);
        }
    }
    if count < 2 {
        return false;
    }
    true
}

fn drop_zero_variance_ld(ref_ld: &DataFrame) -> Result<(DataFrame, Vec<usize>, usize)> {
    let mut keep_cols: Vec<String> = Vec::new();
    let mut keep_idx: Vec<usize> = Vec::new();
    let mut any_zero = false;
    let mut total = 0usize;

    for name in ref_ld.get_column_names() {
        if name == "SNP" {
            keep_cols.push(name.to_string());
            continue;
        }
        let values = extract_f64(ref_ld, name)?;
        if column_is_zero_variance(&values) {
            any_zero = true;
        } else {
            keep_cols.push(name.to_string());
            keep_idx.push(total);
        }
        total += 1;
    }

    anyhow::ensure!(total > 0, "No L2 annotation columns found in reference LD scores");
    if keep_cols.len() == 1 {
        anyhow::bail!("All L2 annotation columns have zero variance.");
    }
    if any_zero {
        println!("  Removing partitioned L2 columns with zero variance.");
    }

    let df = ref_ld.select(&keep_cols)?;
    Ok((df, keep_idx, total))
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

fn filter_strings_by_mask(v: &[String], mask: &[bool]) -> Vec<String> {
    v.iter()
        .zip(mask.iter())
        .filter_map(|(val, keep)| if *keep { Some(val.clone()) } else { None })
        .collect()
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
    if args.return_silly_things {
        println!(
            "WARNING: --return-silly-things is accepted for CLI parity but has no effect in Rust."
        );
    }
    if args.invert_anyway {
        println!("WARNING: --invert-anyway is accepted for CLI parity but has no effect in Rust.");
    }
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
    let z = col("Z").cast(DataType::Float64);
    let sumstats_lf = parse::scan_sumstats(h2_path)?
        .with_columns([
            col("N").cast(DataType::Float64).alias("N"),
            (z.clone() * z).alias("CHI2"),
        ])
        .select([col("SNP"), col("N"), col("CHI2")]);
    // Reference LD keeps all annotation columns (scalar or partitioned).
    let ref_ld =
        load_ld_ref(args.ref_ld.as_deref(), args.ref_ld_chr.as_deref())?.with_columns([all()
            .exclude_cols(["SNP"])
            .as_expr()
            .cast(DataType::Float64)]);
    let w_ld = load_ld(args.w_ld.as_deref(), args.w_ld_chr.as_deref(), "w_l2")?
        .with_columns([col("w_l2").cast(DataType::Float64).alias("w_l2")]);

    let merge_start = Instant::now();
    let sumstats_df = sumstats_lf.collect().context("loading sumstats")?;
    let ref_ld_df = ref_ld.collect().context("loading ref LD scores")?;
    let (ref_ld_df, keep_idx, k_full) = drop_zero_variance_ld(&ref_ld_df)?;
    let merged = smart_merge_on_snp(sumstats_df, ref_ld_df)
        .context("merging sumstats with reference LD scores")?;
    let w_ld_df = w_ld.collect().context("loading weight LD scores")?;
    let merged = smart_merge_on_snp(merged, w_ld_df).context("merging with weight LD scores")?;
    trace!("h2 merge time: {:?}", merge_start.elapsed());

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
    let non_l2_set: std::collections::HashSet<&str> =
        ["SNP", "CHI2", "N", "w_l2"].iter().copied().collect();
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

    let chi2_raw = extract_f64(&merged, "CHI2")?;
    let w_l2_raw = extract_f64(&merged, "w_l2")?;
    let n_vec_raw = extract_f64(&merged, "N")?;
    let ref_l2_raw_k: Vec<Array1<f64>> = l2_cols
        .iter()
        .map(|name| extract_f64(&merged, name))
        .collect::<Result<Vec<_>>>()?;

    let valid_mask: Vec<bool> = (0..chi2_raw.len())
        .map(|i| !chi2_raw[i].is_nan() && !n_vec_raw[i].is_nan())
        .collect();

    let mut chi2 = filter_by_mask(&chi2_raw, &valid_mask);
    let mut w_l2 = filter_by_mask(&w_l2_raw, &valid_mask);
    let mut n_vec = filter_by_mask(&n_vec_raw, &valid_mask);
    let mut ref_l2_k: Vec<Array1<f64>> = ref_l2_raw_k
        .iter()
        .map(|v| filter_by_mask(v, &valid_mask))
        .collect();

    let chisq_max = if k > 1 {
        match args.chisq_max {
            Some(c) => Some(c),
            None => {
                let max_n = n_vec.iter().cloned().fold(0.0f64, f64::max);
                Some((0.001 * max_n).max(80.0))
            }
        }
    } else {
        args.chisq_max
    };

    warn_small_snp_count(chi2.len());

    if let Some(chisq_max) = chisq_max {
        let mask: Vec<bool> = chi2.iter().map(|&c| c < chisq_max).collect();
        let n_removed = mask.iter().filter(|&&b| !b).count();
        if n_removed > 0 {
            println!("  Removed {} SNPs with chi2 > {:.1}", n_removed, chisq_max);
        }
        chi2 = filter_by_mask(&chi2, &mask);
        w_l2 = filter_by_mask(&w_l2, &mask);
        n_vec = filter_by_mask(&n_vec, &mask);
        ref_l2_k = ref_l2_k
            .iter()
            .map(|v| filter_by_mask(v, &mask))
            .collect::<Vec<_>>();
    }
    let n_obs = chi2.len();
    let n_mean = n_vec.mean().unwrap_or(1.0);
    debug!(
        "h2 regression inputs: n_obs={} n_mean={:.3} k={}",
        n_obs, n_mean, k
    );
    trace_array_stats("chi2", &chi2);
    trace_array_stats("w_l2", &w_l2);
    trace_array_stats("N", &n_vec);
    trace_ref_l2_stats(&l2_cols, &ref_l2_k);

    if k > 1 {
        let reg_start = Instant::now();
        if args.two_step.is_some() {
            anyhow::bail!("--two-step is not compatible with partitioned h2");
        }
        let keep_idx_opt = if keep_idx.len() == k_full {
            None
        } else {
            Some(keep_idx.as_slice())
        };
        let m_vec = resolve_m_vec(
            args.m_snps,
            args.ref_ld_chr.as_deref(),
            args.not_m_5_50,
            n_obs,
            k,
            keep_idx_opt,
            Some(k_full),
        );
        let fit = run_h2_partitioned(
            &chi2,
            &ref_l2_k,
            &w_l2,
            &n_vec,
            &l2_cols,
            n_obs,
            n_mean,
            &args,
            Some(&m_vec),
        )?;
        trace!("h2 regression time: {:?}", reg_start.elapsed());

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

    // Scalar h2 (K == 1) — match Python LDSC weights + IRWLS + jackknife.
    let ref_l2 = ref_l2_k.into_iter().next().unwrap();
    let m_snps = resolve_m(
        args.m_snps,
        args.ref_ld_chr.as_deref(),
        args.not_m_5_50,
        n_obs,
    );
    let n_blocks = n_obs.min(args.n_blocks);

    // Priority: --intercept-h2 > --no-intercept > free intercept.
    let fixed_intercept = if let Some(fixed) = args.intercept_h2 {
        Some(fixed)
    } else if args.no_intercept {
        Some(1.0)
    } else {
        None
    };
    let mut two_step = args.two_step;
    if two_step.is_none() && fixed_intercept.is_none() {
        two_step = Some(30.0);
    }
    if let Some(ts) = two_step {
        println!("Using two-step estimator with cutoff at {ts}.");
    }

    let reg_start = Instant::now();
    let res = run_h2_ldsc(
        &chi2,
        &ref_l2,
        &w_l2,
        &n_vec,
        m_snps,
        n_blocks,
        two_step,
        fixed_intercept,
    )?;
    trace!("h2 regression time: {:?}", reg_start.elapsed());

    println!("Total Observed scale h2: {:.4} ({:.4})", res.h2, res.h2_se);
    println!("Lambda GC: {:.4}", res.lambda_gc);
    println!("Mean Chi^2: {:.4}", res.mean_chi2);
    if fixed_intercept.is_some() {
        println!("Intercept: constrained to {:.4}", res.intercept);
    } else {
        println!("Intercept: {:.4} ({:.4})", res.intercept, res.intercept_se);
        if let Some((ratio, ratio_se)) = res.ratio {
            if ratio < 0.0 {
                println!("Ratio < 0 (usually indicates GC correction).");
            } else {
                println!("Ratio: {:.4} ({:.4})", ratio, ratio_se);
            }
        } else {
            println!("Ratio: NA (mean chi^2 < 1)");
        }
    }

    if let (Some(samp_prev), Some(pop_prev)) = (args.samp_prev, args.pop_prev) {
        let c = liability_conversion_factor(samp_prev, pop_prev);
        println!("Liability-scale h2: {:.4}", res.h2 * c);
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

    let z = col("Z").cast(DataType::Float64);
    let sumstats_lf = parse::scan_sumstats(sumstats_path)?
        .with_columns([
            col("N").cast(DataType::Float64).alias("N"),
            (z.clone() * z).alias("CHI2"),
        ])
        .select([col("SNP"), col("N"), col("CHI2")]);
    let ref_ld =
        load_ld_ref(args.ref_ld.as_deref(), args.ref_ld_chr.as_deref())?.with_columns([all()
            .exclude_cols(["SNP"])
            .as_expr()
            .cast(DataType::Float64)]);
    let w_ld = load_ld(args.w_ld.as_deref(), args.w_ld_chr.as_deref(), "w_l2")?
        .with_columns([col("w_l2").cast(DataType::Float64).alias("w_l2")]);

    let sumstats_df = sumstats_lf.collect().context("loading sumstats")?;
    let ref_ld_df = ref_ld.collect().context("loading ref LD scores")?;
    let (ref_ld_df, keep_idx, k_full) = drop_zero_variance_ld(&ref_ld_df)?;
    let merged = smart_merge_on_snp(sumstats_df, ref_ld_df)
        .context("merging sumstats with reference LD scores")?;
    let w_ld_df = w_ld.collect().context("loading weight LD scores")?;
    let merged = smart_merge_on_snp(merged, w_ld_df).context("merging with weight LD scores")?;

    let n_total = merged.height();
    anyhow::ensure!(
        n_total > 0,
        "No SNPs remaining after merging sumstats with LD scores"
    );
    println!(
        "Loaded {} SNPs from '{}' after merging with LD scores",
        n_total, sumstats_path
    );

    let non_l2_set: std::collections::HashSet<&str> = ["SNP", "CHI2", "N", "w_l2"]
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

    let chi2_raw = extract_f64(&merged, "CHI2")?;
    let w_l2_raw = extract_f64(&merged, "w_l2")?;
    let n_vec_raw = extract_f64(&merged, "N")?;
    let ref_l2_raw_k: Vec<Array1<f64>> = l2_cols
        .iter()
        .map(|name| extract_f64(&merged, name))
        .collect::<Result<Vec<_>>>()?;

    let snp_raw: Vec<String> = merged
        .column("SNP")?
        .as_materialized_series()
        .str()?
        .into_iter()
        .map(|opt| opt.unwrap_or("").to_string())
        .collect();

    let valid_mask: Vec<bool> = (0..chi2_raw.len())
        .map(|i| !chi2_raw[i].is_nan() && !n_vec_raw[i].is_nan())
        .collect();

    let mut chi2 = filter_by_mask(&chi2_raw, &valid_mask);
    let mut w_l2 = filter_by_mask(&w_l2_raw, &valid_mask);
    let mut n_vec = filter_by_mask(&n_vec_raw, &valid_mask);
    let mut ref_l2_k: Vec<Array1<f64>> = ref_l2_raw_k
        .iter()
        .map(|v| filter_by_mask(v, &valid_mask))
        .collect();
    let mut keep_snps = filter_strings_by_mask(&snp_raw, &valid_mask);

    warn_small_snp_count(chi2.len());

    let chisq_max = if let Some(c) = args.chisq_max {
        Some(c)
    } else {
        let max_n = n_vec.iter().cloned().fold(0.0f64, f64::max);
        Some((0.001 * max_n).max(80.0))
    };

    if let Some(cmax) = chisq_max {
        let mask: Vec<bool> = chi2.iter().map(|&c| c < cmax).collect();
        let n_removed = mask.iter().filter(|&&b| !b).count();
        if n_removed > 0 {
            println!("  Removed {} SNPs with chi2 > {:.1}", n_removed, cmax);
        }
        chi2 = filter_by_mask(&chi2, &mask);
        w_l2 = filter_by_mask(&w_l2, &mask);
        n_vec = filter_by_mask(&n_vec, &mask);
        ref_l2_k = ref_l2_k
            .iter()
            .map(|v| filter_by_mask(v, &mask))
            .collect::<Vec<_>>();
        keep_snps = filter_strings_by_mask(&keep_snps, &mask);
    }

    let n_obs = chi2.len();
    let n_mean = n_vec.mean().unwrap_or(1.0);

    let keep_idx_opt = if keep_idx.len() == k_full {
        None
    } else {
        Some(keep_idx.as_slice())
    };
    let m_base = resolve_m_vec(
        args.m_snps,
        args.ref_ld_chr.as_deref(),
        args.not_m_5_50,
        n_obs,
        k_base,
        keep_idx_opt,
        Some(k_full),
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
                v.push(val);
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
            None,
            None,
        )
    };
    let m_total: f64 = m_vec.iter().sum();
    let n_blocks = n_obs.min(args.n_blocks);

    let mut x_raw = Array2::<f64>::zeros((n_obs, k));
    for (j, col_v) in ref_l2_k.iter().enumerate() {
        for i in 0..n_obs {
            x_raw[[i, j]] = col_v[i];
        }
    }
    let mut x_tot = Array1::<f64>::zeros(n_obs);
    for i in 0..n_obs {
        let mut sum = 0.0;
        for j in 0..k {
            sum += x_raw[[i, j]];
        }
        x_tot[i] = sum;
    }

    let fixed_intercept = if let Some(fixed) = args.intercept_h2 {
        Some(fixed)
    } else if args.no_intercept {
        Some(1.0)
    } else {
        None
    };

    let y_reg = if let Some(fixed) = fixed_intercept {
        chi2.mapv(|c| c - fixed)
    } else {
        chi2.clone()
    };

    let x_cols = if fixed_intercept.is_some() { k } else { k + 1 };
    let mut x_reg = Array2::<f64>::zeros((n_obs, x_cols));
    for i in 0..n_obs {
        let scale = n_vec[i] / n_mean;
        for j in 0..k {
            x_reg[[i, j]] = x_raw[[i, j]] * scale;
        }
        if fixed_intercept.is_none() {
            x_reg[[i, k]] = 1.0;
        }
    }

    let agg_intercept = fixed_intercept.unwrap_or(1.0);
    let hsq = aggregate(chi2, &x_tot, n_vec, m_total, agg_intercept);
    let initial_w = ldsc_weights(&x_tot, w_l2, n_vec, m_total, hsq, agg_intercept);
    let w_sqrt = initial_w.mapv(|v| v.sqrt());
    let (xw, yw) = weight_xy(&x_reg, &y_reg, &w_sqrt)?;
    let jknife: JackknifeResult =
        jackknife_fast(&xw, &yw, n_blocks, None).context("jackknife for partitioned h2")?;

    let coef = jknife.est.slice(s![..k]).to_owned() / n_mean;
    let h2_per_annot: Vec<f64> = (0..k).map(|j| coef[j] * m_vec[j]).collect();
    let h2_total: f64 = h2_per_annot.iter().sum();
    let intercept = fixed_intercept.unwrap_or_else(|| jknife.est[k]);

    Ok(PartitionedFit {
        est: jknife.est,
        jknife_cov: jknife.jknife_cov,
        delete_values: jknife.delete_values,
        m_vec,
        n_mean,
        h2_per_annot,
        h2_total,
        intercept,
        n_blocks,
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
    m_override: Option<&[f64]>,
) -> Result<PartitionedFit> {
    let fit = fit_h2_partitioned(
        chi2,
        ref_l2_k,
        w_l2,
        n_vec,
        n_obs,
        n_mean,
        args,
        m_override,
    )?;

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

#[allow(dead_code)]
#[derive(Debug, Clone)]
struct HsqFit {
    tot: f64,
    tot_se: f64,
    intercept: f64,
    intercept_se: f64,
    mean_chi2: f64,
    lambda_gc: f64,
    ratio: Option<(f64, f64)>,
    tot_delete_values: Array1<f64>,
}

#[allow(dead_code)]
#[derive(Debug, Clone)]
struct GencovFit {
    tot: f64,
    tot_se: f64,
    intercept: f64,
    intercept_se: f64,
    mean_z1z2: f64,
    tot_delete_values: Array1<f64>,
}

fn mean(arr: &Array1<f64>) -> f64 {
    arr.mean().unwrap_or(f64::NAN)
}

fn sum_ref_l2(ref_l2_k: &[Array1<f64>]) -> Result<Array1<f64>> {
    anyhow::ensure!(!ref_l2_k.is_empty(), "No L2 columns supplied");
    let n = ref_l2_k[0].len();
    let mut out = Array1::<f64>::zeros(n);
    for col in ref_l2_k {
        anyhow::ensure!(col.len() == n, "L2 column length mismatch");
        out += col;
    }
    Ok(out)
}

fn total_from_jknife(
    jknife: &JackknifeResult,
    m_vec: &[f64],
    nbar: f64,
) -> Result<(f64, f64, Array1<f64>)> {
    let k = m_vec.len();
    anyhow::ensure!(
        jknife.est.len() >= k,
        "Jackknife estimate has {} params but K={}",
        jknife.est.len(),
        k
    );
    let coef_cov = jknife.jknife_cov.slice(s![0..k, 0..k]).to_owned() / (nbar * nbar);

    let mut tot = 0.0;
    for (j, m) in m_vec.iter().enumerate().take(k) {
        tot += m * (jknife.est[j] / nbar);
    }

    let mut tot_cov = 0.0;
    for i in 0..k {
        for j in 0..k {
            tot_cov += m_vec[i] * m_vec[j] * coef_cov[[i, j]];
        }
    }
    let tot_se = tot_cov.max(0.0).sqrt();

    let n_blocks = jknife.delete_values.nrows();
    let mut tot_delete = Array1::<f64>::zeros(n_blocks);
    for b in 0..n_blocks {
        let mut sum = 0.0;
        for (j, m) in m_vec.iter().enumerate().take(k) {
            sum += m * jknife.delete_values[[b, j]];
        }
        tot_delete[b] = sum / nbar;
    }

    Ok((tot, tot_se, tot_delete))
}

#[allow(clippy::too_many_arguments)]
fn gencov_weights(
    ld: &Array1<f64>,
    w_ld: &Array1<f64>,
    n1: &Array1<f64>,
    n2: &Array1<f64>,
    m_total: f64,
    h1: f64,
    h2: f64,
    rho_g: f64,
    intercept_gencov: f64,
    intercept_hsq1: f64,
    intercept_hsq2: f64,
) -> Array1<f64> {
    let h1 = h1.clamp(0.0, 1.0);
    let h2 = h2.clamp(0.0, 1.0);
    let rho_g = rho_g.clamp(-1.0, 1.0);

    let mut out = Array1::<f64>::zeros(ld.len());
    for i in 0..ld.len() {
        let ld_i = ld[i].max(1.0);
        let w_i = w_ld[i].max(1.0);
        let a = n1[i] * h1 * ld_i / m_total + intercept_hsq1;
        let b = n2[i] * h2 * ld_i / m_total + intercept_hsq2;
        let c = (n1[i] * n2[i]).sqrt() * rho_g * ld_i / m_total + intercept_gencov;
        let het_w = 1.0 / (a * b + c * c);
        let oc_w = 1.0 / w_i;
        out[i] = het_w * oc_w;
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn run_hsq_ldsc(
    chi2: &Array1<f64>,
    ref_l2_k: &[Array1<f64>],
    w_l2: &Array1<f64>,
    n_vec: &Array1<f64>,
    m_vec: &[f64],
    n_blocks: usize,
    two_step: Option<f64>,
    fixed_intercept: Option<f64>,
) -> Result<HsqFit> {
    let n = chi2.len();
    let k = ref_l2_k.len();
    anyhow::ensure!(k > 0, "run_hsq_ldsc: no L2 columns");
    anyhow::ensure!(
        w_l2.len() == n && n_vec.len() == n,
        "run_hsq_ldsc: length mismatch"
    );
    for col in ref_l2_k {
        anyhow::ensure!(col.len() == n, "run_hsq_ldsc: L2 length mismatch");
    }

    let nbar = mean(n_vec);
    let m_total: f64 = m_vec.iter().sum();
    let ref_l2_tot = sum_ref_l2(ref_l2_k)?;

    let intercept0 = fixed_intercept.unwrap_or(1.0);
    let tot_agg = aggregate(chi2, &ref_l2_tot, n_vec, m_total, intercept0);
    let initial_w = ldsc_weights(&ref_l2_tot, w_l2, n_vec, m_total, tot_agg, intercept0);

    let mut x_scaled = Array2::<f64>::zeros((n, k));
    for (j, col) in ref_l2_k.iter().enumerate() {
        for i in 0..n {
            x_scaled[[i, j]] = n_vec[i] * col[i] / nbar;
        }
    }

    let (y_reg, x_reg) = if let Some(fixed) = fixed_intercept {
        let y_adj = chi2.mapv(|c| c - fixed);
        (y_adj, x_scaled.clone())
    } else {
        let mut x_i = Array2::<f64>::zeros((n, k + 1));
        x_i.slice_mut(s![.., 0..k]).assign(&x_scaled);
        x_i.column_mut(k).fill(1.0);
        (chi2.clone(), x_i)
    };

    if two_step.is_some() && fixed_intercept.is_some() {
        anyhow::bail!("twostep is not compatible with constrained intercept");
    }
    if two_step.is_some() && k > 1 {
        anyhow::bail!("twostep not compatible with partitioned LD Score");
    }

    let jknife = if let Some(twostep) = two_step {
        let mask: Vec<bool> = chi2.iter().map(|&c| c < twostep).collect();
        let n_s1 = mask.iter().filter(|&&b| b).count();

        let mut x1 = Array2::<f64>::zeros((n_s1, x_reg.ncols()));
        let mut y1 = Array1::<f64>::zeros(n_s1);
        let mut w1 = Array1::<f64>::zeros(n_s1);
        let mut n1v = Array1::<f64>::zeros(n_s1);
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
            let hsq = m_total * coef[0] / nbar;
            let intercept = coef[1];
            ldsc_weights(&x1_col0, &w1, &n1v, m_total, hsq, intercept)
        };

        let step1 = irwls_ldsc(&x1, &y1, &initial_w1, update_func1, n_blocks, None)?;
        let step1_int = step1.est[k];

        let y2 = y_reg.mapv(|c| c - step1_int);
        let x2 = x_scaled.clone();
        let update_func2 = |coef: &Array1<f64>| {
            let hsq = m_total * coef[0] / nbar;
            ldsc_weights(&ref_l2_tot, w_l2, n_vec, m_total, hsq, step1_int)
        };

        let seps = update_separators(&step1.separators, &mask);
        let step2 = irwls_ldsc(&x2, &y2, &initial_w, update_func2, n_blocks, Some(&seps))?;

        let num = (&initial_w * &x2.column(0)).sum();
        let denom = (&initial_w * &x2.column(0).mapv(|v| v * v)).sum();
        let c = num / denom;
        combine_twostep(&step1, &step2, c)?
    } else {
        let update_func = |coef: &Array1<f64>| {
            let hsq = m_total * coef[0] / nbar;
            let intercept = fixed_intercept.unwrap_or(coef[k]);
            ldsc_weights(&ref_l2_tot, w_l2, n_vec, m_total, hsq, intercept)
        };
        irwls_ldsc(&x_reg, &y_reg, &initial_w, update_func, n_blocks, None)?
    };

    let (tot, tot_se, tot_delete_values) = total_from_jknife(&jknife, m_vec, nbar)?;

    let (intercept, intercept_se) = if let Some(fixed) = fixed_intercept {
        (fixed, f64::NAN)
    } else {
        (jknife.est[k], jknife.jknife_se[k])
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

    Ok(HsqFit {
        tot,
        tot_se,
        intercept,
        intercept_se,
        mean_chi2,
        lambda_gc,
        ratio,
        tot_delete_values,
    })
}

#[allow(clippy::too_many_arguments)]
fn run_gencov_ldsc(
    z1: &Array1<f64>,
    z2: &Array1<f64>,
    ref_l2_k: &[Array1<f64>],
    w_l2: &Array1<f64>,
    n1: &Array1<f64>,
    n2: &Array1<f64>,
    m_vec: &[f64],
    n_blocks: usize,
    two_step: Option<f64>,
    intercept_gencov: Option<f64>,
    hsq1: &HsqFit,
    hsq2: &HsqFit,
) -> Result<GencovFit> {
    let n = z1.len();
    let k = ref_l2_k.len();
    anyhow::ensure!(k > 0, "run_gencov_ldsc: no L2 columns");
    anyhow::ensure!(
        z2.len() == n && w_l2.len() == n && n1.len() == n && n2.len() == n,
        "run_gencov_ldsc: length mismatch"
    );
    for col in ref_l2_k {
        anyhow::ensure!(col.len() == n, "run_gencov_ldsc: L2 length mismatch");
    }

    let m_total: f64 = m_vec.iter().sum();
    let ref_l2_tot = sum_ref_l2(ref_l2_k)?;
    let sqrt_n1n2 = Array1::from_iter(n1.iter().zip(n2.iter()).map(|(&a, &b)| (a * b).sqrt()));
    let nbar = mean(&sqrt_n1n2);

    let intercept0 = intercept_gencov.unwrap_or(0.0);
    let prod = z1 * z2;
    let tot_agg = aggregate(&prod, &ref_l2_tot, &sqrt_n1n2, m_total, intercept0);
    let initial_w = gencov_weights(
        &ref_l2_tot,
        w_l2,
        n1,
        n2,
        m_total,
        hsq1.tot,
        hsq2.tot,
        tot_agg,
        intercept0,
        hsq1.intercept,
        hsq2.intercept,
    );

    let mut x_scaled = Array2::<f64>::zeros((n, k));
    for (j, col) in ref_l2_k.iter().enumerate() {
        for i in 0..n {
            x_scaled[[i, j]] = sqrt_n1n2[i] * col[i] / nbar;
        }
    }

    let (y_reg, x_reg) = if let Some(fixed) = intercept_gencov {
        let y_adj = prod.mapv(|p| p - fixed);
        (y_adj, x_scaled.clone())
    } else {
        let mut x_i = Array2::<f64>::zeros((n, k + 1));
        x_i.slice_mut(s![.., 0..k]).assign(&x_scaled);
        x_i.column_mut(k).fill(1.0);
        (prod.clone(), x_i)
    };

    if two_step.is_some() && intercept_gencov.is_some() {
        anyhow::bail!("twostep is not compatible with constrained intercept");
    }
    if two_step.is_some() && k > 1 {
        anyhow::bail!("twostep not compatible with partitioned LD Score");
    }

    let jknife = if let Some(twostep) = two_step {
        let mask: Vec<bool> = z1
            .iter()
            .zip(z2.iter())
            .map(|(&a, &b)| a * a < twostep && b * b < twostep)
            .collect();
        let n_s1 = mask.iter().filter(|&&b| b).count();

        let mut x1 = Array2::<f64>::zeros((n_s1, x_reg.ncols()));
        let mut y1 = Array1::<f64>::zeros(n_s1);
        let mut w1 = Array1::<f64>::zeros(n_s1);
        let mut n1v = Array1::<f64>::zeros(n_s1);
        let mut n2v = Array1::<f64>::zeros(n_s1);
        let mut idx = 0usize;
        for i in 0..n {
            if mask[i] {
                x1.row_mut(idx).assign(&x_reg.row(i));
                y1[idx] = y_reg[i];
                w1[idx] = w_l2[i];
                n1v[idx] = n1[i];
                n2v[idx] = n2[i];
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
            let rho_g = m_total * coef[0] / nbar;
            let intercept = coef[1];
            gencov_weights(
                &x1_col0,
                &w1,
                &n1v,
                &n2v,
                m_total,
                hsq1.tot,
                hsq2.tot,
                rho_g,
                intercept,
                hsq1.intercept,
                hsq2.intercept,
            )
        };

        let step1 = irwls_ldsc(&x1, &y1, &initial_w1, update_func1, n_blocks, None)?;
        let step1_int = step1.est[k];

        let y2 = y_reg.mapv(|p| p - step1_int);
        let x2 = x_scaled.clone();
        let update_func2 = |coef: &Array1<f64>| {
            let rho_g = m_total * coef[0] / nbar;
            gencov_weights(
                &ref_l2_tot,
                w_l2,
                n1,
                n2,
                m_total,
                hsq1.tot,
                hsq2.tot,
                rho_g,
                step1_int,
                hsq1.intercept,
                hsq2.intercept,
            )
        };

        let seps = update_separators(&step1.separators, &mask);
        let step2 = irwls_ldsc(&x2, &y2, &initial_w, update_func2, n_blocks, Some(&seps))?;

        let num = (&initial_w * &x2.column(0)).sum();
        let denom = (&initial_w * &x2.column(0).mapv(|v| v * v)).sum();
        let c = num / denom;
        combine_twostep(&step1, &step2, c)?
    } else {
        let update_func = |coef: &Array1<f64>| {
            let rho_g = m_total * coef[0] / nbar;
            let intercept = intercept_gencov.unwrap_or(coef[k]);
            gencov_weights(
                &ref_l2_tot,
                w_l2,
                n1,
                n2,
                m_total,
                hsq1.tot,
                hsq2.tot,
                rho_g,
                intercept,
                hsq1.intercept,
                hsq2.intercept,
            )
        };
        irwls_ldsc(&x_reg, &y_reg, &initial_w, update_func, n_blocks, None)?
    };

    let (tot, tot_se, tot_delete_values) = total_from_jknife(&jknife, m_vec, nbar)?;
    let (intercept, intercept_se) = if let Some(fixed) = intercept_gencov {
        (fixed, f64::NAN)
    } else {
        (jknife.est[k], jknife.jknife_se[k])
    };
    let mean_z1z2 = mean(&prod);

    Ok(GencovFit {
        tot,
        tot_se,
        intercept,
        intercept_se,
        mean_z1z2,
        tot_delete_values,
    })
}

pub fn run_rg(args: RgArgs) -> Result<()> {
    anyhow::ensure!(
        args.rg.len() >= 2,
        "--rg requires at least 2 sumstats files"
    );
    if args.return_silly_things {
        println!(
            "WARNING: --return-silly-things is accepted for CLI parity but has no effect in Rust."
        );
    }
    if args.invert_anyway {
        println!("WARNING: --invert-anyway is accepted for CLI parity but has no effect in Rust.");
    }

    let n_traits = args.rg.len();
    let mut intercept_h2: Vec<Option<f64>> = if args.intercept_h2.is_empty() {
        vec![None; n_traits]
    } else {
        anyhow::ensure!(
            args.intercept_h2.len() == n_traits,
            "--intercept-h2 expects one value per trait ({} values for this run)",
            n_traits
        );
        args.intercept_h2.iter().copied().map(Some).collect()
    };

    let mut intercept_gencov: Vec<Option<f64>> = if args.intercept_gencov.is_empty() {
        vec![None; n_traits]
    } else if args.intercept_gencov.len() == n_traits {
        args.intercept_gencov.iter().copied().map(Some).collect()
    } else if args.intercept_gencov.len() == n_traits - 1 {
        println!(
            "WARNING: --intercept-gencov has {} values; mapping to trait1-vs-others pairs. \
To match Python, provide {} values (first ignored).",
            args.intercept_gencov.len(),
            n_traits
        );
        let mut vals = vec![None; n_traits];
        for (i, v) in args.intercept_gencov.iter().copied().enumerate() {
            vals[i + 1] = Some(v);
        }
        vals
    } else {
        anyhow::bail!(
            "--intercept-gencov expects either {} values (per trait, first ignored) or {} values (trait1 vs others)",
            n_traits,
            n_traits.saturating_sub(1)
        );
    };

    if args.no_intercept {
        intercept_h2 = vec![Some(1.0); n_traits];
        intercept_gencov = vec![Some(0.0); n_traits];
    }

    if !args.samp_prev.is_empty() {
        anyhow::ensure!(
            args.samp_prev.len() == n_traits,
            "--samp-prev expects one value per trait ({} values for this run)",
            n_traits
        );
    }
    if !args.pop_prev.is_empty() {
        anyhow::ensure!(
            args.pop_prev.len() == n_traits,
            "--pop-prev expects one value per trait ({} values for this run)",
            n_traits
        );
    }

    let ref_ld = load_ld_ref(args.ref_ld.as_deref(), args.ref_ld_chr.as_deref())?
        .collect()
        .context("loading ref LD scores")?;
    let w_ld = load_ld(args.w_ld.as_deref(), args.w_ld_chr.as_deref(), "w_l2")?
        .collect()
        .context("loading weight LD scores")?;

    let ref_l2_cols: Vec<String> = ref_ld
        .get_column_names()
        .into_iter()
        .filter(|n| n.as_str() != "SNP")
        .map(|s| s.to_string())
        .collect();
    anyhow::ensure!(
        !ref_l2_cols.is_empty(),
        "No L2 annotation columns found in reference LD scores"
    );
    let k = ref_l2_cols.len();

    let two_step = args.two_step;
    if let Some(ts) = two_step {
        println!("Using two-step estimator with cutoff at {ts}.");
    }

    let file1 = &args.rg[0];
    if !args.no_check_alleles {
        ensure_sumstats_have_alleles(file1)?;
    }

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

    let ss1_df = ss1.collect().context("loading primary sumstats")?;

    for trait_idx in 1..n_traits {
        let file2 = &args.rg[trait_idx];
        println!("Computing rg: {} vs {}", file1, file2);
        if !args.no_check_alleles {
            ensure_sumstats_have_alleles(file2)?;
        }

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

        let merge_start = Instant::now();
        let merged = smart_merge_on_snp(ss1_df.clone(), ss2.collect()?)
            .with_context(|| format!("merging {} vs {}", file1, file2))?;
        trace!("rg merge time: {:?}", merge_start.elapsed());

        let n_obs = merged.height();
        if n_obs == 0 {
            println!("  No overlapping SNPs — skipping pair");
            continue;
        }

        let z1_raw = extract_f64(&merged, "Z1")?;
        let z2_raw = extract_f64(&merged, "Z2")?;
        let w_l2_raw = extract_f64(&merged, "w_l2")?.mapv(|l| l.max(0.0));
        let n1_raw = extract_f64(&merged, "N1")?;
        let n2_raw = extract_f64(&merged, "N2")?;
        let ref_l2_raw_k: Vec<Array1<f64>> = ref_l2_cols
            .iter()
            .map(|name| extract_f64(&merged, name).map(|v| v.mapv(|l| l.max(0.0))))
            .collect::<Result<Vec<_>>>()?;

        let (z1_raw, z2_raw, w_l2_raw, n1_raw, n2_raw, ref_l2_raw_k) = if args.no_check_alleles {
            (z1_raw, z2_raw, w_l2_raw, n1_raw, n2_raw, ref_l2_raw_k)
        } else {
            let align_start = Instant::now();
            let (mask, z2_aligned, n_removed) =
                align_rg_alleles(&merged, &z2_raw).context("aligning alleles")?;
            trace!(
                "rg align time: {:?} (removed {})",
                align_start.elapsed(),
                n_removed
            );
            if n_removed > 0 {
                println!("  Removed {} SNPs with incompatible alleles", n_removed);
            }
            let z1_f = filter_by_mask(&z1_raw, &mask);
            let z2_f = Array1::from_vec(z2_aligned);
            let w_f = filter_by_mask(&w_l2_raw, &mask);
            let n1_f = filter_by_mask(&n1_raw, &mask);
            let n2_f = filter_by_mask(&n2_raw, &mask);
            let ref_f = ref_l2_raw_k
                .iter()
                .map(|v| filter_by_mask(v, &mask))
                .collect::<Vec<_>>();
            (z1_f, z2_f, w_f, n1_f, n2_f, ref_f)
        };

        let mut finite_mask = vec![true; z1_raw.len()];
        for i in 0..z1_raw.len() {
            if !z1_raw[i].is_finite()
                || !z2_raw[i].is_finite()
                || !w_l2_raw[i].is_finite()
                || !n1_raw[i].is_finite()
                || !n2_raw[i].is_finite()
            {
                finite_mask[i] = false;
                continue;
            }
            if ref_l2_raw_k.iter().any(|v| !v[i].is_finite()) {
                finite_mask[i] = false;
            }
        }
        if finite_mask.iter().any(|&v| !v) {
            let removed = finite_mask.iter().filter(|&&v| !v).count();
            println!("  Removed {} SNPs with non-finite values", removed);
        }
        let z1_raw = filter_by_mask(&z1_raw, &finite_mask);
        let z2_raw = filter_by_mask(&z2_raw, &finite_mask);
        let w_l2_raw = filter_by_mask(&w_l2_raw, &finite_mask);
        let n1_raw = filter_by_mask(&n1_raw, &finite_mask);
        let n2_raw = filter_by_mask(&n2_raw, &finite_mask);
        let ref_l2_raw_k: Vec<Array1<f64>> = ref_l2_raw_k
            .into_iter()
            .map(|v| filter_by_mask(&v, &finite_mask))
            .collect();

        let prod_raw = &z1_raw * &z2_raw;
        let (z1, z2, prod, w_l2, n1, n2, ref_l2_k) = if let Some(chisq_max) = args.chisq_max {
            let mask: Vec<bool> = prod_raw.iter().map(|&p| p.abs() < chisq_max).collect();
            let n_removed = mask.iter().filter(|&&b| !b).count();
            if n_removed > 0 {
                println!(
                    "  Removed {} SNPs with |Z1*Z2| >= {:.1}",
                    n_removed, chisq_max
                );
            }
            (
                filter_by_mask(&z1_raw, &mask),
                filter_by_mask(&z2_raw, &mask),
                filter_by_mask(&prod_raw, &mask),
                filter_by_mask(&w_l2_raw, &mask),
                filter_by_mask(&n1_raw, &mask),
                filter_by_mask(&n2_raw, &mask),
                ref_l2_raw_k
                    .iter()
                    .map(|v| filter_by_mask(v, &mask))
                    .collect::<Vec<_>>(),
            )
        } else {
            (
                z1_raw,
                z2_raw,
                prod_raw,
                w_l2_raw,
                n1_raw,
                n2_raw,
                ref_l2_raw_k,
            )
        };
        let n_obs_filtered = prod.len();
        if n_obs_filtered == 0 {
            println!("  No SNPs remaining after filtering — skipping pair");
            continue;
        }

        let n_blocks = n_obs_filtered.min(args.n_blocks);
        let m_vec = resolve_m_vec(
            args.m_snps,
            args.ref_ld_chr.as_deref(),
            args.not_m_5_50,
            n_obs_filtered,
            k,
            None,
            None,
        );
        let n1_mean = mean(&n1);
        let n2_mean = mean(&n2);
        debug!(
            "rg regression inputs: n_obs={} n_obs_filtered={} m_total={:.0} n1_mean={:.3} n2_mean={:.3}",
            n_obs,
            n_obs_filtered,
            m_vec.iter().sum::<f64>(),
            n1_mean,
            n2_mean
        );
        trace_array_stats("rg_prod", &prod);
        trace_array_stats("rg_w_l2", &w_l2);
        trace_array_stats("rg_n1", &n1);
        trace_array_stats("rg_n2", &n2);
        trace_ref_l2_stats(&ref_l2_cols, &ref_l2_k);

        let samp_prev1 = args.samp_prev.first().copied();
        let pop_prev1 = args.pop_prev.first().copied();
        let samp_prev2 = args.samp_prev.get(trait_idx).copied();
        let pop_prev2 = args.pop_prev.get(trait_idx).copied();

        let rg_reg_start = Instant::now();
        let h2_int1 = intercept_h2[0];
        let h2_int2 = intercept_h2[trait_idx];
        let gencov_int = intercept_gencov[trait_idx];

        let hsq1 = run_hsq_ldsc(
            &z1.mapv(|z| (z * z).max(0.0)),
            &ref_l2_k,
            &w_l2,
            &n1,
            &m_vec,
            n_blocks,
            two_step,
            h2_int1,
        )?;
        let hsq2 = run_hsq_ldsc(
            &z2.mapv(|z| (z * z).max(0.0)),
            &ref_l2_k,
            &w_l2,
            &n2,
            &m_vec,
            n_blocks,
            two_step,
            h2_int2,
        )?;

        let gencov = run_gencov_ldsc(
            &z1, &z2, &ref_l2_k, &w_l2, &n1, &n2, &m_vec, n_blocks, two_step, gencov_int, &hsq1,
            &hsq2,
        )?;

        let (rg_ratio, rg_se) = if hsq1.tot > 0.0 && hsq2.tot > 0.0 {
            let rg_ratio = gencov.tot / (hsq1.tot * hsq2.tot).sqrt();
            let denom_delete: Vec<f64> = hsq1
                .tot_delete_values
                .iter()
                .zip(hsq2.tot_delete_values.iter())
                .map(|(&a, &b)| {
                    let prod = a * b;
                    if prod.is_finite() && prod >= 0.0 {
                        prod.sqrt()
                    } else {
                        f64::NAN
                    }
                })
                .collect();
            let numer_delete =
                Array2::from_shape_vec((n_blocks, 1), gencov.tot_delete_values.to_vec())
                    .context("building rg numer delete values")?;
            let denom_delete = Array2::from_shape_vec((n_blocks, 1), denom_delete)
                .context("building rg denom delete values")?;
            let est = Array1::from_vec(vec![rg_ratio]);
            let (_cov, se) = ratio_jackknife(&est, &numer_delete, &denom_delete);
            (rg_ratio, se[0])
        } else {
            (f64::NAN, f64::NAN)
        };
        trace!("rg regression time: {:?}", rg_reg_start.elapsed());

        println!(
            "  gencov = {:.4} ({:.4})  (intercept: {:.4} ± {:.4})",
            gencov.tot, gencov.tot_se, gencov.intercept, gencov.intercept_se
        );
        println!(
            "  h2({}) = {:.4} ({:.4})  (intercept: {:.4} ± {:.4})",
            file1, hsq1.tot, hsq1.tot_se, hsq1.intercept, hsq1.intercept_se
        );
        println!(
            "  h2({}) = {:.4} ({:.4})  (intercept: {:.4} ± {:.4})",
            file2, hsq2.tot, hsq2.tot_se, hsq2.intercept, hsq2.intercept_se
        );
        if args.return_silly_things {
            println!("  rg = {:.4} ({:.4})", rg_ratio, rg_se);
        } else if rg_ratio.is_finite() && !(-1.2..=1.2).contains(&rg_ratio) {
            println!("  rg = nan (nan) (rg out of bounds)");
        } else {
            println!("  rg = {:.4} ({:.4})", rg_ratio, rg_se);
        }

        if let (Some(sp1), Some(pp1), Some(sp2), Some(pp2)) =
            (samp_prev1, pop_prev1, samp_prev2, pop_prev2)
        {
            let c1 = liability_conversion_factor(sp1, pp1);
            let c2 = liability_conversion_factor(sp2, pp2);
            let h2_1_liab = hsq1.tot * c1;
            let h2_2_liab = hsq2.tot * c2;
            let gencov_liab = gencov.tot * (c1 * c2).sqrt();
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
