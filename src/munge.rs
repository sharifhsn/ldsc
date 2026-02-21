/// Summary statistics munging — Polars LazyFrame pipeline for streaming large files.
use anyhow::{Context, Result};
use polars::prelude::*;
use statrs::distribution::{ContinuousCDF, Normal};
use std::fs::File;
use std::io::BufWriter;

use crate::cli::MungeArgs;
use crate::parse;

// ---------------------------------------------------------------------------
// Column-name synonym map
// ---------------------------------------------------------------------------

/// (uppercase_synonym, canonical_name) pairs.
const CNAME_MAP: &[(&str, &str)] = &[
    // SNP identifier
    ("SNP", "SNP"),
    ("MARKERNAME", "SNP"),
    ("SNPID", "SNP"),
    ("RS", "SNP"),
    ("RSID", "SNP"),
    ("RS_NUMBER", "SNP"),
    ("RS_NUMBERS", "SNP"),
    // P-value
    ("P", "P"),
    ("PVALUE", "P"),
    ("P_VALUE", "P"),
    ("PVAL", "P"),
    ("P_VAL", "P"),
    ("GC_PVALUE", "P"),
    // Allele 1 (effect allele)
    ("A1", "A1"),
    ("ALLELE1", "A1"),
    ("ALLELE_1", "A1"),
    ("EFFECT_ALLELE", "A1"),
    ("REFERENCE_ALLELE", "A1"),
    ("INC_ALLELE", "A1"),
    ("EA", "A1"),
    // Allele 2 (non-effect allele)
    ("A2", "A2"),
    ("ALLELE2", "A2"),
    ("ALLELE_2", "A2"),
    ("OTHER_ALLELE", "A2"),
    ("NON_EFFECT_ALLELE", "A2"),
    ("DEC_ALLELE", "A2"),
    ("NEA", "A2"),
    // Sample size
    ("N", "N"),
    ("WEIGHT", "N"), // METAL convention
    // Z-score
    ("Z", "Z"),
    ("ZSCORE", "Z"),
    ("Z-SCORE", "Z"),
    ("GC_ZSCORE", "Z"),
    // Regression coefficients / signed stats
    ("BETA", "BETA"),
    ("B", "BETA"),
    ("EFFECT", "BETA"),
    ("EFFECTS", "BETA"),
    ("OR", "OR"),
    ("LOG_ODDS", "LOG_ODDS"),
    ("SIGNED_SUMSTAT", "SIGNED_SUMSTAT"),
    // Standard error
    ("SE", "SE"),
    ("STDERR", "SE"),
    ("STDERROR", "SE"),
    ("SE_BETA", "SE"),
    // Allele frequency / MAF
    ("FRQ", "FRQ"),
    ("MAF", "FRQ"),
    ("EAF", "FRQ"),
    ("FRQ_U", "FRQ"),
    ("F_U", "FRQ"),
    // Imputation quality
    ("INFO", "INFO"),
    ("IMPINFO", "INFO"),
];

fn cname_lookup(upper: &str) -> Option<&'static str> {
    CNAME_MAP.iter().find(|(k, _)| *k == upper).map(|(_, v)| *v)
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

pub fn run(args: MungeArgs) -> Result<()> {
    if args.no_alleles && args.merge_alleles.is_some() {
        anyhow::bail!("--no-alleles and --merge-alleles are not compatible");
    }
    let lf = parse::scan_sumstats(&args.sumstats)?;
    let lf = apply_ignore(lf, args.ignore.as_deref())?;
    // Apply user overrides before synonym map so non-standard names are canonicalised.
    let lf = apply_col_overrides(lf, &args)?;
    let lf = normalize_columns(lf)?;
    let lf = apply_info_list(lf, args.info_list.as_deref())?;
    let lf = apply_n_override(lf, &args)?;
    let lf = derive_z(lf, args.signed_sumstats.as_deref(), args.a1_inc)?;
    let lf = filter_snps(lf, args.maf, args.n_min, args.info_min, args.no_alleles)?;
    let lf = apply_nstudy_filter(lf, args.nstudy.as_deref(), args.nstudy_min)?;

    let lf = if let Some(ref allele_path) = args.merge_alleles {
        apply_merge_alleles(lf, allele_path)?
    } else {
        lf
    };

    // Drop rows with missing required values (aligns with Python's dropna).
    let lf = drop_missing_required(lf, args.no_alleles)?;

    // Select output columns: --no-alleles omits A1/A2; --keep-maf includes FRQ.
    let lf = {
        let mut output_cols: Vec<Expr> = vec![col("SNP")];
        if !args.no_alleles {
            output_cols.push(col("A1"));
            output_cols.push(col("A2"));
        }
        output_cols.push(col("Z"));
        output_cols.push(col("N"));
        if args.keep_maf {
            output_cols.push(col("FRQ"));
        }
        lf.select(output_cols)
    };

    let mut df = lf.collect()?;

    let n_before = df.height();
    let snp_subset = ["SNP".to_string()];
    df = df
        .unique_stable(Some(snp_subset.as_slice()), UniqueKeepStrategy::First, None)
        .context("deduplicating on SNP column")?;
    let n_dup = n_before - df.height();
    if n_dup > 0 {
        println!("  Removed {} duplicate SNPs", n_dup);
    }

    let out_path = format!("{}.sumstats.gz", args.out);
    write_sumstats_gz(&out_path, &mut df)?;

    println!("Munging complete: {} SNPs -> {}", df.height(), out_path);
    Ok(())
}

// ---------------------------------------------------------------------------
// Ignore columns (--ignore)
// ---------------------------------------------------------------------------

/// Drop user-specified columns before any other processing.
///
/// `ignore_csv` is a comma-separated list of column names to drop.
/// Names are matched case-insensitively against the actual file header.
/// Unknown names are silently skipped.
fn apply_ignore(lf: LazyFrame, ignore_csv: Option<&str>) -> Result<LazyFrame> {
    let csv = match ignore_csv {
        Some(s) => s,
        None => return Ok(lf),
    };

    let to_drop: Vec<&str> = csv.split(',').map(|s| s.trim()).collect();
    let header = lf.clone().limit(0).collect()?;
    let existing: Vec<String> = header
        .get_column_names()
        .iter()
        .map(|s| s.to_string())
        .collect();

    let keep_cols: Vec<Expr> = existing
        .iter()
        .filter(|c| !to_drop.iter().any(|d| c.eq_ignore_ascii_case(d)))
        .map(|c| col(c.as_str()))
        .collect();

    Ok(lf.select(keep_cols))
}

// ---------------------------------------------------------------------------
// User column-name overrides (--snp-col, --n-col, --a1-col, etc.)
// ---------------------------------------------------------------------------

/// Apply explicit column-name overrides before the synonym map runs.
///
/// Each `--xxx-col` flag maps a user-specified column name (case-insensitive)
/// to the canonical internal name.  This lets users handle any non-standard
/// column name that isn't in `CNAME_MAP`.
fn apply_col_overrides(lf: LazyFrame, args: &MungeArgs) -> Result<LazyFrame> {
    // (user_specified_column_name, canonical_name)
    let overrides: &[(Option<&str>, &str)] = &[
        (args.snp_col.as_deref(), "SNP"),
        (args.n_col.as_deref(), "N"),
        (args.n_cas_col.as_deref(), "N_CAS"),
        (args.n_con_col.as_deref(), "N_CON"),
        (args.a1_col.as_deref(), "A1"),
        (args.a2_col.as_deref(), "A2"),
        (args.p_col.as_deref(), "P"),
        (args.frq_col.as_deref(), "FRQ"),
        (args.info_col.as_deref(), "INFO"),
    ];

    // Collect existing columns (no data load needed).
    let header = lf.clone().limit(0).collect()?;
    let existing: Vec<String> = header
        .get_column_names()
        .iter()
        .map(|s| s.to_string())
        .collect();

    let mut old_names: Vec<String> = Vec::new();
    let mut new_names: Vec<String> = Vec::new();

    for (override_opt, canonical) in overrides {
        if let Some(user_name) = override_opt {
            // Case-insensitive match against existing columns.
            if let Some(actual) = existing.iter().find(|e| e.eq_ignore_ascii_case(user_name))
                && actual != canonical
            {
                old_names.push(actual.clone());
                new_names.push(canonical.to_string());
            }
            // If not found: silently skip; the later select() will
            // raise a clear error if a required column is ultimately missing.
        }
    }

    if old_names.is_empty() {
        return Ok(lf);
    }
    Ok(lf.rename(old_names, new_names, false))
}

// ---------------------------------------------------------------------------
// Column normalisation
// ---------------------------------------------------------------------------

/// Rename input columns to their canonical names using CNAME_MAP, case-insensitively.
fn normalize_columns(lf: LazyFrame) -> Result<LazyFrame> {
    let header = lf.clone().limit(0).collect()?;
    let existing = header.get_column_names();

    let mut old_names: Vec<String> = Vec::new();
    let mut new_names: Vec<String> = Vec::new();

    for name in &existing {
        let upper = name.to_uppercase();
        if let Some(canonical) = cname_lookup(&upper) {
            // Rename only when the current name differs from the canonical.
            if *name != canonical {
                old_names.push(name.to_string());
                new_names.push(canonical.to_string());
            }
        }
    }

    if old_names.is_empty() {
        return Ok(lf);
    }
    Ok(lf.rename(old_names, new_names, false))
}

// ---------------------------------------------------------------------------
// Sample-size override
// ---------------------------------------------------------------------------

/// Apply --n / --n-cas / --n-con overrides.
///
/// Priority:
///   1. --n   → overwrite (or create) the N column with a constant value.
///   2. --n-cas + --n-con  → compute N = N-cas + N-con as a constant.
///   3. Neither → leave N column as-is (must be present in the file).
fn apply_n_override(lf: LazyFrame, args: &MungeArgs) -> Result<LazyFrame> {
    // Priority 1: --n constant.
    if let Some(n) = args.n {
        return Ok(lf.with_column(lit(n).alias("N")));
    }
    // Priority 2: --n-cas + --n-con constants.
    if let (Some(n_cas), Some(n_con)) = (args.n_cas, args.n_con) {
        return Ok(lf.with_column(lit(n_cas + n_con).alias("N")));
    }
    // Priority 3: N_CAS + N_CON per-row columns (from --n-cas-col / --n-con-col).
    let header = lf.clone().limit(0).collect()?;
    let cols = header.get_column_names();
    let has = |n: &str| cols.iter().any(|c| *c == n);
    if has("N_CAS") && has("N_CON") {
        return Ok(lf.with_column(
            (col("N_CAS").cast(DataType::Float64) + col("N_CON").cast(DataType::Float64))
                .alias("N"),
        ));
    }
    Ok(lf)
}

// ---------------------------------------------------------------------------
// Z-score derivation
// ---------------------------------------------------------------------------

/// Derive Z-score using one of four strategies (in priority order):
///   1. Z column already present → no-op.
///   2. BETA + SE present → Z = BETA / SE (lazy expression, no data load).
///   3. P present + signed column:
///      a. --signed-sumstats COLNAME,null → use that column with sign(val−null).
///      b. Auto-detect: BETA → sign(BETA), LOG_ODDS → sign(LOG_ODDS), OR → sign(OR − 1).
///      Z = sign × Φ⁻¹(1 − P/2).  Requires materialising the frame.
///   4. `--a1-inc`: A1 always increases → Z = +|Φ⁻¹(1 − P/2)| from P-value.
///   5. None of the above → pass through; final `select` will report the error.
fn derive_z(lf: LazyFrame, signed_sumstats: Option<&str>, a1_inc: bool) -> Result<LazyFrame> {
    let header = lf.clone().limit(0).collect()?;
    let cols = header.get_column_names();
    let has = |n: &str| cols.iter().any(|c| *c == n);

    if has("Z") {
        return Ok(lf);
    }

    // BETA / SE → Z = BETA / SE (lazy expression, no collect needed).
    if has("BETA") && has("SE") {
        return Ok(lf.with_column((col("BETA") / col("SE")).alias("Z")));
    }

    // --a1-inc: A1 always increases → Z = +|Φ⁻¹(1 − P/2)|.
    if a1_inc && has("P") {
        let mut df = lf
            .collect()
            .context("collecting for P→Z (a1-inc) conversion")?;
        let z_col = p_always_positive(&df).context("P→Z conversion with --a1-inc")?;
        df.with_column(z_col.into())
            .context("adding Z column (a1-inc)")?;
        return Ok(df.lazy());
    }

    // P-value + signed column → materialise and compute Z.
    if has("P") {
        let sign_info: Option<(String, f64)> = if let Some(ss) = signed_sumstats {
            // --signed-sumstats COLNAME,null_value
            let parts: Vec<&str> = ss.splitn(2, ',').collect();
            anyhow::ensure!(
                parts.len() == 2,
                "--signed-sumstats must be COLNAME,null_value (e.g. Z,0 or OR,1)"
            );
            let col_upper = parts[0].trim().to_uppercase();
            let null_val: f64 = parts[1]
                .trim()
                .parse()
                .with_context(|| format!("parsing null value from --signed-sumstats '{}'", ss))?;
            // Case-insensitive match against actual columns.
            let actual = cols
                .iter()
                .find(|c| c.to_uppercase() == col_upper)
                .map(|c| c.to_string());
            actual.map(|c| (c, null_val))
        } else if has("BETA") {
            Some(("BETA".to_string(), 0.0))
        } else if has("LOG_ODDS") {
            Some(("LOG_ODDS".to_string(), 0.0))
        } else if has("OR") {
            Some(("OR".to_string(), 1.0))
        } else {
            None
        };

        if let Some((sign_col, null_val)) = sign_info {
            let mut df = lf.collect().context("collecting for P→Z conversion")?;
            let z_col = p_and_sign_to_z(&df, &sign_col, null_val)
                .with_context(|| format!("P→Z using sign column '{}'", sign_col))?;
            df.with_column(z_col.into())
                .context("adding Z column to DataFrame")?;
            return Ok(df.lazy());
        }
    }

    Ok(lf)
}

/// Compute Z = +|Φ⁻¹(1 − P/2)| for each row (--a1-inc: sign always positive).
fn p_always_positive(df: &DataFrame) -> Result<Series> {
    let normal = Normal::new(0.0, 1.0).expect("Normal(0,1)");
    let p_series = df.column("P")?.cast(&DataType::Float64)?;
    let p_ca = p_series.f64().context("P column as f64")?;
    let z_vals: Vec<Option<f64>> = p_ca
        .into_iter()
        .map(|opt_p| {
            let p = opt_p?;
            if !p.is_finite() || p <= 0.0 || p > 1.0 {
                return None;
            }
            let p_clip = p.clamp(1e-300, 1.0 - 1e-15);
            Some(normal.inverse_cdf(1.0 - p_clip / 2.0))
        })
        .collect();
    Ok(Series::new("Z".into(), z_vals))
}

/// Compute Z = sign(value − null_val) × |Φ⁻¹(1 − P/2)| element-wise.
///
/// `null_val` is the "zero-effect" baseline of the signed column:
///   • 0.0 for BETA, LOG_ODDS, Z  →  sign = sign(value)
///   • 1.0 for OR                 →  sign = sign(OR − 1)
fn p_and_sign_to_z(df: &DataFrame, sign_col: &str, null_val: f64) -> Result<Series> {
    let normal = Normal::new(0.0, 1.0).expect("Normal(0,1)");

    let p_series = df.column("P")?.cast(&DataType::Float64)?;
    let p_ca = p_series.f64().context("P column as f64")?;

    let sign_series = df.column(sign_col)?.cast(&DataType::Float64)?;
    let sign_ca = sign_series
        .f64()
        .with_context(|| format!("'{sign_col}' column as f64"))?;

    let z_vals: Vec<Option<f64>> = p_ca
        .into_iter()
        .zip(sign_ca)
        .map(|(opt_p, opt_s)| {
            let p = opt_p?;
            let s = opt_s?;
            if !p.is_finite() || p <= 0.0 || p > 1.0 {
                return None;
            }
            let p_clip = p.clamp(1e-300, 1.0 - 1e-15);
            let abs_z = normal.inverse_cdf(1.0 - p_clip / 2.0);
            // sign = sign(s - null_val)
            let signed = s - null_val;
            let sign = if signed > 0.0 {
                1.0
            } else if signed < 0.0 {
                -1.0
            } else {
                return None; // indeterminate sign
            };
            Some(sign * abs_z)
        })
        .collect();

    Ok(Series::new("Z".into(), z_vals))
}

// ---------------------------------------------------------------------------
// SNP filtering
// ---------------------------------------------------------------------------

/// Apply MAF/N/INFO filters, optionally remove strand-ambiguous SNPs.
/// When `no_alleles` is true the strand-ambiguity check is skipped.
fn filter_snps(
    lf: LazyFrame,
    maf: f64,
    n_min: f64,
    info_min: f64,
    no_alleles: bool,
) -> Result<LazyFrame> {
    let header = lf.clone().limit(0).collect()?;
    let cols = header.get_column_names();
    let has = |n: &str| cols.iter().any(|c| *c == n);

    let mut lf = lf;

    // Python default: if N exists and --n-min is unset, use 90th percentile / 1.5.
    let n_min = if n_min == 0.0 && has("N") {
        compute_default_n_min(&lf)?.unwrap_or(0.0)
    } else {
        n_min
    };

    // Filter by minimum sample size.
    if n_min > 0.0 && has("N") {
        lf = lf.filter(col("N").cast(DataType::Float64).gt_eq(lit(n_min)));
    }

    // Filter by MAF: keep FRQ in [maf, 1−maf].
    if maf > 0.0 && has("FRQ") {
        lf = lf.filter(
            col("FRQ")
                .cast(DataType::Float64)
                .gt_eq(lit(maf))
                .and(col("FRQ").cast(DataType::Float64).lt_eq(lit(1.0 - maf))),
        );
    }

    // Filter by INFO score.
    if info_min > 0.0 && has("INFO") {
        lf = lf.filter(col("INFO").cast(DataType::Float64).gt_eq(lit(info_min)));
    }

    // Remove invalid or strand-ambiguous SNPs unless --no-alleles.
    if !no_alleles && has("A1") && has("A2") {
        lf = lf.with_columns([
            col("A1").str().to_uppercase().alias("A1"),
            col("A2").str().to_uppercase().alias("A2"),
        ]);

        let valid_base = |c: &str| {
            col(c)
                .eq(lit("A"))
                .or(col(c).eq(lit("C")))
                .or(col(c).eq(lit("G")))
                .or(col(c).eq(lit("T")))
        };
        let not_same = col("A1").neq(col("A2"));
        let not_ambig = col("A1")
            .eq(lit("A"))
            .and(col("A2").eq(lit("T")))
            .or(col("A1").eq(lit("T")).and(col("A2").eq(lit("A"))))
            .or(col("A1").eq(lit("C")).and(col("A2").eq(lit("G"))))
            .or(col("A1").eq(lit("G")).and(col("A2").eq(lit("C"))))
            .not();

        let valid = valid_base("A1")
            .and(valid_base("A2"))
            .and(not_same)
            .and(not_ambig);
        lf = lf.filter(valid);
    }

    Ok(lf)
}

fn compute_default_n_min(lf: &LazyFrame) -> Result<Option<f64>> {
    let n_df = lf
        .clone()
        .select([col("N").cast(DataType::Float64)])
        .collect()?;
    let n_series = n_df.column("N")?.f64()?;
    let mut vals: Vec<f64> = n_series.into_iter().flatten().collect();
    if vals.is_empty() {
        return Ok(None);
    }
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let pos = (vals.len() as f64 - 1.0) * 0.9;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    let p90 = if lo == hi {
        vals[lo]
    } else {
        let frac = pos - lo as f64;
        vals[lo] + frac * (vals[hi] - vals[lo])
    };
    Ok(Some(p90 / 1.5))
}

fn drop_missing_required(lf: LazyFrame, no_alleles: bool) -> Result<LazyFrame> {
    let mut predicate = col("SNP")
        .is_not_null()
        .and(col("Z").is_not_null())
        .and(col("N").is_not_null());
    if !no_alleles {
        predicate = predicate
            .and(col("A1").is_not_null())
            .and(col("A2").is_not_null());
    }
    Ok(lf.filter(predicate))
}

fn complement_expr(col_name: &str) -> Expr {
    when(col(col_name).eq(lit("A")))
        .then(lit("T"))
        .when(col(col_name).eq(lit("T")))
        .then(lit("A"))
        .when(col(col_name).eq(lit("C")))
        .then(lit("G"))
        .when(col(col_name).eq(lit("G")))
        .then(lit("C"))
        .otherwise(lit(""))
}

fn apply_merge_alleles(lf: LazyFrame, allele_path: &str) -> Result<LazyFrame> {
    let alleles = parse::scan_sumstats(allele_path)?;
    let alleles = normalize_columns(alleles)?;
    let alleles = alleles.select([
        col("SNP"),
        col("A1").str().to_uppercase().alias("A1_M"),
        col("A2").str().to_uppercase().alias("A2_M"),
    ]);

    let merged = lf.join(
        alleles,
        [col("SNP")],
        [col("SNP")],
        JoinArgs::new(JoinType::Inner),
    );

    let a1 = col("A1");
    let a2 = col("A2");
    let m1 = col("A1_M");
    let m2 = col("A2_M");
    let m1c = complement_expr("A1_M");
    let m2c = complement_expr("A2_M");

    let matches = a1
        .clone()
        .eq(m1.clone())
        .and(a2.clone().eq(m2.clone()))
        .or(a1.clone().eq(m1c.clone()).and(a2.clone().eq(m2c.clone())))
        .or(a1.clone().eq(m2.clone()).and(a2.clone().eq(m1.clone())))
        .or(a1.eq(m2c).and(a2.eq(m1c)));

    let merged = merged.filter(matches);
    let header = merged.clone().limit(0).collect()?;
    let keep: Vec<Expr> = header
        .get_column_names()
        .iter()
        .filter(|c| c.as_str() != "A1_M" && c.as_str() != "A2_M")
        .map(|c| col(c.as_str()))
        .collect();
    Ok(merged.select(keep))
}

// ---------------------------------------------------------------------------
// Info-list helper (--info-list)
// ---------------------------------------------------------------------------

/// Compute the mean of multiple INFO columns and store the result as "INFO".
///
/// Useful when a summary-stats file provides per-population imputation scores
/// (e.g. INFO_EUR, INFO_EAS) and the analyst wants to filter on their mean.
fn apply_info_list(lf: LazyFrame, info_list: Option<&str>) -> Result<LazyFrame> {
    let list = match info_list {
        Some(s) => s,
        None => return Ok(lf),
    };
    let requested: Vec<&str> = list
        .split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .collect();
    if requested.is_empty() {
        return Ok(lf);
    }

    // Match requested names case-insensitively against the actual column names.
    let header = lf.clone().limit(0).collect()?;
    let existing = header.get_column_names();
    let matched: Vec<String> = requested
        .iter()
        .filter_map(|name| {
            existing
                .iter()
                .find(|e| e.to_uppercase() == name.to_uppercase())
                .map(|e| e.to_string())
        })
        .collect();

    anyhow::ensure!(
        !matched.is_empty(),
        "--info-list: none of the listed columns ({}) were found in the file",
        requested.join(", ")
    );

    let n = matched.len() as f64;
    let sum_expr = matched.iter().skip(1).fold(
        col(matched[0].as_str()).cast(DataType::Float64),
        |acc, c| acc + col(c.as_str()).cast(DataType::Float64),
    );
    println!(
        "  --info-list: computing mean INFO from {} columns: {}",
        matched.len(),
        matched.join(", ")
    );
    Ok(lf.with_column((sum_expr / lit(n)).alias("INFO")))
}

// ---------------------------------------------------------------------------
// Nstudy filter helper (--nstudy / --nstudy-min)
// ---------------------------------------------------------------------------

/// Filter SNPs by the minimum number of studies they appear in.
///
/// `nstudy` names the column holding the per-SNP study count;
/// `nstudy_min` is the minimum required value (inclusive).
fn apply_nstudy_filter(
    lf: LazyFrame,
    nstudy: Option<&str>,
    nstudy_min: Option<u64>,
) -> Result<LazyFrame> {
    let (col_name, min_val) = match (nstudy, nstudy_min) {
        (Some(c), Some(m)) => (c, m),
        _ => return Ok(lf),
    };

    let header = lf.clone().limit(0).collect()?;
    let existing = header.get_column_names();
    if let Some(actual) = existing
        .iter()
        .find(|c| c.to_uppercase() == col_name.to_uppercase())
    {
        let actual_owned = actual.to_string();
        println!(
            "  --nstudy-min {}: filtering on column '{}'",
            min_val, actual_owned
        );
        return Ok(lf.filter(
            col(actual_owned.as_str())
                .cast(DataType::Float64)
                .gt_eq(lit(min_val as f64)),
        ));
    }
    println!(
        "  Warning: --nstudy column '{}' not found in file; skipping filter",
        col_name
    );
    Ok(lf)
}

// ---------------------------------------------------------------------------
// Output writer
// ---------------------------------------------------------------------------

/// Write a DataFrame as a gzip-compressed tab-separated file.
fn write_sumstats_gz(path: &str, df: &mut DataFrame) -> Result<()> {
    use flate2::Compression;
    use flate2::write::GzEncoder;

    let file = File::create(path).with_context(|| format!("creating '{}'", path))?;
    let gz = GzEncoder::new(BufWriter::new(file), Compression::fast());
    CsvWriter::new(gz)
        .with_separator(b'\t')
        .finish(df)
        .with_context(|| format!("writing '{}'", path))?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    /// cname_lookup maps known synonyms to their canonical names.
    #[test]
    fn test_cname_lookup_standard() {
        assert_eq!(cname_lookup("SNP"), Some("SNP"));
        assert_eq!(cname_lookup("RSID"), Some("SNP"));
        assert_eq!(cname_lookup("MARKERNAME"), Some("SNP"));
        assert_eq!(cname_lookup("P"), Some("P"));
        assert_eq!(cname_lookup("PVALUE"), Some("P"));
        assert_eq!(cname_lookup("BETA"), Some("BETA"));
        assert_eq!(cname_lookup("B"), Some("BETA"));
        assert_eq!(cname_lookup("SE"), Some("SE"));
        assert_eq!(cname_lookup("STDERR"), Some("SE"));
        assert_eq!(cname_lookup("MAF"), Some("FRQ"));
        assert_eq!(cname_lookup("EAF"), Some("FRQ"));
        assert_eq!(cname_lookup("FRQ_U"), Some("FRQ"));
        assert_eq!(cname_lookup("INFO"), Some("INFO"));
        assert_eq!(cname_lookup("IMPINFO"), Some("INFO"));
        assert_eq!(cname_lookup("Z"), Some("Z"));
        assert_eq!(cname_lookup("ZSCORE"), Some("Z"));
        assert_eq!(cname_lookup("A1"), Some("A1"));
        assert_eq!(cname_lookup("EFFECT_ALLELE"), Some("A1"));
        assert_eq!(cname_lookup("A2"), Some("A2"));
        assert_eq!(cname_lookup("OTHER_ALLELE"), Some("A2"));
        assert_eq!(cname_lookup("N"), Some("N"));
        assert_eq!(cname_lookup("WEIGHT"), Some("N")); // METAL convention
        assert_eq!(cname_lookup("OR"), Some("OR"));
    }

    /// cname_lookup returns None for unknown names.
    #[test]
    fn test_cname_lookup_unknown() {
        assert_eq!(cname_lookup("FOOBAR"), None);
        assert_eq!(cname_lookup(""), None);
        assert_eq!(cname_lookup("CHROM"), None);
        assert_eq!(cname_lookup("EFFECT_SIZE"), None);
    }

    /// cname_lookup is case-sensitive (keys must be UPPER); lowercase misses.
    /// Mirrors Python's clean_header uppercasing before lookup.
    #[test]
    fn test_cname_lookup_case_sensitive() {
        assert_eq!(cname_lookup("snp"), None); // lowercase → no match
        assert_eq!(cname_lookup("SNP"), Some("SNP"));
        assert_eq!(cname_lookup("Zscore"), None);
        assert_eq!(cname_lookup("ZSCORE"), Some("Z"));
    }

    fn read_gz(path: &str) -> String {
        use flate2::read::GzDecoder;
        let file = std::fs::File::open(path).expect("open gz");
        let mut decoder = GzDecoder::new(file);
        let mut s = String::new();
        decoder.read_to_string(&mut s).expect("read gz");
        s
    }

    fn base_args(sumstats: &str, out: &str) -> MungeArgs {
        MungeArgs {
            sumstats: sumstats.to_string(),
            out: out.to_string(),
            merge_alleles: None,
            n_min: 0.0,
            maf: 0.01,
            info_min: 0.9,
            keep_mhc: true,
            n: None,
            n_cas: None,
            n_con: None,
            snp_col: None,
            n_col: None,
            n_cas_col: None,
            n_con_col: None,
            a1_col: None,
            a2_col: None,
            p_col: None,
            frq_col: None,
            info_col: None,
            signed_sumstats: None,
            ignore: None,
            keep_maf: false,
            a1_inc: false,
            no_alleles: false,
            info_list: None,
            nstudy: None,
            nstudy_min: None,
        }
    }

    #[test]
    fn test_merge_alleles_filters_mismatch() {
        let dir = tempfile::tempdir().unwrap();
        let sumstats = dir.path().join("sumstats.txt");
        let merge = dir.path().join("merge.txt");
        let out = dir.path().join("out");

        std::fs::write(
            &sumstats,
            "SNP A1 A2 Z N\nrs1 A G 1.0 100\nrs2 C T 2.0 100\n",
        )
        .unwrap();
        std::fs::write(&merge, "SNP A1 A2\nrs1 A G\nrs2 G C\n").unwrap();

        let mut args = base_args(sumstats.to_str().unwrap(), out.to_str().unwrap());
        args.merge_alleles = Some(merge.to_str().unwrap().to_string());
        run(args).unwrap();

        let content = read_gz(&format!("{}.sumstats.gz", out.to_str().unwrap()));
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 2, "expected header + 1 row");
        assert!(lines[1].starts_with("rs1\t"));
    }

    #[test]
    fn test_invalid_alleles_removed() {
        let dir = tempfile::tempdir().unwrap();
        let sumstats = dir.path().join("sumstats.txt");
        let out = dir.path().join("out");
        std::fs::write(
            &sumstats,
            "SNP A1 A2 Z N\nrs1 A G 1.0 100\nrs2 I D 2.0 100\n",
        )
        .unwrap();
        let args = base_args(sumstats.to_str().unwrap(), out.to_str().unwrap());
        run(args).unwrap();

        let content = read_gz(&format!("{}.sumstats.gz", out.to_str().unwrap()));
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 2, "expected header + 1 row");
        assert!(lines[1].starts_with("rs1\t"));
    }

    #[test]
    fn test_drop_missing_values() {
        let dir = tempfile::tempdir().unwrap();
        let sumstats = dir.path().join("sumstats.txt");
        let out = dir.path().join("out");
        std::fs::write(&sumstats, "SNP A1 A2 Z N\nrs1 A G 1.0 100\nrs2 C T . 100\n").unwrap();
        let args = base_args(sumstats.to_str().unwrap(), out.to_str().unwrap());
        run(args).unwrap();

        let content = read_gz(&format!("{}.sumstats.gz", out.to_str().unwrap()));
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 2, "expected header + 1 row");
        assert!(lines[1].starts_with("rs1\t"));
    }

    #[test]
    fn test_default_n_min_filters_low_n() {
        let dir = tempfile::tempdir().unwrap();
        let sumstats = dir.path().join("sumstats.txt");
        let out = dir.path().join("out");
        std::fs::write(
            &sumstats,
            "SNP A1 A2 Z N\nrs1 A G 1.0 10\nrs2 A G 1.0 20\nrs3 A G 1.0 30\nrs4 A G 1.0 40\nrs5 A G 1.0 100\n",
        )
        .unwrap();
        let args = base_args(sumstats.to_str().unwrap(), out.to_str().unwrap());
        run(args).unwrap();

        let content = read_gz(&format!("{}.sumstats.gz", out.to_str().unwrap()));
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 2, "expected header + 1 row");
        assert!(lines[1].starts_with("rs5\t"));
    }
}
