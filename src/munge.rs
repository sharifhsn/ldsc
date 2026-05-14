/// Summary statistics munging — pure-Rust pipeline (csv + flate2), no polars.
use anyhow::{Context, Result};
use statrs::function::erf::erfc_inv;

use crate::cli::MungeArgs;
use crate::frame::{self, Column, Frame};
use crate::parse;

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
    ("WEIGHT", "N"),
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

pub fn run(args: MungeArgs) -> Result<()> {
    let mut args = args;
    if args.no_alleles && args.merge_alleles.is_some() {
        anyhow::bail!("--no-alleles and --merge-alleles are not compatible");
    }
    if args.daner && args.daner_n {
        anyhow::bail!(
            "--daner and --daner-n are not compatible. Use --daner for sample size from \
             FRQ_A/FRQ_U headers, use --daner-n for values from Nca/Nco columns"
        );
    }

    let mut f = parse::scan_sumstats(&args.sumstats)?;
    apply_ignore(&mut f, args.ignore.as_deref())?;
    apply_daner_overrides(&mut f, &mut args)?;
    apply_col_overrides(&mut f, &args)?;
    normalize_columns(&mut f)?;
    apply_info_list(&mut f, args.info_list.as_deref())?;
    apply_n_override(&mut f, &args)?;
    filter_pvals(&mut f)?;
    derive_z(&mut f, args.signed_sumstats.as_deref(), args.a1_inc)?;
    filter_snps(&mut f, args.maf, args.n_min, args.info_min, args.no_alleles)?;
    apply_nstudy_filter(&mut f, args.nstudy.as_deref(), args.nstudy_min)?;
    if let Some(allele_path) = args.merge_alleles.clone() {
        apply_merge_alleles(&mut f, &allele_path)?;
    }
    drop_missing_required(&mut f)?;

    // Select output columns: --no-alleles omits A1/A2; --keep-maf includes FRQ.
    let mut out_cols: Vec<&str> = vec!["SNP"];
    if !args.no_alleles {
        out_cols.push("A1");
        out_cols.push("A2");
    }
    out_cols.push("Z");
    out_cols.push("N");
    if args.keep_maf {
        out_cols.push("FRQ");
    }
    let mut f = f.select(&out_cols)?;

    let n_before = f.height();
    f = f.unique_first_on("SNP")?;
    let n_dup = n_before - f.height();
    if n_dup > 0 {
        println!("  Removed {} duplicate SNPs", n_dup);
    }

    let out_path = format!("{}.sumstats.gz", args.out);
    frame::write_tsv(&out_path, &f).with_context(|| format!("writing '{}'", out_path))?;

    println!("Munging complete: {} SNPs -> {}", f.height(), out_path);
    Ok(())
}

/// Drop user-specified columns before any other processing.
fn apply_ignore(f: &mut Frame, ignore_csv: Option<&str>) -> Result<()> {
    let csv = match ignore_csv {
        Some(s) => s,
        None => return Ok(()),
    };
    let to_drop: Vec<&str> = csv.split(',').map(|s| s.trim()).collect();
    let drop_names: Vec<String> = f
        .column_names_owned()
        .into_iter()
        .filter(|c| to_drop.iter().any(|d| c.eq_ignore_ascii_case(d)))
        .collect();
    for n in drop_names {
        f.drop_column(&n)?;
    }
    Ok(())
}

fn apply_daner_overrides(f: &mut Frame, args: &mut MungeArgs) -> Result<()> {
    if !args.daner && !args.daner_n {
        return Ok(());
    }
    let cols = f.column_names_owned();
    let find_prefix = |prefix: &str| cols.iter().find(|c| c.starts_with(prefix)).cloned();

    let frq_u = find_prefix("FRQ_U_")
        .with_context(|| "Could not find FRQ_U_* column expected for daner format")?;

    // Drop any other FRQ-synonyms so the daner FRQ_U_* column wins.
    let drop_frq: Vec<String> = cols
        .iter()
        .filter(|name| cname_lookup(&name.to_uppercase()) == Some("FRQ") && name.as_str() != frq_u)
        .cloned()
        .collect();
    for d in drop_frq {
        f.drop_column(&d)?;
    }

    if frq_u != "FRQ" {
        f.rename(&frq_u, "FRQ")?;
    }

    if args.daner {
        let frq_a = find_prefix("FRQ_A_")
            .with_context(|| "Could not find FRQ_A_* column expected for daner format")?;
        let n_con: f64 = frq_u
            .strip_prefix("FRQ_U_")
            .context("FRQ_U_* column missing numeric suffix")?
            .parse()
            .context("Parsing N_con from FRQ_U_* suffix")?;
        let n_cas: f64 = frq_a
            .strip_prefix("FRQ_A_")
            .context("FRQ_A_* column missing numeric suffix")?
            .parse()
            .context("Parsing N_cas from FRQ_A_* suffix")?;
        println!(
            "  --daner: inferred N_cas = {} and N_con = {} from FRQ_[A/U] headers",
            n_cas, n_con
        );
        args.n_cas = Some(n_cas);
        args.n_con = Some(n_con);
    }

    if args.daner_n {
        let nca = cols
            .iter()
            .find(|c| c.as_str() == "Nca")
            .cloned()
            .context("Could not find Nca column expected for daner-n format")?;
        let nco = cols
            .iter()
            .find(|c| c.as_str() == "Nco")
            .cloned()
            .context("Could not find Nco column expected for daner-n format")?;
        if nca != "N_CAS" {
            f.rename(&nca, "N_CAS")?;
        }
        if nco != "N_CON" {
            f.rename(&nco, "N_CON")?;
        }
    }
    Ok(())
}

fn apply_col_overrides(f: &mut Frame, args: &MungeArgs) -> Result<()> {
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
    let existing = f.column_names_owned();
    for (override_opt, canonical) in overrides {
        if let Some(user_name) = override_opt
            && let Some(actual) = existing.iter().find(|e| e.eq_ignore_ascii_case(user_name))
            && actual != *canonical
        {
            f.rename(actual, canonical)?;
        }
    }
    Ok(())
}

fn normalize_columns(f: &mut Frame) -> Result<()> {
    let existing = f.column_names_owned();
    for name in &existing {
        let upper = name.to_uppercase();
        if let Some(canonical) = cname_lookup(&upper)
            && name != canonical
        {
            f.rename(name, canonical)?;
        }
    }
    Ok(())
}

fn apply_info_list(f: &mut Frame, info_list: Option<&str>) -> Result<()> {
    let list = match info_list {
        Some(s) => s,
        None => return Ok(()),
    };
    let requested: Vec<&str> = list
        .split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .collect();
    if requested.is_empty() {
        return Ok(());
    }
    let existing = f.column_names_owned();
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

    println!(
        "  --info-list: computing mean INFO from {} columns: {}",
        matched.len(),
        matched.join(", ")
    );
    let n = matched.len() as f64;
    let casted: Vec<Vec<f64>> = matched
        .iter()
        .map(|c| {
            f.column(c)
                .and_then(|c| c.cast_to_f64())
                .map(|c| c.as_f64().unwrap().to_vec())
        })
        .collect::<Result<Vec<_>>>()?;
    let h = f.height();
    let mut info = vec![0.0; h];
    for col in &casted {
        for (acc, v) in info.iter_mut().zip(col.iter()) {
            *acc += *v;
        }
    }
    for v in info.iter_mut() {
        *v /= n;
    }
    if f.has_column("INFO") {
        f.replace_column("INFO", Column::F64(info))?;
    } else {
        f.add_column("INFO".into(), Column::F64(info))?;
    }
    Ok(())
}

fn apply_n_override(f: &mut Frame, args: &MungeArgs) -> Result<()> {
    let h = f.height();
    if let Some(n) = args.n {
        let v = vec![n; h];
        if f.has_column("N") {
            f.replace_column("N", Column::F64(v))?;
        } else {
            f.add_column("N".into(), Column::F64(v))?;
        }
        return Ok(());
    }
    if let (Some(n_cas), Some(n_con)) = (args.n_cas, args.n_con) {
        let v = vec![n_cas + n_con; h];
        if f.has_column("N") {
            f.replace_column("N", Column::F64(v))?;
        } else {
            f.add_column("N".into(), Column::F64(v))?;
        }
        return Ok(());
    }
    if f.has_column("N_CAS") && f.has_column("N_CON") {
        let n_cas_col = f.column("N_CAS")?.cast_to_f64()?;
        let n_con_col = f.column("N_CON")?.cast_to_f64()?;
        let n_cas = n_cas_col.as_f64()?;
        let n_con = n_con_col.as_f64()?;
        // N_total = N_CAS + N_CON
        let n_total: Vec<f64> = n_cas.iter().zip(n_con.iter()).map(|(a, b)| a + b).collect();
        let max_n = n_total
            .iter()
            .copied()
            .filter(|v| !v.is_nan())
            .fold(f64::NEG_INFINITY, f64::max);
        // P_max = mean of N_CAS/N_total at rows where N_total == max_n.
        let mut sum_p = 0.0;
        let mut count_p: u64 = 0;
        for (i, &nt) in n_total.iter().enumerate() {
            if nt == max_n && nt > 0.0 {
                sum_p += n_cas[i] / nt;
                count_p += 1;
            }
        }
        let p_max = if count_p > 0 {
            sum_p / count_p as f64
        } else {
            1.0
        };
        let n_col: Vec<f64> = n_cas.iter().map(|nc| nc / p_max).collect();
        let _ = f.drop_column("N_CAS").ok();
        let _ = f.drop_column("N_CON").ok();
        if f.has_column("N") {
            f.replace_column("N", Column::F64(n_col))?;
        } else {
            f.add_column("N".into(), Column::F64(n_col))?;
        }
    }
    Ok(())
}

fn filter_pvals(f: &mut Frame) -> Result<()> {
    if !f.has_column("P") {
        return Ok(());
    }
    let casted = f.column("P")?.cast_to_f64()?;
    let p = casted.as_f64()?;
    let mask: Vec<bool> = p
        .iter()
        .map(|v| !v.is_nan() && *v > 0.0 && *v <= 1.0)
        .collect();
    *f = f.filter_rows(&mask)?;
    Ok(())
}

/// Derive Z-score using one of four strategies.
fn derive_z(f: &mut Frame, signed_sumstats: Option<&str>, a1_inc: bool) -> Result<()> {
    let cols = f.column_names_owned();
    let has = |n: &str| cols.iter().any(|c| c == n);

    if has("P") {
        if a1_inc {
            // Z = +|Φ⁻¹(1 − P/2)|
            let p = f.column("P")?.cast_to_f64()?;
            let p = p.as_f64()?;
            let z: Vec<f64> = p.iter().map(|&pv| p_to_abs_z(pv)).collect();
            insert_z(f, z)?;
            return Ok(());
        }
        let sign_info: Option<(String, f64)> = if let Some(ss) = signed_sumstats {
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
            cols.iter()
                .find(|c| c.to_uppercase() == col_upper)
                .map(|c| (c.clone(), null_val))
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
            let p = f.column("P")?.cast_to_f64()?;
            let p = p.as_f64()?;
            let s = f.column(&sign_col)?.cast_to_f64()?;
            let s = s.as_f64()?;
            let z: Vec<f64> = p
                .iter()
                .zip(s.iter())
                .map(|(&pv, &sv)| p_signed_to_z(pv, sv, null_val))
                .collect();
            insert_z(f, z)?;
            return Ok(());
        }
    }
    if has("Z") {
        return Ok(());
    }
    if has("BETA") && has("SE") {
        let beta = f.column("BETA")?.cast_to_f64()?;
        let se = f.column("SE")?.cast_to_f64()?;
        let beta = beta.as_f64()?;
        let se = se.as_f64()?;
        let z: Vec<f64> = beta.iter().zip(se.iter()).map(|(b, s)| b / s).collect();
        insert_z(f, z)?;
    }
    Ok(())
}

fn insert_z(f: &mut Frame, z: Vec<f64>) -> Result<()> {
    if f.has_column("Z") {
        f.replace_column("Z", Column::F64(z))
    } else {
        f.add_column("Z".into(), Column::F64(z))
    }
}

fn p_to_abs_z(p: f64) -> f64 {
    if !p.is_finite() || p <= 0.0 || p > 1.0 {
        return f64::NAN;
    }
    let p_clip = p.clamp(1e-300, 1.0);
    (2.0f64).sqrt() * erfc_inv(p_clip)
}

fn p_signed_to_z(p: f64, signed: f64, null_val: f64) -> f64 {
    if !p.is_finite() || p <= 0.0 || p > 1.0 || signed.is_nan() {
        return f64::NAN;
    }
    let p_clip = p.clamp(1e-300, 1.0);
    let abs_z = (2.0f64).sqrt() * erfc_inv(p_clip);
    let s = signed - null_val;
    if s > 0.0 {
        abs_z
    } else if s < 0.0 {
        -abs_z
    } else {
        f64::NAN
    }
}

fn filter_snps(f: &mut Frame, maf: f64, n_min: f64, info_min: f64, no_alleles: bool) -> Result<()> {
    let h = f.height();
    let mut mask = vec![true; h];

    // Default n_min: 90th-pct N / 1.5 if unset and N present.
    let n_min = if n_min == 0.0 && f.has_column("N") {
        compute_default_n_min(f)?.unwrap_or(0.0)
    } else {
        n_min
    };

    // NaN-aware: a NaN value fails every comparison, so we drop it.
    // (Spelling out `is_nan() || x < min` is equivalent to `!(x >= min)` but
    // avoids the negated-partial-ord clippy lint and reads more clearly.)
    if n_min > 0.0 && f.has_column("N") {
        let n = f.column("N")?.cast_to_f64()?;
        let n = n.as_f64()?;
        for i in 0..h {
            if n[i].is_nan() || n[i] < n_min {
                mask[i] = false;
            }
        }
    }

    if f.has_column("FRQ") {
        let frq = f.column("FRQ")?.cast_to_f64()?;
        let frq = frq.as_f64()?;
        for i in 0..h {
            let v = frq[i];
            if v.is_nan() || !(0.0..=1.0).contains(&v) {
                mask[i] = false;
                continue;
            }
            let m = if v < 0.5 { v } else { 1.0 - v };
            if m.is_nan() || m <= maf {
                mask[i] = false;
            }
        }
    }

    if f.has_column("INFO") {
        let info = f.column("INFO")?.cast_to_f64()?;
        let info = info.as_f64()?;
        for i in 0..h {
            if info[i].is_nan() || info[i] < info_min {
                mask[i] = false;
            }
        }
    }

    if !no_alleles && f.has_column("A1") && f.has_column("A2") {
        f.uppercase_str_column("A1")?;
        f.uppercase_str_column("A2")?;
        let a1 = f.column("A1")?.as_str()?.to_vec();
        let a2 = f.column("A2")?.as_str()?.to_vec();
        for i in 0..h {
            let a1 = a1[i].as_deref().unwrap_or("");
            let a2 = a2[i].as_deref().unwrap_or("");
            if !valid_unambig_pair(a1, a2) {
                mask[i] = false;
            }
        }
    }

    *f = f.filter_rows(&mask)?;
    Ok(())
}

fn valid_unambig_pair(a1: &str, a2: &str) -> bool {
    let valid = |s: &str| matches!(s, "A" | "C" | "G" | "T");
    if !valid(a1) || !valid(a2) || a1 == a2 {
        return false;
    }
    // strand-ambiguous: A-T, T-A, C-G, G-C
    !matches!((a1, a2), ("A", "T") | ("T", "A") | ("C", "G") | ("G", "C"))
}

fn compute_default_n_min(f: &Frame) -> Result<Option<f64>> {
    let n = f.column("N")?.cast_to_f64()?;
    let mut vals: Vec<f64> = n
        .as_f64()?
        .iter()
        .copied()
        .filter(|v| !v.is_nan())
        .collect();
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

fn drop_missing_required(f: &mut Frame) -> Result<()> {
    let h = f.height();
    let mut mask = vec![true; h];
    let names = f.column_names_owned();
    for name in names {
        if name == "INFO" {
            continue;
        }
        match f.column(&name)? {
            Column::F64(v) => {
                for i in 0..h {
                    if v[i].is_nan() {
                        mask[i] = false;
                    }
                }
            }
            Column::Str(v) => {
                for i in 0..h {
                    if v[i].is_none() || v[i].as_deref().map(|s| s.is_empty()).unwrap_or(true) {
                        mask[i] = false;
                    }
                }
            }
        }
    }
    *f = f.filter_rows(&mask)?;
    Ok(())
}

fn apply_nstudy_filter(f: &mut Frame, nstudy: Option<&str>, nstudy_min: Option<u64>) -> Result<()> {
    let (col_name, min_val) = match (nstudy, nstudy_min) {
        (Some(c), Some(m)) => (c, m),
        _ => return Ok(()),
    };
    let existing = f.column_names_owned();
    if let Some(actual) = existing
        .iter()
        .find(|c| c.to_uppercase() == col_name.to_uppercase())
    {
        let actual = actual.clone();
        println!(
            "  --nstudy-min {}: filtering on column '{}'",
            min_val, actual
        );
        let casted = f.column(&actual)?.cast_to_f64()?;
        let v = casted.as_f64()?;
        let mask: Vec<bool> = v
            .iter()
            .map(|&x| !x.is_nan() && x >= min_val as f64)
            .collect();
        *f = f.filter_rows(&mask)?;
    } else {
        println!(
            "  Warning: --nstudy column '{}' not found in file; skipping filter",
            col_name
        );
    }
    Ok(())
}

fn complement(c: &str) -> &'static str {
    match c {
        "A" => "T",
        "T" => "A",
        "C" => "G",
        "G" => "C",
        _ => "",
    }
}

fn apply_merge_alleles(f: &mut Frame, allele_path: &str) -> Result<()> {
    let mut alleles = parse::scan_sumstats(allele_path)?;
    normalize_columns(&mut alleles)?;
    alleles.uppercase_str_column("A1")?;
    alleles.uppercase_str_column("A2")?;
    alleles.rename("A1", "A1_M")?;
    alleles.rename("A2", "A2_M")?;
    let alleles = alleles.select(&["SNP", "A1_M", "A2_M"])?;

    let merged = f.join_inner_on(&alleles, "SNP")?;

    // Allele-match filter.
    let h = merged.height();
    let a1 = merged.column("A1")?.as_str()?;
    let a2 = merged.column("A2")?.as_str()?;
    let m1 = merged.column("A1_M")?.as_str()?;
    let m2 = merged.column("A2_M")?.as_str()?;
    let mut mask = vec![false; h];
    for i in 0..h {
        let a1 = a1[i].as_deref().unwrap_or("");
        let a2 = a2[i].as_deref().unwrap_or("");
        let m1 = m1[i].as_deref().unwrap_or("");
        let m2 = m2[i].as_deref().unwrap_or("");
        let m1c = complement(m1);
        let m2c = complement(m2);
        if (a1 == m1 && a2 == m2)
            || (a1 == m1c && a2 == m2c)
            || (a1 == m2 && a2 == m1)
            || (a1 == m2c && a2 == m1c)
        {
            mask[i] = true;
        }
    }
    let mut out = merged.filter_rows(&mask)?;
    out.drop_column("A1_M")?;
    out.drop_column("A2_M")?;
    *f = out;
    Ok(())
}
