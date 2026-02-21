/// File parsing utilities: Polars LazyFrame readers for `.sumstats` and `.ldscore` files.
use anyhow::{Context, Result};
use ndarray::Array2;
use polars::prelude::*;
use std::fs::File;
use std::path::PathBuf;
use std::sync::{Mutex, OnceLock};

use bzip2_rs::DecoderReader;
use tempfile::Builder as TempBuilder;

static BZ2_TEMPFILES: OnceLock<Mutex<Vec<tempfile::TempPath>>> = OnceLock::new();

fn maybe_decompress_bz2(path: &str) -> Result<PathBuf> {
    if !path.ends_with(".bz2") {
        return Ok(PathBuf::from(path));
    }

    let input = File::open(path).with_context(|| format!("opening bz2 file '{}'", path))?;
    let mut decoder = DecoderReader::new(input);
    let mut tmp = TempBuilder::new()
        .prefix("ldsc_bz2_")
        .suffix(".tmp")
        .tempfile()
        .context("creating temp file for bz2 decompression")?;
    std::io::copy(&mut decoder, &mut tmp)
        .with_context(|| format!("decompressing bz2 file '{}'", path))?;

    let temp_path = tmp.into_temp_path();
    let temp_buf = temp_path.to_path_buf();
    BZ2_TEMPFILES
        .get_or_init(|| Mutex::new(Vec::new()))
        .lock()
        .expect("bz2 tempfiles mutex poisoned")
        .push(temp_path);

    Ok(temp_buf)
}

// ---------------------------------------------------------------------------
// Readers
// ---------------------------------------------------------------------------

/// Scan a `.sumstats[.gz|.bz2]` file as a Polars LazyFrame.
pub fn scan_sumstats(path: &str) -> Result<LazyFrame> {
    let resolved = maybe_decompress_bz2(path)?;
    let resolved_str = resolved.to_string_lossy();
    let lf = LazyCsvReader::new(resolved_str.as_ref().into())
        .with_separator(b'\t')
        .with_has_header(true)
        .finish()?;
    Ok(lf)
}

/// Scan a `.ldscore[.gz|.bz2]` file as a Polars LazyFrame.
pub fn scan_ldscore(path: &str) -> Result<LazyFrame> {
    let resolved = maybe_decompress_bz2(path)?;
    let resolved_str = resolved.to_string_lossy();
    let lf = LazyCsvReader::new(resolved_str.as_ref().into())
        .with_separator(b'\t')
        .with_has_header(true)
        .finish()?;
    Ok(lf)
}

/// Build a per-chromosome path from `prefix`, `chr`, and `suffix`.
///
/// If `prefix` contains `@`, it is used as a placeholder for the chromosome number
/// (matching Python LDSC behaviour: `ld/chr@_scores` → `ld/chr22_scores.l2.ldscore.gz`).
/// Otherwise the chromosome number is appended after the prefix.
fn make_chr_path(prefix: &str, chr: u8, suffix: &str) -> String {
    if prefix.contains('@') {
        format!("{}{}", prefix.replace('@', &chr.to_string()), suffix)
    } else {
        format!("{}{}{}", prefix, chr, suffix)
    }
}

/// Return the chromosomes (1–22) for which `{prefix}{chr}{suffix}` exists.
pub fn get_present_chrs(prefix: &str, suffix: &str) -> Vec<u8> {
    (1u8..=22)
        .filter(|chr| std::path::Path::new(&make_chr_path(prefix, *chr, suffix)).exists())
        .collect()
}

/// Sum M values from per-chromosome `.M` files (one integer per file).
pub fn read_m_total(prefix: &str, suffix: &str) -> Result<f64> {
    let chrs = get_present_chrs(prefix, suffix);
    anyhow::ensure!(
        !chrs.is_empty(),
        "No chromosome M files found for prefix '{}' suffix '{}'",
        prefix,
        suffix
    );

    let mut total = 0.0f64;
    for chr in chrs {
        let path = make_chr_path(prefix, chr, suffix);
        let content =
            std::fs::read_to_string(&path).with_context(|| format!("reading M file '{}'", path))?;
        let m: f64 = content
            .split_whitespace()
            .next()
            .ok_or_else(|| anyhow::anyhow!("empty M file '{}'", path))?
            .parse()
            .with_context(|| format!("parsing M value from '{}'", path))?;
        total += m;
    }
    Ok(total)
}

/// Read per-annotation M values from per-chromosome `.M[_5_50]` files.
///
/// Each file contains K tab-separated values (K=1 for scalar LD scores).
/// Returns a `Vec<f64>` of length K, summed across chromosomes.
pub fn read_m_vec(prefix: &str, suffix: &str) -> Result<Vec<f64>> {
    let chrs = get_present_chrs(prefix, suffix);
    anyhow::ensure!(
        !chrs.is_empty(),
        "No chromosome M files found for prefix '{}' suffix '{}'",
        prefix,
        suffix
    );

    let mut totals: Vec<f64> = Vec::new();
    for chr in chrs {
        let path = make_chr_path(prefix, chr, suffix);
        let content =
            std::fs::read_to_string(&path).with_context(|| format!("reading M file '{}'", path))?;
        let vals: Vec<f64> = content
            .split_whitespace()
            .map(|s| {
                s.parse::<f64>()
                    .with_context(|| format!("parsing M value '{}' in '{}'", s, path))
            })
            .collect::<Result<_>>()?;
        anyhow::ensure!(!vals.is_empty(), "empty M file '{}'", path);
        if totals.is_empty() {
            totals = vals;
        } else {
            anyhow::ensure!(
                totals.len() == vals.len(),
                "M file '{}' has {} values but first chromosome had {} — annotation count mismatch",
                path,
                vals.len(),
                totals.len()
            );
            for (t, v) in totals.iter_mut().zip(vals.iter()) {
                *t += v;
            }
        }
    }
    Ok(totals)
}

/// Concatenate per-chromosome LazyFrames, accepting .gz, .bz2, or plain files.
pub fn concat_chrs_any(prefix: &str, suffixes: &[&str]) -> Result<LazyFrame> {
    for suffix in suffixes {
        let chrs = get_present_chrs(prefix, suffix);
        if chrs.is_empty() {
            continue;
        }

        let frames: Vec<LazyFrame> = chrs
            .iter()
            .map(|chr| {
                let path = make_chr_path(prefix, *chr, suffix);
                let resolved = maybe_decompress_bz2(&path)?;
                let resolved_str = resolved.to_string_lossy();
                LazyCsvReader::new(resolved_str.as_ref().into())
                    .with_separator(b'\t')
                    .with_has_header(true)
                    .finish()
                    .map_err(anyhow::Error::from)
            })
            .collect::<Result<_>>()?;

        return Ok(concat(frames, UnionArgs::default())?);
    }

    anyhow::bail!("No chromosome files found for prefix '{}' (tried {:?})", prefix, suffixes);
}

// ---------------------------------------------------------------------------
// Annotation file reader
// ---------------------------------------------------------------------------

/// Read a partitioned LD score annotation file into a dense `Array2<f64>`.
///
/// Full format (`thin=false`): `CHR SNP BP CM ANNOT1 …` — skips first 4 columns.
/// Thin format (`thin=true`): all columns are annotations.
/// Tries `{prefix}.annot.gz`, `{prefix}.annot.bz2`, then `{prefix}.annot`.
pub fn read_annot(prefix: &str, thin: bool) -> Result<(Array2<f64>, Vec<String>)> {
    let path = {
        let gz = format!("{}.annot.gz", prefix);
        let bz2 = format!("{}.annot.bz2", prefix);
        let plain = format!("{}.annot", prefix);
        if std::path::Path::new(&gz).exists() {
            gz
        } else if std::path::Path::new(&bz2).exists() {
            bz2
        } else if std::path::Path::new(&plain).exists() {
            plain
        } else {
            anyhow::bail!("Annotation file not found: '{}.annot[.gz|.bz2]'", prefix);
        }
    };

    let resolved = maybe_decompress_bz2(&path)?;
    let resolved_str = resolved.to_string_lossy();
    let df = LazyCsvReader::new(resolved_str.as_ref().into())
        .with_separator(b'\t')
        .with_has_header(true)
        .finish()
        .with_context(|| format!("scanning annot file '{}'", path))?
        .collect()
        .with_context(|| format!("reading annot file '{}'", path))?;

    let all_cols = df.get_column_names();

    // Annotation columns start at index 4 for full format, 0 for thin.
    let skip = if thin { 0 } else { 4 };
    anyhow::ensure!(
        all_cols.len() > skip,
        "Annot file '{}' has {} columns; expected > {} ({})",
        path,
        all_cols.len(),
        skip,
        if thin {
            "thin format"
        } else {
            "full format: CHR SNP BP CM + annotations"
        }
    );

    let col_names: Vec<String> = all_cols[skip..].iter().map(|s| s.to_string()).collect();
    let n_annot = col_names.len();
    let n_rows = df.height();

    let mut matrix = Array2::<f64>::zeros((n_rows, n_annot));
    for (j, name) in col_names.iter().enumerate() {
        let s = df
            .column(name)
            .with_context(|| format!("column '{}' in annot file '{}'", name, path))?
            .cast(&DataType::Float64)
            .with_context(|| format!("casting annot column '{}' to f64", name))?;
        let ca = s
            .f64()
            .with_context(|| format!("annot column '{}' as f64", name))?;
        for (i, val) in ca.into_iter().enumerate() {
            matrix[[i, j]] = val.unwrap_or(0.0);
        }
    }

    Ok((matrix, col_names))
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    #[allow(unused_imports)]
    use std::io::Write;

    /// read_m_total should sum values across two single-value chromosome M files.
    /// Mirrors Python Test_M::test_M_loop (sums per-chr M values).
    #[test]
    fn test_read_m_total_two_chrs() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = format!("{}/test", dir.path().to_str().unwrap());

        std::fs::write(format!("{}1.l2.M", prefix), "1000\n").unwrap();
        std::fs::write(format!("{}2.l2.M", prefix), "2000\n").unwrap();

        let total = read_m_total(&prefix, ".l2.M").unwrap();
        assert!((total - 3000.0).abs() < 0.01, "total={}", total);
    }

    /// read_m_total should error when no chromosome M files exist.
    /// Mirrors Python Test_M::test_bad_M.
    #[test]
    fn test_read_m_total_missing_files_errors() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = format!("{}/nofiles", dir.path().to_str().unwrap());
        assert!(read_m_total(&prefix, ".l2.M").is_err());
    }

    /// read_m_vec on single-column M files returns Vec of length 1.
    #[test]
    fn test_read_m_vec_single_col() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = format!("{}/test", dir.path().to_str().unwrap());

        std::fs::write(format!("{}1.l2.M", prefix), "1000\n").unwrap();
        std::fs::write(format!("{}2.l2.M", prefix), "2000\n").unwrap();

        let v = read_m_vec(&prefix, ".l2.M").unwrap();
        assert_eq!(v.len(), 1);
        assert!((v[0] - 3000.0).abs() < 0.01, "v[0]={}", v[0]);
    }

    /// read_m_vec on two-column (partitioned) M files returns Vec of length 2
    /// with correct per-annotation sums.
    #[test]
    fn test_read_m_vec_multi_col() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = format!("{}/test", dir.path().to_str().unwrap());

        // chr1: annot1=100, annot2=200
        std::fs::write(format!("{}1.l2.M", prefix), "100\t200\n").unwrap();
        // chr2: annot1=300, annot2=400
        std::fs::write(format!("{}2.l2.M", prefix), "300\t400\n").unwrap();

        let v = read_m_vec(&prefix, ".l2.M").unwrap();
        assert_eq!(v.len(), 2);
        assert!((v[0] - 400.0).abs() < 0.01, "annot1 sum={}", v[0]); // 100+300
        assert!((v[1] - 600.0).abs() < 0.01, "annot2 sum={}", v[1]); // 200+400
    }

    /// scan_ldscore can read the existing test .ldscore.gz file and returns
    /// a LazyFrame with the expected columns.
    /// Mirrors Python Test_ldscore::test_ldscore.
    #[test]
    fn test_scan_ldscore_gz() {
        let manifest = env!("CARGO_MANIFEST_DIR");
        let path = format!("{}/tests/fixtures/test.l2.ldscore.gz", manifest);
        let lf = scan_ldscore(&path).unwrap();
        let df = lf.collect().unwrap();

        let cols: Vec<&str> = df.get_column_names().iter().map(|s| s.as_str()).collect();
        assert!(cols.contains(&"SNP"), "missing SNP column");
        assert!(cols.contains(&"AL2"), "missing AL2 column");
        assert!(cols.contains(&"BL2"), "missing BL2 column");
        assert_eq!(df.height(), 22, "expected 22 SNPs, got {}", df.height());
    }
}
