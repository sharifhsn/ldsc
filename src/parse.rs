/// File parsing utilities: Polars LazyFrame readers for `.sumstats` and `.ldscore` files.
use anyhow::{Context, Result};
use ndarray::{Array2, s};
use polars::prelude::*;
use std::fs::File;
use std::path::Path;
use std::path::PathBuf;
use std::sync::{Mutex, OnceLock};

use bzip2_rs::DecoderReader;
use flate2::read::GzDecoder;
use tempfile::Builder as TempBuilder;

static BZ2_TEMPFILES: OnceLock<Mutex<Vec<tempfile::TempPath>>> = OnceLock::new();
static WS_TEMPFILES: OnceLock<Mutex<Vec<tempfile::TempPath>>> = OnceLock::new();

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

fn store_temp_path(path: tempfile::TempPath) {
    WS_TEMPFILES
        .get_or_init(|| Mutex::new(Vec::new()))
        .lock()
        .expect("whitespace tempfiles mutex poisoned")
        .push(path);
}

fn first_line_contains_tab(path: &str, gz: bool) -> Result<bool> {
    use std::io::{BufRead, BufReader};
    let file = File::open(path).with_context(|| format!("opening '{}'", path))?;
    let reader: Box<dyn BufRead> = if gz {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    for line in reader.lines() {
        let line = line?;
        if !line.is_empty() {
            return Ok(line.contains('\t'));
        }
    }
    Ok(false)
}

fn normalize_whitespace_to_tsv(path: &str, gz: bool) -> Result<PathBuf> {
    use std::io::{BufRead, BufReader, BufWriter, Write};
    if first_line_contains_tab(path, gz)? {
        return Ok(PathBuf::from(path));
    }

    let file = File::open(path).with_context(|| format!("opening '{}'", path))?;
    let reader: Box<dyn BufRead> = if gz {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut tmp = TempBuilder::new()
        .prefix("ldsc_ws_")
        .suffix(".tmp")
        .tempfile()
        .context("creating temp file for whitespace normalization")?;
    {
        let mut w = BufWriter::new(&mut tmp);
        for line in reader.lines() {
            let line = line?;
            if line.is_empty() {
                writeln!(w)?;
                continue;
            }
            let cols: Vec<&str> = line.split_whitespace().collect();
            if cols.is_empty() {
                writeln!(w)?;
            } else {
                writeln!(w, "{}", cols.join("\t"))?;
            }
        }
    }

    let temp_path = tmp.into_temp_path();
    let temp_buf = temp_path.to_path_buf();
    store_temp_path(temp_path);
    Ok(temp_buf)
}

pub(crate) fn resolve_text_path(path: &str) -> Result<PathBuf> {
    if path.ends_with(".bz2") {
        let tmp = maybe_decompress_bz2(path)?;
        return normalize_whitespace_to_tsv(tmp.to_string_lossy().as_ref(), false);
    }
    if path.ends_with(".gz") {
        return normalize_whitespace_to_tsv(path, true);
    }
    normalize_whitespace_to_tsv(path, false)
}

// ---------------------------------------------------------------------------
// Readers
// ---------------------------------------------------------------------------

/// Scan a `.sumstats[.gz|.bz2]` file as a Polars LazyFrame.
pub fn scan_sumstats(path: &str) -> Result<LazyFrame> {
    let resolved = resolve_text_path(path)?;
    let resolved_str = resolved.to_string_lossy();
    let lf = LazyCsvReader::new(resolved_str.as_ref().into())
        .with_separator(b'\t')
        .with_has_header(true)
        .with_null_values(Some(NullValues::AllColumns(vec!["NA".into(), ".".into()])))
        .finish()?;
    Ok(lf)
}

/// Scan a `.ldscore[.gz|.bz2]` file as a Polars LazyFrame.
pub fn scan_ldscore(path: &str) -> Result<LazyFrame> {
    let resolved = resolve_text_path(path)?;
    let resolved_str = resolved.to_string_lossy();
    let lf = LazyCsvReader::new(resolved_str.as_ref().into())
        .with_separator(b'\t')
        .with_has_header(true)
        .with_null_values(Some(NullValues::AllColumns(vec!["NA".into(), ".".into()])))
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

/// Read per-annotation M values from multiple prefixes (comma-separated in Python).
/// Concatenates the per-prefix M vectors in order.
pub fn read_m_vec_list(prefixes: &[String], suffix: &str) -> Result<Vec<f64>> {
    let mut out: Vec<f64> = Vec::new();
    for prefix in prefixes {
        let mut v = read_m_vec(prefix, suffix)?;
        out.append(&mut v);
    }
    Ok(out)
}

/// Sum total M across multiple prefixes.
pub fn read_m_total_list(prefixes: &[String], suffix: &str) -> Result<f64> {
    let mut total = 0.0f64;
    for prefix in prefixes {
        total += read_m_total(prefix, suffix)?;
    }
    Ok(total)
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
                let resolved = resolve_text_path(&path)?;
                let resolved_str = resolved.to_string_lossy();
                LazyCsvReader::new(resolved_str.as_ref().into())
                    .with_separator(b'\t')
                    .with_has_header(true)
                    .with_null_values(Some(NullValues::AllColumns(vec!["NA".into(), ".".into()])))
                    .finish()
                    .map_err(anyhow::Error::from)
            })
            .collect::<Result<_>>()?;

        return Ok(concat(frames, UnionArgs::default())?);
    }

    anyhow::bail!(
        "No chromosome files found for prefix '{}' (tried {:?})",
        prefix,
        suffixes
    );
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

    let resolved = resolve_text_path(&path)?;
    let resolved_str = resolved.to_string_lossy();
    let df = LazyCsvReader::new(resolved_str.as_ref().into())
        .with_separator(b'\t')
        .with_has_header(true)
        .with_null_values(Some(NullValues::AllColumns(vec!["NA".into(), ".".into()])))
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
// Overlap-annot helpers
// ---------------------------------------------------------------------------

fn get_present_chrs_any(prefix: &str, suffixes: &[&str]) -> Vec<u8> {
    use std::collections::BTreeSet;
    let mut set = BTreeSet::new();
    for suffix in suffixes {
        for chr in get_present_chrs(prefix, suffix) {
            set.insert(chr);
        }
    }
    set.into_iter().collect()
}

fn resolve_frq_path(prefix: &str, chr: Option<u8>) -> Result<String> {
    let base = match chr {
        Some(c) => make_chr_path(prefix, c, ""),
        None => prefix.to_string(),
    };
    let candidates = [".frq.gz", ".frq.bz2", ".frq"];
    for suffix in candidates {
        let path = format!("{}{}", base, suffix);
        if Path::new(&path).exists() {
            return Ok(path);
        }
    }
    anyhow::bail!("Frequency file not found: '{}.frq[.gz|.bz2]'", base);
}

/// Read a .frq file and return a MAF filter mask: keep 0.05 < FRQ < 0.95.
/// Accepts either FRQ or MAF column names (MAF is renamed to FRQ in Python).
pub fn read_frq_mask(path: &str) -> Result<Vec<bool>> {
    let resolved = resolve_text_path(path)?;
    let resolved_str = resolved.to_string_lossy();
    let df = LazyCsvReader::new(resolved_str.as_ref().into())
        .with_separator(b'\t')
        .with_has_header(true)
        .with_null_values(Some(NullValues::AllColumns(vec!["NA".into(), ".".into()])))
        .finish()
        .with_context(|| format!("scanning frq file '{}'", path))?
        .collect()
        .with_context(|| format!("reading frq file '{}'", path))?;

    let col_name = if df.get_column_names().iter().any(|c| c.as_str() == "FRQ") {
        "FRQ"
    } else if df.get_column_names().iter().any(|c| c.as_str() == "MAF") {
        "MAF"
    } else {
        anyhow::bail!("Frequency file '{}' has no FRQ or MAF column", path);
    };

    let s = df
        .column(col_name)
        .with_context(|| format!("column '{}' in frq file '{}'", col_name, path))?
        .cast(&DataType::Float64)
        .with_context(|| format!("casting '{}' to f64 in frq file '{}'", col_name, path))?;
    let ca = s
        .f64()
        .with_context(|| format!("frq column '{}' as f64", col_name))?;
    let mask = ca
        .into_iter()
        .map(|opt| opt.map(|v| v > 0.05 && v < 0.95).unwrap_or(false))
        .collect();
    Ok(mask)
}

/// Compute annotation overlap matrix XᵀX and total SNP count (M_tot).
///
/// When `chr_split` is true, reads per-chromosome `{prefix}{chr}.annot[.gz|.bz2]`.
/// When false, reads single `{prefix}.annot[.gz|.bz2]`.
///
/// `frqfile` / `frqfile_chr` apply the Python MAF filter (0.05 < FRQ < 0.95)
/// before computing XᵀX.
pub fn read_overlap_matrix(
    prefixes: &[String],
    frqfile: Option<&str>,
    frqfile_chr: Option<&str>,
    chr_split: bool,
) -> Result<(Array2<f64>, usize, Vec<String>)> {
    anyhow::ensure!(!prefixes.is_empty(), "No annotation prefixes provided");

    let mut raw_names_per_prefix: Vec<Vec<String>> = Vec::new();
    let mut overlap: Option<Array2<f64>> = None;
    let mut m_tot: usize = 0;

    if chr_split {
        let chrs = get_present_chrs_any(&prefixes[0], &[".annot.gz", ".annot.bz2", ".annot"]);
        anyhow::ensure!(
            !chrs.is_empty(),
            "No .annot files found for prefix '{}'",
            prefixes[0]
        );

        for chr in chrs {
            let mut matrices: Vec<Array2<f64>> = Vec::new();
            let mut n_rows: Option<usize> = None;
            for (idx, prefix) in prefixes.iter().enumerate() {
                let chr_prefix = make_chr_path(prefix, chr, "");
                let (mat, names) = read_annot(&chr_prefix, false)?;
                if raw_names_per_prefix.len() <= idx {
                    raw_names_per_prefix.push(names.clone());
                } else {
                    anyhow::ensure!(
                        raw_names_per_prefix[idx] == names,
                        "Annotation columns for '{}' differ across chromosomes",
                        prefix
                    );
                }
                if let Some(nr) = n_rows {
                    anyhow::ensure!(
                        mat.nrows() == nr,
                        "Annotation rows mismatch across prefixes for chr {}",
                        chr
                    );
                } else {
                    n_rows = Some(mat.nrows());
                }
                matrices.push(mat);
            }

            let n_rows = n_rows.unwrap_or(0);
            let total_cols: usize = matrices.iter().map(|m| m.ncols()).sum();
            let mut stacked = Array2::<f64>::zeros((n_rows, total_cols));
            let mut col_offset = 0usize;
            for mat in matrices {
                let cols = mat.ncols();
                stacked
                    .slice_mut(s![.., col_offset..col_offset + cols])
                    .assign(&mat);
                col_offset += cols;
            }

            let (matrix, m_chr) = if let Some(frq_prefix) = frqfile_chr {
                let frq_path = resolve_frq_path(frq_prefix, Some(chr))?;
                let mask = read_frq_mask(&frq_path)?;
                anyhow::ensure!(
                    mask.len() == n_rows,
                    "FRQ file '{}' has {} rows; expected {}",
                    frq_path,
                    mask.len(),
                    n_rows
                );
                let keep = mask.iter().filter(|&&b| b).count();
                let mut filtered = Array2::<f64>::zeros((keep, total_cols));
                let mut row = 0usize;
                for (i, keep_row) in mask.iter().enumerate() {
                    if *keep_row {
                        filtered.row_mut(row).assign(&stacked.row(i));
                        row += 1;
                    }
                }
                (filtered.t().dot(&filtered), keep)
            } else {
                (stacked.t().dot(&stacked), n_rows)
            };

            m_tot += m_chr;
            if let Some(ref mut acc) = overlap {
                *acc += &matrix;
            } else {
                overlap = Some(matrix);
            }
        }
    } else {
        let mut matrices: Vec<Array2<f64>> = Vec::new();
        let mut n_rows: Option<usize> = None;
        for (idx, prefix) in prefixes.iter().enumerate() {
            let (mat, names) = read_annot(prefix, false)?;
            if raw_names_per_prefix.len() <= idx {
                raw_names_per_prefix.push(names.clone());
            } else {
                anyhow::ensure!(
                    raw_names_per_prefix[idx] == names,
                    "Annotation columns for '{}' differ across files",
                    prefix
                );
            }
            if let Some(nr) = n_rows {
                anyhow::ensure!(
                    mat.nrows() == nr,
                    "Annotation rows mismatch across prefixes"
                );
            } else {
                n_rows = Some(mat.nrows());
            }
            matrices.push(mat);
        }

        let n_rows = n_rows.unwrap_or(0);
        let total_cols: usize = matrices.iter().map(|m| m.ncols()).sum();
        let mut stacked = Array2::<f64>::zeros((n_rows, total_cols));
        let mut col_offset = 0usize;
        for mat in matrices {
            let cols = mat.ncols();
            stacked
                .slice_mut(s![.., col_offset..col_offset + cols])
                .assign(&mat);
            col_offset += cols;
        }

        let (matrix, m_single) = if let Some(frq_prefix) = frqfile {
            let frq_path = resolve_frq_path(frq_prefix, None)?;
            let mask = read_frq_mask(&frq_path)?;
            anyhow::ensure!(
                mask.len() == n_rows,
                "FRQ file '{}' has {} rows; expected {}",
                frq_path,
                mask.len(),
                n_rows
            );
            let keep = mask.iter().filter(|&&b| b).count();
            let mut filtered = Array2::<f64>::zeros((keep, total_cols));
            let mut row = 0usize;
            for (i, keep_row) in mask.iter().enumerate() {
                if *keep_row {
                    filtered.row_mut(row).assign(&stacked.row(i));
                    row += 1;
                }
            }
            (filtered.t().dot(&filtered), keep)
        } else {
            (stacked.t().dot(&stacked), n_rows)
        };

        m_tot += m_single;
        overlap = Some(matrix);
    }

    let overlap = overlap.context("overlap matrix not computed")?;
    let col_names: Vec<String> = if prefixes.len() > 1 {
        raw_names_per_prefix
            .iter()
            .enumerate()
            .flat_map(|(i, names)| names.iter().map(move |n| format!("{}_{}", n, i)))
            .collect()
    } else {
        raw_names_per_prefix.first().cloned().unwrap_or_default()
    };

    Ok((overlap, m_tot, col_names))
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

    #[test]
    fn test_overlap_matrix_no_frq() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("test");
        let annot_path = format!("{}.annot", prefix.to_string_lossy());
        std::fs::write(
            &annot_path,
            "CHR\tBP\tSNP\tCM\tC1\tC2\tC3\n\
1\t1\trs_0\t0\t1\t0\t0\n\
1\t2\trs_1\t0\t0\t1\t1\n\
1\t3\trs_2\t0\t0\t1\t1\n",
        )
        .unwrap();

        let (overlap, m_tot, names) =
            read_overlap_matrix(&[prefix.to_string_lossy().to_string()], None, None, false)
                .unwrap();

        assert_eq!(m_tot, 3);
        assert_eq!(
            names,
            vec!["C1".to_string(), "C2".to_string(), "C3".to_string()]
        );
        assert!((overlap[[0, 0]] - 1.0).abs() < 1e-6);
        assert!((overlap[[1, 1]] - 2.0).abs() < 1e-6);
        assert!((overlap[[2, 2]] - 2.0).abs() < 1e-6);
        assert!((overlap[[1, 2]] - 2.0).abs() < 1e-6);
        assert!((overlap[[2, 1]] - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_overlap_matrix_with_frq() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("test");
        let annot_path = format!("{}.annot", prefix.to_string_lossy());
        let frq_path = format!("{}.frq", prefix.to_string_lossy());

        std::fs::write(
            &annot_path,
            "CHR\tBP\tSNP\tCM\tC1\tC2\tC3\n\
1\t1\trs_0\t0\t1\t0\t0\n\
1\t2\trs_1\t0\t0\t1\t1\n\
1\t3\trs_2\t0\t0\t1\t1\n",
        )
        .unwrap();
        std::fs::write(
            &frq_path,
            "CHR\tSNP\tCM\tBP\tA1\tFRQ\n\
1\trs_0\t0\t1\tA\t0.7\n\
1\trs_1\t0\t2\tA\t0.1\n\
1\trs_2\t0\t3\tA\t0.01\n",
        )
        .unwrap();

        let (overlap, m_tot, _) = read_overlap_matrix(
            &[prefix.to_string_lossy().to_string()],
            Some(prefix.to_string_lossy().as_ref()),
            None,
            false,
        )
        .unwrap();

        assert_eq!(m_tot, 2);
        assert!((overlap[[0, 0]] - 1.0).abs() < 1e-6);
        assert!((overlap[[1, 1]] - 1.0).abs() < 1e-6);
        assert!((overlap[[2, 2]] - 1.0).abs() < 1e-6);
        assert!((overlap[[1, 2]] - 1.0).abs() < 1e-6);
        assert!((overlap[[2, 1]] - 1.0).abs() < 1e-6);
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

    #[test]
    fn test_scan_sumstats_whitespace() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sumstats.txt");
        let content = "\
SNP A1 A2 Z N
rs1 A G 1.0 100
rs2 C T 2.5 200
";
        std::fs::write(&path, content).unwrap();
        let lf = scan_sumstats(path.to_str().unwrap()).unwrap();
        let df = lf.collect().unwrap();
        assert_eq!(df.height(), 2);
        let cols: Vec<String> = df
            .get_column_names()
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert!(cols.contains(&"SNP".to_string()));
        assert!(cols.contains(&"Z".to_string()));
    }

    #[test]
    fn test_scan_ldscore_whitespace() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("ldscore.txt");
        let content = "\
CHR SNP BP L2
1 rs1 100 1.0
1 rs2 200 2.0
";
        std::fs::write(&path, content).unwrap();
        let lf = scan_ldscore(path.to_str().unwrap()).unwrap();
        let df = lf.collect().unwrap();
        assert_eq!(df.height(), 2);
        let cols: Vec<String> = df
            .get_column_names()
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert!(cols.contains(&"SNP".to_string()));
        assert!(cols.contains(&"L2".to_string()));
    }
}
