use crate::la::{MatF, mat_add_in_place, matmul_tn_to};
/// File parsing utilities: Polars LazyFrame readers for `.sumstats` and `.ldscore` files.
use anyhow::{Context, Result};
use faer::{Accum, Par};
use polars::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
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

/// Scan a TSV text file (`.sumstats`, `.ldscore`, etc.) as a Polars LazyFrame.
/// Handles `.gz` and `.bz2` decompression transparently.
pub fn scan_tsv(path: &str) -> Result<LazyFrame> {
    let resolved = resolve_text_path(path)?;
    let resolved_str = resolved.to_string_lossy();
    let lf = LazyCsvReader::new(resolved_str.as_ref().into())
        .with_separator(b'\t')
        .with_has_header(true)
        .with_null_values(Some(NullValues::AllColumns(vec!["NA".into(), ".".into()])))
        .finish()?;
    Ok(lf)
}

/// Alias for backward compatibility.
pub fn scan_sumstats(path: &str) -> Result<LazyFrame> {
    scan_tsv(path)
}

/// Alias for backward compatibility.
pub fn scan_ldscore(path: &str) -> Result<LazyFrame> {
    scan_tsv(path)
}

/// Comma-separated path list → Vec of trimmed, non-empty strings.
pub fn split_paths(raw: &str) -> Vec<String> {
    raw.split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect()
}

/// Build a per-chromosome path; `@` in prefix is replaced with chr number.
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

/// Read per-annotation M values from per-chromosome `.M[_5_50]` files, summed.
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
                scan_tsv(&path)
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

/// Read annotation file into dense `Mat<f64>`. Skips first 4 cols unless `thin`.
pub fn read_annot(prefix: &str, thin: bool) -> Result<(MatF, Vec<String>)> {
    let path = resolve_annot_path(prefix)?;
    read_annot_path(&path, thin)
}

/// Resolve a single annotation file path.
/// Accepts `prefix` with or without `.annot[.gz|.bz2]` suffix.
pub fn resolve_annot_path(prefix: &str) -> Result<String> {
    if prefix.ends_with(".annot") || prefix.ends_with(".annot.gz") || prefix.ends_with(".annot.bz2")
    {
        if std::path::Path::new(prefix).exists() {
            return Ok(prefix.to_string());
        }
        anyhow::bail!("Annotation file not found: '{}'", prefix);
    }

    let gz = format!("{}.annot.gz", prefix);
    let bz2 = format!("{}.annot.bz2", prefix);
    let plain = format!("{}.annot", prefix);
    if std::path::Path::new(&gz).exists() {
        Ok(gz)
    } else if std::path::Path::new(&bz2).exists() {
        Ok(bz2)
    } else if std::path::Path::new(&plain).exists() {
        Ok(plain)
    } else {
        anyhow::bail!("Annotation file not found: '{}.annot[.gz|.bz2]'", prefix);
    }
}

/// Read a partitioned LD score annotation file into a dense `Mat<f64>`.
/// `path` must include the extension (e.g., `.annot.gz`).
pub fn read_annot_path(path: &str, thin: bool) -> Result<(MatF, Vec<String>)> {
    let df = scan_tsv(path)
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

    let mut matrix = MatF::zeros(n_rows, n_annot);
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
            matrix[(i, j)] = val.unwrap_or(0.0);
        }
    }

    Ok((matrix, col_names))
}

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
    let df = scan_tsv(path)
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
pub fn read_overlap_matrix(
    prefixes: &[String],
    frqfile: Option<&str>,
    frqfile_chr: Option<&str>,
    chr_split: bool,
) -> Result<(MatF, usize, Vec<String>)> {
    anyhow::ensure!(!prefixes.is_empty(), "No annotation prefixes provided");

    let mut raw_names_per_prefix: Vec<Vec<String>> = Vec::new();
    let mut overlap: Option<MatF> = None;
    let mut m_tot: usize = 0;

    if chr_split {
        let chrs = get_present_chrs_any(&prefixes[0], &[".annot.gz", ".annot.bz2", ".annot"]);
        anyhow::ensure!(
            !chrs.is_empty(),
            "No .annot files found for prefix '{}'",
            prefixes[0]
        );

        for chr in chrs {
            let mut matrices: Vec<MatF> = Vec::new();
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
            let mut stacked = MatF::zeros(n_rows, total_cols);
            let mut col_offset = 0usize;
            for mat in matrices {
                let cols = mat.ncols();
                for i in 0..n_rows {
                    for j in 0..cols {
                        stacked[(i, col_offset + j)] = mat[(i, j)];
                    }
                }
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
                let mut filtered = MatF::zeros(keep, total_cols);
                let mut row = 0usize;
                for (i, keep_row) in mask.iter().enumerate() {
                    if *keep_row {
                        for j in 0..total_cols {
                            filtered[(row, j)] = stacked[(i, j)];
                        }
                        row += 1;
                    }
                }
                let mut matrix = MatF::zeros(total_cols, total_cols);
                matmul_tn_to(
                    matrix.as_mut(),
                    filtered.as_ref(),
                    filtered.as_ref(),
                    1.0,
                    Accum::Replace,
                    Par::rayon(0),
                );
                (matrix, keep)
            } else {
                let mut matrix = MatF::zeros(total_cols, total_cols);
                matmul_tn_to(
                    matrix.as_mut(),
                    stacked.as_ref(),
                    stacked.as_ref(),
                    1.0,
                    Accum::Replace,
                    Par::rayon(0),
                );
                (matrix, n_rows)
            };

            m_tot += m_chr;
            if let Some(ref mut acc) = overlap {
                mat_add_in_place(acc.as_mut(), matrix.as_ref());
            } else {
                overlap = Some(matrix);
            }
        }
    } else {
        let mut matrices: Vec<MatF> = Vec::new();
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
        let mut stacked = MatF::zeros(n_rows, total_cols);
        let mut col_offset = 0usize;
        for mat in matrices {
            let cols = mat.ncols();
            for i in 0..n_rows {
                for j in 0..cols {
                    stacked[(i, col_offset + j)] = mat[(i, j)];
                }
            }
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
            let mut filtered = MatF::zeros(keep, total_cols);
            let mut row = 0usize;
            for (i, keep_row) in mask.iter().enumerate() {
                if *keep_row {
                    for j in 0..total_cols {
                        filtered[(row, j)] = stacked[(i, j)];
                    }
                    row += 1;
                }
            }
            let mut matrix = MatF::zeros(total_cols, total_cols);
            matmul_tn_to(
                matrix.as_mut(),
                filtered.as_ref(),
                filtered.as_ref(),
                1.0,
                Accum::Replace,
                Par::rayon(0),
            );
            (matrix, keep)
        } else {
            let mut matrix = MatF::zeros(total_cols, total_cols);
            matmul_tn_to(
                matrix.as_mut(),
                stacked.as_ref(),
                stacked.as_ref(),
                1.0,
                Accum::Replace,
                Par::rayon(0),
            );
            (matrix, n_rows)
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

/// Per-SNP metadata from a PLINK .bim file.
/// Columns: CHR  SNP  CM  BP  A1  A2  (tab-separated, no header).
#[derive(Debug, Clone)]
pub struct BimRecord {
    pub chr: u8,
    pub snp: String,
    pub cm: f64,
    pub bp: u32,
    /// 0-based index within the .bed file (position in the full BIM list).
    pub bed_idx: usize,
}

/// Parse a PLINK .bim file.  Returns one `BimRecord` per row.
pub fn parse_bim(path: &str) -> Result<Vec<BimRecord>> {
    let f = File::open(path).with_context(|| format!("opening BIM file '{}'", path))?;
    let reader = BufReader::new(f);
    let mut records = Vec::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("reading BIM line {}", line_no + 1))?;
        let cols: Vec<&str> = line.split('\t').collect();
        anyhow::ensure!(
            cols.len() >= 6,
            "BIM line {}: expected 6 columns, got {}",
            line_no + 1,
            cols.len()
        );

        records.push(BimRecord {
            chr: cols[0].parse::<u8>().unwrap_or(0),
            snp: cols[1].to_string(),
            cm: cols[2].parse::<f64>().unwrap_or(0.0),
            bp: cols[3].parse::<u32>().unwrap_or(0),
            bed_idx: line_no,
        });
    }

    Ok(records)
}
