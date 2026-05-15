use crate::frame::{self, Frame};
use crate::la::{MatF, mat_add_in_place, matmul_tn_to, par_default};
/// File parsing utilities for `.sumstats`, `.ldscore`, `.annot`, `.frq`, `.bim`,
/// `.l2.M[_5_50]`, etc. Pure-Rust readers (csv-style) — no polars.
use anyhow::{Context, Result};
use faer::Accum;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Read a TSV/whitespace-separated text file (`.sumstats`, `.ldscore`, etc.)
/// into a `Frame`. Handles `.gz` and `.bz2` decompression transparently.
pub fn scan_tsv(path: &str) -> Result<Frame> {
    frame::read_tsv(path)
}

/// Alias for backward compatibility.
pub fn scan_sumstats(path: &str) -> Result<Frame> {
    scan_tsv(path)
}

/// Alias for backward compatibility.
pub fn scan_ldscore(path: &str) -> Result<Frame> {
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

/// Concatenate per-chromosome frames, accepting .gz, .bz2, or plain files.
/// Tries each suffix in order; uses the first one that has any chromosome
/// files present.
pub fn concat_chrs_any(prefix: &str, suffixes: &[&str]) -> Result<Frame> {
    for suffix in suffixes {
        let chrs = get_present_chrs(prefix, suffix);
        if chrs.is_empty() {
            continue;
        }

        let frames: Vec<Frame> = chrs
            .iter()
            .map(|chr| {
                let path = make_chr_path(prefix, *chr, suffix);
                scan_tsv(&path)
            })
            .collect::<Result<_>>()?;

        return frame::concat_rows(frames);
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
    let df = scan_tsv(path).with_context(|| format!("scanning annot file '{}'", path))?;

    let all_cols = df.column_names_owned();

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

    let col_names: Vec<String> = all_cols[skip..].to_vec();
    let n_annot = col_names.len();
    let n_rows = df.height();

    let mut matrix = MatF::zeros(n_rows, n_annot);
    for (j, name) in col_names.iter().enumerate() {
        let casted = df
            .column(name)
            .with_context(|| format!("column '{}' in annot file '{}'", name, path))?
            .cast_to_f64()
            .with_context(|| format!("casting annot column '{}' to f64", name))?;
        let v = casted.as_f64()?;
        for (i, &val) in v.iter().enumerate() {
            matrix[(i, j)] = if val.is_nan() { 0.0 } else { val };
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
    let df = scan_tsv(path).with_context(|| format!("scanning frq file '{}'", path))?;

    let col_name = if df.has_column("FRQ") {
        "FRQ"
    } else if df.has_column("MAF") {
        "MAF"
    } else {
        anyhow::bail!("Frequency file '{}' has no FRQ or MAF column", path);
    };

    let casted = df
        .column(col_name)
        .with_context(|| format!("column '{}' in frq file '{}'", col_name, path))?
        .cast_to_f64()
        .with_context(|| format!("casting '{}' to f64 in frq file '{}'", col_name, path))?;
    let mask: Vec<bool> = casted
        .as_f64()?
        .iter()
        .map(|v| !v.is_nan() && *v > 0.05 && *v < 0.95)
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
                    par_default(),
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
                    par_default(),
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
                par_default(),
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
                par_default(),
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
    parse_bim_reader(BufReader::new(f))
}

/// Parse a PLINK .bim file from an in-memory string. Browser-side
/// counterpart to [`parse_bim`].
pub fn parse_bim_str(s: &str) -> Result<Vec<BimRecord>> {
    parse_bim_reader(BufReader::new(s.as_bytes()))
}

/// Parse a PLINK .bim file streamingly from any `BufRead` source.
/// The browser frontend wraps a `Blob`-backed `Read` shim so the BIM
/// never has to live in JS-string heap (which caps at ~256-512 MB
/// depending on browser → ~4-8M SNPs of BIM). Once parsed, the
/// `Vec<BimRecord>` is in wasm linear memory, which has the wider
/// 4 GB cap.
///
/// Per-record memory is ~70-80 bytes once you count the `String`
/// header for the rsID + the heap allocation for typical
/// `rs1234567`-style IDs. So the practical wasm32 cap is ~50M SNPs
/// for one buffered parse — bigger inputs need
/// [`parse_bim_reader_filtered`] to skip rows the caller doesn't
/// keep.
pub fn parse_bim_reader<R: BufRead>(reader: R) -> Result<Vec<BimRecord>> {
    parse_bim_reader_filtered(reader, |_chr| true)
}

/// Streaming BIM parser with an inline `chr` predicate filter. Rows
/// where `keep_chr(chr)` returns `false` are decoded enough to advance
/// `bed_idx` (the per-line counter) but not stored. Browser-side per-
/// chromosome workers use this to consume just their assigned slice of
/// the genome without holding the full `Vec<BimRecord>` for the rest
/// — without it, 4 workers × ~100 MB BIM ≈ 400 MB of redundant parsed
/// state, and at the documented 60M-SNP scale the per-worker fan-out
/// blows past the 4 GB wasm32 cap.
///
/// `bed_idx` is preserved for *every* line (matching the unfiltered
/// path), so the kept records still index correctly into the original
/// BED file's SNP slabs.
pub fn parse_bim_reader_filtered<R: BufRead, F: FnMut(u8) -> bool>(
    reader: R,
    mut keep_chr: F,
) -> Result<Vec<BimRecord>> {
    let mut records = Vec::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("reading BIM line {}", line_no + 1))?;
        let mut cols = line.split('\t');
        let chr_str = cols
            .next()
            .ok_or_else(|| anyhow::anyhow!("BIM line {}: empty", line_no + 1))?;
        let chr = chr_str.parse::<u8>().unwrap_or(0);
        if !keep_chr(chr) {
            continue;
        }
        let snp = cols
            .next()
            .ok_or_else(|| anyhow::anyhow!("BIM line {}: expected SNP column", line_no + 1))?
            .to_string();
        let cm = cols
            .next()
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(0.0);
        let bp = cols.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
        // Columns 5+6 (A1/A2) are still validated by the original
        // 6-col contract — we read them via remaining iterator advance
        // to mirror the strictness of `parse_bim_reader`'s old check.
        let _a1 = cols
            .next()
            .ok_or_else(|| anyhow::anyhow!("BIM line {}: missing A1 col", line_no + 1))?;
        let _a2 = cols
            .next()
            .ok_or_else(|| anyhow::anyhow!("BIM line {}: missing A2 col", line_no + 1))?;
        records.push(BimRecord {
            chr,
            snp,
            cm,
            bp,
            bed_idx: line_no,
        });
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `parse_bim_str` (in-memory) must produce the same records as
    /// `parse_bim` (file). This pins down the contract that the WASM
    /// frontend can swap in `parse_bim_str` without changing semantics.
    #[test]
    fn parse_bim_str_matches_file() {
        let bim = "1\trs1\t0\t100\tA\tG\n\
                   1\trs2\t0.1\t200\tC\tT\n\
                   2\trs3\t0\t300\tG\tA\n";
        let recs = parse_bim_str(bim).expect("parse_bim_str");
        assert_eq!(recs.len(), 3);
        assert_eq!(recs[0].snp, "rs1");
        assert_eq!(recs[0].chr, 1);
        assert_eq!(recs[0].bp, 100);
        assert_eq!(recs[1].cm, 0.1);
        assert_eq!(recs[2].chr, 2);
        // bed_idx is the 0-based row position
        assert_eq!(recs[2].bed_idx, 2);
    }
}
