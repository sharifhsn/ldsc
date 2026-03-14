use crate::la::MatF;
use crate::parse::{self, BimRecord};
use anyhow::{Context, Result};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Write as IoWrite};

/// Count the number of individuals in a PLINK .fam file (one row per individual).
pub fn count_fam(path: &str) -> Result<usize> {
    let f = File::open(path).with_context(|| format!("opening FAM file '{}'", path))?;
    Ok(BufReader::new(f).lines().count())
}

/// Parse (FID, IID) pairs from a PLINK .fam file in order.
pub fn parse_fam(path: &str) -> Result<Vec<(String, String)>> {
    let f = File::open(path).with_context(|| format!("opening FAM file '{}'", path))?;
    let reader = BufReader::new(f);
    let mut ids = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("reading FAM line {}", i + 1))?;
        let cols: Vec<&str> = line.split_whitespace().collect();
        anyhow::ensure!(
            cols.len() >= 2,
            "FAM line {}: expected ≥ 2 columns (FID IID ...), got {}",
            i + 1,
            cols.len()
        );
        ids.push((cols[0].to_string(), cols[1].to_string()));
    }
    Ok(ids)
}

/// Load the 0-based FAM indices for individuals listed in a keep file (FID IID per line).
pub(super) fn load_individual_indices(
    keep_path: &str,
    fam_ids: &[(String, String)],
) -> Result<Vec<isize>> {
    let file =
        File::open(keep_path).with_context(|| format!("opening keep file '{}'", keep_path))?;
    let keep_set: HashSet<(String, String)> = BufReader::new(file)
        .lines()
        .map_while(Result::ok)
        .filter(|l| !l.trim_start().starts_with('#'))
        .filter_map(|l| {
            let mut cols = l.split_whitespace();
            let fid = cols.next()?.to_string();
            let iid = cols.next()?.to_string();
            Some((fid, iid))
        })
        .collect();

    let indices: Vec<isize> = fam_ids
        .iter()
        .enumerate()
        .filter(|(_, (fid, iid))| keep_set.contains(&(fid.clone(), iid.clone())))
        .map(|(i, _)| i as isize)
        .collect();

    anyhow::ensure!(
        !indices.is_empty(),
        "--keep '{}': none of {} individuals matched the FAM ({} individuals)",
        keep_path,
        keep_set.len(),
        fam_ids.len()
    );
    println!(
        "  --keep: retaining {}/{} individuals",
        indices.len(),
        fam_ids.len()
    );
    Ok(indices)
}

/// Write per-SNP LD scores to a gzip TSV.  `col_names` is `["L2"]` or `["{ANNOT}L2", …]`.
pub(super) fn write_ldscore_refs(
    path: &str,
    snps: &[&BimRecord],
    l2: &MatF,
    col_names: &[String],
) -> Result<()> {
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::BufWriter;

    let file = File::create(path).with_context(|| format!("creating output '{}'", path))?;
    let gz = GzEncoder::new(file, Compression::fast());
    let mut writer = BufWriter::new(gz);

    write!(writer, "CHR\tSNP\tBP")?;
    for name in col_names {
        write!(writer, "\t{}", name)?;
    }
    writeln!(writer)?;

    for (i, snp) in snps.iter().enumerate() {
        write!(writer, "{}\t{}\t{}", snp.chr, snp.snp, snp.bp)?;
        for k in 0..col_names.len() {
            write!(writer, "\t{:.3}", l2[(i, k)])?;
        }
        writeln!(writer)?;
    }

    let gz = writer.into_inner().context("flushing gzip buffer")?;
    gz.finish().context("finalising gzip output")?;
    Ok(())
}

/// Write an annotation matrix to a gzip TSV (CHR, SNP, BP, CM, ...).
pub(super) fn write_annot_matrix(
    path: &str,
    snps: &[BimRecord],
    annot: &MatF,
    col_names: &[String],
) -> Result<()> {
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::BufWriter;

    let file = File::create(path).with_context(|| format!("creating output '{}'", path))?;
    let gz = GzEncoder::new(file, Compression::fast());
    let mut writer = BufWriter::new(gz);

    write!(writer, "CHR\tSNP\tBP\tCM")?;
    for name in col_names {
        write!(writer, "\t{}", name)?;
    }
    writeln!(writer)?;

    for (i, snp) in snps.iter().enumerate() {
        write!(writer, "{}\t{}\t{}\t{}", snp.chr, snp.snp, snp.bp, snp.cm)?;
        for k in 0..col_names.len() {
            write!(writer, "\t{}", annot[(i, k)])?;
        }
        writeln!(writer)?;
    }

    let gz = writer.into_inner().context("flushing gzip buffer")?;
    gz.finish().context("finalising gzip output")?;
    Ok(())
}

pub(super) fn format_m_vals(vals: &[f64]) -> String {
    let joined = vals
        .iter()
        .map(|v| format!("{v}"))
        .collect::<Vec<_>>()
        .join("\t");
    format!("{joined}\n")
}

/// Load SNP IDs from a file (one per line; first token used; `#` lines skipped).
pub(super) fn load_snp_set(path: &str) -> Result<HashSet<String>> {
    let file = File::open(path).with_context(|| format!("opening SNP list '{}'", path))?;
    let reader = BufReader::new(file);
    let set: HashSet<String> = reader
        .lines()
        .map_while(Result::ok)
        .filter(|l| !l.trim_start().starts_with('#'))
        .filter_map(|l| {
            l.split_whitespace()
                .next()
                .filter(|s| !s.is_empty())
                .map(|s| s.to_string())
        })
        .collect();
    println!("  Loaded {} SNP IDs from '{}'", set.len(), path);
    Ok(set)
}

/// Load SNP IDs from a file with exactly one column (Python --print-snps behavior).
pub(super) fn load_print_snps(path: &str) -> Result<HashSet<String>> {
    let resolved = parse::resolve_text_path(path)?;
    let file = File::open(&resolved).with_context(|| format!("opening SNP list '{}'", path))?;
    let reader = BufReader::new(file);
    let mut set: HashSet<String> = HashSet::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("reading SNP line {}", i + 1))?;
        if line.trim_start().starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.is_empty() {
            continue;
        }
        if cols.len() != 1 {
            anyhow::bail!("--print-snps must refer to a file with a one column of SNP IDs.");
        }
        set.insert(cols[0].to_string());
    }
    println!("  Loaded {} SNP IDs from '{}'", set.len(), path);
    Ok(set)
}
