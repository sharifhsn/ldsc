/// Annotation file generator — replaces make_annot.py.
///
/// For each SNP in a PLINK .bim file, emits a 0/1 annotation indicating
/// whether the SNP's base-pair position falls within any of a set of genomic
/// regions supplied as a UCSC BED file or derived from a gene set.
///
/// Output format (tab-separated, optionally gzip-compressed):
///   CHR  BP  SNP  CM  <annotation_column>
///
/// This matches the "full annot" format consumed by `ldscore --annot`.
use anyhow::{Context, Result};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use crate::cli::MakeAnnotArgs;
use crate::ldscore::{BimRecord, parse_bim};

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

pub fn run(args: MakeAnnotArgs) -> Result<()> {
    let snps =
        parse_bim(&args.bimfile).with_context(|| format!("reading BIM '{}'", args.bimfile))?;
    println!("Read {} SNPs from '{}'", snps.len(), args.bimfile);

    // Build the 0/1 annotation vector.
    let annotation: Vec<u8> = if let Some(ref bed_path) = args.bed_file {
        annotate_from_bed(&snps, bed_path, args.windowsize, args.nomerge)?
    } else if args.gene_set_file.is_some() && args.gene_coord_file.is_some() {
        annotate_from_gene_set(
            &snps,
            args.gene_set_file.as_deref().unwrap(),
            args.gene_coord_file.as_deref().unwrap(),
            args.windowsize,
        )?
    } else {
        anyhow::bail!(
            "Must provide either --bed-file or both --gene-set-file and --gene-coord-file"
        );
    };

    let col_name = derive_annot_name(&args);
    write_annot_file(&args.annot_file, &snps, &annotation, &col_name)
        .with_context(|| format!("writing annotation file '{}'", args.annot_file))?;

    let n_annotated = annotation.iter().filter(|&&v| v == 1).count();
    println!(
        "Wrote annotation '{}' to '{}': {}/{} SNPs annotated",
        col_name,
        args.annot_file,
        n_annotated,
        snps.len()
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// BED-file based annotation
// ---------------------------------------------------------------------------

/// Annotate SNPs using a UCSC BED file.
fn annotate_from_bed(
    snps: &[BimRecord],
    bed_path: &str,
    windowsize: u32,
    nomerge: bool,
) -> Result<Vec<u8>> {
    let intervals = load_bed_intervals(bed_path, windowsize, nomerge)?;
    println!(
        "  Loaded BED intervals for {} chromosomes from '{}'",
        intervals.len(),
        bed_path
    );
    Ok(snps.iter().map(|s| annotate_snp(s, &intervals)).collect())
}

/// Load BED intervals per chromosome.
///
/// Returns `HashMap<chr_string, sorted_and_merged_intervals>` where intervals
/// are stored as (start_0based_inclusive, end_0based_exclusive) tuples.
/// The `windowsize` is added symmetrically around each interval before merging.
fn load_bed_intervals(
    bed_path: &str,
    windowsize: u32,
    nomerge: bool,
) -> Result<HashMap<String, Vec<(u32, u32)>>> {
    let file = File::open(bed_path).with_context(|| format!("opening BED file '{}'", bed_path))?;
    let reader = BufReader::new(file);
    let mut intervals: HashMap<String, Vec<(u32, u32)>> = HashMap::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("reading BED line {}", i + 1))?;
        let line = line.trim();
        // Skip headers and empty lines.
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 {
            anyhow::bail!(
                "BED line {}: expected at least 3 columns (chrom start end), got {}",
                i + 1,
                cols.len()
            );
        }

        // Strip "chr" prefix so chromosome identifiers match PLINK .bim format.
        let chrom = cols[0].trim_start_matches("chr").to_string();
        let start: u32 = cols[1]
            .parse()
            .with_context(|| format!("BED line {}: invalid start '{}'", i + 1, cols[1]))?;
        let end: u32 = cols[2]
            .parse()
            .with_context(|| format!("BED line {}: invalid end '{}'", i + 1, cols[2]))?;

        // Expand by windowsize (saturating at 0 for start).
        let start = start.saturating_sub(windowsize);
        let end = end.saturating_add(windowsize);

        intervals.entry(chrom).or_default().push((start, end));
    }

    // Sort and optionally merge intervals per chromosome.
    for ivs in intervals.values_mut() {
        ivs.sort_by_key(|&(s, _)| s);
        if !nomerge {
            *ivs = merge_intervals(std::mem::take(ivs));
        }
    }

    Ok(intervals)
}

/// Merge overlapping or adjacent 0-based half-open intervals.
fn merge_intervals(sorted: Vec<(u32, u32)>) -> Vec<(u32, u32)> {
    if sorted.is_empty() {
        return sorted;
    }
    let mut merged: Vec<(u32, u32)> = Vec::with_capacity(sorted.len());
    let mut current = sorted[0];
    for &(s, e) in &sorted[1..] {
        if s <= current.1 {
            // Overlapping or adjacent: extend.
            current.1 = current.1.max(e);
        } else {
            merged.push(current);
            current = (s, e);
        }
    }
    merged.push(current);
    merged
}

/// Return 1 if the SNP's base-pair position (1-based) falls within any interval.
///
/// UCSC BED uses 0-based half-open coordinates [start, end).
/// A 1-based position `bp` is in the interval if: start <= bp−1 < end
/// i.e., start < bp  AND  bp <= end  (when end is the exclusive 0-based bound).
fn annotate_snp(snp: &BimRecord, intervals: &HashMap<String, Vec<(u32, u32)>>) -> u8 {
    let chr_str = snp.chr.to_string();
    if let Some(ivs) = intervals.get(&chr_str)
        && is_in_intervals(ivs, snp.bp)
    {
        return 1;
    }
    0
}

/// Binary-search check: is 1-based `bp` within any sorted 0-based half-open interval?
fn is_in_intervals(intervals: &[(u32, u32)], bp: u32) -> bool {
    // Find the rightmost interval whose start < bp (i.e., start <= bp-1 < bp).
    // partition_point returns the first index where the predicate is false.
    let pos = intervals.partition_point(|&(s, _)| s < bp);
    // The candidate interval is at pos-1 (if it exists).
    if pos > 0 {
        let (_, end) = intervals[pos - 1];
        bp <= end // end is exclusive 0-based → 1-based bp must be ≤ end
    } else {
        false
    }
}

// ---------------------------------------------------------------------------
// Gene-set based annotation
// ---------------------------------------------------------------------------

/// Annotate SNPs using a gene set (list of gene symbols) and a coordinate file.
fn annotate_from_gene_set(
    snps: &[BimRecord],
    gene_set_path: &str,
    coord_path: &str,
    windowsize: u32,
) -> Result<Vec<u8>> {
    // Load target gene symbols.
    let gene_set_file = File::open(gene_set_path)
        .with_context(|| format!("opening gene set file '{}'", gene_set_path))?;
    let target_genes: HashSet<String> = BufReader::new(gene_set_file)
        .lines()
        .map_while(Result::ok)
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .collect();
    println!(
        "  Loaded {} gene symbols from '{}'",
        target_genes.len(),
        gene_set_path
    );

    // Load gene coordinates, keeping only target genes.
    let coord_file = File::open(coord_path)
        .with_context(|| format!("opening gene coordinate file '{}'", coord_path))?;
    let mut intervals: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    let mut genes_found = 0usize;

    for (i, line) in BufReader::new(coord_file).lines().enumerate() {
        let line = line.with_context(|| format!("reading coord file line {}", i + 1))?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.len() < 4 {
            anyhow::bail!(
                "Gene coord file line {}: expected GENE CHR START END, got {} columns",
                i + 1,
                cols.len()
            );
        }

        let gene = cols[0];
        if !target_genes.contains(gene) {
            continue;
        }
        genes_found += 1;

        let chr = cols[1].trim_start_matches("chr").to_string();
        let start: u32 = cols[2]
            .parse()
            .with_context(|| format!("coord file line {}: invalid start '{}'", i + 1, cols[2]))?;
        let end: u32 = cols[3]
            .parse()
            .with_context(|| format!("coord file line {}: invalid end '{}'", i + 1, cols[3]))?;

        let start = start.saturating_sub(windowsize);
        let end = end.saturating_add(windowsize);

        intervals.entry(chr).or_default().push((start, end));
    }

    println!(
        "  Found coordinates for {}/{} target genes",
        genes_found,
        target_genes.len()
    );

    // Sort and merge per chromosome.
    for ivs in intervals.values_mut() {
        ivs.sort_by_key(|&(s, _)| s);
        *ivs = merge_intervals(std::mem::take(ivs));
    }

    Ok(snps.iter().map(|s| annotate_snp(s, &intervals)).collect())
}

// ---------------------------------------------------------------------------
// Output writer
// ---------------------------------------------------------------------------

fn derive_annot_name(args: &MakeAnnotArgs) -> String {
    let source = args
        .bed_file
        .as_deref()
        .or(args.gene_set_file.as_deref())
        .unwrap_or("annot");
    std::path::Path::new(source)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("ANNOT")
        .to_string()
}

/// Write the annotation file.
///
/// Header: `CHR  BP  SNP  CM  <col_name>` (matching Python make_annot.py).
/// If the path ends in `.gz` the output is gzip-compressed.
fn write_annot_file(
    path: &str,
    snps: &[BimRecord],
    annotation: &[u8],
    col_name: &str,
) -> Result<()> {
    if path.ends_with(".gz") {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        let file = File::create(path).with_context(|| format!("creating '{}'", path))?;
        let mut gz = GzEncoder::new(file, Compression::fast());
        write_rows(&mut gz, snps, annotation, col_name)?;
        gz.finish().context("finalizing gzip output")?;
    } else {
        let file = File::create(path).with_context(|| format!("creating '{}'", path))?;
        let mut w = BufWriter::new(file);
        write_rows(&mut w, snps, annotation, col_name)?;
    }
    Ok(())
}

fn write_rows(
    w: &mut impl Write,
    snps: &[BimRecord],
    annotation: &[u8],
    col_name: &str,
) -> Result<()> {
    writeln!(w, "CHR\tBP\tSNP\tCM\t{}", col_name)?;
    for (snp, &ann) in snps.iter().zip(annotation.iter()) {
        writeln!(
            w,
            "{}\t{}\t{}\t{:.6}\t{}",
            snp.chr, snp.bp, snp.snp, snp.cm, ann
        )?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merge_intervals() {
        // Non-overlapping: unchanged.
        let iv = vec![(0, 10), (20, 30)];
        assert_eq!(merge_intervals(iv), vec![(0, 10), (20, 30)]);

        // Overlapping: merged.
        let iv = vec![(0, 15), (10, 25), (30, 40)];
        assert_eq!(merge_intervals(iv), vec![(0, 25), (30, 40)]);

        // Adjacent: merged.
        let iv = vec![(0, 10), (10, 20)];
        assert_eq!(merge_intervals(iv), vec![(0, 20)]);
    }

    #[test]
    fn test_is_in_intervals() {
        // Intervals (0-based half-open): [100, 200) and [300, 400)
        let ivs = vec![(100, 200), (300, 400)];

        // 1-based bp=150: 0-based 149 → in [100,200) → yes (100 < 150 ≤ 200)
        assert!(is_in_intervals(&ivs, 150));

        // bp=200: 0-based 199 → in [100,200) → yes (100 < 200 ≤ 200)
        assert!(is_in_intervals(&ivs, 200));

        // bp=201: 0-based 200 → NOT in [100,200), NOT in [300,400)
        assert!(!is_in_intervals(&ivs, 201));

        // bp=100: 0-based 99 → NOT in [100,200) (100 is not < 100)
        assert!(!is_in_intervals(&ivs, 100));

        // bp=101: in [100,200) → yes
        assert!(is_in_intervals(&ivs, 101));
    }
}
