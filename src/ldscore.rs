/// LD score computation.
///
/// Memory: .bed files are memory-mapped via `bed-reader`; peak RAM is
/// `window_snps × n_indiv × 4 bytes` (≈ 20 MB for 1000G, ≈ 4 GB for biobanks).
///
/// Algorithm: for each SNP i, `L2[i] = Σ_j r²_unbiased(i, j)` where j ranges
/// over the window, using `r²_unbiased = r² − (1−r²)/(n−2)` [Bulik-Sullivan 2015].
/// SNPs are processed in chunks so the inner A^T B product is a single DGEMM.
///
/// Partitioned LD scores (`--annot`): `L2_k[i] = Σ_j r²_unbiased(i,j) × annot[j,k]`
/// yields one column per annotation; `.M`/`.M_5_50` files are tab-separated.
use anyhow::{Context, Result};
use bed_reader::{Bed, ReadOptions};
use ndarray::{Array1, Array2, Axis, ShapeBuilder, s};
use std::collections::{HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, Write as IoWrite};
use std::ops::AddAssign;

use crate::cli::LdscoreArgs;
use crate::parse;

// ---------------------------------------------------------------------------
// BIM / FAM parsing
// ---------------------------------------------------------------------------

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
fn load_individual_indices(keep_path: &str, fam_ids: &[(String, String)]) -> Result<Array1<isize>> {
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

    let indices: Array1<isize> = fam_ids
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

// ---------------------------------------------------------------------------
// Window boundary computation  (getBlockLefts from ldscore/ldscore.py)
// ---------------------------------------------------------------------------

/// Window mode for LD score computation.
#[derive(Debug, Clone, Copy)]
pub enum WindowMode {
    /// Window defined by genetic distance in cM.
    Cm(f64),
    /// Window defined by physical distance in kb.
    Kb(f64),
    /// Window defined by a fixed number of flanking SNPs.
    Snp(usize),
}

/// For each SNP i, compute the index of the leftmost SNP j such that
/// dist(coords[i], coords[j]) <= max_dist.
///
/// `coords` must be non-decreasing (sorted within a chromosome).
/// Returns a `Vec<usize>` where `block_left[i]` is that leftmost index.
fn get_block_lefts_f64(coords: &[f64], max_dist: f64) -> Vec<usize> {
    let m = coords.len();
    let mut block_left = vec![0usize; m];
    let mut j = 0usize;

    for i in 0..m {
        while j < m && (coords[i] - coords[j]).abs() > max_dist {
            j += 1;
        }
        block_left[i] = j;
    }

    block_left
}

/// Fixed-SNP window: block_left[i] = max(0, i - half_window).
fn get_block_lefts_snp(m: usize, half_window: usize) -> Vec<usize> {
    (0..m).map(|i| i.saturating_sub(half_window)).collect()
}

// ---------------------------------------------------------------------------
// Genotype normalisation helpers
// ---------------------------------------------------------------------------

/// Normalize a genotype column in-place: impute NaN with mean, centre, scale to unit variance.
///
/// Returns the minor allele frequency: `MAF = min(mean/2, 1 − mean/2)`.
/// Works on f32; caller casts to f64 into the pre-allocated `b_mat`.
fn normalize_col(col: &mut Array1<f32>, n: usize) -> f64 {
    let (sum, count) = col.iter().fold((0f64, 0usize), |(s, c), &v| {
        if v.is_nan() {
            (s, c)
        } else {
            (s + v as f64, c + 1)
        }
    });
    let avg = if count > 0 { sum / count as f64 } else { 0.0 };
    let freq = (avg / 2.0).clamp(0.0, 1.0);
    let maf = freq.min(1.0 - freq);

    for v in col.iter_mut() {
        if v.is_nan() {
            *v = avg as f32;
        } else {
            *v -= avg as f32;
        }
    }

    let var: f64 = col.iter().map(|&v| (v as f64).powi(2)).sum::<f64>() / n as f64;
    let std = var.sqrt();
    if std > 0.0 {
        let inv_std = (1.0 / std) as f32;
        col.iter_mut().for_each(|v| *v *= inv_std);
    }
    maf
}

// ---------------------------------------------------------------------------
// LD score computation
// ---------------------------------------------------------------------------

/// Compute LD scores for all SNPs in one sequential global pass.
///
/// Returns `(l2, maf)` where `l2` has shape `(m, 1)` (scalar) or `(m, K)`
/// (partitioned with annotation), and `maf` is per-SNP MAF.
#[allow(clippy::too_many_arguments)]
fn compute_ldscore_global(
    all_snps: &[BimRecord],
    bed: &mut Bed,
    n_indiv: usize,
    mode: WindowMode,
    chunk_c: usize,
    annot: Option<&Array2<f64>>,
    iid_indices: Option<&Array1<isize>>,
    pq_exp: Option<f64>,
    yes_really: bool,
) -> Result<(Array2<f64>, Vec<f64>)> {
    let m = all_snps.len();
    if m == 0 {
        let n_annot = annot.map(|a| a.ncols()).unwrap_or(1);
        return Ok((Array2::zeros((0, n_annot)), vec![]));
    }

    // Global block_left: --ld-wind-snps uses integer coordinates 0..m-1.
    let block_left: Vec<usize> = match mode {
        WindowMode::Snp(half) => get_block_lefts_snp(m, half),
        WindowMode::Cm(max_cm) => {
            let coords: Vec<f64> = all_snps.iter().map(|s| s.cm).collect();
            get_block_lefts_f64(&coords, max_cm)
        }
        WindowMode::Kb(max_kb) => {
            let coords: Vec<f64> = all_snps.iter().map(|s| s.bp as f64 / 1000.0).collect();
            get_block_lefts_f64(&coords, max_kb)
        }
    };
    if !yes_really && !block_left.is_empty() && block_left[m - 1] == 0 {
        println!(
            "WARNING: LD window spans the entire chromosome. \
             Use --yes-really to silence this warning."
        );
    }

    let n = n_indiv as f64;
    let mut maf_per_snp = vec![0.0f64; m];

    let max_window_size = (0..m).map(|i| i - block_left[i]).max().unwrap_or(0);

    // Ring buffer: ring_size = max_window + chunk_c guarantees no live slot is overwritten.
    let ring_size = (max_window_size + chunk_c).max(1);
    let mut ring_buf = Array2::<f64>::zeros((n_indiv, ring_size).f()); // F-order: columns contiguous
    let mut ring_next: usize = 0;
    let mut window: VecDeque<(usize, usize)> = VecDeque::new(); // (snp_idx, ring_slot)

    let mut b_mat = Array2::<f64>::zeros((n_indiv, chunk_c)); // C-order for BLAS
    let mut col_buf = Array1::<f32>::zeros(n_indiv);
    let mut a_buf = Array2::<f64>::zeros((n_indiv, max_window_size.max(1)).f()); // F-order

    if annot.is_none() {
        let mut l2_scalar = vec![0.0f64; m];

        for chunk_start in (0..m).step_by(chunk_c) {
            let chunk_end = (chunk_start + chunk_c).min(m);
            let c = chunk_end - chunk_start;

            let a_left = block_left[chunk_start];
            while window
                .front()
                .map(|(idx, _)| *idx < a_left)
                .unwrap_or(false)
            {
                window.pop_front();
            }

            let bed_indices: Array1<isize> = Array1::from_iter(
                all_snps[chunk_start..chunk_end]
                    .iter()
                    .map(|s| s.bed_idx as isize),
            );
            let raw: Array2<f32> = {
                let mut ro = ReadOptions::builder();
                let builder = ro.sid_index(&bed_indices).f32();
                let result = if let Some(iids) = iid_indices {
                    builder.iid_index(iids).read(bed)
                } else {
                    builder.read(bed)
                };
                result
                    .with_context(|| format!("reading BED chunk [{},{})", chunk_start, chunk_end))?
            };

            {
                let mut bv = b_mat.slice_mut(s![.., ..c]);
                for j in 0..c {
                    col_buf.assign(&raw.column(j));
                    let snp_maf = normalize_col(&mut col_buf, n_indiv);
                    maf_per_snp[chunk_start + j] = snp_maf;
                    ndarray::Zip::from(bv.column_mut(j))
                        .and(col_buf.view())
                        .for_each(|d, &s| *d = s as f64);
                }
            }
            let b_slice = b_mat.slice(s![.., ..c]);

            // B×B
            let bb: Array2<f64> = b_slice.t().dot(&b_slice);
            for j in 0..c {
                let j_g = chunk_start + j;
                let pq_j = if let Some(exp) = pq_exp {
                    let p = maf_per_snp[j_g];
                    (p * (1.0 - p)).powf(exp)
                } else {
                    1.0
                };
                l2_scalar[j_g] += pq_j; // diagonal: r(j,j) = 1
                for k in 0..j {
                    let k_g = chunk_start + k;
                    let pq_k = if let Some(exp) = pq_exp {
                        let p = maf_per_snp[k_g];
                        (p * (1.0 - p)).powf(exp)
                    } else {
                        1.0
                    };
                    let r2u = r2_unbiased(bb[[k, j]] / n, n_indiv);
                    l2_scalar[j_g] += r2u * pq_k;
                    l2_scalar[k_g] += r2u * pq_j;
                }
            }

            // A×B
            if !window.is_empty() {
                let w = window.len();
                for (wi, (_, slot)) in window.iter().enumerate() {
                    a_buf.column_mut(wi).assign(&ring_buf.column(*slot));
                }
                let ab: Array2<f64> = a_buf.slice(s![.., ..w]).t().dot(&b_slice);
                for (wi, (k_g, _)) in window.iter().enumerate() {
                    let pq_k = if let Some(exp) = pq_exp {
                        let p = maf_per_snp[*k_g];
                        (p * (1.0 - p)).powf(exp)
                    } else {
                        1.0
                    };
                    for j in 0..c {
                        let j_g = chunk_start + j;
                        let pq_j = if let Some(exp) = pq_exp {
                            let p = maf_per_snp[j_g];
                            (p * (1.0 - p)).powf(exp)
                        } else {
                            1.0
                        };
                        let r2u = r2_unbiased(ab[[wi, j]] / n, n_indiv);
                        l2_scalar[j_g] += r2u * pq_k;
                        l2_scalar[*k_g] += r2u * pq_j;
                    }
                }
            }

            for j in 0..c {
                let slot = ring_next % ring_size;
                ring_buf.column_mut(slot).assign(&b_slice.column(j));
                window.push_back((chunk_start + j, slot));
                ring_next += 1;
            }
        }

        let l2 = Array2::from_shape_vec((m, 1), l2_scalar).expect("shape matches");
        return Ok((l2, maf_per_snp));
    }

    // Partitioned path (annot is Some).
    let annot = annot.unwrap();
    anyhow::ensure!(
        annot.nrows() == m,
        "Annotation matrix has {} rows but BIM has {} SNPs",
        annot.nrows(),
        m
    );
    let n_annot = annot.ncols();
    let mut l2 = Array2::<f64>::zeros((m, n_annot));

    for chunk_start in (0..m).step_by(chunk_c) {
        let chunk_end = (chunk_start + chunk_c).min(m);
        let c = chunk_end - chunk_start;

        let a_left = block_left[chunk_start];
        while window
            .front()
            .map(|(idx, _)| *idx < a_left)
            .unwrap_or(false)
        {
            window.pop_front();
        }

        let bed_indices: Array1<isize> = Array1::from_iter(
            all_snps[chunk_start..chunk_end]
                .iter()
                .map(|s| s.bed_idx as isize),
        );
        let raw: Array2<f32> = {
            let mut ro = ReadOptions::builder();
            let builder = ro.sid_index(&bed_indices).f32();
            let result = if let Some(iids) = iid_indices {
                builder.iid_index(iids).read(bed)
            } else {
                builder.read(bed)
            };
            result.with_context(|| format!("reading BED chunk [{},{})", chunk_start, chunk_end))?
        };

        {
            let mut bv = b_mat.slice_mut(s![.., ..c]);
            for j in 0..c {
                col_buf.assign(&raw.column(j));
                let snp_maf = normalize_col(&mut col_buf, n_indiv);
                maf_per_snp[chunk_start + j] = snp_maf;
                ndarray::Zip::from(bv.column_mut(j))
                    .and(col_buf.view())
                    .for_each(|d, &s| *d = s as f64);
            }
        }
        let b_slice = b_mat.slice(s![.., ..c]);

        let annot_chunk = annot.slice(s![chunk_start..chunk_end, ..]);

        // Scale annotation rows by (p·(1−p))^S for --pq-exp / --per-allele.
        let annot_chunk_eff: Array2<f64> = if let Some(exp) = pq_exp {
            let mut eff = annot_chunk.to_owned();
            for j in 0..c {
                let p = maf_per_snp[chunk_start + j];
                let pq = (p * (1.0 - p)).powf(exp);
                eff.row_mut(j).mapv_inplace(|v| v * pq);
            }
            eff
        } else {
            annot_chunk.to_owned()
        };

        // B×B
        let bb: Array2<f64> = b_slice.t().dot(&b_slice);
        let mut r2u_bb = Array2::<f64>::zeros((c, c));
        for j in 0..c {
            r2u_bb[[j, j]] = 1.0; // diagonal: r(j,j)=1 → r2u=1
            for k in 0..j {
                let r2u = r2_unbiased(bb[[k, j]] / n, n_indiv);
                r2u_bb[[j, k]] = r2u;
                r2u_bb[[k, j]] = r2u;
            }
        }
        let contrib_bb = r2u_bb.dot(&annot_chunk_eff);
        l2.slice_mut(s![chunk_start..chunk_end, ..])
            .add_assign(&contrib_bb);

        // A×B
        if !window.is_empty() {
            let w = window.len();
            let a_left_idx = window.front().map(|(idx, _)| *idx).unwrap_or(0);

            for (wi, (_, slot)) in window.iter().enumerate() {
                a_buf.column_mut(wi).assign(&ring_buf.column(*slot));
            }
            let ab: Array2<f64> = a_buf.slice(s![.., ..w]).t().dot(&b_slice);

            let mut r2u_ab = Array2::<f64>::zeros((w, c));
            for wi in 0..w {
                for j in 0..c {
                    r2u_ab[[wi, j]] = r2_unbiased(ab[[wi, j]] / n, n_indiv);
                }
            }

            let annot_window = annot.slice(s![a_left_idx..chunk_start, ..]);
            let annot_window_eff: Array2<f64> = if let Some(exp) = pq_exp {
                let mut eff = annot_window.to_owned();
                for (wi, (w_g, _)) in window.iter().enumerate() {
                    let p = maf_per_snp[*w_g];
                    let pq = (p * (1.0 - p)).powf(exp);
                    eff.row_mut(wi).mapv_inplace(|v| v * pq);
                }
                eff
            } else {
                annot_window.to_owned()
            };

            l2.slice_mut(s![a_left_idx..chunk_start, ..])
                .add_assign(&r2u_ab.dot(&annot_chunk_eff));
            l2.slice_mut(s![chunk_start..chunk_end, ..])
                .add_assign(&r2u_ab.t().dot(&annot_window_eff));
        }

        for j in 0..c {
            let slot = ring_next % ring_size;
            ring_buf.column_mut(slot).assign(&b_slice.column(j));
            window.push_back((chunk_start + j, slot));
            ring_next += 1;
        }
    }

    Ok((l2, maf_per_snp))
}

fn compute_maf_only(
    all_snps: &[BimRecord],
    bed: &mut Bed,
    n_indiv: usize,
    chunk_c: usize,
    iid_indices: Option<&Array1<isize>>,
) -> Result<Vec<f64>> {
    let m = all_snps.len();
    if m == 0 {
        return Ok(vec![]);
    }

    let mut maf_per_snp = vec![0.0f64; m];
    let mut col_buf = Array1::<f32>::zeros(n_indiv);

    for chunk_start in (0..m).step_by(chunk_c) {
        let chunk_end = (chunk_start + chunk_c).min(m);
        let c = chunk_end - chunk_start;

        let bed_indices: Array1<isize> = Array1::from_iter(
            all_snps[chunk_start..chunk_end]
                .iter()
                .map(|s| s.bed_idx as isize),
        );
        let raw: Array2<f32> = {
            let mut ro = ReadOptions::builder();
            let builder = ro.sid_index(&bed_indices).f32();
            let result = if let Some(iids) = iid_indices {
                builder.iid_index(iids).read(bed)
            } else {
                builder.read(bed)
            };
            result.with_context(|| format!("reading BED chunk [{},{})", chunk_start, chunk_end))?
        };

        for j in 0..c {
            col_buf.assign(&raw.column(j));
            let snp_maf = normalize_col(&mut col_buf, n_indiv);
            maf_per_snp[chunk_start + j] = snp_maf;
        }
    }

    Ok(maf_per_snp)
}

/// Unbiased r² estimator: r² − (1−r²)/(n−2) [Bulik-Sullivan 2015].
#[inline]
fn r2_unbiased(r: f64, n: usize) -> f64 {
    let sq = r * r;
    if n > 2 {
        sq - (1.0 - sq) / (n as f64 - 2.0)
    } else {
        sq
    }
}

// ---------------------------------------------------------------------------
// Output writer
// ---------------------------------------------------------------------------

/// Write per-SNP LD scores to a gzip TSV.  `col_names` is `["L2"]` or `["{ANNOT}L2", …]`.
fn write_ldscore_refs(
    path: &str,
    snps: &[&BimRecord],
    l2: &Array2<f64>,
    col_names: &[String],
) -> Result<()> {
    use flate2::Compression;
    use flate2::write::GzEncoder;

    let file = File::create(path).with_context(|| format!("creating output '{}'", path))?;
    let mut gz = GzEncoder::new(file, Compression::fast());

    write!(gz, "CHR\tSNP\tBP")?;
    for name in col_names {
        write!(gz, "\t{}", name)?;
    }
    writeln!(gz)?;

    for (i, snp) in snps.iter().enumerate() {
        write!(gz, "{}\t{}\t{}", snp.chr, snp.snp, snp.bp)?;
        for k in 0..col_names.len() {
            write!(gz, "\t{:.6}", l2[[i, k]])?;
        }
        writeln!(gz)?;
    }

    gz.finish().context("finalising gzip output")?;
    Ok(())
}

fn format_m_vals(vals: &[f64]) -> String {
    let joined = vals
        .iter()
        .map(|v| format!("{v}"))
        .collect::<Vec<_>>()
        .join("\t");
    format!("{joined}\n")
}

// ---------------------------------------------------------------------------
// SNP set loader (--extract / --print-snps)
// ---------------------------------------------------------------------------

/// Load SNP IDs from a file (one per line; first token used; `#` lines skipped).
fn load_snp_set(path: &str) -> Result<HashSet<String>> {
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

pub fn run(args: LdscoreArgs) -> Result<()> {
    if args.no_print_annot {
        println!(
            "WARNING: --no-print-annot is a no-op in the Rust port; \
             cts-annot always writes an annot file."
        );
    }
    if args.per_allele && args.pq_exp.is_some() {
        anyhow::bail!(
            "Cannot set both --per-allele and --pq-exp (--per-allele is equivalent to --pq-exp 1)."
        );
    }
    let pq_exp = if args.per_allele {
        Some(1.0)
    } else {
        args.pq_exp
    };
    let mode = if let Some(kb) = args.ld_wind_kb {
        WindowMode::Kb(kb)
    } else if let Some(snp) = args.ld_wind_snp {
        WindowMode::Snp(snp)
    } else {
        WindowMode::Cm(args.ld_wind_cm)
    };

    let bim_path = format!("{}.bim", args.bfile);
    let fam_path = format!("{}.fam", args.bfile);
    let bed_path = format!("{}.bed", args.bfile);

    let all_snps_raw =
        parse_bim(&bim_path).with_context(|| format!("parsing BIM '{}'", bim_path))?;
    let n_indiv = count_fam(&fam_path).with_context(|| format!("counting FAM '{}'", fam_path))?;

    println!(
        "Loaded {} SNPs, {} individuals from '{}'",
        all_snps_raw.len(),
        n_indiv,
        args.bfile
    );

    // --extract: pre-filter BIM (also shrinks LD windows).
    let all_snps: Vec<BimRecord> = if let Some(ref extract_path) = args.extract {
        let extract_set = load_snp_set(extract_path)
            .with_context(|| format!("loading --extract file '{}'", extract_path))?;
        let filtered: Vec<BimRecord> = all_snps_raw
            .into_iter()
            .filter(|s| extract_set.contains(&s.snp))
            .collect();
        println!("  After --extract: {} SNPs", filtered.len());
        filtered
    } else {
        all_snps_raw
    };

    // --print-snps: output filter only; all SNPs still in LD windows.
    let print_set: Option<HashSet<String>> = if let Some(ref ps_path) = args.print_snps {
        Some(
            load_snp_set(ps_path)
                .with_context(|| format!("loading --print-snps file '{}'", ps_path))?,
        )
    } else {
        None
    };

    if args.annot.is_some() && args.extract.is_some() {
        println!(
            "WARNING: --annot with --extract is not supported in Python LDSC. \
             Ensure your annot files match the extracted SNP set."
        );
    }

    let mut annot_result: Option<(Array2<f64>, Vec<String>)> = if let Some(ref prefix) = args.annot
    {
        // Collect unique chromosomes in BIM order.
        let mut chrs_seen: HashSet<u8> = HashSet::new();
        let mut chrs: Vec<u8> = Vec::new();
        for snp in &all_snps {
            if chrs_seen.insert(snp.chr) {
                chrs.push(snp.chr);
            }
        }
        chrs.sort_unstable();

        // Read per-chromosome annotation files (prefix{chr}.annot[.gz]).
        let mut mats: Vec<Array2<f64>> = Vec::new();
        let mut col_names: Vec<String> = Vec::new();
        for chr in &chrs {
            let chr_prefix = format!("{}{}", prefix, chr);
            let (mat, names) =
                parse::read_annot(&chr_prefix, args.thin_annot).with_context(|| {
                    format!("reading annotation for chr{} (prefix '{}')", chr, prefix)
                })?;
            if col_names.is_empty() {
                col_names = names;
            }
            mats.push(mat);
        }

        let views: Vec<_> = mats.iter().map(|m| m.view()).collect();
        let combined = ndarray::concatenate(Axis(0), &views)
            .context("concatenating per-chromosome annotation matrices")?;
        anyhow::ensure!(
            combined.nrows() == all_snps.len(),
            "Annotation file has {} rows but BIM has {} SNPs — they must match exactly",
            combined.nrows(),
            all_snps.len()
        );
        println!(
            "Read {} annotations for {} SNPs from '{}*'",
            col_names.len(),
            combined.nrows(),
            prefix
        );
        Some((combined, col_names))
    } else {
        None
    };

    // --keep: subset individuals for LD computation.
    let iid_indices: Option<Array1<isize>> = if let Some(ref keep_path) = args.keep {
        let fam_ids =
            parse_fam(&fam_path).with_context(|| format!("parsing FAM '{}'", fam_path))?;
        Some(
            load_individual_indices(keep_path, &fam_ids)
                .with_context(|| format!("loading keep file '{}'", keep_path))?,
        )
    } else {
        None
    };
    let n_indiv_actual = iid_indices.as_ref().map(|idx| idx.len()).unwrap_or(n_indiv);

    let mut all_snps = all_snps;

    // Optional pre-filter on MAF (Python behavior).
    if args.maf_pre {
        if let Some(thr) = args.maf {
            let mut bed =
                Bed::new(bed_path.as_str()).context("opening BED file for MAF prefilter")?;
            let maf_all = compute_maf_only(
                &all_snps,
                &mut bed,
                n_indiv_actual,
                args.chunk_size,
                iid_indices.as_ref(),
            )
            .context("computing MAF prefilter")?;

            let mut kept_snps: Vec<BimRecord> = Vec::new();
            let mut keep_mask: Vec<bool> = Vec::with_capacity(all_snps.len());
            for (snp, maf) in all_snps.iter().zip(maf_all.iter()) {
                let keep = *maf >= thr;
                keep_mask.push(keep);
                if keep {
                    kept_snps.push(snp.clone());
                }
            }
            if kept_snps.is_empty() {
                anyhow::bail!("--maf-pre removed all SNPs (threshold {})", thr);
            }
            all_snps = kept_snps;

            if let Some((annot, names)) = annot_result.take() {
                let rows: Vec<usize> = keep_mask
                    .iter()
                    .enumerate()
                    .filter_map(|(i, &k)| if k { Some(i) } else { None })
                    .collect();
                let mut filtered = Array2::<f64>::zeros((rows.len(), annot.ncols()));
                for (ri, &src) in rows.iter().enumerate() {
                    filtered.row_mut(ri).assign(&annot.row(src));
                }
                annot_result = Some((filtered, names));
            }

            println!(
                "--maf-pre: kept {} / {} SNPs (MAF >= {})",
                all_snps.len(),
                keep_mask.len(),
                thr
            );
        } else {
            println!("WARNING: --maf-pre has no effect without --maf");
        }
    }

    let annot_ref: Option<&Array2<f64>> = annot_result.as_ref().map(|(m, _)| m);

    let chunk_c = args.chunk_size;

    let mut bed = Bed::new(bed_path.as_str()).context("opening BED file")?;
    let (l2, maf_per_snp) = compute_ldscore_global(
        &all_snps,
        &mut bed,
        n_indiv_actual,
        mode,
        chunk_c,
        annot_ref,
        iid_indices.as_ref(),
        pq_exp,
        args.yes_really,
    )
    .context("computing LD scores")?;

    // bed_idx (original BIM row) != position in all_snps when --extract is active.
    let bed_idx_to_pos: std::collections::HashMap<usize, usize> = all_snps
        .iter()
        .enumerate()
        .map(|(i, s)| (s.bed_idx, i))
        .collect();

    let scale_suffix = pq_exp.map(|exp| format!("_S{}", exp)).unwrap_or_default();
    let col_names: Vec<String> = match &annot_result {
        None => vec![format!("L2{}", scale_suffix)],
        Some((_, names)) => names
            .iter()
            .map(|n| format!("{}L2{}", n, scale_suffix))
            .collect(),
    };

    let pq_per_snp: Option<Vec<f64>> = pq_exp.map(|exp| {
        maf_per_snp
            .iter()
            .map(|&p| (p * (1.0 - p)).powf(exp))
            .collect()
    });

    let mut chrs: Vec<u8> = all_snps
        .iter()
        .map(|s| s.chr)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    chrs.sort();

    for chr in chrs {
        let chr_snps_all: Vec<&BimRecord> = all_snps.iter().filter(|s| s.chr == chr).collect();

        let chr_snps: Vec<&BimRecord> = chr_snps_all
            .iter()
            .copied()
            .filter(|s| {
                let pos = bed_idx_to_pos[&s.bed_idx];
                let maf_ok = args.maf.map(|thr| maf_per_snp[pos] >= thr).unwrap_or(true);
                let print_ok = print_set
                    .as_ref()
                    .map(|set| set.contains(&s.snp))
                    .unwrap_or(true);
                maf_ok && print_ok
            })
            .collect();
        let n_chr = chr_snps.len();

        let chr_l2: Array2<f64> = {
            let mut mat = Array2::<f64>::zeros((n_chr, col_names.len()));
            for (row, s) in chr_snps.iter().enumerate() {
                let pos = bed_idx_to_pos[&s.bed_idx];
                mat.row_mut(row).assign(&l2.row(pos));
            }
            mat
        };

        let out_path = format!("{}{}.l2.ldscore.gz", args.out, chr);
        write_ldscore_refs(&out_path, &chr_snps, &chr_l2, &col_names)
            .with_context(|| format!("writing output for chr {}", chr))?;

        // M: all extracted SNPs; M_5_50: MAF ≥ 0.05; with --annot: sum of annot[:,k].
        let chr_maf_all: Vec<f64> = chr_snps_all
            .iter()
            .map(|s| maf_per_snp[bed_idx_to_pos[&s.bed_idx]])
            .collect();

        let m_vals: Vec<f64> = match &annot_result {
            Some((annot, _)) => (0..col_names.len())
                .map(|k| {
                    chr_snps_all
                        .iter()
                        .map(|s| {
                            let pos = bed_idx_to_pos[&s.bed_idx];
                            let base = annot[[pos, k]];
                            if let Some(pq) = pq_per_snp.as_ref() {
                                base * pq[pos]
                            } else {
                                base
                            }
                        })
                        .sum::<f64>()
                })
                .collect(),
            None => {
                if let Some(pq) = pq_per_snp.as_ref() {
                    vec![
                        chr_snps_all
                            .iter()
                            .map(|s| pq[bed_idx_to_pos[&s.bed_idx]])
                            .sum::<f64>(),
                    ]
                } else {
                    vec![chr_snps_all.len() as f64]
                }
            }
        };

        let m_5_50_vals: Vec<f64> = match &annot_result {
            Some((annot, _)) => (0..col_names.len())
                .map(|k| {
                    chr_snps_all
                        .iter()
                        .zip(chr_maf_all.iter())
                        .filter(|(_, maf)| **maf >= 0.05)
                        .map(|(s, _)| {
                            let pos = bed_idx_to_pos[&s.bed_idx];
                            let base = annot[[pos, k]];
                            if let Some(pq) = pq_per_snp.as_ref() {
                                base * pq[pos]
                            } else {
                                base
                            }
                        })
                        .sum::<f64>()
                })
                .collect(),
            None => {
                if let Some(pq) = pq_per_snp.as_ref() {
                    vec![
                        chr_snps_all
                            .iter()
                            .zip(chr_maf_all.iter())
                            .filter(|(_, maf)| **maf >= 0.05)
                            .map(|(s, _)| pq[bed_idx_to_pos[&s.bed_idx]])
                            .sum::<f64>(),
                    ]
                } else {
                    let m_5_50 = chr_maf_all
                        .iter()
                        .filter(|&&m| (0.05..=0.5).contains(&m))
                        .count();
                    vec![m_5_50 as f64]
                }
            }
        };

        let m_path = format!("{}{}.l2.M", args.out, chr);
        std::fs::write(&m_path, format_m_vals(&m_vals))
            .with_context(|| format!("writing M file '{}'", m_path))?;

        let m_5_50_path = format!("{}{}.l2.M_5_50", args.out, chr);
        std::fs::write(&m_5_50_path, format_m_vals(&m_5_50_vals))
            .with_context(|| format!("writing M_5_50 file '{}'", m_5_50_path))?;

        println!(
            "chr {}: {} SNPs → {} (M={}, M_5_50={})",
            chr,
            n_chr,
            out_path,
            m_vals
                .iter()
                .map(|v| format!("{v}"))
                .collect::<Vec<_>>()
                .join(","),
            m_5_50_vals
                .iter()
                .map(|v| format!("{v}"))
                .collect::<Vec<_>>()
                .join(","),
        );
    }

    Ok(())
}

/// Strand-ambiguous allele pairs (A/T, C/G) — excluded from regression by default.
#[allow(dead_code)]
pub fn is_strand_ambiguous(a1: &str, a2: &str) -> bool {
    matches!((a1, a2), ("A", "T") | ("T", "A") | ("C", "G") | ("G", "C"))
}

/// Return the complement nucleotide.
#[allow(dead_code)]
pub fn complement(allele: &str) -> Option<&'static str> {
    match allele {
        "A" => Some("T"),
        "T" => Some("A"),
        "C" => Some("G"),
        "G" => Some("C"),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -- BIM parser -----------------------------------------------------------

    /// parse_bim should correctly round-trip a 6-column tab-separated line.
    #[test]
    fn test_parse_bim_inline() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "1\trs123\t0.5\t12345\tA\tG").unwrap();
        writeln!(f, "1\trs456\t1.0\t99999\tC\tT").unwrap();

        let recs = parse_bim(f.path().to_str().unwrap()).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].snp, "rs123");
        assert!((recs[0].cm - 0.5).abs() < 1e-9);
        assert_eq!(recs[0].bp, 12345);
        assert_eq!(recs[0].bed_idx, 0);

        assert_eq!(recs[1].snp, "rs456");
        assert_eq!(recs[1].bed_idx, 1);
    }

    // -- Window boundary computation -----------------------------------------

    #[test]
    fn test_block_lefts_f64_basic() {
        let coords = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let bl = get_block_lefts_f64(&coords, 1.5);
        assert_eq!(bl, vec![0, 0, 1, 2, 3]);

        let bl2 = get_block_lefts_f64(&coords, 2.5);
        assert_eq!(bl2, vec![0, 0, 0, 1, 2]);
    }

    #[test]
    fn test_block_lefts_snp_basic() {
        let bl = get_block_lefts_snp(6, 2);
        assert_eq!(bl, vec![0, 0, 0, 1, 2, 3]);
    }

    // -- Column normalisation ------------------------------------------------

    #[test]
    fn test_normalize_col_constant() {
        let mut col = Array1::from_vec(vec![2.0f32, 2.0, 2.0, 2.0]);
        normalize_col(&mut col, 4);
        assert!(
            col.iter().all(|&v| v == 0.0),
            "constant column should be all zeros"
        );
    }

    #[test]
    fn test_normalize_col_unit_variance() {
        let mut col = Array1::from_vec(vec![0.0f32, 1.0, 2.0, 1.0]);
        let n = col.len();
        normalize_col(&mut col, n);

        let mean: f64 = col.iter().map(|&v| v as f64).sum::<f64>() / n as f64;
        let var: f64 = col.iter().map(|&v| (v as f64 - mean).powi(2)).sum::<f64>() / n as f64;
        assert!(mean.abs() < 1e-5, "mean should be ~0, got {}", mean);
        assert!(
            (var - 1.0).abs() < 1e-5,
            "variance should be ~1, got {}",
            var
        );
    }

    // -- Allele helpers ------------------------------------------------------

    #[test]
    fn test_is_strand_ambiguous() {
        assert!(is_strand_ambiguous("A", "T"));
        assert!(is_strand_ambiguous("C", "G"));
        assert!(!is_strand_ambiguous("A", "G"));
    }

    #[test]
    fn test_complement() {
        assert_eq!(complement("A"), Some("T"));
        assert_eq!(complement("X"), None);
    }

    // -- Python test_getBlockLefts exact cases --------------------------------

    /// All SNPs within one giant window → all block_lefts = 0.
    /// Mirrors Python: getBlockLefts(arange(1,6), 5) == zeros(5).
    #[test]
    fn test_block_lefts_all_in_window() {
        let coords = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let bl = get_block_lefts_f64(&coords, 5.0);
        assert_eq!(bl, vec![0, 0, 0, 0, 0]);
    }

    /// Zero-width window → each SNP is its own left boundary.
    /// Mirrors Python: getBlockLefts(arange(1,6), 0) == arange(0,5).
    #[test]
    fn test_block_lefts_zero_window() {
        let coords = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let bl = get_block_lefts_f64(&coords, 0.0);
        assert_eq!(bl, vec![0, 1, 2, 3, 4]);
    }

    /// Irregular spacing with exact boundary conditions.
    /// Mirrors Python: getBlockLefts((1,4,6,7,7,8), 2) == (0,1,1,2,2,2).
    #[test]
    fn test_block_lefts_irregular_spacing() {
        let coords = vec![1.0, 4.0, 6.0, 7.0, 7.0, 8.0];
        let bl = get_block_lefts_f64(&coords, 2.0);
        assert_eq!(bl, vec![0, 1, 1, 2, 2, 2]);
    }

    fn write_plink_small(dir: &std::path::Path) -> String {
        use std::io::Write;
        let prefix = dir.join("toy");
        let bim = prefix.with_extension("bim");
        let fam = prefix.with_extension("fam");
        let bed = prefix.with_extension("bed");

        // .fam: 4 individuals
        let mut f = std::fs::File::create(&fam).unwrap();
        for i in 1..=4 {
            writeln!(f, "F{} I{} 0 0 0 -9", i, i).unwrap();
        }

        // .bim: 2 SNPs
        let mut b = std::fs::File::create(&bim).unwrap();
        writeln!(b, "1\trs1\t0\t100\tA\tG").unwrap();
        writeln!(b, "1\trs2\t0\t200\tA\tG").unwrap();

        // .bed: SNP-major
        let mut bed_f = std::fs::File::create(&bed).unwrap();
        bed_f.write_all(&[0x6C, 0x1B, 0x01]).unwrap();

        // SNP1: all hom-major (MAF=0)
        let snp1 = [0u8, 0u8, 0u8, 0u8];
        bed_f.write_all(&[pack_genotypes(&snp1)]).unwrap();

        // SNP2: genotypes [0,1,1,2] => MAF=0.5
        let snp2 = [0u8, 1u8, 1u8, 2u8];
        bed_f.write_all(&[pack_genotypes(&snp2)]).unwrap();

        prefix.to_string_lossy().to_string()
    }

    fn pack_genotypes(gt: &[u8; 4]) -> u8 {
        // PLINK .bed encoding (2-bit, little-endian per individual):
        // 00 = hom major, 01 = missing, 10 = het, 11 = hom minor
        let code = |g: u8| match g {
            0 => 0b00,
            1 => 0b10,
            2 => 0b11,
            _ => 0b01,
        };
        code(gt[0]) | (code(gt[1]) << 2) | (code(gt[2]) << 4) | (code(gt[3]) << 6)
    }

    fn read_gz(path: &str) -> String {
        use flate2::read::GzDecoder;
        use std::io::Read;
        let file = std::fs::File::open(path).expect("open gz");
        let mut decoder = GzDecoder::new(file);
        let mut s = String::new();
        decoder.read_to_string(&mut s).expect("read gz");
        s
    }

    #[test]
    fn test_maf_pre_changes_m_counts() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = write_plink_small(dir.path());
        let out_dir = tempfile::tempdir().unwrap();
        let out = out_dir.path().join("out").to_string_lossy().to_string();

        let args = LdscoreArgs {
            bfile: prefix,
            out: out.clone(),
            ld_wind_cm: 1.0,
            ld_wind_kb: None,
            ld_wind_snp: Some(1),
            annot: None,
            thin_annot: false,
            extract: None,
            maf: Some(0.1),
            maf_pre: true,
            print_snps: None,
            keep: None,
            per_allele: false,
            pq_exp: None,
            no_print_annot: false,
            chunk_size: 50,
            yes_really: true,
        };
        run(args).unwrap();

        let m_path = format!("{}1.l2.M", out);
        let m = std::fs::read_to_string(m_path).unwrap();
        assert_eq!(m.trim(), "1", "pre-filter should keep 1 SNP");
    }

    #[test]
    fn test_pq_exp_writes_scaled_m_and_suffix() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = write_plink_small(dir.path());
        let out_dir = tempfile::tempdir().unwrap();
        let out = out_dir.path().join("out").to_string_lossy().to_string();

        let args = LdscoreArgs {
            bfile: prefix,
            out: out.clone(),
            ld_wind_cm: 1.0,
            ld_wind_kb: None,
            ld_wind_snp: Some(1),
            annot: None,
            thin_annot: false,
            extract: None,
            maf: None,
            maf_pre: false,
            print_snps: None,
            keep: None,
            per_allele: false,
            pq_exp: Some(1.0),
            no_print_annot: false,
            chunk_size: 50,
            yes_really: true,
        };
        run(args).unwrap();

        let m_path = format!("{}1.l2.M", out);
        let m = std::fs::read_to_string(m_path).unwrap();
        assert!(m.trim().starts_with("0.25"));

        let m_5_50_path = format!("{}1.l2.M_5_50", out);
        let m_5_50 = std::fs::read_to_string(m_5_50_path).unwrap();
        assert!(m_5_50.trim().starts_with("0.25"));

        let ld_path = format!("{}1.l2.ldscore.gz", out);
        let content = read_gz(&ld_path);
        let header = content.lines().next().unwrap_or("");
        assert!(header.contains("L2_S1"));
    }
}
