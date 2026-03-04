use crate::bed::{Bed, ReadOptions};
use crate::la::{
    MatF, MatF32, mat_add_in_place, mat_copy_from, mat_slice, mat_slice_mut, matmul_tn_to,
    matmul_to,
};
/// LD score computation.
///
/// Memory: .bed files are read via the internal BED reader; peak RAM is
/// `window_snps × n_indiv × 4 bytes` (≈ 20 MB for 1000G, ≈ 4 GB for biobanks).
///
/// Algorithm: for each SNP i, `L2[i] = Σ_j r²_unbiased(i, j)` where j ranges
/// over the window, using `r²_unbiased = r² − (1−r²)/(n−2)` [Bulik-Sullivan 2015].
/// SNPs are processed in chunks so the inner A^T B product is a single matmul.
///
/// Partitioned LD scores (`--annot`): `L2_k[i] = Σ_j r²_unbiased(i,j) × annot[j,k]`
/// yields one column per annotation; `.M`/`.M_5_50` files are tab-separated.
use anyhow::{Context, Result};
use faer::{Accum, Par};
use std::collections::{HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write as IoWrite};
use std::time::Instant;
use tracing::trace;

use crate::cli::L2Args;
use crate::cts_annot;
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
fn load_individual_indices(keep_path: &str, fam_ids: &[(String, String)]) -> Result<Vec<isize>> {
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

/// Compute block_left per chromosome to avoid cross-chromosome LD windows.
/// Returns (block_left, any_full_chr_window).
fn get_block_lefts_by_chr(all_snps: &[BimRecord], mode: WindowMode) -> (Vec<usize>, bool) {
    let m = all_snps.len();
    let mut block_left = vec![0usize; m];
    let mut any_full_chr_window = false;

    let mut start = 0usize;
    while start < m {
        let chr = all_snps[start].chr;
        let mut end = start + 1;
        while end < m && all_snps[end].chr == chr {
            end += 1;
        }
        let len = end - start;
        if len == 0 {
            break;
        }

        match mode {
            WindowMode::Snp(half) => {
                if half >= len.saturating_sub(1) {
                    any_full_chr_window = true;
                }
                for i in 0..len {
                    block_left[start + i] = start + i.saturating_sub(half);
                }
            }
            WindowMode::Cm(max_cm) => {
                let coords: Vec<f64> = all_snps[start..end].iter().map(|s| s.cm).collect();
                let local = get_block_lefts_f64(&coords, max_cm);
                if local.last().copied().unwrap_or(0) == 0 {
                    any_full_chr_window = true;
                }
                for (i, left) in local.into_iter().enumerate() {
                    block_left[start + i] = start + left;
                }
            }
            WindowMode::Kb(max_kb) => {
                let coords: Vec<f64> = all_snps[start..end]
                    .iter()
                    .map(|s| s.bp as f64 / 1000.0)
                    .collect();
                let local = get_block_lefts_f64(&coords, max_kb);
                if local.last().copied().unwrap_or(0) == 0 {
                    any_full_chr_window = true;
                }
                for (i, left) in local.into_iter().enumerate() {
                    block_left[start + i] = start + left;
                }
            }
        }

        start = end;
    }

    (block_left, any_full_chr_window)
}

// ---------------------------------------------------------------------------
// Genotype normalisation helpers
// ---------------------------------------------------------------------------

/// Normalize a genotype column in-place: impute NaN with mean, centre, scale to unit variance.
///
/// Returns the minor allele frequency: `MAF = min(mean/2, 1 − mean/2)`.
/// Works on f32; caller casts to f64 into the pre-allocated `b_mat`.
fn normalize_col_f64(col: &mut [f64], n: usize) -> f64 {
    let (sum, count) = col.iter().fold((0f64, 0usize), |(s, c), &v| {
        if v.is_nan() { (s, c) } else { (s + v, c + 1) }
    });
    let avg = if count > 0 { sum / count as f64 } else { 0.0 };
    let freq = (avg / 2.0).clamp(0.0, 1.0);
    let maf = freq.min(1.0 - freq);

    for v in col.iter_mut() {
        if v.is_nan() {
            *v = avg;
        } else {
            *v -= avg;
        }
    }

    let var: f64 = col.iter().map(|&v| v * v).sum::<f64>() / n as f64;
    let std = var.sqrt();
    if std > 0.0 {
        let inv_std = 1.0 / std;
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
    annot: Option<&MatF>,
    iid_indices: Option<&[isize]>,
    pq_exp: Option<f64>,
    yes_really: bool,
) -> Result<(MatF, Vec<f64>)> {
    let m = all_snps.len();
    if m == 0 {
        let n_annot = annot.map(|a| a.ncols()).unwrap_or(1);
        return Ok((MatF::zeros(0, n_annot), vec![]));
    }

    // Compute block_left per chromosome (avoid cross-chromosome LD windows).
    let (block_left, any_full_chr_window) = get_block_lefts_by_chr(all_snps, mode);
    if !yes_really && any_full_chr_window {
        println!(
            "WARNING: LD window spans the entire chromosome. \
             Use --yes-really to silence this warning."
        );
    }

    let n = n_indiv as f64;
    let mut maf_per_snp = vec![0.0f64; m];
    let trace_snp = std::env::var("LDSC_TRACE_SNP").ok();
    let trace_ab_chunks = std::env::var("LDSC_TRACE_AB_CHUNKS").is_ok();
    let trace_idx = trace_snp
        .as_ref()
        .and_then(|snp| all_snps.iter().position(|rec| rec.snp == *snp));
    if let Some(ref snp) = trace_snp
        && trace_idx.is_none()
    {
        trace!("trace_snp '{}' not found in BIM list", snp);
    }
    let mut trace_diag = 0.0f64;
    let mut trace_bb = 0.0f64;
    let mut trace_ab = 0.0f64;

    let block_sizes: Vec<usize> = (0..m)
        .map(|i| {
            let dist = i - block_left[i];
            if dist == 0 {
                0
            } else {
                dist.div_ceil(chunk_c) * chunk_c
            }
        })
        .collect();
    let max_window_size = block_sizes.iter().copied().max().unwrap_or(0);

    // Ring buffer: moderate slack to minimize wrap copies (HPC memory OK).
    // Align to chunk size so wraps happen on chunk boundaries.
    let mut ring_size = (max_window_size + 8 * chunk_c).max(1);
    if chunk_c > 0 {
        ring_size = ring_size.div_ceil(chunk_c) * chunk_c;
    }
    let mut ring_buf = MatF::zeros(n_indiv, ring_size);
    let mut ring_next: usize = 0;
    let mut window: VecDeque<(usize, usize)> = VecDeque::new(); // (snp_idx, ring_slot)

    let mut b_mat = MatF::zeros(n_indiv, chunk_c);
    let mut col_buf_f64 = vec![0.0f64; n_indiv];
    let mut a_buf = MatF::zeros(n_indiv, max_window_size.max(1));
    let mut ab_buf = MatF::zeros(max_window_size.max(1), chunk_c);

    let mut pq_per_snp: Option<Vec<f64>> = pq_exp.map(|_| vec![1.0f64; m]);
    let mut pq_chunk: Vec<f64> = Vec::with_capacity(chunk_c);

    let mut t_read_bed = std::time::Duration::ZERO;
    let mut t_norm = std::time::Duration::ZERO;
    let mut t_bb = std::time::Duration::ZERO;
    let mut t_ab = std::time::Duration::ZERO;
    let mut t_r2u = std::time::Duration::ZERO;
    let mut t_r2u_bb = std::time::Duration::ZERO;
    let mut t_r2u_ab = std::time::Duration::ZERO;
    let mut t_annot = std::time::Duration::ZERO;
    let mut t_annot_bb = std::time::Duration::ZERO;
    let mut t_annot_ab_left = std::time::Duration::ZERO;
    let mut t_annot_ab_right = std::time::Duration::ZERO;
    let mut t_ab_wrap_copy = std::time::Duration::ZERO;
    let mut t_ab_contig = std::time::Duration::ZERO;
    let mut t_ab_wrap_matmul = std::time::Duration::ZERO;
    let mut ab_calls = 0usize;
    let mut ab_contig_calls = 0usize;
    let mut ab_wrap_calls = 0usize;
    let mut ab_sum_w = 0usize;
    let mut ab_sum_c = 0usize;
    let mut ab_max_w = 0usize;
    let mut ab_max_c = 0usize;
    if annot.is_none() {
        let mut l2_scalar = vec![0.0f64; m];

        for chunk_start in (0..m).step_by(chunk_c) {
            let chunk_end = (chunk_start + chunk_c).min(m);
            let c = chunk_end - chunk_start;

            let mut b = block_sizes[chunk_start];
            if b > chunk_start {
                b = chunk_start;
            }
            let a_left = chunk_start - b;
            while window
                .front()
                .map(|(idx, _)| *idx < a_left)
                .unwrap_or(false)
            {
                window.pop_front();
            }

            let bed_indices: Vec<isize> = all_snps[chunk_start..chunk_end]
                .iter()
                .map(|s| s.bed_idx as isize)
                .collect();
            let raw: MatF32 = {
                let mut ro = ReadOptions::builder();
                let mut builder = ro.sid_index(&bed_indices).f32();
                let t = Instant::now();
                let result = if let Some(iids) = iid_indices {
                    builder.iid_index(iids).read(bed)
                } else {
                    builder.read(bed)
                };
                t_read_bed += t.elapsed();
                result
                    .with_context(|| format!("reading BED chunk [{},{})", chunk_start, chunk_end))?
            };

            {
                let t = Instant::now();
                for j in 0..c {
                    for i in 0..n_indiv {
                        col_buf_f64[i] = raw[(i, j)] as f64;
                    }
                    let snp_maf = normalize_col_f64(&mut col_buf_f64, n_indiv);
                    maf_per_snp[chunk_start + j] = snp_maf;
                    if let (Some(exp), Some(ref mut pq_vec)) = (pq_exp, pq_per_snp.as_mut()) {
                        let pq = (snp_maf * (1.0 - snp_maf)).powf(exp);
                        pq_vec[chunk_start + j] = pq;
                    }
                    for i in 0..n_indiv {
                        b_mat[(i, j)] = col_buf_f64[i];
                    }
                }
                t_norm += t.elapsed();
            }
            let b_slice = mat_slice(b_mat.as_ref(), 0..n_indiv, 0..c);

            pq_chunk.clear();
            if let Some(ref pq_vec) = pq_per_snp {
                for j in 0..c {
                    pq_chunk.push(pq_vec[chunk_start + j]);
                }
            } else {
                pq_chunk.resize(c, 1.0);
            }

            // B×B
            let mut bb = MatF::zeros(c, c);
            let t = Instant::now();
            matmul_tn_to(
                bb.as_mut(),
                b_slice,
                b_slice,
                1.0,
                Accum::Replace,
                Par::rayon(0),
            );
            t_bb += t.elapsed();
            let t = Instant::now();
            for j in 0..c {
                let j_g = chunk_start + j;
                let pq_j = pq_chunk[j];
                l2_scalar[j_g] += pq_j; // diagonal: r(j,j) = 1
                if trace_idx == Some(j_g) {
                    trace_diag += pq_j;
                }
                for k in 0..j {
                    let k_g = chunk_start + k;
                    let pq_k = pq_chunk[k];
                    let r2u = r2_unbiased(bb[(k, j)] / n, n_indiv);
                    l2_scalar[j_g] += r2u * pq_k;
                    l2_scalar[k_g] += r2u * pq_j;
                    if trace_idx == Some(j_g) {
                        trace_bb += r2u * pq_k;
                    } else if trace_idx == Some(k_g) {
                        trace_bb += r2u * pq_j;
                    }
                }
            }
            let elapsed = t.elapsed();
            t_r2u += elapsed;
            t_r2u_bb += elapsed;

            // A×B
            if !window.is_empty() {
                let w = window.len();
                let (first_slot, last_slot) = window
                    .front()
                    .and_then(|(_, f)| window.back().map(|(_, l)| (*f, *l)))
                    .unwrap_or((0, 0));
                let contiguous = first_slot <= last_slot && last_slot - first_slot + 1 == w;

                if !contiguous {
                    let t_copy = Instant::now();
                    for (wi, (_, slot)) in window.iter().enumerate() {
                        for i in 0..n_indiv {
                            a_buf[(i, wi)] = ring_buf[(i, *slot)];
                        }
                    }
                    t_ab_wrap_copy += t_copy.elapsed();
                }
                let t = Instant::now();
                if contiguous {
                    ab_contig_calls += 1;
                    let a_view =
                        mat_slice(ring_buf.as_ref(), 0..n_indiv, first_slot..(first_slot + w));
                    matmul_tn_to(
                        mat_slice_mut(ab_buf.as_mut(), 0..w, 0..c),
                        a_view,
                        b_slice,
                        1.0,
                        Accum::Replace,
                        Par::rayon(0),
                    );
                } else {
                    ab_wrap_calls += 1;
                    let a_view = mat_slice(a_buf.as_ref(), 0..n_indiv, 0..w);
                    matmul_tn_to(
                        mat_slice_mut(ab_buf.as_mut(), 0..w, 0..c),
                        a_view,
                        b_slice,
                        1.0,
                        Accum::Replace,
                        Par::rayon(0),
                    );
                }
                let ab_elapsed = t.elapsed();
                if contiguous {
                    t_ab_contig += ab_elapsed;
                } else {
                    t_ab_wrap_matmul += ab_elapsed;
                }
                t_ab += ab_elapsed;
                ab_calls += 1;
                ab_sum_w += w;
                ab_sum_c += c;
                ab_max_w = ab_max_w.max(w);
                ab_max_c = ab_max_c.max(c);
                let t = Instant::now();
                let ab_view = mat_slice(ab_buf.as_ref(), 0..w, 0..c);
                let pq_vec_opt = pq_per_snp.as_ref();
                for (wi, (k_g, _)) in window.iter().enumerate() {
                    let pq_k = if let Some(pq_vec) = pq_vec_opt {
                        pq_vec[*k_g]
                    } else {
                        1.0
                    };
                    for j in 0..c {
                        let j_g = chunk_start + j;
                        let pq_j = pq_chunk[j];
                        let r2u = r2_unbiased(ab_view[(wi, j)] / n, n_indiv);
                        l2_scalar[j_g] += r2u * pq_k;
                        l2_scalar[*k_g] += r2u * pq_j;
                        if trace_idx == Some(j_g) {
                            trace_ab += r2u * pq_k;
                        } else if trace_idx == Some(*k_g) {
                            trace_ab += r2u * pq_j;
                        }
                    }
                }
                let r2u_elapsed = t.elapsed();
                t_r2u += r2u_elapsed;
                t_r2u_ab += r2u_elapsed;

                if trace_ab_chunks {
                    trace!(
                        "l2 ab_chunk chunk_start={} chunk_end={} c={} w={} contig={} did_ab=1 chunk_all_zero=0 window_all_zero=0 ab_ms={:.3} r2u_ab_ms={:.3} annot_left_ms=0.000 annot_right_ms=0.000",
                        chunk_start,
                        chunk_end,
                        c,
                        w,
                        contiguous as u8,
                        ab_elapsed.as_secs_f64() * 1000.0,
                        r2u_elapsed.as_secs_f64() * 1000.0
                    );
                }
            }

            for j in 0..c {
                let slot = ring_next % ring_size;
                for i in 0..n_indiv {
                    ring_buf[(i, slot)] = b_mat[(i, j)];
                }
                window.push_back((chunk_start + j, slot));
                ring_next += 1;
            }
        }

        trace!(
            "l2 timing: bed_read={:?} norm={:?} bb_dot={:?} ab_dot={:?} r2u_bb={:?} r2u_ab={:?} ab_calls={} ab_contig={} ab_wrap={} ab_wrap_copy={:?} ab_contig_matmul={:?} ab_wrap_matmul={:?} ab_avg_w={:.2} ab_avg_c={:.2} ab_max_w={} ab_max_c={}",
            t_read_bed,
            t_norm,
            t_bb,
            t_ab,
            t_r2u_bb,
            t_r2u_ab,
            ab_calls,
            ab_contig_calls,
            ab_wrap_calls,
            t_ab_wrap_copy,
            t_ab_contig,
            t_ab_wrap_matmul,
            if ab_calls > 0 {
                ab_sum_w as f64 / ab_calls as f64
            } else {
                0.0
            },
            if ab_calls > 0 {
                ab_sum_c as f64 / ab_calls as f64
            } else {
                0.0
            },
            ab_max_w,
            ab_max_c
        );
        if let (Some(snp), Some(idx)) = (trace_snp.as_ref(), trace_idx) {
            let total = trace_diag + trace_bb + trace_ab;
            let trace_output = l2_scalar[idx];
            trace!(
                "trace_snp={} idx={} l2_total={:.12} diag={:.12} bb={:.12} ab={:.12} output={:.12}",
                snp,
                idx,
                total,
                trace_diag,
                trace_bb + trace_diag,
                trace_ab,
                trace_output
            );
        }
        let mut l2 = MatF::zeros(m, 1);
        for (i, v) in l2_scalar.into_iter().enumerate() {
            l2[(i, 0)] = v;
        }
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
    let mut l2 = MatF::zeros(m, n_annot);
    let mut r2u_bb = MatF::zeros(chunk_c, chunk_c);
    let mut r2u_ab = MatF::zeros(max_window_size.max(1), chunk_c);

    for chunk_start in (0..m).step_by(chunk_c) {
        let chunk_end = (chunk_start + chunk_c).min(m);
        let c = chunk_end - chunk_start;

        let mut b = block_sizes[chunk_start];
        if b > chunk_start {
            b = chunk_start;
        }
        let a_left = chunk_start - b;
        while window
            .front()
            .map(|(idx, _)| *idx < a_left)
            .unwrap_or(false)
        {
            window.pop_front();
        }

        let bed_indices: Vec<isize> = all_snps[chunk_start..chunk_end]
            .iter()
            .map(|s| s.bed_idx as isize)
            .collect();
        let raw: MatF32 = {
            let mut ro = ReadOptions::builder();
            let mut builder = ro.sid_index(&bed_indices).f32();
            let t = Instant::now();
            let result = if let Some(iids) = iid_indices {
                builder.iid_index(iids).read(bed)
            } else {
                builder.read(bed)
            };
            t_read_bed += t.elapsed();
            result.with_context(|| format!("reading BED chunk [{},{})", chunk_start, chunk_end))?
        };

        {
            let t = Instant::now();
            for j in 0..c {
                for i in 0..n_indiv {
                    col_buf_f64[i] = raw[(i, j)] as f64;
                }
                let snp_maf = normalize_col_f64(&mut col_buf_f64, n_indiv);
                maf_per_snp[chunk_start + j] = snp_maf;
                if let (Some(exp), Some(ref mut pq_vec)) = (pq_exp, pq_per_snp.as_mut()) {
                    let pq = (snp_maf * (1.0 - snp_maf)).powf(exp);
                    pq_vec[chunk_start + j] = pq;
                }
                for i in 0..n_indiv {
                    b_mat[(i, j)] = col_buf_f64[i];
                }
            }
            t_norm += t.elapsed();
        }
        let b_slice = mat_slice(b_mat.as_ref(), 0..n_indiv, 0..c);

        let annot_chunk = mat_slice(annot.as_ref(), chunk_start..chunk_end, 0..n_annot);
        let (annot_chunk_eff, chunk_all_zero) = if let Some(ref pq_vec) = pq_per_snp {
            let mut eff = MatF::zeros(c, n_annot);
            mat_copy_from(eff.as_mut(), annot_chunk);
            for j in 0..c {
                let pq = pq_vec[chunk_start + j];
                for k in 0..n_annot {
                    eff[(j, k)] *= pq;
                }
            }
            let mut all_zero = true;
            for i in 0..c {
                for k in 0..n_annot {
                    if eff[(i, k)] != 0.0 {
                        all_zero = false;
                        break;
                    }
                }
                if !all_zero {
                    break;
                }
            }
            (Some(eff), all_zero)
        } else {
            let mut all_zero = true;
            for i in 0..c {
                for k in 0..n_annot {
                    if annot_chunk[(i, k)] != 0.0 {
                        all_zero = false;
                        break;
                    }
                }
                if !all_zero {
                    break;
                }
            }
            (None, all_zero)
        };

        // B×B
        if !chunk_all_zero {
            let mut bb = MatF::zeros(c, c);
            let t = Instant::now();
            matmul_tn_to(
                bb.as_mut(),
                b_slice,
                b_slice,
                1.0,
                Accum::Replace,
                Par::rayon(0),
            );
            t_bb += t.elapsed();
            let t = Instant::now();
            for j in 0..c {
                for k in 0..c {
                    r2u_bb[(j, k)] = 0.0;
                }
            }
            for j in 0..c {
                r2u_bb[(j, j)] = 1.0; // diagonal: r(j,j)=1 → r2u=1
                for k in 0..j {
                    let r2u = r2_unbiased(bb[(k, j)] / n, n_indiv);
                    r2u_bb[(j, k)] = r2u;
                    r2u_bb[(k, j)] = r2u;
                }
            }
            let elapsed = t.elapsed();
            t_r2u += elapsed;
            t_r2u_bb += elapsed;
            let t = Instant::now();
            let r2u_bb_view = mat_slice(r2u_bb.as_ref(), 0..c, 0..c);
            let mut contrib_bb = MatF::zeros(c, n_annot);
            let contrib_bb = if let Some(ref eff) = annot_chunk_eff {
                matmul_to(
                    contrib_bb.as_mut(),
                    r2u_bb_view,
                    eff.as_ref(),
                    1.0,
                    Accum::Replace,
                    Par::rayon(0),
                );
                contrib_bb
            } else {
                matmul_to(
                    contrib_bb.as_mut(),
                    r2u_bb_view,
                    annot_chunk,
                    1.0,
                    Accum::Replace,
                    Par::rayon(0),
                );
                contrib_bb
            };
            mat_add_in_place(
                mat_slice_mut(l2.as_mut(), chunk_start..chunk_end, 0..n_annot),
                contrib_bb.as_ref(),
            );
            let elapsed = t.elapsed();
            t_annot += elapsed;
            t_annot_bb += elapsed;
        }

        // A×B
        if !window.is_empty() {
            let w = window.len();
            let a_left_idx = window.front().map(|(idx, _)| *idx).unwrap_or(0);
            let (first_slot, last_slot) = window
                .front()
                .and_then(|(_, f)| window.back().map(|(_, l)| (*f, *l)))
                .unwrap_or((0, 0));
            let contiguous = first_slot <= last_slot && last_slot - first_slot + 1 == w;

            if !contiguous {
                let t_copy = Instant::now();
                for (wi, (_, slot)) in window.iter().enumerate() {
                    for i in 0..n_indiv {
                        a_buf[(i, wi)] = ring_buf[(i, *slot)];
                    }
                }
                t_ab_wrap_copy += t_copy.elapsed();
            }
            let annot_window = mat_slice(annot.as_ref(), a_left_idx..chunk_start, 0..n_annot);
            let (annot_window_eff, window_all_zero) = if let Some(ref pq_vec) = pq_per_snp {
                let mut eff = MatF::zeros(w, n_annot);
                mat_copy_from(eff.as_mut(), annot_window);
                for (wi, (w_g, _)) in window.iter().enumerate() {
                    let pq = pq_vec[*w_g];
                    for k in 0..n_annot {
                        eff[(wi, k)] *= pq;
                    }
                }
                let mut all_zero = true;
                for i in 0..w {
                    for k in 0..n_annot {
                        if eff[(i, k)] != 0.0 {
                            all_zero = false;
                            break;
                        }
                    }
                    if !all_zero {
                        break;
                    }
                }
                (Some(eff), all_zero)
            } else {
                let mut all_zero = true;
                for i in 0..w {
                    for k in 0..n_annot {
                        if annot_window[(i, k)] != 0.0 {
                            all_zero = false;
                            break;
                        }
                    }
                    if !all_zero {
                        break;
                    }
                }
                (None, all_zero)
            };

            let mut ab_elapsed = std::time::Duration::ZERO;
            let mut r2u_elapsed = std::time::Duration::ZERO;
            let mut annot_left_elapsed = std::time::Duration::ZERO;
            let mut annot_right_elapsed = std::time::Duration::ZERO;
            let mut did_ab = false;

            if !(chunk_all_zero && window_all_zero) {
                let t = Instant::now();
                if contiguous {
                    ab_contig_calls += 1;
                    let a_view =
                        mat_slice(ring_buf.as_ref(), 0..n_indiv, first_slot..(first_slot + w));
                    matmul_tn_to(
                        mat_slice_mut(ab_buf.as_mut(), 0..w, 0..c),
                        a_view,
                        b_slice,
                        1.0,
                        Accum::Replace,
                        Par::rayon(0),
                    );
                } else {
                    ab_wrap_calls += 1;
                    let a_view = mat_slice(a_buf.as_ref(), 0..n_indiv, 0..w);
                    matmul_tn_to(
                        mat_slice_mut(ab_buf.as_mut(), 0..w, 0..c),
                        a_view,
                        b_slice,
                        1.0,
                        Accum::Replace,
                        Par::rayon(0),
                    );
                }
                ab_elapsed = t.elapsed();
                if contiguous {
                    t_ab_contig += ab_elapsed;
                } else {
                    t_ab_wrap_matmul += ab_elapsed;
                }
                t_ab += ab_elapsed;
                ab_calls += 1;
                ab_sum_w += w;
                ab_sum_c += c;
                ab_max_w = ab_max_w.max(w);
                ab_max_c = ab_max_c.max(c);
                did_ab = true;

                let t = Instant::now();
                let ab_view = mat_slice(ab_buf.as_ref(), 0..w, 0..c);
                for wi in 0..w {
                    for j in 0..c {
                        r2u_ab[(wi, j)] = r2_unbiased(ab_view[(wi, j)] / n, n_indiv);
                    }
                }
                r2u_elapsed = t.elapsed();
                t_r2u += r2u_elapsed;
                t_r2u_ab += r2u_elapsed;

                if !chunk_all_zero {
                    let t = Instant::now();
                    let r2u_ab_view = mat_slice(r2u_ab.as_ref(), 0..w, 0..c);
                    let mut contrib_left = MatF::zeros(w, n_annot);
                    let contrib_left = if let Some(ref eff) = annot_chunk_eff {
                        matmul_to(
                            contrib_left.as_mut(),
                            r2u_ab_view,
                            eff.as_ref(),
                            1.0,
                            Accum::Replace,
                            Par::rayon(0),
                        );
                        contrib_left
                    } else {
                        matmul_to(
                            contrib_left.as_mut(),
                            r2u_ab_view,
                            annot_chunk,
                            1.0,
                            Accum::Replace,
                            Par::rayon(0),
                        );
                        contrib_left
                    };
                    mat_add_in_place(
                        mat_slice_mut(l2.as_mut(), a_left_idx..chunk_start, 0..n_annot),
                        contrib_left.as_ref(),
                    );
                    annot_left_elapsed = t.elapsed();
                    t_annot += annot_left_elapsed;
                    t_annot_ab_left += annot_left_elapsed;
                }

                if !window_all_zero {
                    let t = Instant::now();
                    let r2u_ab_view = mat_slice(r2u_ab.as_ref(), 0..w, 0..c);
                    let mut contrib_right = MatF::zeros(c, n_annot);
                    let contrib_right = if let Some(ref eff) = annot_window_eff {
                        matmul_tn_to(
                            contrib_right.as_mut(),
                            r2u_ab_view,
                            eff.as_ref(),
                            1.0,
                            Accum::Replace,
                            Par::rayon(0),
                        );
                        contrib_right
                    } else {
                        matmul_tn_to(
                            contrib_right.as_mut(),
                            r2u_ab_view,
                            annot_window,
                            1.0,
                            Accum::Replace,
                            Par::rayon(0),
                        );
                        contrib_right
                    };
                    mat_add_in_place(
                        mat_slice_mut(l2.as_mut(), chunk_start..chunk_end, 0..n_annot),
                        contrib_right.as_ref(),
                    );
                    annot_right_elapsed = t.elapsed();
                    t_annot += annot_right_elapsed;
                    t_annot_ab_right += annot_right_elapsed;
                }
            }

            if trace_ab_chunks {
                trace!(
                    "l2 ab_chunk chunk_start={} chunk_end={} c={} w={} contig={} did_ab={} chunk_all_zero={} window_all_zero={} ab_ms={:.3} r2u_ab_ms={:.3} annot_left_ms={:.3} annot_right_ms={:.3}",
                    chunk_start,
                    chunk_end,
                    c,
                    w,
                    contiguous as u8,
                    did_ab as u8,
                    chunk_all_zero as u8,
                    window_all_zero as u8,
                    ab_elapsed.as_secs_f64() * 1000.0,
                    r2u_elapsed.as_secs_f64() * 1000.0,
                    annot_left_elapsed.as_secs_f64() * 1000.0,
                    annot_right_elapsed.as_secs_f64() * 1000.0
                );
            }
        }

        for j in 0..c {
            let slot = ring_next % ring_size;
            for i in 0..n_indiv {
                ring_buf[(i, slot)] = b_mat[(i, j)];
            }
            window.push_back((chunk_start + j, slot));
            ring_next += 1;
        }
    }

    trace!(
        "l2 timing: bed_read={:?} norm={:?} bb_dot={:?} ab_dot={:?} r2u_bb={:?} r2u_ab={:?} annot_total={:?} annot_bb={:?} annot_ab_left={:?} annot_ab_right={:?} ab_calls={} ab_contig={} ab_wrap={} ab_wrap_copy={:?} ab_contig_matmul={:?} ab_wrap_matmul={:?} ab_avg_w={:.2} ab_avg_c={:.2} ab_max_w={} ab_max_c={}",
        t_read_bed,
        t_norm,
        t_bb,
        t_ab,
        t_r2u_bb,
        t_r2u_ab,
        t_annot,
        t_annot_bb,
        t_annot_ab_left,
        t_annot_ab_right,
        ab_calls,
        ab_contig_calls,
        ab_wrap_calls,
        t_ab_wrap_copy,
        t_ab_contig,
        t_ab_wrap_matmul,
        if ab_calls > 0 {
            ab_sum_w as f64 / ab_calls as f64
        } else {
            0.0
        },
        if ab_calls > 0 {
            ab_sum_c as f64 / ab_calls as f64
        } else {
            0.0
        },
        ab_max_w,
        ab_max_c
    );
    Ok((l2, maf_per_snp))
}

fn compute_snp_stats(
    all_snps: &[BimRecord],
    bed: &mut Bed,
    n_indiv: usize,
    chunk_c: usize,
    iid_indices: Option<&[isize]>,
) -> Result<(Vec<f64>, Vec<bool>)> {
    let m = all_snps.len();
    if m == 0 {
        return Ok((vec![], vec![]));
    }

    #[derive(Clone, Copy)]
    struct ByteStats {
        sum: u8,
        count: u8,
        het: u8,
    }

    let lut: [ByteStats; 256] = std::array::from_fn(|b| {
        let byte = b as u8;
        let mut sum = 0u8;
        let mut count = 0u8;
        let mut het = 0u8;
        for k in 0..4 {
            let bits = (byte >> (2 * k)) & 0b11;
            match bits {
                0 => {
                    sum += 2;
                    count += 1;
                }
                1 => {}
                2 => {
                    sum += 1;
                    count += 1;
                    het += 1;
                }
                3 => {
                    count += 1;
                }
                _ => {}
            }
        }
        ByteStats { sum, count, het }
    });

    let bed_indices: Vec<usize> = all_snps.iter().map(|s| s.bed_idx).collect();
    let keep_iids: Option<Vec<usize>> =
        iid_indices.map(|idxs| idxs.iter().map(|&i| i as usize).collect());
    let keep_locs: Option<Vec<(usize, u8)>> = keep_iids.as_ref().map(|idxs| {
        idxs.iter()
            .map(|&iid| (iid / 4, ((iid % 4) * 2) as u8))
            .collect()
    });
    if let Some(ref idxs) = keep_iids {
        debug_assert_eq!(idxs.len(), n_indiv);
    }

    let mut maf_per_snp = vec![0.0f64; m];
    let mut het_miss_ok = vec![true; m];

    let full_bytes = n_indiv / 4;
    let rem = n_indiv % 4;
    let bytes_per_snp = bed.bytes_per_snp();
    let mut file = bed.open_reader()?;
    let mut buf = vec![0u8; bytes_per_snp];
    let mut block_buf: Vec<u8> = Vec::new();
    let chunk_c = chunk_c.max(1);

    let compute_stats = |bytes: &[u8],
                         keep_locs: Option<&Vec<(usize, u8)>>,
                         lut: &[ByteStats; 256],
                         full_bytes: usize,
                         rem: usize|
     -> (u32, u32, u32) {
        let mut sum = 0u32;
        let mut count = 0u32;
        let mut het = 0u32;
        if let Some(locs) = keep_locs {
            for &(byte_idx, shift) in locs {
                let bits = (bytes[byte_idx] >> shift) & 0b11;
                match bits {
                    0 => {
                        sum += 2;
                        count += 1;
                    }
                    1 => {}
                    2 => {
                        sum += 1;
                        count += 1;
                        het += 1;
                    }
                    3 => {
                        count += 1;
                    }
                    _ => {}
                }
            }
        } else {
            for &byte in &bytes[..full_bytes] {
                let stats = lut[byte as usize];
                sum += stats.sum as u32;
                count += stats.count as u32;
                het += stats.het as u32;
            }
            if rem > 0 {
                let byte = bytes[full_bytes];
                for k in 0..rem {
                    let bits = (byte >> (2 * k)) & 0b11;
                    match bits {
                        0 => {
                            sum += 2;
                            count += 1;
                        }
                        1 => {}
                        2 => {
                            sum += 1;
                            count += 1;
                            het += 1;
                        }
                        3 => {
                            count += 1;
                        }
                        _ => {}
                    }
                }
            }
        }
        (sum, count, het)
    };

    for chunk_start in (0..m).step_by(chunk_c) {
        let chunk_end = (chunk_start + chunk_c).min(m);
        let c = chunk_end - chunk_start;
        if c == 0 {
            continue;
        }

        let chunk_bed = &bed_indices[chunk_start..chunk_end];
        let contiguous = chunk_bed
            .iter()
            .enumerate()
            .skip(1)
            .all(|(i, &v)| v == chunk_bed[0] + i);

        if contiguous {
            let start_sid = chunk_bed[0];
            let total_bytes = c * bytes_per_snp;
            block_buf.resize(total_bytes, 0u8);
            let offset = bed.snp_offset(start_sid);
            file.seek(SeekFrom::Start(offset))
                .with_context(|| format!("seeking BED offset {}", offset))?;
            file.read_exact(&mut block_buf)
                .with_context(|| format!("reading BED block [{},{})", chunk_start, chunk_end))?;
            for j in 0..c {
                let start = j * bytes_per_snp;
                let end = start + bytes_per_snp;
                let bytes = &block_buf[start..end];
                let (sum, count, het) =
                    compute_stats(bytes, keep_locs.as_ref(), &lut, full_bytes, rem);
                let count_usize = count as usize;
                let missing = n_indiv.saturating_sub(count_usize);
                let het_miss = het as usize + missing;
                let freq = if count > 0 {
                    sum as f64 / (2.0 * count as f64)
                } else {
                    0.0
                };
                let maf = freq.min(1.0 - freq);
                maf_per_snp[chunk_start + j] = maf;
                het_miss_ok[chunk_start + j] = het_miss < n_indiv;
            }
        } else {
            for (j, &sid) in chunk_bed.iter().enumerate() {
                let offset = bed.snp_offset(sid);
                file.seek(SeekFrom::Start(offset))
                    .with_context(|| format!("seeking BED offset {}", offset))?;
                file.read_exact(&mut buf)
                    .with_context(|| format!("reading BED SNP {}", sid))?;
                let (sum, count, het) =
                    compute_stats(&buf, keep_locs.as_ref(), &lut, full_bytes, rem);
                let count_usize = count as usize;
                let missing = n_indiv.saturating_sub(count_usize);
                let het_miss = het as usize + missing;
                let freq = if count > 0 {
                    sum as f64 / (2.0 * count as f64)
                } else {
                    0.0
                };
                let maf = freq.min(1.0 - freq);
                maf_per_snp[chunk_start + j] = maf;
                het_miss_ok[chunk_start + j] = het_miss < n_indiv;
            }
        }
    }

    Ok((maf_per_snp, het_miss_ok))
}

/// Unbiased r² estimator: r² − (1−r²)/(n−2) [Bulik-Sullivan 2015].
#[inline]
fn r2_unbiased(r: f64, n: usize) -> f64 {
    let sq = r * r;
    let denom = if n > 2 { n as f64 - 2.0 } else { n as f64 };
    sq - (1.0 - sq) / denom
}

// ---------------------------------------------------------------------------
// Output writer
// ---------------------------------------------------------------------------

/// Write per-SNP LD scores to a gzip TSV.  `col_names` is `["L2"]` or `["{ANNOT}L2", …]`.
fn write_ldscore_refs(
    path: &str,
    snps: &[&BimRecord],
    l2: &MatF,
    col_names: &[String],
) -> Result<()> {
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::BufWriter;

    let file = File::create(path).with_context(|| format!("creating output '{}'", path))?;
    let gz = GzEncoder::new(file, Compression::default());
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
fn write_annot_matrix(
    path: &str,
    snps: &[BimRecord],
    annot: &MatF,
    col_names: &[String],
) -> Result<()> {
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::BufWriter;

    let file = File::create(path).with_context(|| format!("creating output '{}'", path))?;
    let gz = GzEncoder::new(file, Compression::default());
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

/// Load SNP IDs from a file with exactly one column (Python --print-snps behavior).
fn load_print_snps(path: &str) -> Result<HashSet<String>> {
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

pub fn run(args: L2Args) -> Result<()> {
    if args.per_allele && args.pq_exp.is_some() {
        anyhow::bail!(
            "Cannot set both --per-allele and --pq-exp (--per-allele is equivalent to --pq-exp 1)."
        );
    }
    if args.annot.is_some() && args.cts_bin.is_some() {
        anyhow::bail!("--annot and --cts-bin are currently incompatible.");
    }
    if args.cts_bin.is_some() && args.extract.is_some() {
        anyhow::bail!("--cts-bin and --extract are currently incompatible.");
    }
    if (args.cts_bin.is_some()) != (args.cts_breaks.is_some()) {
        anyhow::bail!("Must set both or neither of --cts-bin and --cts-breaks.");
    }

    let pq_exp = if args.per_allele {
        Some(1.0)
    } else {
        args.pq_exp
    };
    let mut wind_flags = 0u8;
    if args.ld_wind_cm.is_some() {
        wind_flags += 1;
    }
    if args.ld_wind_kb.is_some() {
        wind_flags += 1;
    }
    if args.ld_wind_snp.is_some() {
        wind_flags += 1;
    }
    if wind_flags > 1 {
        anyhow::bail!("Must specify exactly one --ld-wind option");
    }
    let mode = if let Some(kb) = args.ld_wind_kb {
        WindowMode::Kb(kb)
    } else if let Some(snp) = args.ld_wind_snp {
        WindowMode::Snp(snp)
    } else {
        WindowMode::Cm(args.ld_wind_cm.unwrap_or(1.0))
    };

    let bim_path = format!("{}.bim", args.bfile);
    let fam_path = format!("{}.fam", args.bfile);
    let bed_path = format!("{}.bed", args.bfile);

    let t_bim = Instant::now();
    let all_snps_raw =
        parse_bim(&bim_path).with_context(|| format!("parsing BIM '{}'", bim_path))?;
    trace!("l2 timing: parse_bim {:?}", t_bim.elapsed());

    let t_fam = Instant::now();
    let n_indiv = count_fam(&fam_path).with_context(|| format!("counting FAM '{}'", fam_path))?;
    trace!("l2 timing: count_fam {:?}", t_fam.elapsed());

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
            load_print_snps(ps_path)
                .with_context(|| format!("loading --print-snps file '{}'", ps_path))?,
        )
    } else {
        None
    };

    if args.annot.is_some() && args.extract.is_some() {
        println!(
            "WARNING: --annot with --extract is not supported. \
             Ensure your annot files match the extracted SNP set."
        );
    }

    let t_annot = Instant::now();
    let mut cts_bin_active = false;
    let mut annot_result: Option<(MatF, Vec<String>)> = if let Some(ref prefix) = args.annot {
        let explicit = prefix.ends_with(".annot")
            || prefix.ends_with(".annot.gz")
            || prefix.ends_with(".annot.bz2");
        if explicit {
            let path = parse::resolve_annot_path(prefix)?;
            let (mat, names) = parse::read_annot_path(&path, args.thin_annot)?;
            anyhow::ensure!(
                mat.nrows() == all_snps.len(),
                "Annotation file has {} rows but BIM has {} SNPs — they must match exactly",
                mat.nrows(),
                all_snps.len()
            );
            println!(
                "Read {} annotations for {} SNPs from '{}'",
                names.len(),
                mat.nrows(),
                path
            );
            Some((mat, names))
        } else if let Ok(path) = parse::resolve_annot_path(prefix) {
            let (mat, names) = parse::read_annot_path(&path, args.thin_annot)?;
            anyhow::ensure!(
                mat.nrows() == all_snps.len(),
                "Annotation file has {} rows but BIM has {} SNPs — they must match exactly",
                mat.nrows(),
                all_snps.len()
            );
            println!(
                "Read {} annotations for {} SNPs from '{}'",
                names.len(),
                mat.nrows(),
                path
            );
            Some((mat, names))
        } else {
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
            let mut mats: Vec<MatF> = Vec::new();
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

            let total_rows: usize = mats.iter().map(|m| m.nrows()).sum();
            let n_cols = mats.first().map(|m| m.ncols()).unwrap_or(0);
            let mut combined = MatF::zeros(total_rows, n_cols);
            let mut row_offset = 0usize;
            for mat in &mats {
                let rows = mat.nrows();
                for i in 0..rows {
                    for j in 0..n_cols {
                        combined[(row_offset + i, j)] = mat[(i, j)];
                    }
                }
                row_offset += rows;
            }
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
        }
    } else if let Some(ref cts_bin) = args.cts_bin {
        let breaks = args
            .cts_breaks
            .as_deref()
            .context("--cts-breaks required with --cts-bin")?;
        let snp_ids: Vec<String> = all_snps.iter().map(|s| s.snp.clone()).collect();
        let (mat, names) =
            cts_annot::build_cts_matrix(&snp_ids, cts_bin, breaks, args.cts_names.as_deref())?;
        println!(
            "Read {} annotations for {} SNPs from --cts-bin",
            names.len(),
            mat.nrows()
        );
        cts_bin_active = true;
        Some((mat, names))
    } else {
        None
    };
    trace!("l2 timing: read_annot {:?}", t_annot.elapsed());

    // --keep: subset individuals for LD computation.
    let iid_indices: Option<Vec<isize>> = if let Some(ref keep_path) = args.keep {
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
    let maf_pre = args.maf_pre;

    // Pre-filter SNPs with all-het/missing genotypes (Python behavior).
    let mut maf_prefilter: Option<Vec<f64>> = None;
    if maf_pre {
        let mut bed = Bed::new(bed_path.as_str()).context("opening BED file for SNP prefilter")?;
        let t_stats = Instant::now();
        let (maf_all, het_miss_ok) = compute_snp_stats(
            &all_snps,
            &mut bed,
            n_indiv_actual,
            args.chunk_size,
            iid_indices.as_deref(),
        )
        .context("computing SNP prefilter stats")?;
        trace!("l2 timing: maf_prefilter {:?}", t_stats.elapsed());

        let thr = args.maf.unwrap_or(0.0);
        let mut keep_mask: Vec<bool> = Vec::with_capacity(all_snps.len());
        let mut kept_snps: Vec<BimRecord> = Vec::new();
        for (i, snp) in all_snps.iter().enumerate() {
            let mut keep = het_miss_ok[i];
            keep &= maf_all[i] > thr;
            keep_mask.push(keep);
            if keep {
                kept_snps.push(snp.clone());
            }
        }

        if kept_snps.is_empty() {
            anyhow::bail!("SNP prefilter removed all SNPs");
        }

        if kept_snps.len() != all_snps.len() {
            if let Some((annot, names)) = annot_result.take() {
                let rows: Vec<usize> = keep_mask
                    .iter()
                    .enumerate()
                    .filter_map(|(i, &k)| if k { Some(i) } else { None })
                    .collect();
                let mut filtered = MatF::zeros(rows.len(), annot.ncols());
                for (ri, &src) in rows.iter().enumerate() {
                    for j in 0..annot.ncols() {
                        filtered[(ri, j)] = annot[(src, j)];
                    }
                }
                annot_result = Some((filtered, names));
            }
            all_snps = kept_snps;
        }

        let maf_kept: Vec<f64> = keep_mask
            .iter()
            .zip(maf_all.iter())
            .filter_map(|(k, v)| if *k { Some(*v) } else { None })
            .collect();
        if args.maf.is_some() {
            println!(
                "--maf-pre: kept {} / {} SNPs (MAF > {})",
                all_snps.len(),
                keep_mask.len(),
                thr
            );
        }
        maf_prefilter = Some(maf_kept);
    }

    let mut pq_scaled_annot = false;
    if let (Some(exp), Some(maf_vals)) = (pq_exp, maf_prefilter.as_ref())
        && let Some((annot, names)) = annot_result.take()
    {
        anyhow::ensure!(
            maf_vals.len() == annot.nrows(),
            "MAF length {} does not match annot rows {}",
            maf_vals.len(),
            annot.nrows()
        );
        let mut scaled = annot;
        for (i, maf) in maf_vals.iter().enumerate() {
            let pq = (maf * (1.0 - maf)).powf(exp);
            for j in 0..scaled.ncols() {
                scaled[(i, j)] *= pq;
            }
        }
        annot_result = Some((scaled, names));
        pq_scaled_annot = true;
    }

    let pq_exp_for_compute = if pq_scaled_annot { None } else { pq_exp };

    let annot_ref: Option<&MatF> = annot_result.as_ref().map(|(m, _)| m);

    let chunk_c = args.chunk_size;

    let mut bed = Bed::new(bed_path.as_str()).context("opening BED file")?;
    let t_l2 = Instant::now();
    let (l2, maf_per_snp) = compute_ldscore_global(
        &all_snps,
        &mut bed,
        n_indiv_actual,
        mode,
        chunk_c,
        annot_ref,
        iid_indices.as_deref(),
        pq_exp_for_compute,
        args.yes_really,
    )
    .context("computing LD scores")?;
    trace!("l2 timing: compute_ldscore {:?}", t_l2.elapsed());

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
    let pq_for_m = if pq_scaled_annot {
        None
    } else {
        pq_per_snp.as_ref()
    };

    if cts_bin_active
        && !args.no_print_annot
        && let Some((annot, _)) = annot_result.as_ref()
    {
        let annot_path = format!("{}.annot.gz", args.out);
        write_annot_matrix(&annot_path, &all_snps, annot, &col_names)
            .with_context(|| format!("writing annot matrix '{}'", annot_path))?;
        println!(
            "Writing annot matrix produced by --cts-bin to {}",
            annot_path
        );
    }

    let mut chrs: Vec<u8> = all_snps
        .iter()
        .map(|s| s.chr)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    chrs.sort();

    let out_positions: Vec<usize> = all_snps
        .iter()
        .filter(|s| {
            let pos = bed_idx_to_pos[&s.bed_idx];
            let maf_ok = args.maf.map(|thr| maf_per_snp[pos] > thr).unwrap_or(true);
            let print_ok = print_set
                .as_ref()
                .map(|set| set.contains(&s.snp))
                .unwrap_or(true);
            maf_ok && print_ok
        })
        .map(|s| bed_idx_to_pos[&s.bed_idx])
        .collect();
    if print_set.is_some() && out_positions.is_empty() {
        anyhow::bail!("After merging with --print-snps, no SNPs remain.");
    }

    let extract_rows = |positions: &[usize], n_cols: usize| -> MatF {
        let mut mat = MatF::zeros(positions.len(), n_cols);
        for (row, &pos) in positions.iter().enumerate() {
            for k in 0..n_cols {
                mat[(row, k)] = l2[(pos, k)];
            }
        }
        mat
    };
    let compute_m_vals = |positions: &[usize]| -> (Vec<f64>, Vec<f64>) {
        match &annot_result {
            Some((annot, _)) => {
                let mut m_vals = vec![0.0f64; col_names.len()];
                let mut m_5_50_vals = vec![0.0f64; col_names.len()];
                for &pos in positions {
                    let maf = maf_per_snp[pos];
                    for k in 0..col_names.len() {
                        let mut base = annot[(pos, k)];
                        if let Some(pq) = pq_for_m {
                            base *= pq[pos];
                        }
                        m_vals[k] += base;
                        if maf > 0.05 {
                            m_5_50_vals[k] += base;
                        }
                    }
                }
                (m_vals, m_5_50_vals)
            }
            None => {
                let m_val = if let Some(pq) = pq_per_snp.as_ref() {
                    positions.iter().map(|&pos| pq[pos]).sum::<f64>()
                } else {
                    positions.len() as f64
                };
                let m_5_50_val = if let Some(pq) = pq_per_snp.as_ref() {
                    positions
                        .iter()
                        .filter(|&&pos| maf_per_snp[pos] > 0.05)
                        .map(|&pos| pq[pos])
                        .sum::<f64>()
                } else {
                    positions
                        .iter()
                        .filter(|&&pos| maf_per_snp[pos] > 0.05)
                        .count() as f64
                };
                (vec![m_val], vec![m_5_50_val])
            }
        }
    };

    let t_write = Instant::now();
    let out_snps: Vec<&BimRecord> = all_snps
        .iter()
        .filter(|s| {
            let pos = bed_idx_to_pos[&s.bed_idx];
            let maf_ok = args.maf.map(|thr| maf_per_snp[pos] > thr).unwrap_or(true);
            let print_ok = print_set
                .as_ref()
                .map(|set| set.contains(&s.snp))
                .unwrap_or(true);
            maf_ok && print_ok
        })
        .collect();
    let out_l2 = extract_rows(&out_positions, col_names.len());
    let out_path = format!("{}.l2.ldscore.gz", args.out);
    write_ldscore_refs(&out_path, &out_snps, &out_l2, &col_names)
        .with_context(|| "writing combined LD score output".to_string())?;

    let all_positions: Vec<usize> = (0..all_snps.len()).collect();
    let (m_vals_all, m_5_50_vals_all) = compute_m_vals(&all_positions);
    let m_path = format!("{}.l2.M", args.out);
    std::fs::write(&m_path, format_m_vals(&m_vals_all))
        .with_context(|| format!("writing M file '{}'", m_path))?;
    let m_5_50_path = format!("{}.l2.M_5_50", args.out);
    std::fs::write(&m_5_50_path, format_m_vals(&m_5_50_vals_all))
        .with_context(|| format!("writing M_5_50 file '{}'", m_5_50_path))?;

    for chr in chrs {
        let chr_snps_all: Vec<&BimRecord> = all_snps.iter().filter(|s| s.chr == chr).collect();
        let chr_positions_all: Vec<usize> = chr_snps_all
            .iter()
            .map(|s| bed_idx_to_pos[&s.bed_idx])
            .collect();

        let chr_snps: Vec<&BimRecord> = chr_snps_all
            .iter()
            .copied()
            .filter(|s| {
                let pos = bed_idx_to_pos[&s.bed_idx];
                let maf_ok = args.maf.map(|thr| maf_per_snp[pos] > thr).unwrap_or(true);
                let print_ok = print_set
                    .as_ref()
                    .map(|set| set.contains(&s.snp))
                    .unwrap_or(true);
                maf_ok && print_ok
            })
            .collect();
        let chr_positions_out: Vec<usize> = chr_snps
            .iter()
            .map(|s| bed_idx_to_pos[&s.bed_idx])
            .collect();
        let n_chr = chr_snps.len();

        let chr_l2 = extract_rows(&chr_positions_out, col_names.len());
        let out_path = format!("{}{}.l2.ldscore.gz", args.out, chr);
        write_ldscore_refs(&out_path, &chr_snps, &chr_l2, &col_names)
            .with_context(|| format!("writing output for chr {}", chr))?;

        let (m_vals, m_5_50_vals) = compute_m_vals(&chr_positions_all);
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
    trace!("l2 timing: write_outputs {:?}", t_write.elapsed());

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
        let mut col = vec![2.0f64, 2.0, 2.0, 2.0];
        normalize_col_f64(&mut col, 4);
        assert!(
            col.iter().all(|&v| v == 0.0),
            "constant column should be all zeros"
        );
    }

    #[test]
    fn test_normalize_col_unit_variance() {
        let mut col = vec![0.0f64, 1.0, 2.0, 1.0];
        let n = col.len();
        normalize_col_f64(&mut col, n);

        let mean: f64 = col.iter().sum::<f64>() / n as f64;
        let var: f64 = col.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / n as f64;
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

        let args = L2Args {
            bfile: prefix,
            out: out.clone(),
            ld_wind_cm: Some(1.0),
            ld_wind_kb: None,
            ld_wind_snp: Some(1),
            annot: None,
            cts_bin: None,
            cts_breaks: None,
            cts_names: None,
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

        let args = L2Args {
            bfile: prefix,
            out: out.clone(),
            ld_wind_cm: Some(1.0),
            ld_wind_kb: None,
            ld_wind_snp: Some(1),
            annot: None,
            cts_bin: None,
            cts_breaks: None,
            cts_names: None,
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
