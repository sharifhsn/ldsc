use crate::bed::{Bed, ChunkReader, ReadOptions};
use crate::cli::L2Args;
use crate::cts_annot;
#[cfg(feature = "gpu")]
use crate::gpu::GpuContext;
use crate::la::{
    MatF, MatF32, mat_add_in_place, mat_copy_from, mat_slice, mat_slice_f32, mat_slice_mut,
    mat_slice_mut_f32, matmul_tn_to, matmul_tn_to_f32, matmul_to,
};
use crate::parse;
use anyhow::{Context, Result};
use faer::{Accum, Par};
use std::collections::{HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, Write as IoWrite};

/// GEMM scratch buffers — either f32 (fast) or f64 (default precision).
enum GemmBufs {
    F64 {
        ring_buf: MatF,
        b_mat: MatF,
        a_buf: MatF,
        ab_buf: MatF,
    },
    F32 {
        ring_buf: MatF32,
        b_mat: MatF32,
        a_buf: MatF32,
        ab_buf: MatF32,
    },
}

impl GemmBufs {
    fn new(
        use_f32: bool,
        n_indiv: usize,
        chunk_c: usize,
        ring_size: usize,
        max_window: usize,
    ) -> Self {
        let mw = max_window.max(1);
        if use_f32 {
            Self::F32 {
                ring_buf: MatF32::zeros(n_indiv, ring_size),
                b_mat: MatF32::zeros(n_indiv, chunk_c),
                a_buf: MatF32::zeros(n_indiv, mw),
                ab_buf: MatF32::zeros(mw, chunk_c),
            }
        } else {
            Self::F64 {
                ring_buf: MatF::zeros(n_indiv, ring_size),
                b_mat: MatF::zeros(n_indiv, chunk_c),
                a_buf: MatF::zeros(n_indiv, mw),
                ab_buf: MatF::zeros(mw, chunk_c),
            }
        }
    }
}

/// Extract column-major f32 data from an f32 faer matrix for GPU upload.
#[cfg(feature = "gpu")]
fn mat_to_col_major_f32_from_f32(
    mat: faer::MatRef<'_, f32>,
    nrows: usize,
    ncols: usize,
) -> Vec<f32> {
    let mut data = Vec::with_capacity(nrows * ncols);
    for j in 0..ncols {
        for i in 0..nrows {
            data.push(mat[(i, j)]);
        }
    }
    data
}

/// Extract column-major f32 data from an f64 faer matrix for GPU upload.
#[cfg(feature = "gpu")]
fn mat_to_col_major_f32_from_f64(
    mat: faer::MatRef<'_, f64>,
    nrows: usize,
    ncols: usize,
) -> Vec<f32> {
    let mut data = Vec::with_capacity(nrows * ncols);
    for j in 0..ncols {
        for i in 0..nrows {
            data.push(mat[(i, j)] as f32);
        }
    }
    data
}

/// Write GPU result (row-major f32) back into an f32 faer matrix.
#[cfg(feature = "gpu")]
fn write_gpu_result_f32(
    result: &[f32],
    mut dst: faer::MatMut<'_, f32>,
    nrows: usize,
    ncols: usize,
) {
    for i in 0..nrows {
        for j in 0..ncols {
            dst[(i, j)] = result[i * ncols + j];
        }
    }
}

/// Write GPU result (row-major f32) back into an f64 faer matrix.
#[cfg(feature = "gpu")]
fn write_gpu_result_f64(
    result: &[f32],
    mut dst: faer::MatMut<'_, f64>,
    nrows: usize,
    ncols: usize,
) {
    for i in 0..nrows {
        for j in 0..ncols {
            dst[(i, j)] = result[i * ncols + j] as f64;
        }
    }
}

/// Extract column-major f64 data from an f64 faer matrix for native f64 GPU upload.
#[cfg(feature = "gpu")]
fn mat_to_col_major_f64(mat: faer::MatRef<'_, f64>, nrows: usize, ncols: usize) -> Vec<f64> {
    let mut data = Vec::with_capacity(nrows * ncols);
    for j in 0..ncols {
        for i in 0..nrows {
            data.push(mat[(i, j)]);
        }
    }
    data
}

/// Write native f64 GPU result (row-major f64) back into an f64 faer matrix.
#[cfg(feature = "gpu")]
fn write_gpu_result_f64_native(
    result: &[f64],
    mut dst: faer::MatMut<'_, f64>,
    nrows: usize,
    ncols: usize,
) {
    for i in 0..nrows {
        for j in 0..ncols {
            dst[(i, j)] = result[i * ncols + j];
        }
    }
}

pub use crate::parse::{BimRecord, parse_bim};

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

/// For each SNP, compute the leftmost SNP index within `max_dist`.
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

/// Normalize genotype column in-place (impute NaN→mean, centre, scale). Returns MAF.
/// When pre-computed sum/count/sum_sq are available (fused with copy), does center+scale
/// in a single pass using variance from raw moments: var = E[X²] - E[X]².
fn normalize_col_f64_with_stats(
    col: &mut [f64],
    n: usize,
    sum: f64,
    count: usize,
    sum_sq: f64,
) -> f64 {
    let avg = if count > 0 { sum / count as f64 } else { 0.0 };
    let freq = (avg / 2.0).clamp(0.0, 1.0);
    let maf = freq.min(1.0 - freq);

    // Variance from raw moments: Var(X) = E[X²] - E[X]²
    // After centering, sum of (x-μ)² = Σx² - n*μ² = sum_sq - count*avg²
    // NaN elements become 0 after imputation, contributing avg² each to centered sum.
    let centered_sum_sq = sum_sq - count as f64 * avg * avg;
    // Missing (NaN→0) elements contribute 0² = 0 to centered col, but we used
    // count*avg² above which only covers non-NaN. The NaN positions are set to 0
    // (centered mean), so they contribute nothing extra.
    let var = centered_sum_sq / n as f64;
    let std = var.sqrt();
    let inv_std = if std > 0.0 { 1.0 / std } else { 0.0 };

    // Single fused pass: impute NaN, center, and scale
    if count == col.len() {
        // Fast path: no NaN — tight branchless loop, auto-vectorizes
        for v in col.iter_mut() {
            *v = (*v - avg) * inv_std;
        }
    } else {
        // Slow path: has NaN
        for v in col.iter_mut() {
            if v.is_nan() {
                *v = 0.0;
            } else {
                *v = (*v - avg) * inv_std;
            }
        }
    }
    maf
}

/// f32 variant: normalize in-place using f64 accumulators for precision.
/// When pre-computed sum/count/sum_sq are available (fused with copy), does center+scale
/// in a single pass.
fn normalize_col_f32_with_stats(
    col: &mut [f32],
    n: usize,
    sum: f64,
    count: usize,
    sum_sq: f64,
) -> f32 {
    let avg = if count > 0 { sum / count as f64 } else { 0.0 };
    let freq = (avg / 2.0).clamp(0.0, 1.0);
    let maf = freq.min(1.0 - freq);

    // Variance from raw moments: Var(X) = E[X²] - E[X]²
    let centered_sum_sq = sum_sq - count as f64 * avg * avg;
    let var = centered_sum_sq / n as f64;
    let std = var.sqrt();
    let inv_std_f32 = if std > 0.0 { (1.0 / std) as f32 } else { 0.0f32 };
    let avg_f32 = avg as f32;

    // Single fused pass: impute NaN, center, and scale
    if count == col.len() {
        // Fast path: no NaN — tight branchless loop, auto-vectorizes
        for v in col.iter_mut() {
            *v = (*v - avg_f32) * inv_std_f32;
        }
    } else {
        // Slow path: has NaN
        for v in col.iter_mut() {
            if v.is_nan() {
                *v = 0.0;
            } else {
                *v = (*v - avg_f32) * inv_std_f32;
            }
        }
    }
    maf as f32
}

/// Compute LD scores for all SNPs. Returns `(l2, maf)`.
#[allow(clippy::too_many_arguments, clippy::unnecessary_cast)]
fn compute_ldscore_global(
    all_snps: &[BimRecord],
    bed_path: &str,
    n_indiv: usize,
    mode: WindowMode,
    chunk_c: usize,
    annot: Option<&MatF>,
    iid_indices: Option<&[isize]>,
    pq_exp: Option<f64>,
    yes_really: bool,
    _use_gpu: bool,
    _gpu_tile_cols: Option<usize>,
    _gpu_flex32: bool,
    _gpu_f64: bool,
    use_f32: bool,
    prefetch_bed: bool,
    verbose_timing: bool,
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
    let n_inv = 1.0 / n;
    // Pre-compute r2_unbiased linear constants: r2u(r) = r2u_a * r² + r2u_b
    // r2_unbiased(r, n) = r² - (1 - r²) / (n - 2) = r² * (1 + 1/(n-2)) - 1/(n-2)
    let r2u_denom = if n_indiv > 2 { n - 2.0 } else { n };
    let r2u_denom_inv = 1.0 / r2u_denom;
    let r2u_a = 1.0 + r2u_denom_inv; // coefficient of r²
    let r2u_b = -r2u_denom_inv; // constant term
    // Pre-compute n_inv² · r2u_a so hot loop does 2 muls instead of 3:
    // r2u(val) = val² · n_inv² · r2u_a + r2u_b  (was: r = val·n_inv; r²·r2u_a + r2u_b)
    let n_inv_sq_r2u_a = n_inv * n_inv * r2u_a;
    let mut maf_per_snp = vec![0.0f64; m];

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
    let mut bufs = GemmBufs::new(use_f32, n_indiv, chunk_c, ring_size, max_window_size);
    let mut ring_next: usize = 0;
    let mut window: VecDeque<(usize, usize)> = VecDeque::new(); // (snp_idx, ring_slot)

    let mut pq_per_snp: Option<Vec<f64>> = pq_exp.map(|_| vec![1.0f64; m]);

    // When no annotation is provided, use an all-ones column (scalar LD scores).
    let ones_annot;
    let annot = match annot {
        Some(a) => a,
        None => {
            ones_annot = MatF::from_fn(m, 1, |_, _| 1.0);
            &ones_annot
        }
    };
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

    // Pre-allocate scratch matrices used per chunk (avoids alloc+zero every iteration).
    // All are written with Accum::Replace (full overwrite) so stale data is harmless.
    let mut bb_f64 = if !use_f32 { MatF::zeros(chunk_c, chunk_c) } else { MatF::zeros(0, 0) };
    let mut bb_f32 = if use_f32 { MatF32::zeros(chunk_c, chunk_c) } else { MatF32::zeros(0, 0) };
    let mut contrib_bb = MatF::zeros(chunk_c, n_annot);
    let mut contrib_left = MatF::zeros(max_window_size.max(1), n_annot);
    let mut contrib_right = MatF::zeros(chunk_c, n_annot);

    use std::time::Instant;
    let mut t_bed_read = std::time::Duration::ZERO;
    let mut t_norm = std::time::Duration::ZERO;
    let mut t_bb_dot = std::time::Duration::ZERO;
    let mut t_ab_dot = std::time::Duration::ZERO;
    let mut t_ring_store = std::time::Duration::ZERO;

    #[cfg(feature = "gpu")]
    let gpu_ctx = if _use_gpu {
        match GpuContext::new(_gpu_flex32) {
            Ok(ctx) => {
                if _gpu_f64 && !ctx.capabilities.has_f64 {
                    eprintln!(
                        "GPU: warning: --gpu-f64 requested but f64 arithmetic not supported; \
                         falling back to f32 conversion"
                    );
                }
                Some(ctx)
            }
            Err(e) => {
                eprintln!("GPU: initialization failed ({}), falling back to CPU", e);
                None
            }
        }
    } else {
        None
    };

    // --prefetch-bed: spawn a background reader thread that reads BED chunks one ahead
    // of the compute loop, overlapping I/O with GEMM.
    // Only beneficial when BED reads involve real blocking I/O (cold GPFS/NFS).
    // On local SSD with warm page cache this adds thread contention with rayon and
    // regresses ~10%; leave off for local benchmarks.
    type ChunkMsg = anyhow::Result<(usize, usize, MatF32)>;
    let prefetch_rx: Option<crossbeam_channel::Receiver<ChunkMsg>>;
    let reader_handle: Option<std::thread::JoinHandle<()>>;

    if prefetch_bed {
        struct ChunkSpec {
            chunk_start: usize,
            chunk_end: usize,
            bed_indices: Vec<isize>,
        }
        let chunk_specs: Vec<ChunkSpec> = (0..m)
            .step_by(chunk_c)
            .map(|start| {
                let end = (start + chunk_c).min(m);
                let bed_indices = all_snps[start..end]
                    .iter()
                    .map(|s| s.bed_idx as isize)
                    .collect();
                ChunkSpec { chunk_start: start, chunk_end: end, bed_indices }
            })
            .collect();

        let (filled_tx, filled_rx) = crossbeam_channel::bounded::<ChunkMsg>(1);
        let bed_path_owned = bed_path.to_string();
        let iid_indices_owned: Option<Vec<isize>> = iid_indices.map(|s| s.to_vec());

        let handle = std::thread::spawn(move || {
            let mut bed = match Bed::builder(bed_path_owned.as_str()).build() {
                Ok(b) => b,
                Err(e) => {
                    let _ = filled_tx.send(Err(e));
                    return;
                }
            };
            for spec in chunk_specs {
                let result: anyhow::Result<MatF32> = {
                    let mut ro = ReadOptions::builder();
                    let mut builder = ro.sid_index(&spec.bed_indices).f32();
                    let raw = if let Some(ref iids) = iid_indices_owned {
                        builder.iid_index(iids).read(&mut bed)
                    } else {
                        builder.read(&mut bed)
                    };
                    raw.with_context(|| {
                        format!("reading BED chunk [{},{})", spec.chunk_start, spec.chunk_end)
                    })
                };
                let msg = result.map(|raw| (spec.chunk_start, spec.chunk_end, raw));
                if filled_tx.send(msg).is_err() {
                    break;
                }
            }
        });
        prefetch_rx = Some(filled_rx);
        reader_handle = Some(handle);
    } else {
        prefetch_rx = None;
        reader_handle = None;
    }

    // Sequential path: open BED once with pre-computed ChunkReader; used only when prefetch_bed=false.
    let mut seq_bed = if !prefetch_bed {
        Some(Bed::builder(bed_path).build().context("opening BED file")?)
    } else {
        None
    };
    let mut chunk_reader: Option<ChunkReader<f32>> = if !prefetch_bed {
        let bed = seq_bed.as_ref().unwrap();
        Some(
            ChunkReader::new(bed, iid_indices, chunk_c, true, f32::NAN)
                .context("creating chunk reader")?,
        )
    } else {
        None
    };
    // Check if ALL SNPs are contiguous in the BED file (common: no --extract).
    // If so, use sequential reads (1 seek + streaming reads) instead of seeking per chunk.
    // BufReader::seek(SeekFrom::Start) discards its internal buffer every time,
    // causing 64× read amplification (8MB buffer filled, 124KB used, buffer discarded).
    let all_bed_sequential = m > 0
        && all_snps
            .windows(2)
            .all(|w| w[1].bed_idx == w[0].bed_idx + 1);
    if all_bed_sequential && !prefetch_bed {
        let bed = seq_bed.as_mut().unwrap();
        let cr = chunk_reader.as_mut().unwrap();
        cr.start_sequential(bed, all_snps[0].bed_idx)
            .context("starting sequential BED read")?;
    }

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

        let t0 = Instant::now();
        // Use a reference to avoid moving the owned Mat from prefetch path
        let prefetch_raw: MatF32;
        let raw_ref: &MatF32 = if let Some(ref rx) = prefetch_rx {
            // Prefetch path: receive pre-read chunk from reader thread.
            let (cs, ce, r) = rx.recv().context("prefetch reader closed early")??;
            debug_assert_eq!(cs, chunk_start);
            debug_assert_eq!(ce, chunk_end);
            prefetch_raw = r;
            &prefetch_raw
        } else if all_bed_sequential {
            // Fast path: sequential BED — no seek, BufReader stays warm.
            let bed = seq_bed.as_mut().unwrap();
            let cr = chunk_reader.as_mut().unwrap();
            cr.read_next(bed, c)
                .with_context(|| format!("reading BED chunk [{},{})", chunk_start, chunk_end))?
        } else {
            // Per-chunk contiguous check (--extract with some contiguous ranges)
            let start_bed_idx = all_snps[chunk_start].bed_idx;
            let is_bed_contiguous = all_snps[chunk_end - 1].bed_idx == start_bed_idx + c - 1;
            if is_bed_contiguous {
                let bed = seq_bed.as_mut().unwrap();
                let cr = chunk_reader.as_mut().unwrap();
                cr.read_contiguous(bed, start_bed_idx, c)
                    .with_context(|| {
                        format!("reading BED chunk [{},{})", chunk_start, chunk_end)
                    })?
            } else {
                // Fallback: non-contiguous BED indices (--extract with gaps)
                let bed_indices: Vec<isize> = all_snps[chunk_start..chunk_end]
                    .iter()
                    .map(|s| s.bed_idx as isize)
                    .collect();
                let bed = seq_bed.as_mut().unwrap();
                let mut ro = ReadOptions::builder();
                let mut builder = ro.sid_index(&bed_indices).f32();
                let result = if let Some(iids) = iid_indices {
                    builder.iid_index(iids).read(bed)
                } else {
                    builder.read(bed)
                };
                prefetch_raw = result.with_context(|| {
                    format!("reading BED chunk [{},{})", chunk_start, chunk_end)
                })?;
                &prefetch_raw
            }
        };
        t_bed_read += t0.elapsed();

        let t0 = Instant::now();
        for j in 0..c {
            // Get raw column as contiguous slice (1 bounds check vs 2 per element)
            let raw_col = raw_ref.col(j).try_as_col_major().unwrap().as_slice();
            match bufs {
                GemmBufs::F32 {
                    ref mut b_mat, ..
                } => {
                    // Fused copy + sum + sum_sq: unconditional (NaN check after)
                    let col = b_mat.col_mut(j).try_as_col_major_mut().unwrap().as_slice_mut();
                    let mut sum = 0f64;
                    let mut sum_sq = 0f64;
                    // Tight copy+sum+sum_sq loop — no branch, auto-vectorizes
                    for i in 0..n_indiv {
                        let v = raw_col[i];
                        col[i] = v;
                        let vf = v as f64;
                        sum += vf;
                        sum_sq += vf * vf;
                    }
                    // NaN propagates: if any element was NaN, sum is NaN
                    let (sum, count, sum_sq) = if sum.is_nan() {
                        // Rare slow path: recompute excluding NaN
                        let mut s = 0f64;
                        let mut sq = 0f64;
                        let mut c = 0usize;
                        for &v in col.iter() {
                            if !v.is_nan() {
                                let vf = v as f64;
                                s += vf;
                                sq += vf * vf;
                                c += 1;
                            }
                        }
                        (s, c, sq)
                    } else {
                        (sum, n_indiv, sum_sq)
                    };
                    let snp_maf = normalize_col_f32_with_stats(col, n_indiv, sum, count, sum_sq);
                    maf_per_snp[chunk_start + j] = snp_maf as f64;
                    if let (Some(exp), Some(ref mut pq_vec)) = (pq_exp, pq_per_snp.as_mut()) {
                        let maf_f64 = snp_maf as f64;
                        let pq = (maf_f64 * (1.0 - maf_f64)).powf(exp);
                        pq_vec[chunk_start + j] = pq;
                    }
                }
                GemmBufs::F64 {
                    ref mut b_mat, ..
                } => {
                    // Fused copy + sum + sum_sq: unconditional (NaN check after)
                    let col = b_mat.col_mut(j).try_as_col_major_mut().unwrap().as_slice_mut();
                    let mut sum = 0f64;
                    let mut sum_sq = 0f64;
                    // Tight f32→f64 copy+sum+sum_sq loop — no branch, auto-vectorizes
                    for i in 0..n_indiv {
                        let v = raw_col[i] as f64;
                        col[i] = v;
                        sum += v;
                        sum_sq += v * v;
                    }
                    // NaN propagates: if any element was NaN, sum is NaN
                    let (sum, count, sum_sq) = if sum.is_nan() {
                        // Rare slow path: recompute excluding NaN
                        let mut s = 0f64;
                        let mut sq = 0f64;
                        let mut c = 0usize;
                        for &v in col.iter() {
                            if !v.is_nan() {
                                s += v;
                                sq += v * v;
                                c += 1;
                            }
                        }
                        (s, c, sq)
                    } else {
                        (sum, n_indiv, sum_sq)
                    };
                    let snp_maf = normalize_col_f64_with_stats(col, n_indiv, sum, count, sum_sq);
                    maf_per_snp[chunk_start + j] = snp_maf;
                    if let (Some(exp), Some(ref mut pq_vec)) = (pq_exp, pq_per_snp.as_mut()) {
                        let pq = (snp_maf * (1.0 - snp_maf)).powf(exp);
                        pq_vec[chunk_start + j] = pq;
                    }
                }
            }
        }

        t_norm += t0.elapsed();

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

        // B×B — compute GEMM then extract r2_unbiased into r2u_bb (always f64).
        let t0 = Instant::now();
        if !chunk_all_zero {
            // Inline helper: fill r2u_bb from bb matrix values.
            // Generic over closure to avoid dyn dispatch in inner loop.
            // Uses pre-computed n_inv²·r2u_a to save one multiply per iteration.
            #[inline(always)]
            fn fill_r2u_bb_from<F: Fn(usize, usize) -> f64>(
                r2u_bb: &mut MatF, c: usize, n_inv_sq_r2u_a: f64, r2u_b: f64, bb_val: F,
            ) {
                // Column-major zeroing: inner loop over rows (stride-1)
                for k in 0..c {
                    for j in 0..c {
                        r2u_bb[(j, k)] = 0.0;
                    }
                }
                for j in 0..c {
                    r2u_bb[(j, j)] = 1.0;
                    for k in 0..j {
                        let val = bb_val(k, j);
                        let r2u = val * val * n_inv_sq_r2u_a + r2u_b;
                        r2u_bb[(j, k)] = r2u;
                        r2u_bb[(k, j)] = r2u;
                    }
                }
            }

            match bufs {
                GemmBufs::F32 { ref b_mat, .. } => {
                    let b_slice = mat_slice_f32(b_mat.as_ref(), 0..n_indiv, 0..c);
                    #[allow(unused_mut)]
                    let mut bb_sl = mat_slice_mut_f32(bb_f32.as_mut(), 0..c, 0..c);
                    let mut _did_gpu = false;
                    #[cfg(feature = "gpu")]
                    {
                        if let Some(ref ctx) = gpu_ctx {
                            let b_f32 = mat_to_col_major_f32_from_f32(b_slice, n_indiv, c);
                            let gpu_result = if let Some(tc) = _gpu_tile_cols {
                                if _gpu_flex32 {
                                    ctx.matmul_tn_tiled_flex32(&b_f32, n_indiv, c, &b_f32, c, tc)
                                } else {
                                    ctx.matmul_tn_tiled(&b_f32, n_indiv, c, &b_f32, c, tc)
                                }
                            } else if _gpu_flex32 {
                                ctx.matmul_tn_flex32(&b_f32, n_indiv, c, &b_f32, c)
                            } else {
                                ctx.matmul_tn(&b_f32, n_indiv, c, &b_f32, c)
                            };
                            match gpu_result {
                                Ok(result) => {
                                    write_gpu_result_f32(&result, bb_sl.as_mut(), c, c);
                                    _did_gpu = true;
                                }
                                Err(e) => eprintln!("GPU: f32 B×B matmul failed ({e}), falling back to CPU"),
                            }
                        }
                    }
                    if !_did_gpu {
                        matmul_tn_to_f32(
                            bb_sl,
                            b_slice,
                            b_slice,
                            1.0f32,
                            Accum::Replace,
                            Par::rayon(0),
                        );
                    }
                    fill_r2u_bb_from(&mut r2u_bb, c, n_inv_sq_r2u_a, r2u_b, |k, j| bb_f32[(k, j)] as f64);
                }
                GemmBufs::F64 { ref b_mat, .. } => {
                    let b_slice = mat_slice(b_mat.as_ref(), 0..n_indiv, 0..c);
                    #[allow(unused_mut)]
                    let mut bb_sl = mat_slice_mut(bb_f64.as_mut(), 0..c, 0..c);
                    let mut _did_gpu = false;
                    #[cfg(feature = "gpu")]
                    {
                        if let Some(ref ctx) = gpu_ctx {
                            if _gpu_f64 && ctx.capabilities.has_f64 {
                                // Native f64 path — no precision loss
                                let b_f64 = mat_to_col_major_f64(b_slice, n_indiv, c);
                                let gpu_result = if let Some(tc) = _gpu_tile_cols {
                                    ctx.matmul_tn_tiled_f64(&b_f64, n_indiv, c, &b_f64, c, tc)
                                } else {
                                    ctx.matmul_tn_f64(&b_f64, n_indiv, c, &b_f64, c)
                                };
                                match gpu_result {
                                    Ok(result) => {
                                        write_gpu_result_f64_native(&result, bb_sl.as_mut(), c, c);
                                        _did_gpu = true;
                                    }
                                    Err(e) => eprintln!("GPU: f64 B×B matmul failed ({e}), falling back to CPU"),
                                }
                            } else {
                                // f32 conversion path
                                let b_f32 = mat_to_col_major_f32_from_f64(b_slice, n_indiv, c);
                                let gpu_result = if let Some(tc) = _gpu_tile_cols {
                                    if _gpu_flex32 {
                                        ctx.matmul_tn_tiled_flex32(&b_f32, n_indiv, c, &b_f32, c, tc)
                                    } else {
                                        ctx.matmul_tn_tiled(&b_f32, n_indiv, c, &b_f32, c, tc)
                                    }
                                } else if _gpu_flex32 {
                                    ctx.matmul_tn_flex32(&b_f32, n_indiv, c, &b_f32, c)
                                } else {
                                    ctx.matmul_tn(&b_f32, n_indiv, c, &b_f32, c)
                                };
                                match gpu_result {
                                    Ok(result) => {
                                        write_gpu_result_f64(&result, bb_sl.as_mut(), c, c);
                                        _did_gpu = true;
                                    }
                                    Err(e) => eprintln!("GPU: B×B matmul failed ({e}), falling back to CPU"),
                                }
                            }
                        }
                    }
                    if !_did_gpu {
                        matmul_tn_to(
                            bb_sl,
                            b_slice,
                            b_slice,
                            1.0f64,
                            Accum::Replace,
                            Par::rayon(0),
                        );
                    }
                    fill_r2u_bb_from(&mut r2u_bb, c, n_inv_sq_r2u_a, r2u_b, |k, j| bb_f64[(k, j)]);
                }
            }
            let r2u_bb_view = mat_slice(r2u_bb.as_ref(), 0..c, 0..c);
            // Small matmul (c×c @ c×n_annot) — Par::Seq avoids thread-pool overhead
            if let Some(ref eff) = annot_chunk_eff {
                matmul_to(
                    mat_slice_mut(contrib_bb.as_mut(), 0..c, 0..n_annot),
                    r2u_bb_view,
                    eff.as_ref(),
                    1.0,
                    Accum::Replace,
                    Par::Seq,
                );
            } else {
                matmul_to(
                    mat_slice_mut(contrib_bb.as_mut(), 0..c, 0..n_annot),
                    r2u_bb_view,
                    annot_chunk,
                    1.0,
                    Accum::Replace,
                    Par::Seq,
                );
            }
            mat_add_in_place(
                mat_slice_mut(l2.as_mut(), chunk_start..chunk_end, 0..n_annot),
                mat_slice(contrib_bb.as_ref(), 0..c, 0..n_annot),
            );
        }
        t_bb_dot += t0.elapsed();

        // A×B
        let t0 = Instant::now();
        if !window.is_empty() {
            let w = window.len();
            let a_left_idx = window.front().map(|(idx, _)| *idx).unwrap_or(0);
            let (first_slot, last_slot) = window
                .front()
                .and_then(|(_, f)| window.back().map(|(_, l)| (*f, *l)))
                .unwrap_or((0, 0));
            let contiguous = first_slot <= last_slot && last_slot - first_slot + 1 == w;

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

            if !(chunk_all_zero && window_all_zero) {
                // Compute A^T × B GEMM, then fill r2u_ab (always f64).
                match bufs {
                    GemmBufs::F32 {
                        ref mut ring_buf,
                        ref mut b_mat,
                        ref mut a_buf,
                        ref mut ab_buf,
                        ..
                    } => {
                        if !contiguous {
                            for (wi, (_, slot)) in window.iter().enumerate() {
                                a_buf
                                    .as_mut()
                                    .submatrix_mut(0, wi, n_indiv, 1)
                                    .copy_from(ring_buf.as_ref().submatrix(0, *slot, n_indiv, 1));
                            }
                        }
                        let a_view = if contiguous {
                            mat_slice_f32(
                                ring_buf.as_ref(),
                                0..n_indiv,
                                first_slot..(first_slot + w),
                            )
                        } else {
                            mat_slice_f32(a_buf.as_ref(), 0..n_indiv, 0..w)
                        };
                        let b_sl = mat_slice_f32(b_mat.as_ref(), 0..n_indiv, 0..c);

                        let mut _did_gpu = false;
                        #[cfg(feature = "gpu")]
                        {
                            if let Some(ref ctx) = gpu_ctx {
                                let a_f32 = mat_to_col_major_f32_from_f32(a_view, n_indiv, w);
                                let b_f32 = mat_to_col_major_f32_from_f32(b_sl, n_indiv, c);
                                let gpu_result = if let Some(tc) = _gpu_tile_cols {
                                    if _gpu_flex32 {
                                        ctx.matmul_tn_tiled_flex32(
                                            &a_f32, n_indiv, w, &b_f32, c, tc,
                                        )
                                    } else {
                                        ctx.matmul_tn_tiled(&a_f32, n_indiv, w, &b_f32, c, tc)
                                    }
                                } else if _gpu_flex32 {
                                    ctx.matmul_tn_flex32(&a_f32, n_indiv, w, &b_f32, c)
                                } else {
                                    ctx.matmul_tn(&a_f32, n_indiv, w, &b_f32, c)
                                };
                                match gpu_result {
                                    Ok(result) => {
                                        write_gpu_result_f32(
                                            &result,
                                            mat_slice_mut_f32(ab_buf.as_mut(), 0..w, 0..c),
                                            w,
                                            c,
                                        );
                                        _did_gpu = true;
                                    }
                                    Err(e) => eprintln!("GPU: f32 A×B matmul failed ({e}), falling back to CPU"),
                                }
                            }
                        }
                        if !_did_gpu {
                            matmul_tn_to_f32(
                                mat_slice_mut_f32(ab_buf.as_mut(), 0..w, 0..c),
                                a_view,
                                b_sl,
                                1.0f32,
                                Accum::Replace,
                                Par::rayon(0),
                            );
                        }
                        let ab_view = mat_slice_f32(ab_buf.as_ref(), 0..w, 0..c);
                        // Column-major inner loop with column slices: eliminates faer
                        // tuple-indexing bounds checks (~8.3B checks for full genome).
                        // Pre-computed n_inv²·r2u_a saves one multiply per iteration.
                        for j in 0..c {
                            let ab_col = ab_view.col(j).try_as_col_major().unwrap().as_slice();
                            let r2u_col = r2u_ab.col_mut(j).try_as_col_major_mut().unwrap().as_slice_mut();
                            for (ab_val, r2u_val) in ab_col[..w].iter().zip(r2u_col[..w].iter_mut()) {
                                let val = *ab_val as f64;
                                *r2u_val = val * val * n_inv_sq_r2u_a + r2u_b;
                            }
                        }
                    }
                    GemmBufs::F64 {
                        ref mut ring_buf,
                        ref mut b_mat,
                        ref mut a_buf,
                        ref mut ab_buf,
                        ..
                    } => {
                        if !contiguous {
                            for (wi, (_, slot)) in window.iter().enumerate() {
                                a_buf
                                    .as_mut()
                                    .submatrix_mut(0, wi, n_indiv, 1)
                                    .copy_from(ring_buf.as_ref().submatrix(0, *slot, n_indiv, 1));
                            }
                        }
                        let a_view = if contiguous {
                            mat_slice(ring_buf.as_ref(), 0..n_indiv, first_slot..(first_slot + w))
                        } else {
                            mat_slice(a_buf.as_ref(), 0..n_indiv, 0..w)
                        };
                        let b_sl = mat_slice(b_mat.as_ref(), 0..n_indiv, 0..c);

                        let mut _did_gpu = false;
                        #[cfg(feature = "gpu")]
                        {
                            if let Some(ref ctx) = gpu_ctx {
                                if _gpu_f64 && ctx.capabilities.has_f64 {
                                    // Native f64 path
                                    let a_f64_data = mat_to_col_major_f64(a_view, n_indiv, w);
                                    let b_f64_data = mat_to_col_major_f64(b_sl, n_indiv, c);
                                    let gpu_result = if let Some(tc) = _gpu_tile_cols {
                                        ctx.matmul_tn_tiled_f64(
                                            &a_f64_data, n_indiv, w, &b_f64_data, c, tc,
                                        )
                                    } else {
                                        ctx.matmul_tn_f64(
                                            &a_f64_data, n_indiv, w, &b_f64_data, c,
                                        )
                                    };
                                    match gpu_result {
                                        Ok(result) => {
                                            write_gpu_result_f64_native(
                                                &result,
                                                mat_slice_mut(ab_buf.as_mut(), 0..w, 0..c),
                                                w,
                                                c,
                                            );
                                            _did_gpu = true;
                                        }
                                        Err(e) => eprintln!("GPU: f64 A×B matmul failed ({e}), falling back to CPU"),
                                    }
                                } else {
                                    // f32 conversion path
                                    let a_f32 = mat_to_col_major_f32_from_f64(a_view, n_indiv, w);
                                    let b_f32 = mat_to_col_major_f32_from_f64(b_sl, n_indiv, c);
                                    let gpu_result = if let Some(tc) = _gpu_tile_cols {
                                        if _gpu_flex32 {
                                            ctx.matmul_tn_tiled_flex32(
                                                &a_f32, n_indiv, w, &b_f32, c, tc,
                                            )
                                        } else {
                                            ctx.matmul_tn_tiled(&a_f32, n_indiv, w, &b_f32, c, tc)
                                        }
                                    } else if _gpu_flex32 {
                                        ctx.matmul_tn_flex32(&a_f32, n_indiv, w, &b_f32, c)
                                    } else {
                                        ctx.matmul_tn(&a_f32, n_indiv, w, &b_f32, c)
                                    };
                                    match gpu_result {
                                        Ok(result) => {
                                            write_gpu_result_f64(
                                                &result,
                                                mat_slice_mut(ab_buf.as_mut(), 0..w, 0..c),
                                                w,
                                                c,
                                            );
                                            _did_gpu = true;
                                        }
                                        Err(e) => eprintln!("GPU: A×B matmul failed ({e}), falling back to CPU"),
                                    }
                                }
                            }
                        }
                        if !_did_gpu {
                            matmul_tn_to(
                                mat_slice_mut(ab_buf.as_mut(), 0..w, 0..c),
                                a_view,
                                b_sl,
                                1.0f64,
                                Accum::Replace,
                                Par::rayon(0),
                            );
                        }
                        let ab_view = mat_slice(ab_buf.as_ref(), 0..w, 0..c);
                        // Column-major inner loop with column slices: eliminates faer
                        // tuple-indexing bounds checks. Pre-computed constant saves 1 mul.
                        for j in 0..c {
                            let ab_col = ab_view.col(j).try_as_col_major().unwrap().as_slice();
                            let r2u_col = r2u_ab.col_mut(j).try_as_col_major_mut().unwrap().as_slice_mut();
                            for (ab_val, r2u_val) in ab_col[..w].iter().zip(r2u_col[..w].iter_mut()) {
                                *r2u_val = *ab_val * *ab_val * n_inv_sq_r2u_a + r2u_b;
                            }
                        }
                    }
                }

                // Small matmuls (w×c @ c×n_annot, c×w @ w×n_annot) — Par::Seq
                if !chunk_all_zero {
                    let r2u_ab_view = mat_slice(r2u_ab.as_ref(), 0..w, 0..c);
                    let cl_sl = mat_slice_mut(contrib_left.as_mut(), 0..w, 0..n_annot);
                    if let Some(ref eff) = annot_chunk_eff {
                        matmul_to(
                            cl_sl,
                            r2u_ab_view,
                            eff.as_ref(),
                            1.0,
                            Accum::Replace,
                            Par::Seq,
                        );
                    } else {
                        matmul_to(
                            cl_sl,
                            r2u_ab_view,
                            annot_chunk,
                            1.0,
                            Accum::Replace,
                            Par::Seq,
                        );
                    }
                    mat_add_in_place(
                        mat_slice_mut(l2.as_mut(), a_left_idx..chunk_start, 0..n_annot),
                        mat_slice(contrib_left.as_ref(), 0..w, 0..n_annot),
                    );
                }

                if !window_all_zero {
                    let r2u_ab_view = mat_slice(r2u_ab.as_ref(), 0..w, 0..c);
                    let cr_sl = mat_slice_mut(contrib_right.as_mut(), 0..c, 0..n_annot);
                    if let Some(ref eff) = annot_window_eff {
                        matmul_tn_to(
                            cr_sl,
                            r2u_ab_view,
                            eff.as_ref(),
                            1.0,
                            Accum::Replace,
                            Par::Seq,
                        );
                    } else {
                        matmul_tn_to(
                            cr_sl,
                            r2u_ab_view,
                            annot_window,
                            1.0,
                            Accum::Replace,
                            Par::Seq,
                        );
                    }
                    mat_add_in_place(
                        mat_slice_mut(l2.as_mut(), chunk_start..chunk_end, 0..n_annot),
                        mat_slice(contrib_right.as_ref(), 0..c, 0..n_annot),
                    );
                }
            }
        }

        t_ab_dot += t0.elapsed();

        // Push B columns into ring buffer.
        // When all c slots are contiguous (no wrap), do a single bulk copy.
        let t0 = Instant::now();
        let start_slot = ring_next % ring_size;
        let end_slot = start_slot + c; // may exceed ring_size if wrapping
        let bulk_ok = end_slot <= ring_size;
        match bufs {
            GemmBufs::F32 {
                ref mut ring_buf,
                ref b_mat,
                ..
            } => {
                if bulk_ok {
                    ring_buf
                        .as_mut()
                        .submatrix_mut(0, start_slot, n_indiv, c)
                        .copy_from(b_mat.as_ref().submatrix(0, 0, n_indiv, c));
                    for j in 0..c {
                        window.push_back((chunk_start + j, start_slot + j));
                    }
                } else {
                    for j in 0..c {
                        let slot = (ring_next + j) % ring_size;
                        ring_buf
                            .as_mut()
                            .submatrix_mut(0, slot, n_indiv, 1)
                            .copy_from(b_mat.as_ref().submatrix(0, j, n_indiv, 1));
                        window.push_back((chunk_start + j, slot));
                    }
                }
            }
            GemmBufs::F64 {
                ref mut ring_buf,
                ref b_mat,
                ..
            } => {
                if bulk_ok {
                    ring_buf
                        .as_mut()
                        .submatrix_mut(0, start_slot, n_indiv, c)
                        .copy_from(b_mat.as_ref().submatrix(0, 0, n_indiv, c));
                    for j in 0..c {
                        window.push_back((chunk_start + j, start_slot + j));
                    }
                } else {
                    for j in 0..c {
                        let slot = (ring_next + j) % ring_size;
                        ring_buf
                            .as_mut()
                            .submatrix_mut(0, slot, n_indiv, 1)
                            .copy_from(b_mat.as_ref().submatrix(0, j, n_indiv, 1));
                        window.push_back((chunk_start + j, slot));
                    }
                }
            }
        }
        ring_next += c;
        t_ring_store += t0.elapsed();
    } // end for chunk_start

    if let Some(handle) = reader_handle {
        handle.join().expect("prefetch reader thread panicked");
    }

    if verbose_timing {
        eprintln!(
            "[perf] compute_ldscore: bed_read(stall)={:.3}s norm={:.3}s bb_dot={:.3}s ab_dot={:.3}s ring_store={:.3}s",
            t_bed_read.as_secs_f64(),
            t_norm.as_secs_f64(),
            t_bb_dot.as_secs_f64(),
            t_ab_dot.as_secs_f64(),
            t_ring_store.as_secs_f64(),
        );
    }

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
            bed.read_snp_block(start_sid, c, &mut block_buf)
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
                bed.read_snp_bytes(sid, &mut buf)
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
/// Inlined as `r*r * r2u_a + r2u_b` in hot path for constant n.
#[inline]
#[allow(dead_code)]
fn r2_unbiased(r: f64, n: usize) -> f64 {
    let sq = r * r;
    let denom = if n > 2 { n as f64 - 2.0 } else { n as f64 };
    sq - (1.0 - sq) / denom
}

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

fn format_m_vals(vals: &[f64]) -> String {
    let joined = vals
        .iter()
        .map(|v| format!("{v}"))
        .collect::<Vec<_>>()
        .join("\t");
    format!("{joined}\n")
}

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
    // Skip when --maf and --pq-exp are both absent: the filter removes essentially no SNPs
    // from well-filtered data (1000G), and costs an extra full BED pass (~4s).
    let verbose_timing = args.verbose_timing;
    let t_run_start = std::time::Instant::now();
    let mut maf_prefilter: Option<Vec<f64>> = None;
    if maf_pre && (args.maf.is_some() || pq_exp.is_some()) {
        let mut bed = Bed::builder(bed_path.as_str())
            .build()
            .context("opening BED file for SNP prefilter")?;
        let (maf_all, het_miss_ok) = compute_snp_stats(
            &all_snps,
            &mut bed,
            n_indiv_actual,
            args.chunk_size,
            iid_indices.as_deref(),
        )
        .context("computing SNP prefilter stats")?;

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

    let t_maf_pre = t_run_start.elapsed();
    if verbose_timing {
        eprintln!("[perf] maf_prefilter={:.3}s", t_maf_pre.as_secs_f64());
    }

    let t_compute_start = std::time::Instant::now();
    let (l2, maf_per_snp) = compute_ldscore_global(
        &all_snps,
        bed_path.as_str(),
        n_indiv_actual,
        mode,
        chunk_c,
        annot_ref,
        iid_indices.as_deref(),
        pq_exp_for_compute,
        args.yes_really,
        args.gpu,
        args.gpu_tile_cols,
        args.gpu_flex32,
        args.gpu_f64,
        args.fast_f32,
        args.prefetch_bed,
        args.verbose_timing,
    )
    .context("computing LD scores")?;

    if verbose_timing {
        eprintln!(
            "[perf] compute_ldscore_total={:.3}s",
            t_compute_start.elapsed().as_secs_f64()
        );
    }

    let t_write_start = std::time::Instant::now();
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

    let should_output = |s: &BimRecord| -> bool {
        let pos = bed_idx_to_pos[&s.bed_idx];
        let maf_ok = args.maf.map(|thr| maf_per_snp[pos] > thr).unwrap_or(true);
        let print_ok = print_set
            .as_ref()
            .map(|set| set.contains(&s.snp))
            .unwrap_or(true);
        maf_ok && print_ok
    };

    let out_positions: Vec<usize> = all_snps
        .iter()
        .filter(|s| should_output(s))
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

    let out_snps: Vec<&BimRecord> = all_snps.iter().filter(|s| should_output(s)).collect();
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

    use rayon::prelude::*;
    let chr_messages: Vec<String> = chrs
        .par_iter()
        .map(|chr| -> Result<String> {
            let chr_snps_all: Vec<&BimRecord> =
                all_snps.iter().filter(|s| s.chr == *chr).collect();
            let chr_positions_all: Vec<usize> = chr_snps_all
                .iter()
                .map(|s| bed_idx_to_pos[&s.bed_idx])
                .collect();

            let chr_snps: Vec<&BimRecord> = chr_snps_all
                .iter()
                .copied()
                .filter(|s| should_output(s))
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

            Ok(format!(
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
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    for msg in chr_messages {
        println!("{}", msg);
    }

    if verbose_timing {
        eprintln!(
            "[perf] write_outputs={:.3}s total={:.3}s",
            t_write_start.elapsed().as_secs_f64(),
            t_run_start.elapsed().as_secs_f64(),
        );
    }
    Ok(())
}
