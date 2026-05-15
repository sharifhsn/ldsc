mod compute;
mod io;
mod normalize;
mod snp_stats;
mod window;

pub use crate::parse::{BimRecord, parse_bim, parse_bim_reader, parse_bim_str};
pub use io::{count_fam, count_fam_str, parse_fam, parse_fam_str};
pub use window::WindowMode;

/// Per-chunk progress event emitted from inside the chunked GEMM
/// loop in [`compute_ldscore_global`]. The browser frontend wires
/// this to a `postMessage` so the main-thread UI can show a progress
/// bar; native callers typically pass `|_| {}` and ignore it.
///
/// Note: this is **progress only** (monotonic SNP-coverage counter).
/// L2 values for SNPs in the just-finished chunk are *not* finalized
/// yet — later chunks contribute cross-window r². The full L2 array
/// is returned at the end via [`L2Output`].
#[derive(Debug, Clone, Copy)]
pub struct L2Progress {
    /// Number of chunks processed so far (1-based; 1 after the first
    /// chunk completes, equals `chunks_total` after the last).
    pub chunks_done: usize,
    /// Total number of chunks the chunked-GEMM loop will run for
    /// this call to `compute_ldscore_global`.
    pub chunks_total: usize,
    /// Number of SNPs the loop has covered so far (= chunk_end of
    /// the most-recent chunk).
    pub snps_done: usize,
    /// Total SNPs in this call (= length of the input `all_snps`).
    pub snps_total: usize,
    /// Wall time consumed by *this chunk only* (decode + GEMM +
    /// ring store), in milliseconds. The orchestrator can sum these
    /// across workers for a real-time throughput readout.
    pub chunk_wall_ms: f64,
}

use crate::bed::Bed;
use crate::cli::L2Args;
use crate::cts_annot;
#[cfg(feature = "gpu")]
use crate::gpu::{GpuConfig, GpuContext};
use crate::la::MatF;
use crate::parse;
use anyhow::{Context, Result};
use compute::compute_ldscore_global;
use io::{
    format_m_vals, load_individual_indices, load_print_snps, load_snp_set, write_annot_matrix,
    write_ldscore_refs,
};
use rayon::prelude::*;
use snp_stats::compute_snp_stats;
use std::collections::HashSet;

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
    let t_run_start = web_time::Instant::now();
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

    // --python-compat overrides defaults for bit-identical Python LDSC parity.
    // Effects: chunk-size 50 (Python's default vs ldsc-rs's 200), global-pass
    // (single sequential pass matching Python's cross-chr window bleeding), no
    // SNP-level masking (Python doesn't do it). clap already enforces conflict
    // with --snp-level-masking.
    let chunk_c = if args.python_compat {
        50
    } else {
        args.chunk_size
    };
    let force_python_global_pass = args.python_compat && !args.global_pass;
    if args.python_compat {
        if args.chunk_size != 200 {
            eprintln!(
                "--python-compat overrides --chunk-size {} with c=50 for Python LDSC parity",
                args.chunk_size
            );
        }
        if !args.global_pass {
            eprintln!("--python-compat auto-enables --global-pass for Python LDSC parity");
        }
    }

    let t_maf_pre = t_run_start.elapsed();
    if verbose_timing {
        eprintln!("[perf] maf_prefilter={:.3}s", t_maf_pre.as_secs_f64());
    }

    // --sketch automatically enables f32: CountSketch ±1 entries are exactly
    // representable in f32, so f64 is strictly dominated (same accuracy, ~1.3× slower).
    let use_f32 = args.fast_f32 || args.sketch.is_some();
    if args.sketch.is_some() && !args.fast_f32 {
        println!("  --sketch auto-enables f32 (bit-identical to f64, ~1.3× faster)");
    }

    // GPU + per-chr parallel: each rayon thread allocates full ring/b_mat/a_buf buffers
    // plus competes for the single GPU. At large N this OOMs and serializes GPU work.
    // Auto-enable --global-pass when --gpu is active.
    #[cfg(feature = "gpu")]
    if args.gpu && !args.global_pass {
        eprintln!("GPU: auto-enabling --global-pass (per-chr parallel wastes GPU with contention)");
    }
    #[cfg(feature = "gpu")]
    let force_global_pass = (args.gpu && !args.global_pass) || force_python_global_pass;
    #[cfg(not(feature = "gpu"))]
    let force_global_pass = force_python_global_pass;

    // Create GPU context once (shared across all chromosomes).
    #[cfg(feature = "gpu")]
    let gpu_ctx = if args.gpu {
        match GpuContext::new(args.gpu_flex32) {
            Ok(ctx) => {
                if args.gpu_f64 && !ctx.capabilities.has_f64 {
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
    #[cfg(feature = "gpu")]
    let gpu_config = GpuConfig {
        tile_cols: args.gpu_tile_cols,
        flex32: args.gpu_flex32,
        f64: args.gpu_f64,
    };

    let t_compute_start = web_time::Instant::now();
    let (l2, maf_per_snp) = if args.global_pass || force_global_pass {
        // Legacy single-pass across all chromosomes.
        let bed = Bed::builder(bed_path.as_str())
            .build()
            .context("opening BED file")?;
        let mmap_bed = if args.mmap {
            Some(crate::bed::MmapBed::open(bed_path.as_str()).context("memory-mapping BED file")?)
        } else {
            None
        };
        compute_ldscore_global(
            &all_snps,
            bed,
            mmap_bed,
            n_indiv_actual,
            mode,
            chunk_c,
            annot_ref,
            iid_indices.as_deref(),
            pq_exp_for_compute,
            args.yes_really,
            #[cfg(feature = "gpu")]
            gpu_ctx.as_ref(),
            #[cfg(feature = "gpu")]
            gpu_config,
            use_f32,
            args.verbose_timing,
            args.sketch,
            args.sketch_maf_aware,
            args.snp_level_masking,
            |_| {}, // CLI: no per-chunk progress (terminal stdout would be too noisy)
        )
        .context("computing LD scores")?
    } else {
        // Per-chromosome parallel processing: split SNPs by chromosome,
        // compute each independently via rayon, then concatenate results.
        // Exploits diminishing GEMM thread scaling for ~25% improvement.
        let chr_groups: Vec<(u8, usize, usize)> = {
            let mut groups = Vec::new();
            let mut start = 0usize;
            while start < all_snps.len() {
                let chr = all_snps[start].chr;
                let mut end = start + 1;
                while end < all_snps.len() && all_snps[end].chr == chr {
                    end += 1;
                }
                groups.push((chr, start, end));
                start = end;
            }
            groups
        };
        println!(
            "  Per-chromosome parallel: {} chromosomes (use --global-pass for legacy mode)",
            chr_groups.len()
        );

        let bed_path_str = bed_path.as_str();
        let iid_slice = iid_indices.as_deref();
        let results: Vec<(MatF, Vec<f64>)> = chr_groups
            .par_iter()
            .map(|&(_chr, start, end)| {
                let chr_snps = &all_snps[start..end];
                let chr_annot: Option<MatF> = annot_ref.map(|a| {
                    let mut sub = MatF::zeros(end - start, a.ncols());
                    for i in 0..end - start {
                        for j in 0..a.ncols() {
                            sub[(i, j)] = a[(start + i, j)];
                        }
                    }
                    sub
                });
                // Each chr opens its own Bed (and optional MmapBed) — they
                // run concurrently from rayon tasks so they can't share a
                // single file handle.
                let bed = Bed::builder(bed_path_str)
                    .build()
                    .context("opening BED file")?;
                let mmap_bed = if args.mmap {
                    Some(
                        crate::bed::MmapBed::open(bed_path_str)
                            .context("memory-mapping BED file")?,
                    )
                } else {
                    None
                };
                compute_ldscore_global(
                    chr_snps,
                    bed,
                    mmap_bed,
                    n_indiv_actual,
                    mode,
                    chunk_c,
                    chr_annot.as_ref(),
                    iid_slice,
                    pq_exp_for_compute,
                    args.yes_really,
                    #[cfg(feature = "gpu")]
                    gpu_ctx.as_ref(),
                    #[cfg(feature = "gpu")]
                    gpu_config,
                    use_f32,
                    false, // verbose_timing disabled per-chr (noisy)
                    args.sketch,
                    args.sketch_maf_aware,
                    args.snp_level_masking,
                    |_| {}, // CLI: no per-chr progress (rayon already prints chr msgs)
                )
            })
            .collect::<Result<Vec<_>>>()
            .context("computing LD scores (per-chromosome)")?;

        // Concatenate per-chromosome results into global arrays.
        let n_annot = annot_ref.map(|a| a.ncols()).unwrap_or(1);
        let mut l2 = MatF::zeros(all_snps.len(), n_annot);
        let mut maf = vec![0.0f64; all_snps.len()];
        for (&(_chr, start, _end), (chr_l2, chr_maf)) in chr_groups.iter().zip(results.iter()) {
            for i in 0..chr_l2.nrows() {
                for k in 0..n_annot {
                    l2[(start + i, k)] = chr_l2[(i, k)];
                }
                maf[start + i] = chr_maf[i];
            }
        }
        (l2, maf)
    };

    if verbose_timing {
        eprintln!(
            "[perf] compute_ldscore_total={:.3}s",
            t_compute_start.elapsed().as_secs_f64()
        );
    }

    let t_write_start = web_time::Instant::now();
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

    let chr_messages: Vec<String> = chrs
        .par_iter()
        .map(|chr| -> Result<String> {
            let chr_snps_all: Vec<&BimRecord> = all_snps.iter().filter(|s| s.chr == *chr).collect();
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

    // Print LD Score summary (matches Python's "Summary of LD Scores in ..." output).
    // LDlink's Flask code parses for this string.
    print_ld_score_summary(
        &out_path,
        &out_snps,
        &out_l2,
        &col_names,
        &maf_per_snp,
        &out_positions,
    );

    if verbose_timing {
        eprintln!(
            "[perf] write_outputs={:.3}s total={:.3}s",
            t_write_start.elapsed().as_secs_f64(),
            t_run_start.elapsed().as_secs_f64(),
        );
    }
    Ok(())
}

/// Print LD Score summary statistics matching Python's output format.
/// LDlink parses stdout for "Summary of LD Scores in ..." to extract results.
fn print_ld_score_summary(
    out_path: &str,
    out_snps: &[&BimRecord],
    out_l2: &MatF,
    col_names: &[String],
    maf_per_snp: &[f64],
    out_positions: &[usize],
) {
    let n = out_snps.len();
    if n == 0 {
        return;
    }

    // Collect MAF values for output SNPs
    let maf_vals: Vec<f64> = out_positions.iter().map(|&pos| maf_per_snp[pos]).collect();

    // Build columns: MAF + each L2 column
    let mut all_cols: Vec<(&str, Vec<f64>)> = Vec::new();
    all_cols.push(("MAF", maf_vals));
    for (k, name) in col_names.iter().enumerate() {
        let vals: Vec<f64> = (0..n).map(|i| out_l2[(i, k)]).collect();
        all_cols.push((name, vals));
    }

    println!("\nSummary of LD Scores in {}", out_path);

    // Print header
    let col_width = 8;
    print!("{:>8}", "");
    for (name, _) in &all_cols {
        print!("{:>col_width$}", name);
    }
    println!();

    // Compute and print descriptive stats: mean, std, min, 25%, 50%, 75%, max
    let labels = ["mean", "std", "min", "25%", "50%", "75%", "max"];
    for label in &labels {
        print!("{:>8}", label);
        for (_, vals) in &all_cols {
            let stat = compute_descriptive_stat(vals, label);
            print!("{:>col_width$.4}", stat);
        }
        println!();
    }

    // Print MAF/LD Score Correlation Matrix
    println!();
    println!("MAF/LD Score Correlation Matrix");
    let col_width_corr = 8;
    print!("{:>8}", "");
    for (name, _) in &all_cols {
        print!("{:>col_width_corr$}", name);
    }
    println!();

    for (i, (name_i, vals_i)) in all_cols.iter().enumerate() {
        print!("{:>8}", name_i);
        for (j, (_, vals_j)) in all_cols.iter().enumerate() {
            if i == j {
                print!("{:>col_width_corr$.4}", 1.0);
            } else {
                print!("{:>col_width_corr$.4}", pearson_corr(vals_i, vals_j));
            }
        }
        println!();
    }
}

fn compute_descriptive_stat(vals: &[f64], stat: &str) -> f64 {
    let n = vals.len();
    if n == 0 {
        return f64::NAN;
    }
    match stat {
        "mean" => vals.iter().sum::<f64>() / n as f64,
        "std" => {
            let mean = vals.iter().sum::<f64>() / n as f64;
            let var = vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1).max(1) as f64;
            var.sqrt()
        }
        "min" => vals.iter().copied().fold(f64::INFINITY, f64::min),
        "max" => vals.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        pct => {
            let q = match pct {
                "25%" => 0.25,
                "50%" => 0.50,
                "75%" => 0.75,
                _ => return f64::NAN,
            };
            percentile(vals, q)
        }
    }
}

fn percentile(vals: &[f64], q: f64) -> f64 {
    let mut sorted: Vec<f64> = vals.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n == 1 {
        return sorted[0];
    }
    // Use linear interpolation (same as pandas default)
    let pos = q * (n - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = lo + 1;
    let frac = pos - lo as f64;
    if hi >= n {
        sorted[lo]
    } else {
        sorted[lo] * (1.0 - frac) + sorted[hi] * frac
    }
}

fn pearson_corr(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let mut cov = 0.0;
    let mut vx = 0.0;
    let mut vy = 0.0;
    for i in 0..x.len() {
        let dx = x[i] - mx;
        let dy = y[i] - my;
        cov += dx * dy;
        vx += dx * dx;
        vy += dy * dy;
    }
    if vx == 0.0 || vy == 0.0 {
        return f64::NAN;
    }
    cov / (vx * vy).sqrt()
}

// ── Public compute API for ldsc-web / library callers ─────────────────────────
//
// The `run(L2Args)` orchestrator above does end-to-end disk I/O —
// `--bfile` paths, `--annot` / `--extract` / `--keep` / `--print-snps`
// file loaders, `.l2.ldscore.gz` writers, the Python-style stdout
// summary. It exists because every CLI test fixture and production
// pipeline calls it.
//
// The MVP web demo doesn't have a filesystem, so it can't go through
// `run`. The slim `compute_l2_from_bytes` API below accepts BED / BIM
// / FAM contents already in memory, exposes the headline LD-score
// computation knobs, and returns the L2 / MAF arrays for the caller
// to plot. Annotation / extract / keep / per-allele aren't wired in
// yet — they'll come back as additional fields on `L2Config` when
// the UI grows past v1.

/// Configuration knobs for [`compute_l2_from_bytes`]. Mirrors the
/// subset of `L2Args` that the in-browser demo exposes.
#[derive(Debug, Clone)]
pub struct L2Config {
    /// Window definition (`--ld-wind-kb` / `--ld-wind-snps` / `--ld-wind-cm`).
    pub mode: WindowMode,
    /// `--chunk-size`. Defaults to 200 — Python LDSC is 50.
    pub chunk_size: usize,
    /// `--fast-f32`.
    pub use_f32: bool,
    /// `--sketch d`. `None` runs exact.
    pub sketch: Option<usize>,
    /// `--sketch-maf-aware` (only meaningful when `sketch` is `Some`).
    pub sketch_maf_aware: bool,
    /// `--snp-level-masking`.
    pub snp_level_masking: bool,
    /// `--yes-really` — bypass the "your window covers a whole
    /// chromosome" safety check. In the browser we always have small
    /// inputs by design, so this is fine to default to `true`.
    pub yes_really: bool,
    /// `--pq-exp` (or `Some(1.0)` for `--per-allele`).
    pub pq_exp: Option<f64>,
    /// `--verbose-timing`. Off in the browser; we render our own timing.
    pub verbose_timing: bool,
}

impl Default for L2Config {
    fn default() -> Self {
        Self {
            mode: WindowMode::Kb(1000.0),
            chunk_size: 200,
            use_f32: false,
            sketch: None,
            sketch_maf_aware: false,
            snp_level_masking: false,
            yes_really: true,
            pq_exp: None,
            verbose_timing: false,
        }
    }
}

/// Result of [`compute_l2_from_bytes`]. `l2`, `maf` and `snps` are
/// aligned: `l2[i]` and `maf[i]` correspond to `snps[i]`.
#[derive(Debug, Clone)]
pub struct L2Output {
    pub snps: Vec<BimRecord>,
    pub l2: Vec<f64>,
    pub maf: Vec<f64>,
    pub wall_seconds: f64,
}

/// Compute LD scores from in-memory BED / BIM / FAM contents.
///
/// This is the browser-facing entry point. It performs no I/O: BIM
/// and FAM are passed as `&str` (what `FileReader.readAsText()` hands
/// the Leptos app), BED as a raw byte buffer (`FileReader.readAsArrayBuffer()`).
///
/// The current implementation drives the single-annotation
/// (`K = 1`) path with no `--extract` / `--keep` / `--annot`
/// filtering — that's the headline demo. Annotation support will land
/// as additional `L2Config` fields once the UI grows past v1.
pub fn compute_l2_from_bytes(
    bed_bytes: Vec<u8>,
    bim_text: &str,
    fam_text: &str,
    config: L2Config,
    on_progress: impl FnMut(L2Progress),
) -> Result<L2Output> {
    let snps = crate::parse::parse_bim_str(bim_text)
        .context("parsing BIM text in compute_l2_from_bytes")?;
    anyhow::ensure!(!snps.is_empty(), "compute_l2_from_bytes: BIM has zero SNPs");

    let n_indiv = io::count_fam_str(fam_text);
    anyhow::ensure!(
        n_indiv > 0,
        "compute_l2_from_bytes: FAM has zero individuals"
    );

    let bed = Bed::from_bytes(bed_bytes, n_indiv, snps.len())
        .context("validating BED bytes in compute_l2_from_bytes")?;

    compute_l2_from_bed(bed, snps, n_indiv, config, on_progress)
}

/// Compute LD scores from a pre-opened [`Bed`] + parsed BIM + FAM
/// individual count.
///
/// This is the API ldsc-web's Web Worker calls — the Worker
/// constructs a [`Bed`] backed by a [`crate::bed::BedSource`] that
/// streams chunks out of a JS `Blob` via `FileReaderSync`, so the
/// browser never has to load the whole multi-GB BED into wasm
/// linear memory. Native callers can use the same API with any
/// `Bed` (e.g. `Bed::builder(path).build()`).
pub fn compute_l2_from_bed(
    bed: Bed,
    snps: Vec<BimRecord>,
    n_indiv: usize,
    config: L2Config,
    on_progress: impl FnMut(L2Progress),
) -> Result<L2Output> {
    anyhow::ensure!(!snps.is_empty(), "compute_l2_from_bed: BIM has zero SNPs");
    anyhow::ensure!(n_indiv > 0, "compute_l2_from_bed: FAM has zero individuals");

    let t0 = web_time::Instant::now();
    let (l2_mat, maf) = compute_ldscore_global(
        &snps,
        bed,
        None, // mmap is unavailable / pointless in-browser
        n_indiv,
        config.mode,
        config.chunk_size,
        None, // no --annot in MVP
        None, // no --keep in MVP
        config.pq_exp,
        config.yes_really,
        #[cfg(feature = "gpu")]
        None,
        #[cfg(feature = "gpu")]
        crate::gpu::GpuConfig::default(),
        config.use_f32,
        config.verbose_timing,
        config.sketch,
        config.sketch_maf_aware,
        config.snp_level_masking,
        on_progress,
    )
    .context("computing LD scores in compute_l2_from_bed")?;

    // K = 1: collapse the (n × 1) matrix into a flat Vec<f64>.
    debug_assert_eq!(l2_mat.ncols(), 1);
    let l2: Vec<f64> = (0..l2_mat.nrows()).map(|i| l2_mat[(i, 0)]).collect();

    Ok(L2Output {
        snps,
        l2,
        maf,
        wall_seconds: t0.elapsed().as_secs_f64(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Hand-crafted micro BED: 4 individuals × 2 perfectly-correlated SNPs.
    ///
    /// Genotype counts (count-A1 semantics, what ldsc uses):
    ///   SNP1:  [2, 2, 0, 0]
    ///   SNP2:  [2, 2, 0, 0]
    /// PLINK packs each genotype as 2 bits within a byte, iid 0 in the
    /// least-significant pair. With the LUT `[2, missing, 1, 0]`:
    ///   count=2 → bits `00`
    ///   count=0 → bits `11`
    /// So per SNP the byte is `0b11_11_00_00 = 0xF0`.
    ///
    /// Pure-Rust expectation: SNPs are identical after centering →
    /// r² = 1 between SNP1 and SNP2. With unbiased correction
    /// `r² − (1 − r²)/(N − 2)`, contribution = 1.0 either way. Plus
    /// every SNP includes its self-LD of 1.0. So `l2 = [2.0, 2.0]`
    /// up to float precision, and `maf = [0.5, 0.5]`.
    #[test]
    fn compute_l2_from_bytes_smoke() {
        let mut bed_bytes = vec![0x6c, 0x1b, 0x01];
        bed_bytes.push(0xF0); // SNP 1 slab
        bed_bytes.push(0xF0); // SNP 2 slab

        // BIM: chr 1, SNP id, CM 0, BP 100/200, A1 A2 (tab-separated)
        let bim_text = "1\trs1\t0\t100\tA\tG\n1\trs2\t0\t200\tA\tG\n";
        // FAM: 4 individuals (FID IID PID MID SEX PHEN)
        let fam_text =
            "F1\tI1\t0\t0\t1\t-9\nF2\tI2\t0\t0\t1\t-9\nF3\tI3\t0\t0\t1\t-9\nF4\tI4\t0\t0\t1\t-9\n";

        let cfg = L2Config {
            // 100 kb window — both SNPs land in each other's window (200 - 100 = 100 bp ≤ 100 kb)
            mode: WindowMode::Kb(100.0),
            chunk_size: 2,
            yes_really: true,
            ..L2Config::default()
        };
        let out = compute_l2_from_bytes(bed_bytes, bim_text, fam_text, cfg, |_| {})
            .expect("compute_l2_from_bytes must accept a valid micro BED");

        assert_eq!(out.snps.len(), 2);
        assert_eq!(out.l2.len(), 2);
        assert_eq!(out.maf.len(), 2);
        for &maf in &out.maf {
            assert!((maf - 0.5).abs() < 1e-9, "MAF should be 0.5, got {}", maf);
        }
        for (i, &l2) in out.l2.iter().enumerate() {
            assert!(
                (l2 - 2.0).abs() < 1e-9,
                "SNP {i} LD score should be 2.0 (self + perfect-LD partner), got {l2}"
            );
        }
    }

    /// Configuration mismatches (here: a too-short BED) must produce a
    /// graceful `Result::Err`, not a panic.
    #[test]
    fn compute_l2_from_bytes_rejects_truncated_bed() {
        let bed_bytes = vec![0x6c, 0x1b, 0x01, 0xF0]; // header + only 1 SNP slab
        let bim_text = "1\trs1\t0\t100\tA\tG\n1\trs2\t0\t200\tA\tG\n";
        let fam_text =
            "F1\tI1\t0\t0\t1\t-9\nF2\tI2\t0\t0\t1\t-9\nF3\tI3\t0\t0\t1\t-9\nF4\tI4\t0\t0\t1\t-9\n";

        let cfg = L2Config {
            yes_really: true,
            ..L2Config::default()
        };
        let result = compute_l2_from_bytes(bed_bytes, bim_text, fam_text, cfg, |_| {});
        assert!(
            result.is_err(),
            "truncated BED must be rejected; got {:?}",
            result.as_ref().map(|o| o.l2.len())
        );
    }
}
