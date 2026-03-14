mod compute;
mod io;
mod normalize;
mod snp_stats;
mod window;

pub use crate::parse::{BimRecord, parse_bim};
pub use io::{count_fam, parse_fam};
pub use window::WindowMode;

use crate::bed::Bed;
use crate::cli::L2Args;
use crate::cts_annot;
use crate::la::MatF;
use crate::parse;
use anyhow::{Context, Result};
use compute::compute_ldscore_global;
use io::{
    format_m_vals, load_individual_indices, load_print_snps, load_snp_set, write_annot_matrix,
    write_ldscore_refs,
};
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
    // --keep / --subsample: subset individuals for LD computation.
    let iid_indices: Option<Vec<isize>> = if let Some(ref keep_path) = args.keep {
        let fam_ids =
            parse_fam(&fam_path).with_context(|| format!("parsing FAM '{}'", fam_path))?;
        Some(
            load_individual_indices(keep_path, &fam_ids)
                .with_context(|| format!("loading keep file '{}'", keep_path))?,
        )
    } else if let Some(n_prime) = args.subsample {
        anyhow::ensure!(n_prime > 0, "--subsample must be > 0");
        if n_prime >= n_indiv {
            eprintln!(
                "WARNING: --subsample {} >= n_indiv {}; using all individuals",
                n_prime, n_indiv
            );
            None
        } else {
            // Fisher-Yates partial shuffle to select n_prime indices (deterministic seed 42).
            // Sort after selection so BED reads remain sequential-friendly.
            let mut rng = fastrand::Rng::with_seed(42);
            let mut indices: Vec<usize> = (0..n_indiv).collect();
            for i in 0..n_prime {
                let j = rng.usize(i..n_indiv);
                indices.swap(i, j);
            }
            indices.truncate(n_prime);
            indices.sort_unstable();
            println!(
                "  --subsample: retaining {}/{} individuals",
                n_prime, n_indiv
            );
            Some(indices.iter().map(|&i| i as isize).collect())
        }
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

    // --sketch automatically enables f32: sketch projections produce bit-identical
    // results in f32 vs f64 (Rademacher ±1/√d and CountSketch ±1 are exactly
    // representable), so f64 is strictly dominated (same accuracy, ~1.3× slower).
    let use_f32 = args.fast_f32 || args.sketch.is_some();
    if args.sketch.is_some() && !args.fast_f32 {
        println!("  --sketch auto-enables f32 (bit-identical to f64, ~1.3× faster)");
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
        use_f32,
        args.prefetch_bed,
        args.verbose_timing,
        args.stochastic,
        args.sketch,
        args.sketch_method.as_str(),
        args.mmap,
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

    if verbose_timing {
        eprintln!(
            "[perf] write_outputs={:.3}s total={:.3}s",
            t_write_start.elapsed().as_secs_f64(),
            t_run_start.elapsed().as_secs_f64(),
        );
    }
    Ok(())
}
