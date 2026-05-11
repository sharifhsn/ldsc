//! Bit-packed exact LD score computation (Phase 1 — scalar reference).
//!
//! Compute LD scores directly from packed BED bytes without materializing the
//! N×c f32 genotype matrix or running a dense GEMM. For each SNP pair (j, k) in
//! a per-SNP window, decode genotypes inline, accumulate the centered dot
//! product, and convert to r². This produces *exact* per-SNP windows (matching
//! the LDSC paper's mathematical definition), bit-stable vs the existing
//! `--snp-level-masking` f64 path up to float ULPs.
//!
//! Phase 1 is a scalar reference for correctness validation; wall clock is
//! expected to be 2–4× slower than the current exact path. Phase 2 will add
//! AVX2 byte-level SIMD; Phase 3 will benchmark and decide whether to promote
//! `--exact-bitpacked` as the default exact mode.
//!
//! Design choices:
//! - Iterate per-SNP windows `k in block_left[j]..j` strictly — the iteration
//!   IS the per-SNP mask. No post-GEMM masking needed.
//! - Mean imputation: missing genotype contributes 0 to the centered dot
//!   product (matches the f32 path's NaN→0 normalization).
//! - Parallelism comes from rayon over chromosomes in `super::run`; the inner
//!   loop here is single-threaded to avoid races on `l2[k, c] += …`.
//! - Memory access: prefer `MmapBed::snp_bytes` for zero-copy; fall back to a
//!   `Bed`-backed read into an owned buffer. The plan recommends `--mmap`.

use super::compute::{build_bed_byte_lut, snp_stats_from_bytes};
use super::window::{WindowMode, get_block_lefts_by_chr};
use crate::bed::{Bed, IidPos, MmapBed, precompute_iid_positions, resolve_indices};
use crate::la::MatF;
use crate::parse::BimRecord;
use anyhow::{Context, Result};

/// Per-SNP statistics needed for the bit-packed inner product.
///
/// `mean` and `inv_std` come from the f64-precision (sum, count, sum_sq) over
/// the selected individuals; `inv_std` divides by `n_indiv_sel` (the kept
/// individual count), matching the existing `normalize_col_f64_with_stats`
/// convention. `maf` is the minor allele frequency for the output and pq
/// weighting.
#[derive(Clone, Copy)]
struct SnpStatsBP {
    mean: f64,
    inv_std: f64,
    maf: f64,
}

/// Centered dot product `Σ_i (g_j[i] − μ_j)(g_k[i] − μ_k)` over individuals
/// where both genotypes are non-missing (matches NaN→0 mean imputation).
///
/// Uses a per-pair 4×4 product table indexed by 2-bit genotype codes
/// (0b00=hom-A1=2, 0b01=missing, 0b10=het=1, 0b11=hom-A2=0 under count_a1).
/// Scalar reference — no SIMD, ~4 decode+lookups per byte-pair.
#[inline]
fn pair_centered_dot(
    bytes_j: &[u8],
    bytes_k: &[u8],
    n_indiv_sel: usize,
    iid_positions: &[IidPos],
    all_iids: bool,
    mean_j: f64,
    mean_k: f64,
) -> f64 {
    // Decoded genotype value per 2-bit code; None = missing.
    let g_dec: [Option<f64>; 4] = [
        Some(2.0), // 0b00 = hom A1 = 2  (count_a1=true convention)
        None,      // 0b01 = missing
        Some(1.0), // 0b10 = het = 1
        Some(0.0), // 0b11 = hom A2 = 0
    ];
    // Per-pair table: prod[code_j][code_k] = (g_j − μ_j)(g_k − μ_k) if both
    // non-missing, else 0. 16 entries — cheap to materialize per pair, cached
    // in a register tile during the byte walk.
    let mut prod = [[0.0f64; 4]; 4];
    for (cj, gj_opt) in g_dec.iter().enumerate() {
        for (ck, gk_opt) in g_dec.iter().enumerate() {
            prod[cj][ck] = match (gj_opt, gk_opt) {
                (Some(gj), Some(gk)) => (gj - mean_j) * (gk - mean_k),
                _ => 0.0,
            };
        }
    }

    let mut acc = 0.0f64;
    if all_iids {
        // Full coverage: walk 4 individuals per byte pair.
        let full_bytes = n_indiv_sel / 4;
        let rem = n_indiv_sel % 4;
        for b in 0..full_bytes {
            let bj = bytes_j[b];
            let bk = bytes_k[b];
            acc += prod[(bj & 0b11) as usize][(bk & 0b11) as usize];
            acc += prod[((bj >> 2) & 0b11) as usize][((bk >> 2) & 0b11) as usize];
            acc += prod[((bj >> 4) & 0b11) as usize][((bk >> 4) & 0b11) as usize];
            acc += prod[((bj >> 6) & 0b11) as usize][((bk >> 6) & 0b11) as usize];
        }
        if rem > 0 {
            let bj = bytes_j[full_bytes];
            let bk = bytes_k[full_bytes];
            for r in 0..rem {
                let cj = (bj >> (r * 2)) & 0b11;
                let ck = (bk >> (r * 2)) & 0b11;
                acc += prod[cj as usize][ck as usize];
            }
        }
    } else {
        // Subset of individuals: gather via precomputed (byte_idx, shift) pairs.
        for pos in iid_positions {
            let bj = bytes_j[pos.byte_idx];
            let bk = bytes_k[pos.byte_idx];
            let cj = (bj >> pos.shift) & 0b11;
            let ck = (bk >> pos.shift) & 0b11;
            acc += prod[cj as usize][ck as usize];
        }
    }
    acc
}

/// Compute LD scores for all SNPs via direct iteration over per-SNP windows
/// from packed BED bytes. Returns `(l2, maf_per_snp)`.
///
/// Semantics match `--snp-level-masking` (per-SNP exact windows, the LDSC
/// paper's `ℓ_j = Σ_k r²_{jk}`), not Python LDSC's chunk-level approximation.
#[allow(clippy::too_many_arguments)]
pub(super) fn compute_ldscore_bitpacked(
    all_snps: &[BimRecord],
    bed_path: &str,
    n_indiv_sel: usize,
    mode: WindowMode,
    annot: Option<&MatF>,
    iid_indices: Option<&[isize]>,
    pq_exp: Option<f64>,
    yes_really: bool,
    use_mmap: bool,
) -> Result<(MatF, Vec<f64>)> {
    let m = all_snps.len();
    let n_annot = annot.map(|a| a.ncols()).unwrap_or(1);
    if m == 0 {
        return Ok((MatF::zeros(0, n_annot), vec![]));
    }

    if let Some(a) = annot {
        anyhow::ensure!(
            a.nrows() == m,
            "Annotation matrix has {} rows but BIM has {} SNPs",
            a.nrows(),
            m
        );
    }

    // Per-chromosome window boundaries — strict per-SNP windows enforced by
    // the iteration `k in block_left[j]..j`.
    let (block_left, any_full_chr) = get_block_lefts_by_chr(all_snps, mode);
    if any_full_chr && !yes_really {
        anyhow::bail!(
            "LD window spans an entire chromosome. Pass --yes-really to confirm this is intended."
        );
    }

    // Resolve iid_indices into (byte_idx, shift) positions. Always create the
    // ChunkReader-style iid_positions even when all_iids=true, so a single
    // code path handles both cases. The total FAM count comes from MmapBed or
    // Bed; we open one of them and derive bytes_per_snp / iid_count from it.
    //
    // We always open `Bed` (cheap header check + bytes_per_snp) regardless of
    // --mmap, since MmapBed::open also reads the FAM/BIM files (slow on the
    // synthetic biobank dataset). When --mmap is on we additionally hold a
    // memory map for zero-copy reads; otherwise we use the Bed-backed reader.
    let mut bed = Bed::builder(bed_path)
        .build()
        .with_context(|| format!("opening BED '{}'", bed_path))?;
    let n_indiv_fam = bed.iid_count();
    let bytes_per_snp = bed.bytes_per_snp();
    let iid_idx = resolve_indices(iid_indices, n_indiv_fam)?;
    let n_sel_from_idx = iid_idx.len();
    anyhow::ensure!(
        n_sel_from_idx == n_indiv_sel,
        "internal: selected iid count {} disagrees with iid_indices length {}",
        n_indiv_sel,
        n_sel_from_idx,
    );
    let all_iids = n_indiv_sel == n_indiv_fam
        && iid_idx.iter().enumerate().all(|(i, &v)| v == i);
    let iid_positions = precompute_iid_positions(&iid_idx);

    // Optional memory map for zero-copy byte access — held alive for both
    // the stats pass and the inner-loop scatter, since dropping mid-call
    // would force a Bed-backed fallback.
    let mmap_bed = if use_mmap {
        Some(MmapBed::open(bed_path).with_context(|| format!("mmap'ing BED '{}'", bed_path))?)
    } else {
        None
    };

    // -------- Pass 1: per-SNP statistics ---------------------------------
    //
    // We need (mean, inv_std, maf) for every SNP to compute centered dot
    // products in the inner loop. Read each SNP once; mmap path is zero-copy,
    // Bed-backed path allocates a single bytes_per_snp buffer reused across
    // SNPs.
    let stats: Vec<SnpStatsBP> = {
        let byte_lut = build_bed_byte_lut();
        let mut out = Vec::with_capacity(m);
        let mut buf: Vec<u8> = vec![0u8; bytes_per_snp];
        for snp in all_snps {
            let bytes: &[u8] = if let Some(ref mb) = mmap_bed {
                mb.snp_bytes(snp.bed_idx, 1)
            } else {
                bed.read_snp_bytes(snp.bed_idx, &mut buf)
                    .with_context(|| format!("reading BED SNP {}", snp.bed_idx))?;
                &buf
            };
            let (sum, count, sum_sq) =
                snp_stats_from_bytes(bytes, n_indiv_sel, &iid_positions, all_iids, &byte_lut);
            let mean = if count > 0 {
                sum as f64 / count as f64
            } else {
                0.0
            };
            let freq = (mean / 2.0).clamp(0.0, 1.0);
            let maf = freq.min(1.0 - freq);
            let centered_ss = sum_sq as f64 - count as f64 * mean * mean;
            let var = centered_ss / n_indiv_sel as f64;
            let std_ = var.sqrt();
            let inv_std = if std_ > 0.0 { 1.0 / std_ } else { 0.0 };
            out.push(SnpStatsBP {
                mean,
                inv_std,
                maf,
            });
        }
        out
    };

    let maf_per_snp: Vec<f64> = stats.iter().map(|s| s.maf).collect();

    // Per-SNP pq weight: (maf * (1-maf))^pq_exp, or 1.0 when pq_exp = None.
    let pq: Vec<f64> = match pq_exp {
        Some(exp) => stats
            .iter()
            .map(|s| (s.maf * (1.0 - s.maf)).powf(exp))
            .collect(),
        None => vec![1.0; m],
    };

    // -------- Pass 2: pair iteration -------------------------------------
    //
    // For each j ∈ [0, m), iterate k ∈ [block_left[j], j) and accumulate
    //   l2[j, c] += r²(j, k) × annot_eff[k, c]
    //   l2[k, c] += r²(j, k) × annot_eff[j, c]
    // Diagonal (j == k): r² = 1, l2[j, c] += annot_eff[j, c].
    // where annot_eff[i, c] = (annot[i, c] if Some, else 1.0) × pq[i].
    //
    // r²_unbiased = r²_raw × r2u_a + r2u_b with
    //   r2u_a = 1 + 1/(N−2), r2u_b = −1/(N−2)
    // i.e. r²_unbiased = r²_raw − (1 − r²_raw)/(N − 2).
    let n_f = n_indiv_sel as f64;
    let n_inv = 1.0 / n_f;
    let (r2u_a, r2u_b) = {
        let denom = if n_indiv_sel > 2 {
            n_f - 2.0
        } else {
            n_f.max(1.0)
        };
        (1.0 + 1.0 / denom, -1.0 / denom)
    };

    let mut l2 = MatF::zeros(m, n_annot);

    // Pre-compute annot_eff[j, c] = annot[j, c] × pq[j] (or pq[j] when annot is None).
    // For n_annot == 1 we keep the bulk loop tight by using a flat Vec; for
    // n_annot > 1 we accumulate per-column inside the hot loop.
    let annot_eff: Vec<f64> = match annot {
        Some(a) => {
            let mut v = Vec::with_capacity(m * n_annot);
            for j in 0..m {
                for c in 0..n_annot {
                    v.push(a[(j, c)] * pq[j]);
                }
            }
            v
        }
        None => pq.clone(), // n_annot == 1, annot_eff[j, 0] = pq[j]
    };
    let stride = n_annot;

    // Read-byte helper: zero-copy via mmap when available, else pread into the
    // scratch buffer. We can't use a closure because it would need to borrow
    // both `bed` mutably twice across the j and k fetches — instead we inline
    // the dispatch where needed and own the scratch buffers in the outer
    // scope. To avoid double-buffering issues with mmap (slices can't outlive
    // a re-fetch when the inner loop reads a new k), the mmap path *re-fetches*
    // bytes_j per inner iteration too. That's free (zero-copy) and keeps the
    // borrow checker happy.
    let mut bytes_j_buf: Vec<u8> = vec![0u8; bytes_per_snp];
    let mut bytes_k_buf: Vec<u8> = vec![0u8; bytes_per_snp];

    for j in 0..m {
        let sj = stats[j];
        // Diagonal contribution: r²(j, j) = 1 exactly. The unbiased
        // correction at r² = 1 yields 1 × r2u_a + r2u_b = (1 + 1/(N-2)) − 1/(N-2) = 1,
        // so we use 1.0 directly.
        for c in 0..n_annot {
            l2[(j, c)] += annot_eff[j * stride + c];
        }

        let bl = block_left[j];
        if bl >= j {
            // No left neighbors in window — only the diagonal contributes.
            continue;
        }

        // Fast path: SNP j has zero variance. The GEMM path produces r²_raw=0
        // for every pair (the normalized column is all zeros), so
        // r²_unbiased = r2u_b = −1/(N−2) for every k. Emit that contribution
        // without reading any bytes.
        if sj.inv_std == 0.0 {
            for k in bl..j {
                let r2u = r2u_b;
                for c in 0..n_annot {
                    l2[(j, c)] += r2u * annot_eff[k * stride + c];
                    l2[(k, c)] += r2u * annot_eff[j * stride + c];
                }
            }
            continue;
        }

        // Load bytes for SNP j (non-mmap path only — mmap path slices inline).
        let bed_idx_j = all_snps[j].bed_idx;
        if mmap_bed.is_none() {
            bed.read_snp_bytes(bed_idx_j, &mut bytes_j_buf)
                .with_context(|| format!("reading BED SNP {}", bed_idx_j))?;
        }

        // For each k in [bl, j), compute r²_unbiased and scatter into l2[j] and l2[k].
        for k in bl..j {
            let sk = stats[k];

            // Zero-variance SNP k: same shortcut as the sj zero-variance fast
            // path above. r²_raw = 0 → r²_unbiased = r2u_b.
            let r2u = if sk.inv_std == 0.0 {
                r2u_b
            } else {
                let bed_idx_k = all_snps[k].bed_idx;
                let (bytes_j_view, bytes_k_view): (&[u8], &[u8]) = if let Some(ref mb) = mmap_bed {
                    (mb.snp_bytes(bed_idx_j, 1), mb.snp_bytes(bed_idx_k, 1))
                } else {
                    bed.read_snp_bytes(bed_idx_k, &mut bytes_k_buf)
                        .with_context(|| format!("reading BED SNP {}", bed_idx_k))?;
                    (&bytes_j_buf[..], &bytes_k_buf[..])
                };

                let centered_dot = pair_centered_dot(
                    bytes_j_view,
                    bytes_k_view,
                    n_indiv_sel,
                    &iid_positions,
                    all_iids,
                    sj.mean,
                    sk.mean,
                );
                let r = centered_dot * sj.inv_std * sk.inv_std * n_inv;
                let r2_raw = r * r;
                // Clamp r² to [0, 1] before the LDSC unbiased correction.
                // r² > 1 can only arise from float-noise on near-zero-variance
                // SNPs (or mean-imputation interactions); the dense GEMM path
                // can produce the same artifacts. Clamping keeps the
                // correction sensible.
                let r2_raw = r2_raw.min(1.0);
                r2_raw * r2u_a + r2u_b
            };

            // Scatter into l2[j] and l2[k]. Single-thread inside one chr to
            // avoid races on l2[k] from concurrent j's; parallelism is over
            // chromosomes in the dispatcher.
            for c in 0..n_annot {
                l2[(j, c)] += r2u * annot_eff[k * stride + c];
                l2[(k, c)] += r2u * annot_eff[j * stride + c];
            }
        }
    }

    Ok((l2, maf_per_snp))
}
