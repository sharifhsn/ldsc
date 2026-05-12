//! Bit-packed exact LD score computation.
//!
//! Compute LD scores directly from packed BED bytes without materializing the
//! N×c f32 genotype matrix or running a dense GEMM. For each SNP pair (j, k) in
//! a per-SNP window, decode genotypes inline, accumulate the centered dot
//! product, and convert to r². This produces *exact* per-SNP windows (matching
//! the LDSC paper's mathematical definition), bit-stable vs the existing
//! `--snp-level-masking` f64 path up to float ULPs.
//!
//! Inner kernel strategy (Phase 2): histogram + dot. For each (j, k) pair we
//! compute the joint genotype-code histogram (16 bins, indexed by
//! `(cj << 2) | ck` over 2-bit BED codes), then form the centered dot as a
//! 16-element dot with a per-pair `(g_cj − μ_j)(g_ck − μ_k)` product table.
//! Decoupling the (SIMD-friendly) count from the (f64-precise) product gives
//! exact arithmetic while letting AVX2 dominate the per-individual work.
//!
//! When AVX2 is available the histogram walks 32 BED bytes (= 128 individuals)
//! per iteration: split each byte into 4 positional code vectors via shift+mask,
//! combine into a 32-byte joint-code vector, then 16 cmpeq + movemask + popcnt
//! per joint vector tally each bin in parallel. Phase 1 scalar fallback kept
//! for non-x86_64 targets and the iid-subset path (`--keep`), where the
//! scatter-gather access pattern makes SIMD a wash.
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
use rayon::prelude::*;

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

/// Decoded genotype value per 2-bit BED code; `None` = missing.
const G_DEC: [Option<f64>; 4] = [
    Some(2.0), // 0b00 = hom A1 = 2  (count_a1=true convention)
    None,      // 0b01 = missing
    Some(1.0), // 0b10 = het = 1
    Some(0.0), // 0b11 = hom A2 = 0
];

/// Per-pair 16-entry product table indexed by `(cj << 2) | ck`.
/// `lut[code] = (g_cj − μ_j)(g_ck − μ_k)` if both codes are non-missing,
/// else 0 (matches NaN→0 mean imputation in the GEMM path).
#[inline]
fn build_prod_lut(mean_j: f64, mean_k: f64) -> [f64; 16] {
    let mut lut = [0.0f64; 16];
    for cj in 0..4 {
        for ck in 0..4 {
            lut[(cj << 2) | ck] = match (G_DEC[cj], G_DEC[ck]) {
                (Some(gj), Some(gk)) => (gj - mean_j) * (gk - mean_k),
                _ => 0.0,
            };
        }
    }
    lut
}

/// Centered dot product `Σ_i (g_j[i] − μ_j)(g_k[i] − μ_k)` over individuals
/// where both genotypes are non-missing.
///
/// Dispatches to the AVX2 histogram kernel on x86_64 for the full-coverage
/// case; scalar reference handles non-x86_64 targets and the iid-subset path.
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
    let prod = build_prod_lut(mean_j, mean_k);

    if !all_iids {
        // iid-subset path (--keep): scatter-gather access pattern, scalar.
        return pair_centered_dot_subset(bytes_j, bytes_k, iid_positions, &prod);
    }

    // Dispatch to the SIMD histogram kernel for the dense full-coverage path.
    // x86_64: AVX2 (32 bytes / 128 individuals per chunk).
    // aarch64: NEON (16 bytes / 64 individuals per chunk).
    // Other targets: scalar fallback.
    let counts: [u32; 16] = {
        #[cfg(target_arch = "x86_64")]
        {
            // SAFETY: built with `target-feature=+avx2` for the production
            // musl/glibc targets; native dev builds run on x86_64 CPUs that
            // have had AVX2 since 2013.
            unsafe { pair_joint_histogram_avx2(bytes_j, bytes_k, n_indiv_sel) }
        }
        #[cfg(target_arch = "aarch64")]
        {
            // SAFETY: NEON is a baseline aarch64 feature (Apple Silicon,
            // AWS Graviton, ARMv8+).
            unsafe { pair_joint_histogram_neon(bytes_j, bytes_k, n_indiv_sel) }
        }
        #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
        {
            return pair_centered_dot_scalar_all_iids(
                bytes_j,
                bytes_k,
                n_indiv_sel,
                &prod,
            );
        }
    };
    let mut acc = 0.0f64;
    for c in 0..16 {
        acc += counts[c] as f64 * prod[c];
    }
    acc
}

/// Scalar reference: walks 4 individuals per byte pair through the full
/// individual range. Used as a fallback on non-x86_64/aarch64 targets and
/// kept around for unit-test parity checks. The AVX2/NEON paths inline
/// their own scalar tail handlers for the residual bytes.
#[cfg_attr(
    any(target_arch = "x86_64", target_arch = "aarch64"),
    allow(dead_code)
)]
#[inline]
fn pair_centered_dot_scalar_all_iids(
    bytes_j: &[u8],
    bytes_k: &[u8],
    n_indiv_sel: usize,
    prod: &[f64; 16],
) -> f64 {
    let mut acc = 0.0f64;
    let full_bytes = n_indiv_sel / 4;
    let rem = n_indiv_sel % 4;
    for b in 0..full_bytes {
        let bj = bytes_j[b] as usize;
        let bk = bytes_k[b] as usize;
        acc += prod[((bj & 0b11) << 2) | (bk & 0b11)];
        acc += prod[(((bj >> 2) & 0b11) << 2) | ((bk >> 2) & 0b11)];
        acc += prod[(((bj >> 4) & 0b11) << 2) | ((bk >> 4) & 0b11)];
        acc += prod[(((bj >> 6) & 0b11) << 2) | ((bk >> 6) & 0b11)];
    }
    if rem > 0 {
        let bj = bytes_j[full_bytes] as usize;
        let bk = bytes_k[full_bytes] as usize;
        for r in 0..rem {
            let cj = (bj >> (r * 2)) & 0b11;
            let ck = (bk >> (r * 2)) & 0b11;
            acc += prod[(cj << 2) | ck];
        }
    }
    acc
}

/// Scatter-gather path used when only a subset of individuals is selected
/// (`--keep`). `iid_positions` holds (byte_idx, shift) per kept individual.
#[inline]
fn pair_centered_dot_subset(
    bytes_j: &[u8],
    bytes_k: &[u8],
    iid_positions: &[IidPos],
    prod: &[f64; 16],
) -> f64 {
    let mut acc = 0.0f64;
    for pos in iid_positions {
        let bj = bytes_j[pos.byte_idx];
        let bk = bytes_k[pos.byte_idx];
        let cj = ((bj >> pos.shift) & 0b11) as usize;
        let ck = ((bk >> pos.shift) & 0b11) as usize;
        acc += prod[(cj << 2) | ck];
    }
    acc
}

/// AVX2 joint-code histogram. Returns `counts[code] = #individuals with
/// joint genotype code `(cj << 2) | ck`. Processes 32 BED bytes (= 128
/// individuals) per AVX2 iteration; scalar tail handles the residual.
///
/// # Safety
/// Requires AVX2. We rely on `target-feature=+avx2` in `.cargo/config.toml`
/// for the production targets and on any x86_64 CPU shipped since 2013.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn pair_joint_histogram_avx2(
    bytes_j: &[u8],
    bytes_k: &[u8],
    n_indiv_sel: usize,
) -> [u32; 16] {
    use std::arch::x86_64::*;
    unsafe {
        let mut counts = [0u32; 16];

        let n_full_bytes = n_indiv_sel / 4;
        let rem_indiv = n_indiv_sel % 4;
        let chunks = n_full_bytes / 32;
        let tail_byte_start = chunks * 32;

        let mask_03 = _mm256_set1_epi8(0x03);

        // Hot loop: 32 BED bytes = 128 individuals per AVX2 vector pair.
        for chunk_i in 0..chunks {
            let p_j = bytes_j.as_ptr().add(chunk_i * 32) as *const __m256i;
            let p_k = bytes_k.as_ptr().add(chunk_i * 32) as *const __m256i;
            let ymm_j = _mm256_loadu_si256(p_j);
            let ymm_k = _mm256_loadu_si256(p_k);

            // Extract 4 positional code vectors per byte via shift+mask.
            // _mm256_srli_epi16 shifts 16-bit lanes; the AND-with-0x03 kills
            // any cross-byte contamination so the result holds the 2-bit
            // code in each byte's low bits.
            let cj0 = _mm256_and_si256(ymm_j, mask_03);
            let ck0 = _mm256_and_si256(ymm_k, mask_03);
            let cj1 = _mm256_and_si256(_mm256_srli_epi16::<2>(ymm_j), mask_03);
            let ck1 = _mm256_and_si256(_mm256_srli_epi16::<2>(ymm_k), mask_03);
            let cj2 = _mm256_and_si256(_mm256_srli_epi16::<4>(ymm_j), mask_03);
            let ck2 = _mm256_and_si256(_mm256_srli_epi16::<4>(ymm_k), mask_03);
            let cj3 = _mm256_and_si256(_mm256_srli_epi16::<6>(ymm_j), mask_03);
            let ck3 = _mm256_and_si256(_mm256_srli_epi16::<6>(ymm_k), mask_03);

            // Combine cj/ck into a 4-bit joint code per byte: `(cj << 2) | ck`.
            // Since each code holds values 0..=3 (top 6 bits already zero), the
            // 16-bit-lane left-shift never spills into a neighboring byte.
            let j0 = _mm256_or_si256(_mm256_slli_epi16::<2>(cj0), ck0);
            let j1 = _mm256_or_si256(_mm256_slli_epi16::<2>(cj1), ck1);
            let j2 = _mm256_or_si256(_mm256_slli_epi16::<2>(cj2), ck2);
            let j3 = _mm256_or_si256(_mm256_slli_epi16::<2>(cj3), ck3);

            histogram_16(j0, &mut counts);
            histogram_16(j1, &mut counts);
            histogram_16(j2, &mut counts);
            histogram_16(j3, &mut counts);
        }

        // Scalar tail: process the residual full bytes byte-by-byte.
        for b in tail_byte_start..n_full_bytes {
            let bj = *bytes_j.get_unchecked(b) as usize;
            let bk = *bytes_k.get_unchecked(b) as usize;
            *counts.get_unchecked_mut(((bj & 0b11) << 2) | (bk & 0b11)) += 1;
            *counts.get_unchecked_mut((((bj >> 2) & 0b11) << 2) | ((bk >> 2) & 0b11)) += 1;
            *counts.get_unchecked_mut((((bj >> 4) & 0b11) << 2) | ((bk >> 4) & 0b11)) += 1;
            *counts.get_unchecked_mut((((bj >> 6) & 0b11) << 2) | ((bk >> 6) & 0b11)) += 1;
        }

        // Last partial byte: only the `rem_indiv` low positions are valid.
        if rem_indiv > 0 {
            let bj = *bytes_j.get_unchecked(n_full_bytes) as usize;
            let bk = *bytes_k.get_unchecked(n_full_bytes) as usize;
            for pos in 0..rem_indiv {
                let cj = (bj >> (pos * 2)) & 0b11;
                let ck = (bk >> (pos * 2)) & 0b11;
                *counts.get_unchecked_mut((cj << 2) | ck) += 1;
            }
        }

        counts
    }
}

/// AArch64 NEON joint-code histogram. Same shape as the AVX2 kernel but with
/// 128-bit Q-registers (16 BED bytes = 64 individuals per chunk).
///
/// Critical perf detail: AArch64 has no movemask; the natural per-bin idiom
/// `vshrq_n_u8::<7>(cmpeq) → vaddvq_u8` puts the high-latency vaddvq (4-7
/// cycles on Apple Silicon Firestorm/Avalanche) on the critical path 16
/// times per joint vector, easily losing to LLVM's auto-vectorization of the
/// scalar inner loop. This kernel instead keeps 16 in-register byte
/// accumulators (`bin_counters[16]: uint8x16_t`) and updates them via
/// `vsubq_u8(acc, cmpeq)` — `cmpeq` returns `0xFF` (= −1 mod 256) on match,
/// so the subtract becomes a per-lane increment with throughput ~1 op/cycle
/// and no horizontal-reduce dependency. Periodic flushes (every
/// `BATCH_CHUNKS = 60` chunks, where `60 × 4 positions = 240 < 255`) widen
/// the accumulators to u32 before they overflow.
///
/// # Safety
/// NEON is a baseline aarch64 feature; every aarch64 target we care about
/// (Apple Silicon, AWS Graviton, ARMv8 servers) has it.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn pair_joint_histogram_neon(
    bytes_j: &[u8],
    bytes_k: &[u8],
    n_indiv_sel: usize,
) -> [u32; 16] {
    use std::arch::aarch64::*;
    unsafe {
        let mut counts = [0u32; 16];

        let n_full_bytes = n_indiv_sel / 4;
        let rem_indiv = n_indiv_sel % 4;
        let chunks = n_full_bytes / 16;
        let tail_byte_start = chunks * 16;

        let mask_03 = vdupq_n_u8(0x03);

        // Per-bin in-register accumulators: counts[bin][lane] = matches at
        // lane `lane` for bin `bin`. Each chunk adds at most 4 to any lane
        // (4 positions per byte). u8 caps at 255 → safe for 63 chunks.
        const BATCH_CHUNKS: usize = 60;
        let mut bin_acc = [vdupq_n_u8(0); 16];
        let code_vecs: [uint8x16_t; 16] = std::array::from_fn(|c| vdupq_n_u8(c as u8));

        let mut chunks_in_batch = 0usize;
        for chunk_i in 0..chunks {
            let p_j = bytes_j.as_ptr().add(chunk_i * 16);
            let p_k = bytes_k.as_ptr().add(chunk_i * 16);
            let v_j = vld1q_u8(p_j);
            let v_k = vld1q_u8(p_k);

            // NEON has true byte-level shifts; the AND-with-0x03 only kills
            // residual bits at the high positions, not byte-spillover.
            let cj0 = vandq_u8(v_j, mask_03);
            let ck0 = vandq_u8(v_k, mask_03);
            let cj1 = vandq_u8(vshrq_n_u8::<2>(v_j), mask_03);
            let ck1 = vandq_u8(vshrq_n_u8::<2>(v_k), mask_03);
            let cj2 = vandq_u8(vshrq_n_u8::<4>(v_j), mask_03);
            let ck2 = vandq_u8(vshrq_n_u8::<4>(v_k), mask_03);
            let cj3 = vshrq_n_u8::<6>(v_j); // top 2 bits already zero
            let ck3 = vshrq_n_u8::<6>(v_k);

            let j0 = vorrq_u8(vshlq_n_u8::<2>(cj0), ck0);
            let j1 = vorrq_u8(vshlq_n_u8::<2>(cj1), ck1);
            let j2 = vorrq_u8(vshlq_n_u8::<2>(cj2), ck2);
            let j3 = vorrq_u8(vshlq_n_u8::<2>(cj3), ck3);

            // For each bin, subtract the cmpeq mask from the accumulator.
            // cmpeq → 0xFF for match (= -1 mod 256), 0x00 for non-match;
            // `acc - cmpeq` increments lanes that matched. Manually unrolled
            // so the compiler keeps `bin_acc[*]` in registers across the
            // batch and only spills at the flush.
            macro_rules! upd {
                ($bin:literal) => {{
                    let cv = code_vecs[$bin];
                    let m0 = vceqq_u8(j0, cv);
                    let m1 = vceqq_u8(j1, cv);
                    let m2 = vceqq_u8(j2, cv);
                    let m3 = vceqq_u8(j3, cv);
                    let a0 = vsubq_u8(bin_acc[$bin], m0);
                    let a1 = vsubq_u8(a0, m1);
                    let a2 = vsubq_u8(a1, m2);
                    bin_acc[$bin] = vsubq_u8(a2, m3);
                }};
            }
            upd!(0);
            upd!(1);
            upd!(2);
            upd!(3);
            upd!(4);
            upd!(5);
            upd!(6);
            upd!(7);
            upd!(8);
            upd!(9);
            upd!(10);
            upd!(11);
            upd!(12);
            upd!(13);
            upd!(14);
            upd!(15);

            chunks_in_batch += 1;
            if chunks_in_batch == BATCH_CHUNKS {
                flush_bin_acc_neon(&mut bin_acc, &mut counts);
                chunks_in_batch = 0;
            }
        }
        if chunks_in_batch > 0 {
            flush_bin_acc_neon(&mut bin_acc, &mut counts);
        }

        // Scalar tail: residual full bytes byte-by-byte.
        for b in tail_byte_start..n_full_bytes {
            let bj = *bytes_j.get_unchecked(b) as usize;
            let bk = *bytes_k.get_unchecked(b) as usize;
            *counts.get_unchecked_mut(((bj & 0b11) << 2) | (bk & 0b11)) += 1;
            *counts.get_unchecked_mut((((bj >> 2) & 0b11) << 2) | ((bk >> 2) & 0b11)) += 1;
            *counts.get_unchecked_mut((((bj >> 4) & 0b11) << 2) | ((bk >> 4) & 0b11)) += 1;
            *counts.get_unchecked_mut((((bj >> 6) & 0b11) << 2) | ((bk >> 6) & 0b11)) += 1;
        }

        if rem_indiv > 0 {
            let bj = *bytes_j.get_unchecked(n_full_bytes) as usize;
            let bk = *bytes_k.get_unchecked(n_full_bytes) as usize;
            for pos in 0..rem_indiv {
                let cj = (bj >> (pos * 2)) & 0b11;
                let ck = (bk >> (pos * 2)) & 0b11;
                *counts.get_unchecked_mut((cj << 2) | ck) += 1;
            }
        }

        counts
    }
}

/// Reduce per-bin u8 lane accumulators into the u32 `counts` array, then
/// reset to zero. `vpaddlq_u8` widens 16 u8 lanes pairwise into 8 u16 lanes;
/// `vaddvq_u16` then horizontally reduces to a u16 sum (max 16 × 240 = 3840,
/// fits cleanly).
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
#[inline]
unsafe fn flush_bin_acc_neon(
    bin_acc: &mut [std::arch::aarch64::uint8x16_t; 16],
    counts: &mut [u32; 16],
) {
    use std::arch::aarch64::*;
    for bin in 0..16 {
        let widened = vpaddlq_u8(bin_acc[bin]);
        let sum = vaddvq_u16(widened) as u32;
        counts[bin] += sum;
        bin_acc[bin] = vdupq_n_u8(0);
    }
}

/// For each of 16 possible joint codes, count occurrences in `joint`
/// (32-byte AVX2 vector) and add to `counts`. Issues 16 independent
/// cmpeq → movemask → popcnt chains; modern OoO cores overlap them on
/// the vector and integer pipelines for ~16-20 cycles per call.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn histogram_16(joint: std::arch::x86_64::__m256i, counts: &mut [u32; 16]) {
    use std::arch::x86_64::*;
    unsafe {
        macro_rules! count_one {
            ($code:literal) => {{
                let code_vec = _mm256_set1_epi8($code as i8);
                let mask = _mm256_cmpeq_epi8(joint, code_vec);
                let bits = _mm256_movemask_epi8(mask) as u32;
                *counts.get_unchecked_mut($code) += bits.count_ones();
            }};
        }
        count_one!(0);
        count_one!(1);
        count_one!(2);
        count_one!(3);
        count_one!(4);
        count_one!(5);
        count_one!(6);
        count_one!(7);
        count_one!(8);
        count_one!(9);
        count_one!(10);
        count_one!(11);
        count_one!(12);
        count_one!(13);
        count_one!(14);
        count_one!(15);
    }
}

/// Compute LD scores for all SNPs via direct iteration over per-SNP windows
/// from packed BED bytes. Returns `(l2, maf_per_snp)`.
///
/// Semantics match `--snp-level-masking` (per-SNP exact windows, the LDSC
/// paper's `ℓ_j = Σ_k r²_{jk}`), not Python LDSC's chunk-level approximation.
///
/// `inner_threads` controls the intra-chromosome parallelism budget. The
/// dispatcher picks this based on the chr count vs core count: many chrs →
/// outer-only (pass `1`); few chrs → all cores inside one chr (pass
/// `num_cores`). Passing `0` is treated as `1`.
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
    inner_threads: usize,
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

    // -------- Intra-chromosome parallelism --------------------------------
    //
    // Split the j-range across rayon threads, each with its own `l2_local`
    // accumulator. Without this, the chr-major outer parallelism in `mod.rs`
    // leaves a single-chr run effectively single-threaded, while the GEMM
    // exact path uses faer's internal rayon for the same case. Per-thread
    // buffers avoid races on the scatter into `l2[k]` from neighboring
    // threads' j-ranges.
    //
    // Memory: `num_threads × m × n_annot × 8` bytes peak per chr. For
    // num_threads = 8, m = 80k, n_annot = 1: ~5MB per chr. Acceptable.
    //
    // Load balance: each j has work ∝ `j - block_left[j]`. Window widths are
    // roughly uniform across a chromosome (constant Mb radius), so an even
    // j-range split gives balanced work. Rayon's work-stealing can also
    // rebalance if some threads finish early.
    // Effective threads: caller-supplied budget, capped at `m` to avoid
    // empty-range tasks and floored at 1.
    let effective_threads = inner_threads.max(1).min(m).max(1);
    let chunk = m.div_ceil(effective_threads);

    // Compute per-thread partial l2 accumulators in parallel.
    let block_left_ref = &block_left;
    let stats_ref = &stats;
    let annot_eff_ref = &annot_eff;
    let iid_pos_ref = &iid_positions;
    let mmap_bed_ref = mmap_bed.as_ref();
    let all_snps_ref = all_snps;
    let bed_path_ref = bed_path;

    let partials: Result<Vec<MatF>> = (0..effective_threads)
        .into_par_iter()
        .map(|t| -> Result<MatF> {
            let j_start = t * chunk;
            let j_end = ((t + 1) * chunk).min(m);
            let mut local = MatF::zeros(m, n_annot);
            // Per-thread BED reader for the non-mmap fallback. `Bed` holds a
            // BufReader with seek state, so we can't share one across threads.
            let mut bed_local: Option<Bed> = if mmap_bed_ref.is_none() {
                Some(
                    Bed::builder(bed_path_ref)
                        .build()
                        .with_context(|| format!("opening BED '{}' (thread {})", bed_path_ref, t))?,
                )
            } else {
                None
            };
            let mut bytes_j_buf: Vec<u8> = vec![0u8; bytes_per_snp];
            let mut bytes_k_buf: Vec<u8> = vec![0u8; bytes_per_snp];

            for j in j_start..j_end {
                let sj = stats_ref[j];
                // Diagonal: r²(j, j) = 1 exactly → contributes annot_eff[j, c].
                for c in 0..n_annot {
                    local[(j, c)] += annot_eff_ref[j * stride + c];
                }

                let bl = block_left_ref[j];
                if bl >= j {
                    continue;
                }

                // Zero-variance SNP j: every pair gets r²_unbiased = r2u_b
                // (matches the GEMM path's zero-column behavior). Skip byte
                // reads entirely.
                if sj.inv_std == 0.0 {
                    for k in bl..j {
                        let r2u = r2u_b;
                        for c in 0..n_annot {
                            local[(j, c)] += r2u * annot_eff_ref[k * stride + c];
                            local[(k, c)] += r2u * annot_eff_ref[j * stride + c];
                        }
                    }
                    continue;
                }

                let bed_idx_j = all_snps_ref[j].bed_idx;
                if let Some(ref mut bed) = bed_local {
                    bed.read_snp_bytes(bed_idx_j, &mut bytes_j_buf)
                        .with_context(|| format!("reading BED SNP {}", bed_idx_j))?;
                }

                for k in bl..j {
                    let sk = stats_ref[k];

                    let r2u = if sk.inv_std == 0.0 {
                        r2u_b
                    } else {
                        let bed_idx_k = all_snps_ref[k].bed_idx;
                        let (bytes_j_view, bytes_k_view): (&[u8], &[u8]) =
                            if let Some(mb) = mmap_bed_ref {
                                (mb.snp_bytes(bed_idx_j, 1), mb.snp_bytes(bed_idx_k, 1))
                            } else {
                                let bed = bed_local.as_mut().unwrap();
                                bed.read_snp_bytes(bed_idx_k, &mut bytes_k_buf)
                                    .with_context(|| format!("reading BED SNP {}", bed_idx_k))?;
                                (&bytes_j_buf[..], &bytes_k_buf[..])
                            };

                        let centered_dot = pair_centered_dot(
                            bytes_j_view,
                            bytes_k_view,
                            n_indiv_sel,
                            iid_pos_ref,
                            all_iids,
                            sj.mean,
                            sk.mean,
                        );
                        let r = centered_dot * sj.inv_std * sk.inv_std * n_inv;
                        let r2_raw = (r * r).min(1.0);
                        r2_raw * r2u_a + r2u_b
                    };

                    for c in 0..n_annot {
                        local[(j, c)] += r2u * annot_eff_ref[k * stride + c];
                        local[(k, c)] += r2u * annot_eff_ref[j * stride + c];
                    }
                }
            }
            Ok(local)
        })
        .collect();

    let partials = partials?;

    // Reduce: sum the per-thread MatF buffers into the final l2.
    // We do this sequentially since each partial is m × n_annot f64s and the
    // reduction is memory-bound — rayon parallelism here adds little.
    let mut l2 = MatF::zeros(m, n_annot);
    for partial in &partials {
        for i in 0..m {
            for c in 0..n_annot {
                l2[(i, c)] += partial[(i, c)];
            }
        }
    }

    Ok((l2, maf_per_snp))
}
