use super::normalize::{normalize_col_f32_with_stats, normalize_col_f64_with_stats, sum_sumsq_f32};
use super::window::{WindowMode, get_block_lefts_by_chr};
use crate::bed::{Bed, ChunkReader, IidPos, MmapBed, ReadOptions};
#[cfg(feature = "gpu")]
use crate::gpu::{GpuConfig, GpuContext};
use crate::la::{
    MatF, MatF32, mat_add_in_place, mat_copy_from, mat_slice, mat_slice_f32, mat_slice_mut,
    mat_slice_mut_f32, matmul_tn_to, matmul_tn_to_f32, matmul_to, par_default,
};
use crate::parse::BimRecord;
use anyhow::{Context, Result};
use faer::{Accum, Par};
use rayon::prelude::*;
use std::collections::VecDeque;

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

// ── SendPtr wrapper for raw pointers across rayon threads ──────────────────
/// Wrapper that implements Send+Sync for raw pointers. The caller is
/// responsible for ensuring disjoint access (each thread touches a
/// different column of the matrix).
#[derive(Clone, Copy)]
struct SendPtr<T>(*mut T);
unsafe impl<T> Send for SendPtr<T> {}
unsafe impl<T> Sync for SendPtr<T> {}
impl<T> SendPtr<T> {
    /// Return the raw pointer. Using a method instead of `.0` field access
    /// prevents Rust 2021 disjoint-field capture from stripping the
    /// Send+Sync wrapper in closures.
    #[inline(always)]
    fn ptr(self) -> *mut T {
        self.0
    }
}

// ── CountSketch projection ──────────────────────────────────────────────────
// Hash each individual to a random bucket with a random sign.
// O(N) hash + sign arrays; O(N×c) scatter-add, flat cost in d.

/// CountSketch projection: hash each individual to a bucket with a random sign.
struct CountSketchProj {
    bucket: Vec<u32>, // h(i) ∈ {0,...,d-1} for each individual
    sign_f64: Vec<f64>,
    sign_f32: Vec<f32>,
    d: usize,
}

impl CountSketchProj {
    fn new(n_indiv: usize, d: usize, seed: u64) -> Self {
        let mut rng = fastrand::Rng::with_seed(seed);
        let bucket: Vec<u32> = (0..n_indiv).map(|_| rng.u32(..d as u32)).collect();
        let sign_f64: Vec<f64> = (0..n_indiv)
            .map(|_| if rng.bool() { 1.0 } else { -1.0 })
            .collect();
        let sign_f32: Vec<f32> = sign_f64.iter().map(|&s| s as f32).collect();
        Self {
            bucket,
            sign_f64,
            sign_f32,
            d,
        }
    }

    /// Project b_full_f64[n_indiv x c] into out[d x c] via scatter-add, then
    /// apply ratio estimator renormalization so ||out[:,j]||^2 = n.
    ///
    /// Parallelized over columns (j=0..c) via rayon. Each column's scatter-add
    /// is independent, giving c-way parallelism (c=chunk_size, typically 200).
    fn project_f64(&self, b_full: &MatF, n_indiv: usize, c: usize, n: f64, out: &mut MatF) {
        let d = self.d;
        let bucket = &self.bucket;
        let sign = &self.sign_f64;

        // Get raw pointer + stride for parallel column access.
        // Safety: faer Mat is column-major; columns at ptr + j*col_stride are
        // disjoint memory regions of length >= d. Each rayon task writes only
        // to its own column j, so no data races.
        let out_ptr = out.as_mut().as_ptr_mut();
        let out_col_stride = out.col_stride();
        let src_ptr = b_full.as_ref().as_ptr();
        let src_col_stride = b_full.col_stride();

        // Wrap raw pointers for Send across rayon threads.
        let out_send = SendPtr(out_ptr);
        let src_send = SendPtr(src_ptr as *mut f64);

        (0..c).into_par_iter().for_each(|j| {
            // Safety: each thread accesses a disjoint column of out and reads
            // a disjoint column of b_full.
            unsafe {
                let dst = std::slice::from_raw_parts_mut(
                    out_send.ptr().offset(j as isize * out_col_stride),
                    d,
                );
                let src = std::slice::from_raw_parts(
                    (src_send.ptr() as *const f64).offset(j as isize * src_col_stride),
                    n_indiv,
                );

                // Zero output column.
                for v in dst.iter_mut() {
                    *v = 0.0;
                }

                // Scatter-add: out[bucket[i], j] += sign[i] * src[i]
                for i in 0..n_indiv {
                    *dst.get_unchecked_mut(bucket[i] as usize) += sign[i] * src[i];
                }

                // Ratio estimator renorm: rescale so ||col||^2 = n
                let mut nrm_sq = 0.0f64;
                for &v in dst.iter() {
                    nrm_sq += v * v;
                }
                if nrm_sq > 0.0 {
                    let scale = (n / nrm_sq).sqrt();
                    for v in dst.iter_mut() {
                        *v *= scale;
                    }
                }
            }
        });
    }

    /// f32 variant. Parallelized over columns via rayon.
    fn project_f32(&self, b_full: &MatF32, n_indiv: usize, c: usize, n: f64, out: &mut MatF32) {
        let d = self.d;
        let bucket = &self.bucket;
        let sign = &self.sign_f32;

        let out_ptr = out.as_mut().as_ptr_mut();
        let out_col_stride = out.col_stride();
        let src_ptr = b_full.as_ref().as_ptr();
        let src_col_stride = b_full.col_stride();

        let out_send = SendPtr(out_ptr);
        let src_send = SendPtr(src_ptr as *mut f32);

        (0..c).into_par_iter().for_each(|j| {
            unsafe {
                let dst = std::slice::from_raw_parts_mut(
                    out_send.ptr().offset(j as isize * out_col_stride),
                    d,
                );
                let src = std::slice::from_raw_parts(
                    (src_send.ptr() as *const f32).offset(j as isize * src_col_stride),
                    n_indiv,
                );

                for v in dst.iter_mut() {
                    *v = 0.0;
                }

                for i in 0..n_indiv {
                    *dst.get_unchecked_mut(bucket[i] as usize) += sign[i] * src[i];
                }

                // Ratio estimator renorm (accumulate in f64 for precision).
                let mut nrm_sq = 0.0f64;
                for &v in dst.iter() {
                    nrm_sq += (v as f64) * (v as f64);
                }
                if nrm_sq > 0.0 {
                    let scale = (n / nrm_sq).sqrt() as f32;
                    for v in dst.iter_mut() {
                        *v *= scale;
                    }
                }
            }
        });
    }
}

// ── Fused BED-decode-normalize-project ─────────────────────────────────────
// Avoids materializing the full N×c intermediate matrix (b_full).
// Pass 1: compute per-SNP stats (mean, inv_std, maf) from packed BED bytes.
// Pass 2: tiled decode+normalize+project, accumulating into out[d×c].

// BED decode LUT: bits → genotype value (or NaN for missing).
// The per-(mean, inv_std) normalized-byte LUT builders used by the
// fused CountSketch scatter live in `crate::wasm_simd::scatter` (so
// the wasm SAB-pool dispatch can share them with the native rayon
// path without a circular module dep).

// Byte-level stats LUT: each of the 256 possible BED byte values encodes 4
// genotypes; this LUT gives (sum, count, sum_sq) for the full byte.
#[derive(Clone, Copy)]
struct BedByteStats {
    sum: u8,
    count: u8,
    sum_sq: u8,
}

// Higher-moments LUT used by the MAF-aware bias correction path. Kept
// separate from `BedByteStats` so the existing fast paths don't pay the
// extra-load cost. Sum of G³ per byte ≤ 4·8 = 32 (fits u8). Sum of G⁴ per
// byte ≤ 4·16 = 64 (fits u8).
#[derive(Clone, Copy)]
struct BedByteMoments {
    sum_cube: u8,
    sum_4th: u8,
}

fn build_bed_byte_lut() -> [BedByteStats; 256] {
    std::array::from_fn(|b| {
        let byte = b as u8;
        let (mut sum, mut count, mut sum_sq) = (0u8, 0u8, 0u8);
        for k in 0..4 {
            match (byte >> (2 * k)) & 0b11 {
                0 => {
                    sum += 2;
                    count += 1;
                    sum_sq += 4;
                } // hom A1 = 2
                1 => {} // missing
                2 => {
                    sum += 1;
                    count += 1;
                    sum_sq += 1;
                } // het = 1
                3 => {
                    count += 1;
                } // hom A2 = 0
                _ => {}
            }
        }
        BedByteStats { sum, count, sum_sq }
    })
}

/// Compute (sum, count, sum_sq) of non-missing genotypes for one SNP from raw bytes.
/// Genotypes are 0, 1, 2 (count_a1=true convention).
fn snp_stats_from_bytes(
    snp_bytes: &[u8],
    n_indiv: usize,
    iid_positions: &[IidPos],
    all_iids: bool,
    byte_lut: &[BedByteStats; 256],
) -> (u32, u32, u32) {
    if all_iids {
        let full_bytes = n_indiv / 4;
        let rem = n_indiv % 4;
        let (mut sum, mut count, mut sum_sq) = (0u32, 0u32, 0u32);
        for &byte in &snp_bytes[..full_bytes] {
            let s = byte_lut[byte as usize];
            sum += s.sum as u32;
            count += s.count as u32;
            sum_sq += s.sum_sq as u32;
        }
        if rem > 0 {
            let byte = snp_bytes[full_bytes];
            for k in 0..rem {
                match (byte >> (2 * k)) & 0b11 {
                    0 => {
                        sum += 2;
                        count += 1;
                        sum_sq += 4;
                    }
                    1 => {}
                    2 => {
                        sum += 1;
                        count += 1;
                        sum_sq += 1;
                    }
                    3 => {
                        count += 1;
                    }
                    _ => {}
                }
            }
        }
        (sum, count, sum_sq)
    } else {
        let (mut sum, mut count, mut sum_sq) = (0u32, 0u32, 0u32);
        for pos in iid_positions {
            let byte = snp_bytes[pos.byte_idx];
            match (byte >> pos.shift) & 0b11 {
                0 => {
                    sum += 2;
                    count += 1;
                    sum_sq += 4;
                }
                1 => {}
                2 => {
                    sum += 1;
                    count += 1;
                    sum_sq += 1;
                }
                3 => {
                    count += 1;
                }
                _ => {}
            }
        }
        (sum, count, sum_sq)
    }
}

/// Build the higher-moments LUT (sum of G³ and G⁴ per byte). G ∈ {0, 1, 2}:
/// G³ values are {0, 1, 8}, G⁴ values are {0, 1, 16}.
fn build_bed_byte_moments_lut() -> [BedByteMoments; 256] {
    std::array::from_fn(|b| {
        let byte = b as u8;
        let (mut sum_cube, mut sum_4th) = (0u8, 0u8);
        for k in 0..4 {
            match (byte >> (2 * k)) & 0b11 {
                0 => {
                    sum_cube += 8;
                    sum_4th += 16;
                } // hom A1 = 2 → 2³=8, 2⁴=16
                1 => {} // missing
                2 => {
                    sum_cube += 1;
                    sum_4th += 1;
                } // het = 1 → 1³=1, 1⁴=1
                3 => {} // hom A2 = 0 → 0³=0, 0⁴=0
                _ => {}
            }
        }
        BedByteMoments { sum_cube, sum_4th }
    })
}

/// Compute (sum_cube, sum_4th) of non-missing genotypes for one SNP. Parallels
/// `snp_stats_from_bytes`. Only called when `--sketch-maf-aware` is active so
/// the K=1 sketch path is byte-identical to the pre-MAF-aware code.
fn snp_higher_moments_from_bytes(
    snp_bytes: &[u8],
    n_indiv: usize,
    iid_positions: &[IidPos],
    all_iids: bool,
    moments_lut: &[BedByteMoments; 256],
) -> (u32, u32) {
    if all_iids {
        let full_bytes = n_indiv / 4;
        let rem = n_indiv % 4;
        let (mut sum_cube, mut sum_4th) = (0u32, 0u32);
        for &byte in &snp_bytes[..full_bytes] {
            let m = moments_lut[byte as usize];
            sum_cube += m.sum_cube as u32;
            sum_4th += m.sum_4th as u32;
        }
        if rem > 0 {
            let byte = snp_bytes[full_bytes];
            for k in 0..rem {
                match (byte >> (2 * k)) & 0b11 {
                    0 => {
                        sum_cube += 8;
                        sum_4th += 16;
                    }
                    1 => {}
                    2 => {
                        sum_cube += 1;
                        sum_4th += 1;
                    }
                    3 => {}
                    _ => {}
                }
            }
        }
        (sum_cube, sum_4th)
    } else {
        let (mut sum_cube, mut sum_4th) = (0u32, 0u32);
        for pos in iid_positions {
            let byte = snp_bytes[pos.byte_idx];
            match (byte >> pos.shift) & 0b11 {
                0 => {
                    sum_cube += 8;
                    sum_4th += 16;
                }
                1 => {}
                2 => {
                    sum_cube += 1;
                    sum_4th += 1;
                }
                3 => {}
                _ => {}
            }
        }
        (sum_cube, sum_4th)
    }
}

/// Per-SNP kurtosis quantity `κ_j = K_jj / N² = (Σ_i X_ij⁴) / N²` where X_j is
/// the standardized vector (mean-imputed at missing positions). Derived from
/// raw moments via the centered-fourth-moment identity:
///
///   Σ (G − μ)⁴ = sum_4th − 4μ·sum_cube + 6μ²·sum_sq − 3μ⁴·count
///
/// (using `μ·count = sum` to fold the `−4μ³·sum + μ⁴·count` terms together).
/// Standardization divides by σ⁴ where σ² = centered_ss/N (the normalize.rs
/// convention), giving `K_jj = Σ(G-μ)⁴ × N² / centered_ss²` so
/// `K_jj/N² = Σ(G-μ)⁴ / centered_ss²`. See `docs/countsketch-math-analysis.md` §14.
#[inline]
fn per_snp_kappa(sum: u32, count: u32, sum_sq: u32, sum_cube: u32, sum_4th: u32) -> f64 {
    if count == 0 {
        return 0.0;
    }
    let count_f = count as f64;
    let mean = sum as f64 / count_f;
    let m4 = sum_4th as f64 - 4.0 * mean * sum_cube as f64 + 6.0 * mean * mean * sum_sq as f64
        - 3.0 * mean.powi(4) * count_f;
    let centered_ss = sum_sq as f64 - count_f * mean * mean;
    if centered_ss <= 0.0 {
        return 0.0;
    }
    m4 / (centered_ss * centered_ss)
}

/// Per-SNP κ_j = K_jj/N² computed directly from a decoded f32 column (raw
/// genotype values 0/1/2 with NaN for missing). Used in the non-fused sketch
/// path where raw bytes have already been decoded into a matrix. The result
/// matches `per_snp_kappa` byte-for-byte modulo f64 ordering. `mean` and
/// `centered_ss` are the moments over non-NaN entries (i.e. `sum/count` and
/// `sum_sq - count*mean²` from the existing stats pass).
#[inline]
fn per_snp_kappa_from_col_f32(col: &[f32], mean: f64, centered_ss: f64) -> f64 {
    if centered_ss <= 0.0 {
        return 0.0;
    }
    let mut m4 = 0.0f64;
    for &v in col {
        if !v.is_nan() {
            let d = v as f64 - mean;
            let d2 = d * d;
            m4 += d2 * d2;
        }
    }
    m4 / (centered_ss * centered_ss)
}

/// f64 counterpart of `per_snp_kappa_from_col_f32`.
#[inline]
fn per_snp_kappa_from_col_f64(col: &[f64], mean: f64, centered_ss: f64) -> f64 {
    if centered_ss <= 0.0 {
        return 0.0;
    }
    let mut m4 = 0.0f64;
    for &v in col {
        if !v.is_nan() {
            let d = v - mean;
            let d2 = d * d;
            m4 += d2 * d2;
        }
    }
    m4 / (centered_ss * centered_ss)
}

// ── Fused CountSketch: BED-decode-normalize-scatter-add ─────────────────────
// Avoids materializing the full N×c intermediate matrix AND avoids any GEMM.
// Pass 1: compute per-SNP stats (mean, inv_std, maf) from packed BED bytes.
// Pass 2: parallel over columns — decode each genotype, normalize, scatter-add
//         to the appropriate bucket. O(N×c) total, column-parallel via rayon.
//
// Fused CountSketch has NO GEMM at all — just scatter-adds. The column-parallel
// pattern gives c-way (c=200) parallelism with no SIMD efficiency loss.

/// Fused BED-decode-normalize-CountSketch (f32 variant, primary path since sketch auto-enables f32).
/// Reads raw packed BED bytes, computes per-SNP statistics, then scatter-adds
/// normalized genotypes into d-bucket sketch vectors — all without materializing
/// the full N×c genotype matrix.
///
/// Parallelized over columns (j=0..c) via rayon. Each column's scatter-add is
/// independent, giving c-way parallelism (c=chunk_size, typically 200).
///
/// On return, `out[0..d, 0..c]` contains the projected+renormalized columns,
/// `maf_out[0..c]` contains the per-SNP MAF, and (when `kurt_out` is `Some`)
/// `kurt_out[0..c]` contains the per-SNP kurtosis quantity κ_j = K_jj/N²
/// used by the MAF-aware quadratic bias correction.
#[allow(clippy::too_many_arguments)]
fn countsketch_fused_project_f32(
    raw_bytes: &[u8],
    bytes_per_snp: usize,
    n_indiv: usize,
    c: usize,
    cs: &CountSketchProj,
    n: f64,
    iid_positions: &[IidPos],
    all_iids: bool,
    maf_out: &mut [f64],
    kurt_out: Option<&mut [f64]>,
    out: &mut MatF32,
) {
    let d = cs.d;
    let byte_lut = build_bed_byte_lut();
    let moments_lut = if kurt_out.is_some() {
        Some(build_bed_byte_moments_lut())
    } else {
        None
    };

    // --- Pass 1: per-SNP statistics from packed bytes ---
    struct SnpStats {
        mean: f32,
        inv_std: f32,
        maf: f64,
        kurt: f64,
    }
    let stats: Vec<SnpStats> = (0..c)
        .map(|j| {
            let snp_bytes = &raw_bytes[j * bytes_per_snp..(j + 1) * bytes_per_snp];
            let (sum, count, sum_sq) =
                snp_stats_from_bytes(snp_bytes, n_indiv, iid_positions, all_iids, &byte_lut);
            let mean_f64 = if count > 0 {
                sum as f64 / count as f64
            } else {
                0.0
            };
            let freq = (mean_f64 / 2.0).clamp(0.0, 1.0);
            let maf = freq.min(1.0 - freq);
            let centered_ss = sum_sq as f64 - count as f64 * mean_f64 * mean_f64;
            let var = centered_ss / n_indiv as f64;
            let std = var.sqrt();
            let inv_std = if std > 0.0 {
                (1.0 / std) as f32
            } else {
                0.0f32
            };
            let kurt = if let Some(ref m_lut) = moments_lut {
                let (sum_cube, sum_4th) = snp_higher_moments_from_bytes(
                    snp_bytes,
                    n_indiv,
                    iid_positions,
                    all_iids,
                    m_lut,
                );
                per_snp_kappa(sum, count, sum_sq, sum_cube, sum_4th)
            } else {
                0.0
            };
            SnpStats {
                mean: mean_f64 as f32,
                inv_std,
                maf,
                kurt,
            }
        })
        .collect();

    // Flat parallel-array views of stats[].mean / .inv_std so the wasm
    // SAB-pool dispatch (which carries only raw pointers in `usize`)
    // and the native rayon path can both index per-SNP scalars without
    // touching the `SnpStats` Vec layout.
    let mean_buf: Vec<f32> = stats.iter().map(|s| s.mean).collect();
    let inv_std_buf: Vec<f32> = stats.iter().map(|s| s.inv_std).collect();

    for (j, s) in stats.iter().enumerate() {
        maf_out[j] = s.maf;
    }
    if let Some(kurt_slice) = kurt_out {
        for (j, s) in stats.iter().enumerate() {
            kurt_slice[j] = s.kurt;
        }
    }

    // --- Pass 2: column-parallel fused decode + normalize + scatter-add ---
    // Each column j is independent: read BED bytes for SNP j, decode each
    // individual's genotype, normalize, and scatter-add to out[bucket[i], j].
    let out_ptr = out.as_mut().as_ptr_mut();
    let out_col_stride = out.col_stride();
    let bucket = &cs.bucket;
    let sign = &cs.sign_f32;

    // On wasm32 with `+atomics`, rayon's `into_par_iter` is a no-op
    // (runs serially on the calling thread — no spindle / atomic-wait
    // integration). Dispatch through the manual SAB pool instead so
    // the per-outer scatter splits across all inner Workers. The pool
    // returns `false` (and we fall through) when it isn't initialised
    // or N <= 1; in that case the rayon path below runs the scatter
    // serially on the calling thread.
    #[cfg(all(target_arch = "wasm32", target_feature = "atomics"))]
    {
        let args = crate::wasm_simd::pool::ScatterArgsF32 {
            raw_bytes_ptr: raw_bytes.as_ptr() as usize,
            bytes_per_snp,
            n_indiv,
            c,
            bucket_ptr: bucket.as_ptr() as usize,
            sign_ptr: sign.as_ptr() as usize,
            iid_pos_ptr: iid_positions.as_ptr() as usize,
            iid_pos_len: iid_positions.len(),
            all_iids: if all_iids { 1 } else { 0 },
            n_norm: n,
            out_ptr: out_ptr as usize,
            out_col_stride,
            d,
            mean_ptr: mean_buf.as_ptr() as usize,
            inv_std_ptr: inv_std_buf.as_ptr() as usize,
        };
        // SAFETY: every slice we pulled `.as_ptr()` from above lives
        // until this function returns (Vec / slice borrows owned by
        // the calling stack frame). `parallel_scatter_f32` blocks
        // until all inner Workers finish, so the references are still
        // valid when the dispatch returns. `out` is exclusively owned
        // by us via the `&mut MatF32` argument.
        let dispatched = unsafe { crate::wasm_simd::pool::parallel_scatter_f32(args) };
        if dispatched {
            return;
        }
    }

    // Native rayon (or wasm without an initialised pool) path. Same
    // per-column body as the wasm dispatch — both routes funnel
    // through `scatter::scatter_one_column_f32` so the numerics are
    // identical regardless of where the work runs.
    let out_send = SendPtr(out_ptr);
    let raw_ptr_us = raw_bytes.as_ptr() as usize;
    let bucket_ptr_us = bucket.as_ptr() as usize;
    let sign_ptr_us = sign.as_ptr() as usize;
    let iid_pos_ptr_us = iid_positions.as_ptr() as usize;
    let iid_pos_len = iid_positions.len();
    let mean_ptr_us = mean_buf.as_ptr() as usize;
    let inv_std_ptr_us = inv_std_buf.as_ptr() as usize;

    (0..c).into_par_iter().for_each(|j| {
        // SAFETY: rayon hands each `j` to one task; output column j
        // is non-overlapping with other columns. All read-only views
        // (raw_bytes / bucket / sign / iid_positions / mean_buf /
        // inv_std_buf) outlive the join() implicit in for_each.
        unsafe {
            crate::wasm_simd::scatter::scatter_one_column_f32(
                j,
                raw_ptr_us as *const u8,
                bytes_per_snp,
                n_indiv,
                bucket_ptr_us as *const u32,
                sign_ptr_us as *const f32,
                iid_pos_ptr_us as *const IidPos,
                iid_pos_len,
                all_iids,
                *(mean_ptr_us as *const f32).add(j),
                *(inv_std_ptr_us as *const f32).add(j),
                n,
                out_send.ptr(),
                out_col_stride,
                d,
            );
        }
    });
}

/// Fused BED-decode-normalize-CountSketch (f64 variant).
/// Same algorithm as the f32 variant — see `countsketch_fused_project_f32`.
#[allow(clippy::too_many_arguments)]
fn countsketch_fused_project_f64(
    raw_bytes: &[u8],
    bytes_per_snp: usize,
    n_indiv: usize,
    c: usize,
    cs: &CountSketchProj,
    n: f64,
    iid_positions: &[IidPos],
    all_iids: bool,
    maf_out: &mut [f64],
    kurt_out: Option<&mut [f64]>,
    out: &mut MatF,
) {
    let d = cs.d;
    let byte_lut = build_bed_byte_lut();
    let moments_lut = if kurt_out.is_some() {
        Some(build_bed_byte_moments_lut())
    } else {
        None
    };

    // --- Pass 1: per-SNP statistics from packed bytes ---
    struct SnpStats {
        mean: f64,
        inv_std: f64,
        maf: f64,
        kurt: f64,
    }
    let stats: Vec<SnpStats> = (0..c)
        .map(|j| {
            let snp_bytes = &raw_bytes[j * bytes_per_snp..(j + 1) * bytes_per_snp];
            let (sum, count, sum_sq) =
                snp_stats_from_bytes(snp_bytes, n_indiv, iid_positions, all_iids, &byte_lut);
            let mean = if count > 0 {
                sum as f64 / count as f64
            } else {
                0.0
            };
            let freq = (mean / 2.0).clamp(0.0, 1.0);
            let maf = freq.min(1.0 - freq);
            let centered_ss = sum_sq as f64 - count as f64 * mean * mean;
            let var = centered_ss / n_indiv as f64;
            let std = var.sqrt();
            let inv_std = if std > 0.0 { 1.0 / std } else { 0.0 };
            let kurt = if let Some(ref m_lut) = moments_lut {
                let (sum_cube, sum_4th) = snp_higher_moments_from_bytes(
                    snp_bytes,
                    n_indiv,
                    iid_positions,
                    all_iids,
                    m_lut,
                );
                per_snp_kappa(sum, count, sum_sq, sum_cube, sum_4th)
            } else {
                0.0
            };
            SnpStats {
                mean,
                inv_std,
                maf,
                kurt,
            }
        })
        .collect();

    // Flat parallel-array views of stats[].mean / .inv_std for the
    // pool / rayon dispatch. See the f32 counterpart for rationale.
    let mean_buf: Vec<f64> = stats.iter().map(|s| s.mean).collect();
    let inv_std_buf: Vec<f64> = stats.iter().map(|s| s.inv_std).collect();

    for (j, s) in stats.iter().enumerate() {
        maf_out[j] = s.maf;
    }
    if let Some(kurt_slice) = kurt_out {
        for (j, s) in stats.iter().enumerate() {
            kurt_slice[j] = s.kurt;
        }
    }

    // --- Pass 2: column-parallel fused decode + normalize + scatter-add ---
    let out_ptr = out.as_mut().as_ptr_mut();
    let out_col_stride = out.col_stride();
    let bucket = &cs.bucket;
    let sign = &cs.sign_f64;

    // wasm32+atomics: dispatch through the manual SAB pool. See the
    // f32 sibling for the lifetime / safety reasoning — it applies
    // verbatim here.
    #[cfg(all(target_arch = "wasm32", target_feature = "atomics"))]
    {
        let args = crate::wasm_simd::pool::ScatterArgsF64 {
            raw_bytes_ptr: raw_bytes.as_ptr() as usize,
            bytes_per_snp,
            n_indiv,
            c,
            bucket_ptr: bucket.as_ptr() as usize,
            sign_ptr: sign.as_ptr() as usize,
            iid_pos_ptr: iid_positions.as_ptr() as usize,
            iid_pos_len: iid_positions.len(),
            all_iids: if all_iids { 1 } else { 0 },
            n_norm: n,
            out_ptr: out_ptr as usize,
            out_col_stride,
            d,
            mean_ptr: mean_buf.as_ptr() as usize,
            inv_std_ptr: inv_std_buf.as_ptr() as usize,
        };
        let dispatched = unsafe { crate::wasm_simd::pool::parallel_scatter_f64(args) };
        if dispatched {
            return;
        }
    }

    let out_send = SendPtr(out_ptr);
    let raw_ptr_us = raw_bytes.as_ptr() as usize;
    let bucket_ptr_us = bucket.as_ptr() as usize;
    let sign_ptr_us = sign.as_ptr() as usize;
    let iid_pos_ptr_us = iid_positions.as_ptr() as usize;
    let iid_pos_len = iid_positions.len();
    let mean_ptr_us = mean_buf.as_ptr() as usize;
    let inv_std_ptr_us = inv_std_buf.as_ptr() as usize;

    (0..c).into_par_iter().for_each(|j| unsafe {
        crate::wasm_simd::scatter::scatter_one_column_f64(
            j,
            raw_ptr_us as *const u8,
            bytes_per_snp,
            n_indiv,
            bucket_ptr_us as *const u32,
            sign_ptr_us as *const f64,
            iid_pos_ptr_us as *const IidPos,
            iid_pos_len,
            all_iids,
            *(mean_ptr_us as *const f64).add(j),
            *(inv_std_ptr_us as *const f64).add(j),
            n,
            out_send.ptr(),
            out_col_stride,
            d,
        );
    });
}

/// Unbiased r² estimator: r² − (1−r²)/(n−2) [Bulik-Sullivan 2015].
/// Inlined as `r*r * r2u_a + r2u_b` in hot path for constant n.
#[inline]
#[allow(dead_code)]
pub(super) fn r2_unbiased(r: f64, n: usize) -> f64 {
    let sq = r * r;
    let denom = if n > 2 { n as f64 - 2.0 } else { n as f64 };
    sq - (1.0 - sq) / denom
}

/// Quadratic inversion of the renormalized-cosine bias for CountSketch mode.
///
/// Given the sketched-and-renormalized squared cosine `tilde_rsq = (X̃_j·X̃_k)²/(||X̃_j||²·||X̃_k||²)`,
/// the leading-order bias (from a Taylor expansion of A²/(BC) around the mean) is
/// `E[tilde_rsq] - r² ≈ (1 - r²)(1 - 2r²)/d`. Inverting this as a quadratic in r²:
///   `2r⁴ + (d - 3) r² - (d·tilde_rsq - 1) = 0`
///   `r² = (-(d-3) + √((d-3)² + 8(d·tilde_rsq - 1))) / 4`
///
/// Residual bias is O(1/d²) — ~200× tighter than the linear correction at d=200.
/// See `docs/countsketch-math-analysis.md` for the full derivation and validation.
///
/// After the sketch correction, applies the standard LDSC unbiased r² correction
/// via the precomputed `r2u_a = 1 + 1/(N-2)` and `r2u_b = -1/(N-2)`:
///   `r²_unbiased = r²_sketch_corrected × r2u_a + r2u_b`
///
/// `val` is the un-normalized inner product `X̃_j·X̃_k`; `n_inv_sq = 1/N²` converts
/// it to `tilde_rsq`. `d_f` is the sketch dimension. Fallback to 0.0 in the
/// (vanishingly rare) case of a negative discriminant from extreme sketch noise.
#[inline(always)]
pub(super) fn quadratic_sketch_correction(
    val: f64,
    n_inv_sq: f64,
    d_f: f64,
    r2u_a: f64,
    r2u_b: f64,
) -> f64 {
    let tilde_rsq = val * val * n_inv_sq;
    let dm3 = d_f - 3.0;
    let disc = dm3 * dm3 + 8.0 * (d_f * tilde_rsq - 1.0);
    let rsq_sketch = if disc >= 0.0 {
        (-dm3 + disc.sqrt()) * 0.25
    } else {
        0.0
    };
    rsq_sketch * r2u_a + r2u_b
}

/// MAF-aware variant of `quadratic_sketch_correction`. The asymptotic Taylor
/// derivation dropped a subleading term in `Var(B), Var(C)` involving the
/// per-SNP kurtosis `κ_j = K_jj/N² = E[X_j⁴]/N`. Restoring it gives:
///
///   E[r̃²] − r² ≈ (1−r²)(1−2r²)/d − 2r²(κ_j + κ_k)/d
///
/// which inverts to the same quadratic with the linear coefficient shifted:
///
///   2r⁴ + (d − 3 − 2(κ_j + κ_k)) r² − (d·tilde_rsq − 1) = 0
///
/// Active only when `--sketch-maf-aware` is set. At biobank scale (N=50K) the
/// kurtosis term is O(1/N) ≈ 10⁻⁵ to 10⁻³ depending on MAF, so the practical
/// correction is small; the benchmark cross-product quantifies whether it's
/// measurable. See `docs/countsketch-math-analysis.md` §14.
#[inline(always)]
pub(super) fn quadratic_sketch_correction_maf_aware(
    val: f64,
    n_inv_sq: f64,
    d_f: f64,
    r2u_a: f64,
    r2u_b: f64,
    kurt_j: f64,
    kurt_k: f64,
) -> f64 {
    let tilde_rsq = val * val * n_inv_sq;
    let dm3 = d_f - 3.0 - 2.0 * (kurt_j + kurt_k);
    let disc = dm3 * dm3 + 8.0 * (d_f * tilde_rsq - 1.0);
    let rsq_sketch = if disc >= 0.0 {
        (-dm3 + disc.sqrt()) * 0.25
    } else {
        0.0
    };
    rsq_sketch * r2u_a + r2u_b
}

/// Compute LD scores for all SNPs. Returns `(l2, maf)`.
///
/// The caller is responsible for opening the `Bed` (via
/// [`Bed::builder`] for a file or [`Bed::from_bytes`] for an
/// in-memory buffer) and, on native targets that want zero-copy I/O,
/// optionally pre-opening a matching [`MmapBed`]. The `Bed` is
/// consumed; if the caller wants to compute on multiple chromosomes
/// independently it must open one `Bed` per chromosome.
#[allow(clippy::too_many_arguments, clippy::unnecessary_cast)]
pub(super) fn compute_ldscore_global(
    all_snps: &[BimRecord],
    bed: Bed,
    mmap_bed: Option<MmapBed>,
    n_indiv: usize,
    mode: WindowMode,
    chunk_c: usize,
    annot: Option<&MatF>,
    iid_indices: Option<&[isize]>,
    pq_exp: Option<f64>,
    yes_really: bool,
    #[cfg(feature = "gpu")] gpu_ctx: Option<&GpuContext>,
    #[cfg(feature = "gpu")] gpu_config: GpuConfig,
    use_f32: bool,
    verbose_timing: bool,
    sketch: Option<usize>,
    sketch_maf_aware: bool,
    snp_level_masking: bool,
    on_progress: &mut dyn FnMut(crate::l2::L2Progress),
) -> Result<(MatF, Vec<f64>, crate::l2::L2Perf)> {
    // chunk_c=0 makes `(0..m).step_by(chunk_c)` panic ("step must be
    // non-zero"). Catch the malformed-config case here with a real
    // error instead of a wasm `unreachable executed` trap. The
    // `chunks_total = m.div_ceil(chunk_c.max(1))` line below would
    // otherwise mask this with a defensive-looking `.max(1)` that
    // doesn't actually defend the loop.
    anyhow::ensure!(
        chunk_c >= 1,
        "compute_ldscore_global: chunk_size must be >= 1, got {}",
        chunk_c,
    );
    // `use_mmap` is derived from whether the caller pre-opened a
    // memory-mapped view; that's the only meaningful signal at this
    // layer. Browser callers always pass `mmap_bed = None`.
    let use_mmap = mmap_bed.is_some();
    let m = all_snps.len();
    if m == 0 {
        let n_annot = annot.map(|a| a.ncols()).unwrap_or(1);
        return Ok((
            MatF::zeros(0, n_annot),
            vec![],
            crate::l2::L2Perf {
                m: 0,
                n_indiv,
                ..Default::default()
            },
        ));
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

    // ── Randomized sketching setup ──────────────────────────────────────────
    // Compress individual dimension N→d via random projection P (d×N).
    // All subsequent GEMMs operate on d-dimensional data.
    // Bias correction: E[val̃²] = val²(1+1/d) + N²/d, absorbed into adjusted constants.
    let sketch_dim: Option<usize> = match sketch {
        Some(d) => {
            anyhow::ensure!(
                d >= 3 && d < n_indiv,
                "--sketch {d} is out of bounds (must be 3 ≤ d < n_indiv={n_indiv})"
            );
            Some(d)
        }
        None => None,
    };
    let gemm_n: usize; // row dimension for all GEMM matrices (d when sketching, N otherwise)
    // Constants for the non-sketch path: hot loop is `val² × active_n_inv_sq_r2u_a
    // + active_r2u_b`, which fuses the val→r̂ conversion (val·n_inv) with the
    // LDSC unbiased r² correction `r̂² + (-(1-r̂²)/(N-2))`.
    let active_n_inv_sq_r2u_a: f64;
    let active_r2u_b: f64;
    // Constants for the sketch path: quadratic inversion of the renormalized-cosine
    // bias, followed by the same LDSC unbiased correction. `quad_d_f` is the
    // sketch dimension and `n_inv_sq = 1/N²` converts val² → tilde_rsq.
    let quad_d_f: f64;
    let n_inv_sq: f64 = n_inv * n_inv;

    if let Some(d) = sketch_dim {
        gemm_n = d;
        // Bias correction for the sketched, renormalized squared cosine.
        // After CountSketch projection, each column is rescaled so ||x̃'_j||² = N
        // (the ratio estimator). The resulting squared cosine
        //   r̃'² = (X̃_j^T X̃_k)² / (||X̃_j||² · ||X̃_k||²)
        // has bias `E[r̃'²] - r² ≈ (1 - r²)(1 - 2r²)/d` (Taylor expansion of
        // A²/(BC) around its mean; validated empirically to Monte Carlo noise).
        //
        // We invert this exactly: solve `2r⁴ + (d-3)r² - (d·r̃'² - 1) = 0`
        // for r², giving residual bias O(1/d²). The LDSC unbiased r² correction
        // is then applied on top. See `docs/countsketch-math-analysis.md` for
        // the full derivation, MC validation, and a per-bin analysis on chr22
        // EUR showing the correction shifts high-LD-score SNPs in exactly the
        // predicted direction.
        //
        // The non-sketch constants below are unused in this branch (the
        // quadratic path uses `n_inv_sq`, `quad_d_f`, `r2u_a`, `r2u_b` directly).
        quad_d_f = d as f64;
        active_n_inv_sq_r2u_a = n_inv_sq_r2u_a;
        active_r2u_b = r2u_b;
        if d <= 50 {
            eprintln!(
                "WARNING: --sketch {d} is below the recommended minimum (d ≥ 100). \
                 At d ≤ 50, higher-order bias terms exceed O(1/d) and the quadratic \
                 correction's sqrt amplifies sketch noise; the corrected LD scores \
                 are no better than uncorrected at the LD-score level. Increase d \
                 (cost is essentially flat in d below the GEMM crossover)."
            );
        }
        println!(
            "Sketch mode: projecting N={n_indiv} → d={d} (CountSketch, quadratic bias correction)"
        );
    } else {
        gemm_n = n_indiv;
        active_n_inv_sq_r2u_a = n_inv_sq_r2u_a;
        active_r2u_b = r2u_b;
        quad_d_f = 0.0;
    }

    // CountSketch projection: O(N×c) scatter-add, cost flat in d.
    let count_sketch: Option<CountSketchProj> =
        sketch_dim.map(|d| CountSketchProj::new(n_indiv, d, 42));

    // Fused BED-decode-normalize-CountSketch: single pass over raw BED bytes,
    // scatter-adds directly into the d×c sketch buffer. No intermediate N×c matrix.
    let use_fused_sketch = sketch_dim.is_some();

    let mut maf_per_snp = vec![0.0f64; m];
    // Per-SNP kurtosis κ_j = K_jj/N² (Σ_i X_ij⁴ / N²) for MAF-aware bias correction.
    // Allocated only when --sketch-maf-aware is active; populated in the chunk loop
    // stats pass and consumed at the fill_r2u sites. See compute.rs::per_snp_kappa
    // and docs/countsketch-math-analysis.md §14.
    let mut kurt_per_snp: Vec<f64> = if sketch_maf_aware && sketch_dim.is_some() {
        vec![0.0f64; m]
    } else {
        Vec::new()
    };

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
    let mut bufs = GemmBufs::new(use_f32, gemm_n, chunk_c, ring_size, max_window_size);
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

    // Sketch: temporary full-size buffer for normalization before projection.
    // Fused CountSketch bypasses b_full for contiguous chunks, but non-contiguous
    // chunks (--extract with gaps) fall back to the standard decode+project path
    // which requires b_full. Allocate b_full whenever sketch is active.
    let need_b_full = sketch_dim.is_some();
    let mut b_full_f64: MatF = if need_b_full && !use_f32 {
        MatF::zeros(n_indiv, chunk_c)
    } else {
        MatF::zeros(0, 0)
    };
    let mut b_full_f32: MatF32 = if need_b_full && use_f32 {
        MatF32::zeros(n_indiv, chunk_c)
    } else {
        MatF32::zeros(0, 0)
    };

    let mut r2u_bb = MatF::zeros(chunk_c, chunk_c);
    let mut r2u_ab = MatF::zeros(max_window_size.max(1), chunk_c);

    // Pre-allocate scratch matrices used per chunk (avoids alloc+zero every iteration).
    // All are written with Accum::Replace (full overwrite) so stale data is harmless.
    let mut bb_f64 = if !use_f32 {
        MatF::zeros(chunk_c, chunk_c)
    } else {
        MatF::zeros(0, 0)
    };
    let mut bb_f32 = if use_f32 {
        MatF32::zeros(chunk_c, chunk_c)
    } else {
        MatF32::zeros(0, 0)
    };
    let mut contrib_bb = MatF::zeros(chunk_c, n_annot);
    let mut contrib_left = MatF::zeros(max_window_size.max(1), n_annot);
    let mut contrib_right = MatF::zeros(chunk_c, n_annot);

    use web_time::Instant;
    let mut t_bed_read = std::time::Duration::ZERO;
    let mut t_norm = std::time::Duration::ZERO;
    let mut t_sketch = std::time::Duration::ZERO;
    let mut t_bb_dot = std::time::Duration::ZERO;
    let mut t_ab_dot = std::time::Duration::ZERO;
    let mut t_ring_store = std::time::Duration::ZERO;

    // GPU context and config are passed in from the caller (created once in run()).
    #[cfg(feature = "gpu")]
    let (_gpu_tile_cols, _gpu_flex32, _gpu_f64) =
        (gpu_config.tile_cols, gpu_config.flex32, gpu_config.f64);

    // Sequential path: use the caller-supplied BED with pre-computed
    // ChunkReader. When mmap is enabled (caller pre-opened a `MmapBed`),
    // we still need ChunkReader for metadata (iid_positions, etc.) but
    // I/O goes through MmapBed instead.
    let mut seq_bed = Some(bed);
    let mut chunk_reader: Option<ChunkReader<f32>> = {
        let bed = seq_bed.as_ref().unwrap();
        Some(
            ChunkReader::new(bed, iid_indices, chunk_c, true, f32::NAN)
                .context("creating chunk reader")?,
        )
    };
    // Check if ALL SNPs are contiguous in the BED file (common: no --extract).
    // If so, use sequential reads (1 seek + streaming reads) instead of seeking per chunk.
    // BufReader::seek(SeekFrom::Start) discards its internal buffer every time,
    // causing 64× read amplification (8MB buffer filled, 124KB used, buffer discarded).
    let all_bed_sequential = m > 0
        && all_snps
            .windows(2)
            .all(|w| w[1].bed_idx == w[0].bed_idx + 1);
    if all_bed_sequential && !use_mmap {
        let bed = seq_bed.as_mut().unwrap();
        let cr = chunk_reader.as_mut().unwrap();
        cr.start_sequential(bed, all_snps[0].bed_idx)
            .context("starting sequential BED read")?;
    }

    let chunks_total = m.div_ceil(chunk_c.max(1));
    let mut chunks_done = 0usize;
    for chunk_start in (0..m).step_by(chunk_c) {
        let chunk_started = Instant::now();
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
        // When using the fused sketch path (sequential reader + CountSketch), read raw bytes
        // and skip b_full materialization. Otherwise use the existing decoded-Mat path.
        let sketching = sketch_dim.is_some();
        // Determine if this chunk qualifies for the fused path.
        // Fused requires: (a) contiguous chunk, (b) not prefetch.
        // With mmap: any contiguous chunk qualifies (no sequential mode needed).
        let chunk_is_contiguous = if use_mmap {
            // mmap can access any contiguous range without seek penalty
            let start_bed_idx = all_snps[chunk_start].bed_idx;
            all_snps[chunk_end - 1].bed_idx == start_bed_idx + c - 1
        } else if !all_bed_sequential {
            let start_bed_idx = all_snps[chunk_start].bed_idx;
            all_snps[chunk_end - 1].bed_idx == start_bed_idx + c - 1
        } else {
            all_bed_sequential
        };
        let do_fused = use_fused_sketch && chunk_is_contiguous;

        if do_fused {
            // ── Fused path: read raw bytes, fuse decode+normalize+project ──────────
            // Extract metadata from ChunkReader (immutable borrow)
            let cr_ref = chunk_reader.as_ref().unwrap();
            let bps = cr_ref.bytes_per_snp_val();
            let start_bed_idx = all_snps[chunk_start].bed_idx;

            let raw: std::borrow::Cow<'_, [u8]> = if let Some(ref mb) = mmap_bed {
                // Zero-copy: slice directly from memory-mapped file
                let raw_slice = mb.snp_bytes(start_bed_idx, c);
                // Prefetch next chunk asynchronously
                #[cfg(unix)]
                {
                    let next_start = chunk_start + chunk_c;
                    if next_start < m {
                        let next_end = (next_start + chunk_c).min(m);
                        let next_c = next_end - next_start;
                        let next_bed_idx = all_snps[next_start].bed_idx;
                        mb.advise_willneed(next_bed_idx, next_c);
                    }
                }
                std::borrow::Cow::Borrowed(&raw_slice[..c * bps])
            } else {
                // BufReader path: read into ChunkReader buffer, then copy
                let cr = chunk_reader.as_mut().unwrap();
                let bed = seq_bed.as_mut().unwrap();
                if all_bed_sequential {
                    cr.read_next_raw_only(bed, c).with_context(|| {
                        format!("reading BED chunk [{},{})", chunk_start, chunk_end)
                    })?;
                } else {
                    cr.read_contiguous_raw_only(bed, start_bed_idx, c)
                        .with_context(|| {
                            format!("reading BED chunk [{},{})", chunk_start, chunk_end)
                        })?;
                }
                std::borrow::Cow::Owned(cr.raw_bytes()[..c * bps].to_vec())
            };
            t_bed_read += t0.elapsed();

            let t0 = Instant::now();
            let cr_ref = chunk_reader.as_ref().unwrap();
            let n_i = cr_ref.n_iid();
            let ipos = cr_ref.iid_positions_ref().to_vec();
            let ai = cr_ref.all_iids_flag();
            if let Some(ref cs) = count_sketch {
                // Fused CountSketch: decode + normalize + scatter-add (no GEMM)
                // Split borrows: maf and kurt are disjoint slices of the global vec(s)
                let kurt_slice: Option<&mut [f64]> = if sketch_maf_aware {
                    Some(&mut kurt_per_snp[chunk_start..chunk_start + c])
                } else {
                    None
                };
                match bufs {
                    GemmBufs::F32 { ref mut b_mat, .. } => {
                        countsketch_fused_project_f32(
                            &raw,
                            bps,
                            n_i,
                            c,
                            cs,
                            n,
                            &ipos,
                            ai,
                            &mut maf_per_snp[chunk_start..chunk_start + c],
                            kurt_slice,
                            b_mat,
                        );
                    }
                    GemmBufs::F64 { ref mut b_mat, .. } => {
                        countsketch_fused_project_f64(
                            &raw,
                            bps,
                            n_i,
                            c,
                            cs,
                            n,
                            &ipos,
                            ai,
                            &mut maf_per_snp[chunk_start..chunk_start + c],
                            kurt_slice,
                            b_mat,
                        );
                    }
                }
            }
            // pq_exp: derive from the MAF we just computed
            if let (Some(exp), Some(ref mut pq_vec)) = (pq_exp, pq_per_snp.as_mut()) {
                for j in 0..c {
                    let maf = maf_per_snp[chunk_start + j];
                    pq_vec[chunk_start + j] = (maf * (1.0 - maf)).powf(exp);
                }
            }
            t_sketch += t0.elapsed();
        } else {
            // ── Standard path: read decoded Mat<f32>, normalize, then project ────────
            let prefetch_raw: MatF32;
            let raw_ref: &MatF32 = if all_bed_sequential {
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
                    cr.read_contiguous(bed, start_bed_idx, c).with_context(|| {
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
            // When sketching, normalize into b_full (n_indiv rows), then project into b_mat (d rows).
            // When not sketching, normalize directly into b_mat (n_indiv rows).
            for j in 0..c {
                // Get raw column as contiguous slice (1 bounds check vs 2 per element)
                let raw_col = raw_ref.col(j).try_as_col_major().unwrap().as_slice();
                match bufs {
                    GemmBufs::F32 { ref mut b_mat, .. } => {
                        let col = if sketching {
                            b_full_f32
                                .col_mut(j)
                                .try_as_col_major_mut()
                                .unwrap()
                                .as_slice_mut()
                        } else {
                            b_mat
                                .col_mut(j)
                                .try_as_col_major_mut()
                                .unwrap()
                                .as_slice_mut()
                        };
                        col[..n_indiv].copy_from_slice(&raw_col[..n_indiv]);
                        // SAFETY: we build with +avx2,+fma for the musl target;
                        // native builds have AVX2 on any modern x86_64 CPU.
                        let (sum, sum_sq) = unsafe { sum_sumsq_f32(&raw_col[..n_indiv]) };
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
                        // MAF-aware kurtosis: compute κ_j before normalization
                        // (raw_col still holds 0/1/2/NaN). Cheap O(N) scan.
                        if sketch_maf_aware && sketching {
                            let mean_f64 = if count > 0 { sum / count as f64 } else { 0.0 };
                            let centered_ss = sum_sq - count as f64 * mean_f64 * mean_f64;
                            kurt_per_snp[chunk_start + j] = per_snp_kappa_from_col_f32(
                                &raw_col[..n_indiv],
                                mean_f64,
                                centered_ss,
                            );
                        }
                        let snp_maf =
                            normalize_col_f32_with_stats(col, n_indiv, sum, count, sum_sq);
                        maf_per_snp[chunk_start + j] = snp_maf as f64;
                        if let (Some(exp), Some(ref mut pq_vec)) = (pq_exp, pq_per_snp.as_mut()) {
                            let maf_f64 = snp_maf as f64;
                            let pq = (maf_f64 * (1.0 - maf_f64)).powf(exp);
                            pq_vec[chunk_start + j] = pq;
                        }
                    }
                    GemmBufs::F64 { ref mut b_mat, .. } => {
                        let col = if sketching {
                            b_full_f64
                                .col_mut(j)
                                .try_as_col_major_mut()
                                .unwrap()
                                .as_slice_mut()
                        } else {
                            b_mat
                                .col_mut(j)
                                .try_as_col_major_mut()
                                .unwrap()
                                .as_slice_mut()
                        };
                        // SAFETY: see above.
                        let (sum, sum_sq) = unsafe { sum_sumsq_f32(&raw_col[..n_indiv]) };
                        // Widening copy (f32 → f64) — vectorizes to vcvtps2pd + vmovupd
                        for (dst, &src) in col[..n_indiv].iter_mut().zip(raw_col[..n_indiv].iter())
                        {
                            *dst = src as f64;
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
                        // MAF-aware kurtosis: compute κ_j before normalization
                        // (col still holds 0/1/2/NaN as f64 from the widening copy).
                        if sketch_maf_aware && sketching {
                            let mean_f64 = if count > 0 { sum / count as f64 } else { 0.0 };
                            let centered_ss = sum_sq - count as f64 * mean_f64 * mean_f64;
                            kurt_per_snp[chunk_start + j] =
                                per_snp_kappa_from_col_f64(&col[..n_indiv], mean_f64, centered_ss);
                        }
                        let snp_maf =
                            normalize_col_f64_with_stats(col, n_indiv, sum, count, sum_sq);
                        maf_per_snp[chunk_start + j] = snp_maf;
                        if let (Some(exp), Some(ref mut pq_vec)) = (pq_exp, pq_per_snp.as_mut()) {
                            let pq = (snp_maf * (1.0 - snp_maf)).powf(exp);
                            pq_vec[chunk_start + j] = pq;
                        }
                    }
                }
            }

            t_norm += t0.elapsed();

            // Sketch projection: b_mat (d × c) = CountSketch(b_full)
            // O(N×c) scatter-add, flat cost in d.
            if sketching {
                let t0 = Instant::now();
                let cs = count_sketch.as_ref().unwrap();
                match bufs {
                    GemmBufs::F32 { ref mut b_mat, .. } => {
                        cs.project_f32(&b_full_f32, n_indiv, c, n, b_mat);
                    }
                    GemmBufs::F64 { ref mut b_mat, .. } => {
                        cs.project_f64(&b_full_f64, n_indiv, c, n, b_mat);
                    }
                }
                t_sketch += t0.elapsed();
            }
        } // end standard path else (do_fused was false)

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

        // ---- Exact GEMM path ----

        // B×B — compute GEMM then extract r2_unbiased into r2u_bb (always f64).
        let t0 = Instant::now();
        if !chunk_all_zero {
            // Inline helper: fill r2u_bb from bb matrix values.
            // Generic over closure to avoid dyn dispatch in inner loop.
            // The `apply_correction` closure encodes the per-pair correction
            // (fused linear or quadratic inversion + LDSC unbiased). Indices
            // `(k, j)` are passed through so the MAF-aware correction can look
            // up per-SNP κ for both members of the pair.
            #[inline(always)]
            fn fill_r2u_bb_from<F: Fn(usize, usize) -> f64, G: Fn(usize, usize, f64) -> f64>(
                r2u_bb: &mut MatF,
                c: usize,
                bb_val: F,
                apply_correction: G,
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
                        let r2u = apply_correction(k, j, val);
                        r2u_bb[(j, k)] = r2u;
                        r2u_bb[(k, j)] = r2u;
                    }
                }
            }

            match bufs {
                GemmBufs::F32 { ref b_mat, .. } => {
                    let b_slice = mat_slice_f32(b_mat.as_ref(), 0..gemm_n, 0..c);
                    #[allow(unused_mut)]
                    let mut bb_sl = mat_slice_mut_f32(bb_f32.as_mut(), 0..c, 0..c);
                    let mut _did_gpu = false;
                    #[cfg(feature = "gpu")]
                    {
                        if let Some(ctx) = gpu_ctx {
                            let b_f32 = mat_to_col_major_f32_from_f32(b_slice, gemm_n, c);
                            let gpu_result = if let Some(tc) = _gpu_tile_cols {
                                if _gpu_flex32 {
                                    ctx.matmul_tn_tiled_flex32(&b_f32, gemm_n, c, &b_f32, c, tc)
                                } else {
                                    ctx.matmul_tn_tiled(&b_f32, gemm_n, c, &b_f32, c, tc)
                                }
                            } else if _gpu_flex32 {
                                ctx.matmul_tn_flex32(&b_f32, gemm_n, c, &b_f32, c)
                            } else {
                                ctx.matmul_tn(&b_f32, gemm_n, c, &b_f32, c)
                            };
                            match gpu_result {
                                Ok(result) => {
                                    write_gpu_result_f32(&result, bb_sl.as_mut(), c, c);
                                    _did_gpu = true;
                                }
                                Err(e) => eprintln!(
                                    "GPU: f32 B×B matmul failed ({e}), falling back to CPU"
                                ),
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
                            par_default(),
                        );
                    }
                    if sketch_dim.is_some() {
                        if sketch_maf_aware {
                            let kp = &kurt_per_snp;
                            fill_r2u_bb_from(
                                &mut r2u_bb,
                                c,
                                |k, j| bb_f32[(k, j)] as f64,
                                |k, j, val| {
                                    quadratic_sketch_correction_maf_aware(
                                        val,
                                        n_inv_sq,
                                        quad_d_f,
                                        r2u_a,
                                        r2u_b,
                                        kp[chunk_start + j],
                                        kp[chunk_start + k],
                                    )
                                },
                            );
                        } else {
                            fill_r2u_bb_from(
                                &mut r2u_bb,
                                c,
                                |k, j| bb_f32[(k, j)] as f64,
                                |_k, _j, val| {
                                    quadratic_sketch_correction(
                                        val, n_inv_sq, quad_d_f, r2u_a, r2u_b,
                                    )
                                },
                            );
                        }
                    } else {
                        fill_r2u_bb_from(
                            &mut r2u_bb,
                            c,
                            |k, j| bb_f32[(k, j)] as f64,
                            |_k, _j, val| val * val * active_n_inv_sq_r2u_a + active_r2u_b,
                        );
                    }
                }
                GemmBufs::F64 { ref b_mat, .. } => {
                    let b_slice = mat_slice(b_mat.as_ref(), 0..gemm_n, 0..c);
                    #[allow(unused_mut)]
                    let mut bb_sl = mat_slice_mut(bb_f64.as_mut(), 0..c, 0..c);
                    let mut _did_gpu = false;
                    #[cfg(feature = "gpu")]
                    {
                        if let Some(ctx) = gpu_ctx {
                            if _gpu_f64 && ctx.capabilities.has_f64 {
                                // Native f64 path — no precision loss
                                let b_f64 = mat_to_col_major_f64(b_slice, gemm_n, c);
                                let gpu_result = if let Some(tc) = _gpu_tile_cols {
                                    ctx.matmul_tn_tiled_f64(&b_f64, gemm_n, c, &b_f64, c, tc)
                                } else {
                                    ctx.matmul_tn_f64(&b_f64, gemm_n, c, &b_f64, c)
                                };
                                match gpu_result {
                                    Ok(result) => {
                                        write_gpu_result_f64_native(&result, bb_sl.as_mut(), c, c);
                                        _did_gpu = true;
                                    }
                                    Err(e) => eprintln!(
                                        "GPU: f64 B×B matmul failed ({e}), falling back to CPU"
                                    ),
                                }
                            } else {
                                // f32 conversion path
                                let b_f32 = mat_to_col_major_f32_from_f64(b_slice, gemm_n, c);
                                let gpu_result = if let Some(tc) = _gpu_tile_cols {
                                    if _gpu_flex32 {
                                        ctx.matmul_tn_tiled_flex32(&b_f32, gemm_n, c, &b_f32, c, tc)
                                    } else {
                                        ctx.matmul_tn_tiled(&b_f32, gemm_n, c, &b_f32, c, tc)
                                    }
                                } else if _gpu_flex32 {
                                    ctx.matmul_tn_flex32(&b_f32, gemm_n, c, &b_f32, c)
                                } else {
                                    ctx.matmul_tn(&b_f32, gemm_n, c, &b_f32, c)
                                };
                                match gpu_result {
                                    Ok(result) => {
                                        write_gpu_result_f64(&result, bb_sl.as_mut(), c, c);
                                        _did_gpu = true;
                                    }
                                    Err(e) => eprintln!(
                                        "GPU: B×B matmul failed ({e}), falling back to CPU"
                                    ),
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
                            par_default(),
                        );
                    }
                    if sketch_dim.is_some() {
                        if sketch_maf_aware {
                            let kp = &kurt_per_snp;
                            fill_r2u_bb_from(
                                &mut r2u_bb,
                                c,
                                |k, j| bb_f64[(k, j)],
                                |k, j, val| {
                                    quadratic_sketch_correction_maf_aware(
                                        val,
                                        n_inv_sq,
                                        quad_d_f,
                                        r2u_a,
                                        r2u_b,
                                        kp[chunk_start + j],
                                        kp[chunk_start + k],
                                    )
                                },
                            );
                        } else {
                            fill_r2u_bb_from(
                                &mut r2u_bb,
                                c,
                                |k, j| bb_f64[(k, j)],
                                |_k, _j, val| {
                                    quadratic_sketch_correction(
                                        val, n_inv_sq, quad_d_f, r2u_a, r2u_b,
                                    )
                                },
                            );
                        }
                    } else {
                        fill_r2u_bb_from(
                            &mut r2u_bb,
                            c,
                            |k, j| bb_f64[(k, j)],
                            |_k, _j, val| val * val * active_n_inv_sq_r2u_a + active_r2u_b,
                        );
                    }
                }
            }
            // SNP-level B×B masking: zero r2u_bb entries for pairs outside
            // each other's window. Only triggers at chromosome boundaries
            // where block_left jumps within a chunk.
            if snp_level_masking {
                for j in 1..c {
                    let bl_j = block_left[chunk_start + j];
                    if bl_j <= chunk_start {
                        continue;
                    }
                    let cutoff = bl_j - chunk_start;
                    for k in 0..cutoff.min(j) {
                        r2u_bb[(j, k)] = 0.0;
                        r2u_bb[(k, j)] = 0.0;
                    }
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
                                    .submatrix_mut(0, wi, gemm_n, 1)
                                    .copy_from(ring_buf.as_ref().submatrix(0, *slot, gemm_n, 1));
                            }
                        }
                        let a_view = if contiguous {
                            mat_slice_f32(
                                ring_buf.as_ref(),
                                0..gemm_n,
                                first_slot..(first_slot + w),
                            )
                        } else {
                            mat_slice_f32(a_buf.as_ref(), 0..gemm_n, 0..w)
                        };
                        let b_sl = mat_slice_f32(b_mat.as_ref(), 0..gemm_n, 0..c);

                        let mut _did_gpu = false;
                        #[cfg(feature = "gpu")]
                        {
                            if let Some(ctx) = gpu_ctx {
                                let a_f32 = mat_to_col_major_f32_from_f32(a_view, gemm_n, w);
                                let b_f32 = mat_to_col_major_f32_from_f32(b_sl, gemm_n, c);
                                let gpu_result = if let Some(tc) = _gpu_tile_cols {
                                    if _gpu_flex32 {
                                        ctx.matmul_tn_tiled_flex32(&a_f32, gemm_n, w, &b_f32, c, tc)
                                    } else {
                                        ctx.matmul_tn_tiled(&a_f32, gemm_n, w, &b_f32, c, tc)
                                    }
                                } else if _gpu_flex32 {
                                    ctx.matmul_tn_flex32(&a_f32, gemm_n, w, &b_f32, c)
                                } else {
                                    ctx.matmul_tn(&a_f32, gemm_n, w, &b_f32, c)
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
                                    Err(e) => eprintln!(
                                        "GPU: f32 A×B matmul failed ({e}), falling back to CPU"
                                    ),
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
                                par_default(),
                            );
                        }
                        let ab_view = mat_slice_f32(ab_buf.as_ref(), 0..w, 0..c);
                        // Column-major inner loop with column slices: eliminates faer
                        // tuple-indexing bounds checks (~8.3B checks for full genome).
                        // Branch on sketch mode once outside the inner loop.
                        if sketch_dim.is_some() && sketch_maf_aware {
                            // Precompute per-window-row κ_wi for cache-friendliness.
                            let kurt_a: Vec<f64> = window
                                .iter()
                                .map(|&(g_idx, _)| kurt_per_snp[g_idx])
                                .collect();
                            for j in 0..c {
                                let kurt_b_j = kurt_per_snp[chunk_start + j];
                                let ab_col = ab_view.col(j).try_as_col_major().unwrap().as_slice();
                                let r2u_col = r2u_ab
                                    .col_mut(j)
                                    .try_as_col_major_mut()
                                    .unwrap()
                                    .as_slice_mut();
                                for (wi, (ab_val, r2u_val)) in
                                    ab_col[..w].iter().zip(r2u_col[..w].iter_mut()).enumerate()
                                {
                                    let val = *ab_val as f64;
                                    *r2u_val = quadratic_sketch_correction_maf_aware(
                                        val, n_inv_sq, quad_d_f, r2u_a, r2u_b, kurt_a[wi], kurt_b_j,
                                    );
                                }
                            }
                        } else if sketch_dim.is_some() {
                            for j in 0..c {
                                let ab_col = ab_view.col(j).try_as_col_major().unwrap().as_slice();
                                let r2u_col = r2u_ab
                                    .col_mut(j)
                                    .try_as_col_major_mut()
                                    .unwrap()
                                    .as_slice_mut();
                                for (ab_val, r2u_val) in
                                    ab_col[..w].iter().zip(r2u_col[..w].iter_mut())
                                {
                                    let val = *ab_val as f64;
                                    *r2u_val = quadratic_sketch_correction(
                                        val, n_inv_sq, quad_d_f, r2u_a, r2u_b,
                                    );
                                }
                            }
                        } else {
                            for j in 0..c {
                                let ab_col = ab_view.col(j).try_as_col_major().unwrap().as_slice();
                                let r2u_col = r2u_ab
                                    .col_mut(j)
                                    .try_as_col_major_mut()
                                    .unwrap()
                                    .as_slice_mut();
                                for (ab_val, r2u_val) in
                                    ab_col[..w].iter().zip(r2u_col[..w].iter_mut())
                                {
                                    let val = *ab_val as f64;
                                    *r2u_val = val * val * active_n_inv_sq_r2u_a + active_r2u_b;
                                }
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
                                    .submatrix_mut(0, wi, gemm_n, 1)
                                    .copy_from(ring_buf.as_ref().submatrix(0, *slot, gemm_n, 1));
                            }
                        }
                        let a_view = if contiguous {
                            mat_slice(ring_buf.as_ref(), 0..gemm_n, first_slot..(first_slot + w))
                        } else {
                            mat_slice(a_buf.as_ref(), 0..gemm_n, 0..w)
                        };
                        let b_sl = mat_slice(b_mat.as_ref(), 0..gemm_n, 0..c);

                        let mut _did_gpu = false;
                        #[cfg(feature = "gpu")]
                        {
                            if let Some(ctx) = gpu_ctx {
                                if _gpu_f64 && ctx.capabilities.has_f64 {
                                    // Native f64 path
                                    let a_f64_data = mat_to_col_major_f64(a_view, gemm_n, w);
                                    let b_f64_data = mat_to_col_major_f64(b_sl, gemm_n, c);
                                    let gpu_result = if let Some(tc) = _gpu_tile_cols {
                                        ctx.matmul_tn_tiled_f64(
                                            &a_f64_data,
                                            gemm_n,
                                            w,
                                            &b_f64_data,
                                            c,
                                            tc,
                                        )
                                    } else {
                                        ctx.matmul_tn_f64(&a_f64_data, gemm_n, w, &b_f64_data, c)
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
                                        Err(e) => eprintln!(
                                            "GPU: f64 A×B matmul failed ({e}), falling back to CPU"
                                        ),
                                    }
                                } else {
                                    // f32 conversion path
                                    let a_f32 = mat_to_col_major_f32_from_f64(a_view, gemm_n, w);
                                    let b_f32 = mat_to_col_major_f32_from_f64(b_sl, gemm_n, c);
                                    let gpu_result = if let Some(tc) = _gpu_tile_cols {
                                        if _gpu_flex32 {
                                            ctx.matmul_tn_tiled_flex32(
                                                &a_f32, gemm_n, w, &b_f32, c, tc,
                                            )
                                        } else {
                                            ctx.matmul_tn_tiled(&a_f32, gemm_n, w, &b_f32, c, tc)
                                        }
                                    } else if _gpu_flex32 {
                                        ctx.matmul_tn_flex32(&a_f32, gemm_n, w, &b_f32, c)
                                    } else {
                                        ctx.matmul_tn(&a_f32, gemm_n, w, &b_f32, c)
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
                                        Err(e) => eprintln!(
                                            "GPU: A×B matmul failed ({e}), falling back to CPU"
                                        ),
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
                                par_default(),
                            );
                        }
                        let ab_view = mat_slice(ab_buf.as_ref(), 0..w, 0..c);
                        // Column-major inner loop with column slices: eliminates faer
                        // tuple-indexing bounds checks. Branch on sketch mode once.
                        if sketch_dim.is_some() && sketch_maf_aware {
                            let kurt_a: Vec<f64> = window
                                .iter()
                                .map(|&(g_idx, _)| kurt_per_snp[g_idx])
                                .collect();
                            for j in 0..c {
                                let kurt_b_j = kurt_per_snp[chunk_start + j];
                                let ab_col = ab_view.col(j).try_as_col_major().unwrap().as_slice();
                                let r2u_col = r2u_ab
                                    .col_mut(j)
                                    .try_as_col_major_mut()
                                    .unwrap()
                                    .as_slice_mut();
                                for (wi, (ab_val, r2u_val)) in
                                    ab_col[..w].iter().zip(r2u_col[..w].iter_mut()).enumerate()
                                {
                                    *r2u_val = quadratic_sketch_correction_maf_aware(
                                        *ab_val, n_inv_sq, quad_d_f, r2u_a, r2u_b, kurt_a[wi],
                                        kurt_b_j,
                                    );
                                }
                            }
                        } else if sketch_dim.is_some() {
                            for j in 0..c {
                                let ab_col = ab_view.col(j).try_as_col_major().unwrap().as_slice();
                                let r2u_col = r2u_ab
                                    .col_mut(j)
                                    .try_as_col_major_mut()
                                    .unwrap()
                                    .as_slice_mut();
                                for (ab_val, r2u_val) in
                                    ab_col[..w].iter().zip(r2u_col[..w].iter_mut())
                                {
                                    *r2u_val = quadratic_sketch_correction(
                                        *ab_val, n_inv_sq, quad_d_f, r2u_a, r2u_b,
                                    );
                                }
                            }
                        } else {
                            for j in 0..c {
                                let ab_col = ab_view.col(j).try_as_col_major().unwrap().as_slice();
                                let r2u_col = r2u_ab
                                    .col_mut(j)
                                    .try_as_col_major_mut()
                                    .unwrap()
                                    .as_slice_mut();
                                for (ab_val, r2u_val) in
                                    ab_col[..w].iter().zip(r2u_col[..w].iter_mut())
                                {
                                    *r2u_val =
                                        *ab_val * *ab_val * active_n_inv_sq_r2u_a + active_r2u_b;
                                }
                            }
                        }
                    }
                }

                // SNP-level A×B masking: zero r2u_ab entries where the window
                // SNP is outside the chunk SNP's per-SNP window boundary.
                // Window is sorted by snp_idx; block_left is non-decreasing
                // across the chunk, so the cutoff row can only advance.
                if snp_level_masking {
                    let mut cutoff_row = 0usize;
                    for j in 0..c {
                        let bl_j = block_left[chunk_start + j];
                        while cutoff_row < w
                            && window.get(cutoff_row).is_some_and(|(idx, _)| *idx < bl_j)
                        {
                            cutoff_row += 1;
                        }
                        for wi in 0..cutoff_row {
                            r2u_ab[(wi, j)] = 0.0;
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
                        .submatrix_mut(0, start_slot, gemm_n, c)
                        .copy_from(b_mat.as_ref().submatrix(0, 0, gemm_n, c));
                    for j in 0..c {
                        window.push_back((chunk_start + j, start_slot + j));
                    }
                } else {
                    for j in 0..c {
                        let slot = (ring_next + j) % ring_size;
                        ring_buf
                            .as_mut()
                            .submatrix_mut(0, slot, gemm_n, 1)
                            .copy_from(b_mat.as_ref().submatrix(0, j, gemm_n, 1));
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
                        .submatrix_mut(0, start_slot, gemm_n, c)
                        .copy_from(b_mat.as_ref().submatrix(0, 0, gemm_n, c));
                    for j in 0..c {
                        window.push_back((chunk_start + j, start_slot + j));
                    }
                } else {
                    for j in 0..c {
                        let slot = (ring_next + j) % ring_size;
                        ring_buf
                            .as_mut()
                            .submatrix_mut(0, slot, gemm_n, 1)
                            .copy_from(b_mat.as_ref().submatrix(0, j, gemm_n, 1));
                        window.push_back((chunk_start + j, slot));
                    }
                }
            }
        }
        ring_next += c;
        t_ring_store += t0.elapsed();

        // Fire the per-chunk progress callback. Native callers pass
        // a no-op; the WASM worker turns this into a postMessage
        // back to the main thread for the progress bar. Note: L2
        // contributions for SNPs in this chunk are NOT finalized yet
        // — later chunks add cross-window r² for them. This callback
        // is for *progress* (monotonic SNP-coverage counter), not
        // for streaming finalized output.
        //
        // Closure is `&mut dyn FnMut` (boxed once at the public API
        // boundary in `compute_l2_from_bed_with_progress`) so this
        // ~1900 LOC function only monomorphizes once across all
        // call sites.
        chunks_done += 1;
        on_progress(crate::l2::L2Progress {
            chunks_done,
            chunks_total,
            snps_done: chunk_end,
            snps_total: m,
            chunk_wall_ms: chunk_started.elapsed().as_secs_f64() * 1000.0,
        });
    } // end for chunk_start

    // Always populate the structured per-phase breakdown so callers
    // (browser orchestrator, CLI's verbose_timing path, future
    // benchmarks) can consume it programmatically. Cheap — just
    // fields-from-Durations.
    let perf = crate::l2::L2Perf {
        m,
        n_indiv,
        bed_read_secs: t_bed_read.as_secs_f64(),
        norm_secs: t_norm.as_secs_f64(),
        sketch_secs: sketch_dim.map(|_| t_sketch.as_secs_f64()),
        bb_dot_secs: t_bb_dot.as_secs_f64(),
        ab_dot_secs: t_ab_dot.as_secs_f64(),
        ring_store_secs: t_ring_store.as_secs_f64(),
    };

    if verbose_timing {
        // CLI-only eprintln. The wasm path collects via the L2Perf
        // return value and routes through serde-wasm-bindgen →
        // postMessage → main-thread tracing in worker_client.rs.
        // We intentionally do NOT call `tracing::info!` here: the
        // worker context lacks a tracing subscriber, and adding one
        // (`tracing_wasm::set_as_global_default`) interacts badly
        // with our multi-threaded WASM setup and traps at runtime.
        if let Some(sketch_secs) = perf.sketch_secs {
            eprintln!(
                "[perf] compute_ldscore: bed_read(stall)={:.3}s norm={:.3}s sketch={:.3}s bb_dot={:.3}s ab_dot={:.3}s ring_store={:.3}s",
                perf.bed_read_secs,
                perf.norm_secs,
                sketch_secs,
                perf.bb_dot_secs,
                perf.ab_dot_secs,
                perf.ring_store_secs,
            );
        } else {
            eprintln!(
                "[perf] compute_ldscore: bed_read(stall)={:.3}s norm={:.3}s bb_dot={:.3}s ab_dot={:.3}s ring_store={:.3}s",
                perf.bed_read_secs,
                perf.norm_secs,
                perf.bb_dot_secs,
                perf.ab_dot_secs,
                perf.ring_store_secs,
            );
        }
    }

    Ok((l2, maf_per_snp, perf))
}
