use super::normalize::{normalize_col_f32_with_stats, normalize_col_f64_with_stats, sum_sumsq_f32};
use super::window::{WindowMode, get_block_lefts_by_chr};
use crate::bed::{Bed, ChunkReader, IidPos, ReadOptions};
#[cfg(feature = "gpu")]
use crate::gpu::GpuContext;
use crate::la::{
    MatF, MatF32, mat_add_in_place, mat_copy_from, mat_slice, mat_slice_f32, mat_slice_mut,
    mat_slice_mut_f32, matmul_tn_to, matmul_tn_to_f32, matmul_to,
};
use crate::parse::BimRecord;
use anyhow::{Context, Result};
use faer::{Accum, Par};
use rayon::prelude::*;
use std::collections::VecDeque;

/// Fill a slice with Rademacher random values (+1.0 or -1.0) using bits from `rng`.
fn fill_rademacher_f64(rng: &mut fastrand::Rng, buf: &mut [f64]) {
    for chunk in buf.chunks_mut(64) {
        let bits = rng.u64(..);
        for (j, v) in chunk.iter_mut().enumerate() {
            *v = if (bits >> j) & 1 == 0 { 1.0 } else { -1.0 };
        }
    }
}

/// Fill a slice with Rademacher random values (+1.0 or -1.0) using bits from `rng`.
fn fill_rademacher_f32(rng: &mut fastrand::Rng, buf: &mut [f32]) {
    for chunk in buf.chunks_mut(64) {
        let bits = rng.u64(..);
        for (j, v) in chunk.iter_mut().enumerate() {
            *v = if (bits >> j) & 1 == 0 { 1.0 } else { -1.0 };
        }
    }
}

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
// Replaces the dense d×N Rademacher matrix with O(N) hash + sign arrays.
// Application: O(N×c) scatter-add vs O(d×N×c) for dense GEMM.

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
// BED bits 0b00=2, 0b01=missing, 0b10=1, 0b11=0 (count_a1=true convention).
const GENO_LUT_F64: [f64; 4] = [2.0, f64::NAN, 1.0, 0.0];
const GENO_LUT_F32: [f32; 4] = [2.0, f32::NAN, 1.0, 0.0];

// Byte-level stats LUT: each of the 256 possible BED byte values encodes 4
// genotypes; this LUT gives (sum, count, sum_sq) for the full byte.
#[derive(Clone, Copy)]
struct BedByteStats {
    sum: u8,
    count: u8,
    sum_sq: u8,
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

/// Fused BED-decode-normalize-project (f64).
/// Reads raw packed BED bytes, computes per-SNP statistics, then tiles over
/// individuals to decode+normalize+project without materializing the full N×c
/// intermediate matrix.
///
/// Pass 2 parallelizes over tiles via rayon: each thread decodes a tile of
/// individuals, computes P_tile × B_tile (a d×c partial result), and the
/// partial results are summed via `reduce`.
///
/// On return, `out[0..d, 0..c]` contains the projected+renormalized columns,
/// and `maf_out[0..c]` contains the per-SNP MAF.
#[allow(clippy::too_many_arguments)]
fn sketch_fused_project_f64(
    raw_bytes: &[u8],
    bytes_per_snp: usize,
    n_indiv: usize,
    c: usize,
    d: usize,
    proj: &MatF,
    n: f64,
    iid_positions: &[IidPos],
    all_iids: bool,
    maf_out: &mut [f64],
    out: &mut MatF,
) {
    use rayon::prelude::*;

    let byte_lut = build_bed_byte_lut();

    // --- Pass 1: per-SNP statistics from packed bytes ---
    struct SnpStats {
        mean: f64,
        inv_std: f64,
        maf: f64,
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
            SnpStats { mean, inv_std, maf }
        })
        .collect();

    for (j, s) in stats.iter().enumerate() {
        maf_out[j] = s.maf;
    }

    // --- Pass 2: parallel tiled decode + normalize + project ---
    const TILE: usize = 512;
    let tile_starts: Vec<usize> = (0..n_indiv).step_by(TILE).collect();

    // Each tile produces a d×c partial result; reduce by summing.
    let summed = tile_starts
        .par_iter()
        .map(|&tile_start| {
            let tile_end = (tile_start + TILE).min(n_indiv);
            let tile_n = tile_end - tile_start;

            let mut b_tile = MatF::zeros(tile_n, c);
            let mut partial = MatF::zeros(d, c);

            // Decode + normalize tile into b_tile[tile_n × c]
            for j in 0..c {
                let snp_bytes = &raw_bytes[j * bytes_per_snp..(j + 1) * bytes_per_snp];
                let mean = stats[j].mean;
                let inv_std = stats[j].inv_std;
                let col = b_tile
                    .col_mut(j)
                    .try_as_col_major_mut()
                    .unwrap()
                    .as_slice_mut();
                if all_iids {
                    for (i, v) in col[..tile_n].iter_mut().enumerate() {
                        let iid = tile_start + i;
                        let byte = snp_bytes[iid / 4];
                        let bits = (byte >> ((iid % 4) * 2)) & 0b11;
                        let g = GENO_LUT_F64[bits as usize];
                        *v = if g.is_nan() {
                            0.0
                        } else {
                            (g - mean) * inv_std
                        };
                    }
                } else {
                    for (v, pos) in col[..tile_n]
                        .iter_mut()
                        .zip(iid_positions[tile_start..tile_end].iter())
                    {
                        let byte = snp_bytes[pos.byte_idx];
                        let bits = (byte >> pos.shift) & 0b11;
                        let g = GENO_LUT_F64[bits as usize];
                        *v = if g.is_nan() {
                            0.0
                        } else {
                            (g - mean) * inv_std
                        };
                    }
                }
            }

            // Project: partial[d × c] = P[:, tile_start..tile_end] × b_tile[tile_n × c]
            matmul_to(
                mat_slice_mut(partial.as_mut(), 0..d, 0..c),
                mat_slice(proj.as_ref(), 0..d, tile_start..tile_end),
                mat_slice(b_tile.as_ref(), 0..tile_n, 0..c),
                1.0,
                Accum::Replace,
                Par::Seq,
            );

            partial
        })
        .reduce(
            || MatF::zeros(d, c),
            |mut acc, part| {
                mat_add_in_place(acc.as_mut(), part.as_ref());
                acc
            },
        );

    // Copy summed result into out
    mat_copy_from(
        mat_slice_mut(out.as_mut(), 0..d, 0..c),
        mat_slice(summed.as_ref(), 0..d, 0..c),
    );

    // Ratio estimator renorm: rescale each column so ||out[:,j]||² = n
    for j in 0..c {
        let col = out.col(j).try_as_col_major().unwrap().as_slice();
        let mut nrm_sq = 0.0f64;
        for &v in &col[..d] {
            nrm_sq += v * v;
        }
        if nrm_sq > 0.0 {
            let scale = (n / nrm_sq).sqrt();
            let col_mut = out
                .col_mut(j)
                .try_as_col_major_mut()
                .unwrap()
                .as_slice_mut();
            for v in &mut col_mut[..d] {
                *v *= scale;
            }
        }
    }
}

/// Fused BED-decode-normalize-project (f32 variant).
/// Same parallel tile strategy as the f64 variant — see `sketch_fused_project_f64`.
#[allow(clippy::too_many_arguments)]
fn sketch_fused_project_f32(
    raw_bytes: &[u8],
    bytes_per_snp: usize,
    n_indiv: usize,
    c: usize,
    d: usize,
    proj: &MatF32,
    n: f64,
    iid_positions: &[IidPos],
    all_iids: bool,
    maf_out: &mut [f64],
    out: &mut MatF32,
) {
    use rayon::prelude::*;

    let byte_lut = build_bed_byte_lut();

    // --- Pass 1: per-SNP statistics from packed bytes ---
    struct SnpStats {
        mean: f32,
        inv_std: f32,
        maf: f64,
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
            SnpStats {
                mean: mean_f64 as f32,
                inv_std,
                maf,
            }
        })
        .collect();

    for (j, s) in stats.iter().enumerate() {
        maf_out[j] = s.maf;
    }

    // --- Pass 2: parallel tiled decode + normalize + project ---
    const TILE: usize = 512;
    let tile_starts: Vec<usize> = (0..n_indiv).step_by(TILE).collect();

    let summed = tile_starts
        .par_iter()
        .map(|&tile_start| {
            let tile_end = (tile_start + TILE).min(n_indiv);
            let tile_n = tile_end - tile_start;

            let mut b_tile = MatF32::zeros(tile_n, c);
            let mut partial = MatF32::zeros(d, c);

            for j in 0..c {
                let snp_bytes = &raw_bytes[j * bytes_per_snp..(j + 1) * bytes_per_snp];
                let mean = stats[j].mean;
                let inv_std = stats[j].inv_std;
                let col = b_tile
                    .col_mut(j)
                    .try_as_col_major_mut()
                    .unwrap()
                    .as_slice_mut();
                if all_iids {
                    for (i, v) in col[..tile_n].iter_mut().enumerate() {
                        let iid = tile_start + i;
                        let byte = snp_bytes[iid / 4];
                        let bits = (byte >> ((iid % 4) * 2)) & 0b11;
                        let g = GENO_LUT_F32[bits as usize];
                        *v = if g.is_nan() {
                            0.0
                        } else {
                            (g - mean) * inv_std
                        };
                    }
                } else {
                    for (v, pos) in col[..tile_n]
                        .iter_mut()
                        .zip(iid_positions[tile_start..tile_end].iter())
                    {
                        let byte = snp_bytes[pos.byte_idx];
                        let bits = (byte >> pos.shift) & 0b11;
                        let g = GENO_LUT_F32[bits as usize];
                        *v = if g.is_nan() {
                            0.0
                        } else {
                            (g - mean) * inv_std
                        };
                    }
                }
            }

            faer::linalg::matmul::matmul(
                mat_slice_mut_f32(partial.as_mut(), 0..d, 0..c),
                Accum::Replace,
                mat_slice_f32(proj.as_ref(), 0..d, tile_start..tile_end),
                mat_slice_f32(b_tile.as_ref(), 0..tile_n, 0..c),
                1.0f32,
                Par::Seq,
            );

            partial
        })
        .reduce(
            || MatF32::zeros(d, c),
            |mut acc, part| {
                for j in 0..c {
                    let acc_col = acc
                        .col_mut(j)
                        .try_as_col_major_mut()
                        .unwrap()
                        .as_slice_mut();
                    let part_col = part.col(j).try_as_col_major().unwrap().as_slice();
                    for (a, &p) in acc_col[..d].iter_mut().zip(part_col[..d].iter()) {
                        *a += p;
                    }
                }
                acc
            },
        );

    // Copy summed result into out
    for j in 0..c {
        let out_col = out
            .col_mut(j)
            .try_as_col_major_mut()
            .unwrap()
            .as_slice_mut();
        let sum_col = summed.col(j).try_as_col_major().unwrap().as_slice();
        out_col[..d].copy_from_slice(&sum_col[..d]);
    }

    // Ratio estimator renorm
    for j in 0..c {
        let col = out.col(j).try_as_col_major().unwrap().as_slice();
        let mut nrm_sq = 0.0f64;
        for &v in &col[..d] {
            nrm_sq += (v as f64) * (v as f64);
        }
        if nrm_sq > 0.0 {
            let scale = (n / nrm_sq).sqrt() as f32;
            let col_mut = out
                .col_mut(j)
                .try_as_col_major_mut()
                .unwrap()
                .as_slice_mut();
            for v in &mut col_mut[..d] {
                *v *= scale;
            }
        }
    }
}

// ── Fused CountSketch: BED-decode-normalize-scatter-add ─────────────────────
// Avoids materializing the full N×c intermediate matrix AND avoids any GEMM.
// Pass 1: compute per-SNP stats (mean, inv_std, maf) from packed BED bytes.
// Pass 2: parallel over columns — decode each genotype, normalize, scatter-add
//         to the appropriate bucket. O(N×c) total, column-parallel via rayon.
//
// Unlike fused Rademacher (which fragments one large GEMM into many tiny tile
// GEMMs), fused CountSketch has NO GEMM at all — just scatter-adds. The column-
// parallel pattern gives c-way (c=200) parallelism with no SIMD efficiency loss.

/// Fused BED-decode-normalize-CountSketch (f32 variant, primary path since sketch auto-enables f32).
/// Reads raw packed BED bytes, computes per-SNP statistics, then scatter-adds
/// normalized genotypes into d-bucket sketch vectors — all without materializing
/// the full N×c genotype matrix.
///
/// Parallelized over columns (j=0..c) via rayon. Each column's scatter-add is
/// independent, giving c-way parallelism (c=chunk_size, typically 200).
///
/// On return, `out[0..d, 0..c]` contains the projected+renormalized columns,
/// and `maf_out[0..c]` contains the per-SNP MAF.
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
    out: &mut MatF32,
) {
    let d = cs.d;
    let byte_lut = build_bed_byte_lut();

    // --- Pass 1: per-SNP statistics from packed bytes ---
    struct SnpStats {
        mean: f32,
        inv_std: f32,
        maf: f64,
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
            SnpStats {
                mean: mean_f64 as f32,
                inv_std,
                maf,
            }
        })
        .collect();

    for (j, s) in stats.iter().enumerate() {
        maf_out[j] = s.maf;
    }

    // --- Pass 2: column-parallel fused decode + normalize + scatter-add ---
    // Each column j is independent: read BED bytes for SNP j, decode each
    // individual's genotype, normalize, and scatter-add to out[bucket[i], j].
    let out_ptr = out.as_mut().as_ptr_mut();
    let out_col_stride = out.col_stride();
    let out_send = SendPtr(out_ptr);
    let bucket = &cs.bucket;
    let sign = &cs.sign_f32;

    (0..c).into_par_iter().for_each(|j| {
        let snp_bytes = &raw_bytes[j * bytes_per_snp..(j + 1) * bytes_per_snp];
        let mean = stats[j].mean;
        let inv_std = stats[j].inv_std;

        // Safety: each thread writes to a disjoint column of out.
        unsafe {
            let dst = std::slice::from_raw_parts_mut(
                out_send.ptr().offset(j as isize * out_col_stride),
                d,
            );

            // Zero output column.
            for v in dst.iter_mut() {
                *v = 0.0;
            }

            // Fused decode + normalize + scatter-add
            if all_iids {
                for i in 0..n_indiv {
                    let byte = snp_bytes[i / 4];
                    let bits = (byte >> ((i % 4) * 2)) & 0b11;
                    let g = GENO_LUT_F32[bits as usize];
                    if !g.is_nan() {
                        let val = (g - mean) * inv_std;
                        *dst.get_unchecked_mut(*bucket.get_unchecked(i) as usize) +=
                            *sign.get_unchecked(i) * val;
                    }
                }
            } else {
                for (i, pos) in iid_positions.iter().enumerate() {
                    let byte = snp_bytes[pos.byte_idx];
                    let bits = (byte >> pos.shift) & 0b11;
                    let g = GENO_LUT_F32[bits as usize];
                    if !g.is_nan() {
                        let val = (g - mean) * inv_std;
                        *dst.get_unchecked_mut(*bucket.get_unchecked(i) as usize) +=
                            *sign.get_unchecked(i) * val;
                    }
                }
            }

            // Ratio estimator renorm: rescale so ||col||^2 = n (accumulate in f64)
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
    out: &mut MatF,
) {
    let d = cs.d;
    let byte_lut = build_bed_byte_lut();

    // --- Pass 1: per-SNP statistics from packed bytes ---
    struct SnpStats {
        mean: f64,
        inv_std: f64,
        maf: f64,
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
            SnpStats { mean, inv_std, maf }
        })
        .collect();

    for (j, s) in stats.iter().enumerate() {
        maf_out[j] = s.maf;
    }

    // --- Pass 2: column-parallel fused decode + normalize + scatter-add ---
    let out_ptr = out.as_mut().as_ptr_mut();
    let out_col_stride = out.col_stride();
    let out_send = SendPtr(out_ptr);
    let bucket = &cs.bucket;
    let sign = &cs.sign_f64;

    (0..c).into_par_iter().for_each(|j| {
        let snp_bytes = &raw_bytes[j * bytes_per_snp..(j + 1) * bytes_per_snp];
        let mean = stats[j].mean;
        let inv_std = stats[j].inv_std;

        unsafe {
            let dst = std::slice::from_raw_parts_mut(
                out_send.ptr().offset(j as isize * out_col_stride),
                d,
            );

            for v in dst.iter_mut() {
                *v = 0.0;
            }

            if all_iids {
                for i in 0..n_indiv {
                    let byte = snp_bytes[i / 4];
                    let bits = (byte >> ((i % 4) * 2)) & 0b11;
                    let g = GENO_LUT_F64[bits as usize];
                    if !g.is_nan() {
                        let val = (g - mean) * inv_std;
                        *dst.get_unchecked_mut(*bucket.get_unchecked(i) as usize) +=
                            *sign.get_unchecked(i) * val;
                    }
                }
            } else {
                for (i, pos) in iid_positions.iter().enumerate() {
                    let byte = snp_bytes[pos.byte_idx];
                    let bits = (byte >> pos.shift) & 0b11;
                    let g = GENO_LUT_F64[bits as usize];
                    if !g.is_nan() {
                        let val = (g - mean) * inv_std;
                        *dst.get_unchecked_mut(*bucket.get_unchecked(i) as usize) +=
                            *sign.get_unchecked(i) * val;
                    }
                }
            }

            // Ratio estimator renorm
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

/// Unbiased r² estimator: r² − (1−r²)/(n−2) [Bulik-Sullivan 2015].
/// Inlined as `r*r * r2u_a + r2u_b` in hot path for constant n.
#[inline]
#[allow(dead_code)]
pub(super) fn r2_unbiased(r: f64, n: usize) -> f64 {
    let sq = r * r;
    let denom = if n > 2 { n as f64 - 2.0 } else { n as f64 };
    sq - (1.0 - sq) / denom
}

/// Compute LD scores for all SNPs. Returns `(l2, maf)`.
#[allow(clippy::too_many_arguments, clippy::unnecessary_cast)]
pub(super) fn compute_ldscore_global(
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
    stochastic: Option<usize>,
    sketch: Option<usize>,
    sketch_method: &str,
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

    // ── Randomized sketching setup ──────────────────────────────────────────
    // Compress individual dimension N→d via random projection P (d×N).
    // All subsequent GEMMs operate on d-dimensional data.
    // Bias correction: E[val̃²] = val²(1+1/d) + N²/d, absorbed into adjusted constants.
    let sketch_dim: Option<usize> = sketch.and_then(|d| {
        if d < 3 || d >= n_indiv {
            eprintln!("WARNING: --sketch {d} ignored (must be 3 ≤ d < n_indiv={n_indiv})");
            return None;
        }
        Some(d)
    });
    let gemm_n: usize; // row dimension for all GEMM matrices (d when sketching, N otherwise)
    let active_n_inv_sq_r2u_a: f64;
    let active_r2u_b: f64;
    // Projection matrices (one per precision path, allocated only when sketching)
    let proj_f64: Option<MatF>;
    let proj_f32: Option<MatF32>;

    if let Some(d) = sketch_dim {
        gemm_n = d;
        // Bias correction with column re-normalization (ratio estimator).
        // After projecting, each column is rescaled so ||x̃'_j||² = N.
        // This reduces per-pair variance by ~2× compared to the raw estimator.
        // The squared cosine r̃² = (val̃'/N)² has bias ≈ (2r⁴ - 7r² + 1)/d;
        // the leading linear correction is: r²_corr = r̃² × d/(d-2) - 1/(d-2).
        // Residual r⁴/(d-2) per pair is negligible for typical LD.
        let d_f = d as f64;
        let dm2 = d_f - 2.0;
        active_n_inv_sq_r2u_a = n_inv_sq_r2u_a * d_f / dm2;
        active_r2u_b = r2u_b - r2u_a / dm2;

        // Only allocate the dense Rademacher P matrix when actually needed.
        // CountSketch uses hash arrays instead; fused path uses P at projection time.
        if sketch_method == "countsketch" {
            proj_f64 = None;
            proj_f32 = None;
        } else {
            let inv_sqrt_d = 1.0 / d_f.sqrt();
            let mut rng_proj = fastrand::Rng::with_seed(42);
            if use_f32 {
                let mut p = MatF32::zeros(d, n_indiv);
                for col in 0..n_indiv {
                    let sl = p
                        .col_mut(col)
                        .try_as_col_major_mut()
                        .unwrap()
                        .as_slice_mut();
                    for chunk in sl[..d].chunks_mut(64) {
                        let bits = rng_proj.u64(..);
                        for (j, v) in chunk.iter_mut().enumerate() {
                            *v = if (bits >> j) & 1 == 0 {
                                inv_sqrt_d as f32
                            } else {
                                -(inv_sqrt_d as f32)
                            };
                        }
                    }
                }
                proj_f32 = Some(p);
                proj_f64 = None;
            } else {
                let mut p = MatF::zeros(d, n_indiv);
                for col in 0..n_indiv {
                    let sl = p
                        .col_mut(col)
                        .try_as_col_major_mut()
                        .unwrap()
                        .as_slice_mut();
                    for chunk in sl[..d].chunks_mut(64) {
                        let bits = rng_proj.u64(..);
                        for (j, v) in chunk.iter_mut().enumerate() {
                            *v = if (bits >> j) & 1 == 0 {
                                inv_sqrt_d
                            } else {
                                -inv_sqrt_d
                            };
                        }
                    }
                }
                proj_f64 = Some(p);
                proj_f32 = None;
            }
        }
        println!(
            "Sketch mode: projecting N={} → d={} (bias-corrected, method={})",
            n_indiv, d, sketch_method
        );
    } else {
        gemm_n = n_indiv;
        active_n_inv_sq_r2u_a = n_inv_sq_r2u_a;
        active_r2u_b = r2u_b;
        proj_f64 = None;
        proj_f32 = None;
    }

    // CountSketch projection (allocated only when sketch_method == "countsketch").
    // Uses seed 42, same as Rademacher (only one is active at a time).
    let count_sketch: Option<CountSketchProj> = sketch_dim.and_then(|d| {
        if sketch_method == "countsketch" {
            Some(CountSketchProj::new(n_indiv, d, 42))
        } else {
            None
        }
    });

    // Fused BED-decode-normalize-project for Rademacher: DISABLED.
    // The tiled approach (512-individual tiles with Par::Seq GEMMs + rayon across tiles)
    // saves ~40MB (N×c buffer) but is 2.2× slower than the non-fused path at biobank scale
    // (285s vs 130s). The non-fused path uses a single large P×B GEMM with Par::rayon(0),
    // which is far more efficient due to faer's SIMD panel scheduling.
    //
    // Fused CountSketch: ENABLED. CountSketch is a scatter-add (not a GEMM), so the
    // tile-fragmentation problem doesn't apply. Column-parallel (c-way) scatter-add
    // gives good parallelism without any GEMM overhead.
    let use_fused_sketch = sketch_dim.is_some() && sketch_method == "countsketch";

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

    // Determine if stochastic Hutchinson mode is usable.
    let stochastic_probes: Option<usize> = stochastic.and_then(|t| {
        if t == 0 {
            return None;
        }
        if n_annot > 1 {
            eprintln!(
                "WARNING: --stochastic not supported for partitioned LD scores; using exact GEMM"
            );
            return None;
        }
        if pq_exp.is_some() {
            eprintln!("WARNING: --stochastic not supported with --pq-exp; using exact GEMM");
            return None;
        }
        Some(t)
    });

    // Stochastic Hutchinson scratch buffers — batched probe matrices.
    // Instead of T separate mat-vecs, batch all probes into a single GEMM:
    //   Z (c×T random), Y = B*Z (n×T), Q = B'*Y (c×T), hutch[j] = sum_t Q[j,t]^2
    let mut rng = stochastic_probes.map(|_| fastrand::Rng::with_seed(42));
    let n_probes_alloc = stochastic_probes.unwrap_or(0);
    let mut hutch_z_c_f64 = if n_probes_alloc > 0 && !use_f32 {
        MatF::zeros(chunk_c, n_probes_alloc)
    } else {
        MatF::zeros(0, 0)
    };
    let mut hutch_z_c_f32 = if n_probes_alloc > 0 && use_f32 {
        MatF32::zeros(chunk_c, n_probes_alloc)
    } else {
        MatF32::zeros(0, 0)
    };
    let mut hutch_z_w_f64 = if n_probes_alloc > 0 && !use_f32 {
        MatF::zeros(max_window_size.max(1), n_probes_alloc)
    } else {
        MatF::zeros(0, 0)
    };
    let mut hutch_z_w_f32 = if n_probes_alloc > 0 && use_f32 {
        MatF32::zeros(max_window_size.max(1), n_probes_alloc)
    } else {
        MatF32::zeros(0, 0)
    };
    let mut hutch_y_f64 = if n_probes_alloc > 0 && !use_f32 {
        MatF::zeros(gemm_n, n_probes_alloc)
    } else {
        MatF::zeros(0, 0)
    };
    let mut hutch_y_f32 = if n_probes_alloc > 0 && use_f32 {
        MatF32::zeros(gemm_n, n_probes_alloc)
    } else {
        MatF32::zeros(0, 0)
    };
    let mut hutch_q_f64 = if n_probes_alloc > 0 && !use_f32 {
        MatF::zeros(chunk_c, n_probes_alloc)
    } else {
        MatF::zeros(0, 0)
    };
    let mut hutch_q_f32 = if n_probes_alloc > 0 && use_f32 {
        MatF32::zeros(chunk_c, n_probes_alloc)
    } else {
        MatF32::zeros(0, 0)
    };
    let mut hutch_p_f64 = if n_probes_alloc > 0 && !use_f32 {
        MatF::zeros(max_window_size.max(1), n_probes_alloc)
    } else {
        MatF::zeros(0, 0)
    };
    let mut hutch_p_f32 = if n_probes_alloc > 0 && use_f32 {
        MatF32::zeros(max_window_size.max(1), n_probes_alloc)
    } else {
        MatF32::zeros(0, 0)
    };

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

    use std::time::Instant;
    let mut t_bed_read = std::time::Duration::ZERO;
    let mut t_norm = std::time::Duration::ZERO;
    let mut t_sketch = std::time::Duration::ZERO;
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
                ChunkSpec {
                    chunk_start: start,
                    chunk_end: end,
                    bed_indices,
                }
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
                        format!(
                            "reading BED chunk [{},{})",
                            spec.chunk_start, spec.chunk_end
                        )
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
        // When using the fused sketch path (sequential reader + CountSketch), read raw bytes
        // and skip b_full materialization. Otherwise use the existing decoded-Mat path.
        let sketching = sketch_dim.is_some();
        // Determine if this chunk qualifies for the fused path.
        // Fused requires: (a) sequential mode OR contiguous chunk, (b) not prefetch.
        let chunk_is_contiguous = if !prefetch_bed && !all_bed_sequential {
            let start_bed_idx = all_snps[chunk_start].bed_idx;
            all_snps[chunk_end - 1].bed_idx == start_bed_idx + c - 1
        } else {
            all_bed_sequential
        };
        let do_fused = use_fused_sketch && chunk_is_contiguous;

        if do_fused {
            // ── Fused path: read raw bytes, fuse decode+normalize+project ──────────
            let cr = chunk_reader.as_mut().unwrap();
            let bed = seq_bed.as_mut().unwrap();
            if all_bed_sequential {
                cr.read_next_raw_only(bed, c).with_context(|| {
                    format!("reading BED chunk [{},{})", chunk_start, chunk_end)
                })?;
            } else {
                let start_bed_idx = all_snps[chunk_start].bed_idx;
                cr.read_contiguous_raw_only(bed, start_bed_idx, c)
                    .with_context(|| {
                        format!("reading BED chunk [{},{})", chunk_start, chunk_end)
                    })?;
            }
            t_bed_read += t0.elapsed();

            let t0 = Instant::now();
            let bps = cr.bytes_per_snp_val();
            let raw = cr.raw_bytes()[..c * bps].to_vec(); // copy to avoid borrow conflict
            let n_i = cr.n_iid();
            let ipos = cr.iid_positions_ref().to_vec();
            let ai = cr.all_iids_flag();
            if let Some(ref cs) = count_sketch {
                // Fused CountSketch: decode + normalize + scatter-add (no GEMM)
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
                            b_mat,
                        );
                    }
                }
            } else {
                // Fused Rademacher (currently disabled via use_fused_sketch)
                let d = sketch_dim.unwrap();
                match bufs {
                    GemmBufs::F32 { ref mut b_mat, .. } => {
                        sketch_fused_project_f32(
                            &raw,
                            bps,
                            n_i,
                            c,
                            d,
                            proj_f32.as_ref().unwrap(),
                            n,
                            &ipos,
                            ai,
                            &mut maf_per_snp[chunk_start..chunk_start + c],
                            b_mat,
                        );
                    }
                    GemmBufs::F64 { ref mut b_mat, .. } => {
                        sketch_fused_project_f64(
                            &raw,
                            bps,
                            n_i,
                            c,
                            d,
                            proj_f64.as_ref().unwrap(),
                            n,
                            &ipos,
                            ai,
                            &mut maf_per_snp[chunk_start..chunk_start + c],
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

            // Sketch projection: b_mat (d × c) = proj(b_full)
            // Dispatches to CountSketch (O(N) scatter-add) or Rademacher dense GEMM.
            if sketching {
                let t0 = Instant::now();
                let d = sketch_dim.unwrap();
                if let Some(ref cs) = count_sketch {
                    // CountSketch: O(N×c) scatter-add, includes ratio renorm.
                    match bufs {
                        GemmBufs::F32 { ref mut b_mat, .. } => {
                            cs.project_f32(&b_full_f32, n_indiv, c, n, b_mat);
                        }
                        GemmBufs::F64 { ref mut b_mat, .. } => {
                            cs.project_f64(&b_full_f64, n_indiv, c, n, b_mat);
                        }
                    }
                } else {
                    // Rademacher: dense GEMM + ratio renorm.
                    match bufs {
                        GemmBufs::F32 { ref mut b_mat, .. } => {
                            let p = proj_f32.as_ref().unwrap();
                            faer::linalg::matmul::matmul(
                                mat_slice_mut_f32(b_mat.as_mut(), 0..d, 0..c),
                                Accum::Replace,
                                mat_slice_f32(p.as_ref(), 0..d, 0..n_indiv),
                                mat_slice_f32(b_full_f32.as_ref(), 0..n_indiv, 0..c),
                                1.0f32,
                                Par::rayon(0),
                            );
                            // Re-normalize: ||x̃'_j||² = N (ratio estimator).
                            for j in 0..c {
                                let col = b_mat.col(j).try_as_col_major().unwrap().as_slice();
                                let mut nrm_sq = 0.0f64;
                                for &v in &col[..d] {
                                    nrm_sq += (v as f64) * (v as f64);
                                }
                                if nrm_sq > 0.0 {
                                    let scale = ((n / nrm_sq).sqrt()) as f32;
                                    let col_mut = b_mat
                                        .col_mut(j)
                                        .try_as_col_major_mut()
                                        .unwrap()
                                        .as_slice_mut();
                                    for v in &mut col_mut[..d] {
                                        *v *= scale;
                                    }
                                }
                            }
                        }
                        GemmBufs::F64 { ref mut b_mat, .. } => {
                            let p = proj_f64.as_ref().unwrap();
                            matmul_to(
                                mat_slice_mut(b_mat.as_mut(), 0..d, 0..c),
                                mat_slice(p.as_ref(), 0..d, 0..n_indiv),
                                mat_slice(b_full_f64.as_ref(), 0..n_indiv, 0..c),
                                1.0,
                                Accum::Replace,
                                Par::rayon(0),
                            );
                            // Re-normalize: ||x̃'_j||² = N (ratio estimator)
                            for j in 0..c {
                                let col = b_mat.col(j).try_as_col_major().unwrap().as_slice();
                                let mut nrm_sq = 0.0f64;
                                for &v in &col[..d] {
                                    nrm_sq += v * v;
                                }
                                if nrm_sq > 0.0 {
                                    let scale = (n / nrm_sq).sqrt();
                                    let col_mut = b_mat
                                        .col_mut(j)
                                        .try_as_col_major_mut()
                                        .unwrap()
                                        .as_slice_mut();
                                    for v in &mut col_mut[..d] {
                                        *v *= scale;
                                    }
                                }
                            }
                        }
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

        // ---- Stochastic Hutchinson path (batched GEMM, scalar LD only) ----
        // Instead of T separate mat-vecs, batch all probes into matrix operations:
        //   Z (c×T random), Y = B*Z (n×T), Q = B'*Y = (B'B)*Z (c×T)
        //   hutch[j] = sum_t Q[j,t]^2 / T  ≈  sum_k (B'B)[j,k]^2
        if let Some(n_probes) = stochastic_probes {
            let rng = rng.as_mut().unwrap();
            let n_sq = n * n;
            let inv_t = 1.0 / n_probes as f64;

            /// Sum of squared elements per row of a matrix: result[i] = sum_j mat[i,j]^2
            fn row_sum_sq(mat: faer::MatRef<'_, f64>, nrows: usize, ncols: usize) -> Vec<f64> {
                let mut out = vec![0.0f64; nrows];
                for j in 0..ncols {
                    let col = mat.col(j).try_as_col_major().unwrap().as_slice();
                    for (i, v) in col[..nrows].iter().enumerate() {
                        out[i] += v * v;
                    }
                }
                out
            }

            /// Sum of squared elements per row (f32 input, f64 accumulation).
            fn row_sum_sq_f32(mat: faer::MatRef<'_, f32>, nrows: usize, ncols: usize) -> Vec<f64> {
                let mut out = vec![0.0f64; nrows];
                for j in 0..ncols {
                    let col = mat.col(j).try_as_col_major().unwrap().as_slice();
                    for (i, v) in col[..nrows].iter().enumerate() {
                        let vf = *v as f64;
                        out[i] += vf * vf;
                    }
                }
                out
            }

            // bb via batched Hutchinson
            let t0 = Instant::now();
            let hutch_bb = match bufs {
                GemmBufs::F64 { ref b_mat, .. } => {
                    let b_sl = mat_slice(b_mat.as_ref(), 0..gemm_n, 0..c);
                    // Fill Z (c × T) with Rademacher values
                    for t in 0..n_probes {
                        fill_rademacher_f64(
                            rng,
                            &mut hutch_z_c_f64
                                .col_mut(t)
                                .try_as_col_major_mut()
                                .unwrap()
                                .as_slice_mut()[..c],
                        );
                    }
                    let z_ref = mat_slice(hutch_z_c_f64.as_ref(), 0..c, 0..n_probes);
                    // Y = B * Z  (n × T)
                    matmul_to(
                        mat_slice_mut(hutch_y_f64.as_mut(), 0..gemm_n, 0..n_probes),
                        b_sl,
                        z_ref,
                        1.0,
                        Accum::Replace,
                        Par::rayon(0),
                    );
                    let y_ref = mat_slice(hutch_y_f64.as_ref(), 0..gemm_n, 0..n_probes);
                    // Q = B' * Y = (B'B) * Z  (c × T)
                    matmul_tn_to(
                        mat_slice_mut(hutch_q_f64.as_mut(), 0..c, 0..n_probes),
                        b_sl,
                        y_ref,
                        1.0,
                        Accum::Replace,
                        Par::rayon(0),
                    );
                    row_sum_sq(
                        mat_slice(hutch_q_f64.as_ref(), 0..c, 0..n_probes),
                        c,
                        n_probes,
                    )
                }
                GemmBufs::F32 { ref b_mat, .. } => {
                    let b_sl = mat_slice_f32(b_mat.as_ref(), 0..gemm_n, 0..c);
                    for t in 0..n_probes {
                        fill_rademacher_f32(
                            rng,
                            &mut hutch_z_c_f32
                                .col_mut(t)
                                .try_as_col_major_mut()
                                .unwrap()
                                .as_slice_mut()[..c],
                        );
                    }
                    let z_ref = mat_slice_f32(hutch_z_c_f32.as_ref(), 0..c, 0..n_probes);
                    faer::linalg::matmul::matmul(
                        mat_slice_mut_f32(hutch_y_f32.as_mut(), 0..gemm_n, 0..n_probes),
                        Accum::Replace,
                        b_sl,
                        z_ref,
                        1.0f32,
                        Par::rayon(0),
                    );
                    let y_ref = mat_slice_f32(hutch_y_f32.as_ref(), 0..gemm_n, 0..n_probes);
                    faer::linalg::matmul::matmul(
                        mat_slice_mut_f32(hutch_q_f32.as_mut(), 0..c, 0..n_probes),
                        Accum::Replace,
                        b_sl.transpose(),
                        y_ref,
                        1.0f32,
                        Par::rayon(0),
                    );
                    row_sum_sq_f32(
                        mat_slice_f32(hutch_q_f32.as_ref(), 0..c, 0..n_probes),
                        c,
                        n_probes,
                    )
                }
            };
            // Convert: hutch_bb[j]/T ≈ sum_k G[j,k]^2 (including diagonal n^2)
            for j in 0..c {
                let est = hutch_bb[j] * inv_t;
                l2[(chunk_start + j, 0)] +=
                    1.0 + active_n_inv_sq_r2u_a * (est - n_sq) + (c as f64 - 1.0) * active_r2u_b;
            }
            t_bb_dot += t0.elapsed();

            // ab via batched Hutchinson
            let t0 = Instant::now();
            if !window.is_empty() {
                let w = window.len();
                let (first_slot, last_slot) = window
                    .front()
                    .and_then(|(_, f)| window.back().map(|(_, l)| (*f, *l)))
                    .unwrap_or((0, 0));
                let contiguous = first_slot <= last_slot && last_slot - first_slot + 1 == w;

                match bufs {
                    GemmBufs::F64 {
                        ref ring_buf,
                        ref b_mat,
                        ref mut a_buf,
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

                        // ab left: Z_c (c×T), Y=B*Z_c, P=A'*Y → hutch_ab_left
                        for t in 0..n_probes {
                            fill_rademacher_f64(
                                rng,
                                &mut hutch_z_c_f64
                                    .col_mut(t)
                                    .try_as_col_major_mut()
                                    .unwrap()
                                    .as_slice_mut()[..c],
                            );
                        }
                        let z_ref = mat_slice(hutch_z_c_f64.as_ref(), 0..c, 0..n_probes);
                        matmul_to(
                            mat_slice_mut(hutch_y_f64.as_mut(), 0..gemm_n, 0..n_probes),
                            b_sl,
                            z_ref,
                            1.0,
                            Accum::Replace,
                            Par::rayon(0),
                        );
                        let y_ref = mat_slice(hutch_y_f64.as_ref(), 0..gemm_n, 0..n_probes);
                        matmul_tn_to(
                            mat_slice_mut(hutch_p_f64.as_mut(), 0..w, 0..n_probes),
                            a_view,
                            y_ref,
                            1.0,
                            Accum::Replace,
                            Par::rayon(0),
                        );
                        let hutch_ab_left = row_sum_sq(
                            mat_slice(hutch_p_f64.as_ref(), 0..w, 0..n_probes),
                            w,
                            n_probes,
                        );
                        for (wi, (snp_idx, _)) in window.iter().enumerate() {
                            let est = hutch_ab_left[wi] * inv_t;
                            l2[(*snp_idx, 0)] +=
                                active_n_inv_sq_r2u_a * est + c as f64 * active_r2u_b;
                        }

                        // ab right: Z_w (w×T), Y=A*Z_w, Q=B'*Y → hutch_ab_right
                        for t in 0..n_probes {
                            fill_rademacher_f64(
                                rng,
                                &mut hutch_z_w_f64
                                    .col_mut(t)
                                    .try_as_col_major_mut()
                                    .unwrap()
                                    .as_slice_mut()[..w],
                            );
                        }
                        let z_ref = mat_slice(hutch_z_w_f64.as_ref(), 0..w, 0..n_probes);
                        matmul_to(
                            mat_slice_mut(hutch_y_f64.as_mut(), 0..gemm_n, 0..n_probes),
                            a_view,
                            z_ref,
                            1.0,
                            Accum::Replace,
                            Par::rayon(0),
                        );
                        let y_ref = mat_slice(hutch_y_f64.as_ref(), 0..gemm_n, 0..n_probes);
                        matmul_tn_to(
                            mat_slice_mut(hutch_q_f64.as_mut(), 0..c, 0..n_probes),
                            b_sl,
                            y_ref,
                            1.0,
                            Accum::Replace,
                            Par::rayon(0),
                        );
                        let hutch_ab_right = row_sum_sq(
                            mat_slice(hutch_q_f64.as_ref(), 0..c, 0..n_probes),
                            c,
                            n_probes,
                        );
                        for j in 0..c {
                            let est = hutch_ab_right[j] * inv_t;
                            l2[(chunk_start + j, 0)] +=
                                active_n_inv_sq_r2u_a * est + w as f64 * active_r2u_b;
                        }
                    }
                    GemmBufs::F32 {
                        ref ring_buf,
                        ref b_mat,
                        ref mut a_buf,
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

                        // ab left
                        for t in 0..n_probes {
                            fill_rademacher_f32(
                                rng,
                                &mut hutch_z_c_f32
                                    .col_mut(t)
                                    .try_as_col_major_mut()
                                    .unwrap()
                                    .as_slice_mut()[..c],
                            );
                        }
                        let z_ref = mat_slice_f32(hutch_z_c_f32.as_ref(), 0..c, 0..n_probes);
                        faer::linalg::matmul::matmul(
                            mat_slice_mut_f32(hutch_y_f32.as_mut(), 0..gemm_n, 0..n_probes),
                            Accum::Replace,
                            b_sl,
                            z_ref,
                            1.0f32,
                            Par::rayon(0),
                        );
                        let y_ref = mat_slice_f32(hutch_y_f32.as_ref(), 0..gemm_n, 0..n_probes);
                        faer::linalg::matmul::matmul(
                            mat_slice_mut_f32(hutch_p_f32.as_mut(), 0..w, 0..n_probes),
                            Accum::Replace,
                            a_view.transpose(),
                            y_ref,
                            1.0f32,
                            Par::rayon(0),
                        );
                        let hutch_ab_left = row_sum_sq_f32(
                            mat_slice_f32(hutch_p_f32.as_ref(), 0..w, 0..n_probes),
                            w,
                            n_probes,
                        );
                        for (wi, (snp_idx, _)) in window.iter().enumerate() {
                            let est = hutch_ab_left[wi] * inv_t;
                            l2[(*snp_idx, 0)] +=
                                active_n_inv_sq_r2u_a * est + c as f64 * active_r2u_b;
                        }

                        // ab right
                        for t in 0..n_probes {
                            fill_rademacher_f32(
                                rng,
                                &mut hutch_z_w_f32
                                    .col_mut(t)
                                    .try_as_col_major_mut()
                                    .unwrap()
                                    .as_slice_mut()[..w],
                            );
                        }
                        let z_ref = mat_slice_f32(hutch_z_w_f32.as_ref(), 0..w, 0..n_probes);
                        faer::linalg::matmul::matmul(
                            mat_slice_mut_f32(hutch_y_f32.as_mut(), 0..gemm_n, 0..n_probes),
                            Accum::Replace,
                            a_view,
                            z_ref,
                            1.0f32,
                            Par::rayon(0),
                        );
                        let y_ref = mat_slice_f32(hutch_y_f32.as_ref(), 0..gemm_n, 0..n_probes);
                        faer::linalg::matmul::matmul(
                            mat_slice_mut_f32(hutch_q_f32.as_mut(), 0..c, 0..n_probes),
                            Accum::Replace,
                            b_sl.transpose(),
                            y_ref,
                            1.0f32,
                            Par::rayon(0),
                        );
                        let hutch_ab_right = row_sum_sq_f32(
                            mat_slice_f32(hutch_q_f32.as_ref(), 0..c, 0..n_probes),
                            c,
                            n_probes,
                        );
                        for j in 0..c {
                            let est = hutch_ab_right[j] * inv_t;
                            l2[(chunk_start + j, 0)] +=
                                active_n_inv_sq_r2u_a * est + w as f64 * active_r2u_b;
                        }
                    }
                }
            }
            t_ab_dot += t0.elapsed();
        } else {
            // ---- Exact GEMM path ----

            // B×B — compute GEMM then extract r2_unbiased into r2u_bb (always f64).
            let t0 = Instant::now();
            if !chunk_all_zero {
                // Inline helper: fill r2u_bb from bb matrix values.
                // Generic over closure to avoid dyn dispatch in inner loop.
                // Uses pre-computed n_inv²·r2u_a to save one multiply per iteration.
                #[inline(always)]
                fn fill_r2u_bb_from<F: Fn(usize, usize) -> f64>(
                    r2u_bb: &mut MatF,
                    c: usize,
                    n_inv_sq_r2u_a: f64,
                    r2u_b: f64,
                    bb_val: F,
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
                        let b_slice = mat_slice_f32(b_mat.as_ref(), 0..gemm_n, 0..c);
                        #[allow(unused_mut)]
                        let mut bb_sl = mat_slice_mut_f32(bb_f32.as_mut(), 0..c, 0..c);
                        let mut _did_gpu = false;
                        #[cfg(feature = "gpu")]
                        {
                            if let Some(ref ctx) = gpu_ctx {
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
                                Par::rayon(0),
                            );
                        }
                        fill_r2u_bb_from(
                            &mut r2u_bb,
                            c,
                            active_n_inv_sq_r2u_a,
                            active_r2u_b,
                            |k, j| bb_f32[(k, j)] as f64,
                        );
                    }
                    GemmBufs::F64 { ref b_mat, .. } => {
                        let b_slice = mat_slice(b_mat.as_ref(), 0..gemm_n, 0..c);
                        #[allow(unused_mut)]
                        let mut bb_sl = mat_slice_mut(bb_f64.as_mut(), 0..c, 0..c);
                        let mut _did_gpu = false;
                        #[cfg(feature = "gpu")]
                        {
                            if let Some(ref ctx) = gpu_ctx {
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
                                            write_gpu_result_f64_native(
                                                &result,
                                                bb_sl.as_mut(),
                                                c,
                                                c,
                                            );
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
                                            ctx.matmul_tn_tiled_flex32(
                                                &b_f32, gemm_n, c, &b_f32, c, tc,
                                            )
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
                                Par::rayon(0),
                            );
                        }
                        fill_r2u_bb_from(
                            &mut r2u_bb,
                            c,
                            active_n_inv_sq_r2u_a,
                            active_r2u_b,
                            |k, j| bb_f64[(k, j)],
                        );
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
                                    a_buf.as_mut().submatrix_mut(0, wi, gemm_n, 1).copy_from(
                                        ring_buf.as_ref().submatrix(0, *slot, gemm_n, 1),
                                    );
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
                                if let Some(ref ctx) = gpu_ctx {
                                    let a_f32 = mat_to_col_major_f32_from_f32(a_view, gemm_n, w);
                                    let b_f32 = mat_to_col_major_f32_from_f32(b_sl, gemm_n, c);
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
                                    Par::rayon(0),
                                );
                            }
                            let ab_view = mat_slice_f32(ab_buf.as_ref(), 0..w, 0..c);
                            // Column-major inner loop with column slices: eliminates faer
                            // tuple-indexing bounds checks (~8.3B checks for full genome).
                            // Pre-computed n_inv²·r2u_a saves one multiply per iteration.
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
                        GemmBufs::F64 {
                            ref mut ring_buf,
                            ref mut b_mat,
                            ref mut a_buf,
                            ref mut ab_buf,
                            ..
                        } => {
                            if !contiguous {
                                for (wi, (_, slot)) in window.iter().enumerate() {
                                    a_buf.as_mut().submatrix_mut(0, wi, gemm_n, 1).copy_from(
                                        ring_buf.as_ref().submatrix(0, *slot, gemm_n, 1),
                                    );
                                }
                            }
                            let a_view = if contiguous {
                                mat_slice(
                                    ring_buf.as_ref(),
                                    0..gemm_n,
                                    first_slot..(first_slot + w),
                                )
                            } else {
                                mat_slice(a_buf.as_ref(), 0..gemm_n, 0..w)
                            };
                            let b_sl = mat_slice(b_mat.as_ref(), 0..gemm_n, 0..c);

                            let mut _did_gpu = false;
                            #[cfg(feature = "gpu")]
                            {
                                if let Some(ref ctx) = gpu_ctx {
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
                                            ctx.matmul_tn_f64(
                                                &a_f64_data,
                                                gemm_n,
                                                w,
                                                &b_f64_data,
                                                c,
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
                                            Err(e) => eprintln!(
                                                "GPU: f64 A×B matmul failed ({e}), falling back to CPU"
                                            ),
                                        }
                                    } else {
                                        // f32 conversion path
                                        let a_f32 =
                                            mat_to_col_major_f32_from_f64(a_view, gemm_n, w);
                                        let b_f32 = mat_to_col_major_f32_from_f64(b_sl, gemm_n, c);
                                        let gpu_result = if let Some(tc) = _gpu_tile_cols {
                                            if _gpu_flex32 {
                                                ctx.matmul_tn_tiled_flex32(
                                                    &a_f32, gemm_n, w, &b_f32, c, tc,
                                                )
                                            } else {
                                                ctx.matmul_tn_tiled(
                                                    &a_f32, gemm_n, w, &b_f32, c, tc,
                                                )
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
                                    Par::rayon(0),
                                );
                            }
                            let ab_view = mat_slice(ab_buf.as_ref(), 0..w, 0..c);
                            // Column-major inner loop with column slices: eliminates faer
                            // tuple-indexing bounds checks. Pre-computed constant saves 1 mul.
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
        } // end exact GEMM else branch

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
    } // end for chunk_start

    if let Some(handle) = reader_handle {
        handle.join().expect("prefetch reader thread panicked");
    }

    if verbose_timing {
        if sketch_dim.is_some() {
            eprintln!(
                "[perf] compute_ldscore: bed_read(stall)={:.3}s norm={:.3}s sketch={:.3}s bb_dot={:.3}s ab_dot={:.3}s ring_store={:.3}s",
                t_bed_read.as_secs_f64(),
                t_norm.as_secs_f64(),
                t_sketch.as_secs_f64(),
                t_bb_dot.as_secs_f64(),
                t_ab_dot.as_secs_f64(),
                t_ring_store.as_secs_f64(),
            );
        } else {
            eprintln!(
                "[perf] compute_ldscore: bed_read(stall)={:.3}s norm={:.3}s bb_dot={:.3}s ab_dot={:.3}s ring_store={:.3}s",
                t_bed_read.as_secs_f64(),
                t_norm.as_secs_f64(),
                t_bb_dot.as_secs_f64(),
                t_ab_dot.as_secs_f64(),
                t_ring_store.as_secs_f64(),
            );
        }
    }

    Ok((l2, maf_per_snp))
}
