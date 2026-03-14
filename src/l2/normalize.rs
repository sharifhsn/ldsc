/// AVX2+FMA vectorized sum and sum-of-squares for f32 slices.
/// Returns (sum, sum_sq) as f64. BED values are {0,1,2,NaN} so f32
/// accumulation is exact (max partial sum < 2^23).
///
/// # Safety
/// Requires AVX2 and FMA CPU features to be available.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
pub(super) unsafe fn sum_sumsq_f32(data: &[f32]) -> (f64, f64) {
    use std::arch::x86_64::*;
    unsafe {
        let n = data.len();
        let n8 = n / 8 * 8;
        let ptr = data.as_ptr();

        let mut sum_v = _mm256_setzero_ps();
        let mut sq_v = _mm256_setzero_ps();

        let mut i = 0usize;
        while i < n8 {
            let v = _mm256_loadu_ps(ptr.add(i));
            sum_v = _mm256_add_ps(sum_v, v);
            sq_v = _mm256_fmadd_ps(v, v, sq_v); // v*v + sq_v
            i += 8;
        }

        // Horizontal reduction: ymm → scalar f32
        let lo = _mm256_castps256_ps128(sum_v);
        let hi = _mm256_extractf128_ps(sum_v, 1);
        let sum4 = _mm_add_ps(lo, hi);
        let sum4_hi = _mm_movehl_ps(sum4, sum4);
        let sum2 = _mm_add_ps(sum4, sum4_hi);
        let sum1 = _mm_add_ss(sum2, _mm_shuffle_ps(sum2, sum2, 1));
        let mut sum = _mm_cvtss_f32(sum1) as f64;

        let lo = _mm256_castps256_ps128(sq_v);
        let hi = _mm256_extractf128_ps(sq_v, 1);
        let sq4 = _mm_add_ps(lo, hi);
        let sq4_hi = _mm_movehl_ps(sq4, sq4);
        let sq2 = _mm_add_ps(sq4, sq4_hi);
        let sq1 = _mm_add_ss(sq2, _mm_shuffle_ps(sq2, sq2, 1));
        let mut sum_sq = _mm_cvtss_f32(sq1) as f64;

        // Scalar remainder
        for j in n8..n {
            let v = *ptr.add(j);
            sum += v as f64;
            sum_sq += (v * v) as f64;
        }

        (sum, sum_sq)
    }
}

/// Fallback for non-x86_64 targets.
#[cfg(not(target_arch = "x86_64"))]
pub(super) unsafe fn sum_sumsq_f32(data: &[f32]) -> (f64, f64) {
    let mut sum = 0f64;
    let mut sum_sq = 0f64;
    for &v in data {
        let vf = v as f64;
        sum += vf;
        sum_sq += vf * vf;
    }
    (sum, sum_sq)
}

/// Normalize genotype column in-place (impute NaN→mean, centre, scale). Returns MAF.
/// When pre-computed sum/count/sum_sq are available (fused with copy), does center+scale
/// in a single pass using variance from raw moments: var = E[X²] - E[X]².
pub(super) fn normalize_col_f64_with_stats(
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
pub(super) fn normalize_col_f32_with_stats(
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
    let inv_std_f32 = if std > 0.0 {
        (1.0 / std) as f32
    } else {
        0.0f32
    };
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
