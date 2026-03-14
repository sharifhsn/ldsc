use crate::parse::BimRecord;

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
pub(super) fn get_block_lefts_f64(coords: &[f64], max_dist: f64) -> Vec<usize> {
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
pub(super) fn get_block_lefts_by_chr(
    all_snps: &[BimRecord],
    mode: WindowMode,
) -> (Vec<usize>, bool) {
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

/// For each SNP, compute the rightmost SNP index within `max_dist`.
/// Reverse analogue of `get_block_lefts_f64`.
#[allow(dead_code)]
fn get_block_rights_f64(coords: &[f64], max_dist: f64) -> Vec<usize> {
    let m = coords.len();
    if m == 0 {
        return vec![];
    }
    let mut block_right = vec![m - 1; m];
    let mut j = m - 1;

    for i in (0..m).rev() {
        while j > i && (coords[j] - coords[i]).abs() > max_dist {
            j -= 1;
        }
        block_right[i] = j;
    }

    block_right
}

/// Compute block_right per chromosome (mirror of `get_block_lefts_by_chr`).
/// For each SNP j, block_right[j] = largest k in same chromosome such that
/// |coord_k - coord_j| ≤ max_dist.
#[allow(dead_code)]
pub(super) fn get_block_rights_by_chr(all_snps: &[BimRecord], mode: WindowMode) -> Vec<usize> {
    let m = all_snps.len();
    let mut block_right = vec![0usize; m];

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
                for i in 0..len {
                    block_right[start + i] = start + (i + half).min(len - 1);
                }
            }
            WindowMode::Cm(max_cm) => {
                let coords: Vec<f64> = all_snps[start..end].iter().map(|s| s.cm).collect();
                let local = get_block_rights_f64(&coords, max_cm);
                for (i, right) in local.into_iter().enumerate() {
                    block_right[start + i] = start + right;
                }
            }
            WindowMode::Kb(max_kb) => {
                let coords: Vec<f64> = all_snps[start..end]
                    .iter()
                    .map(|s| s.bp as f64 / 1000.0)
                    .collect();
                let local = get_block_rights_f64(&coords, max_kb);
                for (i, right) in local.into_iter().enumerate() {
                    block_right[start + i] = start + right;
                }
            }
        }

        start = end;
    }

    block_right
}
