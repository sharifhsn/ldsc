use crate::bed::Bed;
use crate::parse::BimRecord;
use anyhow::{Context, Result};

pub(super) fn compute_snp_stats(
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
