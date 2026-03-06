use anyhow::{Context, Result};
use faer::Mat;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

const BED_MAGIC_1: u8 = 0x6c;
const BED_MAGIC_2: u8 = 0x1b;
const BED_MODE_SNP_MAJOR: u8 = 0x01;
const BED_HEADER_LEN: u64 = 3;
const DEFAULT_STREAM_BUFFER_BYTES: usize = 8 * 1024 * 1024;
const DEFAULT_READ_BLOCK_BYTES: usize = 32 * 1024 * 1024;

#[derive(Debug)]
pub struct Bed {
    reader: BufReader<File>,
    iid_count: usize,
    sid_count: usize,
    bytes_per_snp: usize,
    read_block_bytes: usize,
    file_len: u64,
}

pub struct BedBuilder {
    path: PathBuf,
    stream_buffer_bytes: usize,
    read_block_bytes: usize,
}

impl Bed {
    pub fn builder(path: impl AsRef<Path>) -> BedBuilder {
        BedBuilder {
            path: path.as_ref().to_path_buf(),
            stream_buffer_bytes: DEFAULT_STREAM_BUFFER_BYTES,
            read_block_bytes: DEFAULT_READ_BLOCK_BYTES,
        }
    }

    pub fn bytes_per_snp(&self) -> usize {
        self.bytes_per_snp
    }

    pub fn read_snp_bytes(&mut self, sid: usize, buf: &mut [u8]) -> Result<()> {
        anyhow::ensure!(
            sid < self.sid_count,
            "SNP index {} out of range (max {})",
            sid,
            self.sid_count
        );
        anyhow::ensure!(
            buf.len() == self.bytes_per_snp,
            "read_snp_bytes expected buffer of {} bytes, got {}",
            self.bytes_per_snp,
            buf.len()
        );
        let offset = (sid as u64)
            .checked_mul(self.bytes_per_snp as u64)
            .and_then(|v| v.checked_add(BED_HEADER_LEN))
            .context("read_snp_bytes offset overflow")?;
        anyhow::ensure!(
            offset + self.bytes_per_snp as u64 <= self.file_len,
            "read_snp_bytes offset {} out of bounds (len {})",
            offset,
            self.file_len
        );
        self.reader
            .seek(SeekFrom::Start(offset))
            .with_context(|| format!("seeking BED offset {}", offset))?;
        self.reader
            .read_exact(buf)
            .with_context(|| format!("reading BED SNP {}", sid))?;
        Ok(())
    }

    pub fn read_snp_block(&mut self, sid_start: usize, count: usize, buf: &mut [u8]) -> Result<()> {
        anyhow::ensure!(
            sid_start < self.sid_count,
            "SNP index {} out of range (max {})",
            sid_start,
            self.sid_count
        );
        let block_end = sid_start
            .checked_add(count)
            .context("read_snp_block index overflow")?;
        anyhow::ensure!(
            block_end <= self.sid_count,
            "SNP block [{}..{}) out of range (max {})",
            sid_start,
            block_end,
            self.sid_count
        );
        let expected = count
            .checked_mul(self.bytes_per_snp)
            .context("read_snp_block size overflow")?;
        anyhow::ensure!(
            buf.len() == expected,
            "read_snp_block expected buffer of {} bytes, got {}",
            expected,
            buf.len()
        );
        let offset = (sid_start as u64)
            .checked_mul(self.bytes_per_snp as u64)
            .and_then(|v| v.checked_add(BED_HEADER_LEN))
            .context("read_snp_block offset overflow")?;
        anyhow::ensure!(
            offset + expected as u64 <= self.file_len,
            "read_snp_block offset {} out of bounds (len {})",
            offset,
            self.file_len
        );
        self.reader
            .seek(SeekFrom::Start(offset))
            .with_context(|| format!("seeking BED offset {}", offset))?;
        self.reader
            .read_exact(buf)
            .with_context(|| format!("reading BED block start {} count {}", sid_start, count))?;
        Ok(())
    }
}

impl BedBuilder {
    pub fn build(self) -> Result<Bed> {
        let (fam_path, bim_path) = resolve_companion_paths(&self.path)?;
        let iid_count = count_lines(&fam_path)
            .with_context(|| format!("counting FAM lines '{}'", fam_path.display()))?;
        let sid_count = count_lines(&bim_path)
            .with_context(|| format!("counting BIM lines '{}'", bim_path.display()))?;

        let mut f = File::open(&self.path)
            .with_context(|| format!("opening BED file '{}'", self.path.display()))?;
        let mut header = [0u8; 3];
        f.read_exact(&mut header)
            .with_context(|| format!("reading BED header '{}'", self.path.display()))?;
        anyhow::ensure!(
            header[0] == BED_MAGIC_1 && header[1] == BED_MAGIC_2,
            "BED file '{}' has invalid magic bytes",
            self.path.display()
        );
        anyhow::ensure!(
            header[2] == BED_MODE_SNP_MAJOR,
            "BED file '{}' is not SNP-major (mode=1 required)",
            self.path.display()
        );

        let bytes_per_snp_u64 = (iid_count as u64)
            .checked_add(3)
            .context("iid_count overflow")?
            / 4;
        let expected_len = bytes_per_snp_u64
            .checked_mul(sid_count as u64)
            .and_then(|v| v.checked_add(BED_HEADER_LEN))
            .context("BED file length overflow")?;
        let file_len = f
            .metadata()
            .with_context(|| format!("reading BED metadata '{}'", self.path.display()))?
            .len();
        anyhow::ensure!(
            file_len == expected_len,
            "BED file '{}' has invalid length (expected {}, got {})",
            self.path.display(),
            expected_len,
            file_len
        );

        let bytes_per_snp =
            usize::try_from(bytes_per_snp_u64).context("bytes_per_snp does not fit in usize")?;

        let reader = BufReader::with_capacity(self.stream_buffer_bytes, f);

        Ok(Bed {
            reader,
            iid_count,
            sid_count,
            bytes_per_snp,
            read_block_bytes: self.read_block_bytes,
            file_len,
        })
    }
}

pub(crate) struct ReadOptions;

pub(crate) struct ReadOptionsBuilder {
    sid_index: Option<Vec<isize>>,
    iid_index: Option<Vec<isize>>,
    count_a1: bool,
}

pub(crate) struct ReadOptionsTyped<T> {
    sid_index: Option<Vec<isize>>,
    iid_index: Option<Vec<isize>>,
    count_a1: bool,
    missing_value: T,
}

impl ReadOptions {
    pub fn builder() -> ReadOptionsBuilder {
        ReadOptionsBuilder {
            sid_index: None,
            iid_index: None,
            count_a1: true,
        }
    }
}

impl ReadOptionsBuilder {
    pub fn sid_index(&mut self, idx: &[isize]) -> &mut Self {
        self.sid_index = Some(idx.to_vec());
        self
    }

    pub fn f32(&self) -> ReadOptionsTyped<f32> {
        ReadOptionsTyped {
            sid_index: self.sid_index.clone(),
            iid_index: self.iid_index.clone(),
            count_a1: self.count_a1,
            missing_value: f32::NAN,
        }
    }
}

impl<T: BedVal> ReadOptionsTyped<T> {
    pub fn iid_index(&mut self, idx: &[isize]) -> &mut Self {
        self.iid_index = Some(idx.to_vec());
        self
    }

    pub fn read(&self, bed: &mut Bed) -> Result<Mat<T>> {
        let iid_indices = resolve_indices(self.iid_index.as_deref(), bed.iid_count)?;
        let sid_indices = resolve_indices(self.sid_index.as_deref(), bed.sid_count)?;
        let mut out = Mat::from_fn(iid_indices.len(), sid_indices.len(), |_, _| {
            self.missing_value
        });
        if iid_indices.is_empty() || sid_indices.is_empty() {
            return Ok(out);
        }
        self.read_into(bed, &iid_indices, &sid_indices, &mut out)?;
        Ok(out)
    }

    pub fn read_into(
        &self,
        bed: &mut Bed,
        iid_indices: &[usize],
        sid_indices: &[usize],
        out: &mut Mat<T>,
    ) -> Result<()> {
        anyhow::ensure!(
            out.nrows() == iid_indices.len() && out.ncols() == sid_indices.len(),
            "output matrix shape mismatch: expected {}x{}, got {}x{}",
            iid_indices.len(),
            sid_indices.len(),
            out.nrows(),
            out.ncols()
        );

        let iid_positions = precompute_iid_positions(iid_indices);
        let lut = build_lut(self.count_a1, self.missing_value);
        let bytes_per_snp = bed.bytes_per_snp;

        let is_contiguous = sid_indices
            .iter()
            .enumerate()
            .skip(1)
            .all(|(i, &v)| v == sid_indices[0] + i);

        if is_contiguous {
            let mut block_buf: Vec<u8> = Vec::new();
            let max_sids_per_block = (bed.read_block_bytes / bytes_per_snp).max(1);
            let mut col_offset = 0usize;
            while col_offset < sid_indices.len() {
                let block_sids = (sid_indices.len() - col_offset).min(max_sids_per_block);
                let start_sid = sid_indices[col_offset];
                let total_bytes = block_sids * bytes_per_snp;
                block_buf.resize(total_bytes, 0u8);
                bed.read_snp_block(start_sid, block_sids, &mut block_buf)?;
                for j in 0..block_sids {
                    let start = j * bytes_per_snp;
                    let bytes = &block_buf[start..start + bytes_per_snp];
                    decode_column(bytes, &iid_positions, &lut, out, col_offset + j);
                }
                col_offset += block_sids;
            }
            return Ok(());
        }

        let mut buf = vec![0u8; bytes_per_snp];
        for (col_idx, &sid) in sid_indices.iter().enumerate() {
            bed.read_snp_bytes(sid, &mut buf)?;
            decode_column(&buf, &iid_positions, &lut, out, col_idx);
        }

        Ok(())
    }
}

fn resolve_companion_paths(bed_path: &Path) -> Result<(PathBuf, PathBuf)> {
    let bed_str = bed_path.to_string_lossy();
    let base = bed_str.strip_suffix(".bed").unwrap_or(&bed_str);
    let fam = PathBuf::from(format!("{}.fam", base));
    let bim = PathBuf::from(format!("{}.bim", base));
    Ok((fam, bim))
}

fn count_lines(path: &Path) -> Result<usize> {
    let f = File::open(path).with_context(|| format!("opening '{}'", path.display()))?;
    Ok(BufReader::new(f).lines().count())
}

fn resolve_indices(raw: Option<&[isize]>, max: usize) -> Result<Vec<usize>> {
    match raw {
        None => Ok((0..max).collect()),
        Some(list) => {
            let mut out = Vec::with_capacity(list.len());
            for &idx in list {
                let resolved = if idx >= 0 {
                    let u = idx as usize;
                    anyhow::ensure!(u < max, "index {} out of range (max {})", idx, max);
                    u
                } else {
                    let neg = (-idx) as usize;
                    anyhow::ensure!(neg <= max, "index {} out of range (max {})", idx, max);
                    max - neg
                };
                out.push(resolved);
            }
            Ok(out)
        }
    }
}

#[derive(Clone, Copy)]
struct IidPos {
    byte_idx: usize,
    shift: u8,
}

fn precompute_iid_positions(iid_indices: &[usize]) -> Vec<IidPos> {
    iid_indices
        .iter()
        .map(|&iid| IidPos {
            byte_idx: iid / 4,
            shift: ((iid % 4) * 2) as u8,
        })
        .collect()
}

fn build_lut<T: BedVal>(count_a1: bool, missing_value: T) -> [T; 4] {
    if count_a1 {
        [T::from_u8(2), missing_value, T::from_u8(1), T::from_u8(0)]
    } else {
        [T::from_u8(0), missing_value, T::from_u8(1), T::from_u8(2)]
    }
}

fn decode_column<T: BedVal>(
    bytes: &[u8],
    iid_positions: &[IidPos],
    lut: &[T; 4],
    out: &mut Mat<T>,
    col: usize,
) {
    for (row_idx, pos) in iid_positions.iter().enumerate() {
        let byte = bytes[pos.byte_idx];
        let bits = (byte >> pos.shift) & 0b11;
        out[(row_idx, col)] = lut[bits as usize];
    }
}

pub(crate) trait BedVal: Copy {
    fn from_u8(v: u8) -> Self;
}

impl BedVal for f32 {
    fn from_u8(v: u8) -> Self {
        v as f32
    }
}

impl BedVal for f64 {
    fn from_u8(v: u8) -> Self {
        v as f64
    }
}

impl BedVal for i8 {
    fn from_u8(v: u8) -> Self {
        v as i8
    }
}
