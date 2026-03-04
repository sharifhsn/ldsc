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
    #[cfg_attr(not(feature = "bed-polars"), allow(dead_code))]
    path: PathBuf,
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

    #[allow(dead_code)]
    pub fn new(bed_path: &str) -> Result<Self> {
        Bed::builder(bed_path).build()
    }

    #[allow(dead_code)]
    pub fn iid_count(&self) -> usize {
        self.iid_count
    }

    #[allow(dead_code)]
    pub fn sid_count(&self) -> usize {
        self.sid_count
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

    #[allow(dead_code)]
    pub fn read_block_bytes(&self) -> usize {
        self.read_block_bytes
    }
}

impl BedBuilder {
    #[allow(dead_code)]
    pub fn stream_buffer_bytes(mut self, bytes: usize) -> Self {
        self.stream_buffer_bytes = bytes.max(1);
        self
    }

    #[allow(dead_code)]
    pub fn read_block_bytes(mut self, bytes: usize) -> Self {
        self.read_block_bytes = bytes.max(1);
        self
    }

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
            path: self.path,
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
    _marker: std::marker::PhantomData<T>,
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

    #[allow(dead_code)]
    pub fn iid_index(&mut self, idx: &[isize]) -> &mut Self {
        self.iid_index = Some(idx.to_vec());
        self
    }

    #[allow(dead_code)]
    pub fn count_a1(&mut self) -> &mut Self {
        self.count_a1 = true;
        self
    }

    #[allow(dead_code)]
    pub fn count_a2(&mut self) -> &mut Self {
        self.count_a1 = false;
        self
    }

    pub fn f32(&self) -> ReadOptionsTyped<f32> {
        ReadOptionsTyped {
            sid_index: self.sid_index.clone(),
            iid_index: self.iid_index.clone(),
            count_a1: self.count_a1,
            missing_value: f32::NAN,
            _marker: std::marker::PhantomData,
        }
    }

    #[allow(dead_code)]
    pub fn f64(&self) -> ReadOptionsTyped<f64> {
        ReadOptionsTyped {
            sid_index: self.sid_index.clone(),
            iid_index: self.iid_index.clone(),
            count_a1: self.count_a1,
            missing_value: f64::NAN,
            _marker: std::marker::PhantomData,
        }
    }

    #[allow(dead_code)]
    pub fn i8(&self) -> ReadOptionsTyped<i8> {
        ReadOptionsTyped {
            sid_index: self.sid_index.clone(),
            iid_index: self.iid_index.clone(),
            count_a1: self.count_a1,
            missing_value: -127,
            _marker: std::marker::PhantomData,
        }
    }
}

impl<T: BedVal> ReadOptionsTyped<T> {
    #[allow(dead_code)]
    pub fn sid_index(&mut self, idx: &[isize]) -> &mut Self {
        self.sid_index = Some(idx.to_vec());
        self
    }

    pub fn iid_index(&mut self, idx: &[isize]) -> &mut Self {
        self.iid_index = Some(idx.to_vec());
        self
    }

    #[allow(dead_code)]
    pub fn missing_value(&mut self, value: T) -> &mut Self {
        self.missing_value = value;
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

#[cfg(feature = "bed-polars")]
pub mod bed_metadata_polars {
    use super::*;
    use polars::prelude::*;

    pub struct BedMetadataPolars {
        pub bim: DataFrame,
        pub fam: DataFrame,
    }

    impl Bed {
        pub fn metadata_polars(&self) -> Result<BedMetadataPolars> {
            let (fam_path, bim_path) = resolve_companion_paths(&self.path)?;
            let bim = read_bim_polars(&bim_path)?;
            let fam = read_fam_polars(&fam_path)?;
            Ok(BedMetadataPolars { bim, fam })
        }
    }

    fn read_bim_polars(path: &Path) -> Result<DataFrame> {
        let f =
            File::open(path).with_context(|| format!("opening BIM file '{}'", path.display()))?;
        let reader = BufReader::new(f);

        let mut chr: Vec<i32> = Vec::new();
        let mut snp: Vec<String> = Vec::new();
        let mut cm: Vec<f64> = Vec::new();
        let mut bp: Vec<i64> = Vec::new();
        let mut a1: Vec<String> = Vec::new();
        let mut a2: Vec<String> = Vec::new();

        for (line_no, line) in reader.lines().enumerate() {
            let line = line.with_context(|| format!("reading BIM line {}", line_no + 1))?;
            let cols: Vec<&str> = line.split_whitespace().collect();
            anyhow::ensure!(
                cols.len() >= 6,
                "BIM line {}: expected 6 columns, got {}",
                line_no + 1,
                cols.len()
            );
            chr.push(cols[0].parse::<i32>()?);
            snp.push(cols[1].to_string());
            cm.push(cols[2].parse::<f64>()?);
            bp.push(cols[3].parse::<i64>()?);
            a1.push(cols[4].to_string());
            a2.push(cols[5].to_string());
        }

        DataFrame::new(vec![
            Series::new("CHR".into(), chr),
            Series::new("SNP".into(), snp),
            Series::new("CM".into(), cm),
            Series::new("BP".into(), bp),
            Series::new("A1".into(), a1),
            Series::new("A2".into(), a2),
        ])
        .map_err(|e| e.into())
    }

    fn read_fam_polars(path: &Path) -> Result<DataFrame> {
        let f =
            File::open(path).with_context(|| format!("opening FAM file '{}'", path.display()))?;
        let reader = BufReader::new(f);

        let mut fid: Vec<String> = Vec::new();
        let mut iid: Vec<String> = Vec::new();
        let mut father: Vec<String> = Vec::new();
        let mut mother: Vec<String> = Vec::new();
        let mut sex: Vec<i32> = Vec::new();
        let mut pheno: Vec<i32> = Vec::new();

        for (line_no, line) in reader.lines().enumerate() {
            let line = line.with_context(|| format!("reading FAM line {}", line_no + 1))?;
            let cols: Vec<&str> = line.split_whitespace().collect();
            anyhow::ensure!(
                cols.len() >= 6,
                "FAM line {}: expected 6 columns, got {}",
                line_no + 1,
                cols.len()
            );
            fid.push(cols[0].to_string());
            iid.push(cols[1].to_string());
            father.push(cols[2].to_string());
            mother.push(cols[3].to_string());
            sex.push(cols[4].parse::<i32>()?);
            pheno.push(cols[5].parse::<i32>()?);
        }

        DataFrame::new(vec![
            Series::new("FID".into(), fid),
            Series::new("IID".into(), iid),
            Series::new("FATHER".into(), father),
            Series::new("MOTHER".into(), mother),
            Series::new("SEX".into(), sex),
            Series::new("PHENO".into(), pheno),
        ])
        .map_err(|e| e.into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_plink_small(dir: &std::path::Path) -> String {
        let prefix = dir.join("toy");
        let bim = prefix.with_extension("bim");
        let fam = prefix.with_extension("fam");
        let bed = prefix.with_extension("bed");

        // .fam: 4 individuals
        let mut f = std::fs::File::create(&fam).unwrap();
        for i in 1..=4 {
            writeln!(f, "F{} I{} 0 0 0 -9", i, i).unwrap();
        }

        // .bim: 2 SNPs
        let mut b = std::fs::File::create(&bim).unwrap();
        writeln!(b, "1\trs1\t0\t100\tA\tG").unwrap();
        writeln!(b, "1\trs2\t0\t200\tA\tG").unwrap();

        // .bed: SNP-major
        let mut bed_f = std::fs::File::create(&bed).unwrap();
        bed_f.write_all(&[0x6C, 0x1B, 0x01]).unwrap();

        // SNP1: genotypes [0,0,0,0]
        let snp1 = [0u8, 0u8, 0u8, 0u8];
        bed_f.write_all(&[pack_genotypes(&snp1)]).unwrap();

        // SNP2: genotypes [0,1,1,2]
        let snp2 = [0u8, 1u8, 1u8, 2u8];
        bed_f.write_all(&[pack_genotypes(&snp2)]).unwrap();

        prefix.to_string_lossy().to_string()
    }

    fn pack_genotypes(gt: &[u8; 4]) -> u8 {
        // PLINK .bed encoding (2-bit, little-endian per individual):
        // 00 = hom major, 01 = missing, 10 = het, 11 = hom minor
        let code = |g: u8| match g {
            0 => 0b00,
            1 => 0b10,
            2 => 0b11,
            _ => 0b01,
        };
        code(gt[0]) | (code(gt[1]) << 2) | (code(gt[2]) << 4) | (code(gt[3]) << 6)
    }

    #[test]
    fn test_negative_indices() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = write_plink_small(dir.path());
        let mut bed = Bed::new(&format!("{}.bed", prefix)).unwrap();
        let val = ReadOptions::builder()
            .iid_index(&[-1, 0])
            .sid_index(&[-1])
            .f64()
            .read(&mut bed)
            .unwrap();
        assert_eq!(val.nrows(), 2);
        assert_eq!(val.ncols(), 1);
    }

    #[test]
    fn test_count_a1_a2() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = write_plink_small(dir.path());
        let mut bed = Bed::new(&format!("{}.bed", prefix)).unwrap();
        let val_a1 = ReadOptions::builder()
            .sid_index(&[0])
            .iid_index(&[0])
            .f64()
            .read(&mut bed)
            .unwrap();
        let mut bed = Bed::new(&format!("{}.bed", prefix)).unwrap();
        let val_a2 = ReadOptions::builder()
            .count_a2()
            .sid_index(&[0])
            .iid_index(&[0])
            .f64()
            .read(&mut bed)
            .unwrap();
        assert_eq!(val_a1[(0, 0)], 2.0);
        assert_eq!(val_a2[(0, 0)], 0.0);
    }

    #[test]
    fn test_contig_vs_random() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = write_plink_small(dir.path());
        let mut bed = Bed::new(&format!("{}.bed", prefix)).unwrap();
        let contig = ReadOptions::builder()
            .sid_index(&[0, 1])
            .f64()
            .read(&mut bed)
            .unwrap();
        let mut bed = Bed::new(&format!("{}.bed", prefix)).unwrap();
        let random = ReadOptions::builder()
            .sid_index(&[1, 0])
            .f64()
            .read(&mut bed)
            .unwrap();
        assert_eq!(contig[(0, 0)], random[(0, 1)]);
        assert_eq!(contig[(0, 1)], random[(0, 0)]);
    }

    #[test]
    fn test_bad_length() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("bad");
        let bim = prefix.with_extension("bim");
        let fam = prefix.with_extension("fam");
        let bed = prefix.with_extension("bed");

        let mut f = std::fs::File::create(&fam).unwrap();
        writeln!(f, "F1 I1 0 0 0 -9").unwrap();

        let mut b = std::fs::File::create(&bim).unwrap();
        writeln!(b, "1\trs1\t0\t100\tA\tG").unwrap();

        let mut bed_f = std::fs::File::create(&bed).unwrap();
        bed_f.write_all(&[0x6C, 0x1B, 0x01]).unwrap();
        bed_f.write_all(&[0u8]).unwrap(); // should be 1 byte per SNP, ok, but add extra mismatch
        bed_f.write_all(&[0u8]).unwrap();

        let err = Bed::new(bed.to_str().unwrap()).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("invalid length"));
    }

    #[test]
    fn test_bad_header() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("bad_header");
        let bim = prefix.with_extension("bim");
        let fam = prefix.with_extension("fam");
        let bed = prefix.with_extension("bed");

        let mut f = std::fs::File::create(&fam).unwrap();
        writeln!(f, "F1 I1 0 0 0 -9").unwrap();

        let mut b = std::fs::File::create(&bim).unwrap();
        writeln!(b, "1\trs1\t0\t100\tA\tG").unwrap();

        let mut bed_f = std::fs::File::create(&bed).unwrap();
        bed_f.write_all(&[0x00, 0x00, 0x01]).unwrap();
        bed_f.write_all(&[0u8]).unwrap();

        let err = Bed::new(bed.to_str().unwrap()).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("invalid magic bytes"));
    }
}
