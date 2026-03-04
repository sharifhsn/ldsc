use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use faer::Mat;

const BED_MAGIC_1: u8 = 0x6c;
const BED_MAGIC_2: u8 = 0x1b;
const BED_MODE_SNP_MAJOR: u8 = 0x01;
const BED_HEADER_LEN: u64 = 3;

pub struct Bed {
    path: PathBuf,
    iid_count: usize,
    sid_count: usize,
    bytes_per_snp: usize,
}

impl Bed {
    pub fn new(bed_path: &str) -> Result<Self> {
        let path = PathBuf::from(bed_path);
        let (fam_path, bim_path) = resolve_companion_paths(&path)?;
        let iid_count = count_lines(&fam_path)
            .with_context(|| format!("counting FAM lines '{}'", fam_path.display()))?;
        let sid_count = count_lines(&bim_path)
            .with_context(|| format!("counting BIM lines '{}'", bim_path.display()))?;
        let bytes_per_snp = iid_count.div_ceil(4);

        let mut f =
            File::open(&path).with_context(|| format!("opening BED file '{}'", path.display()))?;
        let mut header = [0u8; 3];
        f.read_exact(&mut header)
            .with_context(|| format!("reading BED header '{}'", path.display()))?;
        anyhow::ensure!(
            header[0] == BED_MAGIC_1 && header[1] == BED_MAGIC_2,
            "BED file '{}' has invalid magic bytes",
            path.display()
        );
        anyhow::ensure!(
            header[2] == BED_MODE_SNP_MAJOR,
            "BED file '{}' is not SNP-major (mode=1 required)",
            path.display()
        );

        Ok(Self {
            path,
            iid_count,
            sid_count,
            bytes_per_snp,
        })
    }

    pub(crate) fn open_reader(&self) -> Result<File> {
        File::open(&self.path)
            .with_context(|| format!("opening BED file '{}'", self.path.display()))
    }

    pub(crate) fn bytes_per_snp(&self) -> usize {
        self.bytes_per_snp
    }

    pub(crate) fn snp_offset(&self, sid: usize) -> u64 {
        BED_HEADER_LEN + (sid as u64) * self.bytes_per_snp as u64
    }
}

pub(crate) struct ReadOptions;

pub(crate) struct ReadOptionsBuilder {
    sid_index: Option<Vec<isize>>,
    iid_index: Option<Vec<isize>>,
}

pub(crate) struct ReadOptionsTyped<T> {
    sid_index: Option<Vec<isize>>,
    iid_index: Option<Vec<isize>>,
    _marker: std::marker::PhantomData<T>,
}

impl ReadOptions {
    pub fn builder() -> ReadOptionsBuilder {
        ReadOptionsBuilder {
            sid_index: None,
            iid_index: None,
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

    pub fn f32(&self) -> ReadOptionsTyped<f32> {
        ReadOptionsTyped {
            sid_index: self.sid_index.clone(),
            iid_index: self.iid_index.clone(),
            _marker: std::marker::PhantomData,
        }
    }

    #[allow(dead_code)]
    pub fn f64(&self) -> ReadOptionsTyped<f64> {
        ReadOptionsTyped {
            sid_index: self.sid_index.clone(),
            iid_index: self.iid_index.clone(),
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

    pub fn read(&self, bed: &mut Bed) -> Result<Mat<T>> {
        let iid_indices = resolve_indices(self.iid_index.as_deref(), bed.iid_count)?;
        let sid_indices = resolve_indices(self.sid_index.as_deref(), bed.sid_count)?;
        let mut out = Mat::from_fn(iid_indices.len(), sid_indices.len(), |_, _| T::missing());
        if iid_indices.is_empty() || sid_indices.is_empty() {
            return Ok(out);
        }

        let mut f = File::open(&bed.path)
            .with_context(|| format!("opening BED file '{}'", bed.path.display()))?;
        let mut buf = vec![0u8; bed.bytes_per_snp];

        for (col_idx, &sid) in sid_indices.iter().enumerate() {
            let offset = BED_HEADER_LEN + (sid as u64) * bed.bytes_per_snp as u64;
            f.seek(SeekFrom::Start(offset))
                .with_context(|| format!("seeking BED offset {}", offset))?;
            f.read_exact(&mut buf)
                .with_context(|| format!("reading BED SNP {}", sid))?;

            for (row_idx, &iid) in iid_indices.iter().enumerate() {
                let byte = buf[iid / 4];
                let shift = (iid % 4) * 2;
                let bits = (byte >> shift) & 0b11;
                out[(row_idx, col_idx)] = decode_genotype(bits);
            }
        }

        Ok(out)
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
                if idx < 0 {
                    anyhow::bail!("negative index {} is not supported", idx);
                }
                let idx = idx as usize;
                if idx >= max {
                    anyhow::bail!("index {} out of range (max {})", idx, max);
                }
                out.push(idx);
            }
            Ok(out)
        }
    }
}

fn decode_genotype<T: BedVal>(bits: u8) -> T {
    match bits {
        0 => T::from_u8(2),
        1 => T::missing(),
        2 => T::from_u8(1),
        3 => T::from_u8(0),
        _ => T::missing(),
    }
}

pub(crate) trait BedVal: Copy {
    fn from_u8(v: u8) -> Self;
    fn missing() -> Self;
}

impl BedVal for f32 {
    fn from_u8(v: u8) -> Self {
        v as f32
    }

    fn missing() -> Self {
        f32::NAN
    }
}

impl BedVal for f64 {
    fn from_u8(v: u8) -> Self {
        v as f64
    }

    fn missing() -> Self {
        f64::NAN
    }
}
