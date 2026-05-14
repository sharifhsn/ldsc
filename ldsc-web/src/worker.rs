//! Web Worker entry point.
//!
//! Runs `compute_l2` synchronously inside a Worker so we can use
//! `FileReaderSync` (which is only available in Worker contexts) to
//! stream chunked reads out of a `Blob`-backed BED. This is the path
//! that lifts the ~500 MB - 1 GB main-thread-WASM file-upload ceiling
//! all the way up to "however big a `File` your browser hands you" —
//! only the working set (current chunk + ring buffer + output L2
//! array, ~30 MB total at biobank scale) ever lives in wasm linear
//! memory.
//!
//! Wire-up: the main thread (`l2_panel.rs`) spawns a `Worker` from
//! `assets/worker.js`, posts a `WorkerRequest` (BED `File` + parsed
//! BIM/FAM text + serialised `L2Config`), and listens for a
//! `WorkerResponse` (`l2 / maf / wall_seconds`) coming back. The
//! worker.js shim imports the SAME wasm bundle the main thread
//! loaded and calls [`worker_compute_l2`] below.

use std::io::{Read, Seek, SeekFrom};

use ldsc::bed::Bed;
use ldsc::l2::{L2Config, WindowMode, compute_l2_from_bed, count_fam_str};
use ldsc::parse::parse_bim_str;
use serde::{Deserialize, Serialize};
use wasm_bindgen::JsCast;
use wasm_bindgen::prelude::*;
use web_sys::{Blob, FileReaderSync};

/// `Read + Seek` shim that backs onto a JS [`Blob`] via
/// [`FileReaderSync`] for synchronous chunked reads.
///
/// Only constructible inside a Web Worker — `FileReaderSync` is
/// `undefined` on the main thread (the constructor will throw and
/// `FileReaderSync::new()` will return `Err`).
pub struct BlobBedSource {
    blob: Blob,
    reader: FileReaderSync,
    /// Current read position. `Read::read` advances this; `Seek`
    /// rewrites it directly.
    offset: u64,
    /// Total `blob.size()`, cached so we don't round-trip JS for
    /// every read bound check.
    len: u64,
}

impl BlobBedSource {
    /// Construct from a `Blob` (or any subtype like `File`). Returns
    /// `Err` on the main thread because `FileReaderSync` constructor
    /// throws there.
    pub fn new(blob: Blob) -> Result<Self, JsValue> {
        let reader = FileReaderSync::new()?;
        let len = blob.size() as u64;
        Ok(Self {
            blob,
            reader,
            offset: 0,
            len,
        })
    }
}

impl Read for BlobBedSource {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        if self.offset >= self.len {
            return Ok(0);
        }
        let want = (self.len - self.offset).min(buf.len() as u64) as usize;
        if want == 0 {
            return Ok(0);
        }
        // Blob.slice(start, end) is zero-copy on the JS side — it
        // returns a sub-Blob view, no bytes touched until FileReader
        // actually reads.
        let slice = self
            .blob
            .slice_with_f64_and_f64(self.offset as f64, (self.offset + want as u64) as f64)
            .map_err(|e| {
                std::io::Error::other(format!("Blob.slice({}, {}): {e:?}", self.offset, want))
            })?;
        let array_buf = self.reader.read_as_array_buffer(&slice).map_err(|e| {
            std::io::Error::other(format!("FileReaderSync.readAsArrayBuffer: {e:?}"))
        })?;
        let array_buf: js_sys::ArrayBuffer = array_buf
            .dyn_into()
            .map_err(|_| std::io::Error::other("FileReaderSync returned non-ArrayBuffer"))?;
        let u8_array = js_sys::Uint8Array::new(&array_buf);
        let n = u8_array.byte_length() as usize;
        let n = n.min(buf.len());
        u8_array.copy_to(&mut buf[..n]);
        self.offset += n as u64;
        Ok(n)
    }
}

impl Seek for BlobBedSource {
    fn seek(&mut self, pos: SeekFrom) -> std::io::Result<u64> {
        let new = match pos {
            SeekFrom::Start(n) => n as i64,
            SeekFrom::End(n) => self.len as i64 + n,
            SeekFrom::Current(n) => self.offset as i64 + n,
        };
        if new < 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "BlobBedSource: seek before start",
            ));
        }
        self.offset = new as u64;
        Ok(self.offset)
    }
}

// ── Wire format for postMessage between main thread and worker ────────

/// Window-mode picker as a tagged union friendly to JSON.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", content = "value", rename_all = "lowercase")]
pub enum WireWindowMode {
    Kb(f64),
    Snp(usize),
    Cm(f64),
}

impl From<WireWindowMode> for WindowMode {
    fn from(w: WireWindowMode) -> Self {
        match w {
            WireWindowMode::Kb(v) => WindowMode::Kb(v),
            WireWindowMode::Snp(v) => WindowMode::Snp(v),
            WireWindowMode::Cm(v) => WindowMode::Cm(v),
        }
    }
}

/// Subset of `L2Config` that round-trips as JSON. Mirrors the UI
/// knobs the L2 panel exposes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WireL2Config {
    pub mode: WireWindowMode,
    pub chunk_size: usize,
    pub use_f32: bool,
    pub sketch: Option<usize>,
    pub sketch_maf_aware: bool,
    pub snp_level_masking: bool,
    pub yes_really: bool,
    pub pq_exp: Option<f64>,
    pub verbose_timing: bool,
}

impl From<WireL2Config> for L2Config {
    fn from(w: WireL2Config) -> Self {
        L2Config {
            mode: w.mode.into(),
            chunk_size: w.chunk_size,
            use_f32: w.use_f32,
            sketch: w.sketch,
            sketch_maf_aware: w.sketch_maf_aware,
            snp_level_masking: w.snp_level_masking,
            yes_really: w.yes_really,
            pq_exp: w.pq_exp,
            verbose_timing: w.verbose_timing,
        }
    }
}

/// Compute output we ship back to the main thread. We deliberately
/// do NOT include the parsed BIM (the main thread can re-parse from
/// `bim_text` if it ever needs the full SNP records) — l2 / maf are
/// what charts + stats need.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WireL2Output {
    pub l2: Vec<f64>,
    pub maf: Vec<f64>,
    pub wall_seconds: f64,
    pub n_snps: usize,
}

/// Worker-side compute entry point. Called from `assets/worker.js`
/// after the wasm bundle initialises. Synchronous — the worker
/// thread blocks until compute completes, then posts the result.
///
/// `bed_file` is a `web_sys::File` (subtype of `Blob`); we wrap it
/// in a [`BlobBedSource`] for streaming chunked reads. `bim_text`
/// and `fam_text` are the already-loaded text contents of the BIM /
/// FAM files (small — the BIM is the largest at ~44 MB for full
/// genome 1000G, well within the main-thread upload ceiling).
/// `config_json` is a `WireL2Config` serialised as JSON.
#[wasm_bindgen]
pub fn worker_compute_l2(
    bed_file: web_sys::File,
    bim_text: String,
    fam_text: String,
    config_json: String,
) -> Result<JsValue, JsError> {
    web_sys::console::log_1(
        &format!(
            "worker_compute_l2: BED size={} bim={} fam={}",
            bed_file.size(),
            bim_text.len(),
            fam_text.len()
        )
        .into(),
    );

    let snps = parse_bim_str(&bim_text).map_err(|e| JsError::new(&format!("parse BIM: {e:#}")))?;
    let n_indiv = count_fam_str(&fam_text);

    let bed_blob: Blob = bed_file.into();
    let bed_len = bed_blob.size() as u64;
    let source = BlobBedSource::new(bed_blob)
        .map_err(|e| JsError::new(&format!("BlobBedSource (worker only): {e:?}")))?;
    let bed = Bed::from_source(source, n_indiv, snps.len(), bed_len)
        .map_err(|e| JsError::new(&format!("Bed::from_source: {e:#}")))?;

    let cfg: WireL2Config = serde_json::from_str(&config_json)
        .map_err(|e| JsError::new(&format!("config JSON: {e}")))?;

    let out = compute_l2_from_bed(bed, snps, n_indiv, cfg.into())
        .map_err(|e| JsError::new(&format!("compute_l2_from_bed: {e:#}")))?;

    let wire = WireL2Output {
        n_snps: out.l2.len(),
        l2: out.l2,
        maf: out.maf,
        wall_seconds: out.wall_seconds,
    };
    serde_wasm_bindgen::to_value(&wire).map_err(|e| JsError::new(&format!("output serialise: {e}")))
}
