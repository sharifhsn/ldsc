//! Web Worker entry points.
//!
//! Three wasm-bindgen functions exposed to `assets/worker.js`:
//!
//! * [`worker_scan_bim`] — reads a BIM Blob via [`FileReaderSync`],
//!   returns the unique chromosomes present + each chr's SNP count.
//!   The main-thread orchestrator uses this to partition work across
//!   the per-chromosome compute worker pool.
//!
//! * [`worker_compute_l2_chrs`] — the per-chunk compute entry point.
//!   Streams BIM (filtered to the assigned `chrs_filter`) + BED via
//!   [`FileReaderSync`], runs `compute_l2_from_bed` with a progress
//!   callback that `postMessage`s back to the main thread mid-compute,
//!   returns the L2 / MAF output for the assigned SNPs.
//!
//! Why both: the orchestrator can't know how to split work across
//! workers without the chr summary first, and we don't want every
//! compute worker to re-stream the full BIM just to discover what
//! chrs exist. Scan once, dispatch many.
//!
//! Streaming via [`FileReaderSync`] (Worker-only) means neither file
//! ever lives in main-thread JS heap (cap ~256-512 MB depending on
//! browser) or in wasm linear memory in full — the working set is
//! the current chunk + ring buffer + output L2 array. BED has no
//! practical size limit; BIM is bounded by the parsed
//! `Vec<BimRecord>` size in wasm linear memory (~60 bytes/record →
//! ~60M SNPs at the 4 GB wasm32 cap).

use std::io::{BufReader, Read, Seek, SeekFrom};

use ldsc::bed::Bed;
use ldsc::l2::{
    L2Config, L2Progress, WindowMode, compute_l2_from_bed, count_fam_str, parse_bim_reader,
};
use ldsc::parse::BimRecord;
use serde::{Deserialize, Serialize};
use wasm_bindgen::JsCast;
use wasm_bindgen::prelude::*;
use web_sys::{Blob, DedicatedWorkerGlobalScope, FileReaderSync};

// ── Read shims over JS Blobs via FileReaderSync ─────────────────────

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

/// `Read`-only shim over a JS [`Blob`] for streaming UTF-8 parsing
/// (BIM, FAM, etc.). Sequential — no `Seek` needed because the BIM
/// parser is single-pass.
///
/// Why it exists: BIM is loaded via `FileReader.readAsText` on the
/// main thread today, which forces the entire BIM into JS string
/// heap (~256-512 MB cap depending on browser → ~4-8M SNPs of BIM).
/// Streaming via `FileReaderSync` inside a Worker bypasses that cap;
/// the only ceiling is the parsed `Vec<BimRecord>` in wasm linear
/// memory (~60 bytes/record → ~60M SNPs at the 4 GB wasm32 cap).
///
/// Reads in 8 MB blocks under the hood (the Blob.slice + read happens
/// once per `Read::read` call up to `buf.len()`); wrap in a `BufReader`
/// so callers paying line-by-line don't trigger a JS round-trip per
/// line.
pub struct BlobBimSource {
    blob: Blob,
    reader: FileReaderSync,
    offset: u64,
    len: u64,
}

impl BlobBimSource {
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

impl Read for BlobBimSource {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        if self.offset >= self.len {
            return Ok(0);
        }
        let want = (self.len - self.offset).min(buf.len() as u64) as usize;
        if want == 0 {
            return Ok(0);
        }
        let slice = self
            .blob
            .slice_with_f64_and_f64(self.offset as f64, (self.offset + want as u64) as f64)
            .map_err(|e| {
                std::io::Error::other(format!("BIM Blob.slice({}, {}): {e:?}", self.offset, want))
            })?;
        let array_buf = self.reader.read_as_array_buffer(&slice).map_err(|e| {
            std::io::Error::other(format!("BIM FileReaderSync.readAsArrayBuffer: {e:?}"))
        })?;
        let array_buf: js_sys::ArrayBuffer = array_buf
            .dyn_into()
            .map_err(|_| std::io::Error::other("BIM FileReaderSync returned non-ArrayBuffer"))?;
        let u8_array = js_sys::Uint8Array::new(&array_buf);
        let n = u8_array.byte_length() as usize;
        let n = n.min(buf.len());
        u8_array.copy_to(&mut buf[..n]);
        self.offset += n as u64;
        Ok(n)
    }
}

// ── Wire format for postMessage between main thread and worker ──────

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

/// Result of a [`worker_compute_l2_chrs`] call. The orchestrator
/// concatenates per-worker outputs back together in BIM order.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WireL2Output {
    /// `bed_idx` of each SNP this worker computed, in input order.
    /// The orchestrator uses this to slot results back into the
    /// global L2 array in BIM order even if work was assigned out
    /// of order.
    pub bed_idx: Vec<usize>,
    pub l2: Vec<f64>,
    pub maf: Vec<f64>,
    pub wall_seconds: f64,
    pub n_snps: usize,
}

/// Per-chr summary from [`worker_scan_bim`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WireChrSummary {
    pub chr: u8,
    pub snp_count: usize,
}

// ── Worker entry points ─────────────────────────────────────────────

/// Quick scan over a BIM Blob: returns the chromosomes present plus
/// each chr's SNP count. Used by the orchestrator to plan how to
/// shard work across the compute worker pool.
///
/// Cost: one streaming pass over the BIM. ~100 MB at biobank-50K
/// scale, ~1 second wall on a single thread.
#[wasm_bindgen]
pub fn worker_scan_bim(bim_file: web_sys::File) -> Result<JsValue, JsError> {
    let blob: Blob = bim_file.into();
    let source =
        BlobBimSource::new(blob).map_err(|e| JsError::new(&format!("BlobBimSource: {e:?}")))?;
    let snps = parse_bim_reader(BufReader::new(source))
        .map_err(|e| JsError::new(&format!("parse BIM: {e:#}")))?;

    // Group consecutive same-chr SNPs into summaries. BIM is
    // expected to be sorted by (chr, bp); we rely on that for cheap
    // counting via a single pass.
    let mut summaries: Vec<WireChrSummary> = Vec::new();
    for snp in &snps {
        match summaries.last_mut() {
            Some(last) if last.chr == snp.chr => last.snp_count += 1,
            _ => summaries.push(WireChrSummary {
                chr: snp.chr,
                snp_count: 1,
            }),
        }
    }
    serde_wasm_bindgen::to_value(&summaries)
        .map_err(|e| JsError::new(&format!("scan output serialise: {e}")))
}

/// Per-chunk compute entry point. Called from `assets/worker.js` per
/// assigned shard of chromosomes. Synchronous — the worker thread
/// blocks until compute completes; progress is forwarded mid-compute
/// via `self.postMessage({ kind: 'progress', ... })` from the
/// callback below.
///
/// `chrs_filter_json` is a JSON-encoded `Vec<u8>` (e.g. `"[1,3,5]"`).
/// An empty array means "all chromosomes" — this is the back-compat
/// path that mirrors today's single-worker behavior.
#[wasm_bindgen]
pub fn worker_compute_l2_chrs(
    bed_file: web_sys::File,
    bim_file: web_sys::File,
    fam_text: String,
    chrs_filter_json: String,
    config_json: String,
) -> Result<JsValue, JsError> {
    let chrs_filter: Vec<u8> = serde_json::from_str(&chrs_filter_json)
        .map_err(|e| JsError::new(&format!("chrs_filter JSON: {e}")))?;
    let chrs_set: Option<std::collections::HashSet<u8>> = if chrs_filter.is_empty() {
        None
    } else {
        Some(chrs_filter.iter().copied().collect())
    };

    web_sys::console::log_1(
        &format!(
            "worker_compute_l2_chrs: chrs={chrs_filter:?} bed.size={} bim.size={} fam.len={}",
            bed_file.size(),
            bim_file.size(),
            fam_text.len(),
        )
        .into(),
    );

    // Stream + parse BIM, filter to assigned chrs.
    let bim_blob: Blob = bim_file.into();
    let bim_source =
        BlobBimSource::new(bim_blob).map_err(|e| JsError::new(&format!("BlobBimSource: {e:?}")))?;
    let all_snps = parse_bim_reader(BufReader::new(bim_source))
        .map_err(|e| JsError::new(&format!("parse BIM: {e:#}")))?;
    let total_bim_snps = all_snps.len();

    let snps: Vec<BimRecord> = match chrs_set {
        Some(ref set) => all_snps
            .into_iter()
            .filter(|s| set.contains(&s.chr))
            .collect(),
        None => all_snps,
    };
    if snps.is_empty() {
        // Empty assignment is harmless — return an empty result so
        // the orchestrator can degrade gracefully.
        let wire = WireL2Output {
            bed_idx: Vec::new(),
            l2: Vec::new(),
            maf: Vec::new(),
            wall_seconds: 0.0,
            n_snps: 0,
        };
        return serde_wasm_bindgen::to_value(&wire)
            .map_err(|e| JsError::new(&format!("output serialise: {e}")));
    }
    let bed_idx: Vec<usize> = snps.iter().map(|s| s.bed_idx).collect();

    let n_indiv = count_fam_str(&fam_text);

    // Open BED with `total_bim_snps` (the *full* sid_count) — the
    // BED file size is validated against `iid_count × total_bim_snps`
    // even though we only read a chr's worth of slabs from it.
    let bed_blob: Blob = bed_file.into();
    let bed_len = bed_blob.size() as u64;
    let bed_source =
        BlobBedSource::new(bed_blob).map_err(|e| JsError::new(&format!("BlobBedSource: {e:?}")))?;
    let bed = Bed::from_source(bed_source, n_indiv, total_bim_snps, bed_len)
        .map_err(|e| JsError::new(&format!("Bed::from_source: {e:#}")))?;

    let cfg: WireL2Config = serde_json::from_str(&config_json)
        .map_err(|e| JsError::new(&format!("config JSON: {e}")))?;

    // Progress callback: postMessage back to main on every chunk.
    // js_sys::global() in a Worker context is the
    // DedicatedWorkerGlobalScope; this is the only way to reach
    // self.postMessage from inside synchronous wasm code.
    let global: DedicatedWorkerGlobalScope = js_sys::global().unchecked_into();
    let on_progress = move |p: L2Progress| {
        let obj = js_sys::Object::new();
        let _ = js_sys::Reflect::set(&obj, &"kind".into(), &"progress".into());
        let _ = js_sys::Reflect::set(&obj, &"chunks_done".into(), &(p.chunks_done as f64).into());
        let _ = js_sys::Reflect::set(
            &obj,
            &"chunks_total".into(),
            &(p.chunks_total as f64).into(),
        );
        let _ = js_sys::Reflect::set(&obj, &"snps_done".into(), &(p.snps_done as f64).into());
        let _ = js_sys::Reflect::set(&obj, &"snps_total".into(), &(p.snps_total as f64).into());
        let _ = js_sys::Reflect::set(&obj, &"chunk_wall_ms".into(), &p.chunk_wall_ms.into());
        // Best-effort post; if it fails (worker shutting down) we
        // just drop the update — compute keeps going.
        let _ = global.post_message(&obj);
    };

    let out = compute_l2_from_bed(bed, snps, n_indiv, cfg.into(), on_progress)
        .map_err(|e| JsError::new(&format!("compute_l2_from_bed: {e:#}")))?;

    let wire = WireL2Output {
        n_snps: out.l2.len(),
        bed_idx,
        l2: out.l2,
        maf: out.maf,
        wall_seconds: out.wall_seconds,
    };
    serde_wasm_bindgen::to_value(&wire).map_err(|e| JsError::new(&format!("output serialise: {e}")))
}
