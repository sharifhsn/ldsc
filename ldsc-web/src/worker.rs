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
//! `Vec<BimRecord>` size in wasm linear memory (~75 bytes/record
//! once the rsID `String` heap allocation is counted → ~50M SNPs in
//! ONE worker at the 4 GB wasm32 cap).
//!
//! Per-worker BIM memory is *further* bounded by the per-chr filter
//! — each compute worker only keeps the BimRecords for its assigned
//! chromosomes, via [`parse_bim_reader_filtered`]. So 4 workers split
//! the BIM by chr instead of each parsing the full genome.

use std::cell::Cell;
use std::io::{BufReader, Read, Seek, SeekFrom};

use ldsc::bed::Bed;
use ldsc::l2::{
    L2Config, WindowMode, compute_l2_from_bed_with_progress, count_fam_str, parse_bim_reader,
    parse_bim_reader_filtered,
};
use ldsc::parse::BimRecord;
use serde::{Deserialize, Serialize};
use wasm_bindgen::JsCast;
use wasm_bindgen::prelude::*;
use web_sys::{Blob, DedicatedWorkerGlobalScope, FileReaderSync};

// ── Tracing setup (LANDMINE — DO NOT CALL FROM WORKER) ─────────────
//
// History: this helper was an attempt to enable `tracing::info!` from
// within compute Workers (so the lib's `[perf]` lines would reach the
// DevTools console). It TRAPS with `RuntimeError: unreachable` when
// invoked from a Worker thread on the multi-threaded WASM bundle
// (atomics + shared memory). Root cause is undiagnosed but reliably
// reproducible.
//
// Workaround: worker-side `[perf]` data is returned as `WireL2Perf`
// in the `WireL2Output` postMessage payload and logged on the main
// thread (which has tracing-wasm wired in `main()`). See
// `worker_client.rs` for the main-thread logging.
//
// Function kept (vs deleted) so future devs see this paragraph before
// re-wiring tracing in the worker. If you really need this, find out
// why the multi-threaded wasm bundle traps in `set_as_global_default`
// and `console_error_panic_hook::set_once` FIRST.
#[cfg(target_arch = "wasm32")]
#[allow(dead_code)]
fn ensure_worker_tracing() {
    use std::cell::Cell;
    thread_local! {
        static DONE: Cell<bool> = const { Cell::new(false) };
    }
    DONE.with(|d| {
        if !d.get() {
            // Both gated behind the `debug` cargo feature — production
            // builds drop the panic hook + tracing subscriber (panics
            // already abort via `-C panic=abort` + `panic_immediate_abort`).
            #[cfg(feature = "debug")]
            {
                console_error_panic_hook::set_once();
                tracing_wasm::set_as_global_default();
            }
            d.set(true);
        }
    });
}

// ── Manual SAB worker pool (workstream G) ───────────────────────────
//
// `inner_worker_loop` is the wasm-bindgen entry point each spawned
// inner Worker calls (from `assets/inner_worker.js`). It parks on the
// shared atomics in `ldsc::wasm_simd::pool` and runs assigned slices
// of GEMM dispatched by the outer Worker.
//
// `init_inner_pool` is called by the outer Worker (in
// `worker_compute_l2_chrs` below, on first invocation) once the JS
// side has spawned + initialised all N inner Workers. It tells the
// pool how many workers it has, so dispatches partition output cols
// across the right denominator.

/// Inner Worker entry point. Park forever, processing GEMM slices
/// dispatched by the outer Worker via the shared atomics in
/// `ldsc::wasm_simd::pool`.
///
/// `worker_id` is in `1..N` — worker 0 is the outer Worker, which
/// runs synchronously in `parallel_gemm_tn_f32`.
#[cfg(all(target_arch = "wasm32", target_feature = "atomics"))]
#[wasm_bindgen]
pub fn inner_worker_loop(worker_id: usize) {
    ldsc::wasm_simd::pool::worker_loop(worker_id);
}

/// Outer-side: tell the pool how many inner Workers it has. Called
/// once per outer Worker after the JS side reports all inner Workers
/// are ready.
#[cfg(all(target_arch = "wasm32", target_feature = "atomics"))]
#[wasm_bindgen]
pub fn init_inner_pool(n_workers: usize) {
    ldsc::wasm_simd::pool::init(n_workers);
}

/// Returns the wasm linear memory as a JsValue (a `WebAssembly.Memory`
/// object). The outer Worker passes this to each spawned inner Worker
/// so they can `mod.default({ module_or_path, memory })` against the
/// SAME SAB-backed memory rather than creating a fresh one. Without
/// this, atomic ops in inner Workers would hit a different memory
/// from the outer Worker's and the manual SAB pool would deadlock.
#[cfg(all(target_arch = "wasm32", target_feature = "atomics"))]
#[wasm_bindgen]
pub fn wasm_memory() -> JsValue {
    wasm_bindgen::memory()
}

/// Wrap the BIM source in a BufReader with this much capacity. Default
/// 8 KB triggers ~12,800 JS round-trips for a 1.66M-SNP BIM (~50 µs
/// each in Chrome) — ~640 ms wasted on reader chatter alone. 256 KB
/// drops that to ~400 round-trips, ~20 ms.
const BIM_BUFREADER_CAP: usize = 256 * 1024;

/// Throttle progress posts from the worker side. Compute can fire
/// thousands of chunks/sec at biobank scale; posting every one would
/// queue ~250+ messages/sec on the main thread, each triggering a
/// Leptos re-render of the progress bar. 16 ms is one rAF frame at
/// 60 Hz — fast enough for a smooth bar, slow enough to amortize
/// JS message dispatch.
///
/// We always post the FIRST and LAST chunk regardless, so the bar
/// is responsive at the start and reaches 100% at the end.
const PROGRESS_THROTTLE_MS: f64 = 16.0;

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

/// `Read + Seek` shim serving one chromosome's BED bytes from a
/// pre-loaded `Vec<u8>` in wasm memory, plus the 3-byte BED magic
/// header so `Bed::from_source` validation passes.
///
/// Workstream I attack on the I/O bottleneck: one `await
/// blob.slice(chr_start, chr_end).arrayBuffer()` upfront replaces
/// the N (50-200) `FileReaderSync.readAsArrayBuffer` calls the
/// streaming path would issue during compute. Per the [Emscripten
/// WORKERFS perf issue], `Blob.slice() + readAsArrayBuffer()` has
/// 50-100× overhead vs a typed-array cache on Chrome (per-call
/// Mojo data-pipe copy, browser process → renderer).
///
/// File-relative offsets are preserved (`file_len` is the FULL BED
/// length, not the shard) so we don't have to translate `bed_idx`
/// downstream. Reads outside `[chr_data_start, chr_data_start +
/// chr_bytes.len())` error — the per-chr BIM slice fed to compute
/// only references `bed_idx` in this chr's range, so the
/// downstream `ChunkReader` should never seek there.
///
/// [Emscripten WORKERFS perf issue]: https://github.com/emscripten-core/emscripten/issues/6955
pub struct PreloadedShardSource {
    /// `[0x6c, 0x1b, 0x01]` — the BED magic + SNP-major mode bytes.
    /// Served by `read` when `offset < 3`.
    header: [u8; 3],
    /// Pre-loaded chr-shard bytes (no header). Indexed by
    /// `(offset - chr_data_start_byte)` for reads in the shard range.
    chr_bytes: Vec<u8>,
    /// File-relative byte offset where `chr_bytes[0]` lives in the
    /// original BED. `= 3 + first_bed_idx * bytes_per_snp`.
    chr_data_start_byte: u64,
    /// Current seek position (file-relative).
    offset: u64,
    /// Total original-file length, passed unchanged to
    /// `Bed::from_source` for validation.
    file_len: u64,
}

impl PreloadedShardSource {
    pub fn new(chr_bytes: Vec<u8>, chr_data_start_byte: u64, file_len: u64) -> Self {
        Self {
            header: [0x6c, 0x1b, 0x01],
            chr_bytes,
            chr_data_start_byte,
            offset: 0,
            file_len,
        }
    }

    /// Async-load a byte range from a JS `Blob` into a `Vec<u8>` via
    /// `await blob.slice(start, end).arrayBuffer()`. One JS↔wasm
    /// round-trip per chr instead of the 50-200 FRS calls the
    /// streaming path would issue.
    pub async fn load_blob_range(
        blob: &Blob,
        start_byte: u64,
        end_byte: u64,
    ) -> Result<Vec<u8>, JsValue> {
        let slice = blob.slice_with_f64_and_f64(start_byte as f64, end_byte as f64)?;
        let array_buf_js = wasm_bindgen_futures::JsFuture::from(slice.array_buffer()).await?;
        let array_buf: js_sys::ArrayBuffer = array_buf_js
            .dyn_into()
            .map_err(|_| JsValue::from_str("array_buffer() returned non-ArrayBuffer"))?;
        let u8_array = js_sys::Uint8Array::new(&array_buf);
        let mut bytes = vec![0u8; u8_array.byte_length() as usize];
        u8_array.copy_to(&mut bytes);
        Ok(bytes)
    }
}

impl Read for PreloadedShardSource {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        if self.offset >= self.file_len {
            return Ok(0);
        }
        // Serve the 3-byte BED magic header for reads in [0, 3).
        // `Bed::from_source` reads exactly these 3 bytes once before
        // the BufReader is wrapped, so this only fires during init.
        if self.offset < 3 {
            let from = self.offset as usize;
            let avail = 3 - from;
            let want = avail.min(buf.len());
            buf[..want].copy_from_slice(&self.header[from..from + want]);
            self.offset += want as u64;
            return Ok(want);
        }
        // Serve from the cached chr-shard. `Read::read` is allowed
        // to return short, so we cap at the remaining shard length.
        let shard_end = self.chr_data_start_byte + self.chr_bytes.len() as u64;
        if self.offset < self.chr_data_start_byte || self.offset >= shard_end {
            return Err(std::io::Error::other(format!(
                "PreloadedShardSource: read at offset {} outside cached shard \
                 [{}, {}) — caller (ChunkReader) should never seek outside \
                 the chr's bed_idx range",
                self.offset, self.chr_data_start_byte, shard_end,
            )));
        }
        let local_offset = (self.offset - self.chr_data_start_byte) as usize;
        let avail = (shard_end - self.offset) as usize;
        let want = avail.min(buf.len());
        buf[..want].copy_from_slice(&self.chr_bytes[local_offset..local_offset + want]);
        self.offset += want as u64;
        Ok(want)
    }
}

impl Seek for PreloadedShardSource {
    fn seek(&mut self, pos: SeekFrom) -> std::io::Result<u64> {
        let new = match pos {
            SeekFrom::Start(n) => n as i64,
            SeekFrom::End(n) => self.file_len as i64 + n,
            SeekFrom::Current(n) => self.offset as i64 + n,
        };
        if new < 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "PreloadedShardSource: seek before start",
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
            verbose_timing: false,
        }
    }
}

/// Serializable mirror of [`ldsc::l2::L2Perf`]. The lib's struct is
/// serde-free; this is the conversion type so the per-phase
/// breakdown survives the `serde-wasm-bindgen` round-trip from
/// worker to main.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct WireL2Perf {
    pub m: usize,
    pub n_indiv: usize,
    pub bed_read_secs: f64,
    pub norm_secs: f64,
    pub sketch_secs: Option<f64>,
    pub bb_dot_secs: f64,
    pub ab_dot_secs: f64,
    pub ring_store_secs: f64,
}

impl From<ldsc::l2::L2Perf> for WireL2Perf {
    fn from(p: ldsc::l2::L2Perf) -> Self {
        Self {
            m: p.m,
            n_indiv: p.n_indiv,
            bed_read_secs: p.bed_read_secs,
            norm_secs: p.norm_secs,
            sketch_secs: p.sketch_secs,
            bb_dot_secs: p.bb_dot_secs,
            ab_dot_secs: p.ab_dot_secs,
            ring_store_secs: p.ring_store_secs,
        }
    }
}

impl WireL2Perf {
    /// Sum per-chr breakdowns into a per-worker total. Called once
    /// per chr inside `worker_compute_l2_chrs`'s chr-loop.
    pub fn accumulate(&mut self, other: &WireL2Perf) {
        self.m += other.m;
        self.n_indiv = other.n_indiv; // identical across chrs within a worker
        self.bed_read_secs += other.bed_read_secs;
        self.norm_secs += other.norm_secs;
        match (&mut self.sketch_secs, other.sketch_secs) {
            (Some(acc), Some(s)) => *acc += s,
            (None, Some(s)) => self.sketch_secs = Some(s),
            _ => {}
        }
        self.bb_dot_secs += other.bb_dot_secs;
        self.ab_dot_secs += other.ab_dot_secs;
        self.ring_store_secs += other.ring_store_secs;
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
    /// Per-phase wall-time breakdown summed across all chrs this
    /// worker handled. Logged by the main-thread orchestrator's
    /// `compute_l2_chrs_done` handler via `tracing::info` (which IS
    /// wired on main, unlike workers).
    #[serde(default)]
    pub perf: WireL2Perf,
}

/// Per-chr summary from [`worker_scan_bim`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WireChrSummary {
    pub chr: u8,
    pub snp_count: usize,
}

/// Serialisable mirror of [`ldsc::h2::H2Result`] — the canonical
/// 5-tuple working stat geneticists quote in papers:
///   h² (SE), intercept (SE), ratio (SE), mean χ², λ_GC
///
/// Plus diagnostic counters: how many sumstats SNPs joined against
/// the LD scores, and the wall-time spent inside the regression.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WireH2Result {
    pub h2: f64,
    pub h2_se: f64,
    pub intercept: f64,
    pub intercept_se: f64,
    pub mean_chi2: f64,
    pub lambda_gc: f64,
    pub ratio: Option<f64>,
    pub ratio_se: Option<f64>,
    /// Count of SNPs that overlapped between the LD-score input and
    /// the sumstats file (after χ² > max_chi2 filter, if any).
    pub n_snps_used: usize,
    /// Total SNPs in the sumstats file (before overlap with LD).
    pub n_snps_sumstats: usize,
    /// Total SNPs in the LD-score input (= total reference panel size
    /// `M` used by the regression).
    pub m_snps_ref: usize,
    /// Mean N (sample size) across joined SNPs.
    pub mean_n: f64,
    /// Wall time inside the worker, measured with `web_time::Instant`
    /// (not `std::time::Instant` — the latter panics on
    /// wasm32-unknown-unknown; see Workstream L diagnosis).
    pub wall_seconds: f64,
}

impl From<ldsc::h2::H2Result> for WireH2Result {
    fn from(r: ldsc::h2::H2Result) -> Self {
        let (ratio, ratio_se) = match r.ratio {
            Some((r, se)) => (Some(r), Some(se)),
            None => (None, None),
        };
        Self {
            h2: r.h2,
            h2_se: r.h2_se,
            intercept: r.intercept,
            intercept_se: r.intercept_se,
            mean_chi2: r.mean_chi2,
            lambda_gc: r.lambda_gc,
            ratio,
            ratio_se,
            n_snps_used: 0,
            n_snps_sumstats: 0,
            m_snps_ref: 0,
            mean_n: 0.0,
            wall_seconds: 0.0,
        }
    }
}

// ── Worker entry points ─────────────────────────────────────────────

/// Quick scan over a BIM Blob: returns the chromosomes present plus
/// each chr's SNP count. Used by the orchestrator to plan how to
/// shard work across the compute worker pool.
///
/// Cost: one streaming pass over the BIM. ~100 MB at biobank-50K
/// scale, ~1 second wall on a single thread. The full
/// `Vec<BimRecord>` is built during the parse — the cost is bounded
/// by the lifetime of this function call (Rust drops it on return),
/// so it doesn't compound across worker fan-out.
#[wasm_bindgen]
pub fn worker_scan_bim(bim_file: web_sys::File) -> Result<JsValue, JsError> {
    let blob: Blob = bim_file.into();
    let source =
        BlobBimSource::new(blob).map_err(|e| JsError::new(&format!("BlobBimSource: {e:?}")))?;
    let snps = parse_bim_reader(BufReader::with_capacity(BIM_BUFREADER_CAP, source))
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
/// Async — returns a JS `Promise<WireL2Output>`. The JS caller in
/// `worker.js` `await`s it. Async because each chr's BED bytes are
/// pre-loaded via `await blob.slice(...).arrayBuffer()` before
/// `compute_l2_from_bed_with_progress` runs (workstream I).
#[wasm_bindgen]
pub async fn worker_compute_l2_chrs(
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

    // We need TWO things from the BIM that the filtered parse can't
    // give us in one normal pass:
    //
    //   - `total_bim_snps` (the *full* sid_count, used by Bed::from_source
    //     to validate the BED file size against `iid × total + 3`).
    //   - the per-chr SNP records for our assigned chrs.
    //
    // The cheap way: stream the BIM once with a counter that ticks
    // every line and a filter that decides what to keep. The
    // predicate closes over `&mut total` so we get both in one pass
    // without ever materializing the full Vec<BimRecord>.
    let bim_blob: Blob = bim_file.into();
    let bim_source =
        BlobBimSource::new(bim_blob).map_err(|e| JsError::new(&format!("BlobBimSource: {e:?}")))?;
    let mut total_bim_snps: usize = 0;
    let snps: Vec<BimRecord> = parse_bim_reader_filtered(
        BufReader::with_capacity(BIM_BUFREADER_CAP, bim_source),
        |chr| {
            total_bim_snps += 1;
            match chrs_set {
                Some(ref set) => set.contains(&chr),
                None => true,
            }
        },
    )
    .map_err(|e| JsError::new(&format!("parse BIM: {e:#}")))?;

    if snps.is_empty() {
        // Empty assignment is harmless — return an empty result so
        // the orchestrator can degrade gracefully.
        let wire = WireL2Output {
            bed_idx: Vec::new(),
            l2: Vec::new(),
            maf: Vec::new(),
            wall_seconds: 0.0,
            n_snps: 0,
            perf: WireL2Perf::default(),
        };
        return serde_wasm_bindgen::to_value(&wire)
            .map_err(|e| JsError::new(&format!("output serialise: {e}")));
    }
    let bed_idx: Vec<usize> = snps.iter().map(|s| s.bed_idx).collect();
    let n_snps_in = snps.len();
    let total_assigned = snps.len();

    // Group consecutive same-chr SNPs. The BIM filter preserves the
    // input order, and the BIM is sorted by (chr, bp), so each chr
    // appears as a contiguous run within `snps`. We must call
    // `compute_l2_from_bed` ONCE PER CHR rather than once for the
    // whole worker shard, otherwise the chunked GEMM's chunk
    // boundaries straddle chr transitions and compute spurious
    // cross-chr r² inside the within-chunk B^T·B (~12% mean L2
    // inflation observed on full 1000G in default chunked mode).
    // This mirrors the native CLI's per-chr rayon dispatch in
    // `src/l2/mod.rs::run` (the non-`--global-pass` path).
    let mut chr_ranges: Vec<(usize, usize)> = Vec::new();
    {
        let mut start = 0usize;
        while start < snps.len() {
            let chr = snps[start].chr;
            let mut end = start + 1;
            while end < snps.len() && snps[end].chr == chr {
                end += 1;
            }
            chr_ranges.push((start, end));
            start = end;
        }
    }

    let n_indiv = count_fam_str(&fam_text);
    let bed_blob: Blob = bed_file.into();
    let bed_len = bed_blob.size() as u64;
    // Pre-compute bytes_per_snp once so the per-chr loop can
    // calculate each chr's byte range without re-deriving it.
    // PLINK BED layout: 2 bits per individual, rounded UP to the
    // nearest byte. Matches `Bed::from_source`'s computation.
    let bytes_per_snp: u64 = (n_indiv as u64).saturating_add(3) / 4;

    let cfg_wire: WireL2Config = serde_json::from_str(&config_json)
        .map_err(|e| JsError::new(&format!("config JSON: {e}")))?;
    let cfg: L2Config = cfg_wire.into();

    let global: DedicatedWorkerGlobalScope = js_sys::global()
        .dyn_into()
        .map_err(|_| JsError::new("worker_compute_l2_chrs called outside a DedicatedWorker"))?;

    // Throttle + cumulative-progress state lives in `Cell`s shared
    // across per-chr calls. Each per-chr call's L2Progress reports
    // chunks_done relative to its own chr; we add `snps_done_base`
    // (cumulative SNPs from prior chrs in this worker) before
    // posting upstream.
    let last_post_ms: Cell<f64> = Cell::new(0.0);
    let snps_done_base: Cell<usize> = Cell::new(0);

    let mut accum_l2: Vec<f64> = Vec::with_capacity(snps.len());
    let mut accum_maf: Vec<f64> = Vec::with_capacity(snps.len());
    let mut total_wall_seconds = 0.0f64;
    let mut accum_perf: WireL2Perf = WireL2Perf::default();

    // Move snps into a temp Vec we can drain by chr-range. Avoids
    // cloning each chr's slice.
    let mut snps_remaining = snps;

    for &(start, end) in &chr_ranges {
        let chr_n = end - start;
        // Drain the first `chr_n` records (consumes the head of
        // snps_remaining; preserves order).
        let chr_snps: Vec<BimRecord> = snps_remaining.drain(..chr_n).collect();

        // Check whether the chr's SNPs are contiguous in the BED
        // file (i.e. `bed_idx` increases by exactly 1 per SNP). The
        // common case is YES — no `--extract` means BIM rows map
        // straight to BED rows, and we filter BIM by chr in
        // BIM-sorted order. If contiguous, we can pre-load the
        // chr's bytes in ONE `await blob.slice(...).arrayBuffer()`
        // call and serve all subsequent reads from RAM. If NOT
        // contiguous (currently only possible if a future
        // `--extract`-like feature lands), fall back to the
        // streaming `BlobBedSource` path.
        let first_bed_idx = chr_snps[0].bed_idx;
        let last_bed_idx = chr_snps[chr_n - 1].bed_idx;
        let contiguous = (last_bed_idx - first_bed_idx + 1) == chr_n;

        let bed_chr = if contiguous {
            // Async pre-load this chr's BED bytes. One JS↔wasm
            // round-trip per chr instead of ~50-200 (8 MB BufReader
            // serves 50-200 chunks per chr at biobank c=200). The
            // browser fulfils `arrayBuffer()` as a single contiguous
            // Mojo data-pipe copy, which is meaningfully faster
            // than N back-to-back `FileReaderSync` slice reads.
            let chr_data_start = 3u64 + first_bed_idx as u64 * bytes_per_snp;
            let chr_data_end = 3u64 + (last_bed_idx as u64 + 1) * bytes_per_snp;
            let chr_bytes =
                PreloadedShardSource::load_blob_range(&bed_blob, chr_data_start, chr_data_end)
                    .await
                    .map_err(|e| {
                        JsError::new(&format!(
                            "PreloadedShardSource::load_blob_range (chr {start}..{end}): {e:?}",
                        ))
                    })?;
            let source = PreloadedShardSource::new(chr_bytes, chr_data_start, bed_len);
            Bed::from_source(source, n_indiv, total_bim_snps, bed_len).map_err(|e| {
                JsError::new(&format!(
                    "Bed::from_source (preloaded shard, chr {start}..{end}): {e:#}",
                ))
            })?
        } else {
            // Fall back to the streaming BlobBedSource path. Cheap
            // Blob clone (refcounted by the browser); each per-SNP
            // read still pays the FileReaderSync round-trip cost.
            let bed_source_chr = BlobBedSource::new(bed_blob.clone())
                .map_err(|e| JsError::new(&format!("BlobBedSource (chr {start}..{end}): {e:?}")))?;
            Bed::from_source(bed_source_chr, n_indiv, total_bim_snps, bed_len).map_err(|e| {
                JsError::new(&format!("Bed::from_source (chr {start}..{end}): {e:#}"))
            })?
        };

        // Per-chr progress callback wraps the shared throttle +
        // cumulative-offset state. `is_first` / `is_last` semantics
        // are computed against the cumulative counter, not the
        // per-chr one, so the bar advances smoothly across chr
        // boundaries within the worker.
        let on_progress = |p: ldsc::l2::L2Progress| {
            let cumulative = snps_done_base.get() + p.snps_done;
            let now_ms = js_sys::Date::now();
            let is_first = cumulative == p.snps_done && p.chunks_done == 1;
            let is_last = cumulative == total_assigned;
            let is_due = now_ms - last_post_ms.get() >= PROGRESS_THROTTLE_MS;
            if !(is_first || is_last || is_due) {
                return;
            }
            last_post_ms.set(now_ms);
            let obj = js_sys::Object::new();
            let _ = js_sys::Reflect::set(&obj, &"kind".into(), &"progress".into());
            let _ = js_sys::Reflect::set(&obj, &"snps_done".into(), &(cumulative as f64).into());
            if let Err(e) = global.post_message(&obj) {
                web_sys::console::warn_1(
                    &format!("worker progress postMessage failed: {e:?}").into(),
                );
            }
        };

        let chr_out =
            compute_l2_from_bed_with_progress(bed_chr, chr_snps, n_indiv, cfg.clone(), on_progress)
                .map_err(|e| {
                    JsError::new(&format!("compute_l2_from_bed (chr {start}..{end}): {e:#}"))
                })?;

        // Length invariant per chr.
        if chr_out.l2.len() != chr_n || chr_out.maf.len() != chr_n {
            return Err(JsError::new(&format!(
                "internal: chr {start}..{end} length mismatch (l2={}, maf={}, expected={})",
                chr_out.l2.len(),
                chr_out.maf.len(),
                chr_n,
            )));
        }

        accum_l2.extend(chr_out.l2);
        accum_maf.extend(chr_out.maf);
        total_wall_seconds += chr_out.wall_seconds;
        accum_perf.accumulate(&WireL2Perf::from(chr_out.perf));
        snps_done_base.set(snps_done_base.get() + chr_n);
    }

    // Pool-level length invariant — guards against future drift in
    // the lib's L2Output shape that would silently truncate via the
    // triple-zip in worker_client::assemble_partials.
    if bed_idx.len() != accum_l2.len() || accum_l2.len() != accum_maf.len() {
        return Err(JsError::new(&format!(
            "internal: length mismatch in WireL2Output (bed_idx={}, l2={}, maf={}, input snps={})",
            bed_idx.len(),
            accum_l2.len(),
            accum_maf.len(),
            n_snps_in,
        )));
    }
    let out = ldsc::l2::L2Output {
        snps: Vec::new(), // unused on the wire side
        l2: accum_l2,
        maf: accum_maf,
        wall_seconds: total_wall_seconds,
        // Pool-level L2Output uses a default L2Perf because the
        // per-chr perfs were already accumulated into `accum_perf`
        // above (which goes out via `WireL2Output.perf` below).
        perf: ldsc::l2::L2Perf::default(),
    };

    let wire = WireL2Output {
        n_snps: out.l2.len(),
        bed_idx,
        l2: out.l2,
        maf: out.maf,
        wall_seconds: out.wall_seconds,
        perf: accum_perf,
    };
    serde_wasm_bindgen::to_value(&wire).map_err(|e| JsError::new(&format!("output serialise: {e}")))
}

/// Decode a sumstats payload: if it starts with the gzip magic
/// (`0x1f 0x8b`), inflate it via `flate2`; otherwise interpret as
/// UTF-8. Returns the decoded text plus the original compressed size
/// (`Some(z)` if gz, `None` if plain).
///
/// Separated from [`worker_compute_h2`] so it's testable on native —
/// the actual h² entry is wasm-only because it speaks `web_sys` /
/// `Blob`, but this helper is pure-Rust.
pub fn decode_sumstats_bytes(
    bytes: Vec<u8>,
) -> std::result::Result<(String, Option<usize>), String> {
    let is_gzipped = bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b;
    if is_gzipped {
        use std::io::Read;
        let compressed_size = bytes.len();
        let mut decoder = flate2::read::GzDecoder::new(&bytes[..]);
        let mut text = String::new();
        decoder
            .read_to_string(&mut text)
            .map_err(|e| format!("decompress sumstats.gz: {e}"))?;
        Ok((text, Some(compressed_size)))
    } else {
        let text =
            String::from_utf8(bytes).map_err(|e| format!("sumstats not valid UTF-8: {e}"))?;
        Ok((text, None))
    }
}

/// Run LDSC's h² regression in a Web Worker against:
///
/// * `l2` — LD scores per BIM row (output of [`worker_compute_l2_chrs`]).
/// * `maf` — MAF per BIM row, same length as `l2`. Used to compute
///   `M_5_50` (= count of SNPs with MAF in [0.05, 0.5]) — the
///   reference-panel size LDSC's h² normalisation expects. This is
///   the value written to the `.l2.M_5_50` files by the CLI and
///   loaded by `ldsc h2 --ref-ld-chr …`; passing the full BIM count
///   instead overestimates h² by roughly `M_full / M_5_50` (≈ 1.3×
///   on 1000G EUR).
/// * `bed_idx` — `bed_idx` from the same L2 output. Maps L2 array
///   positions back to BIM rows so we can join sumstats by SNP rsID.
///   Same length as `l2`.
/// * `bim_file` — the user's `.bim`. Re-streamed here to get the
///   per-SNP rsID strings (those aren't returned by the L2 worker).
/// * `sumstats_file` — a `.sumstats` Blob, either plain text or
///   gzipped. The header on the bytes (`0x1f 0x8b`) decides; LDSC's
///   `munge_sumstats` writes gzipped by default. Columns:
///   `SNP \t Z \t N` with a header line. Whitespace-tolerant.
/// * `n_blocks` — block jackknife block count. LDSC default = 200.
/// * `two_step` — chi² cutoff for the two-step estimator (None ⇒
///   single-step regression). LDSC default = 30.0.
///
/// Returns [`WireH2Result`] = the canonical 5-tuple plus join
/// diagnostics. All timing uses [`web_time::Instant`] — `std::time::Instant`
/// panics on `wasm32-unknown-unknown` ("time not implemented"), which
/// was the original undiagnosed cause of the rolled-back
/// `worker_h2_demo` trap (see Workstream L probe sweep).
///
/// `async` purely so we can `await blob.text()` on the sumstats Blob
/// — that's the simplest path to get the file contents into Rust
/// without a `FileReaderSync` round-trip per chunk. For typical
/// `.sumstats` files (a few MB; ~1.2M SNPs) this loads in <100 ms.
#[wasm_bindgen]
pub async fn worker_compute_h2(
    l2: Vec<f64>,
    maf: Vec<f64>,
    bed_idx: Vec<u32>,
    bim_file: web_sys::File,
    sumstats_file: web_sys::File,
    n_blocks: usize,
    two_step: Option<f64>,
) -> Result<JsValue, JsError> {
    // CRITICAL: NEVER import `std::time::Instant` here. It panics on
    // wasm32-unknown-unknown with "time not implemented on this
    // platform". Use `web_time::Instant` (already a dep) which routes
    // to `performance.now()` and works inside Web Workers.
    let t_start = web_time::Instant::now();

    if l2.len() != bed_idx.len() || l2.len() != maf.len() {
        return Err(JsError::new(&format!(
            "worker_compute_h2: length mismatch (l2={}, maf={}, bed_idx={})",
            l2.len(),
            maf.len(),
            bed_idx.len(),
        )));
    }
    if l2.is_empty() {
        return Err(JsError::new(
            "worker_compute_h2: L2 input is empty — nothing to regress against",
        ));
    }
    // M_5_50: the count of common SNPs (MAF in [0.05, 0.5]) in the
    // reference panel. LDSC normalises h² by this M, NOT the full BIM
    // count — passing the full count overestimates h² by `M_full /
    // M_5_50` (≈ 1.3× on 1000G EUR). Computed inline from the maf
    // vector the L2 worker already returns.
    let m_5_50: usize = maf.iter().filter(|&&m| m >= 0.05).count();
    if m_5_50 == 0 {
        return Err(JsError::new(
            "worker_compute_h2: no SNPs have MAF ≥ 0.05 — h² regression \
             requires common variants for the M_5_50 normalisation",
        ));
    }

    web_sys::console::log_1(
        &format!(
            "worker_compute_h2: l2.len={} bim.size={} sumstats.size={}",
            l2.len(),
            bim_file.size(),
            sumstats_file.size(),
        )
        .into(),
    );

    // Stream BIM via `FileReaderSync` (Worker-only) and build a map
    // from rsID → position in the `l2` / `bed_idx` arrays.
    //
    // The full BIM has the same row count as the reference panel; we
    // only need the subset whose `bed_idx` matches our `bed_idx` set.
    // Doing one streaming pass with a HashSet of expected bed_idx is
    // cheaper than parsing every row and discarding the rest at
    // biobank scale (1.66M SNPs).
    use std::collections::HashMap;
    let expected_bed_idx: std::collections::HashSet<u32> = bed_idx.iter().copied().collect();
    let bim_blob: Blob = bim_file.into();
    let bim_source =
        BlobBimSource::new(bim_blob).map_err(|e| JsError::new(&format!("BlobBimSource: {e:?}")))?;
    let bim_snps: Vec<BimRecord> = parse_bim_reader_filtered(
        BufReader::with_capacity(BIM_BUFREADER_CAP, bim_source),
        |_chr| true, // need all rows so bed_idx assignment is correct
    )
    .map_err(|e| JsError::new(&format!("parse BIM: {e:#}")))?;
    // Build rsID → position-in-l2 map. `bed_idx_to_l2_pos[bed_idx]`
    // tells us which slot of the `l2` array a given BIM-row maps to.
    let mut bed_idx_to_l2_pos: HashMap<u32, usize> = HashMap::with_capacity(bed_idx.len());
    for (pos, &bi) in bed_idx.iter().enumerate() {
        bed_idx_to_l2_pos.insert(bi, pos);
    }
    let mut rsid_to_l2_pos: HashMap<String, usize> = HashMap::with_capacity(bed_idx.len());
    for snp in &bim_snps {
        let bi = snp.bed_idx as u32;
        if expected_bed_idx.contains(&bi)
            && let Some(&pos) = bed_idx_to_l2_pos.get(&bi)
        {
            rsid_to_l2_pos.insert(snp.snp.clone(), pos);
        }
    }
    web_sys::console::log_1(
        &format!(
            "worker_compute_h2: rsid→pos map built ({} entries)",
            rsid_to_l2_pos.len(),
        )
        .into(),
    );

    // Async-load the sumstats file as raw bytes (need bytes, not a
    // pre-decoded String, so we can detect + decompress gzip).
    // `Blob.arrayBuffer()` returns a Promise<ArrayBuffer> on the
    // worker-thread; await it then copy into a Vec<u8>.
    //
    // We detect gzip by the two-byte magic (`0x1f 0x8b`) at the head
    // of the blob, NOT by filename — users sometimes rename
    // `.sumstats.gz` to `.sumstats` (or vice versa), and `munge_sumstats`
    // canonically writes gzipped. Magic-byte sniffing always wins
    // over file-extension parsing.
    let sumstats_blob: Blob = sumstats_file.into();
    let array_buf_js = wasm_bindgen_futures::JsFuture::from(sumstats_blob.array_buffer())
        .await
        .map_err(|e| JsError::new(&format!("sumstats Blob.arrayBuffer(): {e:?}")))?;
    let array_buf: js_sys::ArrayBuffer = array_buf_js
        .dyn_into()
        .map_err(|_| JsError::new("sumstats Blob.arrayBuffer() returned non-ArrayBuffer"))?;
    let u8_array = js_sys::Uint8Array::new(&array_buf);
    let mut sumstats_bytes = vec![0u8; u8_array.byte_length() as usize];
    u8_array.copy_to(&mut sumstats_bytes);

    let (sumstats_text, compressed_size) = decode_sumstats_bytes(sumstats_bytes)
        .map_err(|e| JsError::new(&format!("sumstats: {e}")))?;
    if let Some(z) = compressed_size {
        web_sys::console::log_1(
            &format!(
                "worker_compute_h2: sumstats loaded ({} compressed bytes → {} chars)",
                z,
                sumstats_text.len(),
            )
            .into(),
        );
    } else {
        web_sys::console::log_1(
            &format!(
                "worker_compute_h2: sumstats loaded ({} chars)",
                sumstats_text.len(),
            )
            .into(),
        );
    }

    // Parse `.sumstats`: header line (SNP, Z, N — order-independent),
    // then tab-separated rows. We accept arbitrary whitespace
    // separators to match LDSC's loose parser.
    let mut snp_col: Option<usize> = None;
    let mut z_col: Option<usize> = None;
    let mut n_col: Option<usize> = None;
    let mut lines = sumstats_text.lines();
    let header = lines
        .next()
        .ok_or_else(|| JsError::new("sumstats: empty file"))?;
    for (i, col) in header.split_whitespace().enumerate() {
        match col.to_ascii_uppercase().as_str() {
            "SNP" | "RSID" => snp_col = Some(i),
            "Z" => z_col = Some(i),
            "N" => n_col = Some(i),
            _ => {}
        }
    }
    let snp_col =
        snp_col.ok_or_else(|| JsError::new("sumstats header missing required column SNP"))?;
    let z_col = z_col.ok_or_else(|| JsError::new("sumstats header missing required column Z"))?;
    let n_col = n_col.ok_or_else(|| JsError::new("sumstats header missing required column N"))?;

    let mut chi2_vec: Vec<f64> = Vec::new();
    let mut l2_vec: Vec<f64> = Vec::new();
    let mut n_vec: Vec<f64> = Vec::new();
    let mut n_sumstats_total: usize = 0;
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split_whitespace().collect();
        let needed = snp_col.max(z_col).max(n_col);
        if cols.len() <= needed {
            continue;
        }
        n_sumstats_total += 1;
        let rsid = cols[snp_col];
        let Some(&pos) = rsid_to_l2_pos.get(rsid) else {
            continue;
        };
        let Ok(z) = cols[z_col].parse::<f64>() else {
            continue;
        };
        let Ok(n) = cols[n_col].parse::<f64>() else {
            continue;
        };
        if !z.is_finite() || !n.is_finite() || n <= 0.0 {
            continue;
        }
        chi2_vec.push(z * z);
        l2_vec.push(l2[pos]);
        n_vec.push(n);
    }

    let n_joined = chi2_vec.len();
    if n_joined < 200 {
        return Err(JsError::new(&format!(
            "worker_compute_h2: only {n_joined} SNPs overlapped between LD scores and sumstats — \
             need at least 200 for jackknife. Check that the sumstats SNP column matches the \
             rsIDs in your BIM file."
        )));
    }

    let mean_n: f64 = n_vec.iter().sum::<f64>() / (n_vec.len() as f64);
    web_sys::console::log_1(
        &format!(
            "worker_compute_h2: joined {} / {} sumstats SNPs (mean N = {:.0}); running regression",
            n_joined, n_sumstats_total, mean_n,
        )
        .into(),
    );

    // Call into the lib's regression. `m_snps` is `M_5_50` (common-
    // variant count from the reference panel) — matches what the CLI
    // reads from `.l2.M_5_50` files. Using the full BIM count instead
    // overestimates h² by `M_full / M_5_50` (verified vs CLI on
    // 1000G_eur + simulated sumstats: identical intercept / mean χ² /
    // λ_GC / ratio across browser and CLI, h² recovers to within
    // M-ratio precision once M_5_50 is used).
    use ldsc::la::col_from_vec;
    let chi2 = col_from_vec(chi2_vec);
    let l2_col = col_from_vec(l2_vec);
    let n_vec_col = col_from_vec(n_vec.clone());
    let n_blocks = n_blocks.max(2).min(n_joined / 2);
    let result = ldsc::h2::run_h2_ldsc(
        &chi2,
        &l2_col,
        &l2_col, // weights LD = ref LD (single-annotation regression)
        &n_vec_col,
        m_5_50 as f64,
        n_blocks,
        two_step,
        None, // unconstrained intercept
    )
    .map_err(|e| JsError::new(&format!("run_h2_ldsc: {e:#}")))?;

    let mut wire = WireH2Result::from(result);
    wire.n_snps_used = n_joined;
    wire.n_snps_sumstats = n_sumstats_total;
    wire.m_snps_ref = m_5_50;
    wire.mean_n = mean_n;
    wire.wall_seconds = t_start.elapsed().as_secs_f64();

    web_sys::console::log_1(
        &format!(
            "worker_compute_h2: done in {:.3}s — h²={:.4} (±{:.4}), intercept={:.4} (±{:.4}), \
             mean χ²={:.4}, λ_GC={:.4}",
            wire.wall_seconds,
            wire.h2,
            wire.h2_se,
            wire.intercept,
            wire.intercept_se,
            wire.mean_chi2,
            wire.lambda_gc,
        )
        .into(),
    );

    serde_wasm_bindgen::to_value(&wire)
        .map_err(|e| JsError::new(&format!("h2 output serialise: {e}")))
}

#[cfg(test)]
mod tests {
    use super::decode_sumstats_bytes;

    /// Plain-text path: returns the bytes verbatim, with `None` for
    /// the compressed-size sentinel.
    #[test]
    fn decode_plain_text_sumstats() {
        let plain = b"SNP\tZ\tN\nrs1\t1.5\t1000\nrs2\t-0.3\t1000\n";
        let (text, compressed) = decode_sumstats_bytes(plain.to_vec()).unwrap();
        assert_eq!(compressed, None);
        assert_eq!(text, std::str::from_utf8(plain).unwrap());
    }

    /// Gzip path: round-trip through flate2 produces the same text
    /// as the plain-text version, with `Some(z)` for the compressed
    /// size. Compresses identical content + verifies decompression.
    #[test]
    fn decode_gzipped_sumstats_roundtrips() {
        use flate2::{Compression, write::GzEncoder};
        use std::io::Write;

        let plain = b"SNP\tZ\tN\nrs1\t1.5\t1000\nrs2\t-0.3\t1000\nrs3\t2.1\t1000\n";
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(plain).unwrap();
        let gz_bytes = encoder.finish().unwrap();

        // Smoke check that the gz wrapping actually fired (first two
        // bytes are the magic).
        assert_eq!(gz_bytes[0], 0x1f);
        assert_eq!(gz_bytes[1], 0x8b);

        let (text, compressed) = decode_sumstats_bytes(gz_bytes.clone()).unwrap();
        assert_eq!(compressed, Some(gz_bytes.len()));
        assert_eq!(text.as_bytes(), &plain[..]);
    }

    /// Truncated gz payload surfaces as an error string, NOT a panic.
    /// The h² UI path catches the error and routes it to `on_error`.
    #[test]
    fn decode_truncated_gzipped_errors() {
        // Gzip magic but no body — decoder errors part-way through.
        let bad = vec![0x1f, 0x8b, 0x08, 0x00];
        let err = decode_sumstats_bytes(bad).unwrap_err();
        assert!(
            err.contains("sumstats.gz"),
            "expected decompress error string, got: {err}"
        );
    }

    /// Empty input: shorter than the 2-byte gz magic check, so
    /// falls through to the plain-text branch and returns an empty
    /// string.
    #[test]
    fn decode_empty_bytes() {
        let (text, compressed) = decode_sumstats_bytes(Vec::new()).unwrap();
        assert_eq!(compressed, None);
        assert_eq!(text, "");
    }

    /// Non-UTF-8 bytes (and no gz magic) surface as an error, not a
    /// panic.
    #[test]
    fn decode_invalid_utf8_errors() {
        let bad = vec![0xfe, 0xfe, 0xfe];
        let err = decode_sumstats_bytes(bad).unwrap_err();
        assert!(err.contains("UTF-8"), "expected utf-8 error, got: {err}");
    }
}
