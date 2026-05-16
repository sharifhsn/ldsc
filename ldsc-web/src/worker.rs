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
            console_error_panic_hook::set_once();
            tracing_wasm::set_as_global_default();
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

        // Fresh Bed handle per chr — Bed is consumed by
        // compute_l2_from_bed. BlobBedSource just wraps a Blob
        // reference (cheap clone in JS, refcounted by the browser).
        let bed_source_chr = BlobBedSource::new(bed_blob.clone())
            .map_err(|e| JsError::new(&format!("BlobBedSource (chr {start}..{end}): {e:?}")))?;
        let bed_chr = Bed::from_source(bed_source_chr, n_indiv, total_bim_snps, bed_len)
            .map_err(|e| JsError::new(&format!("Bed::from_source (chr {start}..{end}): {e:#}")))?;

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
