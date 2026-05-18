//! Main-thread Web Worker pool client.
//!
//! Orchestrates a per-chromosome worker pool against the user's
//! uploaded BED + BIM + FAM. Single entry point [`spawn_compute_l2_pool`]:
//!
//! 1. Spawn a one-off scan worker, post the BIM File handle. Worker
//!    streams BIM via `FileReaderSync`, returns the chr summary
//!    (`Vec<{chr, snp_count}>`).
//! 2. Partition the chr list across `min(distinct_chrs,
//!    navigator.hardwareConcurrency, MAX_POOL)` compute workers,
//!    largest-first / smallest-residual to balance load. Empty
//!    buckets are dropped before workers are spawned.
//! 3. Spawn the compute workers in parallel. Each worker streams
//!    the BIM (filtered to its assigned chrs) + BED via
//!    `BlobBedSource`, runs `compute_l2_from_bed_with_progress`
//!    with a callback that `postMessage`s back mid-compute (throttled
//!    to ~1 per 16 ms, see worker.rs). Each worker is monitored by
//!    a per-pool watchdog setInterval that aborts the run if no
//!    progress for [`WATCHDOG_TIMEOUT_MS`].
//! 4. Aggregate per-worker results into one [`WorkerComputeResult`]
//!    with SNPs concatenated in BIM order (sorted by `bed_idx`).
//!    Wall-time is `max` across workers (parallel), progress events
//!    fold into one global `(snps_done, snps_total)` counter.
//!
//! Lifecycle: workers are terminated as soon as their final
//! `compute_l2_chrs_done` message arrives. On any error path
//! (worker.onerror, parse failure, timeout) all sibling workers
//! are terminated and `on_error` fires. A `Rc<Cell<bool>>` "completed"
//! flag enforces that `on_done` and `on_error` are mutually
//! exclusive — without it, a `worker.terminate()` cleanup that
//! triggers an `onerror` could fire `on_error` after `on_done`
//! has already fired.

use std::cell::{Cell, RefCell};
use std::rc::Rc;

use js_sys::Reflect;
use wasm_bindgen::JsCast;
use wasm_bindgen::prelude::*;
use web_sys::{MessageEvent, Url, Worker, WorkerOptions, WorkerType};

/// Hard cap on pool size, regardless of `navigator.hardwareConcurrency`.
///
/// Each outer worker is its own JS Web Worker with its own SAB-backed
/// wasm linear memory (~150 MB at biobank N=50K, chunk=200, f32; less
/// at 1000G scale). 8 outer × ~150 MB = ~1.2 GB renderer SAB total —
/// fits comfortably under Chrome's per-page address-space budget but
/// closer to it than the previous 4-worker cap.
///
/// Why bump it: the fused CountSketch scatter phase is **single-
/// threaded per outer worker** (rayon is a no-op on wasm32 with our
/// current toolchain — no spindle/atomic-wait), so each outer worker
/// runs scatter serially across its chr-shard. On biobank d=200 the
/// scatter is 60s/worker out of ~70s wall — 4× scatter parallelism
/// hides 14 cores on M5 Pro. Bumping to 8 outers (2× more scatter
/// parallelism) drops biobank wall by an estimated 25-30s without
/// touching the kernel.
///
/// `inner_threads_for_pool` will reduce per-outer rayon threads to
/// keep `outer × inner ≤ navigator.hardwareConcurrency` and avoid
/// CPU oversubscription. On HC=18 (M5 Pro) with 8 outers, we get
/// 2 inner threads each (16 total = no oversubscription).
const MAX_POOL: usize = 8;

/// Per-worker watchdog: if no progress / no done message arrives in
/// this window, abort the pool. Catches hung workers (silent wasm
/// panic that doesn't trigger `onerror`, browser tab throttled into
/// the background, structured-clone of a >2 GB File rejected).
/// 120 s is generous — biobank-scale exact f64 chr1 is ~30 s on a
/// single thread.
const WATCHDOG_TIMEOUT_MS: f64 = 120_000.0;

/// How often the watchdog polls last-activity timestamps.
const WATCHDOG_POLL_MS: i32 = 5_000;

// Type aliases for the once-cell callbacks the orchestrator passes
// around. Boxed so we can erase the concrete generic type.
type DoneCb = Rc<RefCell<Option<Box<dyn FnOnce(WorkerComputeResult)>>>>;
type ErrorCb = Rc<RefCell<Option<Box<dyn FnOnce(String)>>>>;
type ProgressCb = Rc<RefCell<dyn FnMut(PoolProgress)>>;
type AbortCb = Rc<dyn Fn(String)>;

/// Result-of-pool payload sent to `on_done`. Same shape as the
/// single-worker result was; SNPs are in BIM order.
#[derive(Clone, Debug)]
pub struct WorkerComputeResult {
    /// Per-SNP `bed_idx` (original BIM row index) in the same order
    /// as `l2`/`maf`. Used by the h² worker to map SNPs back to BIM
    /// rsIDs for sumstats join. Sorted ascending by construction
    /// (BIM order) since `assemble_partials` sorts by `bed_idx`.
    pub bed_idx: Vec<u32>,
    pub l2: Vec<f64>,
    pub maf: Vec<f64>,
    /// Wall time = max across workers (parallel; this is the
    /// user-perceived elapsed compute, not the sum of CPU work).
    pub wall_seconds: f64,
    pub n_snps: usize,
    /// Per-outer-worker timing breakdown, in worker index order
    /// (`worker[0]`, `worker[1]`, …). Used by the L2 panel's results
    /// card to surface the I/O / sketch / GEMM split without having
    /// to dig into DevTools. Same data the worker_client.rs
    /// `tracing::info!` line logs to console — both are derived
    /// from each worker's `WireL2Output.perf`.
    pub per_worker_perf: Vec<crate::worker::WireL2Perf>,
    /// Each outer worker's `wall_seconds` (max-of-chrs inside that
    /// worker), in the same order as `per_worker_perf`. The pool
    /// wall is `max(per_worker_wall_seconds)`; per-worker numbers
    /// expose the long-tail (typically chr1's outer).
    pub per_worker_wall_seconds: Vec<f64>,
}

/// Mid-compute progress event aggregated across the pool.
#[derive(Clone, Copy, Debug)]
pub struct PoolProgress {
    pub snps_done: usize,
    pub snps_total: usize,
}

/// Find the trunk-emitted wasm bundle paths from the current
/// document. Trunk inserts `<link rel="modulepreload"
/// href="/ldsc/ldsc-web-{hash}.js">` and `<link rel="preload"
/// href="/ldsc/ldsc-web-{hash}_bg.wasm">` after the build, both
/// hash-suffixed; we look those up by suffix match.
pub fn discover_wasm_urls() -> Result<(String, String), JsValue> {
    let win = web_sys::window().ok_or_else(|| JsValue::from_str("no window"))?;
    let doc = win
        .document()
        .ok_or_else(|| JsValue::from_str("no document"))?;
    let head = doc.head().ok_or_else(|| JsValue::from_str("no <head>"))?;
    let children = head.children();

    let mut js_url: Option<String> = None;
    let mut bg_url: Option<String> = None;
    for i in 0..children.length() {
        let Some(el) = children.item(i) else {
            continue;
        };
        if el.tag_name().eq_ignore_ascii_case("link") {
            let href = el.get_attribute("href").unwrap_or_default();
            if href.contains("ldsc-web-") {
                if href.ends_with(".js") {
                    js_url = Some(href);
                } else if href.ends_with("_bg.wasm") {
                    bg_url = Some(href);
                }
            }
        }
    }

    Ok((
        js_url.ok_or_else(|| JsValue::from_str("ldsc-web .js bundle not found in <head> links"))?,
        bg_url.ok_or_else(|| JsValue::from_str("ldsc-web _bg.wasm not found in <head> links"))?,
    ))
}

/// `navigator.hardwareConcurrency`, clamped into `[1, MAX_POOL]`.
/// Falls back to 2 if the browser doesn't report a value (Safari in
/// some cross-origin-isolated contexts; older mobile browsers).
fn pool_cap() -> usize {
    let nav_hwc = web_sys::window()
        .map(|w| w.navigator().hardware_concurrency() as usize)
        .unwrap_or(2);
    nav_hwc.clamp(1, MAX_POOL)
}

/// Hard cap on rayon threads inside each compute Worker. 4 because
/// the wasm32-simd128 GEMM kernel saturates IPC on a single core
/// pretty quickly, and the marginal value of the 5th+ rayon thread
/// is below the cost of one more nested Web Worker (each carries a
/// SAB-shared wasm linear memory + JS heap).
const MAX_INNER_THREADS: usize = 4;

/// Returns the number of rayon threads to spawn inside each compute
/// Worker. Total active cores = `outer_pool_size * inner_threads`,
/// which we cap at the device's `navigator.hardwareConcurrency` to
/// avoid CPU oversubscription (the OS scheduler handles minor
/// over-spawn fine, but pinning all cores hurts UI responsiveness).
///
/// Examples (assuming MAX_POOL=4, MAX_INNER_THREADS=4):
///   HC=2  → outer=2, inner=1, total=2
///   HC=4  → outer=4, inner=1, total=4
///   HC=8  → outer=4, inner=2, total=8
///   HC=16 → outer=4, inner=4, total=16
fn inner_threads_for_pool(outer_pool_size: usize) -> usize {
    let nav_hwc = web_sys::window()
        .map(|w| w.navigator().hardware_concurrency() as usize)
        .unwrap_or(2);
    (nav_hwc / outer_pool_size.max(1)).clamp(1, MAX_INNER_THREADS)
}

/// Partition `(chr, snp_count)` pairs across at most `n_workers`
/// non-empty buckets, largest-first / smallest-residual. Returns
/// only the buckets that received at least one chr; the orchestrator
/// uses the result length as the actual worker count.
///
/// Greedy LPT (longest processing time first) is within 4/3 of
/// optimal — fine for a 22-chr / ≤4-worker partition.
fn partition_chrs(summaries: &[(u8, usize)], n_workers: usize) -> Vec<Vec<u8>> {
    let n_workers = n_workers.max(1).min(summaries.len().max(1));
    let mut indexed: Vec<(u8, usize)> = summaries.to_vec();
    // Sort by SNP count descending — std::cmp::Reverse turns the
    // ascending sort_by_key into a descending one and silences the
    // newer nightly's `clippy::unnecessary_sort_by` lint.
    indexed.sort_by_key(|t| std::cmp::Reverse(t.1));

    let mut buckets: Vec<(usize, Vec<u8>)> = (0..n_workers).map(|_| (0, Vec::new())).collect();
    for (chr, n) in indexed {
        let (idx, _) = buckets
            .iter()
            .enumerate()
            .min_by_key(|(_, (load, _))| *load)
            .unwrap();
        buckets[idx].0 += n;
        buckets[idx].1.push(chr);
    }
    buckets
        .into_iter()
        .filter_map(|(_, chrs)| if chrs.is_empty() { None } else { Some(chrs) })
        .collect()
}

/// Spawn a per-chr Web Worker pool to compute L2 against the supplied
/// inputs. Calls `on_progress` periodically with running totals,
/// `on_done` exactly once with the assembled output on success, or
/// `on_error` exactly once on any failure (mutually exclusive with
/// `on_done`).
pub fn spawn_compute_l2_pool(
    bed_file: web_sys::File,
    bim_file: web_sys::File,
    fam_text: String,
    config_json: String,
    on_progress: impl FnMut(PoolProgress) + 'static,
    on_done: impl FnOnce(WorkerComputeResult) + 'static,
    on_error: impl FnOnce(String) + 'static,
) -> Result<(), JsValue> {
    let (wasm_js_url, wasm_bg_url) = discover_wasm_urls()?;
    let wasm_js_url_abs = absolutize(&wasm_js_url)?;
    let wasm_bg_url_abs = absolutize(&wasm_bg_url)?;
    tracing::info!(
        "spawn_compute_l2_pool: js={wasm_js_url} bg={wasm_bg_url} bed.size={} bim.size={} fam.len={}",
        bed_file.size(),
        bim_file.size(),
        fam_text.len(),
    );

    // ── Shared state across all worker callbacks ────────────────────
    //
    // `on_progress` fires many times (per chunk per worker, throttled
    // worker-side to ~1 per 16 ms). `on_done` and `on_error` fire
    // *at most once each, mutually exclusive*. The single-shot
    // callbacks live in `Option<Box<...>>` inside RefCells so the
    // firing path can `take()` them by value (FnOnce).
    //
    // The `completed` flag is the source of truth for "we are done"
    // — both the success and failure paths consult it via
    // `replace(true)` before firing their callback. This closes the
    // race where the success path takes `on_done`, then a stray
    // `worker.onerror` (e.g. from the terminate() teardown of one
    // of the sibling workers) takes `on_error` and fires it as well.
    let completed: Rc<Cell<bool>> = Rc::new(Cell::new(false));
    let on_progress: ProgressCb = Rc::new(RefCell::new(on_progress));
    let on_done: DoneCb = Rc::new(RefCell::new(Some(Box::new(on_done))));
    let on_error: ErrorCb = Rc::new(RefCell::new(Some(Box::new(on_error))));

    // Pool of compute workers (filled after the scan completes).
    let pool: Rc<RefCell<Vec<Rc<Worker>>>> = Rc::new(RefCell::new(Vec::new()));
    // Per-worker partial outputs, keyed by worker index.
    let partials: Rc<RefCell<Vec<Option<crate::worker::WireL2Output>>>> =
        Rc::new(RefCell::new(Vec::new()));
    // Running progress tally (sum across workers).
    let snps_done_per_worker: Rc<RefCell<Vec<usize>>> = Rc::new(RefCell::new(Vec::new()));
    let snps_total: Rc<RefCell<usize>> = Rc::new(RefCell::new(0));
    // Latest wall_seconds reported by each worker (we report max).
    let wall_seconds_per_worker: Rc<RefCell<Vec<f64>>> = Rc::new(RefCell::new(Vec::new()));
    // Watchdog: per-worker last-activity timestamp (Date.now()).
    let last_activity_ms: Rc<RefCell<Vec<f64>>> = Rc::new(RefCell::new(Vec::new()));
    // Watchdog setInterval handle (so we can clear it on completion).
    let watchdog_handle: Rc<Cell<Option<i32>>> = Rc::new(Cell::new(None));

    // Helper: fire on_error and tear down everything (idempotent).
    // Both abort and the success path consult `completed` first so
    // exactly one of `on_done` / `on_error` fires per pool run.
    let abort_pool: AbortCb = {
        let pool = pool.clone();
        let on_error = on_error.clone();
        let on_done = on_done.clone();
        let completed = completed.clone();
        let watchdog_handle = watchdog_handle.clone();
        Rc::new(move |msg: String| {
            if completed.replace(true) {
                return;
            }
            clear_watchdog(&watchdog_handle);
            terminate_pool(&pool);
            // Drop the (unused) on_done so the FnOnce doesn't outlive
            // the pool indefinitely.
            let _ = on_done.borrow_mut().take();
            if let Some(cb) = on_error.borrow_mut().take() {
                cb(msg);
            }
        })
    };

    // Pre-clone the shared inputs that need to live across multiple
    // postMessage calls. SendWrapper is a no-op marker on wasm32
    // (single-threaded) — it just makes !Send web_sys types fit
    // through Leptos signal storage.
    let bed_file = Rc::new(send_wrapper::SendWrapper::new(bed_file));
    let bim_file = Rc::new(send_wrapper::SendWrapper::new(bim_file));
    let fam_text = Rc::new(fam_text);
    let config_json = Rc::new(config_json);

    // ── Step 1: spawn the scan worker ───────────────────────────────
    let scan_worker = Rc::new(new_module_worker()?);

    let scan_worker_for_handler = scan_worker.clone();
    let abort_for_scan = abort_pool.clone();
    let pool_for_scan = pool.clone();
    let partials_for_scan = partials.clone();
    let snps_done_for_scan = snps_done_per_worker.clone();
    let snps_total_for_scan = snps_total.clone();
    let wall_for_scan = wall_seconds_per_worker.clone();
    let activity_for_scan = last_activity_ms.clone();
    let watchdog_for_scan = watchdog_handle.clone();
    let on_progress_for_scan = on_progress.clone();
    let on_done_for_scan = on_done.clone();
    let on_error_for_scan = on_error.clone();
    let completed_for_scan = completed.clone();
    let abort_for_compute = abort_pool.clone();
    let bed_file_for_scan = bed_file.clone();
    let bim_file_for_scan = bim_file.clone();
    let fam_text_for_scan = fam_text.clone();
    let config_json_for_scan = config_json.clone();
    let wasm_js_url_for_scan = wasm_js_url_abs.clone();
    let wasm_bg_url_for_scan = wasm_bg_url_abs.clone();

    let scan_on_message: Closure<dyn FnMut(MessageEvent)> = Closure::new(
        move |ev: MessageEvent| {
            let data = ev.data();
            let kind = string_field(&data, "kind");
            match kind.as_str() {
                "ready" => { /* scan worker booted */ }
                "scan_bim_done" => {
                    let summaries_js = match Reflect::get(&data, &"summaries".into()) {
                        Ok(v) => v,
                        Err(e) => {
                            abort_for_scan(format!("scan response missing .summaries: {e:?}"));
                            return;
                        }
                    };
                    let summaries: Vec<crate::worker::WireChrSummary> =
                        match serde_wasm_bindgen::from_value(summaries_js) {
                            Ok(v) => v,
                            Err(e) => {
                                abort_for_scan(format!("scan summaries deserialise: {e}"));
                                return;
                            }
                        };
                    if summaries.is_empty() {
                        abort_for_scan("BIM contains no SNPs (scan returned empty summary)".into());
                        return;
                    }
                    let total: usize = summaries.iter().map(|s| s.snp_count).sum();
                    *snps_total_for_scan.borrow_mut() = total;
                    tracing::info!(
                        "scan_bim_done: {} chrs, {} total SNPs",
                        summaries.len(),
                        total
                    );

                    terminate_worker(&scan_worker_for_handler);

                    let n_workers_max = pool_cap().min(summaries.len()).max(1);
                    let pairs: Vec<(u8, usize)> =
                        summaries.iter().map(|s| (s.chr, s.snp_count)).collect();
                    let assignments = partition_chrs(&pairs, n_workers_max);
                    let n_workers = assignments.len();
                    let inner_threads = inner_threads_for_pool(n_workers);
                    tracing::info!(
                        "Spawning {n_workers} compute workers × {inner_threads} rayon threads; assignments = {assignments:?}"
                    );

                    *partials_for_scan.borrow_mut() = vec![None; n_workers];
                    *snps_done_for_scan.borrow_mut() = vec![0usize; n_workers];
                    *wall_for_scan.borrow_mut() = vec![0.0f64; n_workers];

                    // Initialize watchdog timestamps to now (so we don't
                    // immediately abort on a slow first chunk).
                    let now_ms = js_sys::Date::now();
                    *activity_for_scan.borrow_mut() = vec![now_ms; n_workers];

                    if let Err(e) = start_watchdog(
                        &watchdog_for_scan,
                        activity_for_scan.clone(),
                        abort_for_compute.clone(),
                    ) {
                        abort_for_compute(format!("watchdog setInterval: {e:?}"));
                        return;
                    }

                    for (worker_idx, chrs) in assignments.into_iter().enumerate() {
                        if let Err(e) = spawn_compute_worker(
                            worker_idx,
                            chrs,
                            inner_threads,
                            wasm_js_url_for_scan.clone(),
                            wasm_bg_url_for_scan.clone(),
                            bed_file_for_scan.clone(),
                            bim_file_for_scan.clone(),
                            fam_text_for_scan.clone(),
                            config_json_for_scan.clone(),
                            pool_for_scan.clone(),
                            partials_for_scan.clone(),
                            snps_done_for_scan.clone(),
                            snps_total_for_scan.clone(),
                            wall_for_scan.clone(),
                            activity_for_scan.clone(),
                            watchdog_for_scan.clone(),
                            on_progress_for_scan.clone(),
                            on_done_for_scan.clone(),
                            on_error_for_scan.clone(),
                            completed_for_scan.clone(),
                            abort_for_compute.clone(),
                        ) {
                            abort_for_compute(format!("spawn_compute_worker[{worker_idx}]: {e:?}"));
                            return;
                        }
                    }
                }
                "error" => {
                    let msg = string_field(&data, "error");
                    abort_for_scan(format!("scan worker: {msg}"));
                }
                other => {
                    tracing::warn!("scan worker: unknown kind '{other}'");
                }
            }
        },
    );
    scan_worker.set_onmessage(Some(scan_on_message.as_ref().unchecked_ref()));
    scan_on_message.forget();

    let abort_for_scan_err = abort_pool.clone();
    let scan_on_err: Closure<dyn FnMut(JsValue)> = Closure::new(move |ev: JsValue| {
        abort_for_scan_err(format!("scan worker.onerror: {ev:?}"));
    });
    scan_worker.set_onerror(Some(scan_on_err.as_ref().unchecked_ref()));
    scan_on_err.forget();

    // Init + scan messages to the scan worker. The scan worker just
    // streams BIM and returns the chr summary — no GEMM, so no rayon
    // pool needed (innerThreads = 1 → worker.js skips initThreadPool).
    post_init(&scan_worker, &wasm_js_url_abs, &wasm_bg_url_abs, 1)?;
    let scan_msg = js_sys::Object::new();
    Reflect::set(&scan_msg, &"kind".into(), &"scan_bim".into())?;
    Reflect::set(&scan_msg, &"bimFile".into(), file_jsvalue(&bim_file))?;
    scan_worker.post_message(&scan_msg)?;

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn spawn_compute_worker(
    worker_idx: usize,
    chrs: Vec<u8>,
    inner_threads: usize,
    wasm_js_url: JsValue,
    wasm_bg_url: JsValue,
    bed_file: Rc<send_wrapper::SendWrapper<web_sys::File>>,
    bim_file: Rc<send_wrapper::SendWrapper<web_sys::File>>,
    fam_text: Rc<String>,
    config_json: Rc<String>,
    pool: Rc<RefCell<Vec<Rc<Worker>>>>,
    partials: Rc<RefCell<Vec<Option<crate::worker::WireL2Output>>>>,
    snps_done_per_worker: Rc<RefCell<Vec<usize>>>,
    snps_total: Rc<RefCell<usize>>,
    wall_seconds_per_worker: Rc<RefCell<Vec<f64>>>,
    last_activity_ms: Rc<RefCell<Vec<f64>>>,
    watchdog_handle: Rc<Cell<Option<i32>>>,
    on_progress: ProgressCb,
    on_done: DoneCb,
    on_error: ErrorCb,
    completed: Rc<Cell<bool>>,
    abort_pool: AbortCb,
) -> Result<(), JsValue> {
    let worker = Rc::new(new_module_worker()?);
    pool.borrow_mut().push(worker.clone());

    let chrs_filter_json = serde_json::to_string(&chrs)
        .map_err(|e| JsValue::from_str(&format!("chrs serialise: {e}")))?;

    let worker_for_handler = worker.clone();
    let pool_for_handler = pool.clone();
    let partials_for_handler = partials.clone();
    let snps_done_for_handler = snps_done_per_worker.clone();
    let snps_total_for_handler = snps_total.clone();
    let wall_for_handler = wall_seconds_per_worker.clone();
    let activity_for_handler = last_activity_ms.clone();
    let watchdog_for_handler = watchdog_handle.clone();
    let on_progress_for_handler = on_progress.clone();
    let on_done_for_handler = on_done.clone();
    let on_error_for_handler = on_error.clone();
    let completed_for_handler = completed.clone();
    let abort_for_handler = abort_pool.clone();

    let on_msg: Closure<dyn FnMut(MessageEvent)> = Closure::new(move |ev: MessageEvent| {
        // The "completed" flag is also checked inside the success +
        // abort paths below, but a quick guard here saves us
        // unnecessary work on late messages from a worker we already
        // tore down.
        if completed_for_handler.get() {
            return;
        }
        let data = ev.data();
        let kind = string_field(&data, "kind");
        match kind.as_str() {
            "ready" => { /* compute worker booted */ }
            "progress" => {
                let snps_done = number_field(&data, "snps_done") as usize;
                snps_done_for_handler.borrow_mut()[worker_idx] = snps_done;
                activity_for_handler.borrow_mut()[worker_idx] = js_sys::Date::now();
                let total_done: usize = snps_done_for_handler.borrow().iter().sum();
                let total = *snps_total_for_handler.borrow();
                if let Ok(mut cb) = on_progress_for_handler.try_borrow_mut() {
                    cb(PoolProgress {
                        snps_done: total_done,
                        snps_total: total,
                    });
                }
            }
            "compute_l2_chrs_done" => {
                let result_js = match Reflect::get(&data, &"result".into()) {
                    Ok(v) => v,
                    Err(e) => {
                        abort_for_handler(format!("worker[{worker_idx}] missing .result: {e:?}"));
                        return;
                    }
                };
                let partial: crate::worker::WireL2Output =
                    match serde_wasm_bindgen::from_value(result_js) {
                        Ok(v) => v,
                        Err(e) => {
                            abort_for_handler(format!(
                                "worker[{worker_idx}] result deserialise: {e}"
                            ));
                            return;
                        }
                    };
                wall_for_handler.borrow_mut()[worker_idx] = partial.wall_seconds;
                snps_done_for_handler.borrow_mut()[worker_idx] = partial.n_snps;
                activity_for_handler.borrow_mut()[worker_idx] = js_sys::Date::now();
                // Per-worker [perf] breakdown via main-thread tracing.
                // Compute workers themselves can't safely use
                // tracing-wasm or console_error_panic_hook (they trap
                // at runtime in our multi-threaded WASM setup), so we
                // route per-phase totals back through the postMessage
                // result and log them here on main where the tracing
                // subscriber from `main()` is installed.
                let p = &partial.perf;
                let sketch_blurb = match p.sketch_secs {
                    Some(s) => format!(" sketch={s:.3}s"),
                    None => String::new(),
                };
                tracing::info!(
                    "[perf] worker[{worker_idx}] m={} n_indiv={} wall={:.3}s bed_read(stall)={:.3}s norm={:.3}s{sketch_blurb} bb_dot={:.3}s ab_dot={:.3}s ring_store={:.3}s",
                    p.m,
                    p.n_indiv,
                    partial.wall_seconds,
                    p.bed_read_secs,
                    p.norm_secs,
                    p.bb_dot_secs,
                    p.ab_dot_secs,
                    p.ring_store_secs,
                );
                partials_for_handler.borrow_mut()[worker_idx] = Some(partial);
                terminate_worker(&worker_for_handler);

                // All workers done? Aggregate + fire on_done.
                let all_done = partials_for_handler.borrow().iter().all(|p| p.is_some());
                if all_done {
                    if completed_for_handler.replace(true) {
                        return;
                    }
                    clear_watchdog(&watchdog_for_handler);
                    let assembled = assemble_partials(&partials_for_handler.borrow());
                    let per_worker_wall_seconds: Vec<f64> = wall_for_handler.borrow().clone();
                    let max_wall = per_worker_wall_seconds
                        .iter()
                        .copied()
                        .fold(0.0f64, f64::max);
                    // Snapshot each worker's perf in worker-index
                    // order. Workers that didn't report (shouldn't
                    // happen at this point — `all_done` checked
                    // above) get a default zero-filled WireL2Perf.
                    let per_worker_perf: Vec<crate::worker::WireL2Perf> = partials_for_handler
                        .borrow()
                        .iter()
                        .map(|opt| opt.as_ref().map(|o| o.perf.clone()).unwrap_or_default())
                        .collect();
                    let result = WorkerComputeResult {
                        bed_idx: assembled.0,
                        l2: assembled.1,
                        maf: assembled.2,
                        wall_seconds: max_wall,
                        n_snps: assembled.3,
                        per_worker_perf,
                        per_worker_wall_seconds,
                    };
                    // Defensive teardown of any stragglers (already
                    // terminated above) + drop on_error so the
                    // FnOnce doesn't outlive the pool.
                    terminate_pool(&pool_for_handler);
                    let _ = on_error_for_handler.borrow_mut().take();
                    if let Some(cb) = on_done_for_handler.borrow_mut().take() {
                        cb(result);
                    }
                }
            }
            "error" => {
                let msg = string_field(&data, "error");
                abort_for_handler(format!("worker[{worker_idx}]: {msg}"));
            }
            other => {
                tracing::warn!("worker[{worker_idx}]: unknown kind '{other}'");
            }
        }
    });
    worker.set_onmessage(Some(on_msg.as_ref().unchecked_ref()));
    on_msg.forget();

    let abort_for_err = abort_pool.clone();
    let on_err: Closure<dyn FnMut(JsValue)> = Closure::new(move |ev: JsValue| {
        abort_for_err(format!("worker[{worker_idx}].onerror: {ev:?}"));
    });
    worker.set_onerror(Some(on_err.as_ref().unchecked_ref()));
    on_err.forget();

    // Init + compute messages. Compute workers each spawn
    // `inner_threads` rayon workers (via wasm-bindgen-rayon's
    // `initThreadPool`) before processing their assigned chr-shard.
    post_init(&worker, &wasm_js_url, &wasm_bg_url, inner_threads)?;
    let compute_msg = js_sys::Object::new();
    Reflect::set(&compute_msg, &"kind".into(), &"compute_l2_chrs".into())?;
    Reflect::set(&compute_msg, &"bedFile".into(), file_jsvalue(&bed_file))?;
    Reflect::set(&compute_msg, &"bimFile".into(), file_jsvalue(&bim_file))?;
    Reflect::set(&compute_msg, &"famText".into(), &(*fam_text).clone().into())?;
    Reflect::set(
        &compute_msg,
        &"chrsFilterJson".into(),
        &chrs_filter_json.into(),
    )?;
    Reflect::set(
        &compute_msg,
        &"configJson".into(),
        &(*config_json).clone().into(),
    )?;
    worker.post_message(&compute_msg)?;
    Ok(())
}

/// Concatenate per-worker partial outputs into a global bed_idx /
/// L2 / MAF triple, sorted by `bed_idx` so all three follow BIM
/// order. The h² flow needs `bed_idx` to look up rsIDs in the BIM
/// for the sumstats join.
///
/// Length consistency `bed_idx.len() == l2.len() == maf.len()` is
/// enforced upstream by `worker_compute_l2_chrs` (`worker.rs`); the
/// `debug_assert_eq!` here is a redundant defense in depth so we
/// catch any future drift in dev builds even if the wasm-side
/// assertion is bypassed.
fn assemble_partials(
    partials: &[Option<crate::worker::WireL2Output>],
) -> (Vec<u32>, Vec<f64>, Vec<f64>, usize) {
    let n_total: usize = partials
        .iter()
        .filter_map(|p| p.as_ref().map(|o| o.n_snps))
        .sum();

    let mut tuples: Vec<(u32, f64, f64)> = Vec::with_capacity(n_total);
    for partial in partials.iter().flatten() {
        debug_assert_eq!(partial.bed_idx.len(), partial.l2.len());
        debug_assert_eq!(partial.l2.len(), partial.maf.len());
        for ((b, l), m) in partial
            .bed_idx
            .iter()
            .zip(partial.l2.iter())
            .zip(partial.maf.iter())
        {
            tuples.push((*b as u32, *l, *m));
        }
    }
    tuples.sort_by_key(|t| t.0);

    let mut bed_idx = Vec::with_capacity(n_total);
    let mut l2 = Vec::with_capacity(n_total);
    let mut maf = Vec::with_capacity(n_total);
    for (b, l, m) in tuples {
        bed_idx.push(b);
        l2.push(l);
        maf.push(m);
    }
    (bed_idx, l2, maf, n_total)
}

// ── Watchdog ────────────────────────────────────────────────────────

fn start_watchdog(
    watchdog_handle: &Rc<Cell<Option<i32>>>,
    last_activity_ms: Rc<RefCell<Vec<f64>>>,
    abort_pool: AbortCb,
) -> Result<(), JsValue> {
    let win = web_sys::window().ok_or_else(|| JsValue::from_str("no window"))?;
    let cb: Closure<dyn FnMut()> = Closure::new(move || {
        let now_ms = js_sys::Date::now();
        let stale: Vec<usize> = last_activity_ms
            .borrow()
            .iter()
            .enumerate()
            .filter_map(|(i, &t)| {
                if now_ms - t > WATCHDOG_TIMEOUT_MS {
                    Some(i)
                } else {
                    None
                }
            })
            .collect();
        if !stale.is_empty() {
            abort_pool(format!(
                "watchdog: workers {:?} silent for >{}s",
                stale,
                (WATCHDOG_TIMEOUT_MS / 1000.0) as u32
            ));
        }
    });
    let handle = win.set_interval_with_callback_and_timeout_and_arguments_0(
        cb.as_ref().unchecked_ref(),
        WATCHDOG_POLL_MS,
    )?;
    cb.forget();
    watchdog_handle.set(Some(handle));
    Ok(())
}

fn clear_watchdog(watchdog_handle: &Rc<Cell<Option<i32>>>) {
    if let Some(handle) = watchdog_handle.take()
        && let Some(win) = web_sys::window()
    {
        win.clear_interval_with_handle(handle);
    }
}

// ── Worker lifecycle helpers ────────────────────────────────────────

fn new_module_worker() -> Result<Worker, JsValue> {
    let opts = WorkerOptions::new();
    opts.set_type(WorkerType::Module);
    Worker::new_with_options("worker.js", &opts)
}

/// Cleanly tear down a single worker: detach the message + error
/// handlers (so JS releases its references to the leaked Closures),
/// then terminate. Idempotent; safe to call on an already-terminated
/// worker.
fn terminate_worker(worker: &Worker) {
    worker.set_onmessage(None);
    worker.set_onerror(None);
    worker.terminate();
}

/// Tear down every worker currently in the pool.
fn terminate_pool(pool: &Rc<RefCell<Vec<Rc<Worker>>>>) {
    for w in pool.borrow().iter() {
        terminate_worker(w);
    }
}

fn post_init(
    worker: &Worker,
    wasm_js_url: &JsValue,
    wasm_bg_url: &JsValue,
    inner_threads: usize,
) -> Result<(), JsValue> {
    let init = js_sys::Object::new();
    Reflect::set(&init, &"kind".into(), &"init".into())?;
    Reflect::set(&init, &"wasmJsUrl".into(), wasm_js_url)?;
    Reflect::set(&init, &"wasmBgUrl".into(), wasm_bg_url)?;
    // worker.js will skip pool init if `innerThreads <= 1` (the scan
    // worker, and any single-core fallback).
    Reflect::set(
        &init,
        &"innerThreads".into(),
        &JsValue::from_f64(inner_threads as f64),
    )?;
    worker.post_message(&init)?;
    Ok(())
}

fn absolutize(url: &str) -> Result<JsValue, JsValue> {
    let base = web_sys::window()
        .ok_or_else(|| JsValue::from_str("no window"))?
        .location()
        .href()?;
    Ok(JsValue::from_str(&Url::new_with_base(url, &base)?.href()))
}

/// `&File` derefs to `&JsValue` (web_sys deref chain through Blob /
/// Object), so `Reflect::set(.., file_jsvalue(&rc), ..)` works
/// without an explicit clone.
fn file_jsvalue(rc: &Rc<send_wrapper::SendWrapper<web_sys::File>>) -> &JsValue {
    rc.as_ref().as_ref()
}

/// Read a `.kind` (or other) string field from a postMessage payload,
/// returning an empty string if missing or non-string. Less verbose
/// than `Reflect::get + as_string + unwrap_or_default` at every site.
fn string_field(data: &JsValue, key: &str) -> String {
    Reflect::get(data, &JsValue::from_str(key))
        .ok()
        .and_then(|v| v.as_string())
        .unwrap_or_default()
}

/// Read a numeric field (postMessage `f64`) with `0` as the default.
fn number_field(data: &JsValue, key: &str) -> f64 {
    Reflect::get(data, &JsValue::from_str(key))
        .ok()
        .and_then(|v| v.as_f64())
        .unwrap_or(0.0)
}

// ── h² compute (Workstream L) ───────────────────────────────────────
//
// Spawns a single sibling Worker that reads the BIM (for SNP rsID
// lookup) + sumstats (for χ² / N per SNP), joins against the
// already-computed LD scores from the L2 panel, and runs LDSC's
// IRWLS regression with block jackknife SE. Returns the canonical
// 5-tuple stat geneticists quote in papers: h² (SE), intercept (SE),
// ratio, mean χ², λ_GC.
//
// Lifecycle:
//   1. Spawn one fresh Worker (no rayon pool; the h² regression is
//      sequential and tiny — biobank-scale wall is < 1 s).
//   2. Init wasm bundle with innerThreads=1.
//   3. Post `compute_h2` with the L2 array + bed_idx (from prior L2
//      run) + the BIM and sumstats File handles + jackknife config.
//   4. On `compute_h2_done`, hand back the H2Result and terminate.
//   5. On `error`, fire `on_error` with the message and terminate.
//
// Why one-shot vs reusing the L2 pool: each L2 worker is terminated
// after its compute completes, and the regression has different
// arg-shape needs (no BED, but yes sumstats). Cleaner to spin up a
// dedicated worker than retrofit the pool dispatch.
pub fn spawn_compute_h2(
    l2: Vec<f64>,
    maf: Vec<f64>,
    bed_idx: Vec<u32>,
    bim_file: web_sys::File,
    sumstats_file: web_sys::File,
    n_blocks: usize,
    two_step: Option<f64>,
    on_done: impl FnOnce(crate::worker::WireH2Result) + 'static,
    on_error: impl FnOnce(String) + 'static,
) -> Result<(), JsValue> {
    let (wasm_js_url, wasm_bg_url) = discover_wasm_urls()?;
    let wasm_js_url = absolutize(&wasm_js_url)?;
    let wasm_bg_url = absolutize(&wasm_bg_url)?;

    let worker = Rc::new(new_module_worker()?);
    let on_done: Rc<RefCell<Option<Box<dyn FnOnce(crate::worker::WireH2Result)>>>> =
        Rc::new(RefCell::new(Some(Box::new(on_done))));
    let on_error: Rc<RefCell<Option<Box<dyn FnOnce(String)>>>> =
        Rc::new(RefCell::new(Some(Box::new(on_error))));
    // Pack the L2 + MAF + bed_idx + file handles into a closure that
    // fires on `ready`. Doing it this way means the heavy payload
    // doesn't outlive the worker's first dispatch.
    let payload: Rc<
        RefCell<Option<(Vec<f64>, Vec<f64>, Vec<u32>, web_sys::File, web_sys::File)>>,
    > = Rc::new(RefCell::new(Some((l2, maf, bed_idx, bim_file, sumstats_file))));

    let worker_for_handler = worker.clone();
    let on_done_for_handler = on_done.clone();
    let on_error_for_handler = on_error.clone();
    let payload_for_handler = payload.clone();

    let on_msg: Closure<dyn FnMut(MessageEvent)> = Closure::new(move |ev: MessageEvent| {
        let data = ev.data();
        let kind = string_field(&data, "kind");
        match kind.as_str() {
            "ready" => {
                let Some((l2v, mafv, bedv, bim, sumstats)) =
                    payload_for_handler.borrow_mut().take()
                else {
                    return;
                };
                // Marshal L2 + MAF + bed_idx as typed arrays (proven
                // safe by the Workstream L probe sweep).
                let l2_ta = js_sys::Float64Array::new_with_length(l2v.len() as u32);
                l2_ta.copy_from(&l2v);
                let maf_ta = js_sys::Float64Array::new_with_length(mafv.len() as u32);
                maf_ta.copy_from(&mafv);
                let bed_ta = js_sys::Uint32Array::new_with_length(bedv.len() as u32);
                bed_ta.copy_from(&bedv);
                let msg = js_sys::Object::new();
                let _ = Reflect::set(&msg, &"kind".into(), &"compute_h2".into());
                let _ = Reflect::set(&msg, &"l2".into(), &l2_ta);
                let _ = Reflect::set(&msg, &"maf".into(), &maf_ta);
                let _ = Reflect::set(&msg, &"bedIdx".into(), &bed_ta);
                let _ = Reflect::set(&msg, &"bimFile".into(), bim.as_ref());
                let _ = Reflect::set(&msg, &"sumstatsFile".into(), sumstats.as_ref());
                let _ = Reflect::set(&msg, &"nBlocks".into(), &(n_blocks as f64).into());
                if let Some(t) = two_step {
                    let _ = Reflect::set(&msg, &"twoStep".into(), &t.into());
                }
                if let Err(e) = worker_for_handler.post_message(&msg) {
                    if let Some(cb) = on_error_for_handler.borrow_mut().take() {
                        cb(format!("compute_h2 postMessage: {e:?}"));
                    }
                    terminate_worker(&worker_for_handler);
                }
            }
            "compute_h2_done" => {
                let result_js = match Reflect::get(&data, &"result".into()) {
                    Ok(v) => v,
                    Err(e) => {
                        if let Some(cb) = on_error_for_handler.borrow_mut().take() {
                            cb(format!("compute_h2 missing .result: {e:?}"));
                        }
                        terminate_worker(&worker_for_handler);
                        return;
                    }
                };
                let parsed: crate::worker::WireH2Result =
                    match serde_wasm_bindgen::from_value(result_js) {
                        Ok(v) => v,
                        Err(e) => {
                            if let Some(cb) = on_error_for_handler.borrow_mut().take() {
                                cb(format!("compute_h2 result deserialise: {e}"));
                            }
                            terminate_worker(&worker_for_handler);
                            return;
                        }
                    };
                if let Some(cb) = on_done_for_handler.borrow_mut().take() {
                    cb(parsed);
                }
                terminate_worker(&worker_for_handler);
            }
            "error" => {
                let err = string_field(&data, "error");
                if let Some(cb) = on_error_for_handler.borrow_mut().take() {
                    cb(format!("compute_h2 worker error: {err}"));
                }
                terminate_worker(&worker_for_handler);
            }
            _ => {
                tracing::warn!("[compute_h2] unexpected message kind: {kind}");
            }
        }
    });
    worker.set_onmessage(Some(on_msg.as_ref().unchecked_ref()));
    on_msg.forget();

    let worker_for_err = worker.clone();
    let on_error_for_err = on_error.clone();
    let on_err: Closure<dyn FnMut(JsValue)> = Closure::new(move |ev: JsValue| {
        if let Some(cb) = on_error_for_err.borrow_mut().take() {
            cb(format!("compute_h2 worker.onerror: {ev:?}"));
        }
        terminate_worker(&worker_for_err);
    });
    worker.set_onerror(Some(on_err.as_ref().unchecked_ref()));
    on_err.forget();

    // h² is sequential — no rayon pool needed.
    post_init(&worker, &wasm_js_url, &wasm_bg_url, 1)?;
    Ok(())
}
