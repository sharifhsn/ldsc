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
//!    biggest-first / smallest-residual to balance load.
//! 3. Spawn the compute workers in parallel. Each worker streams
//!    the BIM (filtered to its assigned chrs), opens the BED via
//!    `BlobBedSource`, runs `compute_l2_from_bed` with a progress
//!    callback that `postMessage`s back mid-compute.
//! 4. Aggregate per-worker results into one [`WorkerComputeResult`]
//!    with SNPs concatenated in BIM order (sorted by `bed_idx`).
//!    Total wall-time is `max` across workers (parallel),
//!    progress events fold into one global `(snps_done, snps_total)`
//!    counter.
//!
//! Lifecycle: workers are terminated as soon as their final
//! `compute_l2_chrs_done` message arrives. On any error path
//! (worker.onerror, parse failure, partial result), all sibling
//! workers are terminated and `on_error` fires.

use std::cell::RefCell;
use std::rc::Rc;

use js_sys::Reflect;
use wasm_bindgen::JsCast;
use wasm_bindgen::prelude::*;
use web_sys::{MessageEvent, Url, Worker, WorkerOptions, WorkerType};

/// Hard cap on pool size, regardless of `navigator.hardwareConcurrency`.
/// Each worker holds a chunk × N × precision genotype slab in wasm
/// linear memory (~40 MB at biobank N=50K, chunk=200, f32) — 4 workers
/// keeps that bounded around 160 MB renderer-side. Bump if profiling
/// shows the chr-shard distribution is the bottleneck.
const MAX_POOL: usize = 4;

// Shared type aliases for the once-cell callbacks the orchestrator
// passes around. Boxed so we can erase the concrete generic type
// (which would otherwise infect every Rc clone).
type DoneCb = Rc<RefCell<Option<Box<dyn FnOnce(WorkerComputeResult)>>>>;
type ErrorCb = Rc<RefCell<Option<Box<dyn FnOnce(String)>>>>;
type ProgressCb = Rc<RefCell<dyn FnMut(PoolProgress)>>;

/// Result-of-pool payload sent to `on_done`. Same shape as the
/// single-worker `WorkerComputeResult` was; SNPs are in BIM order.
#[derive(Clone, Debug)]
pub struct WorkerComputeResult {
    pub l2: Vec<f64>,
    pub maf: Vec<f64>,
    /// Wall time = max across workers (parallel; this is the
    /// user-perceived elapsed compute, not the sum of CPU work).
    pub wall_seconds: f64,
    pub n_snps: usize,
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
    let doc = web_sys::window().unwrap().document().unwrap();
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

/// `navigator.hardwareConcurrency`, capped at [`MAX_POOL`].
fn pool_cap() -> usize {
    let nav_hwc = web_sys::window()
        .map(|w| w.navigator())
        .and_then(|n| Reflect::get(&n, &"hardwareConcurrency".into()).ok())
        .and_then(|v| v.as_f64())
        .map(|v| v as usize)
        .unwrap_or(2);
    nav_hwc.clamp(1, MAX_POOL)
}

/// Partition `(chr, snp_count)` pairs across `n_workers` lists,
/// largest-first / smallest-residual. Returns one Vec<u8> of chrs
/// per worker; index in the outer Vec = worker index.
///
/// Greedy LPT (longest processing time first) gives load balance
/// within ~4/3 of optimal — fine for a 22-chr / 4-worker partition.
fn partition_chrs(summaries: &[(u8, usize)], n_workers: usize) -> Vec<Vec<u8>> {
    let n_workers = n_workers.max(1);
    let mut indexed: Vec<(u8, usize)> = summaries.to_vec();
    // Largest snp_count first.
    indexed.sort_by(|a, b| b.1.cmp(&a.1));

    let mut buckets: Vec<(usize, Vec<u8>)> = (0..n_workers).map(|_| (0, Vec::new())).collect();
    for (chr, n) in indexed {
        // Assign to the bucket with the smallest current load.
        let (idx, _) = buckets
            .iter()
            .enumerate()
            .min_by_key(|(_, (load, _))| *load)
            .unwrap();
        buckets[idx].0 += n;
        buckets[idx].1.push(chr);
    }
    buckets.into_iter().map(|(_, chrs)| chrs).collect()
}

/// Spawn a per-chr Web Worker pool to compute L2 against the supplied
/// inputs. Calls `on_progress` periodically with running totals,
/// `on_done` once with the assembled output, or `on_error` once on
/// any failure. Each callback is mutually exclusive: `on_done` and
/// `on_error` both fire at most once each.
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
    // `on_progress` runs many times (per chunk per worker); `on_done`
    // and `on_error` run once each. We Rc<RefCell<Option<_>>> the
    // FnOnce callbacks so they can be `take()`n on the firing message
    // without violating ownership.
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

    // Helper: tear down everything and fire `on_error` (idempotent).
    let abort_pool = {
        let pool = pool.clone();
        let on_error = on_error.clone();
        let on_done = on_done.clone();
        Rc::new(move |msg: String| {
            for w in pool.borrow().iter() {
                w.terminate();
            }
            // First error wins; subsequent errors are dropped.
            let _ = on_done.borrow_mut().take();
            if let Some(cb) = on_error.borrow_mut().take() {
                cb(msg);
            }
        })
    };

    // ── Step 1: spawn the scan worker ───────────────────────────────
    let scan_worker = new_module_worker()?;
    let scan_worker = Rc::new(scan_worker);

    // Closure factory: given a JsValue, produce a "compute" boot
    // sequence that fans out to N compute workers.
    let bed_file = Rc::new(send_wrapper::SendWrapper::new(bed_file));
    let bim_file = Rc::new(send_wrapper::SendWrapper::new(bim_file));
    let fam_text = Rc::new(fam_text);
    let config_json = Rc::new(config_json);

    let scan_worker_for_handler = scan_worker.clone();
    let abort_pool_for_scan = abort_pool.clone();
    let pool_for_scan = pool.clone();
    let partials_for_scan = partials.clone();
    let snps_done_for_scan = snps_done_per_worker.clone();
    let snps_total_for_scan = snps_total.clone();
    let wall_for_scan = wall_seconds_per_worker.clone();
    let on_progress_for_scan = on_progress.clone();
    let on_done_for_scan = on_done.clone();
    let abort_for_compute = abort_pool.clone();
    let bed_file_for_scan = bed_file.clone();
    let bim_file_for_scan = bim_file.clone();
    let fam_text_for_scan = fam_text.clone();
    let config_json_for_scan = config_json.clone();
    let wasm_js_url_for_scan = wasm_js_url_abs.clone();
    let wasm_bg_url_for_scan = wasm_bg_url_abs.clone();

    let scan_on_message: Closure<dyn FnMut(MessageEvent)> =
        Closure::new(move |ev: MessageEvent| {
            let data = ev.data();
            let kind = Reflect::get(&data, &"kind".into())
                .ok()
                .and_then(|v| v.as_string())
                .unwrap_or_default();
            match kind.as_str() {
                "ready" => { /* scan worker booted */ }
                "scan_bim_done" => {
                    let summaries_js = match Reflect::get(&data, &"summaries".into()) {
                        Ok(v) => v,
                        Err(e) => {
                            abort_pool_for_scan(format!("scan response missing .summaries: {e:?}"));
                            return;
                        }
                    };
                    let summaries: Vec<crate::worker::WireChrSummary> =
                        match serde_wasm_bindgen::from_value(summaries_js) {
                            Ok(v) => v,
                            Err(e) => {
                                abort_pool_for_scan(format!("scan summaries deserialise: {e}"));
                                return;
                            }
                        };
                    if summaries.is_empty() {
                        abort_pool_for_scan(
                            "BIM contains no SNPs (scan returned empty summary)".into(),
                        );
                        return;
                    }
                    let total: usize = summaries.iter().map(|s| s.snp_count).sum();
                    *snps_total_for_scan.borrow_mut() = total;
                    tracing::info!(
                        "scan_bim_done: {} chrs, {} total SNPs",
                        summaries.len(),
                        total
                    );

                    // Done with scan worker.
                    scan_worker_for_handler.terminate();

                    // Plan the partition.
                    let n_workers = pool_cap().min(summaries.len()).max(1);
                    let pairs: Vec<(u8, usize)> =
                        summaries.iter().map(|s| (s.chr, s.snp_count)).collect();
                    let assignments = partition_chrs(&pairs, n_workers);
                    let active_workers = assignments.iter().filter(|a| !a.is_empty()).count();
                    tracing::info!(
                        "Spawning {} compute workers; assignments = {:?}",
                        active_workers,
                        assignments
                    );

                    *partials_for_scan.borrow_mut() = vec![None; assignments.len()];
                    *snps_done_for_scan.borrow_mut() = vec![0usize; assignments.len()];
                    *wall_for_scan.borrow_mut() = vec![0.0f64; assignments.len()];

                    // Spawn each compute worker.
                    for (worker_idx, chrs) in assignments.into_iter().enumerate() {
                        if chrs.is_empty() {
                            // Nothing assigned; record an empty partial.
                            partials_for_scan.borrow_mut()[worker_idx] =
                                Some(crate::worker::WireL2Output {
                                    bed_idx: Vec::new(),
                                    l2: Vec::new(),
                                    maf: Vec::new(),
                                    wall_seconds: 0.0,
                                    n_snps: 0,
                                });
                            continue;
                        }
                        if let Err(e) = spawn_compute_worker(
                            worker_idx,
                            chrs,
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
                            on_progress_for_scan.clone(),
                            on_done_for_scan.clone(),
                            abort_for_compute.clone(),
                        ) {
                            abort_for_compute(format!("spawn_compute_worker[{worker_idx}]: {e:?}"));
                            return;
                        }
                    }
                }
                "error" => {
                    let msg = Reflect::get(&data, &"error".into())
                        .ok()
                        .and_then(|v| v.as_string())
                        .unwrap_or_else(|| "<no error message>".to_string());
                    abort_pool_for_scan(format!("scan worker: {msg}"));
                }
                other => {
                    tracing::warn!("scan worker: unknown kind '{other}'");
                }
            }
        });
    scan_worker.set_onmessage(Some(scan_on_message.as_ref().unchecked_ref()));
    scan_on_message.forget();

    let abort_for_scan_err = abort_pool.clone();
    let scan_on_err: Closure<dyn FnMut(JsValue)> = Closure::new(move |ev: JsValue| {
        abort_for_scan_err(format!("scan worker.onerror: {ev:?}"));
    });
    scan_worker.set_onerror(Some(scan_on_err.as_ref().unchecked_ref()));
    scan_on_err.forget();

    // Init + scan messages to the scan worker.
    post_init(&scan_worker, &wasm_js_url_abs, &wasm_bg_url_abs)?;
    let scan_msg = js_sys::Object::new();
    Reflect::set(&scan_msg, &"kind".into(), &"scan_bim".into())?;
    Reflect::set(&scan_msg, &"bimFile".into(), bim_file.as_ref().as_ref())?;
    scan_worker.post_message(&scan_msg)?;

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn spawn_compute_worker(
    worker_idx: usize,
    chrs: Vec<u8>,
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
    on_progress: ProgressCb,
    on_done: DoneCb,
    abort_pool: Rc<dyn Fn(String)>,
) -> Result<(), JsValue> {
    let worker = Rc::new(new_module_worker()?);
    pool.borrow_mut().push(worker.clone());

    let chrs_filter_json = serde_json::to_string(&chrs)
        .map_err(|e| JsValue::from_str(&format!("chrs serialise: {e}")))?;

    // Per-worker on_message handler.
    let worker_for_handler = worker.clone();
    let pool_for_handler = pool.clone();
    let partials_for_handler = partials.clone();
    let snps_done_for_handler = snps_done_per_worker.clone();
    let snps_total_for_handler = snps_total.clone();
    let wall_for_handler = wall_seconds_per_worker.clone();
    let on_progress_for_handler = on_progress.clone();
    let on_done_for_handler = on_done.clone();
    let abort_for_handler = abort_pool.clone();

    let on_msg: Closure<dyn FnMut(MessageEvent)> = Closure::new(move |ev: MessageEvent| {
        let data = ev.data();
        let kind = Reflect::get(&data, &"kind".into())
            .ok()
            .and_then(|v| v.as_string())
            .unwrap_or_default();
        match kind.as_str() {
            "ready" => { /* compute worker booted */ }
            "progress" => {
                let snps_done = Reflect::get(&data, &"snps_done".into())
                    .ok()
                    .and_then(|v| v.as_f64())
                    .map(|v| v as usize)
                    .unwrap_or(0);
                snps_done_for_handler.borrow_mut()[worker_idx] = snps_done;
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
                partials_for_handler.borrow_mut()[worker_idx] = Some(partial);
                worker_for_handler.terminate();

                // All workers done? Aggregate + fire on_done.
                let all_done = partials_for_handler.borrow().iter().all(|p| p.is_some());
                if all_done {
                    let assembled = assemble_partials(&partials_for_handler.borrow());
                    let max_wall = wall_for_handler
                        .borrow()
                        .iter()
                        .copied()
                        .fold(0.0f64, f64::max);
                    let result = WorkerComputeResult {
                        l2: assembled.0,
                        maf: assembled.1,
                        wall_seconds: max_wall,
                        n_snps: assembled.2,
                    };
                    // Tear down any stragglers (defensive — they're
                    // already terminated above).
                    for w in pool_for_handler.borrow().iter() {
                        w.terminate();
                    }
                    if let Some(cb) = on_done_for_handler.borrow_mut().take() {
                        cb(result);
                    }
                }
            }
            "error" => {
                let msg = Reflect::get(&data, &"error".into())
                    .ok()
                    .and_then(|v| v.as_string())
                    .unwrap_or_else(|| "<no error message>".to_string());
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

    // Init + compute messages.
    post_init(&worker, &wasm_js_url, &wasm_bg_url)?;
    let compute_msg = js_sys::Object::new();
    Reflect::set(&compute_msg, &"kind".into(), &"compute_l2_chrs".into())?;
    Reflect::set(&compute_msg, &"bedFile".into(), bed_file.as_ref().as_ref())?;
    Reflect::set(&compute_msg, &"bimFile".into(), bim_file.as_ref().as_ref())?;
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

/// Concatenate per-worker partial outputs into a global L2 / MAF
/// pair, sorted by `bed_idx` so the result follows BIM order.
fn assemble_partials(
    partials: &[Option<crate::worker::WireL2Output>],
) -> (Vec<f64>, Vec<f64>, usize) {
    // Total SNP count.
    let n_total: usize = partials
        .iter()
        .filter_map(|p| p.as_ref().map(|o| o.n_snps))
        .sum();

    // Gather all (bed_idx, l2, maf) tuples, sort.
    let mut tuples: Vec<(usize, f64, f64)> = Vec::with_capacity(n_total);
    for partial in partials.iter().flatten() {
        for ((b, l), m) in partial
            .bed_idx
            .iter()
            .zip(partial.l2.iter())
            .zip(partial.maf.iter())
        {
            tuples.push((*b, *l, *m));
        }
    }
    tuples.sort_by_key(|t| t.0);

    let mut l2 = Vec::with_capacity(n_total);
    let mut maf = Vec::with_capacity(n_total);
    for (_, l, m) in tuples {
        l2.push(l);
        maf.push(m);
    }
    (l2, maf, n_total)
}

fn new_module_worker() -> Result<Worker, JsValue> {
    let opts = WorkerOptions::new();
    opts.set_type(WorkerType::Module);
    Worker::new_with_options("worker.js", &opts)
}

fn post_init(worker: &Worker, wasm_js_url: &JsValue, wasm_bg_url: &JsValue) -> Result<(), JsValue> {
    let init = js_sys::Object::new();
    Reflect::set(&init, &"kind".into(), &"init".into())?;
    Reflect::set(&init, &"wasmJsUrl".into(), wasm_js_url)?;
    Reflect::set(&init, &"wasmBgUrl".into(), wasm_bg_url)?;
    worker.post_message(&init)?;
    Ok(())
}

fn absolutize(url: &str) -> Result<JsValue, JsValue> {
    let base = web_sys::window().unwrap().location().href()?;
    Ok(JsValue::from_str(&Url::new_with_base(url, &base)?.href()))
}
