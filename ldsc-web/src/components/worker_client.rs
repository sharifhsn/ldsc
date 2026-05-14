//! Main-thread Web Worker client.
//!
//! Spawns the worker.js shim, hands it the wasm bundle URLs (which
//! are hash-suffixed by trunk on every build, so the main thread is
//! the natural source of truth), and orchestrates the request /
//! response message flow for `compute_l2`.
//!
//! Lifecycle is one-shot per call: spawn → init → compute → response
//! → terminate. Re-spawning per Run click is fine because worker
//! startup is sub-100 ms (just imports the wasm bundle the browser
//! already cached after the main thread initial load).

use std::cell::RefCell;
use std::rc::Rc;

use js_sys::Reflect;
use wasm_bindgen::JsCast;
use wasm_bindgen::prelude::*;
use web_sys::{Blob, BlobPropertyBag, MessageEvent, Url, Worker, WorkerOptions, WorkerType};

/// Result-of-compute payload sent from worker → main. Mirrors
/// [`crate::worker::WireL2Output`].
#[derive(Clone, Debug)]
pub struct WorkerComputeResult {
    pub l2: Vec<f64>,
    pub maf: Vec<f64>,
    pub wall_seconds: f64,
    pub n_snps: usize,
}

/// Find the trunk-emitted wasm bundle paths from the current
/// document. Trunk inserts `<link rel="modulepreload"
/// href="/ldsc/ldsc-web-{hash}.js">` and `<link rel="preload"
/// href="/ldsc/ldsc-web-{hash}_bg.wasm">` after the build, both
/// hash-suffixed; we look those up by suffix match.
pub fn discover_wasm_urls() -> Result<(String, String), JsValue> {
    let doc = web_sys::window().unwrap().document().unwrap();
    let head = doc.head().ok_or_else(|| JsValue::from_str("no <head>"))?;
    // `query_selector_all` lives on the ParentNode mixin, which
    // web-sys 0.3 doesn't expose directly on Document or
    // HtmlHeadElement. Walk children manually instead — there are
    // only ~10 link elements in the head.
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

/// Spawn the `worker.js` shim and run `worker_compute_l2` against
/// the supplied inputs. `on_done` fires on success, `on_error` on
/// any failure (worker can't load, wasm panic in worker, etc.).
///
/// Both callbacks are invoked exactly once. The worker is
/// terminated as soon as either fires.
pub fn spawn_compute_l2(
    bed_file: web_sys::File,
    bim_text: String,
    fam_text: String,
    config_json: String,
    on_done: impl FnOnce(WorkerComputeResult) + 'static,
    on_error: impl FnOnce(String) + 'static,
) -> Result<(), JsValue> {
    let (wasm_js_url, wasm_bg_url) = discover_wasm_urls()?;
    tracing::info!(
        "spawn_compute_l2: js={wasm_js_url} bg={wasm_bg_url} bed.size={}",
        bed_file.size()
    );

    let worker_opts = WorkerOptions::new();
    worker_opts.set_type(WorkerType::Module);
    let worker = Worker::new_with_options("worker.js", &worker_opts)?;
    let worker = Rc::new(worker);

    // Wrap the callbacks in `Option<...>` inside a `RefCell` so the
    // message handler can take them by value on the firing message
    // (FnOnce) without violating ownership.
    let on_done = Rc::new(RefCell::new(Some(on_done)));
    let on_error = Rc::new(RefCell::new(Some(on_error)));

    let worker_for_handler = worker.clone();
    let on_done_for_handler = on_done.clone();
    let on_error_for_handler = on_error.clone();
    let on_message: Closure<dyn FnMut(MessageEvent)> = Closure::new(move |ev: MessageEvent| {
        let data = ev.data();
        let kind = Reflect::get(&data, &JsValue::from_str("kind"))
            .ok()
            .and_then(|v| v.as_string())
            .unwrap_or_default();

        match kind.as_str() {
            "ready" => {
                tracing::info!("worker ready");
                // Init done — no-op; we already queued the compute
                // message immediately after init below.
            }
            "compute_l2_done" => {
                let result = match parse_compute_result(&data) {
                    Ok(r) => r,
                    Err(e) => {
                        if let Some(cb) = on_error_for_handler.borrow_mut().take() {
                            cb(format!("worker response parse: {e}"));
                        }
                        worker_for_handler.terminate();
                        return;
                    }
                };
                if let Some(cb) = on_done_for_handler.borrow_mut().take() {
                    cb(result);
                }
                worker_for_handler.terminate();
            }
            "error" => {
                let msg = Reflect::get(&data, &JsValue::from_str("error"))
                    .ok()
                    .and_then(|v| v.as_string())
                    .unwrap_or_else(|| "<no error message>".to_string());
                if let Some(cb) = on_error_for_handler.borrow_mut().take() {
                    cb(msg);
                }
                worker_for_handler.terminate();
            }
            _ => {
                tracing::warn!("worker: unknown message kind '{kind}'");
            }
        }
    });
    worker.set_onmessage(Some(on_message.as_ref().unchecked_ref()));
    // Leak the closure — the worker holds the only reference and
    // we want it alive for as long as the worker runs. The
    // worker.terminate() call inside the handler implicitly drops
    // the message channel.
    on_message.forget();

    let on_error_for_err = on_error.clone();
    let worker_for_err = worker.clone();
    let on_worker_error: Closure<dyn FnMut(JsValue)> = Closure::new(move |ev: JsValue| {
        let msg = format!("worker.onerror: {ev:?}");
        if let Some(cb) = on_error_for_err.borrow_mut().take() {
            cb(msg);
        }
        worker_for_err.terminate();
    });
    worker.set_onerror(Some(on_worker_error.as_ref().unchecked_ref()));
    on_worker_error.forget();

    // 1. Init message: tell worker where the wasm bundle lives.
    //    Use absolute URLs so the worker can `import()` them
    //    regardless of its own base path.
    let init_obj = js_sys::Object::new();
    Reflect::set(&init_obj, &"kind".into(), &"init".into())?;
    Reflect::set(&init_obj, &"wasmJsUrl".into(), &absolutize(&wasm_js_url)?)?;
    Reflect::set(&init_obj, &"wasmBgUrl".into(), &absolutize(&wasm_bg_url)?)?;
    worker.post_message(&init_obj)?;

    // 2. Compute message: hand over the BED File + parsed BIM/FAM
    //    text + config JSON.
    let compute_obj = js_sys::Object::new();
    Reflect::set(&compute_obj, &"kind".into(), &"compute_l2".into())?;
    Reflect::set(&compute_obj, &"bedFile".into(), &bed_file)?;
    Reflect::set(&compute_obj, &"bimText".into(), &bim_text.into())?;
    Reflect::set(&compute_obj, &"famText".into(), &fam_text.into())?;
    Reflect::set(&compute_obj, &"configJson".into(), &config_json.into())?;
    worker.post_message(&compute_obj)?;

    Ok(())
}

fn absolutize(url: &str) -> Result<JsValue, JsValue> {
    let base = web_sys::window().unwrap().location().href()?;
    Ok(JsValue::from_str(&Url::new_with_base(url, &base)?.href()))
}

fn parse_compute_result(data: &JsValue) -> Result<WorkerComputeResult, String> {
    let result = Reflect::get(data, &JsValue::from_str("result"))
        .map_err(|e| format!("missing .result: {e:?}"))?;
    let l2_js = Reflect::get(&result, &JsValue::from_str("l2"))
        .map_err(|e| format!("missing .l2: {e:?}"))?;
    let maf_js = Reflect::get(&result, &JsValue::from_str("maf"))
        .map_err(|e| format!("missing .maf: {e:?}"))?;
    let wall = Reflect::get(&result, &JsValue::from_str("wall_seconds"))
        .ok()
        .and_then(|v| v.as_f64())
        .ok_or_else(|| "wall_seconds not a number".to_string())?;
    let n_snps = Reflect::get(&result, &JsValue::from_str("n_snps"))
        .ok()
        .and_then(|v| v.as_f64())
        .map(|v| v as usize)
        .ok_or_else(|| "n_snps not a number".to_string())?;

    let l2 = js_sys::Array::from(&l2_js)
        .iter()
        .map(|v| v.as_f64().unwrap_or(f64::NAN))
        .collect();
    let maf = js_sys::Array::from(&maf_js)
        .iter()
        .map(|v| v.as_f64().unwrap_or(f64::NAN))
        .collect();

    Ok(WorkerComputeResult {
        l2,
        maf,
        wall_seconds: wall,
        n_snps,
    })
}

// Two unused imports we may need later; suppress warnings while the
// inline-blob worker option (commented out above) isn't in use.
#[allow(dead_code)]
fn _suppress_unused() {
    let _ = Blob::new_with_str_sequence;
    let _ = BlobPropertyBag::new;
}
