// ldsc-web Web Worker entry point.
//
// Loads the same wasm bundle the main thread uses, then dispatches
// inbound messages to the appropriate wasm export. Three active
// message kinds:
//
//   init         : { wasmJsUrl, wasmBgUrl }
//                  → load wasm bundle, post { kind: 'ready' }
//
//   scan_bim     : { bimFile }
//                  → call wasm.worker_scan_bim, post
//                    { kind: 'scan_bim_done', summaries: [{chr, snp_count}] }
//
//   compute_l2_chrs : { bedFile, bimFile, famText, chrsFilterJson, configJson }
//                  → call wasm.worker_compute_l2_chrs, post
//                    { kind: 'compute_l2_chrs_done', result: WireL2Output }
//                    Mid-compute, the wasm callback posts
//                    { kind: 'progress', snps_done } directly via
//                    DedicatedWorkerGlobalScope.postMessage. We
//                    don't intercept those — they flow straight through.
//
// Unknown / typo'd `kind` values short-circuit to a structured
// error rather than silent drop.
//
// We use ES module workers (`new Worker(url, { type: 'module' })`)
// so we can `import()` the trunk-emitted ES module bundle directly.
// All modern browsers support module workers (Chrome 80+, Firefox
// 114+, Safari 15+).

// Source-of-truth set of accepted message kinds. Keep in sync with
// the Rust dispatch table in ldsc-web/src/components/worker_client.rs.
const KNOWN_KINDS = new Set(['init', 'scan_bim', 'compute_l2_chrs']);

let wasm = null;
let initPromise = null;

self.onmessage = async (e) => {
    const data = e.data;
    if (!data || typeof data !== 'object') {
        self.postMessage({
            kind: 'error',
            error: `worker.js: expected an object payload, got ${typeof data}`,
        });
        return;
    }
    const { kind } = data;
    if (!KNOWN_KINDS.has(kind)) {
        self.postMessage({
            kind: 'error',
            error: `worker.js: unknown message kind ${JSON.stringify(kind)}`,
        });
        return;
    }

    try {
        if (kind === 'init') {
            // First message: load the wasm bundle. The main thread
            // tells us where to find it (paths are hash-suffixed).
            // If a previous init failed, `initPromise` is a rejected
            // Promise — resetting it on failure (in the catch below)
            // means the orchestrator can retry by sending another
            // init without spawning a new worker.
            const { wasmJsUrl, wasmBgUrl } = data;
            if (!initPromise) {
                initPromise = (async () => {
                    const mod = await import(wasmJsUrl);
                    await mod.default(wasmBgUrl);
                    wasm = mod;
                    self.postMessage({ kind: 'ready' });
                })();
            }
            await initPromise;
            return;
        }

        // Make sure init is done. The main thread always sends an
        // init message before any compute/scan, but we await
        // defensively in case messages get reordered.
        if (initPromise) await initPromise;
        if (!wasm) {
            self.postMessage({ kind: 'error', error: 'wasm not initialised' });
            return;
        }

        if (kind === 'scan_bim') {
            const { bimFile } = data;
            const t0 = performance.now();
            const summaries = wasm.worker_scan_bim(bimFile);
            const dt = performance.now() - t0;
            self.postMessage({ kind: 'scan_bim_done', summaries, wallSecondsJs: dt / 1000 });
            return;
        }

        if (kind === 'compute_l2_chrs') {
            const { bedFile, bimFile, famText, chrsFilterJson, configJson } = data;
            const t0 = performance.now();
            const result = wasm.worker_compute_l2_chrs(
                bedFile, bimFile, famText, chrsFilterJson, configJson,
            );
            const dt = performance.now() - t0;
            self.postMessage({ kind: 'compute_l2_chrs_done', result, wallSecondsJs: dt / 1000 });
            return;
        }

        // Unreachable — KNOWN_KINDS guard above ensures we never
        // get here, but keep a defensive postMessage in case the
        // set drifts out of sync with the dispatch table.
        self.postMessage({ kind: 'error', error: `worker.js: unhandled kind '${kind}'` });
    } catch (err) {
        // If init failed, reset `initPromise` so the next init
        // message can retry. Otherwise (compute/scan failure), the
        // wasm bundle is fine and we just forward the exception.
        if (kind === 'init') {
            initPromise = null;
            wasm = null;
        }
        // Forward any exception (including wasm-side panics that
        // propagated out as JsError → JS exception) back to the main
        // thread so the UI can render it instead of dying silently.
        self.postMessage({
            kind: 'error',
            error: (err && err.stack) ? err.stack : String(err),
        });
    }
};
