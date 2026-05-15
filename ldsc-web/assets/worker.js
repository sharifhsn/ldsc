// ldsc-web Web Worker entry point.
//
// Loads the same wasm bundle the main thread uses, then dispatches
// inbound messages to the appropriate wasm export. There are three
// active message kinds:
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
//                    { kind: 'progress', chunks_done, chunks_total, snps_done, snps_total, chunk_wall_ms }
//                    directly via DedicatedWorkerGlobalScope.postMessage.
//                    We don't intercept those — they flow straight through.
//
// We use ES module workers (`new Worker(url, { type: 'module' })`)
// so we can `import()` the trunk-emitted ES module bundle directly.
// All modern browsers support module workers (Chrome 80+, Firefox
// 114+, Safari 15+).

let wasm = null;
let initPromise = null;

self.onmessage = async (e) => {
    try {
        const { kind } = e.data || {};

        if (kind === 'init') {
            // First message: load the wasm bundle. The main thread
            // tells us where to find it (paths are hash-suffixed).
            const { wasmJsUrl, wasmBgUrl } = e.data;
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
            const { bimFile } = e.data;
            const t0 = performance.now();
            const summaries = wasm.worker_scan_bim(bimFile);
            const dt = performance.now() - t0;
            self.postMessage({ kind: 'scan_bim_done', summaries, walkSecondsJs: dt / 1000 });
            return;
        }

        if (kind === 'compute_l2_chrs') {
            const { bedFile, bimFile, famText, chrsFilterJson, configJson } = e.data;
            const t0 = performance.now();
            const result = wasm.worker_compute_l2_chrs(
                bedFile, bimFile, famText, chrsFilterJson, configJson,
            );
            const dt = performance.now() - t0;
            self.postMessage({ kind: 'compute_l2_chrs_done', result, walkSecondsJs: dt / 1000 });
            return;
        }

        self.postMessage({ kind: 'error', error: `unknown message kind: ${kind}` });
    } catch (err) {
        // Forward any exception (including wasm-side panics that
        // propagated out as JsError → JS exception) back to the main
        // thread so the UI can render it instead of dying silently.
        self.postMessage({
            kind: 'error',
            error: (err && err.stack) ? err.stack : String(err),
        });
    }
};
