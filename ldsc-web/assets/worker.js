// ldsc-web Web Worker entry point.
//
// Loads the same wasm bundle the main thread uses, then waits for an
// init message that hands us the (hash-suffixed) bundle URL plus the
// computation payload. Why route the URL through postMessage rather
// than hard-coding it: trunk content-hashes the JS+wasm filenames on
// every build, so a static URL would 404 after every rebuild. The
// main thread already knows the active hashes (it loaded itself from
// them) and is the natural source of truth.
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

        if (kind === 'compute_l2') {
            // Make sure init is done. The main thread always sends
            // an init message before compute_l2, but we await
            // defensively in case messages get reordered.
            if (initPromise) await initPromise;
            if (!wasm) {
                self.postMessage({ kind: 'error', error: 'wasm not initialised' });
                return;
            }

            const { bedFile, bimText, famText, configJson } = e.data;
            const t0 = performance.now();
            const result = wasm.worker_compute_l2(bedFile, bimText, famText, configJson);
            const dt = performance.now() - t0;
            self.postMessage({ kind: 'compute_l2_done', result, walkSecondsJs: dt / 1000 });
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
