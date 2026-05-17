// ldsc-web Web Worker entry point.
//
// Loads the same wasm bundle the main thread uses, then dispatches
// inbound messages to the appropriate wasm export. Three active
// message kinds:
//
//   init         : { wasmJsUrl, wasmBgUrl, innerThreads }
//                  → load wasm bundle, init wasm-bindgen-rayon thread
//                    pool with `innerThreads` rayon workers (skipped
//                    when innerThreads ≤ 1), post { kind: 'ready' }.
//                    The pool MUST be initialised here, inside the
//                    Worker — initialising from main thread fails
//                    silently when crossbeam mutex contention later
//                    lowers to `Atomics.wait` (see plan workstream
//                    C.3 retrospective). Each Worker has its own wasm
//                    instance + memory + rayon pool; pools don't
//                    share memory across outer Workers.
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
            // First message: load the wasm bundle + init the rayon
            // thread pool. The main thread tells us where to find the
            // bundle (paths are hash-suffixed) and how many rayon
            // workers to spawn inside this outer Worker.
            //
            // If a previous init failed, `initPromise` is a rejected
            // Promise — resetting it on failure (in the catch below)
            // means the orchestrator can retry by sending another
            // init without spawning a new outer worker.
            const { wasmJsUrl, wasmBgUrl, innerThreads } = data;
            if (!initPromise) {
                initPromise = (async () => {
                    const mod = await import(wasmJsUrl);
                    await mod.default(wasmBgUrl);
                    wasm = mod;
                    // Manual SAB-backed parallel GEMM pool (workstream G):
                    // spawn `innerThreads - 1` nested Workers (this outer
                    // Worker counts as worker 0, hence the -1). Each
                    // nested Worker shares this outer Worker's wasm
                    // linear memory via SharedArrayBuffer (the
                    // `--shared-memory --import-memory` link-args make
                    // wasm.memory a SAB-backed WebAssembly.Memory; we
                    // pass it to each inner Worker's wasm init below).
                    //
                    // After all inner Workers report `ready`, we tell
                    // the wasm pool how many workers it has via
                    // `init_inner_pool(N)`; subsequent dispatches in
                    // `wasm_simd::pool::parallel_gemm_tn_f32` then
                    // partition output cols across N total workers.
                    if (typeof innerThreads === 'number' && innerThreads > 1) {
                        const innerCount = innerThreads - 1;
                        // Resolve inner_worker.js relative to OUR own
                        // location (same dir, served by trunk under
                        // /ldsc/).
                        const innerJsUrl = new URL('./inner_worker.js', self.location.href).href;
                        const innerWorkers = [];
                        const readyPromises = [];
                        for (let i = 0; i < innerCount; i++) {
                            const w = new Worker(innerJsUrl, { type: 'module' });
                            innerWorkers.push(w);
                            const workerId = i + 1; // outer is worker 0
                            readyPromises.push(new Promise((resolve, reject) => {
                                w.onmessage = (e) => {
                                    if (e.data?.kind === 'ready') {
                                        resolve();
                                    } else if (e.data?.kind === 'error') {
                                        reject(new Error(`inner_worker[${workerId}]: ${e.data.error}`));
                                    }
                                };
                                w.onerror = (ev) => {
                                    reject(new Error(`inner_worker[${workerId}].onerror: ${ev.message || ev}`));
                                };
                            }));
                            w.postMessage({
                                kind: 'init',
                                wasmJsUrl,
                                wasmBgUrl,
                                memory: mod.wasm_memory(),
                                workerId,
                            });
                        }
                        await Promise.all(readyPromises);
                        // Pool is now ready: tell the wasm side.
                        mod.init_inner_pool(innerThreads);
                        // Keep references so they aren't GC'd. (See
                        // wasm-bindgen-rayon's workerHelpers.js for
                        // the same pattern + Firefox bug rationale.)
                        self.__inner_workers = innerWorkers;
                    }
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
            // `worker_compute_l2_chrs` is now `async fn` on the Rust
            // side (workstream I — async pre-load of chr-shard bytes
            // via `await blob.slice(...).arrayBuffer()` per chr). The
            // wasm export therefore returns a Promise; await it.
            const result = await wasm.worker_compute_l2_chrs(
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
