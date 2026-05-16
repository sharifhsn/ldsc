// Inner Web Worker for the manual SAB-backed parallel GEMM pool
// (workstream G). One nested Worker spawned per rayon thread inside
// each outer compute Worker. Shares the outer Worker's wasm linear
// memory via SharedArrayBuffer (the `--shared-memory --import-memory`
// link-args make the wasm declare its memory as `(import ...
// shared)`, and the JS side passes the outer Worker's
// `WebAssembly.Memory` to `mod.default(...)` here).
//
// Lifecycle:
//
//   1. Outer worker spawns this file via `new Worker('inner_worker.js',
//      { type: 'module' })` and posts a single 'init' message with:
//
//        { wasmJsUrl, wasmBgUrl, memory, workerId }
//
//      where `memory` is the outer Worker's `WebAssembly.Memory`
//      object (shared SAB-backed) and `workerId` is the inner
//      Worker's index in `1..N` (worker 0 is the outer Worker).
//
//   2. We import the wasm bundle, instantiate it with the SHARED
//      memory (so atomics in this Worker hit the same SAB the outer
//      Worker writes), post a 'ready' acknowledgement, then call
//      `wasm.inner_worker_loop(workerId)` which parks forever in the
//      shared-atomics dispatch loop in `ldsc::wasm_simd::pool`.
//
// We use ES module workers (`type: 'module'`) so we can dynamic-import
// the trunk-emitted JS bundle directly. All modern browsers support
// module workers (Chrome 80+, Firefox 114+, Safari 15+).

self.onmessage = async (e) => {
    const data = e.data;
    if (!data || data.kind !== 'init') {
        self.postMessage({
            kind: 'error',
            error: `inner_worker.js: expected {kind:'init', ...}, got ${JSON.stringify(data)}`,
        });
        return;
    }
    try {
        const { wasmJsUrl, wasmBgUrl, memory, workerId } = data;
        const mod = await import(wasmJsUrl);
        // wasm-bindgen's __wbg_init accepts an init descriptor where
        // we can override the `memory` import — that's how we bind
        // this nested Worker's wasm instance to the outer Worker's
        // SAB-backed memory.
        //
        // **CRITICAL:** pass `thread_stack_size` so that
        // `__wbindgen_start(stack_size)` allocates a UNIQUE stack +
        // thread-local-storage region for this inner Worker in the
        // SHARED linear memory. Without this, every inner Worker
        // shares the outer Worker's TLS region — concurrent
        // thread_local!() access from multiple inner Workers
        // silently corrupts memory and crashes on large workloads
        // (e.g. full 1000G; small workloads like chr22 don't trigger
        // it because the rare TLS access doesn't race). Must be a
        // multiple of the wasm page size (65 536). 2 MB matches
        // wasm-bindgen-rayon's WASM_BINDGEN_THREADS_STACK_SIZE
        // default and is more than enough for our compute path
        // (deepest call stack is ~10 frames in the GEMM kernel).
        const THREAD_STACK_SIZE = 2 * 1024 * 1024; // 2 MB
        await mod.default({
            module_or_path: wasmBgUrl,
            memory,
            thread_stack_size: THREAD_STACK_SIZE,
        });
        self.postMessage({ kind: 'ready', workerId });
        // Park forever (returns only on exit signal). All work happens
        // via shared atomics — no further postMessages between outer
        // and inner workers during compute.
        mod.inner_worker_loop(workerId);
    } catch (err) {
        self.postMessage({
            kind: 'error',
            error: (err && err.stack) ? err.stack : String(err),
        });
    }
};
