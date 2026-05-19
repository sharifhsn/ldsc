//! The L2 panel — the headline demo surface.
//!
//! Layout (LDLink card idiom):
//! - Upload card: three file rows (.bed / .bim / .fam) with loaded badges.
//! - Configuration card: mode preset + headline knobs + advanced.
//! - Run button.
//! - Results card: timing + summary stats + plotters charts +
//!   downloadable .l2.ldscore.gz.
//!
//! All compute happens by calling `ldsc::l2::compute_l2_from_bytes`
//! synchronously on the main thread. For the MVP that's fine on
//! 1000G-sized inputs (single-second runs); biobank-scale will need
//! the wasm-bindgen-rayon Web Worker pool from Workstream C.3.

use leptos::ev::MouseEvent;
use leptos::prelude::*;

use super::charts::render_l2_vs_maf;
use super::file_upload::{FileUploadRow, LoadedFile, ReadAs};
use super::worker_client::{
    PoolProgress, WorkerComputeResult, spawn_compute_h2, spawn_compute_l2_pool,
};
use crate::worker::{WireH2Result, WireL2Config, WireWindowMode};

/// Mode preset radio (mirrors the "Mode" section of the CLI README's
/// recommendations table). Selecting one auto-toggles the underlying
/// flags in `effective_config`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ModePreset {
    /// Default-recommended: CountSketch d=1000 (truth-cluster: h² within
    /// ~0.003 of per-SNP exact across N=503 to N=100K, Pearson r ≥ 0.993
    /// universally). See preprint §"Optimal d as a function of N".
    Sketch,
    /// Chunked exact (Rust default; chunk_size=200 — slightly inflates
    /// L2 vs Python by mean ~0.2 / max ~6).
    Default,
    /// Bit-identical to Python LDSC. chunk_size=50, global_pass.
    PythonCompat,
    /// LDSC-paper-canonical per-SNP exact windows. Combine with Sketch
    /// for the headline `--sketch 1000 --snp-level-masking` mode.
    PerSnpExact,
}

impl ModePreset {
    fn label(self) -> &'static str {
        match self {
            ModePreset::Sketch => "Sketch (recommended)",
            ModePreset::Default => "Default (chunked exact)",
            ModePreset::PythonCompat => "Python-compat (bit-identical to Python LDSC)",
            ModePreset::PerSnpExact => "Per-SNP exact (LDSC paper definition)",
        }
    }
    fn description(self) -> &'static str {
        match self {
            ModePreset::Sketch => {
                "CountSketch projection on the individual axis. Fast and accurate; \
                 d=1000 + per-SNP masking puts h² within ~0.003 of per-SNP exact \
                 across all panel sizes (1000G → biobank, N=503 to N=100K). \
                 d=5000 for r ≥ 0.999 if per-SNP accuracy matters."
            }
            ModePreset::Default => {
                "Rust default. Chunk size 200 trades a tiny accuracy hit for ~4× \
                 faster GEMM than Python LDSC's chunk size 50."
            }
            ModePreset::PythonCompat => {
                "Forces chunk size 50 + global pass. Output bit-identical to Python LDSC."
            }
            ModePreset::PerSnpExact => {
                "Per-SNP exact windows after r² masking. Matches the LDSC paper's \
                 mathematical definition. Combine with Sketch for the headline mode."
            }
        }
    }
}

#[derive(Default, Clone)]
struct RunStats {
    n_snps: usize,
    mean_l2: f64,
    maf_l2_corr: f64,
    wall_seconds: f64,
}

impl RunStats {
    fn from_output(out: &WorkerComputeResult) -> Self {
        Self {
            n_snps: out.l2.len(),
            mean_l2: mean(&out.l2),
            maf_l2_corr: pearson(&out.maf, &out.l2),
            wall_seconds: out.wall_seconds,
        }
    }
}

#[component]
pub fn L2Panel() -> impl IntoView {
    // ── Upload state (one signal per file) ─────────────────────────────
    let bed = RwSignal::new(LoadedFile::default());
    let bim = RwSignal::new(LoadedFile::default());
    let fam = RwSignal::new(LoadedFile::default());
    // .sumstats is optional — only needed if the user wants to run
    // the h² regression after L2 finishes. Format = LDSC's munged
    // sumstats text (header + tab-sep SNP/Z/N rows). Uploaded as a
    // raw File handle so the h² worker can `await blob.text()`
    // without main-thread copies.
    let sumstats = RwSignal::new(LoadedFile::default());

    // ── Config state ───────────────────────────────────────────────────
    let preset = RwSignal::new(ModePreset::Sketch);
    // d=1000 is the empirically-validated universal sweet spot from
    // the N×d sweep in `preprint/data/dn_sweep_full.csv`: across
    // N ∈ {503, 10K, 20K, 50K, 100K}, d=1000 with --snp-level-masking
    // gives Pearson r ≥ 0.993 vs exact and h² within ~0.003 of
    // per-SNP exact at all sizes, at <15s wall even at N=100K.
    let sketch_d = RwSignal::new(1000u32);
    let snp_level_masking = RwSignal::new(true); // canonical per the README headline
    // f32 is on by default — ~1.85× faster than f64 with negligible
    // accuracy hit at chr22 demo scale. Sketch modes implicitly use
    // f32 anyway (CountSketch ±1 is exactly representable), so this
    // toggle is only meaningful for the chunked-exact / Python-compat
    // / per-SNP-exact paths.
    let fast_f32 = RwSignal::new(true);
    let ld_wind_kb = RwSignal::new(1000.0f64);

    // ── Run state ──────────────────────────────────────────────────────
    let is_running = RwSignal::new(false);
    let result: RwSignal<Option<WorkerComputeResult>> = RwSignal::new(None);
    let error: RwSignal<Option<String>> = RwSignal::new(None);
    // Live progress: `(snps_done, snps_total)`. `(0, 0)` = unknown
    // (scan worker hasn't reported yet). Updated mid-compute by
    // progress messages from each per-chr worker.
    let progress: RwSignal<(usize, usize)> = RwSignal::new((0, 0));

    // Per-bed_idx mapping from L2 array position back to BIM row.
    // Held alongside `result` so the h² worker can pass it through
    // unchanged. Reset to `None` whenever a new L2 run kicks off.
    let bed_idx_for_h2: RwSignal<Option<Vec<u32>>> = RwSignal::new(None);

    // h² state — separate from L2's `is_running` because they run
    // sequentially but as independent worker jobs.
    let h2_running = RwSignal::new(false);
    let h2_result: RwSignal<Option<WireH2Result>> = RwSignal::new(None);
    let h2_error: RwSignal<Option<String>> = RwSignal::new(None);

    // Re-render the MAF/L2 scatter whenever a new L2 result lands.
    // The histogram was dropped as part of the demo refocus
    // (Workstream K — geneticists never publish L2 histograms; the
    // percentiles in the stat grid already convey the distribution
    // shape, and the MAF/L2 scatter is the one plot with documented
    // diagnostic precedent in the LDSC ecosystem).
    Effect::new(move |_| {
        if let Some(out) = result.get() {
            // Defer one tick so the canvas element exists in the DOM.
            // Guard against detachment (user switched modules / route
            // mid-render); the canvas may be gone by the time the
            // closure runs.
            request_animation_frame(move || {
                let doc = match web_sys::window().and_then(|w| w.document()) {
                    Some(d) => d,
                    None => return,
                };
                if doc.get_element_by_id("chart-scatter").is_some()
                    && let Err(e) = render_l2_vs_maf("chart-scatter", &out.maf, &out.l2)
                {
                    tracing::error!("scatter render: {e:?}");
                }
            });
        }
    });

    let all_loaded = move || {
        bed.with(|f| f.file_handle.is_some())
            && bim.with(|f| f.file_handle.is_some())
            && fam.with(|f| f.text.is_some())
    };

    let on_run = move |ev: MouseEvent| {
        ev.prevent_default();
        if is_running.get() {
            return;
        }
        if !all_loaded() {
            error.set(Some(
                "Please upload .bed, .bim, AND .fam before running.".to_string(),
            ));
            return;
        }
        error.set(None);
        result.set(None);
        progress.set((0, 0));
        is_running.set(true);
        // h² is gated on L2 success; reset its state too.
        h2_result.set(None);
        h2_error.set(None);
        bed_idx_for_h2.set(None);

        let preset_now = preset.get();
        let sketch_d_now = sketch_d.get();
        let snp_mask_now = snp_level_masking.get();
        let f32_now = fast_f32.get();
        let kb_now = ld_wind_kb.get();

        // BED + BIM both stay as `web_sys::File` handles — never
        // read into main-thread JS heap or wasm memory. The
        // per-chromosome worker pool streams both via FileReaderSync
        // inside each Worker. FAM is tiny (~1 MB at biobank scale)
        // and stays as a String passed through postMessage.
        //
        // We re-check each handle here (in addition to the
        // `all_loaded()` button-disable gate) because that gate is
        // client-side and a determined user could bypass it via
        // devtools. Surfacing as an `error.set(...)` is the safe
        // failure mode; the previous `.expect("bed loaded")` would
        // panic the wasm module.
        let Some(bed_file) = bed.with(|f| f.file_handle.clone()).map(|h| h.take()) else {
            error.set(Some("BED file not loaded".into()));
            is_running.set(false);
            return;
        };
        let Some(bim_file) = bim.with(|f| f.file_handle.clone()).map(|h| h.take()) else {
            error.set(Some("BIM file not loaded".into()));
            is_running.set(false);
            return;
        };
        let Some(fam_text) = fam.with(|f| f.text.clone()) else {
            error.set(Some("FAM file not loaded".into()));
            is_running.set(false);
            return;
        };

        let cfg = build_config(preset_now, sketch_d_now, snp_mask_now, f32_now, kb_now);
        let cfg_json = match serde_json::to_string(&cfg) {
            Ok(s) => s,
            Err(e) => {
                error.set(Some(format!("config JSON serialize: {e}")));
                is_running.set(false);
                return;
            }
        };

        tracing::info!(
            "Spawning L2 worker pool: BED size={} BIM size={} FAM len={}",
            bed_file.size(),
            bim_file.size(),
            fam_text.len()
        );

        let on_progress = move |p: PoolProgress| {
            progress.set((p.snps_done, p.snps_total));
        };
        let on_done = move |out: WorkerComputeResult| {
            tracing::info!("L2 done: {} SNPs in {:.2}s", out.n_snps, out.wall_seconds);
            // Stash bed_idx separately so the h² worker can request
            // it without making the rest of the UI thread through a
            // bigger payload.
            bed_idx_for_h2.set(Some(out.bed_idx.clone()));
            result.set(Some(out));
            is_running.set(false);
        };
        let on_err = move |msg: String| {
            tracing::error!("L2 worker pool failed: {msg}");
            error.set(Some(msg));
            is_running.set(false);
        };

        if let Err(e) = spawn_compute_l2_pool(
            bed_file,
            bim_file,
            fam_text,
            cfg_json,
            on_progress,
            on_done,
            on_err,
        ) {
            tracing::error!("spawn_compute_l2_pool failed: {e:?}");
            error.set(Some(format!("worker pool spawn: {e:?}")));
            is_running.set(false);
        }
    };

    // h² run handler. Fires the dedicated h² worker against the
    // already-computed L2 + bed_idx + the uploaded BIM and sumstats.
    // Gated on (L2 done && sumstats loaded) by the button's
    // `disabled` predicate below.
    //
    // `#[allow(unused_variables)]` on the binding: the Leptos `view!`
    // macro's nested-closure capture path defers `on_run_h2`'s
    // capture to a per-render closure, which rustc's static analyser
    // misses. The button works (verified end-to-end). Suppress the
    // false positive instead of an underscore prefix (which would
    // misleadingly imply the binding is dead).
    #[allow(unused_variables)]
    let on_run_h2 = move |ev: MouseEvent| {
        ev.prevent_default();
        if h2_running.get() {
            return;
        }
        let Some(l2_out) = result.get() else {
            h2_error.set(Some("Run L2 first.".into()));
            return;
        };
        let Some(bed_idx) = bed_idx_for_h2.get() else {
            h2_error.set(Some("L2 result is missing bed_idx — re-run L2.".into()));
            return;
        };
        let Some(bim_handle) = bim.with(|f| f.file_handle.clone()).map(|h| h.take()) else {
            h2_error.set(Some("BIM file not loaded.".into()));
            return;
        };
        let Some(sumstats_handle) = sumstats.with(|f| f.file_handle.clone()).map(|h| h.take())
        else {
            h2_error.set(Some(
                "Upload a .sumstats file (LDSC munge_sumstats output) to run h².".into(),
            ));
            return;
        };

        h2_error.set(None);
        h2_result.set(None);
        h2_running.set(true);

        let on_done_h2 = move |r: WireH2Result| {
            tracing::info!(
                "[h²] done: h²={:.4} (±{:.4}), intercept={:.4} (±{:.4}), mean χ²={:.4}, λ_GC={:.4}, n={}",
                r.h2,
                r.h2_se,
                r.intercept,
                r.intercept_se,
                r.mean_chi2,
                r.lambda_gc,
                r.n_snps_used,
            );
            h2_result.set(Some(r));
            h2_running.set(false);
        };
        let on_err_h2 = move |msg: String| {
            tracing::error!("[h²] failed: {msg}");
            h2_error.set(Some(msg));
            h2_running.set(false);
        };

        // LDSC defaults: 200 jackknife blocks, two-step at χ² = 30
        // (the standard cutoff in the original LDSC paper). Bypass
        // both knobs in v1 — they're exposed in the CLI but not in
        // the demo UI.
        if let Err(e) = spawn_compute_h2(
            l2_out.l2.clone(),
            l2_out.maf.clone(),
            bed_idx,
            bim_handle,
            sumstats_handle,
            200,
            Some(30.0),
            on_done_h2,
            on_err_h2,
        ) {
            tracing::error!("spawn_compute_h2 failed: {e:?}");
            h2_error.set(Some(format!("h² worker spawn: {e:?}")));
            h2_running.set(false);
        }
    };

    view! {
        // ── Upload card ────────────────────────────────────────────────
        <div class="card">
            <div class="card-header">
                "L2 — LD scores from .bed/.bim/.fam"
                <span class="subtitle">"Upload your PLINK genotype trio (.sumstats optional, for h²)"</span>
            </div>
            <div class="card-body">
                <FileUploadRow ext=".bed" read_as=ReadAs::BlobHandle file=bed />
                <FileUploadRow ext=".bim" read_as=ReadAs::BlobHandle file=bim />
                <FileUploadRow ext=".fam" read_as=ReadAs::Text file=fam />
                <FileUploadRow
                    ext=".sumstats"
                    accept=".sumstats,.sumstats.gz"
                    read_as=ReadAs::BlobHandle
                    file=sumstats
                />
                <p class="text-muted mt-2 mb-0" style="font-size: 0.825rem;">
                    "Files stay in your browser — no data is uploaded anywhere. \
                    Compute runs locally via WebAssembly. The .sumstats file is optional \
                    (uncompressed or .gz — both work; LDSC's munge_sumstats writes .gz \
                    by default); upload it to also run h² heritability after L2 finishes."
                </p>
            </div>
        </div>

        // ── Configuration card ─────────────────────────────────────────
        <div class="card">
            <div class="card-header">
                "Mode"
                <span class="subtitle">"Pick a speed / accuracy trade-off"</span>
            </div>
            <div class="card-body">
                <ModeRadios preset=preset />
                {move || {
                    if matches!(preset.get(), ModePreset::Sketch) {
                        view! { <SketchSlider sketch_d=sketch_d fam=fam /> }.into_any()
                    } else {
                        view! { <span></span> }.into_any()
                    }
                }}

                <hr class="my-3" />

                <div class="row g-3">
                    <div class="col-sm-6">
                        <div class="form-check">
                            <input class="form-check-input" type="checkbox" id="opt-snp-mask"
                                prop:checked=move || snp_level_masking.get()
                                on:change=move |_| snp_level_masking.update(|b| *b = !*b) />
                            <label class="form-check-label" for="opt-snp-mask">
                                <code>"--snp-level-masking"</code>
                                <span class="text-muted ms-2" style="font-size: 0.825rem;">
                                    "(per-SNP exact windows; LDSC paper math)"
                                </span>
                            </label>
                        </div>
                    </div>
                    <div class="col-sm-6">
                        <div class="form-check">
                            <input class="form-check-input" type="checkbox" id="opt-fast-f32"
                                prop:checked=move || fast_f32.get()
                                on:change=move |_| fast_f32.update(|b| *b = !*b) />
                            <label class="form-check-label" for="opt-fast-f32">
                                <code>"--fast-f32"</code>
                                <span class="text-muted ms-2" style="font-size: 0.825rem;">
                                    "(default — f32 GEMM, ~1.85× faster than f64; sketch modes implicitly use f32)"
                                </span>
                            </label>
                        </div>
                    </div>
                </div>

                <div class="row g-3 mt-2">
                    <div class="col-sm-6">
                        <label class="form-label" for="opt-ld-wind">
                            <code>"--ld-wind-kb"</code>
                        </label>
                        <input class="form-control" type="number" id="opt-ld-wind"
                            min="10" max="100000" step="100"
                            prop:value=move || ld_wind_kb.get()
                            on:input=move |ev| {
                                if let Ok(v) = event_target_value(&ev).parse::<f64>() {
                                    ld_wind_kb.set(v);
                                }
                            } />
                    </div>
                </div>
            </div>
        </div>

        // ── Run button ─────────────────────────────────────────────────
        <div class="d-flex justify-content-end mb-3">
            <button class="btn btn-primary btn-lg"
                disabled=move || is_running.get() || !all_loaded()
                on:click=on_run>
                {move || if is_running.get() { "Running…" } else { "Run" }}
            </button>
        </div>

        // ── Progress bar (visible while running) ───────────────────────
        {move || {
            if !is_running.get() {
                return view! { <span></span> }.into_any();
            }
            let (done, total) = progress.get();
            let pct = if total > 0 {
                (done as f64 / total as f64 * 100.0).clamp(0.0, 100.0)
            } else {
                0.0
            };
            let label = if total == 0 {
                "Scanning BIM…".to_string()
            } else {
                format!("{} / {} SNPs ({:.1}%)", format_int(done), format_int(total), pct)
            };
            view! {
                <div class="mb-3">
                    <div class="progress" role="progressbar" style="height: 1.25rem;">
                        <div class="progress-bar progress-bar-striped progress-bar-animated"
                             style=move || format!("width: {:.1}%; background-color: #2a71a5;", pct)>
                        </div>
                    </div>
                    <div class="text-muted text-end" style="font-size: 0.825rem;">
                        {label}
                    </div>
                </div>
            }.into_any()
        }}

        // ── Error box ──────────────────────────────────────────────────
        {move || error.get().map(|msg| view! {
            <div class="error-box">{msg}</div>
        })}

        // ── Results card ───────────────────────────────────────────────
        {move || result.get().map(|out| {
            let stats = RunStats::from_output(&out);
            let perf_view = perf_breakdown_view(&out);
            let (corr_class, corr_note) = classify_maf_l2_corr(stats.maf_l2_corr);
            view! {
                <div class="card">
                    <div class="card-header">"L2 — LD Scores"</div>
                    <div class="card-body">
                        // Tight stat grid: the four numbers Python LDSC's
                        // CLI summary actually prints (wall, SNPs, mean L2,
                        // MAF/L2 r). Median / max / mean MAF dropped — the
                        // research summary in docs/perf-log.md notes only
                        // mean L2 + MAF/L2 r are routinely cited.
                        <div class="stat-grid">
                            <div class="stat">
                                <div class="label">"Wall time"</div>
                                <div class="value">{format!("{:.2} s", stats.wall_seconds)}</div>
                            </div>
                            <div class="stat">
                                <div class="label">"SNPs"</div>
                                <div class="value">{format_int(stats.n_snps)}</div>
                            </div>
                            <div class="stat">
                                <div class="label">"Mean L2"</div>
                                <div class="value">{format!("{:.3}", stats.mean_l2)}</div>
                            </div>
                            <div class="stat">
                                <div class="label">"MAF / L2 r"</div>
                                <div class="value">{format!("{:.3}", stats.maf_l2_corr)}</div>
                                <div class=corr_class style="font-size: 0.7rem; margin-top: 0.15rem;">
                                    {corr_note}
                                </div>
                            </div>
                        </div>

                        <canvas id="chart-scatter" class="chart" width="880" height="320"></canvas>
                        <p class="text-muted mt-1 mb-0" style="font-size: 0.75rem;">
                            "QC: MAF and L2 should be positively correlated (LDSC wiki: r ≈ 0.27 in EUR panels). \
                             Negative r → check for MAF filter / panel issues."
                        </p>

                        <details class="mt-3">
                            <summary class="text-muted" style="font-size: 0.825rem;">
                                "Compute timing breakdown (where the wall went)"
                            </summary>
                            <div class="mt-2">{perf_view}</div>
                        </details>
                    </div>
                </div>
            }
        })}

        // ── h² heritability card (gated on L2 result + sumstats) ─────
        //
        // Renders below the L2 results card whenever L2 has finished.
        // The button is disabled until the user uploads a .sumstats
        // file, with inline help to point them at LDSC's munge_sumstats.
        {move || {
            // Show the card only after the user has finished an L2 run.
            // (Before that, the .sumstats picker is still visible in
            // the upload card, but there's nothing to regress against.)
            if result.get().is_none() {
                return view! { <span></span> }.into_any();
            }
            let h2_ready = sumstats.with(|f| f.file_handle.is_some());
            view! {
                <div class="card">
                    <div class="card-header">
                        "h² — Heritability"
                        <span class="subtitle">
                            "LDSC regression with block jackknife SE; 5-tuple matches the CLI output"
                        </span>
                    </div>
                    <div class="card-body">
                        {move || if !h2_ready {
                            view! {
                                <p class="text-muted mb-0" style="font-size: 0.9rem;">
                                    "Upload a "<code>".sumstats"</code>" or "
                                    <code>".sumstats.gz"</code>" file (output of "
                                    <code>"ldsc munge-sumstats"</code>") in the upload card \
                                    above to run heritability regression against the LD scores \
                                    you just computed. Need columns "<code>"SNP"</code>", "
                                    <code>"Z"</code>", "<code>"N"</code>"."
                                </p>
                            }.into_any()
                        } else {
                            view! {
                                <div class="d-flex justify-content-end mb-3">
                                    <button class="btn btn-primary"
                                        disabled=move || h2_running.get()
                                        on:click=on_run_h2>
                                        {move || if h2_running.get() {
                                            "Running h²…"
                                        } else {
                                            "Run h²"
                                        }}
                                    </button>
                                </div>
                            }.into_any()
                        }}

                        {move || h2_error.get().map(|msg| view! {
                            <div class="error-box">{msg}</div>
                        })}

                        {move || h2_result.get().map(|r| view! {
                            <H2ResultCard r=r />
                        })}
                    </div>
                </div>
            }.into_any()
        }}
    }
}

/// Renders the 5-tuple stat-geneticists quote in papers, plus the
/// MAF/L2 r badge equivalent (sign-and-magnitude checks on h² and
/// intercept).
#[component]
fn H2ResultCard(r: WireH2Result) -> impl IntoView {
    // Quick sanity classifiers — h² should be in [0, 1] for a sane
    // run on a single complex trait; LDSC commonly returns slightly
    // negative values for noisy / low-h² traits (those are valid
    // point estimates that the SE swallows). Intercept ≈ 1 means
    // no confounding; > 1.05 hints at stratification.
    let h2_class = if (0.0..=1.0).contains(&r.h2) {
        "h2-good"
    } else {
        "h2-warn"
    };
    let intercept_class = if (0.95..=1.05).contains(&r.intercept) {
        "h2-good"
    } else if (0.90..=1.10).contains(&r.intercept) {
        "h2-warn"
    } else {
        "h2-bad"
    };
    let ratio_view = match (r.ratio, r.ratio_se) {
        (Some(rv), Some(se)) => view! {
            <div class="stat">
                <div class="label">"Ratio"</div>
                <div class="value">{format!("{:.3}", rv)}</div>
                <div class="text-muted" style="font-size: 0.7rem; margin-top: 0.15rem;">
                    {format!("± {:.3} SE", se)}
                </div>
            </div>
        }
        .into_any(),
        _ => view! {
            <div class="stat">
                <div class="label">"Ratio"</div>
                <div class="value">"—"</div>
                <div class="text-muted" style="font-size: 0.7rem; margin-top: 0.15rem;">
                    "NA (mean χ² ≈ 1)"
                </div>
            </div>
        }
        .into_any(),
    };
    let snps_used = r.n_snps_used;
    let snps_sumstats = r.n_snps_sumstats;
    let m_snps_ref = r.m_snps_ref;
    let mean_n = r.mean_n;
    let wall = r.wall_seconds;
    view! {
        <div class="stat-grid">
            <div class="stat">
                <div class="label">"h²"</div>
                <div class="value">{format!("{:.4}", r.h2)}</div>
                <div class=h2_class style="font-size: 0.7rem; margin-top: 0.15rem;">
                    {format!("± {:.4} SE", r.h2_se)}
                </div>
            </div>
            <div class="stat">
                <div class="label">"Intercept"</div>
                <div class="value">{format!("{:.4}", r.intercept)}</div>
                <div class=intercept_class style="font-size: 0.7rem; margin-top: 0.15rem;">
                    {format!("± {:.4} SE", r.intercept_se)}
                </div>
            </div>
            {ratio_view}
            <div class="stat">
                <div class="label">"Mean χ²"</div>
                <div class="value">{format!("{:.4}", r.mean_chi2)}</div>
            </div>
            <div class="stat">
                <div class="label">"λ_GC"</div>
                <div class="value">{format!("{:.4}", r.lambda_gc)}</div>
            </div>
        </div>

        <p class="text-muted mt-3 mb-0" style="font-size: 0.825rem;">
            {format!(
                "Joined {} of {} sumstats SNPs against {} LD-score reference SNPs (mean N = {:.0}). \
                 Regression completed in {:.3} s.",
                format_int(snps_used),
                format_int(snps_sumstats),
                format_int(m_snps_ref),
                mean_n,
                wall,
            )}
        </p>

        <details class="mt-3">
            <summary class="text-muted" style="font-size: 0.825rem;">
                "How to read this"
            </summary>
            <div class="mt-2 text-muted" style="font-size: 0.825rem;">
                <p class="mb-2">
                    <strong>"h²"</strong>" — narrow-sense heritability captured by common SNPs in
                    the reference panel. For most complex traits this is roughly 0.1–0.5; values
                    near 0 (or slightly negative) are LDSC's natural point estimate for
                    low-heritability traits — the SE is the honest measure of uncertainty."
                </p>
                <p class="mb-2">
                    <strong>"Intercept"</strong>" — should be ≈ 1.0 under the LDSC null (no
                    confounding). Values > 1.05 suggest population stratification or cryptic
                    relatedness; the "<strong>"Ratio"</strong>" = (intercept − 1) / (mean χ² − 1)
                    is the share of inflation attributable to confounding (Bulik-Sullivan 2015)."
                </p>
                <p class="mb-0">
                    <strong>"Mean χ², λ_GC"</strong>" — pre-regression diagnostics. λ_GC > 1 is
                    expected for polygenic traits and does not by itself indicate confounding."
                </p>
            </div>
        </details>
    }
}

/// Categorise the MAF/L2 Pearson correlation against the LDSC
/// wiki's expectation (positive, r ≈ 0.27 in EUR panels). Returns
/// (css-class, short-note) for the UI badge under the stat tile.
fn classify_maf_l2_corr(r: f64) -> (&'static str, &'static str) {
    if r >= 0.10 {
        ("h2-good", "✓ positive (typical EUR ≈ 0.27)")
    } else if r >= 0.0 {
        ("h2-warn", "weak (typical EUR ≈ 0.27)")
    } else {
        ("h2-bad", "⚠ negative — check MAF filter / panel")
    }
}

/// Per-phase wall-time breakdown for the long-tail outer (the worker
/// whose `wall_seconds` matches the overall pool wall). Pre-H this
/// is dominated by `bed_read(stall)`; Workstream H+ should shift the
/// dominant phase toward `sketch` / GEMM. Surfaces the data the
/// `[perf]` console log already prints, so users can see where the
/// wall is going without opening DevTools.
fn perf_breakdown_view(out: &WorkerComputeResult) -> impl IntoView + use<> {
    if out.per_worker_perf.is_empty() {
        return view! { <div></div> }.into_any();
    }
    // Find the long-tail worker (max wall) and use it for the
    // breakdown. This is the worker that defines the overall pool
    // wall, so its phases are the ones that matter for optimisation.
    let (long_tail_idx, long_tail_wall) = out
        .per_worker_wall_seconds
        .iter()
        .copied()
        .enumerate()
        .fold(
            (0usize, 0.0f64),
            |acc, (i, w)| if w > acc.1 { (i, w) } else { acc },
        );
    let p = &out.per_worker_perf[long_tail_idx];

    // Each row: (label, seconds, css-class). Display only non-zero
    // rows so sketch-only runs don't show empty GEMM rows.
    let mut rows: Vec<(String, f64, &'static str)> = Vec::new();
    rows.push((
        "BED read (FileReaderSync stall)".to_string(),
        p.bed_read_secs,
        "phase-io",
    ));
    rows.push(("Normalize".to_string(), p.norm_secs, "phase-cpu"));
    if let Some(s) = p.sketch_secs {
        rows.push(("Sketch (scatter-add)".to_string(), s, "phase-cpu"));
    }
    if p.bb_dot_secs > 0.0 {
        rows.push((
            "GEMM B^T·B (within-chunk)".to_string(),
            p.bb_dot_secs,
            "phase-cpu",
        ));
    }
    if p.ab_dot_secs > 0.0 {
        rows.push((
            "GEMM A^T·B (cross-window)".to_string(),
            p.ab_dot_secs,
            "phase-cpu",
        ));
    }
    if p.ring_store_secs > 0.0 {
        rows.push(("Ring buffer".to_string(), p.ring_store_secs, "phase-cpu"));
    }

    // Bar widths normalised to the long-tail worker's wall.
    let max_bar = long_tail_wall.max(1e-6);
    let total_accounted: f64 = rows.iter().map(|r| r.1).sum();
    let other = (long_tail_wall - total_accounted).max(0.0);
    if other > 0.01 {
        rows.push(("Other / dispatch".to_string(), other, "phase-other"));
    }

    let n_workers = out.per_worker_wall_seconds.len();
    let row_views = rows
        .into_iter()
        .map(|(label, secs, css)| {
            let pct = (secs / max_bar * 100.0).clamp(0.0, 100.0);
            let bar_color = match css {
                "phase-io" => "#d97706",  // amber — I/O stall
                "phase-cpu" => "#2a71a5", // blue — CPU work
                _ => "#9ca3af",           // grey — other
            };
            view! {
                <div class="perf-row">
                    <div class="perf-label">{label}</div>
                    <div class="perf-bar-wrap">
                        <div class="perf-bar"
                             style=format!("width: {:.1}%; background-color: {};",
                                           pct, bar_color)>
                        </div>
                    </div>
                    <div class="perf-value">
                        {format!("{:.2} s ({:.0}%)", secs, secs / max_bar * 100.0)}
                    </div>
                </div>
            }
        })
        .collect::<Vec<_>>();

    view! {
        <div class="perf-breakdown mt-4">
            <div class="perf-header">
                <strong>"Where the wall went — worker[" {long_tail_idx} "] (long tail of "
                {n_workers} ", chr handled = " {p.m} " SNPs)"</strong>
                <span class="text-muted ms-2" style="font-size: 0.825rem;">
                    "Pool wall = max across workers = "
                    {format!("{:.2}s", out.wall_seconds)}
                </span>
            </div>
            <div class="perf-rows mt-2">
                {row_views}
            </div>
            <div class="text-muted mt-1" style="font-size: 0.75rem;">
                "Amber = file I/O stall. Blue = CPU work (SIMD compute)."
            </div>
        </div>
    }
    .into_any()
}

#[component]
fn ModeRadios(preset: RwSignal<ModePreset>) -> impl IntoView {
    let radio = move |m: ModePreset| {
        let id = format!("mode-{:?}", m);
        let id_for_label = id.clone();
        let active = move || preset.get() == m;
        view! {
            <div class="form-check mb-2">
                <input class="form-check-input"
                    type="radio"
                    name="mode-preset"
                    id=id
                    prop:checked=active
                    on:change=move |_| preset.set(m) />
                <label class="form-check-label" for=id_for_label>
                    <strong>{m.label()}</strong>
                    <div class="text-muted" style="font-size: 0.825rem;">
                        {m.description()}
                    </div>
                </label>
            </div>
        }
    };

    view! {
        <div>
            {radio(ModePreset::Sketch)}
            {radio(ModePreset::Default)}
            {radio(ModePreset::PythonCompat)}
            {radio(ModePreset::PerSnpExact)}
        </div>
    }
}

#[component]
fn SketchSlider(sketch_d: RwSignal<u32>, fam: RwSignal<LoadedFile>) -> impl IntoView {
    // Inline these helpers — Leptos's `#[component]` macro reparents
    // the body in a way that snake_case module-level fns can't always
    // be resolved from. Inlining sidesteps that wart cleanly.
    //
    // **N×d sweep, M5 Pro, --snp-level-masking, BMI Yengo 2018**
    // (`preprint/data/dn_sweep_full.csv`): wall, Pearson r vs exact,
    // and |Δh²| vs exact at each (N, d) cell:
    //
    //   d      N=503   N=10K   N=20K   N=50K   N=100K   r at d  |Δh²| at d=1000
    //   200    2.3s    3.2s    3.4s    7.1s    12.4s    ~0.967  0.003-0.007
    //   500    3.0s    10.5s   4.8s    7.7s    13.8s    ~0.987  0.000-0.006
    //   1000   ---     7.4s    6.4s    9.7s    14.4s    ~0.994  0.000-0.003
    //   1600   ---     9.6s    13.3s   12.1s   17.5s    ~0.996  0.000-0.002
    //   5000   ---     23.9s   29.4s   27.4s   32.1s    ~0.999  0.000-0.002
    //
    // **d=1000 is the universal truth-cluster sweet spot** (Pearson
    // r ≥ 0.993 and h² within one regression-SE of per-SNP exact
    // across all N). Wall is <15s even at N=100K. d=200 is faster
    // but drifts h² by 0.007-0.013 from truth (not in the cluster).
    // For strict per-SNP accuracy (e.g. partitioned heritability)
    // bump to d=5000 for r ≥ 0.999.
    let pick_d = |_n_indiv: usize| -> u32 { 1000 };
    let n_from_fam =
        |text: &str| -> usize { text.lines().filter(|l| !l.trim().is_empty()).count() };
    let recommended = move || {
        fam.with(|f| {
            f.text.as_ref().map(|t| {
                let n = n_from_fam(t);
                (n, pick_d(n))
            })
        })
    };
    view! {
        <div class="mt-2 mb-3" style="margin-left: 1.75rem;">
            <label class="form-label" for="sketch-d">
                <code>"--sketch d"</code>
                " = "
                <strong>{move || sketch_d.get()}</strong>
            </label>
            <div class="text-muted" style="font-size: 0.78rem; margin-bottom: 0.5rem;">
                "Conceptually: "
                <code>"d"</code>
                " is the number of \"synthetic individuals\" we compress \
                 your N real individuals into. CountSketch randomly hashes \
                 each individual into one of "
                <code>"d"</code>
                " buckets with a random ±1 sign, then sums them. Pairwise \
                 SNP inner products (LD scores) survive with per-pair noise \
                 ∝ 1/d. Compute scales with N, accuracy scales with d alone \
                 — so the right "
                <code>"d"</code>
                " is the same regardless of panel size."
            </div>
            <input class="form-range" type="range" id="sketch-d"
                min="100" max="3200" step="100"
                prop:value=move || sketch_d.get() as f64
                on:input=move |ev| {
                    if let Ok(v) = event_target_value(&ev).parse::<u32>() {
                        sketch_d.set(v);
                    }
                } />
            <div class="d-flex justify-content-between text-muted" style="font-size: 0.75rem;">
                <span>"100"</span><span>"500"</span><span>"1000 (rec.)"</span><span>"2000"</span><span>"3200"</span>
            </div>
            {move || match recommended() {
                Some((n, opt)) => view! {
                    <div class="mt-2">
                        <button type="button" class="btn btn-sm btn-outline-primary"
                            on:click=move |_| sketch_d.set(opt)>
                            "⚡ Speed-optimize for your data (d=" {opt} ", N=" {n} ")"
                        </button>
                        <div class="text-muted mt-1" style="font-size: 0.75rem;">
                            "Sets d to the empirically-validated truth-cluster sweet spot \
                             from the N×d sweep: h² within ~0.003 of per-SNP exact across \
                             N=503 to N=100K, Pearson r ≥ 0.993 universally, wall < 15s \
                             even at N=100K. d=200 would be faster (~2 s saved) but drifts \
                             h² by ~0.01 from truth."
                        </div>
                    </div>
                }.into_any(),
                None => view! {
                    <div class="text-muted mt-2" style="font-size: 0.75rem;">
                        "Upload a .fam file to enable the speed-optimize button."
                    </div>
                }.into_any(),
            }}
            <div class="text-muted mt-2" style="font-size: 0.825rem;">
                "Larger d = more accurate, slower. d=1000 is the empirically-validated \
                 universal sweet spot (truth-cluster across all sample sizes). \
                 d=5000 for r ≥ 0.999 if you need per-SNP accuracy. \
                 d ≤ 50 is unstable and is rejected by the library."
            </div>
        </div>
    }
}

/// Translate the UI mode + flags into a JSON-serialisable wire
/// config. We send this across the worker postMessage boundary, so
/// it has to be `Serialize` (= `WireL2Config` from `worker::*`).
fn build_config(
    preset: ModePreset,
    sketch_d: u32,
    snp_mask: bool,
    fast_f32: bool,
    ld_wind_kb: f64,
) -> WireL2Config {
    let mut cfg = WireL2Config {
        mode: WireWindowMode::Kb(ld_wind_kb),
        chunk_size: 200,
        use_f32: fast_f32,
        sketch: None,
        sketch_maf_aware: false,
        snp_level_masking: snp_mask,
        yes_really: true,
        pq_exp: None,
    };
    match preset {
        ModePreset::Sketch => {
            cfg.sketch = Some(sketch_d as usize);
            // Sketch implies f32 internally; the toggle still works
            // as a no-op explicit hint.
            cfg.use_f32 = true;
        }
        ModePreset::Default => {
            // chunk_size=200 (default), no sketch.
        }
        ModePreset::PythonCompat => {
            cfg.chunk_size = 50;
            cfg.snp_level_masking = false;
        }
        ModePreset::PerSnpExact => {
            cfg.snp_level_masking = true;
        }
    }
    cfg
}

fn mean(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return f64::NAN;
    }
    let n = xs.len() as f64;
    xs.iter().filter(|x| x.is_finite()).sum::<f64>() / n
}

fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len().min(y.len());
    if n < 2 {
        return f64::NAN;
    }
    let mx = mean(&x[..n]);
    let my = mean(&y[..n]);
    let mut cov = 0.0;
    let mut vx = 0.0;
    let mut vy = 0.0;
    for i in 0..n {
        let dx = x[i] - mx;
        let dy = y[i] - my;
        cov += dx * dy;
        vx += dx * dx;
        vy += dy * dy;
    }
    if vx == 0.0 || vy == 0.0 {
        return f64::NAN;
    }
    cov / (vx * vy).sqrt()
}

fn format_int(n: usize) -> String {
    let s = n.to_string();
    let chars: Vec<char> = s.chars().rev().collect();
    let mut out = String::new();
    for (i, c) in chars.iter().enumerate() {
        if i > 0 && i % 3 == 0 {
            out.push(' ');
        }
        out.push(*c);
    }
    out.chars().rev().collect()
}

fn request_animation_frame(f: impl FnOnce() + 'static) {
    use wasm_bindgen::JsCast;
    use wasm_bindgen::closure::Closure;
    let cb = Closure::once_into_js(f);
    let _ = web_sys::window()
        .expect("window")
        .request_animation_frame(cb.unchecked_ref());
}
