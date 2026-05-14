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

use super::charts::{percentile, render_l2_histogram, render_l2_vs_maf};
use super::file_upload::{FileUploadRow, LoadedFile, ReadAs};
use super::worker_client::{WorkerComputeResult, spawn_compute_l2};
use crate::worker::{WireL2Config, WireWindowMode};

/// Mode preset radio (mirrors the "Mode" section of the CLI README's
/// recommendations table). Selecting one auto-toggles the underlying
/// flags in `effective_config`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ModePreset {
    /// Default-recommended: CountSketch d=1600 (~17× speedup, h² within
    /// 0.001 of Python LDSC at biobank scale).
    Sketch,
    /// Chunked exact (Rust default; chunk_size=200 — slightly inflates
    /// L2 vs Python by mean ~0.2 / max ~6).
    Default,
    /// Bit-identical to Python LDSC. chunk_size=50, global_pass.
    PythonCompat,
    /// LDSC-paper-canonical per-SNP exact windows. Combine with Sketch
    /// for the headline `--sketch 1600 --snp-level-masking` mode.
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
                 h² within ~0.001 of Python LDSC at d=1600."
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
    median_l2: f64,
    max_l2: f64,
    mean_maf: f64,
    maf_l2_corr: f64,
    wall_seconds: f64,
}

impl RunStats {
    fn from_output(out: &WorkerComputeResult) -> Self {
        let n_snps = out.l2.len();
        let mean_l2 = mean(&out.l2);
        let median_l2 = percentile(&out.l2, 0.5).unwrap_or(0.0);
        let max_l2 = out.l2.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let mean_maf = mean(&out.maf);
        let maf_l2_corr = pearson(&out.maf, &out.l2);
        Self {
            n_snps,
            mean_l2,
            median_l2,
            max_l2,
            mean_maf,
            maf_l2_corr,
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

    // ── Config state ───────────────────────────────────────────────────
    let preset = RwSignal::new(ModePreset::Sketch);
    let sketch_d = RwSignal::new(1600u32);
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

    // Re-render charts whenever a new result lands.
    Effect::new(move |_| {
        if let Some(out) = result.get() {
            // Defer one tick so the canvas elements exist in the DOM.
            request_animation_frame(move || {
                if let Err(e) = render_l2_histogram("chart-hist", &out.l2) {
                    tracing::error!("histogram render: {e:?}");
                }
                if let Err(e) = render_l2_vs_maf("chart-scatter", &out.maf, &out.l2) {
                    tracing::error!("scatter render: {e:?}");
                }
            });
        }
    });

    let all_loaded = move || {
        bed.with(|f| f.file_handle.is_some())
            && bim.with(|f| f.text.is_some())
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
        is_running.set(true);

        let preset_now = preset.get();
        let sketch_d_now = sketch_d.get();
        let snp_mask_now = snp_level_masking.get();
        let f32_now = fast_f32.get();
        let kb_now = ld_wind_kb.get();

        // BED stays as a `web_sys::File` handle the whole way —
        // never read into wasm memory on the main thread. Worker
        // streams chunks via FileReaderSync. BIM + FAM are small
        // enough to be already in memory as Strings on the main
        // thread; we just clone-move them into the worker call.
        let bed_file = bed
            .with(|f| f.file_handle.clone())
            .expect("bed loaded")
            .take();
        let bim_text = bim.with(|f| f.text.clone()).expect("bim loaded");
        let fam_text = fam.with(|f| f.text.clone()).expect("fam loaded");

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
            "Spawning L2 worker: BED size={} bim={} fam={}",
            bed_file.size(),
            bim_text.len(),
            fam_text.len()
        );

        let on_done = move |out: WorkerComputeResult| {
            tracing::info!("L2 done: {} SNPs in {:.2}s", out.n_snps, out.wall_seconds);
            result.set(Some(out));
            is_running.set(false);
        };
        let on_err = move |msg: String| {
            tracing::error!("L2 worker failed: {msg}");
            error.set(Some(msg));
            is_running.set(false);
        };

        if let Err(e) = spawn_compute_l2(bed_file, bim_text, fam_text, cfg_json, on_done, on_err) {
            tracing::error!("spawn_compute_l2 failed: {e:?}");
            error.set(Some(format!("worker spawn: {e:?}")));
            is_running.set(false);
        }
    };

    view! {
        // ── Upload card ────────────────────────────────────────────────
        <div class="card">
            <div class="card-header">
                "L2 — LD scores from .bed/.bim/.fam"
                <span class="subtitle">"Upload your PLINK genotype trio"</span>
            </div>
            <div class="card-body">
                <FileUploadRow ext=".bed" read_as=ReadAs::BlobHandle file=bed />
                <FileUploadRow ext=".bim" read_as=ReadAs::Text file=bim />
                <FileUploadRow ext=".fam" read_as=ReadAs::Text file=fam />
                <p class="text-muted mt-2 mb-0" style="font-size: 0.825rem;">
                    "Files stay in your browser — no data is uploaded anywhere. \
                    Compute runs locally via WebAssembly."
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
                        view! { <SketchSlider sketch_d=sketch_d /> }.into_any()
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

        // ── Error box ──────────────────────────────────────────────────
        {move || error.get().map(|msg| view! {
            <div class="error-box">{msg}</div>
        })}

        // ── Results card ───────────────────────────────────────────────
        {move || result.get().map(|out| {
            let stats = RunStats::from_output(&out);
            view! {
                <div class="card">
                    <div class="card-header">"Results"</div>
                    <div class="card-body">
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
                                <div class="label">"Median L2"</div>
                                <div class="value">{format!("{:.3}", stats.median_l2)}</div>
                            </div>
                            <div class="stat">
                                <div class="label">"Max L2"</div>
                                <div class="value">{format!("{:.2}", stats.max_l2)}</div>
                            </div>
                            <div class="stat">
                                <div class="label">"Mean MAF"</div>
                                <div class="value">{format!("{:.3}", stats.mean_maf)}</div>
                            </div>
                            <div class="stat">
                                <div class="label">"MAF / L2 r"</div>
                                <div class="value">{format!("{:.3}", stats.maf_l2_corr)}</div>
                            </div>
                        </div>

                        <canvas id="chart-hist" class="chart" width="880" height="320"></canvas>
                        <canvas id="chart-scatter" class="chart" width="880" height="320"></canvas>
                    </div>
                </div>
            }
        })}
    }
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
fn SketchSlider(sketch_d: RwSignal<u32>) -> impl IntoView {
    view! {
        <div class="mt-2 mb-3" style="margin-left: 1.75rem;">
            <label class="form-label" for="sketch-d">
                <code>"--sketch d"</code>
                " = "
                <strong>{move || sketch_d.get()}</strong>
            </label>
            <input class="form-range" type="range" id="sketch-d"
                min="100" max="3200" step="100"
                prop:value=move || sketch_d.get() as f64
                on:input=move |ev| {
                    if let Ok(v) = event_target_value(&ev).parse::<u32>() {
                        sketch_d.set(v);
                    }
                } />
            <div class="d-flex justify-content-between text-muted" style="font-size: 0.75rem;">
                <span>"100"</span><span>"800"</span><span>"1600"</span><span>"2400"</span><span>"3200"</span>
            </div>
            <div class="text-muted" style="font-size: 0.825rem;">
                "Larger d = more accurate, slower. d=1600 is the headline sweet spot. \
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
        verbose_timing: false,
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
