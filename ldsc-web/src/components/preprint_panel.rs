//! Preprint reader — iframes the standalone HTML rendering of
//! `preprint/main.typ` (built by `preprint/scripts/build_web_artifact.sh`
//! and vendored at `ldsc-web/assets/preprint.html`).
//!
//! The vendored HTML is a single self-contained doc: one big inline
//! SVG with all 20 pages stacked, plus a `.tsel` semantic-text
//! overlay that makes the rendered text selectable and Ctrl-F
//! searchable in any modern browser. No JS bridge, no typst.ts
//! renderer WASM, no dynamic imports — just an `<iframe src=...>`
//! pointed at the pre-rendered file.
//!
//! Architecture trade-offs (decided 2026-05-19):
//!   - svg_html: 8.1 MB raw, 2.17 MB gzipped on the wire.
//!   - typst.ts + vector IR alternative: 3.16 MB combined gzipped
//!     (1 MB renderer WASM + ~2 MB IR + 83 KB JS shim), MORE
//!     complex, MORE moving parts. We picked svg_html.
//!
//! Mobile note: the iframe has a fixed CSS height (`--preprint-h`
//! in palette.css, currently 80vh) and scrolls internally. iOS
//! Safari renders the iframe lazily, so the giant SVG inside
//! doesn't get rasterised until scrolled into the iframe's
//! visible region — no memory-ceiling concerns at 612 × 16632
//! page extent.

use leptos::prelude::*;

#[component]
pub fn PreprintPanel() -> impl IntoView {
    view! {
        <div class="card">
            <div class="card-header d-flex justify-content-between align-items-center">
                <span>
                    "Preprint — ldsc-rs: a fast Rust reimplementation of LD Score Regression"
                </span>
                // Download fallback. `download` attr asks the
                // browser to save rather than navigate; the
                // explicit filename keeps the saved file from
                // being called the generic "preprint.pdf" if the
                // user is bulk-downloading from multiple sources.
                <a href="preprint.pdf"
                   download="ldsc-rs-preprint.pdf"
                   class="btn btn-sm btn-outline-primary"
                   title="Download the high-fidelity PDF (2 MB)">
                    "Download PDF"
                </a>
            </div>
            // `p-0` strips the card-body's default padding so the
            // iframe goes edge-to-edge inside the card — looks
            // cleaner than a framed reader with extra inner gutter.
            <div class="card-body p-0">
                <iframe src="preprint.html"
                        class="preprint-iframe"
                        title="ldsc-rs preprint (rendered with typst.ts svg_html)"
                        // sandbox: allow-same-origin lets in-doc
                        // cross-refs resolve and lets the embedded
                        // .tsel selection styling inherit. allow-
                        // scripts is needed for the tiny ~40-LOC
                        // drag-selection shim injected by
                        // preprint/scripts/build_web_artifact.sh —
                        // typst.ts svg_html puts each text run in
                        // its own <foreignObject>, and browsers
                        // refuse to extend a Range across that
                        // boundary during mouse drag, so the shim
                        // drives selection manually via
                        // caretPositionFromPoint. The injected
                        // script is our own and the iframe loads
                        // a same-origin static asset (preprint.html)
                        // — no third-party JS runs here.
                        sandbox="allow-same-origin allow-scripts"
                        loading="lazy"></iframe>
            </div>
        </div>
    }
}
