//! Preprint reader — iframes the canonical PDF directly and lets
//! the browser's built-in PDF viewer handle rendering, selection,
//! scrolling, zoom, and Ctrl-F search.
//!
//! Why the browser's PDF viewer instead of a custom HTML rendering?
//! Each major browser ships a battle-tested PDF reader:
//!   - Chrome / Edge: PDFium (Google's C++ engine, hardened by
//!     years of fuzzing).
//!   - Firefox: pdf.js (Mozilla's JS engine, used by hundreds of
//!     millions of Firefox users daily).
//!   - Safari: PDFKit (Apple's Cocoa framework, same engine as
//!     Preview.app).
//! Drag-selection, range-select-across-pages, copy-paste, find-on-
//! page — all native, all rock-solid, zero JS we maintain.
//!
//! An earlier iteration of this component iframed a typst.ts
//! `svg_html` rendering (preprint.html, 8.1 MB raw / 2.17 MB
//! gzipped) and tried to compensate for typst.ts wrapping every
//! text run in its own SVG <foreignObject> (6,946 of them!) with
//! a hand-rolled selection shim. The shim worked but was custom
//! code we owned; the PDF iframe is one HTML attribute change
//! and the selection problem disappears entirely. The PDF is the
//! same 2 MB payload either way (it was already the download
//! fallback under the SVG approach), so this is a strict
//! simplification — less JS, less CSS, less surface area, more
//! reliable selection.
//!
//! The browser PDF viewer renders its own toolbar (page count,
//! zoom, download, print) which replaces the seamless dark-themed
//! card-body look. For a "Preprint" tab this is arguably more
//! appropriate UX — users immediately recognise it as a paper
//! and know how to navigate it.

use leptos::prelude::*;

#[component]
pub fn PreprintPanel() -> impl IntoView {
    view! {
        <div class="card">
            <div class="card-header d-flex justify-content-between align-items-center">
                <span>
                    "Preprint — ldsc-rs: a fast Rust reimplementation of LD Score Regression"
                </span>
                // Explicit download link with a friendlier filename
                // than the generic "preprint.pdf". The browser PDF
                // viewer's own toolbar also has a Save button, but
                // some mobile browsers (iOS Safari) hide the
                // viewer chrome and this gives users a guaranteed
                // path to the file.
                <a href="preprint.pdf"
                   download="ldsc-rs-preprint.pdf"
                   class="btn btn-sm btn-outline-primary"
                   title="Download the PDF (2 MB)">
                    "Download PDF"
                </a>
            </div>
            // `p-0` strips the card-body's default padding so the
            // iframe goes edge-to-edge inside the card.
            <div class="card-body p-0">
                // type="application/pdf" hints to the browser to
                // use the built-in PDF plugin even if it might
                // otherwise content-sniff. No sandbox attribute:
                // the browser's PDF viewer is trusted code (not
                // arbitrary HTML/JS); sandboxing it can break
                // PDF features (link navigation, page anchors)
                // on some browsers.
                <iframe src="preprint.pdf"
                        type="application/pdf"
                        class="preprint-iframe"
                        title="ldsc-rs preprint"
                        loading="lazy"></iframe>
            </div>
        </div>
    }
}
