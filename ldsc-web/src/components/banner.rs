//! Project banner — the ldsc-web equivalent of LDLink's "official
//! website of the United States government" strip. Same visual slot
//! (thin grey bar, full-width, sits above the dark nav), but shows
//! project links instead of gov chrome since we're not NIH.

use leptos::prelude::*;

#[component]
pub fn Banner() -> impl IntoView {
    view! {
        <div class="project-banner">
            <div class="d-flex justify-content-between align-items-center"
                 style="max-width: 1200px; margin: 0 auto;">
                <span>
                    <strong>"Open source"</strong>
                    " · "
                    <a href="https://github.com/sharifhsn/ldsc"
                       target="_blank"
                       rel="noopener noreferrer">
                        "sharifhsn/ldsc on GitHub"
                    </a>
                    " · "
                    <a href="https://github.com/sharifhsn/ldsc#readme"
                       target="_blank"
                       rel="noopener noreferrer">"README"</a>
                </span>
                <span style="color: #777; font-size: 0.78rem;">
                    "Pure-Rust WASM port. No data leaves your browser."
                </span>
            </div>
        </div>
    }
}
