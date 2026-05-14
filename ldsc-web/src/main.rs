//! ldsc-web — Leptos CSR frontend for the in-browser ldsc demo.
//!
//! Mount point is `<div id="app-root">` in `index.html`. The
//! placeholder text inside that div ("Loading ldsc-rs…") gets
//! replaced when Leptos hydrates / mounts on top.
//!
//! Visual design follows LDLink (CBIIT/nci-webtools-dceg-linkage):
//! a thin grey project banner (replaces the gov banner), a dark
//! `#4a4a4a` top nav with module tabs (gold underline on active),
//! a centered card column for content, and a sticky `#2a71a5`
//! footer. See the `assets/palette.css` overrides applied on top of
//! Bootstrap 5.

mod components;

use leptos::mount::mount_to;
use leptos::prelude::document;
use leptos::prelude::*;
use wasm_bindgen::JsCast;

use components::{Banner, Footer, L2Panel, NavBar};

/// Module currently shown in the main panel area. h2 / rg are stubbed
/// in v1 — the L2 panel is the headline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Module {
    L2,
    H2,
    Rg,
}

/// Root component. Lays out the LDLink-style chrome (banner / nav /
/// main / footer) and routes the active module to a panel.
#[component]
fn App() -> impl IntoView {
    let module = RwSignal::new(Module::L2);

    view! {
        <Banner />
        <NavBar module=module />
        <main>
            {move || match module.get() {
                Module::L2 => view! { <L2Panel /> }.into_any(),
                Module::H2 => view! {
                    <div class="card">
                        <div class="card-header">"H2 — Heritability regression"</div>
                        <div class="card-body">
                            <p class="text-muted mb-0">
                                "h² estimation from sumstats + LD scores is planned for v2.
                                The CLI fully supports it — see "
                                <code>"ldsc h2 --h2 …"</code>
                                ". For now, please use the CLI for h² / rg analyses."
                            </p>
                        </div>
                    </div>
                }.into_any(),
                Module::Rg => view! {
                    <div class="card">
                        <div class="card-header">"Rg — Genetic correlation"</div>
                        <div class="card-body">
                            <p class="text-muted mb-0">
                                "Genetic correlation across two GWAS sumstats is planned
                                for v2. See "
                                <code>"ldsc rg --rg … --ref-ld-chr … --w-ld-chr …"</code>
                                " in the CLI for now."
                            </p>
                        </div>
                    </div>
                }.into_any(),
            }}
        </main>
        <Footer />
    }
}

fn main() {
    // Funnel Rust panics through console.error with a stack trace
    // (otherwise they become a generic "Unreachable executed").
    console_error_panic_hook::set_once();

    // Lightweight tracing → console.log. Trust level configurable via
    // `?log=debug` later; default to info.
    tracing_wasm::set_as_global_default();

    let mount_target = document()
        .get_element_by_id("app-root")
        .expect("index.html must have <div id=\"app-root\">")
        .unchecked_into::<web_sys::HtmlElement>();

    // Clear the "Loading ldsc-rs…" placeholder before mount — Leptos
    // appends new children rather than replacing existing ones.
    mount_target.set_inner_html("");

    // Keep the mount alive for the lifetime of the page; the Leptos
    // 0.8 API hands back an UnmountHandle and complains if we drop
    // it implicitly.
    mount_to(mount_target, App).forget();
}
