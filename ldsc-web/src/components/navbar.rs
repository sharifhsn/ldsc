//! Top navigation bar — LDLink's dark `#4a4a4a` module switcher
//! (active tab gets a 4-px gold underline).
//!
//! Module routing is single-page; the parent passes an `RwSignal<Module>`
//! that we mutate on click.

use leptos::ev::MouseEvent;
use leptos::prelude::*;

use crate::Module;

#[component]
pub fn NavBar(module: RwSignal<Module>) -> impl IntoView {
    let tab = move |label: &'static str, m: Module, enabled: bool| -> AnyView {
        let active = move || module.get() == m;
        let on_click = move |ev: MouseEvent| {
            ev.prevent_default();
            if enabled {
                module.set(m);
            }
        };

        view! {
            <a href="#"
               class="nav-tab"
               class:active=active
               class:disabled=move || !enabled
               on:click=on_click>
                {label}
            </a>
        }
        .into_any()
    };

    view! {
        <nav class="app-nav">
            <div class="d-flex justify-content-between align-items-center"
                 style="max-width: 1200px; margin: 0 auto; padding: 0 1.25rem;">
                <a href="#" class="wordmark" on:click=|ev| ev.prevent_default()>
                    "ldsc-rs"
                    <span class="v">" · v0.5.0 in WASM"</span>
                </a>
                <div>
                    {tab("L2", Module::L2, true)}
                    {tab("H2", Module::H2, true)}
                    {tab("Rg", Module::Rg, true)}
                    {tab("Preprint", Module::Preprint, true)}
                </div>
            </div>
        </nav>
    }
}
