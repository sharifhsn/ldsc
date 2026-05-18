//! `<input type="file">` row with the LDLink-styled chrome plus a
//! "loaded" indicator. Reads the file into a `Vec<u8>` (BED) or
//! `String` (BIM/FAM) via gloo-file's Promise-based readers.
//!
//! The loaded bytes are shared back to the parent via a callback —
//! the parent owns the `RwSignal<Option<Vec<u8>>>` (or `String`).

use gloo_file::{Blob, futures::read_as_bytes, futures::read_as_text};
use leptos::ev::Event;
use leptos::prelude::document;
use leptos::prelude::*;
use leptos::task::spawn_local;
use wasm_bindgen::JsCast;
use web_sys::{HtmlInputElement, MouseEvent};

/// What kind of payload to read out of the picked file.
#[derive(Clone, Copy)]
pub enum ReadAs {
    /// Read as raw bytes — small files where it's cheap. Currently
    /// unused; the BED uses `BlobHandle` instead so that multi-GB
    /// files don't have to live in wasm linear memory.
    #[allow(dead_code)]
    Bytes,
    /// Read as UTF-8 text — for `.bim` (~44 MB at full-genome 1000G)
    /// and `.fam` (~9 KB-1.3 MB at biobank).
    Text,
    /// Don't read at all — just hold the `web_sys::File` handle so
    /// the eventual consumer (Web Worker) can stream chunks itself
    /// via `FileReaderSync`. This is what the `.bed` row uses; it
    /// makes the upload essentially free regardless of file size.
    BlobHandle,
}

/// One uploaded file's loaded state. The size is tracked separately
/// from the payload so the UI can show "200 MB" without re-measuring.
#[derive(Clone, Default)]
pub struct LoadedFile {
    pub name: String,
    pub size_bytes: usize,
    pub bytes: Option<Vec<u8>>,
    pub text: Option<String>,
    /// Raw `web_sys::File` handle for streaming — set by
    /// [`ReadAs::BlobHandle`]. Wrapped in `SendWrapper` because
    /// `web_sys::File` is `!Send` but Leptos's `RwSignal` storage
    /// internally tracks Send-ness.
    pub file_handle: Option<send_wrapper::SendWrapper<web_sys::File>>,
}

impl LoadedFile {
    pub fn is_loaded(&self) -> bool {
        self.bytes.is_some() || self.text.is_some() || self.file_handle.is_some()
    }
}

#[component]
pub fn FileUploadRow(
    /// Visible extension label shown in the row chrome (e.g. `".bed"`).
    /// Used as a default for the picker's `accept` attribute when
    /// `accept` is not supplied.
    ext: &'static str,
    read_as: ReadAs,
    /// Parent-owned signal that receives the loaded file.
    file: RwSignal<LoadedFile>,
    /// Optional comma-separated `accept` value for the underlying
    /// `<input type="file">`. Use to broaden what the OS picker shows
    /// vs the displayed label — e.g. `ext=".sumstats"` (label) with
    /// `accept=".sumstats,.sumstats.gz"` (picker accepts both). When
    /// `None`, the picker accepts exactly `ext`.
    #[prop(optional)]
    accept: Option<&'static str>,
) -> impl IntoView {
    let accept_attr: &'static str = accept.unwrap_or(ext);
    let input_id = format!(
        "file-upload-{}",
        ext.trim_start_matches('.').split(',').next().unwrap_or(ext),
    );
    let input_id_for_click = input_id.clone();
    let input_id_for_label = input_id.clone();

    let on_change = move |ev: Event| {
        let input: HtmlInputElement = ev
            .target()
            .and_then(|t| t.dyn_into::<HtmlInputElement>().ok())
            .expect("file picker change event target");
        let Some(file_list) = input.files() else {
            return;
        };
        let Some(picked) = file_list.get(0) else {
            return;
        };

        let name = picked.name();
        let size_bytes = picked.size() as usize;

        match read_as {
            ReadAs::Bytes => {
                let blob: Blob = picked.into();
                spawn_local(async move {
                    match read_as_bytes(&blob).await {
                        Ok(bytes) => file.set(LoadedFile {
                            name,
                            size_bytes,
                            bytes: Some(bytes),
                            text: None,
                            file_handle: None,
                        }),
                        Err(e) => tracing::error!("read_as_bytes failed: {e:?}"),
                    }
                });
            }
            ReadAs::Text => {
                let blob: Blob = picked.into();
                spawn_local(async move {
                    match read_as_text(&blob).await {
                        Ok(text) => file.set(LoadedFile {
                            name,
                            size_bytes,
                            bytes: None,
                            text: Some(text),
                            file_handle: None,
                        }),
                        Err(e) => tracing::error!("read_as_text failed: {e:?}"),
                    }
                });
            }
            ReadAs::BlobHandle => {
                // Just keep the File handle. No bytes touched here —
                // the consumer (Web Worker) will stream chunks via
                // FileReaderSync as compute progresses.
                file.set(LoadedFile {
                    name,
                    size_bytes,
                    bytes: None,
                    text: None,
                    file_handle: Some(send_wrapper::SendWrapper::new(picked)),
                });
            }
        }
    };

    let row_class = move || {
        if file.with(|f| f.is_loaded()) {
            "file-row loaded"
        } else {
            "file-row"
        }
    };

    let filename_view = move || {
        file.with(|f| {
            if f.is_loaded() {
                view! { <span class="filename">{f.name.clone()}</span> }.into_any()
            } else {
                view! {
                    <span class="filename placeholder">"Click to choose a file…"</span>
                }
                .into_any()
            }
        })
    };

    let size_view = move || {
        file.with(|f| {
            if f.is_loaded() {
                view! { <span class="size">{format_bytes(f.size_bytes)}</span> }.into_any()
            } else {
                view! { <span class="size">""</span> }.into_any()
            }
        })
    };

    let check_view = move || {
        file.with(|f| {
            if f.is_loaded() {
                view! { <span class="check">"✓ loaded"</span> }.into_any()
            } else {
                view! { <span class="check" style="visibility:hidden">"-"</span> }.into_any()
            }
        })
    };

    // Make the entire row clickable as a label proxy — the actual
    // <input> is hidden but still functional.
    let on_row_click = move |_ev: MouseEvent| {
        if let Some(el) = document().get_element_by_id(&input_id_for_click)
            && let Ok(input) = el.dyn_into::<HtmlInputElement>()
        {
            input.click();
        }
    };

    view! {
        <div class=row_class on:click=on_row_click style="cursor: pointer;">
            <span class="ext">{ext}</span>
            {filename_view}
            {size_view}
            {check_view}
            <label for=input_id_for_label.clone() style="display: none;"></label>
            <input
                id=input_id
                type="file"
                accept=accept_attr
                on:change=on_change
                on:click=|ev: MouseEvent| ev.stop_propagation()
                style="display: none;"
            />
        </div>
    }
}

fn format_bytes(n: usize) -> String {
    if n < 1024 {
        format!("{n} B")
    } else if n < 1024 * 1024 {
        format!("{:.1} KB", n as f64 / 1024.0)
    } else if n < 1024 * 1024 * 1024 {
        format!("{:.1} MB", n as f64 / (1024.0 * 1024.0))
    } else {
        format!("{:.2} GB", n as f64 / (1024.0 * 1024.0 * 1024.0))
    }
}
