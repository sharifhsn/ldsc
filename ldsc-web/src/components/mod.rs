//! Leptos components for the ldsc-web UI.
//!
//! Each component owns one visual region of the LDLink-inspired
//! chrome (banner / nav / footer) or one application surface (L2
//! panel + the helper file_upload / charts). Components communicate
//! via `RwSignal`s passed by the root `App` component in `main.rs`.

mod banner;
mod charts;
mod file_upload;
mod footer;
mod l2_panel;
mod navbar;
mod worker_client;

pub use banner::Banner;
pub use footer::Footer;
pub use l2_panel::L2Panel;
pub use navbar::NavBar;
