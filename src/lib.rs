//! ldsc — Rust reimplementation of LD Score Regression (Bulik-Sullivan et al.).
//!
//! On `wasm32-unknown-unknown` with `+atomics`, we use the unstable
//! `core::arch::wasm32::memory_atomic_wait32` / `memory_atomic_notify`
//! intrinsics to implement a manual SAB-backed worker pool (see
//! `wasm_simd::pool`). The intrinsics are gated behind the
//! `stdarch_wasm_atomic_wait` feature, opted into here.
#![cfg_attr(
    all(target_arch = "wasm32", target_feature = "atomics"),
    feature(stdarch_wasm_atomic_wait)
)]

//!
//! This crate is consumed in three ways:
//!
//! 1. **As a binary** — `cargo install ldsc` produces the `ldsc` CLI built
//!    from `src/main.rs`. That binary is a thin shim around the entry
//!    points re-exported here (`run_munge`, `run_l2`, `run_h2`, `run_rg`,
//!    `run_make_annot`, `run_cts_annot`).
//!
//! 2. **As a Rust library** — downstream Rust callers can take a
//!    fully-parsed [`cli::L2Args`] / [`cli::H2Args`] / etc. and run a
//!    subcommand programmatically by calling the same `run_*` functions
//!    used by the CLI.
//!
//! 3. **As a WASM library** — the `ldsc-web` Leptos frontend depends on
//!    this crate built for `wasm32-unknown-unknown`. The browser path
//!    bypasses the disk-backed `run_*` orchestrators and calls the
//!    lower-level building blocks ([`bed`], [`parse`], [`l2::compute`],
//!    [`regressions`], [`munge`]) against in-memory inputs. The
//!    in-memory loaders (`from_bytes` / `from_str` variants) live next
//!    to the disk-backed loaders in those modules.
//!
//! All modules are re-exported here. Anything marked `pub(crate)` inside
//! a module is intentionally not part of the public API and may change
//! between releases.

pub mod bed;
pub mod cli;
pub mod cts_annot;
pub mod frame;
#[cfg(feature = "gpu")]
pub mod gpu;
pub mod h2;
pub mod irwls;
pub mod jackknife;
pub mod l2;
pub mod la;
pub mod make_annot;
pub mod munge;
pub mod parse;
pub mod regressions;
pub mod wasm_simd;

// Convenience re-exports so callers can `use ldsc::run_l2;` rather than
// reaching into submodules.
pub use cts_annot::run as run_cts_annot;
pub use l2::run as run_l2;
pub use make_annot::run as run_make_annot;
pub use munge::run as run_munge;
pub use regressions::{run_h2, run_rg};
