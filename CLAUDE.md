# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Rust rewrite of [Bulik-Sullivan et al.'s LDSC](https://github.com/bulik/ldsc) (LD Score Regression), a tool for estimating heritability and genetic correlation from GWAS summary statistics. Six subcommands: `munge-sumstats`, `l2`, `h2`, `rg`, `make-annot`, `cts-annot`. The goal is numerical parity with the Python implementation at 18× speedup (23× with `--fast-f32`).

## Build & Test Commands

```bash
cargo build --release          # Always use release; debug is too slow
cargo check                    # Quick type-check
cargo clippy --release -- -D warnings   # Lint (CI enforces zero warnings)
cargo fmt --check              # Format check
cargo test --release           # Unit tests (no inline tests currently; this runs doc tests)
bash test_rust.sh --build      # Integration smoke tests (needs data/1000G_phase3_common_norel.{bed,bim,fam})
```

Production/benchmark build (musl + mimalloc, AVX2+FMA via `.cargo/config.toml`):
```bash
cargo build --release --features mimalloc --target x86_64-unknown-linux-musl
```
**Assembly analysis must use the musl target** — this is the production binary shipped to AWS HPC.
The native (glibc) build may produce different codegen. Always use `--target x86_64-unknown-linux-musl`
when inspecting assembly or benchmarking.

GPU build (optional, requires CUDA toolkit):
```bash
cargo build --release --features gpu
```

The release binary is at `target/release/ldsc` (native) or `target/x86_64-unknown-linux-musl/release/ldsc` (musl).

## CI

Defined in `.github/workflows/ci.yml`: runs `cargo check`, `cargo fmt --check`, `cargo clippy --release -- -D warnings`, and `cargo test --release` on every push/PR to main.

## Architecture

### Module Structure

All source lives in `src/`. The binary is a clap-derive CLI dispatcher (`main.rs` → `cli.rs`) that routes to subcommand modules:

- **`cli.rs`** — Pure clap derive structs (`MungeArgs`, `L2Args`, `H2Args`, `RgArgs`, `MakeAnnotArgs`, `CtsAnnotArgs`). No logic.
- **`main.rs`** — CLI dispatch + global thread pool setup. Injects implicit subcommand for backward-compatible `--h2`/`--l2` flag usage.
- **`munge.rs`** — Polars LazyFrame pipeline for GWAS summary statistics preprocessing. Streams input without loading entire file into RAM.
- **`l2.rs`** (~1800 lines, largest module) — LD score computation with ring-buffer genotype store and chunked GEMM. Sequential global pass across all chromosomes matching Python's cross-chromosome boundary behavior.
- **`regressions.rs`** (~2500 lines) — h2/rg regression drivers. Scalar and partitioned paths, two-step estimator, liability-scale conversion.
- **`h2.rs`** — IRWLS regression + block jackknife logic used by regressions.rs.
- **`irwls.rs`** / **`jackknife.rs`** — Iteratively Reweighted Least Squares and parallel block jackknife (rayon `par_iter` over 200 leave-one-out blocks).
- **`parse.rs`** — File I/O: `scan_tsv()`, `read_annot()`, `read_m_vec()`, BIM/FAM parsing, per-chromosome file concatenation.
- **`la.rs`** — Linear algebra type aliases and helpers wrapping `faer`. `MatF = Mat<f64>`, `ColF = Mat<f64>` (n×1).
- **`bed.rs`** — Custom PLINK BED reader (no external crate). Builder pattern: `Bed::builder(path).build()`, then `ReadOptions::builder().sid_index(vec).read(&bed)`.
- **`gpu.rs`** — Optional (`#[cfg(feature = "gpu")]`) CUDA matmul via CubeCL.
- **`make_annot.rs`** / **`cts_annot.rs`** — Annotation file generators.

### Key Data Flow (l2 subcommand)

1. Parse BIM → `Vec<BimRecord>`, apply `--extract` filter
2. Compute `block_left` window boundaries from coordinates
3. Sequential chunked pass: read genotype chunks from BED → normalize → ring-buffer GEMM (`B^T·B` within-chunk + `A^T·B` cross-window) → accumulate unbiased r² into L2 array
4. Group by chromosome, write per-chr `.l2.ldscore.gz` + `.l2.M` / `.l2.M_5_50`

### Key Data Flow (h2/rg subcommands)

1. Load sumstats + LD scores via Polars LazyFrame, inner-join on SNP
2. Detect annotation columns (K=1 scalar, K>1 partitioned)
3. IRWLS regression with block jackknife SE (200 blocks, rayon-parallelized)
4. Optional two-step: estimate intercept on low-chi² SNPs, fix and re-estimate

## Dependencies (key rationale)

- **`faer 0.24`** — All linear algebra (matmul, SVD). Replaced ndarray+OpenBLAS entirely. Parallelism via `Par::rayon(0)`.
- **`polars 0.53`** — Lazy CSV streaming for munge and LD score file loading.
- **`rayon 1`** — Parallel jackknife blocks and faer internal parallelism.
- **`cubecl 0.10.0-pre.2`/`cubek-matmul 0.10.0-pre.2`** - optional for GPU LD score calculation
- No bed-reader crate, no ndarray, no BLAS, no thiserror.

## Polars 0.53 Gotchas

- `LazyFrame::schema()` does not exist; use `lf.clone().limit(0).collect()` to inspect columns.
- `Series::new` requires `PlSmallStr` name: `"name".into()`.
- `DataFrame::with_column` takes `Column` not `Series`: `series.into()`.
- `Expr` does not implement `Not`; use `.not()` method.
- `get_column_names()` returns `Vec<PlSmallStr>`; use `.as_str()` to compare with `&str`.

## BED Reader API (src/bed.rs)

```rust
let bed = Bed::builder(path).build()?;
let mut ro = ReadOptions::builder();  // bind first to avoid temp borrow
ro.sid_index(vec_of_isize);           // indices are Vec<isize>
ro.iid_index(vec_of_isize);
let mat: Mat<f32> = ro.f32().read(&bed)?;
```

## Testing Against Python

Python LDSC reference: `uv run --project ldsc_py python ldsc_py/ldsc.py ...`

Parity/benchmark scripts in `scripts/`:
- `verify_parity_l2.sh` — LD score parity + perf on 200k extract
- `verify_parity_h2.sh` — Adversarial h2 correctness checks
- `verify_parity_munge.sh` — Munge output comparison
- `bench_l2.sh` — Quick l2 benchmark (5k SNPs default)
- `check_l2_tiny_py_vs_rust.sh` — Uses `--extract` on full 1.66M BED
- `aws-bench.sh` — Build, push to ECR, run on AWS Batch (c7a.4xlarge EPYC). **Use this for verified perf numbers; local timings are noisy.**

## Important Behavioral Notes

- `--extract` pre-filters BIM before LD computation (affects windows); `--print-snps` post-filters output only.
- `--maf` in l2 is a pre-filter (before LD computation), matching Python.
- `bed_idx` (original BIM row) differs from `pos` (index in filtered `all_snps`) when `--extract` is active.
- L2 uses a global sequential pass across all chromosomes (not per-chromosome parallel), matching Python's cross-chromosome window bleeding behavior.
- All CM columns in the 1000G reference data are 0 — never use `--ld-wind-cm` for benchmarking with this data.
- `make-annot` BED intervals are 0-based half-open; SNP BP=b matches [start,end) iff start < b <= end.

## Data Files

- Full genome: `data/1000G_phase3_common_norel.{bed,bim,fam}` (1,664,852 SNPs, 2,490 individuals)
- Benchmark subsets: `data/bench_5k.*`, `data/bench_200k.*`, plus chr22 variants

## Important Docs, Read When Needed
- `generated_docs/`: All dependencies have documentation locally saved here, prefer checking those over web searching for documentation. This has a lot of text so be careful when reading to not hammer your context window.
- `docs/ldsc_paper.md`/`docs/ldsc_supplementary_note.md`: Original LD Score paper and math.
- `docs/perf-log.md`: Active log of attempts to improve performance, should be edited after running a perf test.
- `docs/hpc.md`: Example of an HPC that `ldsc` is intended to run on, in this case UPenn.
- `ldsc.wiki/`: The wiki for the Bulik's LDSC implementation, has example workflows.
- `ldsc_py/`: Bulik's LDSC implementation.

## GPU

There's an experimental GPU option, code is available in `src/gpu.rs`, based on the CubeCL library for now using CUDA backend.
Documentation at `docs/gpu-feasibility.md`.

## Optimization Workflow
- **Verification**: Start with verifying that the optimization makes sense logically.
- **Implementation**: Then implement it in code.
- **Adversarial Pass**: Then run an adversarial pass over your code to ensure there are no bugs. If you see any, fix them and run the adversarial pass until there are none.
- **Local Parity/Perf**: Test with a smaller dataset parity with Python, and check for any obvious perf wins.
- **HPC Perf Test**: Run a proper perf benchmark with hyperfine on the AWS HPC.
- **Document**: Document all findings in `docs/perf-log.md`.
- **Commit**: Commit all changes with an informative commit message in Conventional Commits syntax.ß


## Release Process

```bash
cargo release patch            # dry run
cargo release patch --execute  # tag + publish
```
