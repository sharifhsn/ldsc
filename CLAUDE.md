# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Rust rewrite of [Bulik-Sullivan et al.'s LDSC](https://github.com/bulik/ldsc) (LD Score Regression), a tool for estimating heritability and genetic correlation from GWAS summary statistics. Six subcommands: `munge-sumstats`, `l2`, `h2`, `rg`, `make-annot`, `cts-annot`. Numerically exact parity with Python (`max_abs_diff=0` in f64 mode). Current performance vs Python on AWS EPYC 7R13:

- **1000G (N=2,490)**: 37.7× exact, 101× with `--sketch 50`
- **Biobank (N=50K)**: 20× with `--sketch 50` (best), 1.84× exact-f32

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
- **`l2/`** — LD score computation, split into submodules:
  - **`mod.rs`** (~570 lines) — `run()` orchestrator: arg validation, BIM/FAM loading, `--extract`/`--annot`/`--keep` wiring, output writing.
  - **`compute.rs`** (~1500 lines) — `compute_ldscore_global`: ring-buffer GEMM loop (sequential global pass, scalar + partitioned + sketch paths). Also `GemmBufs`, GPU compat helpers.
  - **`window.rs`** — `WindowMode` enum + `get_block_lefts_*`: per-chromosome window boundary computation.
  - **`normalize.rs`** — `normalize_col_f{32,64}_with_stats`: impute NaN→mean, centre, unit-variance. AVX2+FMA `sum_sumsq_f32`.
  - **`snp_stats.rs`** — `compute_snp_stats`: fast BED scan for MAF + het/missing prefilter.
  - **`io.rs`** — FAM parsers (`count_fam`, `parse_fam`), SNP-set loaders, gzip TSV writers.
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
- `verify_parity_all.sh` — Comprehensive 8-mode parity sweep (exact/sketch × f64/f32 vs Python, chr22 data)
- `aws-bench.sh` — Build, push to ECR, run on AWS Batch (c6a.4xlarge EPYC). **Use this for verified perf numbers; local timings are noisy.**
- `biobank-bench.sh` — AWS Batch biobank 50K benchmark (all 8 mode×precision variants)
- `make_synthetic_biobank.py` — Generate synthetic N=50K BED from 1000G (~21× replication of N=2,490, 1% noise)

## Important Behavioral Notes

- `--extract` pre-filters BIM before LD computation (affects windows); `--print-snps` post-filters output only.
- `--maf` in l2 is a pre-filter (before LD computation), matching Python.
- `bed_idx` (original BIM row) differs from `pos` (index in filtered `all_snps`) when `--extract` is active.
- L2 uses a global sequential pass across all chromosomes (not per-chromosome parallel), matching Python's cross-chromosome window bleeding behavior.
- All CM columns in the 1000G reference data are 0 — never use `--ld-wind-cm` for benchmarking with this data.
- `make-annot` BED intervals are 0-based half-open; SNP BP=b matches [start,end) iff start < b <= end.
- `--sketch` projection matrix P uses deterministic seed 42 (`fastrand::Rng::with_seed(42)`). Same `--sketch d` on same data always produces identical output.

## Data Files

- Full genome: `data/1000G_phase3_common_norel.{bed,bim,fam}` (1,664,852 SNPs, 2,490 individuals)
- Biobank synthetic: `data/biobank_50k.{bed,bim,fam}` (1,664,852 SNPs, 50,000 individuals, ~20GB BED)
- Benchmark subsets: `data/bench_5k.*`, `data/bench_200k.*`, plus chr22 variants
- S3 mirror: `s3://ldsc-bench-data-270497617191/` (all datasets including biobank_50k)

## Important Docs, Read When Needed
- `generated_docs/`: All dependencies have documentation locally saved here, prefer checking those over web searching for documentation. This has a lot of text so be careful when reading to not hammer your context window.
- `docs/ldsc_paper.md`/`docs/ldsc_supplementary_note.md`: Original LD Score paper and math.
- `docs/perf-log.md`: Active log of attempts to improve performance, should be edited after running a perf test.
- `docs/hpc.md`: Example of an HPC that `ldsc` is intended to run on, in this case UPenn.
- `ldsc.wiki/`: The wiki for the Bulik's LDSC implementation, has example workflows.
- `ldsc_py/`: Bulik's LDSC implementation.

## Approximate / Fast Modes

Three opt-in modes trade precision for speed in the `l2` subcommand:

- **`--fast-f32`** — Use f32 for genotype storage and GEMM. ~1.85× speedup (bandwidth-bound). Numerically close to f64; safe for downstream h2/rg. Combinable with all other modes.

- **`--sketch <d>`** — Randomized dimensionality reduction via CountSketch before GEMM. Fused BED-decode-normalize-scatter-add kernel, O(N×c) regardless of d. Cost is flat in d until d≈√N, so d=200 has same speed as d=50 but much better accuracy.
  - Accuracy depends on d: d=50 r≈0.81, d=100 r≈0.85, d=200 r≈0.93, d=500 r≈0.97
  - **Automatically enables f32** — CountSketch ±1 entries are exactly representable in f32. `--fast-f32` is redundant with `--sketch`.
  - Deterministic (seed=42). Same `--sketch d` on same data always produces identical output.


### Performance Summary (AWS EPYC 7R13 c6a.4xlarge, 16 vCPU, 1.66M SNPs)

**1000G reference (N=2,490)** — Python baseline: 1548.5s

| Mode | Time | vs Python |
|------|------|-----------|
| exact | 41.1s | 37.7× |
| --sketch 50 | 15.3s | 101× |
| --sketch 200 | 25.4s | 61× |

**Biobank scale (N=50,000)** — Python baseline: ~1548s (extrapolated)

| Mode | Time | vs exact-f64 |
|------|------|-------------|
| exact-f64 | 665.9s | 1.0× |
| exact-f32 | 361.9s | 1.84× |
| **--sketch 50** | **33.1s** | **20.1×** |
| **--sketch 200** | **33.8s** | **19.7×** |
| --sketch 500 | 36.2s | 18.4× |
| --sketch 1000 | 39.7s | 16.8× |
| subsample-5k+sketch50 | 24.9s | 26.7× |

Key insight: Fused CountSketch eliminates the N×c intermediate buffer entirely — it reads packed BED bytes and scatter-adds directly into the d×c sketch buffer. Cost is O(N×c) independent of d, so d=200 is "free" accuracy vs d=50.

## GPU

There's an experimental GPU option, code is available in `src/gpu.rs`, based on the CubeCL library for now using CUDA backend.
Documentation at `docs/gpu-feasibility.md`.

### CUDA Version Compatibility

`cudarc` (CubeCL's CUDA backend) eagerly loads driver symbols at runtime matching the **compile-time** CUDA version. If the binary is built for CUDA 13.1 but runs on a system with an older driver, it panics on missing symbols.

**Control the target CUDA version at build time:**
```bash
CUDARC_CUDA_VERSION=12040 cargo build --release --features gpu   # CUDA 12.4 (AWS default)
CUDARC_CUDA_VERSION=11080 cargo build --release --features gpu   # CUDA 11.8 (older HPC)
```

Without `CUDARC_CUDA_VERSION`, cudarc auto-detects via `nvcc --version` (good for local dev). In Docker builds without `nvcc`, it falls back to CUDA 13.1 (latest) which will crash on older drivers.

**Target versions by environment:**
| Target | `CUDARC_CUDA_VERSION` | Driver | Notes |
|--------|----------------------|--------|-------|
| Local dev (RTX 2070) | auto (omit) | 590 | nvcc auto-detects |
| AWS g5/g6 (AL2023 AMI) | `12040` | 580 | Safe default |
| AWS g5/g6 (AL2 AMI) | `12040` | 550 | Safe default |
| PMACS P100 | `11080` | 418-470 | No tensor cores |

Check target CUDA version: `nvidia-smi` (top-right shows "CUDA Version").

## Optimization Workflow
- **Verification**: Start with verifying that the optimization makes sense logically.
- **Implementation**: Then implement it in code.
- **Adversarial Pass**: Then run an adversarial pass over your code to ensure there are no bugs. If you see any, fix them and run the adversarial pass until there are none.
- **Local Perf Gate**: Run `bash scripts/local-perf-gate.sh --full` and compare against known baselines below. **Do NOT submit AWS jobs until local numbers look sane.** For N-dependent changes (anything touching GEMM parallelism, tiling, or buffer allocation), also test with `--biobank` if data is available.
- **HPC Perf Test**: Run a proper perf benchmark with hyperfine on the AWS HPC. Always confirm with the user before launching — each biobank run costs ~$1-2 and takes 2+ hours.
- **Document**: Document all findings in `docs/perf-log.md`.
- **Commit**: Commit all changes with an informative commit message in Conventional Commits syntax.

### Local Perf Baselines (bench_200k, N=2490, `--ld-wind-kb 1000`)

Reference timings (Ryzen 5 5600X, 2026-03-14):

| Mode | Time | vs exact-f64 |
|------|------|-------------|
| exact-f64 | 13.6s | 1.0× |
| exact-f32 | 12.0s | 1.13× |
| sketch-50 (auto-f32) | 3.1s | 4.4× |
| sketch-100 | 2.7s | 5.1× |
| sketch-200 | 3.3s | 4.2× |
| countsketch-50 | 1.7s | 8.0× |
| countsketch-100 | 3.6s | 3.8× |

**Red flags**: any mode >30% slower than these baselines, or sketch slower than exact.

### Perf Regression Lessons Learned

- **Fused CountSketch is O(N×c), not O(d×N×c).** Its scatter-add kernel bypasses GEMM entirely; d=200 has the same cost as d=50. A Rademacher-style fused tile kernel (512-individual tiles, many small `Par::Seq` GEMMs) caused a 2.2× regression (285s vs 130s) — avoid fragmenting the GEMM.
- **Parallelism is N-dependent.** A change that's neutral at N=2,490 can be catastrophic at N=50,000. Always think about how parallelism scales with N before submitting HPC jobs.
- **Diff the hot path against the last validated binary before any AWS run.** The sketch projection path in `compute.rs` (search for "Sketch projection") must go through the CountSketch scatter-add path, not any GEMM.


## Release Process

```bash
cargo release patch            # dry run
cargo release patch --execute  # tag + publish
```
