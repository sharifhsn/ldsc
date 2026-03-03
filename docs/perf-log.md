# Performance Log

This log captures **post-change** performance measurements and parity checks.
Each entry should list the dataset, command, and key timings so we can track
how changes affect runtime.

## 2026-03-03
- Change: low-risk l2 optimizations (buffer reuse, contiguous BED indexing, `general_mat_mul`).
- Dataset: 1000G.EUR.QC, extract 500k SNPs (`/Users/sharif/Code/ldsc/perf/l2/extract_500000.snps`).
- Command: `ldsc l2 --bfile /Users/sharif/Code/ldsc/data/1000G.EUR.QC --out /Users/sharif/Code/ldsc/perf/l2/rust_l2_after --ld-wind-cm 1 --extract /Users/sharif/Code/ldsc/perf/l2/extract_500000.snps`
- Rust timing: `real 22.39s` (from `/Users/sharif/Code/ldsc/perf/l2/rust_l2_after.time`).
- Previous baseline (pre-opt): `real 27.19s` (from `/Users/sharif/Code/ldsc/perf/l2/rust_l2.time`).
- Speedup: ~17.7%.
- Parity check (Python vs Rust, 500k extract): **FAILED**.
  - Script: `EXTRACT_N=500000 /Users/sharif/Code/ldsc/scripts/verify_parity_l2.sh`
  - Result: `max_abs_diff=0.001` at `rs4275489` (tolerance 1e-3 exceeded).

## 2026-03-03 (later)
- Change: reverted low-risk l2 optimizations (removed buffer reuse, contiguous BED index fast path, and `general_mat_mul` usage).
- Parity check (Python vs Rust, 500k extract): **FAILED**.
  - Script: `EXTRACT_N=500000 /Users/sharif/Code/ldsc/scripts/verify_parity_l2.sh`
  - Result: `max_abs_diff=0.001` at `rs4275489` (tolerance 1e-3 exceeded).
  - Implication: parity drift persists without those optimizations.

## 2026-03-03 (later, option D)
- Change: added f64 default normalization with opt-in `--fast-f32`.
- Dataset: 1000G.EUR.QC (full dataset, no extract).
- Command (f64 default): `ldsc l2 --bfile /Users/sharif/Code/ldsc/data/1000G.EUR.QC --out /Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64 --ld-wind-cm 1`
- Rust timing (f64): `real 298.17s` (from `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64.time`).
- Command (fast f32): `ldsc l2 --bfile /Users/sharif/Code/ldsc/data/1000G.EUR.QC --out /Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f32 --ld-wind-cm 1 --fast-f32`
- Rust timing (fast f32): `real 329.21s` (from `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f32.time`).
- Result: `--fast-f32` is ~10.4% slower on this run.

## 2026-03-03 (baseline: f64-only)
- Change: removed `--fast-f32`; L2 normalization is f64-only.
- Dataset: 1000G.EUR.QC (full dataset, no extract).
- Command: `ldsc l2 --bfile /Users/sharif/Code/ldsc/data/1000G.EUR.QC --out /Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace --ld-wind-cm 1`
- Trace log: `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.stdout`
- Timing (from trace): `maf_prefilter=52.005s`, `compute_ldscore=215.833s` with
  `bed_read=54.602s`, `norm=11.225s`, `bb_dot=9.473s`, `ab_dot=93.586s`, `r2u=17.356s`,
  and `write_outputs=6.830s`.
- Wall time: `real 276.58s` (from `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.stderr`).
- Output fingerprints: `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.sha256` (69 files).

## 2026-03-03 (ab_dot optimizations)
- Change: preallocated `ab` GEMM buffer, contiguous ring-buffer fast path for `A` windows, and precomputed pq weights.
- Dataset: 1000G.EUR.QC (full dataset, no extract).
- Command: `ldsc l2 --bfile /Users/sharif/Code/ldsc/data/1000G.EUR.QC --out /Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace --ld-wind-cm 1`
- Trace log: `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.stdout`
- Timing (from trace): `maf_prefilter=49.605s`, `compute_ldscore=207.446s` with
  `bed_read=54.576s`, `norm=11.607s`, `bb_dot=9.268s`, `ab_dot=90.946s`, `r2u=34.483s`,
  and `write_outputs=7.102s`.
- Wall time: `real 267.37s` (from `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.stderr`).
- Output parity: `shasum -a 256 -c /Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.sha256` → OK.
