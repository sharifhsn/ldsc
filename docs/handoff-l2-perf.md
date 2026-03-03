# Handoff: L2 Perf/Parity Work (2026-03-03)

## Location
- Repo: `/Users/sharif/Code/ldsc`

## Current Goal
- Aggressively optimize `l2` while **preserving exact parity**.
- Parity check rule: **full 1000G run** and **SHA256 match** against baseline outputs.

## Perf Directory Minimal Keep List (enforced)
Only keep the following files:
- `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.sha256`
- `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.stdout`
- `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.stderr`

All other `perf/` subdirectories and files are removed.

## Baseline (f64-only, parity-correct)

**Baseline outputs (hash manifest)**
- `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.sha256`

**Baseline trace run**
```
ldsc l2 --bfile /Users/sharif/Code/ldsc/data/1000G.EUR.QC \
  --out /Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace \
  --ld-wind-cm 1
```

**Baseline trace/timing**
- Trace log: `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.stdout`
- Wall time: `real 276.58s` (from `.stderr`)
- Breakdown:
  - `maf_prefilter=52.005s`
  - `compute_ldscore=215.833s`
    - `bed_read=54.602s`
    - `norm=11.225s`
    - `bb_dot=9.473s`
    - `ab_dot=93.586s`
    - `r2u=17.356s`
  - `write_outputs=6.830s`

## Latest changes (ab_dot optimizations) — parity preserved

**Changes implemented**
- Preallocated `ab` GEMM buffer (`ab_buf`) and reused it.
- Contiguous ring-buffer fast path (skip `A` copy when window slots are contiguous).
- Precomputed `pq` weights per chunk/window to avoid `powf` in inner loops.
- Added `general_mat_mul` for explicit GEMM into preallocated buffer.

**Files changed**
- `/Users/sharif/Code/ldsc/src/l2.rs`
- `/Users/sharif/Code/ldsc/docs/perf-log.md`

**Perf after changes**
- Trace log: `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.stdout`
- Wall time: `real 267.37s`
- Breakdown:
  - `maf_prefilter=49.605s`
  - `compute_ldscore=207.446s`
    - `bed_read=54.576s`
    - `norm=11.607s`
    - `bb_dot=9.268s`
    - `ab_dot=90.946s`
    - `r2u=34.483s` (increased)
  - `write_outputs=7.102s`

**Parity validation (new rule)**
```
cd /Users/sharif/Code/ldsc/perf/l2
shasum -a 256 -c rust_l2_full_f64_trace.sha256
```
- Result: OK  
- Log: `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.sha256.check`

**Perf log updated**
- `/Users/sharif/Code/ldsc/docs/perf-log.md`

## Removed / Reverted
- Removed `--fast-f32` (f32 path). L2 now always uses f64 normalization (matches Python).

## How to Continue (workstation)

1. Pull latest repo state.
2. Run full 1000G trace:
```
/usr/bin/time -p env RUST_LOG=ldsc=trace ldsc l2 \
  --bfile /Users/sharif/Code/ldsc/data/1000G.EUR.QC \
  --out /Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace \
  --ld-wind-cm 1 \
  > /Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.stdout \
  2> /Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.stderr
```
3. Verify parity via hash:
```
cd /Users/sharif/Code/ldsc/perf/l2
shasum -a 256 -c rust_l2_full_f64_trace.sha256
```
4. Record any new changes and timings in `/Users/sharif/Code/ldsc/docs/perf-log.md`.

## Adversarial Workflow (Parity-Safe Perf Iteration)
Use this every time you change performance-sensitive code:

1. **Start clean**: ensure only the minimal keep list exists under `perf/`.
   - `ls /Users/sharif/Code/ldsc/perf/l2` should show exactly the three keep files.
2. **Run full 1000G trace** (Rust only, no Python):
```
/usr/bin/time -p env RUST_LOG=ldsc=trace ldsc l2 \
  --bfile /Users/sharif/Code/ldsc/data/1000G.EUR.QC \
  --out /Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace \
  --ld-wind-cm 1 \
  > /Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.stdout \
  2> /Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.stderr
```
3. **Parity check**: compare outputs against the baseline hash manifest.
```
cd /Users/sharif/Code/ldsc/perf/l2
shasum -a 256 -c rust_l2_full_f64_trace.sha256
```
   - If any file mismatches: treat as parity regression, revert or investigate before proceeding.
4. **Update perf log**: record new timing breakdowns and wall time in
   `/Users/sharif/Code/ldsc/docs/perf-log.md`.

Notes:
- The hash manifest is the source of truth for parity (not decimal rounding).
- If you need to regenerate the manifest (rare), do it once and document why.

## Known Hotspots
- `ab_dot` remains dominant (~90–94s).
- `r2u` time doubled after last change; investigate cache/loop ordering.

## Next Ideas (not implemented)
- Optimize `r2u` loop (vectorized/in-place transform to reduce passes).
- Explore ring buffer layout to avoid extra copies and minimize `r2u` cache misses.
- Consider multi-threading only in `r2u` if BLAS is single-threaded.
