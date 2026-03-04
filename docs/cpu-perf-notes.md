# CPU Performance Notes (No Algorithm Changes)

This note captures **potential** CPU-only performance wins that do **not** alter
algorithm design or numerical results. The items below are **not prioritized**;
we should only pursue them if profiling shows they matter for common workloads.

## 1) Jackknife Allocation Elimination
**Current behavior:** Each jackknife block builds `x_jack`, `y_jack`, `w_jack`
by allocating/copying slices and concatenating. This creates multiple full-size
allocations per block.

**Why it might help:**
- Large `n` and ~200 blocks can create heavy allocator and memory bandwidth
  pressure.
- Reusing per-thread buffers could reduce copies and heap churn.

**Why it might not:**
- Per-block IRWLS uses SVD (LAPACK) which can dominate runtime.
- README notes I/O + merge joins may dominate overall wall time.

**Expectation:** Potentially meaningful only if jackknife is a large fraction of
runtime. Otherwise likely a low single-digit % win.

## 2) LD Score Chunk Loop Allocations
**Current behavior:** Per chunk, allocates `bed_indices`, `raw`, and intermediate
matrices (`bb`, `ab`, `r2u_*`).

**Why it might help:**
- Many small chunks can create allocation overhead and cache churn.

**Why it might not:**
- `faer` matmul calls usually dominate for non-trivial sizes.
- `.bed` I/O can be more expensive than allocation.

**Expectation:** Likely a small (single-digit %) win unless chunk sizes are tiny
and compute is very light.

## 3) Matmul Thread Tuning (No Code Changes)
**Current behavior:** `faer` handles heavy GEMM. Rayon thread count is controlled
via CLI.

**Why it might help:**
- A heuristic for thread count based on `n_indiv` and window size can reduce
  oversubscription overhead.
- Small matrices can benefit from forcing sequential matmul.

**Expectation:** Can be material for some systems, but depends on hardware.

## Decision Gate
Do not implement the above without profiling that shows:
- The targeted area is a significant share of wall time.
- Expected speedup is meaningful for common workloads.

## HPC-Only / Whole-Dataset Window Ideas (Potentially High Impact)
These are **not suitable to validate on a weak dev CPU**; they should be tested
on the HPC where bandwidth and core counts match production.

### A) GEMM Backend Swap for `--yes-really`
**Idea:** Route `ab_dot`/`bb_dot` matmul to an alternate high-performance GEMM
backend or faer’s fastest kernels when `--yes-really` is used.

**Why it might help:** Whole-dataset windows are dominated by GEMM; better kernels
can yield large wins.

**Parity risk:** Low (same math, different FP ordering). Expect tiny float diffs only.

### B) Block/Tiling Parallelism for Whole-Dataset Windows
**Idea:** Parallelize across large correlation tiles (explicit blocking) instead of
small chunk windows, using bigger tiles for better cache reuse.

**Why it might help:** Reduces overhead and improves cache/bandwidth utilization.

**Parity risk:** Low, but reordering can change floating-point summation order slightly.

### C) Aggressive Chunk/Ring Tuning for `--yes-really`
**Idea:** Increase chunk sizes and ring size only when `--yes-really` is set.

**Why it might help:** Whole-dataset windows benefit from fewer, larger GEMMs.

**Parity risk:** None (same math).

### D) Pre-scale/Normalize Packed Buffers
**Idea:** Precompute normalized genotype blocks into contiguous buffers to avoid
repeated scaling inside inner loops.

**Why it might help:** Reduces memory traffic inside GEMM hot path.

**Parity risk:** None if done deterministically.

## Optional Fast-F32 Path (Not Parity-Safe)
**Idea:** Provide a compile-time `fast-f32` feature for `l2` to compute in f32
and cast outputs to f64 for writing.

**Why it might help:** Faster GEMM on some CPUs, lower memory bandwidth.

**Parity risk:** Medium/high. Results will differ vs f64; acceptable only for a
separate “fast” mode. Observed full-1000G output deltas:
`mean_abs_diff=3.03e-4`, `rmse=5.79e-4`, `max_abs_diff=0.008`,
`max_rel_diff=6.51e-4` (relative to f64).

## Profiling Guidance
Measure:
- Total time in `jackknife()` vs `ldscore()`.
- Within jackknife: time in IRWLS vs time copying/allocating arrays.
- Within `ldscore`: matmul time vs allocation + bed I/O.
