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
- BLAS `dot` calls usually dominate for non-trivial sizes.
- `.bed` I/O can be more expensive than allocation.

**Expectation:** Likely a small (single-digit %) win unless chunk sizes are tiny
and compute is very light.

## 3) BLAS Backend and Thread Tuning (No Code Changes)
**Current behavior:** BLAS is used for heavy GEMM. OpenBLAS thread count is set
from CLI.

**Why it might help:**
- System-tuned BLAS on HPC can be faster than a generic static build.
- A heuristic for thread count based on `n_indiv` and window size can reduce
  oversubscription overhead.

**Expectation:** Can be material for some systems, but depends on hardware.

## Decision Gate
Do not implement the above without profiling that shows:
- The targeted area is a significant share of wall time.
- Expected speedup is meaningful for common workloads.

## Profiling Guidance
Measure:
- Total time in `jackknife()` vs `ldscore()`.
- Within jackknife: time in IRWLS vs time copying/allocating arrays.
- Within `ldscore`: BLAS time vs allocation + bed I/O.
