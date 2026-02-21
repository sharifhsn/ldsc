# GPU Feasibility Note (Future Plan)

## Context
This repository already achieves major speedups on CPU by combining:
- A ring-buffer genotype store with a single DGEMM per chunk (`AᵀB`, `BᵀB`).
- Tuned BLAS thread counts to avoid oversubscription on small windows.
- A global sequential LD pass to preserve cross-chromosome boundary effects.
- Parallel block jackknife for h2/rg variance estimation.

These choices are documented in `README.md` and implemented primarily in `src/ldscore.rs`,
`src/jackknife.rs`, and `src/irwls.rs`.

## Why GPU Is Deferred
- No evidence yet that users want or can use GPU acceleration.
- No local GPU available to validate correctness or performance.
- For 1000G-scale workloads, PCIe transfer and kernel launch overheads may erase gains.

## Where GPU Could Help
GPU acceleration is most promising for **large-n** datasets (biobank-scale) where:
- The ring buffer can remain resident on the GPU.
- Only the new `B` chunk is streamed each step.
- GEMM dominates runtime more than I/O or normalization.

## Candidate Approaches
### CUDA + cuBLAS (primary)
- Use cuBLAS for `AᵀB` and `BᵀB`.
- Keep `A` (ring buffer) on the GPU; stream `B` from host each chunk.
- Start by keeping r²_unbiased + L2 accumulation on CPU; move it to GPU only if
  PCIe transfer becomes a bottleneck.

### WGPU/Vulkan (not recommended)
- WGPU lacks a mature BLAS stack and reliable f64 support across devices.
- Compute shaders would require custom GEMM kernels and extensive tuning.
- Expected performance gains are unlikely to exceed OpenBLAS for typical LDSC sizes.

## Precision
- Default to `f32` for GPU path; allow optional `f64` if hardware supports it.
- Require explicit validation of numeric drift against the CPU baseline.

## Go/No-Go Criteria
- **Go** if end-to-end `ldscore` speedup is **>= 2x** on large datasets and the
  numeric deltas are within acceptable thresholds.
- **No-go** if GEMM gains are negated by data transfer or accuracy loss is too high.

## Future Plan (If Revisited)
1. Profile CPU runtime to quantify GEMM vs I/O/normalization cost at biobank scale.
2. Build a minimal CUDA POC for GEMM only (cuBLAS), leaving accumulation on CPU.
3. Validate output correctness against CPU on a fixed dataset.
4. Measure end-to-end speedup and determine whether to proceed.

## Revisit Triggers
- Evidence that users want GPU acceleration (e.g., issues, surveys, HPC usage).
- Access to a test GPU to validate correctness and performance.
