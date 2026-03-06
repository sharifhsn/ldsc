# GPU Feasibility Analysis for LD Score Computation

## Workload Profile

The hot inner loop of `compute_ldscore_global` is a dense GEMM:
`A^T × B` where A is `(n_indiv × window)` and B is `(n_indiv × chunk)`.

Per chunk iteration:
1. **CPU**: Read BED chunk from disk, decode genotypes
2. **CPU**: Normalize columns (impute NaN→mean, centre, scale) — branchy, NaN-checking
3. **GEMM**: `B^T × B` (within-chunk LD) and `A^T × B` (window × chunk LD)
4. **CPU**: Unbiased r² transform, annotation accumulation into l2 output

Only step 3 benefits from GPU. Steps 1-2 and 4 remain on CPU.

## Memory Working Set

Ring buffer + scratch matrices at different scales:

| Component | 1000G (n=2.5k) | Biobank (n=500k) f64 | Biobank f16 |
|-----------|----------------|----------------------|-------------|
| Ring buffer (n × ring_size) | 27 MB | 5.2 GB | **1.3 GB** |
| B chunk (n × 200) | 4 MB | 800 MB | **200 MB** |
| a_buf scratch (n × max_window) | 14 MB | 2.6 GB | **660 MB** |
| **Total** | ~50 MB | ~10 GB | **~2.5 GB** |

GPU VRAM constraints:
- P100: 16 GB — fits biobank f64 for --ld-wind-kb 1000, tight
- A100 40 GB: comfortable for f64, ample for f16
- Whole-chr no-window (n=500k, window=137k): ring = 137 GB f64 / 34 GB f16 — does not fit any single GPU

## PCIe Transfer Overhead

Each chunk iteration transfers the B matrix to GPU:

| Precision | B chunk size (n=500k, c=200) | PCIe Gen3 x16 (~12 GB/s) | PCIe Gen4 x16 (~25 GB/s) |
|-----------|------------------------------|---------------------------|---------------------------|
| f64 | 800 MB | 67 ms | 32 ms |
| f32 | 400 MB | 33 ms | 16 ms |
| f16 | 200 MB | 17 ms | 8 ms |

The ring buffer (A matrix) stays resident on GPU — only the new B chunk transfers per iteration.
Result matrix (window × chunk) is tiny (~1 MB) — transfer back is negligible.

## Compute Time Per Chunk

GEMM: (n=500k, window=700, chunk=200), FLOPs ≈ 2 × 500k × 700 × 200 = 140 GFLOP

| Hardware | f64 | f32 | f16 (tensor core) |
|----------|-----|-----|-------------------|
| 40-core Xeon (~2.5 TFLOPS f64) | 56 ms | 28 ms | N/A (Cascade Lake) |
| P100 (4.7 TFLOPS f64) | 30 ms | 15 ms | N/A (no tensor cores) |
| A100 (19.5 TFLOPS f64) | 7 ms | 2 ms | **0.4 ms** |
| H100 (67 TFLOPS f64) | 2 ms | 0.7 ms | **0.14 ms** |

## End-to-End Per-Chunk Estimate (biobank, --ld-wind-kb 1000)

| Step | CPU-only (40-core Xeon) | GPU (A100, f16 matmul) |
|------|------------------------|------------------------|
| BED read + normalize | 50-100 ms | 50-100 ms (CPU) |
| PCIe transfer B | — | 8-32 ms |
| GEMM (A^T×B + B^T×B) | 50-100 ms | **~1 ms** |
| r2u + accumulate | 5-10 ms | 5-10 ms (CPU) |
| **Total per chunk** | **~120-210 ms** | **~70-145 ms** |
| **Speedup** | baseline | **~1.5-2×** |

BED read + normalize dominates in both cases. GPU makes the matmul essentially free but doesn't help the I/O-bound steps.

## When GPU Clearly Wins

**Whole-chromosome no-window at large n_indiv**: window grows to 100k+ SNPs, matmul becomes O(n × M × chunk) per step. CPU takes minutes per chunk; GPU (with tiled out-of-core if ring exceeds VRAM) takes seconds.

**Very large n_indiv (>1M)**: matmul fraction of total time grows, PCIe fraction shrinks proportionally. The compute advantage dominates.

**Reduced precision (f16/bf16)**: 4× memory reduction enables:
- Everything fits in GPU VRAM
- PCIe transfer time halved vs f32, quartered vs f64
- Tensor core throughput 16× over f64
- Our fast-f32 results (max_abs_diff=0.001) suggest reduced precision is viable for LD scores

## When GPU Doesn't Help

**Small n_indiv (n=2,490 like 1000G)**: matmul is tiny, PCIe overhead exceeds compute time. GPU would be slower than CPU.

**--ld-wind-kb 1000 (realistic window)**: BED I/O + normalization dominate regardless. GPU gives ~1.5-2× total speedup at best — not worth the engineering complexity.

**GPFS I/O bound**: at biobank scale, reading 100+ GB BED from network filesystem may dominate all compute.

## f16 Precision Analysis

f16 has ~3.3 decimal digits of precision, range ±65504.

The viable approach is **mixed precision**:
- Store genotype matrices in f16 (saves memory + bandwidth)
- Matmul in f16 with **f32 accumulation** (exactly what tensor cores do)
- r2u computation and l2 accumulation stay in f64

The matmul computes correlations `r = (X^T Y) / n`. These get squared and bias-corrected into LD scores that only need ~3-4 significant digits for downstream h2/rg regression. Our fast-f32 results (f32 matmul + f64 accumulation) show max_abs_diff=0.001, supporting that reduced precision in the matmul is acceptable.

Going from f32→f16 matmul adds more error, but the key insight is that f16's value isn't the compute speedup — it's the **4× memory footprint reduction** that:
1. Makes biobank-scale data fit in GPU VRAM
2. Halves PCIe transfer time vs f32
3. Doubles effective CPU memory bandwidth (relevant even without GPU)

## PMACS HPC Specifics

Target: Dell C6420 nodes (80 cores, 256-512 GB RAM) + 2× Tesla P100 GPU nodes.

**P100 verdict**: marginal benefit. f64 matmul only ~2× faster than 40-core Xeon. No tensor cores for f16. PCIe Gen3 transfer overhead significant. The **22-parallel-LSF-jobs strategy** (one per chromosome on CPU nodes) gives comparable or better wall-time reduction with zero code changes.

**CPU node scaling** (from perf-log measurements on 6-core Ryzen 5600X):
- Parallel scaling saturated at ~6 threads on desktop (2.8× over single-threaded)
- 40-core Xeon with higher aggregate memory bandwidth should scale further (est. 8-15×)
- NUMA awareness matters: `numactl --localalloc` recommended for dual-socket nodes

## Practical GPU Implementation (if pursued)

1. Keep BED read + normalize on CPU
2. Ring buffer lives on GPU; stream B chunks via pinned memory + async PCIe copy
3. `A^T × B` and `B^T × B` via cuBLAS (f16 input, f32 accumulation on tensor cores)
4. Stream r2u result matrix back to CPU for annotation accumulation in f64
5. Double-buffer: while GPU computes chunk k, CPU reads + normalizes chunk k+1

Would require: cuBLAS bindings, GPU memory management, async pipeline. Significant engineering effort.

## Recommendation

**For PMACS (P100 GPUs)**: don't pursue GPU. Use parallel LSF jobs across CPU nodes.

**For future A100/H100 clusters**: GPU becomes compelling at biobank scale (n>100k) with f16 matmul + f32 accumulation, especially for whole-chromosome windows. The 4× memory reduction from f16 is more valuable than the raw compute speedup.

**Biggest practical wins available now (no GPU needed)**:
- 22 parallel LSF jobs (one per chromosome): reduces wall time from sequential ~38 min to ~5 min bottleneck
- `--rayon-threads` tuned to node core count: better intra-matmul parallelism on 40-core nodes
- fast-f32 compile flag: 1.4× speedup with 0.001 max drift (already implemented)
