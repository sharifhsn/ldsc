# Performance Deep-Dive

This document covers algorithmic complexity, scaling analysis, and downstream regression
impact for each `l2` mode. For headline benchmarks, see the main [README](../README.md#performance).

## Algorithmic complexity

This section documents the per-chunk and total complexity of each `l2` mode. The variables are:

| Symbol | Meaning | Typical value |
|--------|---------|---------------|
| M | Total SNPs | 1.66M (full genome) |
| N | Individuals in reference panel | 2,490 (1000G) to 500K+ (UKB) |
| c | Chunk size (`--chunk-size`) | 200 |
| w | Max window size in SNPs | ~1,000–4,000 (depends on `--ld-wind-*`) |
| d | Sketch dimension (`--sketch d`) | 50–10,000 |

The main loop iterates over M/c chunks. Each chunk reads c SNP columns, computes the within-chunk
product BB = B^T B (c × c) and the cross-window product AB = A^T B (w × c), then accumulates
unbiased r² values into the L2 array.

**Note on w scaling:** With distance-based windows (`--ld-wind-kb`, `--ld-wind-cm`), w is
proportional to local SNP density. For a fixed genome, denser SNP panels (larger M) produce
proportionally larger w, making total compute **O(M² × N)** rather than O(M × N × w) with
fixed w. At 1.66M SNPs with `--ld-wind-kb 1000`, w ≈ 550; at 10M SNPs, w ≈ 3,300.
Only `--ld-wind-snps` gives a fixed w independent of density.

### Exact mode (default)

```
Per chunk:
  BED decode + normalize .... O(N × c)       decode 2-bit genotypes, center, scale
  BB = Bᵀ·B (within-chunk) .. O(N × c²)      c×c matmul, N inner dimension
  AB = Aᵀ·B (cross-window) .. O(N × w × c)   w×c matmul, N inner dimension
  Ring store ................. O(N × c)       copy c columns into ring buffer

Total: O(M × N × (c + w))

Memory: O(N × (w + c))   ring buffer + scratch
BED I/O: O(M × N / 4) bytes   sequential scan of packed 2-bit genotypes
```

At 1000G scale (N=2,490), the BB and AB matmuls are small and fast. At biobank scale
(N=50K), GEMM dominates: BB alone is 50K × 200² = 2B FLOPs per chunk, ×8,324 chunks.

### `--fast-f32`

Same algorithm as exact, but genotype storage and GEMM use f32 instead of f64. This halves
memory bandwidth (the bottleneck at large N) for a ~1.85× speedup at N=50K. Complexity is
identical; the constant factor changes.

### `--sketch d` (fused CountSketch)

CountSketch replaces a dense O(d × N × c) projection GEMM with a hash-based scatter-add that
is O(N × c) regardless of d. Each individual i is assigned a random bucket h(i) ∈ {0,...,d−1}
and sign σ(i) ∈ {±1}. The projected value for bucket k is the sum of σ(i) × x(i) over all
individuals i with h(i) = k. This is the only sketch method shipped (the previously available
Rademacher dense projection was removed in 5209844 — at every realistic workload it was either
unnecessary at small N or 3.6× slower than CountSketch at biobank N).

The fused kernel reads raw packed BED bytes and performs decode + normalize + scatter-add in a
single pass, eliminating the N × c intermediate buffer entirely.

```
Per chunk:
  Pass 1 — SNP statistics ... O(N × c / 4)   byte-level LUT over packed BED
  Pass 2 — fused scatter-add  O(N × c)        decode genotype, normalize, scatter to bucket
  Ratio estimator renorm .... O(d × c)        rescale columns so ||col||² = N
  BB = B̃ᵀ·B̃ (within-chunk) .. O(d × c²)      c×c matmul, d inner dimension
  AB = Ãᵀ·B̃ (cross-window) .. O(d × w × c)   w×c matmul, d inner dimension
  Ring store ................. O(d × c)

Total: O(M × (N + d × (c + w)))
     = O(M × N)  when N >> d × (c + w)  (scatter-add dominates)

Memory: O(N) for bucket/sign arrays, O(d × (w + c)) for ring buffer + scratch
        NO N×c buffer — the key memory advantage
BED I/O: O(M × N / 4) bytes   same full scan, but reads packed bytes (not decoded f32)
```

**Performance is constant in d**: the scatter-add is O(N × c) regardless of d. The output
buffer (d × c) fits in L1 cache for any d ≤ ~2000. Only the downstream BB (O(d × c²)) and AB
(O(d × w × c)) scale with d, and these are negligible until d² approaches N:

| d | Scatter-add FLOPs/chunk | BB + AB FLOPs/chunk (w=2000) | Observed local time (N=2490) |
|---|------------------------|-----------------------------|-----------------------------|
| 50 | 500K | 2M + 20M = 22M | 619ms |
| 200 | 500K | 8M + 80M = 88M | 644ms |
| 500 | 500K | 50M + 200M = 250M | 760ms |
| 1000 | 500K | 200M + 400M = 600M | 1506ms |
| 2000 | 500K | 800M + 800M = 1.6B | 2258ms |

The crossover where BB+AB cost matches scatter-add cost is approximately **d ≈ √N**:

```
BB + AB = O(d × c × (c + w))
scatter = O(N × c)
crossover: d × (c + w) ≈ N  →  d ≈ N / (c + w) ≈ N / 2200

For practical purposes (c + w ≈ 2200):
  N = 2,490   → d_crossover ≈ 1
  N = 50,000  → d_crossover ≈ 23
  N = 500,000 → d_crossover ≈ 227
```

At biobank scale (N=50K), d up to ~1000 is essentially free — the scatter-add over 50K
individuals dwarfs the d²-scaling GEMM. At UKB scale (N=500K), d up to ~5000 would be free.

**Recommended d values** (LD score accuracy, measured on bench_200k N=2,490):

| d | Pearson r vs exact | Median relative error | Speed vs d=50 |
|---|-------------------|----------------------|---------------|
| 50 | 0.74 | 39% | 1.0× |
| 100 | 0.87 | 27% | ~1.0× |
| 200 | 0.92 | 20% | ~1.0× |
| **500** | **0.97** | **11%** | **~1.0×** |
| 1000 | 0.99 | 7% | ~1.5× slower |
| 2000 | 0.99 | 5% | ~3× slower |

**The best d depends on your downstream use case** — see
[Downstream regression impact](#downstream-regression-impact) below. For LD scores alone
(visualization, QC), d=500 is sufficient. For downstream h2/rg regression, d ≥ 5000 is
recommended to avoid attenuation bias.

Accuracy depends only on d, not on N: the variance of CountSketch inner products is
Var[⟨x̃,ỹ⟩] ≈ (N² + ⟨x,y⟩²) / d for unit-normalized columns. This ratio is independent of N
after normalization, so d=500 gives the same r ≈ 0.97 whether N = 2,490 or N = 500,000.

### `--snp-level-masking`

Default LDSC (Python and `ldsc-rs` without this flag) uses chunk-level window eviction:
`block_left[chunk_start]` is applied to all c SNPs in a chunk, so SNPs late in a chunk
accumulate r² from SNPs outside their *true* per-SNP window. The LDSC paper defines LD scores
with exact per-SNP windows — chunking is an implementation detail not in the math.

`--snp-level-masking` corrects this post-GEMM at negligible cost (O(w + c) per chunk via a
monotonic cutoff scan, exploiting that `block_left` is non-decreasing across the chunk):

- B×B mask: zeroes within-chunk pairs that cross chromosome boundaries
- A×B mask: zeroes `r2u_ab[(wi, j)]` where window SNP `wi` falls outside chunk SNP `j`'s
  per-SNP `block_left`

Impact (1.66M SNPs, `--ld-wind-kb 1000`): 99.99% of SNPs get different L2 scores (avg ~11.3%
lower), and h² estimates rise ~15-16% on real GWAS (BMI: 0.1045 → 0.1209; SCZ: 0.3292 → 0.3825).
The effect grows for narrower windows: `--ld-wind-kb 100` produces +47.7% h². Default off for
exact Python parity.

### Summary table

| Mode | Compute per chunk | Memory | Performance scales with d? |
|------|-------------------|--------|---------------------------|
| Exact f64 | O(N × c × (c + w)) | O(N × (w + c)) | n/a |
| Exact f32 (`--fast-f32`) | Same, 2× bandwidth | Same, half bytes | n/a |
| **Fused CountSketch (`--sketch d`)** | **O(N × c + d × c × (c + w))** | **O(N + d × (w + c))** | **No** (until d ≈ √N) |
| `--snp-level-masking` (any mode) | + O(w + c) post-GEMM mask | unchanged | n/a |

## Downstream regression impact

Sketch LD scores introduce measurement error in the regressor (L2), which causes
**errors-in-variables attenuation bias** in h2 regression: the slope (h2) is systematically
underestimated and the intercept is inflated. This effect is much larger than the per-SNP
LD score error would suggest — a Pearson r of 0.95 on LD scores can still produce ~8%
bias in the h2 estimate.

The following table shows downstream h2 regression results using LD scores computed at
biobank scale (N=50,000, chr22, 24,624 SNPs) with simulated phenotypes (true h2=0.5).
The exact-f64 result serves as the ground truth.

| Mode | h2 estimate (SE) | h2 % error | Intercept (SE) | LD score r |
|------|-----------------|------------|----------------|------------|
| **exact-f64** | 0.4526 (0.0048) | baseline | 3.02 (0.20) | 1.0 |
| **exact-f32** | **0.4526 (0.0048)** | **0.0%** | **3.02 (0.20)** | 1.000000 |
| CS d=50 | 0.3499 (0.0099) | -22.7% | 15.66 (0.17) | 0.825 |
| CS d=100 | 0.4153 (0.0072) | -8.2% | 12.37 (0.20) | 0.911 |
| CS d=200 | 0.4142 (0.0089) | -8.5% | 8.80 (0.24) | 0.954 |
| CS d=500 | 0.4138 (0.0065) | -8.6% | 5.92 (0.23) | 0.981 |
| CS d=1000 | 0.4365 (0.0058) | -3.6% | 4.53 (0.23) | 0.991 |
| CS d=2000 | 0.4383 (0.0054) | -3.2% | 3.80 (0.22) | 0.995 |
| CS d=5000 | 0.4439 (0.0050) | -1.9% | 3.48 (0.21) | 0.998 |
| **CS d=10000** | **0.4525 (0.0049)** | **-0.02%** | **3.22 (0.21)** | 0.999 |

Key observations:

- **exact-f32 produces identical h2 to exact-f64.** Every digit matches. This is the
  recommended mode when exact downstream regression is needed — 1.84× faster at zero
  accuracy cost.

- **CountSketch attenuation persists to high d.** Even at d=500 (LD score r=0.98), h2
  is 8.6% low. The inflection point where h2 error drops below one SE is around d=1000.

- **d=5000 gives ~2% h2 bias** and is still ~3× faster than exact-f32 at biobank scale.
  **d=10000 essentially recovers the exact answer** (0.02% error) while remaining ~2× faster
  than exact-f32.

- **The intercept absorbs sketch noise.** This is expected: noisy LD scores act as
  measurement error in the regressor, causing classical attenuation bias (slope biased
  toward zero, intercept biased upward).

### Recommendation by use case

| Use case | Recommended mode | Speedup vs exact-f64 |
|----------|-----------------|---------------------|
| Exact h2/rg needed | `--fast-f32` | 1.84× |
| h2 within ~2% | `--sketch 5000` | ~5× |
| h2 within ~4% | `--sketch 1000` | ~9× |
| LD scores only (QC, visualization) | `--sketch 200` | ~20× |
| Quick screening | `--sketch 50` | ~20× |

## Scaling to denser SNP panels

With `--ld-wind-kb` (the standard), compute is O(M² × N) because w grows proportionally
to SNP density. Both Python and Rust have the same asymptotic scaling, but Rust's constant
factor is ~40-85× better due to SIMD GEMM, ring-buffer reuse, and parallel jackknife.

CountSketch breaks this scaling: its scatter-add is O(M × N) (independent of w), and
only the downstream BB+AB GEMM is O(M × d × w). At d << N, CountSketch speedup over exact
*grows* with denser SNP panels.

| Scenario | M | w (approx) | Rust exact-f32 (est.) | Rust CS-5000 (est.) | Python (est.) |
|---|---|---|---|---|---|
| 1000G SNPs | 1.66M | 550 | 6 min | ~55 s | ~8.6 hr |
| Dense array (5M) | 5M | 1,700 | ~55 min | ~5 min | ~3.2 days |
| Imputed common (10M) | 10M | 3,300 | ~3.6 hr | ~17 min | ~12.8 days |

Estimates assume linear scaling in M and M² scaling in w for exact modes; CountSketch
estimates assume O(M) scatter-add dominates. All at N=50,000.

## Why Python is slow

The original Python implementation is bottlenecked by three independent factors:

1. **GIL-blocked jackknife.** `jackknife.py` runs 200 leave-one-block-out IRWLS refits sequentially.
   Each refit is a `scipy.linalg.lstsq` call that releases the GIL, but Python loop overhead and
   NumPy's per-call allocation dominate at this problem size.

2. **Per-SNP NumPy allocation in the LD score loop.** `ldscore.py` calls `np.dot` in a Python-level
   loop with fresh array views on each of the ~33,000 chunks for a 1M-SNP genome. Python's boxing
   overhead and NumPy's internal allocation path are not amortised.

3. **Sequential LD computation.** The GIL prevents genuine thread-level parallelism in the
   correlation loop.

## What the Rust implementation does differently

### 1. Ring-buffer genotype store (`l2/compute.rs`)

Python allocates a new `rfuncA` matrix every chunk. Rust pre-allocates a single F-order
`MatF` of shape `(n_indiv, ring_size)` where `ring_size = max_window + chunk_c`. SNP columns
are written into successive ring slots modulo `ring_size`; evicted slots are reused. This eliminates
~33,000 heap allocations for a 1M-SNP genome and improves cache locality because each active column
occupies a contiguous 8-byte stride in memory.

### 2. Single matmul per chunk

For each chunk of B SNPs the computation is:

```
BB = Bᵀ · B          (chunk × chunk, unbiased r²)
AB = Aᵀ · B          (window × chunk, unbiased r²)
```

Both are single `faer` matmul calls. The window matrix `A` is assembled from ring slots into a
pre-allocated column-major `a_buf` so columns are contiguous in memory and the matmul kernel can
stride through them without gather operations.

### 3. Threading control

Small matrix multiplications benefit from fewer threads; the `--rayon-threads` flag controls the
global Rayon thread count used by jackknife and matrix ops.

### 4. Global sequential pass — no cross-chromosome boundary artefact

Python processes all chromosomes as a single ordered dataset. With `--ld-wind-snps`, the last 100
SNPs of chromosome k and the first 100 of chromosome k+1 are within each other's windows. The 1000G
reference panel contains five continental populations, creating population-stratification-driven
Pearson r up to 0.38 across chromosome boundaries. Earlier versions of the Rust code ran per-chromosome
in parallel, which zeroed out these cross-boundary contributions and produced L2 values 1–2 units too
low for boundary SNPs. The current implementation mirrors Python: a single global pass over all SNPs
in BIM order, with per-chromosome files written from the global L2 array after the fact.

### 5. Parallel block jackknife (`jackknife.rs`)

The 200 leave-one-block-out IRWLS refits are independent. Rayon's `into_par_iter` distributes them
across all available cores. Each refit allocates two `faer` matrices and one SVD call; the total
wall time for h2 and rg is dominated by the file I/O and merge join, not the jackknife.

### 6. Polars LazyFrame for munge (`munge.rs`)

`munge_sumstats.py` uses pandas, which loads the entire file into RAM before filtering. The Rust
implementation uses Polars `LazyCsvReader`, which pushes column selection, renaming, and filter
predicates into a query plan that streams the file in chunks. For large GWAS files (> 1 M SNPs) the
peak RAM is proportional to the output size, not the input size.
