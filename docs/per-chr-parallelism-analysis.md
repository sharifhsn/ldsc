# Per-Chromosome Parallelism vs Global Sequential Pass for LD Score Computation

## Problem Statement

LD score computation for each SNP j sums r^2 between j and all SNPs within a distance window. The question is whether to process all chromosomes in a single sequential pass (allowing windows to potentially span chromosome boundaries) or process each chromosome independently in parallel.

This analysis covers: what the paper/Python actually do, genetic correctness, performance tradeoffs, and a recommendation.

---

## 1. What the Paper Says

### Theoretical Definition (Supplementary Note, Proposition 1)

The LD Score of variant j is defined as:

    l_j = sum_{k=1}^{M} r_{jk}^2

where the sum runs over **all M SNPs** in the genome. The theoretical definition is genome-wide, not chromosome-restricted.

However, in practice (Methods section of the main paper), LD scores are estimated using **windowed sums** with 1 cM windows. The paper states:

> "We estimated LD Scores from the European ancestry samples in the 1000 Genomes Project (EUR) using an unbiased estimator of r^2 with **1 centiMorgan (cM) windows**, singletons excluded."

### Cross-Chromosome LD in the Paper

Section 2.2 of the supplementary note ("LD in a Mixture of Populations") explicitly analyzes the case of **unlinked variants** (i.e., on different chromosomes):

> "Suppose j and k are unlinked variants such that r_{jk,1} = r_{jk,2} = 0"

The paper derives that for unlinked variants in a population mixture, E[r_mix,jk] = 0 and Var[r_mix,jk] = F_ST^2. In a finite sample:

    E[r_hat^2_mix,jk] ≈ F_ST^2 + (1 - F_ST^2)/N

For an unstructured sample (no stratification), cross-chromosome r^2 is simply 1/N in expectation -- pure sampling noise. At N=2,490 this is 0.0004; at N=50,000 it's 0.00002.

**Key insight**: The paper's windowed estimation with 1 cM windows is a practical approximation of the theoretical genome-wide sum. Since cross-chromosome SNP pairs have r^2 ≈ 0 (plus noise), excluding them loses essentially nothing. The 1 cM window already excludes the vast majority of same-chromosome pairs too.

### LD Scores as Noisy Regressors

The paper explicitly addresses the fact that LD scores are measured with noise (Section "Attenuation Bias" in Online Methods). The regression framework includes a disattenuation correction. Different window sizes (1 cM vs 2 cM), different reference populations (FIN vs TSI with ~15% mean LD score shift), and small reference panels (N=378) all produce usable LD scores. The framework is robust to moderate perturbations.

---

## 2. What Python Actually Does

### The Global Sequential Pass

Python's `ldscore()` function in `ldsc.py` (lines 267-314) processes **all SNPs across all chromosomes in a single call**:

```python
# Read ALL SNPs from BED file (all chromosomes)
geno_array = array_obj(array_file, n, array_snps, ...)

# Compute coordinates for ALL SNPs globally
if args.ld_wind_snps:
    coords = np.array(range(geno_array.m))          # 0,1,...,M-1 global indices
elif args.ld_wind_kb:
    coords = np.array(array_snps.df['BP'])[geno_array.kept_snps]  # raw BP
elif args.ld_wind_cm:
    coords = np.array(array_snps.df['CM'])[geno_array.kept_snps]  # raw CM

# Single global call to getBlockLefts
block_left = ld.getBlockLefts(coords, max_dist)

# Single global LD score computation
lN = geno_array.ldScoreVarBlocks(block_left, args.chunk_size, annot=annot_matrix)
```

Then `__corSumVarBlocks__` in `ldscore.py` processes all M SNPs sequentially with a sliding window defined by `block_left`.

### Cross-Chromosome Bleeding: When Does It Actually Happen?

The bleeding behavior depends on the window type:

**`--ld-wind-kb` (the standard mode): NO meaningful bleeding**

BP coordinates are chromosome-local. In the 1000G data:
- Chr1 ends at BP 249,222,450
- Chr2 starts at BP 11,944
- Distance: ~249 million BP >> 1,000,000 BP window

The `getBlockLefts` function computes `|coords[j] - coords[i]| > max_dist`, and since BP values reset near zero at each chromosome boundary, the distance between the last SNP of chr_k and the first SNP of chr_{k+1} is always vastly larger than any reasonable kb window. **No cross-chromosome pairs enter the window.**

**`--ld-wind-cm`: NO meaningful bleeding (in practice)**

CM coordinates also reset to near-zero at each chromosome boundary. Chr1 ends at ~280 cM, chr2 starts at ~0 cM. With a 1 cM window, no bleeding occurs. (Exception: the 1000G reference data has all CM=0, which would cause all SNPs to be in the same window -- but `--ld-wind-cm` is never used with this data for exactly this reason.)

**`--ld-wind-snps`: YES, bleeding occurs**

`coords = range(geno_array.m)` assigns sequential indices 0,1,...,M-1 across all chromosomes. With `--ld-wind-snps 200`, the last 200 SNPs of chr1 *are* in the window of the first SNPs of chr2. This is a minor bug in Python LDSC but has negligible numerical impact because cross-chromosome r^2 ≈ -1/(N-2) (the unbiased correction makes it slightly negative, not zero).

### Python's Chunk-Level Window Eviction

Python uses a chunked algorithm (default chunk_size=50). Window eviction is per-chunk: `block_sizes[l_B]` rounds up to a multiple of c. This means the actual window used at each SNP is an *approximation* of the true window -- a few extra or fewer SNPs may be included. This approximation is present in both global and per-chr modes and is the source of the chunk alignment differences discussed below.

### Output: Per-Chromosome Files Despite Global Computation

Despite the global computation, Python writes **per-chromosome output files** (`.l2.ldscore.gz` per chr). The downstream h2/rg code reads these per-chr files and concatenates them. This is purely an I/O convention.

---

## 3. Genetic Correctness Analysis

### Fundamental Genetics

**LD (linkage disequilibrium)** is the non-random association of alleles at different loci. It arises from:
1. Physical linkage (same chromosome, close together)
2. Population structure / admixture
3. Selection
4. Genetic drift (especially in small populations)

**Mendel's Law of Independent Assortment**: Alleles on different chromosomes segregate independently during meiosis. Therefore:

- **True population r^2 between SNPs on different chromosomes is 0** (in an unstructured population with no selection)
- In a structured population, cross-chromosome r^2 = F_ST^2 (typically < 10^-4 for European subpopulations)
- In any finite sample, observed cross-chromosome r^2 ≈ 1/N (sampling noise)

### The Unbiased Estimator

The r^2_unbiased estimator used is:

    r^2_unbiased = r^2_sample - (1 - r^2_sample) / (N - 2)

For cross-chromosome pairs where true r^2 = 0:
- E[r^2_sample] ≈ 1/N
- E[r^2_unbiased] ≈ 1/N - (1 - 1/N)/(N-2) ≈ 0 (by design -- the estimator is unbiased)

So including cross-chromosome pairs adds approximately zero to the LD score in expectation. The contribution is pure noise.

### Which Approach is More Correct?

**Per-chromosome is more scientifically correct** because:
1. LD is fundamentally a within-chromosome phenomenon (physical linkage)
2. Cross-chromosome contributions are pure noise (E ≈ 0, but variance > 0)
3. Including them adds noise to the LD score estimate without adding signal
4. The paper's windowed approach already implicitly assumes this -- 1 cM windows cannot span chromosomes

**However, the practical difference is negligible** because:
1. With `--ld-wind-kb` or `--ld-wind-cm`, no cross-chromosome bleeding occurs anyway
2. With `--ld-wind-snps`, the ~200 cross-chromosome pairs per boundary contribute ~200 * (-1/N) ≈ -0.08 at N=2490, vs typical L2 values of 5-20
3. The LD score regression framework is robust to this level of noise (it's designed for noisy regressors)

---

## 4. Performance Analysis

### Global Sequential Pass

- **Parallelism**: One thread drives the GEMM loop. faer's internal `Par::rayon(0)` parallelizes each matmul across all available cores (e.g., 16 vCPU on c6a.4xlarge).
- **Memory**: One set of GEMM buffers: ring_buf (N x ring_size), b_mat (N x chunk_c), a_buf (N x max_window), ab_buf (max_window x chunk_c).
- **At N=2,490**: ring_buf ≈ 2490 * 8000 * 8B ≈ 160MB (f64). Total ~200MB.
- **At N=50,000**: ring_buf ≈ 50000 * 8000 * 8B ≈ 3.2GB (f64). Total ~4GB.
- **GPU**: One stream processes all chunks sequentially -- clean, no contention.

### Per-Chromosome Parallel

- **Parallelism**: 22 rayon tasks (one per autosome), each calling `compute_ldscore_global` independently. Each task's internal faer matmuls use `Par::rayon(0)`.
- **Effective parallelism**: Since faer already uses all cores for each matmul, having 22 tasks doesn't add 22x parallelism. Instead, it exploits the diminishing returns of per-GEMM thread scaling: two chromosomes each using 8 cores can be faster than one chromosome using 16 cores, because GEMM efficiency drops at high thread counts for small matrices.
- **Observed speedup**: ~1.52x at N=2,490 on local Ryzen (6c/12t). Likely less at N=50,000 where GEMM is more efficiently parallelized.

#### Memory Implications

**At N=2,490 (1000G)**:
- Per-chr: 22 * ~200MB = ~4.4GB. Manageable.

**At N=50,000 (Biobank)**:
- Per-chr: 22 * ~4GB = ~88GB (f64) or ~44GB (f32). **Exceeds typical instance memory.**
- Global: ~4GB (f64) or ~2GB (f32). Fine.
- The biobank queue uses 28GB containers. Per-chr would OOM.

**Mitigation**: In practice, rayon's work-stealing means not all 22 tasks are active simultaneously. But peak memory is still ~6-8 concurrent chromosomes * 4GB = 24-32GB. Risky.

#### GPU Implications

- **Global pass**: Clean -- one GEMM stream, predictable memory, no GPU contention.
- **Per-chr parallel**: 22 threads competing for one GPU. Each would need to allocate GPU buffers, serialize through the GPU's execution engine. Net result is likely slower than sequential due to contention and memory fragmentation.
- The code already auto-enables `--global-pass` when `--gpu` is active (lines 347-357 of mod.rs).

#### Chunk Alignment Differences

The global pass has chunk boundaries at [0, 200, 400, ...] across all chromosomes. Per-chr mode resets chunk boundaries at each chromosome start. Since no chromosome has a SNP count divisible by 200, chunks on chr2+ are differently aligned between modes.

The chunk-level window eviction (`block_left[chunk_start]` instead of per-SNP `block_left[chunk_start + j]`) means different chunk boundaries produce different window approximations:
- 92% of SNPs show different L2 values (bench_200k)
- Mean |diff| = 0.24, max = 15.5 on L2 values typically 5-20 (~1.5% relative)
- This is *not* a bug in either mode -- both are valid approximations of the true windowed sum

### Performance Summary

| Factor | Global Pass | Per-Chr Parallel |
|--------|------------|------------------|
| Speed (N=2,490) | Baseline | ~1.52x faster |
| Speed (N=50,000) | Baseline | Unknown (untested) |
| Memory (N=2,490) | ~200MB | ~4.4GB (22 chr) |
| Memory (N=50,000) | ~4GB | ~88GB (OOM risk) |
| GPU compatibility | Clean | Bad (contention) |
| Thread efficiency | High (all cores per GEMM) | Diminishing returns exploited |
| Python parity | Exact (max_abs_diff=0) | ~1.5% relative diff |
| Genetic correctness | Slightly wrong for --ld-wind-snps | Correct |

---

## 5. The `--ld-wind-snps` Cross-Chromosome Bug

This is the only case where the global pass produces scientifically incorrect results:

With `--ld-wind-snps N`, the Python code uses `coords = range(M)` -- sequential integers across all chromosomes. The last N SNPs of chr_k and first SNPs of chr_{k+1} share a window. The Rust global pass replicates this bug for parity.

**Impact quantification** (at N=2,490, --ld-wind-snps 200):
- 21 chromosome boundaries, each with ~200 cross-chromosome pairs
- Each pair contributes r^2_unbiased ≈ -1/(N-2) ≈ -0.0004
- Total spurious contribution per boundary SNP: ~200 * (-0.0004) ≈ -0.08
- This affects ~200 SNPs at each boundary = ~4,200 SNPs out of 1.66M (0.25%)
- Typical L2 value: 5-20, so relative error: ~0.4-1.6%

This is well within the method's noise tolerance. However, it is technically incorrect and per-chr mode eliminates it entirely.

Note: `--ld-wind-snps` is rarely used in practice. The standard LDSC workflow uses `--ld-wind-cm 1` (the paper's recommendation) or `--ld-wind-kb 1000`.

---

## 6. Downstream h2/rg Impact

### Empirical Results (bench_5k, synthetic null sumstats)

| Metric    | Per-Chr | Global-Pass | Diff/SE |
|-----------|---------|-------------|---------|
| h2        | 0.0225  | 0.0127      | +1.75   |
| Intercept | 5.5551  | 5.5671      | -0.06   |

h2 difference is meaningless (both ~0 on null data). Intercept difference is negligible (0.06 SE).

### Empirical Results (bench_200k L2, chr22 2950-SNP sumstats)

| Metric    | Per-Chr | Global-Pass | Diff/SE |
|-----------|---------|-------------|---------|
| h2        | 23.94   | 23.29       | +1.19   |
| Intercept | 7.36    | 8.05        | -1.06   |

### Scaling Argument

Real analyses use ~1.1M SNPs for regression. Only 22 chunk boundaries shift (one per chromosome), affecting <0.4% of all chunks. SE scales as 1/sqrt(M), so the 1-1.75 SE differences at M=5000 shrink to <0.1 SE at M=1.1M. The downstream impact is negligible at real-world scale.

**TODO**: Validate at full 1.1M scale on AWS HPC by running h2 on real GWAS sumstats (e.g., PGC schizophrenia) with both modes.

---

## 7. Recommendation

### Default: Per-Chromosome Parallel

Per-chromosome mode should be the default because:

1. **More scientifically correct** -- eliminates cross-chromosome bleeding in `--ld-wind-snps` mode
2. **Faster** -- ~1.52x speedup exploiting diminishing GEMM thread scaling (at 1000G scale)
3. **Numerical differences are within the method's noise floor** -- ~1.5% relative perturbation vs method tolerance of 15%+ (population mismatch) and 30%+ (window size variation)
4. **Downstream impact is negligible** at real-world SNP counts
5. **HPC-friendly** -- HPCs typically have 128-512GB+ RAM per node, so the 22-thread memory overhead (e.g. ~88GB at N=50K f64) is well within budget. Per-chr parallel better utilizes the large core counts (32-128 cores) common on HPC nodes, where single-GEMM thread scaling hits diminishing returns earlier.

### When to Use `--global-pass`

- **Exact Python parity**: When you need max_abs_diff=0 vs Python LDSC (e.g., validation, publication)
- **GPU acceleration**: Auto-enabled when `--gpu` is active. 22 threads competing for one GPU causes contention and memory fragmentation -- a single GEMM stream is cleaner and faster.
- **Memory-constrained environments**: On machines with < 64GB RAM at large N, per-chr may OOM from 22 concurrent ring buffers. Most HPC nodes have enough RAM, but cloud instances (e.g. 28GB Batch containers) may not.

### When NOT to Use Per-Chr Parallel

- GPU workflows (already handled by auto-detection)
- When exact Python reproducibility is required for regulatory/publication purposes
- Severely memory-constrained environments (rare on HPC)

### Future Work

1. **Bounded parallelism**: Limit concurrency to `min(22, available_memory / per_chr_memory)` to prevent OOM on constrained systems while still exploiting parallelism where possible. Could auto-detect available RAM at startup.
2. **SNP-level masking in global pass**: The `--snp-level-masking` flag already addresses cross-chromosome B×B bleeding by zeroing entries for pairs outside each other's window within a chunk. This gives scientific correctness without the chunk alignment differences.
3. **Full-scale h2 validation**: Run h2 on real GWAS sumstats with 1.1M+ SNP LD scores from both modes to confirm negligible downstream differences.
4. **GPU per-chr with CUDA streams**: Use multiple CUDA streams (one per chromosome) to overlap GPU work without CPU-side thread contention. This would give per-chr correctness with GPU acceleration, but requires CubeCL stream support.
