# LDSC — Rust Rewrite

[![CI](https://github.com/sharifhsn/ldsc/actions/workflows/ci.yml/badge.svg)](https://github.com/sharifhsn/ldsc/actions)
[![crates.io](https://img.shields.io/crates/v/ldsc.svg)](https://crates.io/crates/ldsc)
[![License: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](LICENSE)
[![MSRV: 1.85](https://img.shields.io/badge/rustc-1.85%2B-orange.svg)](https://blog.rust-lang.org/2025/02/20/Rust-1.85.0.html)

A compiled, statically-typed rewrite of [Bulik-Sullivan et al.'s LDSC](https://github.com/bulik/ldsc) in Rust.
Implements six subcommands — `munge-sumstats`, `l2`, `h2`, `rg`, `make-annot`, `cts-annot` — with
numerically identical output and a **~38× speedup** on LD score computation
(exact mode, 1000G N=2,490; up to **101×** with `--sketch`). Approximate modes
(`--sketch`, `--subsample`) trade per-SNP precision for additional throughput.

---

## Install

Fastest (no Rust required):

```bash
docker run --rm ghcr.io/sharifhsn/ldsc:latest --help
```

Local install options:

- Prebuilt binaries from GitHub Releases (see “Prebuilt Binaries” below).
- Cargo install (requires Rust): `cargo install ldsc`

## Quick Start

The typical LDSC workflow — preprocess summary statistics, then estimate heritability or genetic
correlation — mirrors the [upstream wiki tutorial](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation).

**Step 1: Download pre-computed European LD scores** (skip `l2` for European GWAS)

```bash
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
tar -jxvf eur_w_ld_chr.tar.bz2   # inner .l2.ldscore.gz files are already gzip-compressed
bunzip2 w_hm3.snplist.bz2
```

**Step 2: Pre-process summary statistics**

```bash
ldsc munge-sumstats \
  --sumstats my_gwas.txt \
  --N 50000 \
  --merge-alleles w_hm3.snplist \
  --out my_trait
```

**Step 3a: Estimate heritability**

```bash
ldsc h2 \
  --h2 my_trait.sumstats.gz \
  --ref-ld-chr eur_w_ld_chr/ \
  --w-ld-chr eur_w_ld_chr/ \
  --out my_trait_h2
```

**Step 3b: Estimate genetic correlation**

```bash
ldsc rg \
  --rg trait1.sumstats.gz,trait2.sumstats.gz \
  --ref-ld-chr eur_w_ld_chr/ \
  --w-ld-chr eur_w_ld_chr/ \
  --out trait1_vs_trait2
```

---

## Usage

### munge-sumstats

Pre-processes GWAS summary statistics into the `.sumstats.gz` format consumed by `h2` and `rg`.
Input summary statistics may be plain, `.gz`, or `.bz2`, and can be tab- or whitespace-delimited.

```bash
ldsc munge-sumstats \
  --sumstats my_gwas.txt.gz \
  --out output_prefix \
  [--merge-alleles w_hm3.snplist] \
  [--signed-sumstats BETA,0] \
  [--N 50000] \
  [--info-min 0.9] \
  [--maf-min 0.01]
```

Key flags: `--signed-sumstats COLNAME,null_value` tells the tool which column carries effect direction and what the
null value is (e.g. `BETA,0`, `OR,1`, `Z,0`). Without this flag the tool auto-detects from BETA/LOG_ODDS/OR/Z columns.
`--a1-inc` skips the signed column and treats all Z-scores as positive (A1 is always the risk allele).
`--merge-alleles` enforces allele concordance (mismatches are removed), matching Python behavior.
Use `--daner` or `--daner-n` for Ripke daner formats (infers N from FRQ_[A/U] headers or Nca/Nco columns).

### l2

Computes LD scores from a PLINK binary file set (`.bed/.bim/.fam`).
Annotation inputs (`.annot`) may be plain, `.gz`, or `.bz2`.

> **Tip for European GWAS:** Pre-computed 1000G phase 3 LD scores are available from the
> [Broad LDSCORE page](https://data.broadinstitute.org/alkesgroup/LDSCORE/). Download
> `eur_w_ld_chr.tar.bz2`; after `tar -jxvf`, the inner `.l2.ldscore.gz` files are already
> gzip-compressed and work directly with `ldsc`. Non-European populations require computing
> your own LD scores from an appropriate reference panel.

```bash
ldsc l2 \
  --bfile /path/to/1000G_EUR \
  --out out/eur \
  --ld-wind-cm 1.0 \
  [--annot annotations/BaselineLD.] \
  [--extract snplist.txt] \
  [--maf 0.01] \
  [--keep keep_individuals.txt] \
  [--per-allele] \
  [--pq-exp 1.0]
```

`ldsc l2` warns if the LD window spans an entire chromosome; use `--yes-really` to silence.

Window flags are mutually exclusive: `--ld-wind-cm` (genetic distance, default 1.0), `--ld-wind-kb`
(physical distance), or `--ld-wind-snps` (fixed flanking SNP count).

Partitioned LD scores with `--annot prefix`: accepts either a single `{prefix}.annot[.gz|.bz2]` file or
per-chromosome `{prefix}{chr}.annot[.gz|.bz2]` files for each chromosome present in the BIM. Outputs one L2 column
per annotation and corresponding `.l2.M` / `.l2.M_5_50` files.

`--per-allele` is equivalent to `--pq-exp 1` (weights each r² by p·(1−p)). Use `--pq-exp S` to
apply (p·(1−p))^S weighting; output columns and `.M` files receive a `_S{S}` suffix.
`--no-print-annot` suppresses the `.annot.gz` output produced by `--cts-bin`.

### h2

Estimates SNP heritability.
LD score inputs may be plain, `.gz`, or `.bz2`.

```bash
ldsc h2 \
  --h2 trait.sumstats.gz \
  --ref-ld-chr eur_w_ld_chr/ \
  --w-ld-chr eur_w_ld_chr/ \
  --out results
```

`--ref-ld-chr prefix` appends the chromosome number then `.l2.ldscore.gz`. So
`--ref-ld-chr eur_w_ld_chr/` reads `eur_w_ld_chr/1.l2.ldscore.gz` … `eur_w_ld_chr/22.l2.ldscore.gz`.
If the chromosome number falls in the middle of the filename, use `@` as a placeholder:
`--ref-ld-chr ld/chr@_scores` → `ld/chr1_scores.l2.ldscore.gz`, etc.
The same convention applies to `--w-ld-chr`.
You may pass a comma-separated list to `--ref-ld` / `--ref-ld-chr` (Python behavior);
`--w-ld` / `--w-ld-chr` must point to a single fileset.

Common options: `--no-intercept`, `--intercept-h2 VALUE`, `--two-step 30`, `--chisq-max 80`,
`--samp-prev 0.1 --pop-prev 0.01` (liability-scale conversion),
`--print-coefficients` (partitioned h2: per-annotation τ and enrichment).

**Overlapping annotations:** use `--overlap-annot` with `--frqfile-chr prefix` (or `--frqfile` for
single filesets) to match Python’s overlap-adjusted results. When enabled, LDSC writes
`<out>.results` with overlap-aware proportion/enrichment columns.

**Cell-type-specific h2:** use `--h2-cts` and `--ref-ld-chr-cts` (see the LDSC wiki for `.ldcts`
format). Output is written to `<out>.cell_type_results.txt`. Add `--print-all-cts` to report
coefficients for all CTS LD score prefixes in each line.

### rg

Estimates genetic correlations across all pairs from a list of summary statistic files.
LD score inputs may be plain, `.gz`, or `.bz2`.

```bash
ldsc rg \
  --rg trait1.sumstats.gz,trait2.sumstats.gz,trait3.sumstats.gz \
  --ref-ld-chr eur_w_ld_chr/ \
  --w-ld-chr eur_w_ld_chr/ \
  --out results
```

`--ref-ld-chr` / `--w-ld-chr` follow the same prefix convention as `h2` (see above).

Common options: `--no-intercept`, `--intercept-h2 1,1,1` (one per trait), `--intercept-gencov 0.0,0.0` (per-pair), `--two-step 30`,
`--samp-prev` / `--pop-prev` (comma-separated, one value per input file).

### make-annot

Generates 0/1 annotation files from a UCSC BED file or a gene set.

```bash
# From a BED file:
ldsc make-annot \
  --bimfile my_data.bim \
  --bed-file regions.bed \
  --annot-file output.annot.gz \
  --windowsize 100000

# From a gene set:
ldsc make-annot \
  --bimfile my_data.bim \
  --gene-set-file immune_genes.txt \
  --gene-coord-file ENSG_coord.txt \
  --annot-file output.annot.gz \
  --windowsize 100000
```

### cts-annot

Bins one or more continuous annotations into categories and writes a `.annot` file
compatible with `l2 --annot` (Python `--cts-bin` preprocessing).

```bash
ldsc cts-annot \
  --bimfile my_data.bim \
  --cts-bin DAF.txt,DIST.txt \
  --cts-breaks 0.1,0.25,0.4x10,100,1000 \
  --cts-names DAF,DIST_TO_GENE \
  --annot-file cts.annot.gz
```

---

## Installation Details

Native builds require Rust ≥ 1.85. The Rust implementation uses `faer` for dense linear algebra.

### Fast f32 mode

Pass `--fast-f32` to `ldsc l2` to run core matmuls in f32 while accumulating
in f64. This halves the memory bandwidth for GEMM and can be significantly faster on CPUs
with 256-bit SIMD.

Observed speedup varies with panel size (hyperfine-validated, AWS EPYC 7R13):
- **N=2,490 (1000G):** ~1.3× faster vs f64 exact
- **N=50,000 (biobank):** ~1.85× faster vs f64 exact (bandwidth-limited, f32 halves per-chunk BED footprint)
- Per-SNP LD score deltas vs f64: `max_abs_diff ≈ 0.008`
- **Downstream h2/rg regression: identical to f64** — h2 estimate, SE, intercept, and ratio
  all match to displayed precision. The per-SNP noise is far below the regression's sensitivity
  threshold. **This is the recommended mode for biobank-scale data when exact h2/rg is needed.**

**Note:** `--sketch` automatically enables f32 (see below), so `--fast-f32` is only needed
for exact mode.

```bash
# f32 matmul (runtime flag, no rebuild needed)
ldsc l2 --bfile … --out … --fast-f32

# default f64 (parity-safe)
ldsc l2 --bfile … --out …
```


### Random projection sketch (`--sketch`, experimental)

Pass `--sketch d` to `ldsc l2` to compress the individual dimension from N to d via
a random projection before all GEMM operations. This reduces the inner dimension of
every matrix multiply, trading precision for speed. After projection, each column is
re-normalized (ratio estimator) to reduce per-pair variance by ~2×.

**`--sketch` automatically enables f32.** Sketch projections in f64 and f32 produce
bit-identical output because the projection entries (±1/√d for Rademacher, ±1 for
CountSketch) are exactly representable in f32. There is zero accuracy cost and ~1.3×
speed gain, so f64 sketch is strictly dominated. You do not need to pass `--fast-f32`.

**This is an approximation.** All SNPs share the same random projection matrix, so
errors are correlated across SNPs. The bias is corrected in expectation via adjusted
r²_unbiased constants, but a single run can show systematic shifts. Requires d ≥ 3.

| d (dim) | Speedup (1.66M SNPs, N=2,490) | Pearson r vs exact | Median \|rel\| error | Recommended for |
|---------|-------------------------------|-------------------|---------------------|-----------------|
| 25      | 2.88× (14.3s)                 | ~0.73             | ~97%                | quick screening |
| 50      | **~101×** (15.3s)             | ~0.81             | ~52%                | rough estimates, biobank speed |
| 100     | ~82× (~19s est.)              | ~0.85             | ~13%                | moderate accuracy |
| 200     | ~61× (25.4s)                  | ~0.93             | ~6%                 | recommended default |
| 500     | ~49× (31.5s)                  | ~0.97             | ~3%                 | high accuracy   |

Accuracy measured at N=2,490 (1000G). At biobank scale (N=50K), CountSketch is recommended
over Rademacher — see [Projection methods](#projection-methods---sketch-method) below.

**Important: downstream h2/rg regression accuracy.** Sketch LD scores introduce measurement
error that causes attenuation bias in h2 regression — the h2 estimate is systematically low
and the intercept is inflated. This effect is much larger than the per-SNP LD score error
suggests. See [Downstream regression impact](#downstream-regression-impact) for detailed
benchmarks. If exact h2/rg is needed, use `--fast-f32` (exact, 1.84× faster) or CountSketch
at d ≥ 5000.

```bash
# Sketch (faster, approximate — d=200 is a good default for LD scores)
ldsc l2 --bfile … --out … --sketch 200

# High-accuracy sketch for downstream h2/rg (d=5000, ~2% h2 bias)
ldsc l2 --bfile … --out … --sketch 5000 --sketch-method countsketch

# Exact (default, numerically identical to Python)
ldsc l2 --bfile … --out …
```

#### Projection methods (`--sketch-method`)

Two projection methods are available via `--sketch-method`:

- **`rademacher`** (default) — Dense ±1/√d random matrix. Projection cost is O(d×N×c) via
  GEMM. Leverages faer's SIMD-optimized matmul for high throughput. Slightly more accurate
  than CountSketch at the same d.

- **`countsketch`** — Hash-based projection where each individual is assigned a random bucket
  b ∈ {0,...,d−1} and sign σ ∈ {±1}. Projection cost is O(N×c) — a scatter-add instead of a
  full GEMM. Faster than Rademacher at large N (fewer FLOPs), but slightly noisier because each
  individual contributes to exactly one bucket instead of all d dimensions.

| Method | d=50 median error | d=200 median error | Projection cost | Best for |
|--------|------------------|--------------------|-----------------|----------|
| Rademacher | 20.6% | 6.1% | O(d×N×c) | small-to-medium N, best accuracy |
| CountSketch | 21.8% | 7.5% | O(N×c) | large N (biobank), fastest projection |

To match Rademacher d=50 accuracy with CountSketch, use d=100–200.

```bash
# CountSketch projection (faster at large N)
ldsc l2 --bfile … --out … --sketch 100 --sketch-method countsketch
```

#### Sketch limitations

- **Rademacher:** avoid d > 500 — cache thrashing from the d×N projection matrix causes
  severe regression at d=1000.
- **CountSketch:** d up to √N is free (scatter-add is O(N×c), independent of d). Above √N,
  the downstream d×d matmul begins to matter. Even at d=10000 (20% of N=50K), CountSketch
  is still ~2× faster than exact-f32 and ~4× faster than exact-f64.
- **Downstream h2/rg bias:** sketch LD scores cause errors-in-variables attenuation in
  regression. See [Downstream regression impact](#downstream-regression-impact) for
  benchmarks. Use d ≥ 5000 for h2/rg, or `--fast-f32` for exact regression.
- Incompatible with `--gpu`.
- Results are deterministic (fixed PRNG seed 42).
- Works with partitioned annotations.

### Downstream regression impact

Sketch LD scores cause **errors-in-variables attenuation bias** in h2 regression: the h2
estimate is systematically low and the intercept is inflated. exact-f32 has zero downstream
error. For full benchmarks, see [docs/performance-deep-dive.md](docs/performance-deep-dive.md#downstream-regression-impact).

| Use case | Recommended mode | Speedup vs exact-f64 |
|----------|-----------------|---------------------|
| Exact h2/rg needed | `--fast-f32` | 1.84× |
| h2 within ~2% | `--sketch 5000 --sketch-method countsketch` | ~5× |
| h2 within ~4% | `--sketch 1000 --sketch-method countsketch` | ~9× |
| LD scores only (QC, visualization) | `--sketch 200 --sketch-method countsketch` | ~20× |

### Individual subsampling (`--subsample`, experimental)

Pass `--subsample N'` to `ldsc l2` to randomly select N' individuals from the reference
panel before LD computation. This reduces both BED I/O and GEMM cost proportionally —
runtime scales with N' instead of N.

**This is an approximation.** LD scores from a subsample are noisier than from the full
panel. For biobank-scale data (N > 10K), N' = 2,000–5,000 gives accurate LD scores for
common variants (MAF > 5%). The subsample is deterministic (seed 42) and sorted to preserve
sequential BED read order. Cannot be combined with `--keep`.

```bash
# Subsample 5000 individuals from a 50K panel
ldsc l2 --bfile biobank_50k --out out --subsample 5000

# Full panel (default)
ldsc l2 --bfile biobank_50k --out out
```

### Exact per-SNP window boundaries (`--snp-level-masking`)

By default, this tool reproduces Python LDSC's output exactly — including a known approximation
in how window boundaries are applied during LD score computation.

**The approximation:** Python LDSC processes SNPs in chunks and evicts window entries once per
chunk using the left boundary of the first SNP in each chunk (`block_left[chunk_start]`). This
means later SNPs in the same chunk can include r² contributions from SNPs that fall outside
their true window. The Python codebase comments on this explicitly, calling it a vectorization
approximation. The [Bulik-Sullivan et al. 2015 paper](https://www.nature.com/articles/ng.3211)
defines LD scores with exact per-SNP window boundaries — chunking is not mentioned.

**The fix:** `--snp-level-masking` applies the paper's definition exactly. After each GEMM,
r² entries for SNP pairs outside each SNP's true window are zeroed before accumulation into
LD scores. The implementation uses a monotonic scan (O(window + chunk)), so the overhead is
negligible — the GEMM itself is 74% of runtime and is unchanged.

**Impact (full 1000G, 1.66M SNPs, `--ld-wind-kb 1000`):**

| Metric | Value |
|--------|-------|
| SNPs with different L2 scores | 99.99% |
| Mean relative L2 difference | 11.3% (masked scores are lower) |
| h2 change — BMI | +15.7% (0.1045 → 0.1209) |
| h2 change — SCZ | +16.2% (0.3292 → 0.3825) |

The effect is larger at narrower windows:

| Window | Mean L2 rel diff | h2 change (BMI) |
|--------|-----------------|-----------------|
| 100 kb | 36.4% | +47.7% |
| 500 kb | 15.8% | +20.3% |
| 1000 kb | 11.3% | +15.7% |
| 2000 kb | 7.7% | +11.5% |

The h2 increase is driven by the regressor: lower LD scores produce a steeper regression
slope, which translates to higher heritability estimates. Whether the paper-correct LD scores
produce more accurate h2 estimates in practice is an open empirical question — the LDSC
framework is tolerant of LD score perturbations, and the Python approximation has been used
in thousands of published analyses.

```bash
# Paper-correct LD scores (exact per-SNP window boundaries)
ldsc l2 --bfile data/1000G_phase3_common_norel --out ld_scores/snp_exact \
    --ld-wind-kb 1000 --snp-level-masking

# Default (Python-identical, for reproducibility with published results)
ldsc l2 --bfile data/1000G_phase3_common_norel --out ld_scores/python_compat \
    --ld-wind-kb 1000
```

**Note:** The default (no flag) is Python-identical and should be used when comparing against
published LD score files or replicating results from the Python tool. Use `--snp-level-masking`
when computing new reference panels where theoretical correctness is preferred.

### BED prefetch for networked storage (`--prefetch-bed`)

On HPC clusters with networked filesystems (GPFS, NFS, Lustre), BED file reads travel over
the network and can block for tens of milliseconds per chunk. `--prefetch-bed` reads the next
BED chunk on a background thread while the compute thread runs GEMM, overlapping I/O latency
with computation.

```bash
# Recommended for GPFS/NFS (HPC) — first run or cold filesystem cache
ldsc l2 --bfile … --out … --prefetch-bed

# Default (no flag) — always correct, optimal for local SSD
ldsc l2 --bfile … --out …
```

**When to use it:**

| Storage | Cache state | Use `--prefetch-bed`? |
|---------|-------------|----------------------|
| Local SSD | Warm (repeated runs) | **No** — regresses ~10% |
| Local SSD | Cold (first run) | Neutral to slight benefit |
| GPFS / NFS / Lustre | Cold (typical HPC job) | **Yes** — hides network latency |
| GPFS / NFS / Lustre | Warm (same job, 2nd pass) | Probably neutral |

**Why it regresses on local SSD:** With a warm OS page cache the BED file is already in RAM,
so `read()` is just a `memcpy` from kernel pages — no blocking I/O. The reader thread still
uses CPU for 2-bit genotype decoding and competes with the rayon GEMM thread pool. On a
6-core / 12-thread machine, adding a 13th thread preempts GEMM workers and slows the whole run.

**Why it helps on GPFS:** Each `read()` call blocks waiting for data over InfiniBand. While the
reader thread waits, the CPU is free for GEMM — genuine parallelism with no resource contention.
Benefit scales with GPFS cache miss rate (highest on the first job run after data is staged).

**Note for PMACS users:** The cluster runs CentOS 7 (kernel 3.10). `io_uring` (kernel 5.1+) is
not available; `--prefetch-bed` uses standard POSIX `pread` on a background `std::thread`.
Output is always bit-identical to the default path.

### Memory-mapped BED I/O (`--mmap`)

For HPC deployments with GPFS, Lustre, or other networked filesystems, `--mmap` uses
memory-mapped I/O instead of buffered reads. This provides:

- **Zero-copy access** for the fused CountSketch path (eliminates ~20GB of memcpy at biobank scale)
- **OS-managed readahead** via `MADV_SEQUENTIAL` — GPFS can prefetch across storage nodes in parallel
- **Async prefetch** via `MADV_WILLNEED` on the next chunk — replaces `--prefetch-bed` without thread contention
- **No seek invalidation** — unlike `BufReader`, mmap'd pages stay resident once faulted

```bash
# Recommended for GPFS/Lustre HPC
ldsc l2 --bfile … --out … --mmap

# Can combine with any mode
ldsc l2 --bfile … --out … --mmap --sketch 200 --sketch-method countsketch
```

**Note:** On local SSD with warm cache, `--mmap` regresses ~15% due to page fault overhead.
Use the default (no flag) for local storage. `--mmap` is designed for networked filesystems
where it replaces both `--prefetch-bed` and buffered reads.

### Optional GPU acceleration (experimental)

Build with `--features gpu` to enable CUDA-accelerated matrix multiplication via
[CubeCL](https://github.com/tracel-ai/cubecl). Requires a CUDA toolkit at build time.

```bash
cargo build --release --features gpu
ldsc l2 --bfile … --out … --gpu              # f32 compute (default)
ldsc l2 --bfile … --out … --gpu --gpu-f64    # native f64 compute
ldsc l2 --bfile … --out … --gpu --gpu-flex32 # half-precision compute, f32 accumulation
```

Precision options:
- `--gpu`: Default f32 compute
- `--gpu-f64`: Native f64 on GPU (slower but numerically exact)
- `--gpu-flex32`: Half-precision compute with f32 accumulation (fastest, slight accuracy loss)
- `--gpu-tile-cols N`: Split large window matrices into VRAM-fitting tiles

At 1000G scale (n=2,490), GPU is transfer-bound and slower than CPU. GPU acceleration
targets biobank-scale cohorts (n >= 50k) where each chunk's GEMM is large enough for
compute to dominate PCIe transfer.

Default build (CPU only):

```bash
cargo build --release
```

### Prebuilt Binaries

Releases include Linux, macOS, and Windows archives that contain `ldsc`, `LICENSE`, and `README.md`.

```bash
# Linux (x86_64)
curl -L -o ldsc_linux-x86_64.tar.gz \
  https://github.com/sharifhsn/ldsc/releases/latest/download/ldsc_linux-x86_64.tar.gz
tar -xzf ldsc_linux-x86_64.tar.gz
./ldsc --help

# macOS (Apple Silicon)
curl -L -o ldsc_macos-aarch64.tar.gz \
  https://github.com/sharifhsn/ldsc/releases/latest/download/ldsc_macos-aarch64.tar.gz
tar -xzf ldsc_macos-aarch64.tar.gz
./ldsc --help

# Windows (x86_64)
# Download the zip from the release page and extract:
# https://github.com/sharifhsn/ldsc/releases/latest/download/ldsc_windows-x86_64.zip

```

### Docker

Images are published to the GitHub Container Registry on every push to `main` and for each version tag.

```bash
docker pull ghcr.io/sharifhsn/ldsc:latest

# Run with local data mounted
docker run --rm \
  -v /path/to/data:/data \
  ghcr.io/sharifhsn/ldsc:latest \
  h2 --h2         /data/trait.sumstats.gz \
     --ref-ld-chr /data/eur_w_ld_chr/ \
     --w-ld-chr   /data/eur_w_ld_chr/ \
     --out        /data/results
```

Version tags (`v1.2.3`) produce `:1.2.3`, `:1.2`, and `:latest`. Pushes to `main` produce a `:main`
tag and a short-SHA tag (`:sha-XXXXXXX`).

### Building from source

Requires a Rust toolchain (≥ 1.85; edition 2024 features used).

```bash
cargo build --release
# binary: target/release/ldsc
```

The release profile sets `opt-level = 3`, `lto = "thin"`, `codegen-units = 1`.

### Static binary with mimalloc (recommended for HPC)

For a fully static binary that runs on any Linux (including CentOS 7 / RHEL 7 HPC clusters),
build with musl and the `mimalloc` feature:

```bash
# Requires the musl target (one-time setup)
rustup target add x86_64-unknown-linux-musl

cargo build --release --features mimalloc --target x86_64-unknown-linux-musl
strip target/x86_64-unknown-linux-musl/release/ldsc
# Copy the single binary to your cluster — no dependencies needed
```

This produces a static-pie binary with zero shared library dependencies. The `mimalloc`
feature replaces musl's default allocator with Microsoft's [mimalloc](https://github.com/microsoft/mimalloc),
which eliminates the ~12% performance penalty musl's single-threaded allocator incurs under
rayon parallelism. On AWS EPYC benchmarks, musl + mimalloc is actually **4.5% faster** than
the default glibc build (69.9s vs 73.2s on full 1000G).

Without `--features mimalloc`, a musl build works but is ~12% slower due to allocator contention.

### Runtime tuning (optional)

The following flags are available for performance tuning:

- `--rayon-threads N`: Rayon thread count for jackknife in `h2`/`rg`.
- `--polars-threads N`: Polars thread count for CSV streaming in `munge-sumstats`.
- `--prefetch-bed`: Background BED reader thread for `l2`. Beneficial on networked filesystems (GPFS/NFS); hurts on local SSD with warm cache. See [BED prefetch](#bed-prefetch-for-networked-storage---prefetch-bed) for full guidance.
- `--mmap`: Memory-mapped BED I/O for `l2`. Recommended for GPFS/Lustre HPC; replaces `--prefetch-bed`. See [mmap](#memory-mapped-bed-io---mmap) for details.

---

## Performance

### LD score computation (`l2`)

Benchmarks on AWS c6a.4xlarge (AMD EPYC 7R13, 16 vCPU) using 1000 Genomes Phase 3
(1,664,852 SNPs, n = 2,490 individuals). Measured with `hyperfine` (1 warmup + 3 timed runs).
Static musl binary with mimalloc, AVX2+FMA target features.

#### 1000G reference panel (N = 2,490)

| Mode | Full genome wall time | vs Python | Accuracy |
|------|----------------------|-----------|----------|
| Python | 25 min 49 s | 1.0× | reference |
| **Rust f64** (default) | **41.1 s** | **~38×** | exact (`max_abs_diff = 0`) |
| **Rust f32** (`--fast-f32`) | **~32 s** | **~48×** | `max_abs_diff ≈ 0.008` |
| **Rust sketch** (`--sketch 200`) | **25.4 s** | **~61×** | Pearson r ≈ 0.93 vs exact |
| **Rust sketch** (`--sketch 50`) | **15.3 s** | **~101×** | Pearson r ≈ 0.81 vs exact |

`--ld-wind-kb 1000`, `--chunk-size 200`. `--sketch` automatically enables f32 (bit-identical
to f64 for sketch, ~1.3× faster). The `--fast-f32` and `--sketch` paths trade exact parity for speed; the default f64 path is numerically identical to Python across
all 1,664,852 SNPs.

#### Biobank scale (N = 50,000)

At biobank-scale N, GEMM dominates and larger panels expose more parallelism. Python runtime is
extrapolated assuming linear O(N) scaling; Rust runtimes are hyperfine-measured (AWS EPYC 7R13,
same hardware as above).

| Mode | Wall time | vs exact-f64 | Notes |
|------|-----------|--------------|-------|
| exact-f64 | 665.9 s | 1.0× | numerically exact |
| exact-f32 (`--fast-f32`) | 361.9 s | **1.84×** | halved BED bandwidth |
| `--sketch 50` (Rademacher) | 118.4 s | 5.6× | r ≈ 0.81 |
| `--sketch 200` (Rademacher) | 151.8 s | 4.4× | r ≈ 0.93 |
| **`--sketch 50 --sketch-method countsketch`** | **33.1 s** | **20.1×** | r ≈ 0.81 |
| **`--sketch 200 --sketch-method countsketch`** | **33.8 s** | **19.7×** | r ≈ 0.93 |
| `--sketch 500 --sketch-method countsketch` | 36.2 s | 18.4× | r ≈ 0.97 |
| `--sketch 1000 --sketch-method countsketch` | 39.7 s | 16.8× | r ≈ 0.99 |
| `--sketch 5000 --sketch-method countsketch` | ~55 s\* | ~12×\* | r ≈ 0.998, h2 ~2% low |
| `--sketch 10000 --sketch-method countsketch` | ~75 s\* | ~9×\* | r ≈ 0.999, h2 exact |
| `--subsample 5000 --sketch 50` | 24.9 s | 26.7× | fastest (compounds two approx.) |

\*d=5000 and d=10000 timings estimated from local scaling (1.7× and 2.3× vs d=200); AWS
times will differ but the relative ordering is stable. Even at d=10000 (20% of N), CountSketch
is still **~5× faster than exact-f32** and **~9× faster than exact-f64**.

Fused CountSketch reads packed BED bytes and scatter-adds directly into the d×c sketch buffer,
eliminating the N×c intermediate entirely. Cost is O(N×c) independent of d, so d=200 has the
same speed as d=50 but much better accuracy. A 28 GB container is recommended for N=50K.

Speedup varies with window size (200k-SNP extract, same machine):

| Window | Python | Rust f64 | Speedup |
|--------|--------|----------|---------|
| `--ld-wind-kb 100` | 44.2 s | 5.1 s | **8.7×** |
| `--ld-wind-kb 500` | 48.4 s | 6.2 s | **7.8×** |
| `--ld-wind-kb 1000` | 53.7 s | 8.6 s | **6.2×** |
| `--ld-wind-kb 2000` | 61.8 s | 12.9 s | **4.8×** |

#### Scaling to larger reference panels

Runtime scales with both M (SNP count) and N (individual count). The ring-buffer algorithm
keeps memory bounded by the LD window size, not the total SNP count.

**Scaling with M (fixed N=2,490, exact-f64):**

| M (SNPs) | Est. wall time (AWS EPYC 16v) | BED size | Peak memory |
|----------|-------------------------------|----------|-------------|
| 1.66M | 41 s *(measured)* | 1 GB | ~100 MB |
| 5M | ~124 s | 3 GB | ~100 MB |
| 10M | ~247 s | 6 GB | ~100 MB |
| 50M | ~1,235 s (~21 min) | 30 GB | ~100 MB |

**Scaling with N (fixed M=1.66M, full-genome):**

| N (individuals) | Mode | Wall time | vs Python |
|-----------------|------|-----------|-----------|
| 2,490 (1000G) | exact-f64 | 41.1 s *(measured)* | ~38× |
| 2,490 (1000G) | `--sketch 50` | 15.3 s *(measured)* | ~101× |
| 50,000 (biobank) | exact-f64 | 665.9 s *(measured)* | — |
| 50,000 (biobank) | countsketch-200 | 33.8 s *(measured)* | — |
| 50,000 (biobank) | countsketch-5000 | ~55 s *(estimated)* | — |
| 50,000 (biobank) | countsketch-10000 | ~75 s *(estimated)* | — |

At biobank N, GEMM cost is O(N × w × c) per chunk; CountSketch reduces this to O(N×c) scatter-add,
independent of d. Rademacher sketch reduces inner dimension from N to d but is still GEMM-bound.
BED I/O is sequential and throughput-bound; `--fast-f32` halves BED read bytes per chunk.

Additional UKBB I/O benchmarks on this machine (Apple M4, 10 CPU cores, 24 GB RAM, macOS 26.3
build 25D125). These highlight I/O-heavy workflows and the impact of the Rust pipeline’s faster
parsing and joins.

Dataset: `/Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.sample8m.tsv`
(~497 MB) for `munge-sumstats`; and `/Users/sharif/Code/ldsc/data/UKBB.ALL.ldscore/UKBB.EUR.l2.ldscore.gz`
for `h2`/`rg` (381,831 SNPs after merge).
Quick local checks for `l2` default to a 50k SNP extract via `scripts/bench_l2_py3_vs_rust.sh`.

| Workflow | Rust | Python | Speedup |
|---------|------|--------|---------|
| munge-sumstats | 3.74 s | 62.65 s | **16.75×** |
| h2 | 0.90 s | 7.81 s | **8.68×** |
| rg (two traits) | 2.93 s | 28.09 s | **9.59×** |

Reference HPC hardware (Penn PMACS) for scaling tests:
- 19 Dell C6420 quad node systems (76 compute nodes, 6,080 cores total, CentOS 7.8, 80 CPU cores per node with hyper-threading, 256 GB or 512 GB RAM per node, 56 GB/s EDR or 100 GB/s FDR InfiniBand to the filesystem, 10 Gb/s Ethernet).
- 1 Dell R940 big memory system (1.5 TB RAM, 96 CPU cores, 10 Gb/s Ethernet, 100 GB/s FDR InfiniBand).
- 2 GPU nodes (1× Nvidia Tesla P100, 512 GB RAM, 88 CPU cores, 10 Gb/s Ethernet, 100 GB/s FDR InfiniBand).
- 4.2 PB IBM Spectrum Scale (GPFS) disk storage (2 tiers, no backup).
- 1.3 PB mirrored archive tape storage.
- LSF job scheduling system.

---

## Differences from Python

### Command structure

Python LDSC consists of three separate scripts; this crate consolidates them into subcommands of a
single `ldsc` binary:

| Python | Rust |
|--------|------|
| `python munge_sumstats.py --sumstats … --out …` | `ldsc munge-sumstats --sumstats … --out …` |
| `python ldsc.py --l2 --bfile … --out …` | `ldsc l2 --bfile … --out …` |
| `python ldsc.py --h2 … --ref-ld-chr …` | `ldsc h2 --h2 … --ref-ld-chr …` |
| `python ldsc.py --rg … --ref-ld-chr …` | `ldsc rg --rg … --ref-ld-chr …` |
| `python make_annot.py --bimfile … --bed-file …` | `ldsc make-annot --bimfile … --bed-file …` |
| `python ldsc.py --cts-bin …` | `ldsc cts-annot …` |

Python's `--l2` flag (LD score estimation mode) becomes the `l2` subcommand. The `--h2` and
`--rg` flags (regression modes) become `h2` and `rg` subcommands.

### Flag compatibility

Python flag names are supported directly.

### Behavioural differences

- **`--maf` in l2**: default now matches Python (MAF prefilter before LD computation).
- **`--n-min` default**: when `--n-min` is 0, Rust now matches Python (90th percentile / 1.5).
- **`--yes-really`**: Rust warns when the LD window spans a whole chromosome and
  `--yes-really` is not set (Python errors).
- **`--chunksize`**: Python requires explicit chunking for large files; Rust uses Polars LazyFrame
  streaming and ignores chunk size for munge.
- **`--return-silly-things` / `--invert-anyway`**: accepted flags for CLI parity; Rust never clips
  results and always uses a least-squares solver (warnings emitted).
- **`--no-print-annot`**: only affects `--cts-bin` output; suppresses the `.annot.gz` file.
- **`--cts-bin` workflow**: supported directly by `ldsc l2` (also available as a separate
  preprocessor via `ldsc cts-annot`).

---

## No-op Flags (Warned)

The following Python flags are accepted for CLI parity but do not change behavior in Rust:

- `h2/rg --return-silly-things`
- `h2/rg --invert-anyway`

---

## Further Reading

- **[Performance deep-dive](docs/performance-deep-dive.md)** — algorithmic complexity for each
  mode, scaling analysis for dense SNP panels (O(M² × N) with distance-based windows),
  downstream h2/rg regression accuracy by sketch dimension, and why Python is slow.
- **[Architecture & source map](docs/architecture.md)** — module-level code map, key data-flow
  invariants, and dependency rationale.

---

## Maintainers

### Release process

Releases are cut with `cargo-release` and tagged as `vX.Y.Z`. Tag pushes trigger the release
workflow, which builds and uploads platform tarballs to GitHub Releases.

```bash
cargo release patch
cargo release patch --execute
```

### Building the image locally

Requires Docker with BuildKit (default since Docker 23):

```bash
docker build -t ldsc .
```

The multi-stage `Dockerfile` uses [cargo-chef](https://github.com/LukeMathWalker/cargo-chef) to
cache dependency compilation in a separate layer, so incremental rebuilds only recompile changed
source files.
