# LDSC — Rust Rewrite

[![CI](https://github.com/sharifhsn/ldsc/actions/workflows/ci.yml/badge.svg)](https://github.com/sharifhsn/ldsc/actions)
[![crates.io](https://img.shields.io/crates/v/ldsc.svg)](https://crates.io/crates/ldsc)
[![License: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](LICENSE)
[![MSRV: 1.85](https://img.shields.io/badge/rustc-1.85%2B-orange.svg)](https://blog.rust-lang.org/2025/02/20/Rust-1.85.0.html)

A compiled, statically-typed rewrite of [Bulik-Sullivan et al.'s LDSC](https://github.com/bulik/ldsc) in Rust.
Implements six subcommands — `munge-sumstats`, `l2`, `h2`, `rg`, `make-annot`, `cts-annot` — with
numerically identical output and a **~38× speedup** on LD score computation
(exact mode, 1000G N=2,490; up to **101×** with `--sketch`). Approximate modes
(`--sketch`, `--stochastic`, `--subsample`) trade per-SNP precision for additional throughput.

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

### Fast f32 mode (experimental)

Pass `--fast-f32` to `ldsc l2` to run core matmuls in f32 while accumulating
in f64. This halves the memory bandwidth for GEMM and can be significantly faster on CPUs
with 256-bit SIMD. It is **not parity-safe** compared to the default f64 path.

Observed speedup varies with panel size (hyperfine-validated, AWS EPYC 7R13):
- **N=2,490 (1000G):** ~1.3× faster vs f64 exact
- **N=50,000 (biobank):** ~1.85× faster vs f64 exact (bandwidth-limited, f32 halves per-chunk BED footprint)
- Output deltas vs f64: `max_abs_diff ≈ 0.008` (well within h2/rg regression noise)

**Note:** `--sketch` automatically enables f32 (see below), so `--fast-f32` is only needed
for exact or `--stochastic` modes.

```bash
# f32 matmul (runtime flag, no rebuild needed)
ldsc l2 --bfile … --out … --fast-f32

# default f64 (parity-safe)
ldsc l2 --bfile … --out …
```

### Stochastic trace estimation (`--stochastic`, experimental)

Pass `--stochastic T` to `ldsc l2` to use Hutchinson's stochastic trace estimator
with T random Rademacher probe vectors instead of exact GEMM. This approximates
`diag(R²)` without forming the full correlation matrix, trading precision for speed.

**This is an approximation.** Per-SNP LD scores will differ from the exact (default)
path. The expected relative error per SNP is ~sqrt(2/T). Downstream h2/rg estimates
are typically robust to this noise, but users should validate on their specific
application before relying on stochastic scores.

| T (probes) | Mean relative error | Median |rel| error | Wall time (1.66M SNPs) |
|-----------|--------------------|--------------------|----------------------|
| 50        | ~2%                | ~7%                | 36.2s (13% faster)   |
| 100       | ~1%                | ~5%                | not recommended*     |
| 1000      | ~0.4%              | ~1.5%              | slower than exact    |

\*T=100 currently exhibits a memory performance regression on some hardware. Use T=50
for the speed benefit, or omit the flag for exact computation.

```bash
# Stochastic (faster, approximate)
ldsc l2 --bfile … --out … --stochastic 50

# Exact (default, numerically identical to Python)
ldsc l2 --bfile … --out …
```

**Limitations:**
- Scalar (non-partitioned) LD scores only. If `--annot` is provided with multiple
  annotation columns, the flag is ignored and exact GEMM is used.
- Incompatible with `--gpu`.
- Results are deterministic (fixed PRNG seed) but not identical across runs with
  different T values.

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

Accuracy measured at N=2,490 (1000G). At biobank scale (N=50K), sketch speedup is bounded by
BED read bandwidth: `--sketch 50` achieves **5.3×** over exact-f64, and is the recommended
biobank mode.

```bash
# Sketch (faster, approximate — d=200 is a good default)
ldsc l2 --bfile … --out … --sketch 200

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
  the downstream d×d matmul begins to matter. d=500–1000 recommended.
- Incompatible with `--gpu` and `--stochastic`.
- Results are deterministic (fixed PRNG seed 42).
- Works with partitioned annotations (unlike `--stochastic`).

### Individual subsampling (`--subsample`, experimental)

Pass `--subsample N'` to `ldsc l2` to randomly select N' individuals from the reference
panel before LD computation. This reduces both BED I/O and GEMM cost proportionally —
runtime scales as O(N'²) instead of O(N²).

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
| **Rust stochastic** (`--stochastic 50`) | **36.2 s** | **~43×** | ~7% median per-SNP error |
| **Rust sketch** (`--sketch 200`) | **25.4 s** | **~61×** | Pearson r ≈ 0.93 vs exact |
| **Rust sketch** (`--sketch 50`) | **15.3 s** | **~101×** | Pearson r ≈ 0.81 vs exact |

`--ld-wind-kb 1000`, `--chunk-size 200`. `--sketch` automatically enables f32 (bit-identical
to f64 for sketch, ~1.3× faster). The `--fast-f32`, `--stochastic`, and `--sketch` paths
trade exact parity for speed; the default f64 path is numerically identical to Python across
all 1,664,852 SNPs.

#### Biobank scale (N = 50,000)

At biobank-scale N, GEMM dominates and larger panels expose more parallelism. Python runtime is
extrapolated from the 1000G linear O(N²) scaling; Rust runtimes are hyperfine-measured (AWS EPYC 7R13,
same hardware as above).

| Mode | Wall time | vs exact-f64 | Notes |
|------|-----------|--------------|-------|
| exact-f64 | 665.9 s | 1.0× | numerically exact |
| exact-f32 (`--fast-f32`) | 361.9 s | **1.84×** | halved BED bandwidth |
| `--sketch 50` (Rademacher) | 118.4 s | 5.6× | r ≈ 0.81 |
| `--sketch 200` (Rademacher) | 151.8 s | 4.4× | r ≈ 0.93 |
| **`--sketch 50 --sketch-method countsketch`** | **33.1 s** | **20.1×** | r ≈ 0.81 |
| **`--sketch 200 --sketch-method countsketch`** | **33.8 s** | **19.7×** | r ≈ 0.93, recommended |
| `--sketch 500 --sketch-method countsketch` | 36.2 s | 18.4× | r ≈ 0.97 |
| `--subsample 5000 --sketch 50` | 24.9 s | 26.7× | fastest (compounds two approx.) |

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

At biobank N, GEMM cost is O(N²) per window; CountSketch reduces this to O(N×c) scatter-add,
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

## Performance Deep-Dive

### Algorithmic complexity

This section documents the per-chunk and total complexity of each `l2` mode. The variables are:

| Symbol | Meaning | Typical value |
|--------|---------|---------------|
| M | Total SNPs | 1.66M (full genome) |
| N | Individuals in reference panel | 2,490 (1000G) to 500K+ (UKB) |
| c | Chunk size (`--chunk-size`) | 200 |
| w | Max window size in SNPs | ~1,000–4,000 (depends on `--ld-wind-*`) |
| d | Sketch dimension (`--sketch d`) | 50–1,000 |

The main loop iterates over M/c chunks. Each chunk reads c SNP columns, computes the within-chunk
product BB = B^T B (c × c) and the cross-window product AB = A^T B (w × c), then accumulates
unbiased r² values into the L2 array.

#### Exact mode (default)

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

#### `--fast-f32`

Same algorithm as exact, but genotype storage and GEMM use f32 instead of f64. This halves
memory bandwidth (the bottleneck at large N) for a ~1.85× speedup at N=50K. Complexity is
identical; the constant factor changes.

#### `--sketch d` with Rademacher projection (default sketch method)

Rademacher sketch pre-multiplies each decoded chunk by a fixed random matrix P (d × N) with
entries ±1/√d, reducing the individual dimension from N to d before all downstream GEMM ops.

```
Per chunk:
  BED decode + normalize .... O(N × c)        same as exact
  Projection P·B ............ O(d × N × c)    dense GEMM (the dominant cost)
  Ratio estimator renorm .... O(d × c)        rescale columns so ||col||² = N
  BB = B̃ᵀ·B̃ (within-chunk) .. O(d × c²)      c×c matmul, d inner dimension
  AB = Ãᵀ·B̃ (cross-window) .. O(d × w × c)   w×c matmul, d inner dimension
  Ring store ................. O(d × c)       copy c columns (d rows, not N)

Total: O(M × (d × N + d × (c + w)))
     = O(M × d × N)  when N >> c + w  (projection dominates)

Memory: O(d × N) for projection matrix P, O(d × (w + c)) for ring buffer + scratch,
        O(N × c) for decoded genotype buffer (b_full)
BED I/O: O(M × N / 4) bytes   same full scan
```

The projection GEMM (d × N × c per chunk) is the bottleneck. At d=50, N=50K, c=200: 500M
FLOPs/chunk. faer's SIMD matmul with Par::rayon(0) parallelizes this across all cores.

Rademacher speedup over exact: the downstream GEMM (BB, AB) shrinks from O(N) to O(d) inner
dimension, giving a theoretical N/d speedup on those phases. But the projection itself is
O(d × N × c) — same order as the exact BB — so total speedup is bounded by ~2× from GEMM
reduction alone. The real win is at large N where the O(N × c²) exact BB is the bottleneck.

**Performance scales with d**: higher d = more projection FLOPs = slower. At d=200, projection
is 4× more expensive than d=50. This is why Rademacher sketch shows increasing runtimes with d:

| d | Projection FLOPs/chunk (N=50K) | Observed biobank time |
|---|-------------------------------|----------------------|
| 50 | 500M | 129.5s |
| 100 | 1B | 138.6s |
| 200 | 2B | 159.3s |

#### `--sketch d --sketch-method countsketch` (fused CountSketch)

CountSketch replaces the dense projection GEMM with a hash-based scatter-add. Each individual
i is assigned a random bucket h(i) ∈ {0,...,d−1} and sign σ(i) ∈ {±1}. The projected
value for bucket k is the sum of σ(i) × x(i) over all individuals i with h(i) = k.

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

**Recommended d values** (accuracy vs. speed, measured on bench_200k N=2,490):

| d | Pearson r vs exact | Median relative error | Speed vs d=50 |
|---|-------------------|----------------------|---------------|
| 50 | 0.74 | 39% | 1.0× |
| 100 | 0.87 | 27% | ~1.0× |
| 200 | 0.92 | 20% | ~1.0× |
| **500** | **0.97** | **11%** | **~1.0×** |
| 1000 | 0.99 | 7% | ~1.5× slower |
| 2000 | 0.99 | 5% | ~3× slower |

**d=500 is the recommended default** — it achieves r=0.97 at no performance cost relative to
d=50 (scatter-add dominates at any realistic N). d=1000 is the "high accuracy" option for
users needing r > 0.98, with a modest ~50% slowdown from BB+AB at small N.

Accuracy depends only on d, not on N: the variance of CountSketch inner products is
Var[⟨x̃,ỹ⟩] ≈ (N² + ⟨x,y⟩²) / d for unit-normalized columns. This ratio is independent of N
after normalization, so d=500 gives the same r ≈ 0.97 whether N = 2,490 or N = 500,000.

#### `--subsample N'`

Subsampling selects N' individuals before any computation. This reduces the effective N in all
formulas above:

```
Exact + subsample:      O(M × N' × (c + w))     compute
Sketch + subsample:     O(M × (d × N' + ...))    compute
BED I/O:                O(M × N / 4) bytes        still reads full packed bytes
Decode:                 O(M × N')                 only decodes N' positions
```

BED I/O remains O(M × N / 4) because the PLINK format is SNP-major — the full byte range for
each SNP must be read even if only N' individuals are extracted. But decode and all GEMM ops
scale with N', giving near-linear speedup:

| N' | Theoretical speedup (GEMM) | Observed (biobank N=50K, AWS) |
|----|---------------------------|------------------------------|
| 50,000 (full) | 1× | 665s (exact-f64) |
| 10,000 | 25× | 90s |
| 5,000 | 100× | 51s |
| 2,000 | 625× | I/O-bound |

At N'=5K, LD score variance for common variants (MAF > 5%) increases by only ~10× relative to
N=50K, which is negligible for downstream h2/rg regression (LD scores are inputs to a second
regression that averages over ~1M SNPs).

#### `--stochastic T`

Hutchinson's trace estimator uses T random probe vectors instead of forming the full correlation
matrix. Each probe requires an O(N × M_window) matrix-vector product.

```
Per chunk:
  BED decode + normalize .... O(N × c)
  T probe MVPs .............. O(T × N × c)    T matrix-vector products
  L2 accumulation ........... O(T × c)

Total: O(M × T × N × c)

Memory: O(N × T) for probe buffer + O(N × c) for genotype chunk
```

Compared to exact (O(M × N × c × (c + w))), stochastic is faster when T < c + w. At
T=50, c=200, w~2000: T << c + w, so stochastic should be ~44× faster in theory. In practice,
the sequential MVP loop causes cache pressure at T > 50, and the method is only compatible with
scalar (non-partitioned) LD scores.

#### Summary table

| Mode | Compute per chunk | Memory | Performance scales with d? |
|------|-------------------|--------|---------------------------|
| Exact f64 | O(N × c × (c + w)) | O(N × (w + c)) | n/a |
| Exact f32 | Same, 2× bandwidth | Same, half bytes | n/a |
| Rademacher sketch | O(d × N × c) | O(d × N + d × (w + c)) | Yes — linear in d |
| **Fused CountSketch** | **O(N × c + d × c × (c + w))** | **O(N + d × (w + c))** | **No** (until d ≈ √N) |
| Subsample N' (exact) | O(N' × c × (c + w)) | O(N' × (w + c)) | n/a |
| Stochastic T | O(T × N × c) | O(N × (T + c)) | n/a (scales with T) |

### Why Python is slow

The original Python implementation is bottlenecked by three independent factors:

1. **GIL-blocked jackknife.** `jackknife.py` runs 200 leave-one-block-out IRWLS refits sequentially.
   Each refit is a `scipy.linalg.lstsq` call that releases the GIL, but Python loop overhead and
   NumPy's per-call allocation dominate at this problem size.

2. **Per-SNP NumPy allocation in the LD score loop.** `ldscore.py` calls `np.dot` in a Python-level
   loop with fresh array views on each of the ~33,000 chunks for a 1M-SNP genome. Python's boxing
   overhead and NumPy's internal allocation path are not amortised.

3. **Sequential LD computation.** The GIL prevents genuine thread-level parallelism in the
   correlation loop.

### What the Rust implementation does differently

#### 1. Ring-buffer genotype store (`l2/compute.rs`)

Python allocates a new `rfuncA` matrix every chunk. Rust pre-allocates a single F-order
`MatF` of shape `(n_indiv, ring_size)` where `ring_size = max_window + chunk_c`. SNP columns
are written into successive ring slots modulo `ring_size`; evicted slots are reused. This eliminates
~33,000 heap allocations for a 1M-SNP genome and improves cache locality because each active column
occupies a contiguous 8-byte stride in memory.

#### 2. Single matmul per chunk

For each chunk of B SNPs the computation is:

```
BB = Bᵀ · B          (chunk × chunk, unbiased r²)
AB = Aᵀ · B          (window × chunk, unbiased r²)
```

Both are single `faer` matmul calls. The window matrix `A` is assembled from ring slots into a
pre-allocated column-major `a_buf` so columns are contiguous in memory and the matmul kernel can
stride through them without gather operations.

#### 3. Threading control

Small matrix multiplications benefit from fewer threads; the `--rayon-threads` flag controls the
global Rayon thread count used by jackknife and matrix ops.

#### 4. Global sequential pass — no cross-chromosome boundary artefact

Python processes all chromosomes as a single ordered dataset. With `--ld-wind-snps`, the last 100
SNPs of chromosome k and the first 100 of chromosome k+1 are within each other's windows. The 1000G
reference panel contains five continental populations, creating population-stratification-driven
Pearson r up to 0.38 across chromosome boundaries. Earlier versions of the Rust code ran per-chromosome
in parallel, which zeroed out these cross-boundary contributions and produced L2 values 1–2 units too
low for boundary SNPs. The current implementation mirrors Python: a single global pass over all SNPs
in BIM order, with per-chromosome files written from the global L2 array after the fact.

#### 5. Parallel block jackknife (`jackknife.rs`)

The 200 leave-one-block-out IRWLS refits are independent. Rayon's `into_par_iter` distributes them
across all available cores. Each refit allocates two `faer` matrices and one SVD call; the total
wall time for h2 and rg is dominated by the file I/O and merge join, not the jackknife.

#### 6. Polars LazyFrame for munge (`munge.rs`)

`munge_sumstats.py` uses pandas, which loads the entire file into RAM before filtering. The Rust
implementation uses Polars `LazyCsvReader`, which pushes column selection, renaming, and filter
predicates into a query plan that streams the file in chunks. For large GWAS files (> 1 M SNPs) the
peak RAM is proportional to the output size, not the input size.

---

## Source Map (for LLMs)

```
src/
├── main.rs          Clap dispatch — parses CLI, calls into subcommand modules.
│
├── cli.rs           All argument structs (MungeArgs, L2Args, H2Args, RgArgs,
│                    MakeAnnotArgs, CtsAnnotArgs). No logic — pure clap derive macros.
│
├── gpu.rs           [cfg(feature = "gpu")] CUDA matmul via CubeCL/cubek-matmul.
│                    · GpuContext — holds ComputeClient + Strategy
│                    · matmul_tn — A^T × B single-shot GPU GEMM
│                    · matmul_tn_tiled — tiled variant for VRAM-limited windows
│
├── parse.rs         File I/O helpers:
│                    · scan_sumstats / scan_ldscore  → Polars LazyFrame
│                    · concat_chrs(prefix, suffix)   → concat per-chr files
│                    · read_m_total / read_m_vec      → .l2.M files
│                    · read_annot(prefix, thin)       → MatF + col names
│
├── munge.rs         munge-sumstats pipeline (Polars LazyFrame, no data loaded until
│                    collect). Internal functions:
│                    · apply_ignore, apply_col_overrides, normalize_columns
│                    · apply_info_list, apply_n_override
│                    · derive_z (BETA/SE → Z; P + sign → Z; --a1-inc)
│                    · filter_snps, apply_nstudy_filter
│                    · write_sumstats_gz (gzip TSV output)
│
├── l2/              LD score computation (split into submodules):
│   ├── mod.rs       · run — orchestrates BIM read, --extract / --annot / --keep,
│   │                  calls compute_ldscore_global, writes per-chr .l2.ldscore.gz
│   │                  and .l2.M / .l2.M_5_50 files
│   ├── compute.rs   · GemmBufs — enum holding f32 or f64 scratch buffers (--fast-f32)
│   │                · compute_ldscore_global — ring-buffer GEMM loop (sequential,
│   │                  scalar and partitioned + sketch + stochastic paths)
│   │                · r2_unbiased — r² − (1−r²)/(n−2)
│   ├── window.rs    · WindowMode — Cm / Kb / Snp enum
│   │                · get_block_lefts_f64, get_block_lefts_by_chr — window boundaries
│   ├── normalize.rs · normalize_col_f{32,64}_with_stats — impute NaN → mean,
│   │                  centre, unit-variance; AVX2+FMA sum_sumsq_f32
│   ├── snp_stats.rs · compute_snp_stats — fast BED scan for MAF + het/missing
│   └── io.rs        · parse_bim, count_fam, parse_fam — PLINK file parsers
│                    · load_individual_indices — --keep FID/IID → isize indices
│                    · write_ldscore_refs — gzip TSV output
│                    · load_snp_set — HashSet<String> from --extract / --print-snps
│
├── irwls.rs         Iteratively Re-Weighted Least Squares.
│                    · IrwlsResult — est + optional jackknife fields
│                    · irwls(x, y, weights, n_iter) — pre-alloc xw/yw, SVD solve,
│                      reweight on fitted values; zero-alloc inner loop
│
├── jackknife.rs     Block jackknife variance estimation.
│                    · jackknife(x, y, weights, n_blocks, n_iter) →
│                      full-data IRWLS + n_blocks parallel leave-one-out refits
│                      (rayon par_iter) → pseudo-values → SE + covariance matrix
│
└── regressions.rs   h2 and rg regression drivers.
                     · run_h2 — loads sumstats + LD scores, inner-joins on SNP,
                       detects K annotation columns, dispatches to scalar (K=1)
                       or partitioned (K>1) path; supports --two-step, --no-intercept,
                       --intercept-h2, --print-coefficients, liability-scale output
                     · run_h2_partitioned — K-column design matrix, per-annotation
                       enrichment, resolves M vector from per-annotation M files
                     · run_h2_scalar — shared by standalone h2 and rg univariate sub-fits
                     · run_rg — iterates trait pairs; gencov regression + univariate h2
                       per trait; --two-step, --intercept-gencov, --no-intercept,
                       liability-scale rg; prints summary table
                     · load_ld_ref / load_ld — LazyFrame readers for ref and weight LD
                     · resolve_m / resolve_m_vec — reads .l2.M_5_50 or falls back to n_obs
                     · liability_scale_h2 — observed → liability scale conversion
                     · print_jackknife_diagnostics — --print-cov / --print-delete-vals

make_annot.rs        BED → 0/1 annotation generator.
                     · annotate_from_bed — loads BED intervals per chromosome,
                       sorts and merges, binary-search annotation per SNP
                     · annotate_from_gene_set — gene symbols → coordinate lookup →
                       same interval merge/annotate pipeline
                     · write_annot_file — CHR BP SNP CM ANNOT TSV (optional .gz)
```

### Key data-flow invariants

- `l2 --annot prefix` reads `{prefix}{chr}.annot[.gz]` for every chromosome found in
  the BIM (not a single `prefix.annot.gz` file).
- `--extract` filters the BIM *before* window computation; `--print-snps` filters only the output.
- `bed_idx` (original BIM row index) differs from `pos` (index in the filtered `all_snps` slice)
  when `--extract` is active; `bed_idx_to_pos` in `run()` maps between them.
- `--keep` passes `iid_indices: Option<&[isize]>` to the internal BED reader; `n_indiv_actual` (not
  the FAM total) is used for normalization and the r²-unbiased correction.
- The ring buffer `ring_size = max_window + chunk_c` guarantees no live window slot is overwritten
  before it has been consumed in the A×B product.
- `rg --no-intercept` propagates `fixed_intercept = Some(1.0)` to both the gencov regression and
  each univariate `run_h2_scalar` call, matching Python's behaviour.

### Dependency rationale

| Crate | Version | Role |
|-------|---------|------|
| internal `bed` module | - | minimal PLINK .bed reader tailored to LDSC |
| `polars` | 0.53 | lazy CSV streaming (munge + LD score file loading) |
| `faer` | 0.24 | dense matrix algebra; SVD for IRWLS |
| `rayon` | 1 | data-parallel jackknife blocks |
| `statrs` | 0.18 | Normal CDF/quantile for P→Z conversion |
| `clap` | 4 | derive-macro CLI argument parsing |
| `anyhow` | 1 | error propagation |
| `flate2` | 1 | gzip output for .sumstats.gz and .ldscore.gz |
| `cubecl` | 0.10 | (optional, `gpu` feature) multi-backend GPU compute |
| `cubek-matmul` | 0.2 | (optional, `gpu` feature) autotuned GPU matmul |
| `fastrand` | 2 | Rademacher random generation for `--stochastic` and `--sketch` modes |
| `mimalloc` | 0.1 | (optional, `mimalloc` feature) fast allocator for musl builds |

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
