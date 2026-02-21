# LDSC — Rust Rewrite

[![CI](https://github.com/sharifhsn/ldsc/actions/workflows/ci.yml/badge.svg)](https://github.com/sharifhsn/ldsc/actions)
[![crates.io](https://img.shields.io/crates/v/ldsc.svg)](https://crates.io/crates/ldsc)
[![License: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](LICENSE)
[![MSRV: 1.85](https://img.shields.io/badge/rustc-1.85%2B-orange.svg)](https://blog.rust-lang.org/2025/02/20/Rust-1.85.0.html)

A compiled, statically-typed rewrite of [Bulik-Sullivan et al.'s LDSC](https://github.com/bulik/ldsc) in Rust.
Implements the same five subcommands — `munge-sumstats`, `ldscore`, `h2`, `rg`, `make-annot` — with
identical numerical output and a 7× speedup on LD score computation.

---

## Quick Start

The typical LDSC workflow — preprocess summary statistics, then estimate heritability or genetic
correlation — mirrors the [upstream wiki tutorial](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation).

**Step 1: Download pre-computed European LD scores** (skip `ldscore` for European GWAS)

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
  --n 50000 \
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

## Installation

```bash
cargo install ldsc
```

Requires Rust ≥ 1.85 and a C toolchain with Fortran support to build the default
statically-linked OpenBLAS. If you opt into the system BLAS feature, you only
need the system OpenBLAS development package.

On Debian/Ubuntu (default static build):

```bash
sudo apt-get install cmake gfortran libgfortran-dev
```

---

## Docker

Pre-built images are published to the [GitHub Container Registry](https://ghcr.io/sharifhsn/ldsc)
on every push to `main` and for each version tag:

```bash
# Latest stable release
docker pull ghcr.io/sharifhsn/ldsc:latest

# Specific version
docker pull ghcr.io/sharifhsn/ldsc:0.1.0

# Run with local data mounted
docker run --rm \
  -v /path/to/data:/data \
  ghcr.io/sharifhsn/ldsc:latest \
  h2 --h2         /data/trait.sumstats.gz \
     --ref-ld-chr /data/eur_w_ld_chr/ \
     --w-ld-chr   /data/eur_w_ld_chr/ \
     --out        /data/results
```

Version tags (`v1.2.3`) produce `:1.2.3`, `:1.2`, and `:latest` image tags.
Pushes to `main` produce a `:main` tag and a short-SHA tag (`:sha-XXXXXXX`).

### Building the image locally

Requires Docker with BuildKit (default since Docker 23):

```bash
docker build -t ldsc .
```

The multi-stage `Dockerfile` uses [cargo-chef](https://github.com/LukeMathWalker/cargo-chef) to
cache dependency compilation in a separate layer, so incremental rebuilds only recompile changed
source files. The runtime image is `debian:bookworm-slim` plus `libgfortran5` for the
statically-linked OpenBLAS.

---

## Building from source

Requires a Rust toolchain (≥ 1.85; edition 2024 features used). OpenBLAS is linked statically — no
runtime library installation needed.

```bash
cargo build --release
# binary: target/release/ldsc
```

The release profile sets `opt-level = 3`, `lto = "thin"`, `codegen-units = 1`.

### BLAS configuration

By default this crate builds **OpenBLAS from source** and links it statically
(`blas-openblas-static`). For CI or HPC environments that prefer a system BLAS,
use the system feature instead:

```bash
# Debian/Ubuntu system OpenBLAS
sudo apt-get install libopenblas-dev pkg-config
cargo build --release --no-default-features --features blas-openblas-system
```

The system feature skips the OpenBLAS source build and links via `pkg-config`.
Keep the default static build if you want a self-contained binary.

### Runtime tuning (optional)

The following global flags are available for performance tuning but are **not
heavily battle-tested**. Use them only when needed:

- `--blas-threads N`: OpenBLAS thread count (default 4; affects all subcommands).
- `--rayon-threads N`: Rayon thread count for jackknife in `h2`/`rg`.
- `--polars-threads N`: Polars thread count for CSV streaming in `munge-sumstats`.

`ldsc --version` prints the compiled BLAS backend (e.g., `openblas-static`).

---

## Usage

### munge-sumstats

Pre-processes GWAS summary statistics into the `.sumstats.gz` format consumed by `h2` and `rg`.

```bash
ldsc munge-sumstats \
  --sumstats my_gwas.txt.gz \
  --out output_prefix \
  [--merge-alleles w_hm3.snplist] \
  [--signed-sumstats BETA,0] \
  [--n 50000] \
  [--info-min 0.9] \
  [--maf 0.01]
```

Key flags: `--signed-sumstats COLNAME,null_value` tells the tool which column carries effect direction and what the
null value is (e.g. `BETA,0`, `OR,1`, `Z,0`). Without this flag the tool auto-detects from BETA/LOG_ODDS/OR/Z columns.
`--a1-inc` skips the signed column and treats all Z-scores as positive (A1 is always the risk allele).

### ldscore

Computes LD scores from a PLINK binary file set (`.bed/.bim/.fam`).

> **Tip for European GWAS:** Pre-computed 1000G phase 3 LD scores are available from the
> [Broad LDSCORE page](https://data.broadinstitute.org/alkesgroup/LDSCORE/). Download
> `eur_w_ld_chr.tar.bz2`; after `tar -jxvf`, the inner `.l2.ldscore.gz` files are already
> gzip-compressed and work directly with `ldsc`. Non-European populations require computing
> your own LD scores from an appropriate reference panel.

```bash
ldsc ldscore \
  --bfile /path/to/1000G_EUR \
  --out out/eur \
  --ld-wind-cm 1.0 \
  [--annot annotations/BaselineLD.] \
  [--extract snplist.txt] \
  [--keep keep_individuals.txt] \
  [--per-allele] \
  [--blas-threads 4]
```

Window flags are mutually exclusive: `--ld-wind-cm` (genetic distance, default 1.0), `--ld-wind-kb`
(physical distance), or `--ld-wind-snp` (fixed flanking SNP count).

Partitioned LD scores with `--annot prefix`: expects `{prefix}{chr}.annot[.gz]` for each chromosome present in
the BIM, outputs one L2 column per annotation and corresponding `.l2.M` / `.l2.M_5_50` files.

### h2

Estimates SNP heritability.

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

Common options: `--no-intercept`, `--intercept-h2 VALUE`, `--two-step 30`, `--chisq-max 80`,
`--samp-prev 0.1 --pop-prev 0.01` (liability-scale conversion),
`--print-coefficients` (partitioned h2: per-annotation τ and enrichment).

> **Partitioned h2 note:** The standard baseline model workflow from the upstream wiki requires
> `--overlap-annot` (for overlapping annotation categories) and `--frqfile-chr` (for per-category
> M_5_50 computation). These flags are not yet implemented; see [Unimplemented Features](#unimplemented-features).

### rg

Estimates genetic correlations across all pairs from a list of summary statistic files.

```bash
ldsc rg \
  --rg trait1.sumstats.gz,trait2.sumstats.gz,trait3.sumstats.gz \
  --ref-ld-chr eur_w_ld_chr/ \
  --w-ld-chr eur_w_ld_chr/ \
  --out results
```

`--ref-ld-chr` / `--w-ld-chr` follow the same prefix convention as `h2` (see above).

Common options: `--no-intercept`, `--intercept-gencov 0.0,0.0` (per-pair), `--two-step 30`,
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

---

## Performance

Benchmarks against the original Python on a 16-core desktop (AMD Ryzen 9 5950X) using 1000 Genomes
Phase 3 (n = 2,504 individuals, `--ld-wind-snp 100`).

| Dataset | SNPs | Rust | Python | Speedup |
|---------|------|------|--------|---------|
| chr22 | 24,624 | 1 s | 7 s | **7.0×** |
| 20% genome | 333,000 | 12 s | 88 s | **7.3×** |

Correctness: all 1,664,851 SNPs in the full 1000G genome verified to match Python within 0.001
(max diff 0.000508, median 0.000250) after fixing the four algorithmic bugs described below.

---

## Differences from Python

### Command structure

Python LDSC consists of three separate scripts; this crate consolidates them into subcommands of a
single `ldsc` binary:

| Python | Rust |
|--------|------|
| `python munge_sumstats.py --sumstats … --out …` | `ldsc munge-sumstats --sumstats … --out …` |
| `python ldsc.py --l2 --bfile … --out …` | `ldsc ldscore --bfile … --out …` |
| `python ldsc.py --h2 … --ref-ld-chr …` | `ldsc h2 --h2 … --ref-ld-chr …` |
| `python ldsc.py --rg … --ref-ld-chr …` | `ldsc rg --rg … --ref-ld-chr …` |
| `python make_annot.py --bimfile … --bed-file …` | `ldsc make-annot --bimfile … --bed-file …` |

Python's `--l2` flag (LD score estimation mode) becomes the `ldscore` subcommand. The `--h2` and
`--rg` flags (regression modes) become `h2` and `rg` subcommands.

### Flag renames

| Python | Rust | Note |
|--------|------|------|
| `--ld-wind-snps` | `--ld-wind-snp` | trailing `s` removed |
| `--N` / `--N-col` | `--n` / `--n-col` | lowercase |
| `--M` | `--m-snps` | renamed |
| `--snp` / `--a1` / `--a2` / `--p` / `--frq` / `--info` | `--snp-col` / `--a1-col` / ... | `--*-col` suffix |
| `--maf-min` | `--maf` | renamed |

### Behavioural differences

- **`--maf` in ldscore**: applied as a post-filter on output (Python applies it pre-computation).
  For typical thresholds (≥ 0.01) this makes no numerical difference.
- **`--n-min` default**: Python defaults to 90th-percentile / 2; Rust defaults to 0 (disabled).
- **`--yes-really`**: Python refuses `--ld-wind-cm` spanning an entire chromosome without this flag;
  Rust has no such guard.
- **`--chunksize`**: Python requires explicit chunking for large files; Rust uses Polars LazyFrame
  streaming and ignores chunk size for munge.
- **`--return-silly-things` / `--invert-anyway`**: accepted flags for CLI parity; Rust never clips
  results and always uses a least-squares solver.

---

## Unimplemented Features

The following Python flags have no Rust equivalent and are not planned:

| Subcommand | Flag | Reason |
|------------|------|--------|
| `munge-sumstats` | `--daner` / `--daner-n` | Ripke daner format; niche use case |
| `ldscore` | `--pq-exp` | Generalised per-allele weighting; rare |
| `ldscore` | `--cts-bin` / `--cts-breaks` / `--cts-names` | Continuous variable partitioning |
| `ldscore` | `--no-print-annot` | Print-suppression flag; trivial |
| `h2` | `--h2-cts` / `--ref-ld-chr-cts` / `--print-all-cts` | Cell-type-specific h2 |
| `h2` / `rg` | `--overlap-annot` / `--frqfile-chr` | Overlapping annotation categories |
| `rg` | `--intercept-h2` | Fix per-trait h2 intercepts in rg |
| All | `.bz2` input files | Only `.gz` supported |

---

## Performance Deep-Dive

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

#### 1. Ring-buffer genotype store (`ldscore.rs`)

Python allocates a new `rfuncA` matrix every chunk. Rust pre-allocates a single F-order
`Array2<f64>` of shape `(n_indiv, ring_size)` where `ring_size = max_window + chunk_c`. SNP columns
are written into successive ring slots modulo `ring_size`; evicted slots are reused. This eliminates
~33,000 heap allocations for a 1M-SNP genome and improves cache locality because each active column
occupies a contiguous 8-byte stride in memory.

#### 2. Single DGEMM per chunk

For each chunk of B SNPs the computation is:

```
BB = Bᵀ · B          (chunk × chunk, unbiased r²)
AB = Aᵀ · B          (window × chunk, unbiased r²)
```

Both are single `ndarray::dot` calls dispatched to OpenBLAS DGEMM. The window matrix `A` is
assembled from ring slots into a pre-allocated F-order `a_buf` so columns are contiguous in memory
and the DGEMM kernel can stride through them without gather operations.

#### 3. Tuned BLAS thread count

OpenBLAS defaults to using all available cores, which creates thread-spawning overhead that
outweighs the parallelism benefit for the small matrices in 1000G-scale LD computation
(n ≈ 2,500, window ≈ 200). The Rust binary calls `openblas_set_num_threads(4)` at startup
(overridable with `--blas-threads`). This is optimal for 1000G; biobank data (n > 10,000) may
benefit from higher values.

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
across all available cores. Each refit allocates two ndarray views and one LAPACK SVD call; the total
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
├── cli.rs           All argument structs (MungeArgs, LdscoreArgs, H2Args, RgArgs,
│                    MakeAnnotArgs). No logic — pure clap derive macros.
│
├── parse.rs         File I/O helpers:
│                    · scan_sumstats / scan_ldscore  → Polars LazyFrame
│                    · concat_chrs(prefix, suffix)   → concat per-chr files
│                    · read_m_total / read_m_vec      → .l2.M files
│                    · read_annot(prefix, thin)       → Array2<f64> + col names
│
├── munge.rs         munge-sumstats pipeline (Polars LazyFrame, no data loaded until
│                    collect). Internal functions:
│                    · apply_ignore, apply_col_overrides, normalize_columns
│                    · apply_info_list, apply_n_override
│                    · derive_z (BETA/SE → Z; P + sign → Z; --a1-inc)
│                    · filter_snps, apply_nstudy_filter
│                    · write_sumstats_gz (gzip TSV output)
│
├── ldscore.rs       LD score computation. Key types and functions:
│                    · BimRecord — CHR/SNP/CM/BP/bed_idx struct
│                    · parse_bim, count_fam, parse_fam — PLINK file parsers
│                    · load_individual_indices — --keep FID/IID file → isize indices
│                    · WindowMode — Cm / Kb / Snp enum
│                    · get_block_lefts_f64, get_block_lefts_snp — window boundaries
│                    · normalize_col — impute NaN → mean, centre, unit-variance (f32)
│                    · compute_ldscore_global — ring-buffer DGEMM loop (sequential,
│                      scalar and partitioned paths share the same pre-alloc buffers)
│                    · r2_unbiased — r² − (1−r²)/(n−2)
│                    · write_ldscore_refs — gzip TSV output
│                    · load_snp_set — HashSet<String> from --extract / --print-snps
│                    · run — orchestrates BIM read, --extract / --annot / --keep,
│                      calls compute_ldscore_global, writes per-chr .l2.ldscore.gz
│                      and .l2.M / .l2.M_5_50 files
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

- `ldscore --annot prefix` reads `{prefix}{chr}.annot[.gz]` for every chromosome found in
  the BIM (not a single `prefix.annot.gz` file).
- `--extract` filters the BIM *before* window computation; `--print-snps` filters only the output.
- `bed_idx` (original BIM row index) differs from `pos` (index in the filtered `all_snps` slice)
  when `--extract` is active; `bed_idx_to_pos` in `run()` maps between them.
- `--keep` passes `iid_indices: Option<&Array1<isize>>` to `bed-reader`; `n_indiv_actual` (not
  the FAM total) is used for normalization and the r²-unbiased correction.
- The ring buffer `ring_size = max_window + chunk_c` guarantees no live window slot is overwritten
  before it has been consumed in the A×B product.
- `rg --no-intercept` propagates `fixed_intercept = Some(1.0)` to both the gencov regression and
  each univariate `run_h2_scalar` call, matching Python's behaviour.

### Dependency rationale

| Crate | Version | Role |
|-------|---------|------|
| `bed-reader` | 1 | mmap-based PLINK .bed reading; only touched pages loaded |
| `polars` | 0.53 | lazy CSV streaming (munge + LD score file loading) |
| `ndarray` + `ndarray-linalg` | 0.16 + 0.17 | dense matrix algebra; SVD for IRWLS |
| `blas-src` + OpenBLAS | 0.10 | statically-linked BLAS for DGEMM |
| `rayon` | 1 | data-parallel jackknife blocks |
| `statrs` | 0.18 | Normal CDF/quantile for P→Z conversion |
| `clap` | 4 | derive-macro CLI argument parsing |
| `anyhow` / `thiserror` | 1 / 2 | error propagation |
| `flate2` | 1 | gzip output for .sumstats.gz and .ldscore.gz |
