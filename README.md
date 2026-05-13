# LDSC — Rust Rewrite

[![CI](https://github.com/sharifhsn/ldsc/actions/workflows/ci.yml/badge.svg)](https://github.com/sharifhsn/ldsc/actions)
[![crates.io](https://img.shields.io/crates/v/ldsc.svg)](https://crates.io/crates/ldsc)
[![License: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](LICENSE)
[![MSRV: 1.85](https://img.shields.io/badge/rustc-1.85%2B-orange.svg)](https://blog.rust-lang.org/2025/02/20/Rust-1.85.0.html)

**ldsc-rs** is a from-scratch Rust reimplementation of [Bulik-Sullivan et al.'s
LDSC](https://github.com/bulik/ldsc) — the standard tool for SNP-heritability
(h²) and genetic-correlation (rg) estimation from GWAS summary statistics. All
six original subcommands (`munge-sumstats`, `l2`, `h2`, `rg`, `make-annot`,
`cts-annot`) with numerically identical default output and a **~38× speedup**
on LD-score computation (1000 Genomes reference, N=2,490). CLI is compatible
with Python LDSC's `--l2`/`--h2`/`--rg` syntax for drop-in replacement in
pipelines like [NCI LDlink](https://ldlink.nih.gov/).

LDSC computes per-SNP **LD scores** ℓ\_j = Σ\_k r²\_{jk} (the sum of squared
correlations between SNP j and other SNPs in a window), then regresses GWAS
chi-squares on those LD scores. The slope is proportional to SNP heritability;
two GWAS regressed jointly recover their genetic correlation. The compute
bottleneck is the LD-score step at biobank scale: a custom 50K-individual
reference panel takes 109 minutes in Python LDSC (measured on r6a.8xlarge,
32 vCPU, 240 GiB) and 11 minutes in Rust `--fast-f32` exact mode; the standard
CountSketch shortcut cuts the time further but introduces systematic bias into
downstream h².

## What's different from Python LDSC

1. **Drop-in replacement at the default.** Output is numerically identical to
   Python LDSC when invoked with `--python-compat` (`max_abs_diff = 0` across
   18,627 chr22 1000G EUR SNPs verified bit-identical). The h² regression is
   bit-identical to Python on identical LD-score inputs. ~38× faster than
   Python on LD-score computation at 1000G; ~18× faster at biobank (Rust
   `--snp-level-masking --fast-f32` vs measured Python); ~9× on h²/rg.

2. **Principled CountSketch (`--sketch d`).** The fused
   BED-decode → normalize → scatter-add kernel eliminates the N×c intermediate
   matrix entirely. Bias is corrected by inverting the renormalized-cosine
   bias `(1−r²)(1−2r²)/d` exactly (residual O(1/d²)) — derivation and
   Monte-Carlo validation in
   [`docs/countsketch-math-analysis.md`](docs/countsketch-math-analysis.md).
   At biobank N=50K, sketch is 17–24× faster than `--fast-f32` exact.

3. **`--sketch 1600` reproduces both truth clusters at biobank speed.** At
   d=1600 the masking flag becomes a switch between the two h² definitions
   used in the field, both within 0.001 h² on three real GWAS at biobank:
   - `--sketch 1600` (plain) → matches **Python LDSC** h² (chunked-window
     algorithm, the canonical implementation everyone has used since 2014).
   - `--sketch 1600 --snp-level-masking` → matches **per-SNP exact** h² (the
     LDSC paper's mathematical definition of ℓ\_j = Σ\_k r²\_jk).

   Both run in 21 s — that's **311× faster than Python LDSC at biobank** (21 s
   vs measured 6,541 s on r6a.8xlarge) and matches whichever truth you need.

   And with `--snp-level-masking` on, sketch+mask is **~12% more accurate
   than Python LDSC** at biobank: Python's chunked-window approximation has
   systematic 10–14% h² bias vs the paper's per-SNP exact definition (verified:
   Python h² = Rust chunk-exact h² to 4 decimals, both with Δ −0.015 / −0.020
   / −0.047 vs per-SNP exact). Sketch+mask gets to the paper definition that
   Python's algorithm structurally can't.

   Cross-validated against GCTA, Python LDSC, chunk-exact, and per-SNP exact
   — full validation in [`docs/perf-log.md`](docs/perf-log.md) 2026-05-12 /
   2026-05-13 entries.

## Where this sits in the LD-score-regression tool landscape

Three implementations exist for computing LD scores from a reference panel —
producing slightly different LD-score values that all agree on h² to ~1%
on standard EUR-on-EUR setups:

| Tool | Windowing | Bias correction | Speed (1000G N=2,490) |
|---|---|---|---|
| **GCTA** `--ld-meanrsq` | Block-with-overlap (two-pass averaging) | `--ld-score-adj` | C++ baseline |
| **Python LDSC** | 50-SNP chunked (undocumented approximation) | unbiased r² | 25 min 49 s |
| **ldsc-rs** | 200-SNP chunked (default) → per-SNP exact (`--snp-level-masking`) | unbiased r² + quadratic CountSketch correction | 41 s exact, ~25 s `--sketch 200` |

GCTA was the tool the original LDSC paper used to compute the published 1000G
reference scores; Python LDSC is the canonical implementation everyone has
replicated since 2014; ldsc-rs adds biobank-scale CountSketch + the
`--sketch d --snp-level-masking` combo that reaches the LDSC paper's exact
mathematical definition at sketch speed. For h² and rg regression itself,
ldsc-rs is bit-identical to Python LDSC given identical LD-score inputs.

For h² estimation specifically (rather than LD-score generation), other
tools occupy adjacent niches: BOLT-LMM and REGENIE estimate h² directly from
individual-level genotypes via mixed models (more powerful, more compute);
ldsc-rs and Python LDSC operate on summary statistics only (no genotypes
needed at h² time).

---

## Install

Fastest (no Rust required):

```bash
docker run --rm ghcr.io/sharifhsn/ldsc:latest --help
```

Local install options:

- Prebuilt binaries from GitHub Releases (see [Prebuilt Binaries](#prebuilt-binaries) below).
- Cargo install (requires Rust ≥ 1.85): `cargo install ldsc`

## Quick Start

Two common workflows. Both produce numerically identical h²/rg output to
Python LDSC.

### A) Replicating a published EUR analysis (uses pre-computed LD scores)

Skip LD-score computation by downloading the standard 1000G EUR reference:

```bash
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
tar -jxvf eur_w_ld_chr.tar.bz2
bunzip2 w_hm3.snplist.bz2

# Munge sumstats → standardized format
ldsc munge-sumstats --sumstats my_gwas.txt --N 50000 \
                    --merge-alleles w_hm3.snplist --out my_trait

# h² (single trait) or rg (pair)
ldsc h2 --h2 my_trait.sumstats.gz --ref-ld-chr eur_w_ld_chr/ \
        --w-ld-chr eur_w_ld_chr/ --out my_trait_h2

ldsc rg --rg trait1.sumstats.gz,trait2.sumstats.gz \
        --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ \
        --out trait1_vs_trait2
```

### B) Computing your own LD scores at biobank scale

For a biobank-scale reference panel (N up to 50K+), the recommended high-
accuracy fast mode is `--sketch d --snp-level-masking`. At d=1600 this
matches per-SNP exact h² within 0.001 at ~17× the speed of exact:

```bash
# High-accuracy biobank LD scores: ~21s on N=50K, full 1.66M SNPs
ldsc l2 --bfile my_biobank_panel \
        --ld-wind-kb 1000 --sketch 1600 --snp-level-masking --mmap \
        --out my_ld_scores

# Then h² with your custom LD scores
ldsc h2 --h2 my_trait.sumstats.gz \
        --ref-ld-chr my_ld_scores --w-ld-chr my_ld_scores \
        --out my_trait_h2
```

For 1000G-scale references (N≈500-2,500), the speed difference between sketch
and exact is small enough that exact is preferable: `ldsc l2 --bfile … --out …
--ld-wind-kb 1000 --snp-level-masking --fast-f32`.

## Choosing your `l2` mode

Speedups below are biobank wall-clock (N=50,000, 1.66M SNPs, AWS c6a.4xlarge
unless noted). The Python LDSC baseline is **6,541 s (109 min)** measured on
r6a.8xlarge (32 vCPU, 240 GiB) at `--chunk-size 200 --ld-wind-kb 1000`.

**At d=1600, one flag toggles which "truth" the sketch reproduces** (verified
on 3 real GWAS at biobank, all matches within 0.001 h²):
- `--sketch 1600` (plain) → matches Python LDSC h² (chunked-window algorithm)
- `--sketch 1600 --snp-level-masking` → matches per-SNP exact h² (LDSC paper definition)

Both run in ~21 s — masking is a post-GEMM cutoff scan with negligible cost.

| Goal | Use | Time | vs Python | h² match (3 real GWAS) |
|---|---|---:|---:|---|
| Match Python LDSC bit-for-bit (1000G replication) | `--python-compat` | — | 1.0× (1000G) | bit-identical |
| Exact LD scores, paper-canonical math | `--snp-level-masking --fast-f32` | 361 s | **18×** | per-SNP exact (truth) |
| **Reproduce Python LDSC's h² at biobank** (within 0.001) | **`--sketch 1600 --mmap`** | **21 s** | **311×** | Python within 0.001 |
| **Reach paper-canonical h² at biobank** (within 0.001) | **`--sketch 1600 --snp-level-masking --mmap`** | **21 s** | **311×** | per-SNP exact within 0.001 |
| Faster, ~0.003 h² shift on Height | `--sketch 1000 --snp-level-masking --mmap` | 19 s | 344× | per-SNP exact within 0.003 |
| Fastest LD scores (QC, visualization, screening) | `--sketch 200 --mmap` | 15 s | 436× | per-SNP exact within 0.013 |

See [Performance](#performance) for measured timings; full cross-method h²
validation tables are in `docs/perf-log.md` 2026-05-12 / 2026-05-13 entries.

---

## Differences from Python

For evaluation: what changes if you switch from Python LDSC to ldsc-rs?

### Output

- **At the default**: identical to Python LDSC for h²/rg with bit-identical
  LD scores when you pass `--python-compat` to `l2`.
- **Without `--python-compat`**: ldsc-rs `l2` defaults to `--chunk-size 200`
  (vs Python's 50) for ~4× faster GEMM. Chunked windowing changes per-SNP LD
  scores by ~0.2 mean / ~6 max on chr22 EUR (Pearson r=0.9997). Downstream
  h² impact is <1% on EUR-on-EUR setups, but grows to 10–14% on biobank-scale
  multi-population references — see the [`--snp-level-masking`
  section](#--snp-level-masking--paper-canonical-per-snp-windows).
- **`--sketch d`** adds the CountSketch approximation; combinable with
  `--snp-level-masking` to recover paper-exact h². See
  [Choosing your `l2` mode](#choosing-your-l2-mode).

### Command structure

Python LDSC consists of three separate scripts; ldsc-rs consolidates them
into subcommands of a single `ldsc` binary, while still accepting the
original Python flag syntax:

| Python | Rust |
|--------|------|
| `python munge_sumstats.py --sumstats … --out …` | `ldsc munge-sumstats --sumstats … --out …` |
| `python ldsc.py --l2 --bfile … --out …` | `ldsc l2 --bfile … --out …` (or `ldsc --l2 …`) |
| `python ldsc.py --h2 … --ref-ld-chr …` | `ldsc h2 --h2 … --ref-ld-chr …` (or `ldsc --h2 …`) |
| `python ldsc.py --rg … --ref-ld-chr …` | `ldsc rg --rg … --ref-ld-chr …` (or `ldsc --rg …`) |
| `python make_annot.py --bimfile … --bed-file …` | `ldsc make-annot --bimfile … --bed-file …` |
| `python ldsc.py --cts-bin …` | `ldsc cts-annot …` |

### Behavioural differences

- **`--maf` in l2**: default matches Python (MAF prefilter before LD computation).
- **`--n-min` default**: when `--n-min` is 0, matches Python (90th percentile / 1.5).
- **`--yes-really`**: Rust warns when the LD window spans a whole chromosome and
  `--yes-really` is not set (Python errors).
- **`--chunksize` in munge**: Polars LazyFrame streaming ignores chunk size for munge.
- **`--cts-bin` workflow**: supported directly by `ldsc l2` and via the separate
  `ldsc cts-annot` preprocessor.
- **`--return-silly-things` / `--invert-anyway`**: accepted for CLI parity but no-ops;
  ldsc-rs never clips results and always uses a least-squares solver (warning emitted).
- **`--no-print-annot`**: only affects `--cts-bin` output (suppresses `.annot.gz`).

### Drop-in via argv[0] symlinks

If the binary is invoked via a symlink whose filename contains "munge"
(e.g. `munge_sumstats.py`), it auto-routes to the `munge-sumstats`
subcommand. Same for `--l2` / `--h2` / `--rg` flag-style invocation:

```bash
ln -s /usr/local/bin/ldsc /usr/local/bin/munge_sumstats.py
munge_sumstats.py --sumstats file.txt --out munged  # works
```

stdout output matches the Python format for downstream parsing (the
"Total Observed scale h2: ...", Lambda GC, Mean Chi^2, Intercept, Ratio
lines that pipelines like [NCI LDlink](https://ldlink.nih.gov/) parse).

---

## Subcommand reference

### munge-sumstats

Pre-processes GWAS summary statistics into the `.sumstats.gz` format consumed by `h2` and `rg`.
Input may be plain, `.gz`, or `.bz2`, and tab- or whitespace-delimited.

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

See [Speed/accuracy modes](#speedaccuracy-modes-for-l2) below for `--sketch`, `--snp-level-masking`,
`--fast-f32`, `--mmap`, and `--gpu`.

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

## Speed/accuracy modes for `l2`

Five flags affect the speed/accuracy tradeoff in `l2`. They compose freely.

### `--fast-f32` — exact, ~1.85× faster

Run core matmuls in f32 with f64 accumulation. Halves GEMM bandwidth.

- N=2,490 (1000G): ~1.3× faster vs f64 exact
- N=50,000 (biobank): ~1.85× faster vs f64 exact
- Per-SNP LD score deltas vs f64: `max_abs_diff ≈ 0.008`
- **Downstream h2/rg: identical to f64** to displayed precision

`--sketch` automatically enables f32; `--fast-f32` is only needed for exact mode.

### `--sketch d` — CountSketch random projection

Compress the individual dimension from N to d via a CountSketch
projection (a hash-based ±1 random matrix) before all GEMM operations.
Fused decode-normalize-scatter-add kernel; cost is **O(N×c) regardless
of d** until d ≈ √N. After projection, columns are renormalized (ratio
estimator) so ‖x̃‖² = N exactly. Bias `(1−r²)(1−2r²)/d` is corrected by
quadratic inversion (residual O(1/d²)). See
[`docs/countsketch-math-analysis.md`](docs/countsketch-math-analysis.md)
for derivation and Monte-Carlo validation.

LD-score level accuracy as a function of d (measured on a 200k-SNP extract,
N=2,490):

| d | LD-score Pearson r vs chunk-exact |
|---:|:---:|
| 100 | 0.87 |
| **200** (default) | **0.92** |
| 500 | 0.97 |
| 1000 | 0.99 |

`--sketch` automatically enables f32 (CountSketch ±1 entries are exactly
representable in f32). Deterministic (seed=42). Requires `3 ≤ d < N`
(ldsc-rs errors if out of bounds); warns if `d ≤ 50` (Taylor truncation
+ sqrt amplification breaks the bias correction). Incompatible with
`--gpu`. Works with partitioned annotations.

**Plain `--sketch d` carries downstream h² attenuation bias** even with the
quadratic correction — the residual error comes not from the sketch itself
but from the chunked-window approximation that the sketch path inherits.
Combine with `--snp-level-masking` to fix this.

#### Choosing d for your reference panel

The optimal d depends on N. Two facts to understand:

1. **Cost scales with d slowly until the GEMM crossover** `d* ≈ √N`.
   `T(d) ≈ T_scatter + T_GEMM × d`. Below `d*` doubling d costs almost
   nothing; above `d*` cost grows linearly. Empirically: `d* ≈ 260` at
   1000G N=2,490; `d* ≈ 2,545` at biobank N=50K. So you have a wide
   "free range" of d at biobank scale and a narrow one at 1000G scale.

2. **Sketch attenuation bias on h² grows with N at fixed d.** Sketch
   per-pair noise variance is ~2/d regardless of N, but at larger N
   the true LD-score signal is tighter, so the relative noise impact
   on the regression is larger. Larger N → need larger d for the same
   h² accuracy.

Recommendations (combine with `--snp-level-masking` for the h² rows):

| N (reference panel) | Use | Speedup vs exact-f32 |
|---:|---|---:|
| ~500 – 2,500 (1000G-scale) | `--snp-level-masking --fast-f32` (no sketch) | 1.0× (exact is already fast at this N) |
| ~10,000 | `--sketch 1000 --snp-level-masking` | ~5× (extrapolated) |
| **~50,000 (biobank)** | **`--sketch 1600 --snp-level-masking`** | **~17× (measured)** |
| ~500,000 (UKBB-scale) | `--sketch 5000 --snp-level-masking` (untested) | ~30-50× (extrapolated) |

The biobank row is the well-tested one (full validation in
[`docs/perf-log.md`](docs/perf-log.md) 2026-05-12 entry). UKBB-scale d
follows roughly `d ≈ √N`. For LD-score-only workflows (visualization,
QC — no h² downstream), `--sketch 200` is fine at any N: LD-score
Pearson r ≥ 0.92 regardless of N (per the small-N table above; the per-
pair variance bound `2/d` is N-independent for normalized columns).

### `--snp-level-masking` — paper-canonical per-SNP windows

The LDSC paper defines ℓ\_j = Σ\_k r²\_{jk} over an exact per-SNP window.
Python LDSC and ldsc-rs both default to a chunked approximation that uses
one common window for all SNPs in a chunk (faster GEMM, slightly inflated
LD scores at chunk edges).

`--snp-level-masking` zeroes r² entries outside each SNP's true per-SNP
window after the GEMM. Negligible cost (O(window + chunk) monotonic scan).

The h² impact depends strongly on (LD-reference N, population match):

| Setup | Windowing gap (chunked vs per-SNP exact) |
|-------|-----------------------------------------:|
| 1000G EUR, N=503 (typical h² use case) | 1–2% h² |
| Biobank synthetic, multi-pop, N=50,000 | 10–14% h² |

The gap grows with N because higher-N r² estimates are tighter, so the
"spurious" pairs at chunk boundaries carry real LD signal not noise. See
[`docs/perf-log.md`](docs/perf-log.md) 2026-05-12 entry for the full
cross-method h² breakdown.

### `--sketch d --snp-level-masking` — the recommended biobank mode

At d=1600, the masking flag toggles which truth cluster the sketch
reproduces — both within 0.001 h² across BMI 2010 / BMI 2018 / Height
2018 at biobank scale:

```bash
# Match Python LDSC chunked h² (the canonical implementation), 311× faster
ldsc l2 --bfile … --out … --ld-wind-kb 1000 --sketch 1600 --mmap            # ~21s

# Match per-SNP exact h² (the LDSC paper's mathematical definition)
ldsc l2 --bfile … --out … --ld-wind-kb 1000 --sketch 1600 --snp-level-masking --mmap  # ~21s

# 1000G: smaller speed advantage; exact (~4× slower) is usually preferred
ldsc l2 --bfile … --out … --ld-wind-kb 1000 --snp-level-masking --fast-f32
```

The two flag combinations run at the same speed (masking is a post-GEMM
cutoff scan with negligible cost). Pick which truth you want to reproduce:

- **Python LDSC h²** if you're replicating published results that used Python LDSC's
  chunked-window approximation. `--sketch 1600` matches within 0.001.
- **Per-SNP exact h²** (LDSC paper's mathematical definition) if you want the most
  defensible h² estimate. `--sketch 1600 --snp-level-masking` matches within 0.001,
  and is ~12% MORE ACCURATE than Python LDSC at biobank scale (Python's chunked-
  window approximation has 10–14% h² bias vs the paper definition).

Cross-validated against per-SNP exact, GCTA, chunk-exact, and Python LDSC at biobank
N=50K and at 1000G N=503. At d=1000 the combo is slightly faster (~19 s, ~19×) but
Height shifts ~0.003 h², still within the regression SE of ~0.008. At 1000G scale
(N=503), the d=480 combo reaches GCTA-quality LD scores (Pearson r=0.994 vs per-SNP
exact, lowest max-error of any method).

### `--mmap` — memory-mapped BED I/O

For HPC deployments with GPFS, Lustre, or other networked filesystems,
`--mmap` uses memory-mapped I/O instead of buffered reads:

- **Zero-copy access** for the fused CountSketch path (eliminates ~20GB of memcpy at biobank scale)
- **OS-managed readahead** via `MADV_SEQUENTIAL` — GPFS can prefetch across storage nodes in parallel
- **Async prefetch** via `MADV_WILLNEED` on the next chunk — no reader thread contention
- **No seek invalidation** — unlike `BufReader`, mmap'd pages stay resident once faulted

On local SSD with warm cache, `--mmap` regresses ~15% due to page fault overhead.
Use the default (no flag) for local storage. `--mmap` is designed for networked filesystems.

### `--gpu` (experimental, behind feature flag)

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
targets biobank-scale cohorts (n ≥ 50K) where each chunk's GEMM is large enough for
compute to dominate PCIe transfer. Incompatible with `--sketch`.

---

## Performance

Benchmarks on AWS c6a.4xlarge (AMD EPYC 7R13, 16 vCPU). Static musl binary
with mimalloc, AVX2+FMA target features. 1000G timings are hyperfine
(1 warmup + 3 timed runs); biobank V0.5.0 sketch+masking timings are
single-shot AWS Batch on-demand.

### LD-score computation (`l2`)

#### 1000G reference panel (N = 2,490, 1.66M SNPs)

| Mode | Wall time | vs Python | Accuracy |
|------|----------:|----------:|----------|
| Python LDSC | 25 min 49 s | 1.0× | reference |
| **Rust f64** (default) | **41.1 s** | **~38×** | exact (`max_abs_diff = 0` on chr22 EUR with `--python-compat`) |
| Rust f32 (`--fast-f32`) | ~32 s | ~48× | `max_abs_diff ≈ 0.008` chr22 EUR |
| Rust sketch (`--sketch 200`) | 25.4 s | ~61× | Pearson r ≈ 0.93 vs chunk-exact |
| Rust sketch (`--sketch 50`) | 15.3 s | ~101× | Pearson r ≈ 0.81 *(QC only)* |

#### Biobank scale (N = 50,000, 1.66M SNPs, h² accuracy on real GWAS)

Python LDSC at biobank measured on r6a.8xlarge (32 vCPU, 240 GiB, same
chunk-size 200 + ld-wind-kb 1000 as Rust); Rust modes measured on
c6a.4xlarge (16 vCPU, 28 GiB).

Python LDSC and Rust chunked-exact produce **bit-identical** h² values
at biobank (both use the same chunked-window approximation; the algorithm
matches, only the language differs). Both diverge from per-SNP exact —
the LDSC paper's mathematical definition — by 10–14% h² because of
chunked-window bias.

| Mode | Wall time | vs Python | h² \|Δ\| vs per-SNP exact (BMI 2010 / BMI 2018 / Height 2018) | h² \|Δ\| vs Python LDSC |
|------|----------:|--------:|---|---|
| Python LDSC | 6,541 s (109 min) | 1.0× | 0.015 / 0.020 / 0.047 (chunked-window bias) | 0 (reference) |
| `--fast-f32` (chunked-exact, ≡ Python h²) | 358 s | 18× | 0.015 / 0.020 / 0.047 | **0 / 0 / 0** (bit-identical to Python) |
| `--snp-level-masking --fast-f32` (per-SNP exact, "math truth") | 361 s | **18×** | 0 (truth) | 0.015 / 0.020 / 0.047 |
| **`--sketch 1600`** (matches Python h² within 0.001) | **21 s** | **311×** | 0.015 / 0.020 / 0.047 | **0.0001 / 0.0001 / 0.0008** |
| **`--sketch 1600 --snp-level-masking`** (matches per-SNP exact within 0.001) | **21 s** | **311×** | **0.0002 / 0.0003 / 0.0005** | 0.015 / 0.020 / 0.047 |
| `--sketch 1000` (chunked) | 18 s | 363× | 0.015 / 0.020 / 0.044 | 0.0004 / 0.0005 / 0.0031 |
| `--sketch 200` (default, fastest) | 15 s | 436× | 0.013 / 0.013 / 0.009 | 0.002 / 0.008 / 0.038 |

Cross-method h² validated on BMI 2010, BMI 2018, Height 2018 against
per-SNP exact, GCTA, Python LDSC, and chunk-exact references. Full tables
in [`docs/perf-log.md`](docs/perf-log.md) 2026-05-12 entry.

Fused CountSketch reads packed BED bytes and scatter-adds directly into
the d×c sketch buffer, eliminating the N×c intermediate entirely. Cost is
O(N×c) independent of d — d=200 has the same scatter time as d=1000;
larger d only matters when d² × c starts contributing meaningfully (≈ d=2000+
at biobank). A 28 GB container is recommended for N=50K.

### Window-size sensitivity (200k-SNP extract)

| Window | Python | Rust f64 | Speedup |
|--------|--------|----------|---------|
| `--ld-wind-kb 100` | 44.2 s | 5.1 s | **8.7×** |
| `--ld-wind-kb 500` | 48.4 s | 6.2 s | **7.8×** |
| `--ld-wind-kb 1000` | 53.7 s | 8.6 s | **6.2×** |
| `--ld-wind-kb 2000` | 61.8 s | 12.9 s | **4.8×** |

### Scaling

The ring-buffer algorithm keeps memory bounded by the LD window size, not
total SNP count.

**Scaling with M (fixed N=2,490, exact-f64):**

| M (SNPs) | Est. wall time (AWS EPYC 16v) | BED size | Peak memory |
|----------|-------------------------------|----------|-------------|
| 1.66M | 41 s *(measured)* | 1 GB | ~100 MB |
| 5M | ~124 s | 3 GB | ~100 MB |
| 10M | ~247 s | 6 GB | ~100 MB |
| 50M | ~1,235 s (~21 min) | 30 GB | ~100 MB |

**Scaling with N (fixed M=1.66M, full-genome):**

| N (individuals) | Mode | Wall time |
|-----------------|------|-----------|
| 2,490 (1000G) | exact-f64 | 41.1 s *(measured)* |
| 2,490 (1000G) | `--sketch 200` | 25.4 s *(measured)* |
| 50,000 (biobank) | exact-f32 (`--fast-f32 --global-pass`) | 361 s *(measured)* |
| 50,000 (biobank) | `--sketch 1600 --snp-level-masking --mmap` | 21 s *(measured)* |
| 50,000 (biobank) | `--sketch 200 --mmap` | 15 s *(measured)* |

At biobank N, GEMM cost is O(N × w × c) per chunk; CountSketch reduces this
to O(N×c) scatter-add, independent of d. BED I/O is sequential and
throughput-bound; `--fast-f32` halves BED read bytes per chunk.

### Other subcommands (UKBB I/O, Apple M4 / 10-core)

| Workflow | Rust | Python | Speedup |
|---------|------|--------|---------|
| munge-sumstats | 3.74 s | 62.65 s | **16.75×** |
| h2 | 0.90 s | 7.81 s | **8.68×** |
| rg (two traits) | 2.93 s | 28.09 s | **9.59×** |

(Dataset: 497 MB Pan-UKBB sumstats input; UKBB.EUR LD scores, 381,831 SNPs
after merge.)

---

## Build and install details

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

- `--rayon-threads N`: Rayon thread count for jackknife in `h2`/`rg`.
- `--polars-threads N`: Polars thread count for CSV streaming in `munge-sumstats`.
- `--mmap`: Memory-mapped BED I/O for `l2`. Recommended for GPFS/Lustre HPC.

---

## Further Reading

- **[CountSketch math analysis](docs/countsketch-math-analysis.md)** —
  rigorous derivation of the CountSketch bias on the renormalized cosine,
  Monte Carlo validation, and the §15 sketch+masking finding.
- **[Performance deep-dive](docs/performance-deep-dive.md)** — algorithmic
  complexity for each mode, scaling analysis for dense SNP panels
  (O(M² × N) with distance-based windows), downstream h2/rg regression
  accuracy by sketch dimension, and why Python is slow.
- **[Performance log](docs/perf-log.md)** — append-only log of every perf
  experiment; the 2026-05-12 entry has the cross-method h² validation
  (GCTA + Python LDSC + chunk-exact + per-SNP exact + sketch+masking).
- **[Architecture & source map](docs/architecture.md)** — module-level code
  map, key data-flow invariants, and dependency rationale.
- **[GCTA source audit](docs/gcta-source-audit.md)** — what GCTA's
  `--ld-meanrsq` actually does, since GCTA is the tool the original LDSC
  paper used to compute the published 1000G reference scores.

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

## How to Cite

If you use ldsc-rs in your research, please cite:

> Haason, S. & Khan, Y. (2026). ldsc-rs: Exact and approximate LD Score Regression at biobank scale. *bioRxiv*. <!-- TODO: add DOI after submission -->

```bibtex
@article{haason2026ldscrs,
  author = {Haason, Sharif and Khan, Yousef},
  title = {ldsc-rs: Exact and approximate {LD} {Score} Regression at biobank scale},
  journal = {bioRxiv},
  year = {2026},
}
```
