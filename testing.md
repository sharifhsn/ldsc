# Testing Notes: Rust LDSC vs Python (Pan-UKBB sample)

Date: 2026-03-01

This document captures the exact steps, scripts, environment setup, and current findings from the Rust vs Python LDSC performance/accuracy comparison using Pan-UKBB data on a Macbook Air.

## Key Artifacts

- Pan-UKBB phenotype file:
  - `/Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.tsv.bgz`
- UKBB LD scores (downloaded):
  - `/Users/sharif/Code/ldsc/data/UKBB.ALL.ldscore/`
- Preparation script (Polars + uv script deps):
  - `/Users/sharif/Code/ldsc/scripts/panukbb_prepare_sumstats.py`
- Benchmark script:
  - `/Users/sharif/Code/ldsc/scripts/bench_ukbb_py3_vs_rust.sh`
- Pan-UKBB summary doc:
  - `/Users/sharif/Code/ldsc/docs/pan-ukbb-summary.md`

## Pan-UKBB Specifics That Matter

- The phenotype file uses `neglog10_pval_*` columns (not raw P). We convert `P = 10^(-neglog10_pval)`.
- SNP ID scheme in the phenotype file: `chr`, `pos`, `ref`, `alt` → we construct `SNP = chr:pos:ref:alt`.
- Sample size (EUR) from phenotype manifest metadata:
  - `n_cases_EUR = 367192` (used as N in prep step).

## Python 3.9 Environment with uv (per ldsc_py3/environment3.yml)

The Python implementation is run from `/Users/sharif/Code/ldsc/ldsc_py3`. It expects Python 3.9+ and pinned libs.

### Create clean venv
```bash
cd /Users/sharif/Code/ldsc/ldsc_py3
uv venv --clear --python 3.9
```

### Install build helpers and deps
```bash
uv pip install setuptools==65.5.1 wheel Cython==0.29.36
uv pip install numpy==1.21.5
uv pip install --no-build-isolation \
  pandas==1.3.3 scipy==1.7.3 bitarray==2.0.0 \
  python-dateutil==2.8.2 pytz==2022.4 six==1.16.0
```

### Verify
```bash
uv run python -V
# Expect: Python 3.9.6
```

## Prepare LDSC-ready TSV from Pan-UKBB

Uses Polars for speed (streaming, progress logging). Output is a plain TSV (not gz) to avoid Python LDSC binary/text mode issues.

```bash
uv run /Users/sharif/Code/ldsc/scripts/panukbb_prepare_sumstats.py \
  --input /Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.tsv.bgz \
  --output /Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.ldsc.tsv \
  --mode pop \
  --pop EUR \
  --n 367192 \
  --drop-low-confidence \
  --log-every 2000000
```

Produces:
- `/Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.ldsc.tsv`

## Create a Larger Sample (avoid <200k SNP warning)

The goal is a larger slice that still runs on a Macbook Air but yields
>200k SNPs after LD merge, so the Python warning does not appear.

```bash
head -n 4000001 \
  /Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.ldsc.tsv \
  > /Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.sample.tsv
```

## Munge Sumstats (Python LDSC)

```bash
cd /Users/sharif/Code/ldsc/ldsc_py3
uv run python /Users/sharif/Code/ldsc/ldsc_py3/munge_sumstats.py \
  --sumstats /Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.sample.tsv \
  --out /Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.sample
```

Outputs:
- `/Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.sample.sumstats.gz`

Observed in munge:
- Input rows: 4,000,000
- Remaining after QC: 3,128,466

## Benchmark Script

The benchmark script handles:
- Python env setup (if missing)
- `cargo clippy --release` and `cargo build --release`
- timing of Python vs Rust runs

Run:
```bash
SUMSTATS=/Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.sample.sumstats.gz \
POP=EUR \
ID_SCHEME=chrpos \
LDMS=none \
TWO_STEP=30 \
/Users/sharif/Code/ldsc/scripts/bench_ukbb_py3_vs_rust.sh
```

## Correctness Parity Script (H2)

Runs adversarial correctness checks against Python LDSC using small synthetic datasets
(unsorted LD order, varying N with L2 < 1 clamp, and chisq-max boundary).

```bash
/Users/sharif/Code/ldsc/scripts/verify_parity_h2.sh
```

Notes:
- Uses small `--n-blocks` values to avoid Python errors on tiny datasets.
- Compares key summary numbers with a small tolerance (5e-4).

## Future Work: Fuzzing

Goal: stress Python vs Rust parity with adversarial inputs beyond the hand-crafted cases.

Ideas:
- Random SNP ordering, duplicated SNP IDs, and missing SNPs between sumstats/LD scores.
- NaN/inf injection in Z/N/L2 columns and boundary values (e.g., chisq_max edges).
- Property-based checks in Rust (e.g., proptest) with scripted parity runs against Python.
- Longer-term: `cargo-fuzz` harnesses for merge/parsing logic and regression stability.

### Important Differences: Python vs Rust flags
- Python `ldsc.py --ref-ld` expects a **prefix** (no `.l2.ldscore.gz`). The script strips suffixes for Python automatically.
- Rust uses the full file path.
- Rust now reads `M` from `*.l2.M_5_50` for single-file `--ref-ld` and passes `--M` explicitly (script logic).
- `TWO_STEP=30` is passed to both.

## Results (Larger Sample)

From `perf/ukbb/*.time`:
- Python: `real 4.93s`
- Rust: `real 2.78s`

Accuracy:
- Python (LDSC):
  - h2 = 0.1393 (SE 0.0161)
  - intercept = 1.2134
- Rust:
  - h2 = 0.1076 (SE 0.0149)
  - intercept = 1.3345

Both used:
- 201,908 SNPs after LD merge
- M = 7,009,236 (from `UKBB.EUR.l2.M_5_50`)
- two-step estimator with cutoff 30

### Conclusion
Performance comparison works and the <200k warning is gone, but outputs are still **not identical**. Remaining difference likely from two-step/intercept implementation details. Next steps:

1. Force the Python intercept into Rust via `--intercept-h2` for a strict alignment test.
2. Audit Rust two-step implementation against Python (`ldscore/regressions.py`) to achieve bitwise-equivalent output.

## How to Scale Up Sample Size Gradually

Repeat the sample creation + munge + benchmark with larger slices:

```bash
# Example: 1,000,000 SNPs + header
head -n 1000001 /Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.ldsc.tsv \
  > /Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.sample.tsv

cd /Users/sharif/Code/ldsc/ldsc_py3
uv run python /Users/sharif/Code/ldsc/ldsc_py3/munge_sumstats.py \
  --sumstats /Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.sample.tsv \
  --out /Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.sample

SUMSTATS=/Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.sample.sumstats.gz \
POP=EUR ID_SCHEME=chrpos LDMS=none TWO_STEP=30 \
/Users/sharif/Code/ldsc/scripts/bench_ukbb_py3_vs_rust.sh
```

## Known Issues / Fixes Applied

- Python LDSC expects `--ref-ld` prefix (no `.l2.ldscore.gz`), so the script strips suffixes for Python only.
- Rust single-file `--ref-ld` lacked M sidecar auto-detection; the benchmark script now reads `*.l2.M_5_50` and passes `--M`.
- Python 3.9 venv is required for ldsc_py3 dependencies.
- Python `read_csv` no longer emits "compression has no effect" warnings for gzipped sumstats.

## Files Generated During This Process

- `/Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.ldsc.tsv`
- `/Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.sample.tsv`
- `/Users/sharif/Code/ldsc/data/biomarkers-30600-both_sexes-irnt.sample.sumstats.gz`
- `/Users/sharif/Code/ldsc/perf/ukbb/py3_h2.time`
- `/Users/sharif/Code/ldsc/perf/ukbb/rust_h2.time`
