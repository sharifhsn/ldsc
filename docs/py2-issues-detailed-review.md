# LDSC (Python 2 repo) Issues — Detailed Review (Comments + Recency Weighted)

This document summarizes the most discussed and/or recently active issues from the Python 2 LDSC repository (bulik/ldsc). It is intended to capture **what problems users actually hit**, the **common root causes**, and the **typical fixes** so you don’t need to read the full threads.

Selection heuristic:
- High comment count (discussion intensity)
- Recency (recently updated threads are prioritized even if comment count is lower)

## Cross‑Cutting Themes

### 1) Installation / environment breakage
- Python 2 is deprecated; conda channels no longer provide python=2.7 on newer architectures (notably macOS arm64).
- Many failures are due to mismatched or missing dependencies (numpy/pandas/bitarray) or confusion about running py2 code in py3.

### 2) Munging failures due to input schema problems
- Most munge errors trace to unexpected headers, non‑numeric values in numeric columns, or mis‑specified `--snp` / `--a1` / `--a2` flags.
- Several issues reflect misunderstandings around **A1/A2 semantics** and the direction of signed summary stats.

### 3) Reference LD / annotation file path problems
- Cell‑type analysis issues are frequently due to **relative paths in .ldcts** or missing LD score files.
- Index errors in `--ref-ld-chr` usually trace to missing or mis‑named LD score files.

### 4) Interpretability of partitioned h2 outputs
- Negative Prop._h2 for continuous annotations is not necessarily a bug; the meaningful quantity is the coefficient.

### 5) Non‑EUR / build38 reference resources
- People struggle to obtain the right HapMap3 list, LD scores, or weights for non‑EUR and hg38 datasets.

---

## Issue‑by‑Issue Summaries

### Issue #26 — munge_sumstats.py (ImportError: cannot import name sumstats)
**Problem**: Fresh users got `ImportError: cannot import name sumstats` when running `munge_sumstats.py` from the tutorial.
**Root cause**: Incomplete/incorrect repo checkout (missing `ldscore/sumstats.py` or `__init__.py`).
**Common fix**: Re‑clone the repo; ensure you run `ldsc.py -h` from the root and that `ldscore/` contains expected files.
**Notes**: Thread also includes early discussion about minor rg differences due to code updates and APoE exclusion in Alzheimer’s example.

### Issue #37 — median value issue for munge_sumstats.py
**Problem**: `ValueError: median value of SIGNED_SUMSTATS ... should be close to 0`.
**Root cause**: Signed statistic column not actually centered at the null (often unstandardized effect sizes). LSDC expects signed stats consistent with null value (0 or 1).
**Fix**: Standardize betas (e.g., divide by SD) or use a correctly signed Z column; verify the signed stat and null value in `--signed-sumstats`.

### Issue #369 — module parse has no attribute (py3 port confusion)
**Problem**: `AttributeError: module 'parse' has no attribute 'sumstats'` when running in Python 3.
**Root cause**: Broken/partial py3 port or mixed imports (local `parse.py` shadowing `ldscore.parse`).
**Fix**: Use py2 environment, or fully convert code to py3 and fix imports; ensure `import ldscore.parse as ps` is used consistently.
**Notes**: Many users attempt ad‑hoc fixes with pip‑installed “parse”/“regressions” modules, which is incorrect.

### Issue #425 — Is `--merge-alleles` required for rg?
**Problem**: rg fails with `KeyError: Incompatible alleles ...` despite intersecting SNPs.
**Root cause**: Indels or allele mismatch despite SNP intersection. ldsc’s allele matching fails on indels and strand ambiguous variants.
**Fix**: Remove indels; run `munge_sumstats.py` with `--merge-alleles` consistently. Ensure allele harmonization is correct.
**Extra**: Thread also discusses sample overlap. LDSC does not estimate overlap; intercept can be computed but that’s not the same as overlap.

### Issue #374 — Could not find SNP column
**Problem**: `ValueError: Could not find SNP column`.
**Root cause**: `--snp` argument doesn’t match the actual header or includes quotes/special characters; header mismatch.
**Fix**: Ensure the column name matches exactly and is not quoted; verify with `head` and the “Interpreting column names” section.

### Issue #43 — fatal error with median beta in munge
**Problem**: median beta far from 0 triggers error.
**Root cause**: Summary stats are not centered at the null, often because the file is filtered to significant SNPs or betas aren’t standardized.
**Fix**: Use full GWAS, convert to Z scores or use correct signed stat; use `--a1-inc` if appropriate. Check allele columns for missing values.

### Issue #366 — IndexError: list index out of range (ref LD)
**Problem**: `IndexError` in ldscore_fromlist when loading ref LD.
**Root cause**: ldsc cannot find any `*.l2.ldscore.gz` files for the given prefix; likely wrong path or missing files.
**Fix**: Verify `--ref-ld-chr` path and that files exist for chromosomes. Users resolved by downloading correct eur_w_ld_chr files.

### Issue #66 — Error converting summary statistics (No objects to concatenate)
**Problem**: `ValueError: No objects to concatenate` during munge.
**Root cause**: Data parsing yields empty chunks, often due to header/data mismatch or invalid values in numeric columns.
**Fix**: Validate header length vs row length; ensure numeric columns have no stray strings; test smaller subsets.

### Issue #145 — munge_sumstats takes hours / appears stuck
**Problem**: `munge_sumstats.py` extremely slow on large files.
**Root cause**: Default chunk size (5,000,000) causes heavy memory pressure and slow pandas operations.
**Fix**: Reduce `--chunksize` (e.g., 500,000). This commonly reduces runtime from hours to minutes.

### Issue #447 — S‑LDSC with `--ref-ld-chr-cts` questions
**Problem**: Confusion about output (only coefficients) and interpretation of control annotations.
**Key clarifications**: `h2‑cts` reports coefficients only; prop.h2/enrichment for CTS is not provided without code changes. Control annotations are included as covariates to adjust the model; coefficients are conditional on them.
**Takeaway**: This is about model interpretation, not a bug.

### Issue #349 — cell type specific analyses index error
**Problem**: `IndexError: list index out of range` when running `--h2-cts`.
**Root cause**: Relative paths in `.ldcts` not resolving because `ldsc.py` run from a different directory.
**Fix**: Run ldsc from the directory containing the `.ldcts`, or change the `.ldcts` paths to absolute.

### Issue #178 — Anaconda3 does not work
**Problem**: Users assume LDSC works in Python 3; conda env fails; confusion about py2/py3.
**Root cause**: py2 dependencies and numpy pinning; Python 2 deprecation.
**Fix**: Use a Python 2 environment (historically), or port to Python 3 / use a maintained fork.

### Issue #147 — `Cahoy.1.l2.ldscore` cannot be opened
**Problem**: Cell‑type analyses fail with missing LD score file errors.
**Root cause**: `.ldcts` uses relative paths; user runs from wrong location or moved files.
**Fix**: Provide correct `--ref-ld-chr` list or run from directory expected by `.ldcts`.

### Issue #432 — python 2.7 no longer available in conda
**Problem**: conda can’t install python=2.7 on arm64; `PackagesNotFoundError`.
**Root cause**: Py2 is deprecated; many channels removed it.
**Fix**: Use a Python 3 port/fork, or run on Linux x86 where older environments still exist.

### Issue #438 — Negative Prop._h2 in partitioned h2
**Problem**: Users see significant negative Prop._h2 values.
**Root cause**: Prop._h2 is not meaningful for continuous annotations; coefficients are the interpretable metric.
**Fix**: Use `--print-coefficients` and interpret coefficients. Ensure `--overlap-annot` is used for baseline models.

### Issue #470 — LDSC on build38 data
**Problem**: Users want to run LDSC on hg38 without liftover.
**Root cause**: Official ref LD scores and annotations are hg19. Users must rebuild LD scores, weights, and frequency files for hg38.
**Fix**: Build hg38 annotations and LD scores from hg38 reference genotypes; use published hg38 resources when available.

### Issue #407 — Suitable w_hm3.snplist for AFR
**Problem**: Users need non‑EUR hapmap3 lists and weights.
**Root cause**: Official `w_hm3.snplist` is typically EUR‑centric; non‑EUR populations need custom maps and weights.
**Fix**: Use HapMap3 MAP files for target populations, lift over as needed, and compute matching LD scores/weights.

### Issue #439 — Filtering of munge_sumstats / A1/A2 semantics
**Problem**: Confusion about why ambiguous SNPs are removed and whether A1 is ref or effect allele.
**Resolution**: Ambiguous SNPs are removed to avoid strand errors. A1 is assumed to be the reference allele for signed stats; beta must be oriented accordingly. Users must ensure A1 alignment with LD reference assumptions.

### Issue #322 — HapMap3 variants / genome build
**Problem**: Users want updated HapMap3 lists and correct build mapping.
**Root cause**: HapMap3 is old (hg18). Liftover loses variants; tools use custom subsets.
**Fix**: Liftover with care or use existing curated lists from large pipelines (pan‑UKBB, etc.). No definitive upstream guidance provided.

### Issue #373 — PyPI version vs GitHub version mismatch
**Problem**: PyPI lists a py3 version (2.0.1) not present on GitHub.
**Root cause**: PRs for py3 support were not merged into the main repo.
**Fix**: Use maintained forks or external ports (e.g., cbiit/ldsc).

### Issue #493 — make_annot.py duplicate output rows
**Problem**: Annotation output has duplicate rows when multiple SNPs share the same BP.
**Root cause**: make_annot assumes per‑chromosome input and unique BP positions; collisions create duplicates on merge.
**Fix**: Split by chromosome and ensure unique positions; this is expected behavior per maintainer.

### Issue #190 — make_annot produces all‑zero annotations
**Problem**: Gene‑set annotations are all zero.
**Root cause**: Gene coordinate file mismatch (chrom naming, build mismatch, or coordinate format) or windowsize too small.
**Fix**: Verify gene coordinate file matches build and chromosome naming; test with known positive controls; consider alternative coordinate sources.

### Issue #329 — `ValueError: could not convert string to float: OR`
**Problem**: Numeric columns contain string values (e.g., OR labels), causing conversion failures.
**Root cause**: Dirty data or unexpected non‑numeric values in numeric columns.
**Fix**: Validate column types before munging; remove or correct invalid values; test smaller subsets.

### Issue #133 — A1/A2 allele interpretation confusion
**Problem**: Documentation says A1 is effect allele, but logs interpret A1 as reference for signed stats, causing confusion.
**Resolution**: A1 is the allele to which the signed stat corresponds. If your effect allele is alt, supply A1=alt and ensure sign matches. The term “reference allele” is overloaded and must be interpreted consistently with signed stats.

---

## Practical Takeaways for a Rust Rewrite
- Provide explicit, user‑facing validation of paths (LD scores, `.ldcts`) and file existence to avoid “IndexError” confusion.
- Clearly document A1/A2 expectations and signed stats alignment; consider auto‑detection warnings for mis‑centered signed stats.
- Expose `--chunksize` (or equivalent) in any Rust munging pipeline to avoid slow processing of huge files.
- Improve error messages for missing LD score files and invalid numeric columns.
- Support better diagnostics for allele mismatches and indel handling during rg.

---

If you want, I can also produce a grouped “taxonomy” version of this doc (installation, munging, LD scores, h2/rg, annotations, reference data) with recommendations for Rust parity fixes.
