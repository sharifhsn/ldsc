# Feature Parity Notes (Python LDSC vs Rust)

This document records parity status for features that historically differed from the
Python LDSC implementation.

Run date: 2026-02-21

## Implemented Parity

### Overlapping annotation partitioned h2 (`--overlap-annot` + `--frqfile[-chr]`)

**Python behavior**
1. Reads `.annot` matrices (possibly multiple prefixes) across all SNPs.
2. Applies MAF filtering using `.frq` files (keep `0.05 < FRQ < 0.95`).
3. Computes the annotation overlap matrix `XᵀX`.
4. Uses the overlap matrix to compute `Prop._h2`, enrichment, and related stats.

**Rust status**
- Implemented in `ldsc h2` when `--overlap-annot` is set.
- Writes `<out>.results` with overlap-aware proportions and enrichment.
- Uses `.frq` filtering identical to Python; requires `--frqfile` or `--frqfile-chr`
  unless `--not-m-5-50` is set.

### Cell-type specific heritability (`--h2-cts`, `--ref-ld-chr-cts`)

**Python behavior**
- Reads `.ldcts` (label + comma-separated LD score prefixes).
- For each line: merges baseline LD + CTS LD, runs regression, and reports the first
  coefficient.
- Writes `<out>.cell_type_results.txt`; optional `--print-all-cts` prints all CTS coefficients.

**Rust status**
- Implemented in `ldsc h2` via `--h2-cts`, `--ref-ld-chr-cts`, and `--print-all-cts`.
- Uses the same default `chisq_max` heuristic as Python for partitioned regressions.

### Continuous-annotation binning (`--cts-bin` workflow)

**Python behavior**
- Builds a one-hot annotation matrix by binning continuous values, then feeds it to LD score
  estimation.

**Rust status**
- Implemented as a separate preprocessing step: `ldsc cts-annot`.
- This generates a `.annot` file that is then used with `ldsc ldscore --annot`.
- Keeping it separate avoids changes to the LD score hot path.

## Additional Parity Flags

These previously missing flags are now implemented:
1. `ldscore --pq-exp` (generalized per-allele weighting). Output column names and `.M` files
   include the `_S{exp}` suffix, matching Python.
2. `munge-sumstats --daner` / `--daner-n` (Ripke daner formats). `--daner` infers `N_cas/N_con`
   from FRQ headers; `--daner-n` consumes `Nca/Nco` columns.
3. `rg --intercept-h2` (fix per-trait h2 intercepts in rg).
4. `ldscore --no-print-annot` is accepted for CLI parity and emits a warning (Rust does not
   generate `.annot` from `ldscore`).
