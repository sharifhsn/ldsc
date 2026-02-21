# Feature Parity Notes (Python LDSC vs Rust)

This document captures why certain Python LDSC features remain unimplemented in the Rust
port and the expected effort/impact if we choose to add them later.

Run date: 2026-02-21

## Unimplemented Features and Rationale

### Overlapping annotation partitioned h2 (`--overlap-annot` + `--frqfile[-chr]`)

**Python behavior**
1. Reads the `.annot` matrices (possibly multiple prefixes) for all SNPs.
2. Applies MAF filtering using `.frq` files (MAF > 0.05).
3. Computes the annotation overlap matrix `XᵀX` across all categories.
4. Uses that overlap matrix to compute `Prop._h2`, enrichment, and related statistics.

**Rust status**
- Not implemented. Rust currently uses `.l2.M` / `.l2.M_5_50` and does not read
  per-annotation overlap matrices.

**Why it is non-trivial**
- Requires reading all `.annot` matrices for each chromosome and computing `XᵀX`
  (dense matrix) across all annotations. This is heavy I/O and memory, especially
  for baseline models with many categories.
- Requires integrating per-chromosome `.frq` handling and MAF filtering identical to
  Python for reproducible overlap counts.
- Requires extending the regression output to mirror Python’s overlap-specific results.

**Impact**
- Moderate-to-high complexity and potentially large performance regressions for large
  annotation sets. Needs careful design to avoid prohibitive memory use.

### Cell-type specific heritability (`--h2-cts`, `--ref-ld-chr-cts`)

**Python behavior**
- Reads a `.ldcts` file containing multiple cell-type annotations.
- For each line, merges baseline LD scores + cell-type LD scores and runs an h2 regression.
- Outputs a `.cell_type_results.txt` file with coefficient stats.

**Rust status**
- Not implemented.

**Why it is non-trivial**
- Requires per-line assembly of design matrices and repeated regression (many runs).
- Requires cross-checks for SNP coverage between baseline and CTS LD scores.
- Depends on the overlap/annotation machinery for correct reporting in standard
  workflows (often used alongside `--overlap-annot`).

**Impact**
- Significant incremental complexity; also increases runtime due to many regressions.

## Decision

These features are deliberately left unimplemented for now because:
1. They require heavy additional I/O and dense matrix operations not present in the current
   Rust pipeline.
2. Correctness relies on subtle details of MAF filtering and annotation alignment.
3. Implementing them risks large performance regressions without careful optimization.

If priority changes, the recommended path is:
1. Implement `.frq` parsing + MAF filtering in `parse.rs`.
2. Implement `annot` overlap matrix computation in a streaming or blockwise fashion.
3. Extend `h2` output to match Python’s overlap and CTS results.
