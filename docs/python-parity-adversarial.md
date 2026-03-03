# Adversarial Parity Review: Python LDSC (ldsc_py2/ldsc_py3) vs Rust Rewrite

This document is an **adversarial parity pass** comparing the Rust rewrite to the original Python LDSC codebases located at:
- `/Users/sharif/Code/ldsc/ldsc_py2`
- `/Users/sharif/Code/ldsc/ldsc_py3`

It incorporates **known parity gaps** from existing docs and **additional gaps discovered** by inspecting Python sources.

Sources referenced:
- Python CLI + logic: `/Users/sharif/Code/ldsc/ldsc_py3/ldsc.py`, `/Users/sharif/Code/ldsc/ldsc_py3/ldscore/sumstats.py`, `/Users/sharif/Code/ldsc/ldsc_py3/ldscore/regressions.py`, `/Users/sharif/Code/ldsc/ldsc_py3/munge_sumstats.py`
- Rust CLI + logic: `/Users/sharif/Code/ldsc/src/cli.rs`, `/Users/sharif/Code/ldsc/src/regressions.rs`, `/Users/sharif/Code/ldsc/src/munge.rs`, `/Users/sharif/Code/ldsc/src/l2.rs`, `/Users/sharif/Code/ldsc/src/jackknife.rs`
- Existing parity notes: `/Users/sharif/Code/ldsc/docs/feature-parity.md`

---

## Summary of Major Parity Gaps (High Risk)

### 1) `rg` Pairing Semantics (Python: trait1 vs all; Rust: adjacent pairs)
**Python behavior**: In `sumstats.estimate_rg`, Python computes rg between the **first trait and every other trait** (`p1` vs `p2..pn`).
- See: `/Users/sharif/Code/ldsc/ldsc_py3/ldscore/sumstats.py` (`estimate_rg` iterates `rg_paths[1:n_pheno]`).

**Rust behavior**: Rust loops over `args.rg.windows(2)` and computes **only adjacent pairs** (A‑B, B‑C, C‑D...).
- See: `/Users/sharif/Code/ldsc/src/regressions.rs` (`for (pair_idx, pair) in args.rg.windows(2)`)

**Impact**: Results differ fundamentally for any `--rg` list of length >2. This is a **functional mismatch**.

---

### 2) Default Two‑Step Behavior for h2/rg
**Python behavior**: If scalar regression (`n_annot == 1`) and neither `--two-step` nor intercepts are set, Python **defaults to `two_step = 30`** for both h2 and rg.
- See: `/Users/sharif/Code/ldsc/ldsc_py3/ldscore/sumstats.py` (`estimate_h2` and `estimate_rg` set `args.two_step = 30` when conditions met).

**Rust behavior**: Rust uses two‑step **only when explicitly provided**.
- See: `/Users/sharif/Code/ldsc/src/regressions.rs`

**Impact**: Default numerical behavior diverges from Python for typical scalar h2/rg runs.

---

### 3) Partitioned h2 Uses Old Weights Path in Python, IRWLS in Rust
**Python behavior**: When `n_annot > 1`, `old_weights=True` is passed to `Hsq`, which triggers a **different fitting path** (weighted least squares + jackknife) and bypasses IRWLS.
- See: `/Users/sharif/Code/ldsc/ldsc_py3/ldscore/sumstats.py` (`old_weights=True` when `n_annot > 1`).
- See: `/Users/sharif/Code/ldsc/ldsc_py3/ldscore/regressions.py` (`old_weights` path uses `LstsqJackknifeFast`).

**Rust behavior**: Partitioned h2 always uses IRWLS (`jackknife::jackknife` → `irwls`).
- See: `/Users/sharif/Code/ldsc/src/regressions.rs`

**Impact**: Potentially significant differences in partitioned h2 estimates and SEs.

---

### 4) N‑Scaling in Design Matrix (Per‑SNP N)
**Python behavior**: In LDSC regressions, the design matrix is scaled by `N / Nbar` per SNP (`x = N * L2 / Nbar`).
- See: `/Users/sharif/Code/ldsc/ldsc_py3/ldscore/regressions.py` (`Nbar` scaling in `LD_Score_Regression`).

**Rust behavior**: Design matrix uses raw `L2` only; `N` is used in weights, and final h2 uses `N_mean`.
- See: `/Users/sharif/Code/ldsc/src/regressions.rs` (`x[[i,0]] = ref_l2[i]`)

**Impact**: If N varies across SNPs, results will diverge.

---

### 5) `rg` Partitioned LD Scores (Multi‑annotation)
**Python behavior**: `rg` supports multiple LD score columns (`ref_ld_cnames`) as predictors.
- See: `/Users/sharif/Code/ldsc/ldsc_py3/ldscore/sumstats.py` and `regressions.py`.

**Rust behavior**: `rg` loads only a scalar `L2` column and ignores partitioned LD scores.
- See: `/Users/sharif/Code/ldsc/src/regressions.rs` (`load_ld` selects only `L2`).

**Impact**: Partitioned rg is not supported in Rust.

---

### 6) `n_blocks` Capping
**Python behavior**: `n_blocks = min(n_snp, args.n_blocks)` for both h2 and rg.
- See: `/Users/sharif/Code/ldsc/ldsc_py3/ldscore/sumstats.py`.

**Rust behavior**: Uses `args.n_blocks` directly. If `n_blocks > n_obs`, block size becomes 0 and can break jackknife logic.
- See: `/Users/sharif/Code/ldsc/src/jackknife.rs` (block size = `n / n_blocks`).

**Impact**: Potential runtime errors or invalid SEs on small datasets.

---

## Additional Parity Gaps (Medium/Low Risk)

### 7) Sumstats delimiter handling (Python supports whitespace; Rust tab‑only)
**Python**: `munge_sumstats.py` reads with `sep=r'\s+'` (whitespace‑delimited).
**Rust**: `scan_sumstats` uses separator `\t` only.
- Python: `/Users/sharif/Code/ldsc/ldsc_py3/munge_sumstats.py`
- Rust: `/Users/sharif/Code/ldsc/src/parse.rs`

**Impact**: Whitespace‑delimited sumstats parse in Python but likely fail in Rust.

---

### 8) Out‑of‑bounds P/INFO/FRQ handling
**Python**: Explicitly filters and warns on out‑of‑bounds `P`, `INFO`, `FRQ`.
**Rust**: No explicit P/INFO/FRQ out‑of‑range checks; may accept invalid values silently.
- Python: `/Users/sharif/Code/ldsc/ldsc_py3/munge_sumstats.py` (`filter_pvals`, `filter_info`, `filter_frq`).
- Rust: `/Users/sharif/Code/ldsc/src/munge.rs` (no corresponding checks).

---

### 9) Default `chisq_max` for partitioned h2
**Python**: For partitioned h2, default `chisq_max = max(0.001 * Nmax, 80)` if not specified.
**Rust**: No default for partitioned h2 (except in `h2-cts`).

---

### 10) `--maf` default behavior (pre‑filter)
**Python**: `--maf` filters SNPs before LD computation.
**Rust**: Default now matches Python (prefilter).

---

### 11) `--annot` + `--extract` incompatibility
**Python**: Raises error if `--annot` and `--extract` are used together.
**Rust**: Only warns and proceeds.

---

### 12) Output artifacts (.log, .cov, .delete)
**Python**: Writes `.log`, `.cov`, `.delete`, `.part_delete` files.
**Rust**: Prints covariance/delete values to stdout and does not generate these files.

---

### 13) `--pickle` output format
**Python**: Supports `--pickle` for LD score output.
**Rust**: Not implemented.

---

### 14) LD score output columns
**Known gap from parity doc**: Python `.l2.ldscore` includes `CM` and `MAF`; Rust outputs only `CHR SNP BP` + L2 columns.
- See: `/Users/sharif/Code/ldsc/docs/feature-parity.md`

---

## Known Parity Gaps Already Documented (For Completeness)

From `/Users/sharif/Code/ldsc/docs/feature-parity.md`:
- LD score output format lacks `CM`/`MAF` columns.
- `l2 --no-print-annot` is accepted but no‑op.
- Missing explicit filtering for indels and out‑of‑range P values in `munge`.

---

## Additional Potential Differences Worth Verifying

These are not fully proven but merit targeted validation:

1. **Allele mismatch filtering logic**: Rust `align_rg_alleles` uses direct base complement logic; Python uses `MATCH_ALLELES` lookup table. Likely equivalent but should be validated on edge cases with allele case/encoding differences.
2. **IRWLS vs LstsqJackknifeFast** in partitioned h2 may cause significant SE differences.
3. **N scaling**: differences grow with heterogeneous N across SNPs.
4. **Whitespace input**: Rust claims whitespace support in docs but parser is tab‑only.

---

## Recommended Parity Fixes (If Exact Parity Is Required)

1. Change `rg` pairing to match Python (trait1 vs all others). Consider also adding a flag for adjacent pairs if needed.
2. Apply Python default `--two-step 30` for scalar h2/rg when no intercept is fixed.
3. Add Python’s `old_weights` path for partitioned h2.
4. Scale design matrix by `N / Nbar` as in Python regression.
5. Cap `n_blocks` at `n_obs`.
6. Support whitespace‑delimited sumstats input or auto‑detect delimiter.
7. Restore `.log`, `.cov`, `.delete` outputs for compatibility.
8. Add optional `--pickle` output support if required by downstream tooling.

---

## Closing Notes
This pass intentionally assumes **parity must be exact**. Several high‑impact gaps were found, especially in `rg` pairing logic, default two‑step behavior, partitioned weighting, and N‑scaling. Those should be prioritized if strict parity is required.
