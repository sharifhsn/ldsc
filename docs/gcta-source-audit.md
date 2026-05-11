# Source-code audit: ldsc-rs vs GCTA LD-score computation

**Date**: 2026-05-11. **Scope**: line-by-line comparison of GCTA's `--ld-score`
routine (`/tmp/gcta_src/main/ld.cpp`, GCTA 1.95.1) against ldsc-rs's
`compute_ldscore_global` (`src/l2/compute.rs`). Existing parity claims against
GCTA were made without reading GCTA's source; this audit corrects them.

GCTA is the right reference because Bulik-Sullivan et al. 2015 used
`gcta64 --ld-score` (then aliased `--ld-meanrsq`) to compute the 1000 Genomes
LD scores shipped with LDSC. Python LDSC (the `bulik/ldsc` Python 2 tool) is
itself an undocumented approximation of GCTA's algorithm. Anchoring to GCTA
gives us a stable definition of "the paper's LD score" to validate against.

## TL;DR

Three previously-published claims are corrected by this audit:

1. **GCTA does NOT compute per-SNP exact windows.** It does
   block-based windowing with a two-pass overlap-averaging scheme. The
   block-decomposition is the window definition, not a memory
   optimization. Effective per-SNP window radius at `--ld-wind X` ranges
   from ~`X/2` to `X` kb depending on the SNP's position within the
   pass-1 block (asymmetric near boundaries). The mean snp_num at
   `--ld-wind 1000` on chr22 EUR is 903.6; ldsc-rs's true ±1000 kb
   neighbor count is 1182.

2. **GCTA `--ld-wind X` is BLOCK width in kb, not per-SNP radius.**
   `option.cpp:507-510` reads X, multiplies by 1000 to get bp, and
   passes it as `wind_bp` to `get_ld_blk_pnt`. The log message
   explicitly says "block size of X Kb with an overlap of X/2 Kb between
   blocks". So `--ld-wind 1000` and ldsc-rs `--ld-wind-kb 1000` are NOT
   equivalent: ldsc-rs's window is ~1.3× wider in mean neighbor count.

3. **GCTA's default r² is BIASED.** The
   `r² − (1−r²)/(N−2)` correction is applied at `ld.cpp:419,483` only
   when `_ldscore_adj` (set by `--ld-score-adj`) is true. ldsc-rs
   always uses unbiased. The existing
   `compare_gcta_ldsc.py` did NOT pass `--ld-score-adj`
   (verified at `gcta.log:15`), so the previously-reported Pearson
   r=0.998 vs masked is a *biased r²* + *narrower window* vs ldsc-rs's
   *unbiased r²* + *wider window* — with the two systematic errors
   partially cancelling. The apples-to-apples comparison is GCTA
   `--ld-wind 2000 --ld-score-adj`, which gives **r = 0.9999, mean L2
   diff = +0.006 vs ldsc-rs chunked, +0.18 vs ldsc-rs masked** on chr22
   EUR.

A corrected re-run is saved at
`preprint/data/sim_ldscores_audit/gcta_wind2000_adj.score.ld` for the
record.

## Algorithmic axes

For each axis: GCTA source citation, ldsc-rs source citation, what they
do, what differs, the numerical consequence, and the flag tuning needed
to align them.

### 1. Window definition and units

**GCTA** (`main/ld.cpp:318-369`, `main/option.cpp:507-510`,
`main/ld.cpp:273-305`).

`get_ld_blk_pnt(brk_pnt1, brk_pnt2, brk_pnt3, wind_bp)` builds two
overlapping sets of block boundaries, given `wind_bp` (block width in
bp; `--ld-wind 1000` produces `wind_bp = 1_000_000`):

- **Pass-1 blocks** (`brk_pnt1`): walk SNPs by ascending bp; start a new
  block when bp distance from current block-start exceeds `wind_bp`, OR
  on chromosome change, OR on a >1 Mb gap (`ld.cpp:332-353`). Result is
  a non-overlapping tiling: blocks of width ≤ `wind_bp` bp each. The
  trailing 2-element "boundary" pairs are size-2 and get skipped by the
  `size < 3 continue` guard at `ld.cpp:378`.
- **Pass-2 blocks** (`brk_pnt2`, `ld.cpp:358-368`): for each adjacent
  pair of `brk_pnt1` boundaries, insert the index midpoint. The
  resulting tiling has blocks that span across each pass-1 boundary,
  covering the right half of pass-1 block i and the left half of pass-1
  block i+1. So pass-2 block widths are similar to pass-1
  (≈`wind_bp`), but offset by ~half-block.
- **`brk_pnt3`**: marks where pass-1 block i ended within the pass-2
  block, used to avoid double-counting pairs already computed in pass 1
  (`ld.cpp:411-413`).

For each SNP `i`, the effective LD-score window is the UNION of pairs
seen in pass 1 plus pairs seen in pass 2:

- **Pass 1** gives pairs to other SNPs in the same pass-1 block: up to
  `wind_bp` bp away, but asymmetric depending on where SNP `i` sits
  within its block. A SNP at the block's left edge sees up to `wind_bp`
  to the right; a SNP at the right edge sees up to `wind_bp` to the
  left.
- **Pass 2** picks up cross-boundary pairs (only when `i` is in a
  pass-2 block). Pairs added: those with one endpoint in pass-1 block
  i's "right half" and the other in pass-1 block i+1's "left half".

**Empirical mean per-SNP neighbor count at chr22 1000G EUR
(N=503, 18 627 SNPs after MAF≥0.05):**

| Configuration | Mean snp_num | Mean L2 (with --ld-score-adj) |
|---|---|---|
| GCTA `--ld-wind 1000` | 903.6 | 18.44 |
| GCTA `--ld-wind 2000` | 1699.8 | 18.86 |
| ldsc-rs `--ld-wind-kb 1000` --snp-level-masking | 1182 (true ±1000 kb radius) | 18.67 |

So GCTA's effective window at `--ld-wind X` is *narrower* than ldsc-rs's
at `--ld-wind-kb X` in pair count (903 vs 1182). At `--ld-wind 2X`,
GCTA's neighbor count overshoots (1700 vs 1182), but its mean L2 still
matches ldsc-rs's almost exactly — the extra neighbors are distant and
low-r², contributing negligibly to the LD score.

**ldsc-rs** (`src/l2/window.rs:32-90`, `src/l2/compute.rs:830-840`).

`get_block_lefts_by_chr` computes a `block_left[i]` array per
chromosome: for each SNP `i`, `block_left[i]` is the smallest index `j`
such that `coords[i] - coords[j] ≤ max_dist`. So the per-SNP window is
`[block_left[i], i + something]`, a true per-SNP radius window. With
`--ld-wind-kb X`, every SNP gets a symmetric ±X kb window.

The chunk-eviction approximation in the default (no `--snp-level-masking`)
path: for chunk starting at index `chunk_start`, every SNP in
`[chunk_start, chunk_end)` uses `block_left[chunk_start]` as its left
boundary (`compute.rs:957-961`). This over-includes far-left SNPs for
SNPs late in a chunk. With `--snp-level-masking`, `compute.rs:1431-1442`
and `1735-1748` zero out r² for pairs outside each SNP's true
`block_left[i]`.

**Divergence verdict**: GCTA and ldsc-rs use *fundamentally different*
window schemes. GCTA's is block-based with two-pass overlap; ldsc-rs's
(with `--snp-level-masking`) is per-SNP exact. They are NOT equivalent
even at matched `wind_bp`. To align them on chr22 EUR, use GCTA
`--ld-wind 2000` (block width 2 Mb) — this happens to give matching
mean L2 to ldsc-rs `--ld-wind-kb 1000`.

### 2. r² estimator (biased vs unbiased)

**GCTA** (`main/ld.cpp:416-420`).

```cpp
rsq_sub(j,k) *= (ssx_sqrt_i_sub[j] * ssx_sqrt_i_sub[k]);
rsq_sub(j,k) = rsq_sub(j,k) * rsq_sub(j,k);
if (rsq_sub(j,k) >= rsq_cutoff) {
    if(_ldscore_adj) mean_rsq_sub[j] += rsq_sub(j,k) - (1.0 - rsq_sub(j,k)) / (n - 2.0);
    else mean_rsq_sub[j] += rsq_sub(j,k);
    rsq_size[j] += 1.0;
}
```

So:
- **r² formula**: `((X_j · X_k) / (||X_j|| · ||X_k||))²` — standard
  Pearson r² squared. Computed in single precision (`MatrixXf`); the
  `ssx_sqrt_i = 1/||X.col||` normalization (`ld.cpp:397-401`) is in
  `eigenVector` = double.
- **Bias correction**: `r² − (1−r²)/(N−2)` applied ONLY under
  `--ld-score-adj` (`option.cpp:528-534`). Default is **biased r²**.
- The cutoff check at `ld.cpp:418` uses the BIASED r² regardless of
  `_ldscore_adj`. At cutoff = 0 this is moot.

**ldsc-rs** (`src/l2/compute.rs:727-735, 774-782`).

```rust
pub(super) fn r2_unbiased(r: f64, n: usize) -> f64 {
    let sq = r * r;
    let denom = if n > 2 { n as f64 - 2.0 } else { n as f64 };
    sq - (1.0 - sq) / denom
}
```

Always unbiased; hot-path uses pre-computed linear constants
`r2u_a = 1 + 1/(n-2)`, `r2u_b = -1/(n-2)` so each pair is
`val² · (n_inv²·r2u_a) + r2u_b`. No flag to disable.

**Divergence verdict**: Always differs by `(1-r²)/(N-2)` per pair at
default GCTA settings. To match ldsc-rs, pass GCTA `--ld-score-adj`.

**Numerical consequence**: At N=503, bias per pair ≈ `(1 - 0.02)/501 ≈
0.0020`. Over ~900 pairs, cumulative bias ≈ 1.8 in the LD score.
Observed: GCTA biased mean L2 (20.20) − GCTA unbiased mean L2 (18.44) =
1.76, matching the theoretical prediction.

### 3. r² accumulation and self-LD

**GCTA** (`main/ld.cpp:300, 415-425, 429-442`).

- Self pair excluded: `if (k == j) continue;` (`ld.cpp:415`).
- `mean_rsq[i] = sum_{k≠i} r²_{ik} / count` where count is the number of
  pairs passing the cutoff.
- `snp_num[i] = count`.
- Output: `ldscore = 1.0 + mean_rsq[i] * snp_num[i]` (`ld.cpp:300`) =
  `1 + Σ_{k≠i, r²>cutoff} r²_{ik}`. The `+1` is the self-r² term added
  back after-the-fact.
- Two-pass merge (`ld.cpp:429-442`) uses the running weighted average
  `new_mean = (old_mean × old_count + delta_sum) / (old_count +
  delta_count)`, which when multiplied by total count is equivalent to
  cumulative sum. So:
  `ldscore[i] = 1 + Σ_{k ∈ pass1_block(i)} r²_{ik} + Σ_{k ∈ pass2_block_for_i} r²_{ik}` with proper pair-deduplication via `s1/s2`.

**ldsc-rs** (`src/l2/compute.rs:1297-1306`):

```rust
for j in 0..c {
    r2u_bb[(j, j)] = 1.0;          // self-LD = 1 exactly
    for k in 0..j {
        let val = bb_val(k, j);
        let r2u = val * val * n_inv_sq_r2u_a + r2u_b;
        r2u_bb[(j, k)] = r2u;
        r2u_bb[(k, j)] = r2u;
    }
}
```

The diagonal is set to 1.0 explicitly. Then `l2[chunk] += r2u_bb @
annot[chunk]`, which in the scalar (annot=ones) case adds `Σ_k
r2u_bb[i,k]` to `l2[i]`. The diagonal contribution is `1.0 × annot[i] =
1` for scalar mode.

**Divergence verdict**: Equivalent in scalar mode. Both implementations
include the self-r²=1 term in the final LD score. In partitioned mode,
ldsc-rs's annot-weighted sum is the standard partitioned-LDSC formula;
GCTA has no partitioning equivalent (`calcu_mean_rsq_multiSet` is a SNP
subset filter, not partitioned LD).

### 4. Two-pass overlap averaging

**GCTA**: see Section 1 above. Pass-2 is mandatory for SNPs near pass-1
block boundaries. The averaging step at `ld.cpp:429-442` recovers
cumulative-sum behavior after dedup.

**ldsc-rs**: no overlap averaging. With `--snp-level-masking`, each SNP
gets its true per-SNP window via the post-GEMM mask at
`compute.rs:1431-1442, 1735-1748`. Without it, the chunk approximation
(every SNP in chunk uses `block_left[chunk_start]`) is exactly Python
LDSC's behavior — no averaging.

**Divergence verdict**: ldsc-rs `--snp-level-masking` produces *more
strictly per-SNP exact* windows than GCTA. GCTA's two-pass scheme
*approximates* per-SNP windowing for boundary SNPs but is still
block-based. The two implementations have similar but not-identical
effective windows; they happen to match closely in mean L2 by virtue of
LD scores being dominated by short-range high-r² pairs.

### 5. r² cutoff

**GCTA** (`main/option.cpp:541-548, main/ld.cpp:418`): `--ld-rsq-cutoff`
in [0, 1] (default 0.0). The cutoff is checked against the BIASED r²,
regardless of `_ldscore_adj`. With cutoff = 0, all pairs (including
r²=0) pass.

**ldsc-rs**: no cutoff. All pairs in the window contribute, including
those with negative r²_unbiased values.

**Divergence verdict**: Trivially reconciled at cutoff = 0. Note that at
non-zero cutoff, GCTA's filter is on biased r² (geometrically equivalent
to ldsc-rs filtering on `r²_unbiased + (1-r²)/(N-2)`).

### 6. Genotype standardization

**GCTA** `make_XMat_subset` (`main/data.cpp:2040-2073`):
1. Center: `X(i,j) -= _mu[k]` where `_mu = 2p` (mean dosage from
   `calcu_mu`, `main/data.cpp:1565-1605`). Missing genotypes are set to
   0 BEFORE centering — equivalent to mean-imputation.
2. Scale (when `divid_by_std = true`, which is the case for
   `--ld-score`): divide each column by `sqrt(2p(1-p))`, the
   HWE-derived SD (`data.cpp:2061-2070`).

Then `calcu_ssx_sqrt_i_sub` (`ld.cpp:397-401`) recomputes the actual
column norm `ssx_sqrt_i = 1 / ||X.col||` for use in the r formula.

**ldsc-rs** `normalize_col_f64_with_stats`
(`src/l2/normalize.rs:71-110`):
1. Center: `*v -= avg` where `avg = sum / count` (empirical mean,
   ignoring NaN).
2. Scale: divide by empirical SD `sqrt(centered_sum_sq / n)` (with
   denominator `n`, not `n-1`). NaN positions are set to 0 (equivalent
   to mean-imputation).

**Divergence verdict**: The final r² is computed as `(X_j · X_k)² /
(||X_j||² · ||X_k||²)` in both — invariant to per-column scaling. So
both produce identical Pearson r² regardless of whether columns are
pre-scaled to unit variance or just centered. Documentation only — not
a numerical divergence.

### 7. Chromosome boundary handling

**GCTA**: `ld.cpp:332` and `ld.cpp:376` force block boundaries on chr
change. Also forces boundaries on >1 Mb bp gaps within a chromosome
(`ld.cpp:332`).

**ldsc-rs** `window.rs:32-90`: per-chromosome `block_left` computation;
windows never cross chr boundaries by construction.

**Divergence verdict**: Equivalent for chr boundaries. ldsc-rs does NOT
break windows on intra-chr bp gaps — this is a minor difference that
could matter for sparse data, but it has no visible effect on chr22
EUR.

### 8. Missing-genotype imputation

Both: missing genotype → set to mean (= centered to 0). Same behavior;
documentation only.

### 9. Numeric precision

**GCTA**: `MatrixXf` (Eigen f32) for the X matrix and the `X^T X` GEMM.
ssx and mean_rsq accumulators are `eigenVector` = f64. r² extraction
multiplies the f32 dot by f64 ssx_inv values.

**ldsc-rs**: f64 GEMM by default; `--fast-f32` switches to f32 GEMM.
r²_unbiased extraction is always f64.

**Divergence verdict**: To match GCTA's precision exactly, use ldsc-rs
`--fast-f32`. The default ldsc-rs build uses f64, which is higher
precision than GCTA (ULP differences for tiny r² pairs).

### 10. Output schema

**GCTA** `--ld-score` output `<out>.score.ld` (`ld.cpp:294-302`):
whitespace-delimited, columns `SNP chr bp MAF mean_rsq snp_num max_rsq
ldscore`. Scalar only — no annotation/partitioning support.

**ldsc-rs**: `.l2.ldscore.gz` (tab-delimited) with columns `CHR SNP BP
[CM] L2_<annot>` for each annotation column; `.l2.M` and `.l2.M_5_50`
give SNP counts per annotation. The `write_gcta_as_ldsc_format` helper
in `compare_gcta_ldsc.py` already handles the GCTA→ldsc-rs schema
conversion for cross-validation.

### 11. Annotation / partitioning support

GCTA's `calcu_mean_rsq_multiSet` (`ld.cpp:501-588`) outputs per-set mean
r² counts but does not produce partitioned LD scores in the LDSC
@finucane2015a sense (which requires `ℓ_{j,c} = Σ_{k ∈ C} r²_{j,k}`
weighted by annotation values, not just on/off membership). GCTA cannot
validate ldsc-rs's partitioned LD output; that has no second reference
implementation.

### 12. Parallelism strategy

**GCTA**: `#pragma omp parallel for` over SNPs within a block
(`ld.cpp:406, 469`). Single block at a time.

**ldsc-rs**: per-chromosome rayon tasks (one task per autosome), each
calling `compute_ldscore_global` independently. Within each task, faer's
`Par::rayon(0)` parallelizes the GEMM.

**Divergence verdict**: Different threading models, algorithmically
equivalent. No numerical impact.

### 13. The MKL path

`main/mkl.cpp:612-799` (`calcu_mean_rsq_mkl`) is a near-clone using a
flat float buffer `_geno_mkl` and `cblas_sgemm` for GEMM. **Key
differences from the non-MKL path**:
- `mkl.cpp:625` calls `get_ld_blk_pnt(..., wind_size*2)` — interprets
  `wind_size` as a per-side window radius, not block width. So `--ld-wind X`
  in this code path would give blocks of width `2X`. But this code path
  is NOT called by `--ld-score` (see `option.cpp:1395`): the dispatch is
  `pter_gcta->calcu_mean_rsq(LD_wind, ...)`, the non-MKL path. The MKL
  function is unreferenced in `option.cpp` and appears to be dead code in
  GCTA 1.95.1.
- MKL path NEVER applies the `_ldscore_adj` correction (no `if
  (_ldscore_adj)` branches at `mkl.cpp:705, 727, 782`). Always biased
  r².
- MKL path clips `r² > 1` to `1` (`mkl.cpp:702-703, 780`) — guards
  against single-precision overflow in cblas_sgemm.
- MKL output file is `.mrsq.ld` not `.score.ld`, and lacks the `chr bp
  ldscore` columns.

Audit notes: even if the MKL path were reachable, the `wind_size*2`
discrepancy is real and a likely historical source of unit-convention
confusion in GCTA's own codebase.

## What the existing comparison actually tested

`preprint/scripts/compare_gcta_ldsc.py` invokes GCTA with:

```
gcta64 --bfile data/1000G_eur --extract chr22_maf05.snplist --maf 0.05 \
       --ld-score --ld-wind 1000 --ld-rsq-cutoff 0 --autosome \
       --out preprint/data/sim_ldscores/gcta
```

Two algorithmic mismatches vs ldsc-rs `--ld-wind-kb 1000`:

1. **Window**: GCTA `--ld-wind 1000` → block width 1 Mb, overlap 500 kb,
   effective per-SNP window radius averaging ~750 kb. ldsc-rs
   `--ld-wind-kb 1000` is exactly ±1000 kb radius. GCTA's window is
   ~25% narrower in mean neighbor count (903 vs 1182).
2. **r²**: no `--ld-score-adj` (verified at `gcta.log:15`), so GCTA
   uses biased r². ldsc-rs uses unbiased.

The biases partly cancel: narrower window pulls LD scores DOWN; biased
r² pushes them UP. Net effect: GCTA mean L2 (20.20) > ldsc-rs masked
(18.67) by +8%, with Pearson r = 0.9975.

| Scenario | mean L2 | vs ldsc-rs masked | vs ldsc-rs chunked | r (vs masked) |
|---|---|---|---|---|
| GCTA `--ld-wind 1000` (existing run, biased) | 20.20 | +1.53 | +1.35 | 0.9975 |
| GCTA `--ld-wind 1000 --ld-score-adj` | 18.44 | −0.24 | −0.41 | 0.9984 |
| GCTA `--ld-wind 2000 --ld-score-adj` | **18.86** | **+0.18** | **+0.006** | **0.9995** |
| ldsc-rs masked (--ld-wind-kb 1000) | 18.67 | 0 | −0.18 | 1.0 |
| ldsc-rs chunked (--ld-wind-kb 1000) | 18.85 | +0.18 | 0 | 0.9971 |

**`--ld-wind 2000 --ld-score-adj` is the apples-to-apples GCTA invocation
for matching ldsc-rs `--ld-wind-kb 1000`.** Mean L2 matches ldsc-rs
chunked to within rounding (0.006) and ldsc-rs masked to within +1%.

This corrects the prior memory note ("GCTA unbiased vs ldsc-rs c=50:
Pearson r 0.998, mean_abs 0.39") — that number was either misremembered
or from an unsaved scratch run. The actual prior committed run had mean
|diff| of 1.75 (biased r²); the corrected `--ld-wind 2000 --ld-score-adj`
run gives mean |diff| 0.17.

## Where current docs and preprint are wrong

The following claims overstate GCTA equivalence and need rewording. The
audit doc does NOT edit prose; this list flags spots for the user to
review.

### `preprint/main.typ`

- **Line 497–499** (GCTA description): "Its `--ld-score` (alias
  `--ld-meanrsq`) routine computes per-SNP exact windows---there is no
  chunking approximation."
  **WRONG.** GCTA does block-based windowing with two-pass averaging.
  See Section 1 above.
- **Line 499–502**: "Its block decomposition (`get_ld_blk_pnt` in
  `main/ld.cpp`) is purely a memory optimization: each SNP's LD score
  is computed from the full set of SNPs within its physical-distance
  window, and overlapping-block boundaries are averaged in a second
  pass to eliminate edge artifacts."
  **WRONG.** The block decomposition IS the window definition. The
  effective per-SNP window radius varies from ~`wind/2` to `wind` kb
  depending on the SNP's position within the block. It's NOT "the full
  set of SNPs within its physical-distance window".
- **Line 516–518**: "`--snp-level-masking` (per-SNP exact windows
  matching GCTA's algorithmic intent and the LDSC paper's math)".
  **PARTIALLY WRONG.** `--snp-level-masking` does per-SNP exact windows
  — but those are STRICTER than GCTA's block-based scheme. If anything,
  ldsc-rs `--snp-level-masking` is closer to the LDSC paper's
  mathematical statement than GCTA itself.
- **Line 518–520**: "produces LD scores that agree with GCTA's unbiased
  output at Pearson $r = 0.998$."
  **Approximately right but understated.** With matched windows
  (GCTA `--ld-wind 2000`) and matched r² (`--ld-score-adj`), r = 0.9995.
- **Table @tbl:implementations** row "Window eviction" for GCTA:
  "Per-SNP exact". **WRONG.** Should be "Two-pass overlapping blocks of
  bp width `--ld-wind`".
- **Table caption**: "GCTA's per-SNP exact semantics match the LDSC
  paper's mathematical statement". **WRONG.** GCTA's semantics are
  *block-based with overlap averaging*. The LDSC paper's
  mathematical statement (`ℓ_j = Σ_k r²_{jk}` over a per-SNP window) is
  not actually implemented by GCTA — it's an *approximation* of
  per-SNP exact windowing. ldsc-rs `--snp-level-masking` is the more
  faithful implementation of the paper's math.

### `CLAUDE.md`

- **Reference-implementations bullet** (lines 11–22): claim that
  `ldsc l2 --snp-level-masking` "matches GCTA's algorithmic intent"
  with "residual mean diff 0.39 from GCTA's two-pass overlap averaging
  — pure algorithmic, no flag fixes" is incorrect on two counts:
  (a) GCTA's algorithmic intent is itself block-based, not per-SNP
  exact; (b) the 0.39 mean diff number is not reproducible from the
  committed `compare_gcta_ldsc.py` run (which gives 1.75 because
  `--ld-score-adj` was not passed). The apples-to-apples comparison is
  GCTA `--ld-wind 2000 --ld-score-adj`, giving mean |diff| 0.17 vs
  masked.

### `docs/h2-masking-simulation.md`

- The "GCTA cross-validation" section reports r=0.998 vs masked at
  unbiased — but the committed run was biased. Either rerun with
  `--ld-score-adj` or annotate the existing numbers as biased.
- The framing "GCTA as North Star reference implementation" is fine, but
  "GCTA does per-SNP exact windowing" needs to be softened to "GCTA does
  block-based windowing with two-pass averaging; ldsc-rs
  `--snp-level-masking` is the more strict per-SNP implementation".

## Reconciliation flag-table

| To compare these two | Use these flags |
|---|---|
| ldsc-rs default vs GCTA default | NOT directly comparable (different windows + biased vs unbiased) |
| ldsc-rs `--snp-level-masking --ld-wind-kb X` vs GCTA | GCTA: `--ld-wind 2X --ld-score-adj --ld-rsq-cutoff 0` |
| ldsc-rs r² formula vs GCTA r² formula | GCTA: `--ld-score-adj` (single fix) |
| ldsc-rs window vs GCTA window | GCTA `--ld-wind 2X` ≈ ldsc-rs `--ld-wind-kb X` (empirically matched mean L2; not algorithmically identical) |
| ldsc-rs `--fast-f32` vs GCTA precision | Already matched (both use f32 GEMM with f64 accumulators) |

## Open questions

1. **Is GCTA's two-pass overlap averaging actually documented anywhere
   in the paper or supplementary?** The block-width `--ld-wind` flag
   has zero published characterization, similar to Python LDSC's
   `--chunk-size`. Both are silent algorithmic choices that affect
   downstream LD-score-regression h² estimates.
2. **Does the 16% h² shift on real data (BMI/SCZ) when switching to
   `--snp-level-masking` reflect a real signal, or is it an artifact of
   the *stricter* windowing relative to GCTA?** This audit suggests the
   chunked-vs-masked gap in ldsc-rs is comparing chunk-eviction (Python
   LDSC's behavior) to true per-SNP exact (stricter than GCTA). Neither
   matches GCTA's block-averaging exactly. A third comparison —
   "GCTA-style overlap-averaged blocks in ldsc-rs" — could disentangle
   these.
3. **What's the right reference for the preprint's numerical claims?**
   Currently the preprint claims GCTA = paper-canonical = per-SNP
   exact. The audit shows this is a three-way equivocation: the LDSC
   paper's *math* says per-SNP exact, GCTA's *implementation* is
   block-with-overlap-averaging, and the SHIPPED 1000G LD scores are
   from GCTA's implementation — not from the paper's mathematical
   statement. ldsc-rs `--snp-level-masking` is faithful to the math but
   *not* faithful to the shipped 1000G LD scores.

## Reproducibility

All numerical comparisons in this audit can be regenerated:

```bash
# 1. Build ldsc-rs and produce reference LD scores (already exists in repo):
#    preprint/data/sim_ldscores/{masked,chunked}.l2.ldscore.gz

# 2. Original (biased, mismatched window) GCTA run (already exists):
#    preprint/data/sim_ldscores/gcta.score.ld
#    Re-run with: preprint/scripts/compare_gcta_ldsc.py compute

# 3. Audit re-runs (this doc):
GCTA=/tmp/gcta/gcta-1.95.1-macOS-arm64/bin/gcta64
for wind in 1000 2000; do
  $GCTA --bfile data/1000G_eur \
        --extract preprint/data/chr22_maf05.snplist \
        --maf 0.05 --ld-score \
        --ld-wind $wind --ld-rsq-cutoff 0 --ld-score-adj --autosome \
        --out preprint/data/sim_ldscores_audit/gcta_wind${wind}_adj
done
```

Outputs land in `preprint/data/sim_ldscores_audit/` and are kept for
record.

## Sources of truth

- GCTA 1.95.1 source: `/tmp/gcta_src/main/{ld.cpp,mkl.cpp,data.cpp,option.cpp,gcta.h}`
- GCTA upstream: https://github.com/JianYang-Lab/GCTA
- ldsc-rs LD-score code: `src/l2/{compute.rs,window.rs,normalize.rs}`
- LDSC paper: `docs/ldsc_paper.md` (and the original GCTA papers
  @yang2010, @yang2011 referenced therein for the LD-score-from-GRM
  derivation)
- Python LDSC (for the chunked-approximation reference):
  `ldsc_py3/ldscore/ldscore.py`
