# Adversarial Review Notes (2026-03-26)

## 1. Unbiased r² formula is wrong/inconsistent (line 589)
Paper writes r̂²_unbiased = (Nr̂² - 1)/(N-1). Code uses r² - (1-r²)/(N-2).
Need to either show algebraic equivalence or use the form matching the code,
with citation to specific equation in Bulik-Sullivan supplementary note.
**Status: DONE** — Fixed formula to match code (adjusted R², denom N-2), added
derivation from supp note Eq 1.7, noted O(1/N²) difference.

## 2. Introduction cites no related work besides NCI LDscore
No mention of BOLT-LMM, GCTA-LDMS, LDpred, SumHer, PLINK --ld-window, etc.
Need paragraph positioning ldsc-rs relative to other LD computation tools.
**Status: TODO**

## 3. CountSketch accuracy model r(d) ≈ 1 - C/d is empirical only
No derivation connecting CountSketch variance bounds to Pearson r of LD score
vector. Where does constant 8 come from? Why does accuracy saturate above d≈1000?
**Status: TODO**

## 4. Biobank synthetic dataset weakly justified
Replicating 1000G 21× creates rank-deficient matrix (true rank ~2490).
Doesn't address how this affects: CountSketch accuracy, cost model, or
whether 16% h² attenuation from chunking would replicate on real data.
**Status: TODO**

## 5. Per-SNP masking section doesn't address whether masked or chunked is more accurate
16% h² increase — but is the masked estimate less biased, or just different?
Does LDSC regression theory assume exact or approximate windows?
Thousands of published h² estimates used inflated LD scores.
**Status: TODO**

## 6. Methods section on LD score computation is too thin
Missing: GEMM decomposition, ring-buffer pseudocode, normalization procedure,
how unbiased estimator is applied to GEMM output.
**Status: DONE** — Added chunked GEMM decomposition (B^TB within-chunk, A^TB
cross-chunk), genotype normalization description, window boundary algorithm
with block_left rounding, and SNP-level masking mechanism.

## 7. ldsc_py_opt comparison is a single data point with no breakdown
No description of what optimizations were applied or where the ceiling comes from.
Should show why remaining gap is structural.
**Status: TODO**

## 8. No error bars on many benchmarks
Table 1 Rust exact (41.1s) has no ±. Sketch times are local estimates scaled to AWS.
Table 2 sketch modes have no ±. Clarify which are direct measurements.
**Status: TODO**

## 9. Abstract doesn't mention SNP-level masking finding
16% h² change from undocumented approximation is arguably the most impactful
finding but is absent from abstract.
**Status: TODO**

## 10. No Supplementary Information
All 26 sketch measurements, per-chromosome comparisons, full h2/rg logs
would strengthen reproducibility.
**Status: TODO**
