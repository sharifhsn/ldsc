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
**Status: PARTIAL** — full §"Three implementations: GCTA, Python LDSC,
ldsc-rs" added to `preprint/main.typ` with the algorithmic-choice table
@tbl:implementations and supporting cross-validation in
[docs/h2-masking-simulation.md](../docs/h2-masking-simulation.md). GCTA
positioned as the canonical reference (it's the tool the original LDSC
paper used to compute the 1000G LD scores). BOLT-LMM, LDpred, SumHer
positioning still missing — those are h² estimators rather than LD-score
generators, so a single positioning paragraph in the Introduction (rather
than a numerical comparison) would close this fully.

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
**Status: DONE** (for the algorithmic part; real-data mechanism still open).

The masking section in `preprint/main.typ` now explicitly states that the
algorithmic correction accounts for only ~1% of the 15-16% real-data effect,
with the remainder open for future work. The added §"chunk-size knob and
implementation parity" presents the sensitivity sweep (h²_true=0.5: ranges
+10.5% to +11.6% across c=25..1000 — full chunk-size effect ~1.1pp of
relative bias). The added §"Three implementations" frames `--snp-level-masking`
as the paper-canonical mode matching GCTA, not as a correction to a bug.

Simulation: controlled study with known true h² shows both estimators are
biased *upward*; masking adds ~1% to that bias (paired Δ at h²=0.5: +1.2% of
chunked, 95% CI [+0.43%, +0.79%]). All four LD-score generators tested
(Python LDSC, ldsc-rs chunked, ldsc-rs masked, GCTA) produce h² estimates
within 1% of each other — bias is in the regression machinery at small N,
not in any LD-score generator. See
[docs/h2-masking-simulation.md](../docs/h2-masking-simulation.md) for the
4-way triangulation and `preprint/data/h2_simulation_chunk_sweep.csv` for
the chunk-size sweep raw data.

Still TODO (genuinely open): localize what *does* explain the real-data
+15% shift on BMI/SCZ. Candidates: population-structure interactions,
finite-sample biases compounding with M=1.66M (vs simulation's M=18.6K),
sparse polygenicity, or LDSC-model misspecification. Listed as follow-ups
in the doc.

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
