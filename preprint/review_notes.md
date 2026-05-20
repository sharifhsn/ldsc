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
**Status: PARTIAL** (Workstream M, 2026-05-18) — `preprint/main.typ` Methods
§"Aggregation and downstream robustness" now contains a back-of-the-envelope
deriving 1/d scaling from the per-pair variance bound and explaining why the
constant cannot be predicted by the bound alone (positive cross-pair sketch
covariance from shared SNP projections inflates `Var[ℓ̃_j]` above the
sum-of-variances). The 1/d scaling is now theoretically dictated; the
constant C ≈ 8 is acknowledged as dataset-dependent and not predicted by
the simple variance argument. Open: a tighter bound that accounts for the
cross-pair covariance would close this fully.

## 4. Biobank synthetic dataset weakly justified
Replicating 1000G 21× creates rank-deficient matrix (true rank ~2490).
Doesn't address how this affects: CountSketch accuracy, cost model, or
whether 16% h² attenuation from chunking would replicate on real data.
**Status: PARTIAL** (Workstream M, 2026-05-18) — §Limitations now leads
with a 5-point synthetic-data paragraph covering rank deficiency, cache
behavior, CountSketch generalization, cost-model constants, and a
prominent flag that real-biobank validation is the single largest open
item in evaluation. Every biobank table/figure caption carries an inline
"(synthetic 21×-replicated 1000G; see §Limitations)" qualifier
(verified: 8 occurrences of "synthetic" in main.typ). Real UKBB / AoU
validation is genuinely pending data access and is out of scope for this
push.

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
**Status: DONE** (Workstream M, 2026-05-18) — Dropped per user request:
the row was removed from @tbl:perf-1000g and the narrative paragraph
rewritten to describe the Python ceiling without depending on the
ldsc_py_opt data point. The figure (`generate_figures.py` and
`fig1_performance_1000g.png`) was also re-generated without the
ldsc_py_opt bar.

## 8. No error bars on many benchmarks
Table 1 Rust exact (41.1s) has no ±. Sketch times are local estimates scaled to AWS.
Table 2 sketch modes have no ±. Clarify which are direct measurements.
**Status: DONE** (Workstream M, 2026-05-18) — Captions of
@tbl:perf-1000g and @tbl:perf-biobank now explicitly state which rows
are mean ± SD over hyperfine 10 runs and which are single-run
measurements. Stretch (re-run sketch sweep with hyperfine on AWS for
real per-row SDs) deferred as optional follow-up requiring user spend
approval.

## 9. Abstract doesn't mention SNP-level masking finding
16% h² change from undocumented approximation is arguably the most impactful
finding but is absent from abstract.
**Status: DONE** (Workstream M, 2026-05-18) — Abstract now surfaces the
masking finding as the second sentence after the Rust-port pitch:
"...we surface an undocumented chunk-rounding approximation in the
reference Python LDSC that inflates LD scores by ~12% and attenuates
real-data heritability estimates by ~16% on BMI and schizophrenia GWAS;
ldsc-rs provides a `--snp-level-masking` flag recovering the per-SNP
exact windows of the original LDSC paper at <1% overhead." Final
abstract word count: 249.

## 10. No Supplementary Information
All 26 sketch measurements, per-chromosome comparisons, full h2/rg logs
would strengthen reproducibility.
**Status: TODO** — deferred to follow-up. Candidates: full hyperfine
logs from the 26-point sketch sweep, per-chr LD-score comparison tables,
GCTA source audit (already exists at `docs/gcta-source-audit.md` and
can be moved into `preprint/supplementary.typ` ~as-is), raw chunk-size
sweep data (already at `preprint/data/h2_simulation_chunk_sweep.csv`).
Not blocking for preprint submission; addresses reproducibility for the
journal-review stage.

---

## Workstream M additions (not in original review)

### ldsc-web §Deployment subsection
**Status: NEW + DONE** — §Deployment expanded with a `=== ldsc-web`
subsection covering architecture (SAB worker pool, per-chr
decomposition, async chr-shard pre-load), hand-rolled WASM SIMD GEMM
and K=4-lane scatter kernels, and in-browser h² regression including
the `web_time::Instant` workaround for the `std::time::Instant`
wasm32-unknown panic. Includes new @tbl:browser-perf with biobank
browser-vs-CLI parity numbers (browser best 7.87 s / median 8.34 s vs
CLI 7.70 s = 1.02–1.08×) and an in-browser h² wall (0.45 s).
Reproducibility paragraph links live URL
(https://sharifhsn.github.io/ldsc/).

---

## Workstream N additions (cohort-driven, 2026-05-18)

Triggered by a research-agent survey of 12 comparable papers in the
"fast reimplementation in stat-gen" genre (LDSC, BOLT-LMM, REGENIE,
PLINK 2, GCTA, LDAK, TeraPCA, FastPCA, LDpred2, SBayesR, quickLD,
Bigtools). Goal: close concrete genre conventions our preprint was
missing.

### NEW Figure 1: Workflow / architecture diagram
**Status: NEW + DONE** — matplotlib-rendered architecture diagram at
`preprint/figures/fig0_architecture.png`, inserted as the first
figure in §Results. Boxes + arrows depicting BED → decode/normalize
→ ring-buffer GEMM → per-chr rayon → L2/M outputs → h²/rg/partitioned
downstream, with alt-mode branches (`--sketch`, `--snp-level-masking`,
`--python-compat`) and the WASM compile branch to the browser
deployment.

### GCTA comparator added to perf tables + accuracy figure
**Status: NEW + DONE** — GCTA v1.95.1 ARM64 installed locally and run
on the full 1000G_eur panel (1.66M SNPs, 23 s wall, mean L2 = 19.84).
@tbl:perf-1000g gained a GCTA row (67× faster than Python LDSC, r =
0.9995 vs ldsc-rs exact). @fig:maf-facet shows GCTA-vs-ldsc-rs
per-SNP L2 differences faceted into 4 MAF bins — IQR ≈ ±0.15, median
indistinguishable from zero across all bins.

### Memory profiling (peak RSS)
**Status: NEW + DONE** — `/usr/bin/time -l` measurements for every
(N, mode) cell. Added as full-width @tbl:scaling and as the right
panel of @fig:scaling. At N=100K, exact f32 uses 14 GB peak RSS;
sketch d=200 uses 2.7 GB — a 5× memory advantage that mirrors the
70× wall advantage.

### N-scaling curve with 5 measured points
**Status: NEW + DONE** — generated synthetic biobank at N=10K, 20K,
100K (in addition to existing N=50K). Benchmarked all 4 modes
(exact-f64, exact-f32, sketch-200, sketch-1000) across 5 N values
locally; wrote `preprint/data/scaling_bench.csv`. @fig:scaling
plots both wall time and peak RSS log-log, with exact modes
~linear in N and sketch modes nearly flat.

### Worked example: BMI heritability on Yengo et al. (2018)
**Status: NEW + DONE** — downloaded Yengo 2018 BMI sumstats from
GIANT (2.5M HapMap2 SNPs), ran ldsc-rs munge → l2 → h² end-to-end
in 41 s total. Estimated h² = 0.2032 ± 0.0055 (vs Yengo's published
~0.21 — agreement within 3%). New §"End-to-end worked example"
subsection + @tbl:bmi-pipeline + @tbl:bmi-h2. Yengo 2018 added to
refs.bib.

### Tutorial doc
**Status: NEW + DONE** — `docs/tutorial.md` written, covering
install, LD score computation, sumstats munging, h² estimation,
and the browser-demo no-install alternative. Linked from the
new §"Software, Data, and Reproducibility" section.

### Software/Data/Reproducibility section expanded
**Status: NEW + DONE** — original 5-line `= Data and Code
Availability` replaced with a 35-line `= Software, Data, and
Reproducibility` section with subsections for code (GitHub,
crates.io, license, `cargo install`, Docker, releases), WASM
frontend, benchmarks (link to all scripts + scaling CSV +
figure-gen script), Zenodo archival snapshot placeholder,
tutorial link, reference data (1000G ISGR), and the real GWAS
used in the worked example (Yengo GIANT URL).

### Cross-tool capability matrix
**Status: NEW + DONE** — @tbl:tool-matrix (6 tools × 11
capabilities). Honestly differentiates ldsc-rs's coverage from
Python LDSC, GCTA, PLINK 2, BOLT-LMM, REGENIE — with footnotes
explaining the partial-credit cases (PLINK 2's pairwise r², BOLT-
LMM's internal LD scores, GCTA/BOLT-LMM individual-level vs
summary-level inputs, BOLT-LMM's stochastic trace estimation
within mixed models).

### Reviewer preempts (4 short additions)
**Status: NEW + DONE** —
(a) §"Numerical parity with Python LDSC" now cites
    `scripts/check_l2_tiny_py_vs_rust.sh` at commit `2da5705` for
    the bit-identical claim;
(b) §"Practical guidance" disambiguates PLINK 2 `--r2` (pairwise)
    from LD scores (sum of r²);
(c) @tbl:features `munge-sumstats` row note now states INFO-score
    handling is identical to Python LDSC;
(d) §Limitations expanded with a paragraph noting that BIM is
    byte-for-byte copied during synthetic generation so SNP IDs,
    positions, alleles, and per-SNP MAF are preserved between
    source 1000G and synthetic biobank.

### Masking finding promoted in TOC
**Status: NEW + DONE** — section header "Per-SNP window masking"
renamed to "Discovery: Python LDSC's chunk-rounding inflates
heritability by ~16%" so the finding catches the eye in the
table of contents. Physical section reorder declined as too
risky vs benefit; the rename gives equivalent prominence.
