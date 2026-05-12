# CountSketch in ldsc-rs: mathematical analysis

**Date**: 2026-05-11. **Scope**: rigorous derivation and Monte Carlo
validation of the CountSketch estimator used in `--sketch d`, including
the renormalization step, the linear bias correction, and the LD-score-level
behavior. Identifies three concrete defects in the current implementation
and outlines four alternative optimization paths.

The implementation under analysis is in `src/l2/compute.rs`:
- `CountSketchProj::project_f64` / `project_f32` (lines 178-302)
- Fused BED-decode + project kernels (lines 304-725)
- Bias-correction constants `active_n_inv_sq_r2u_a` / `active_r2u_b` (lines 799-818)
- Hot-loop r² extraction (lines 1297-1306 for B×B, 1599-1611 / 1716-1727 for A×B)

All claims here are validated by Monte Carlo simulation (N=2000, 30 000 trials
per data point, full reproduction script at end of doc).

---

## 1. Problem setup

Let $X \in \mathbb{R}^{N \times M}$ be the genotype matrix with columns
standardized to mean 0 and $\|X_j\|^2 = N$ (so empirical variance 1).

The quantities of interest:
- empirical correlation $\hat r_{jk} := X_j^\top X_k / N$
- empirical squared correlation $\hat r_{jk}^2$
- LD score $\ell_j := \sum_{k \in W(j)} \hat r_{jk}^2$ over a window

CountSketch is parameterized by a hash $h:[N] \to [d]$ and sign
$s:[N] \to \{\pm 1\}$, drawn uniformly and independently. The sketch
matrix $S \in \mathbb{R}^{d \times N}$ has entries $S_{a,i} = s_i \cdot \mathbb{1}[h(i) = a]$.

Sketched columns: $\tilde X_j = S X_j \in \mathbb{R}^d$, i.e.,
$\tilde X_{j,a} = \sum_{i: h(i)=a} s_i X_{ij}$.

In code: `CountSketchProj::project_f64` performs exactly this scatter-add.

## 2. Base estimator: $\mathbb{E}[\tilde X_j^\top \tilde X_k] = X_j^\top X_k$ exactly

### Lemma 1 (unbiased dot product)

Expand:
$$\tilde X_j^\top \tilde X_k = \sum_{a=1}^d \tilde X_{j,a} \tilde X_{k,a}
= \sum_{i,i'} s_i s_{i'} \mathbb{1}[h(i)=h(i')] X_{ij} X_{i'k}.$$

For $i \neq i'$, $\mathbb{E}[s_i s_{i'}] = 0$ by independence and zero-mean.
For $i = i'$, $s_i^2 = 1$ and $\mathbb{1}[h(i)=h(i)] = 1$, so the diagonal
contributes $\sum_i X_{ij} X_{ik} = X_j^\top X_k$.

$\mathbb{E}[\tilde X_j^\top \tilde X_k] = X_j^\top X_k. \quad\blacksquare$

This is exact, not asymptotic — no Taylor expansion or large-$N$ limit
needed.

### Lemma 2 (variance)

The fluctuation is the off-diagonal piece:
$$\tilde X_j^\top \tilde X_k - X_j^\top X_k = \sum_{i \neq i'} s_i s_{i'}
\mathbb{1}[h(i)=h(i')] X_{ij} X_{i'k}.$$

Taking expectation of the square, sign-pairing kills all terms except
those where $\{i, i'\}$ and $\{l, l'\}$ pair up exactly. Two cases survive:
- $(l, l') = (i, i')$: contributes $\frac{1}{d} X_{ij}^2 X_{i'k}^2$
- $(l, l') = (i', i)$: contributes $\frac{1}{d} X_{ij} X_{ik} X_{i'j} X_{i'k}$

Summing over $i \neq i'$ and applying
$\sum_{i \neq i'} a_i b_{i'} = (\sum a_i)(\sum b_{i'}) - \sum a_i b_i$:

$$\mathrm{Var}[\tilde X_j^\top \tilde X_k] = \frac{1}{d}\left(
\|X_j\|^2 \|X_k\|^2 + (X_j^\top X_k)^2 - 2 K_{jk}\right)$$

where $K_{jk} := \sum_i X_{ij}^2 X_{ik}^2$ is the **cross-fourth-moment**.

For standardized columns ($\|X_j\|^2 = N$, $X_j^\top X_k = N r_{jk}$):

$$\mathrm{Var}[\tilde r_{jk}^{(0)}] = \frac{1 + r_{jk}^2 - 2K_{jk}/N^2}{d}
\approx \frac{1 + r_{jk}^2}{d}$$

where $\tilde r_{jk}^{(0)} := \tilde X_j^\top \tilde X_k / N$ is the
"raw" sketched correlation (no renormalization).

### Data dependence of the variance

The leading-order $(1+r^2)/d$ bound is only correct when $K_{jk}/N^2 \ll 1$.
For real genotype data:

- Standardized binomial(2, $p$) at MAF $p = 0.5$: $K_{jj}/N = 1$ exactly.
  $K_{jk}/N^2 = O(1/N)$, negligible at biobank scale.
- At MAF $p$ small: standardized values are
  $\{-\sqrt{2p/(1-p)}, (1-2p)/\sqrt{2p(1-p)}, \sqrt{2(1-p)/p}\}$
  with probabilities $\{(1-p)^2, 2p(1-p), p^2\}$. The fourth moment
  $\mathbb{E}[X^4] = (1-2p)^2 \cdot (2p(1-p)/(2p(1-p)))^2 + p^2 \cdot 4(1-p)^2/p^2 + (1-p)^2 \cdot 4p^2/(1-p)^2 = O(1/p)$ at small $p$.

  So at MAF=0.01, $K_{jj}/N \approx 50$, ten times larger than at MAF=0.5.

**Implication**: the simple $(1+r^2)/d$ variance bound is loose for rare
SNPs. The actual variance of sketched LD scores at low MAF is larger
than the asymptotic prediction. This is not currently surfaced in the
ldsc-rs docs or preprint.

## 3. The squared estimator's bias

Squaring an unbiased estimator introduces bias equal to its variance:
$$\mathbb{E}[(\tilde r_{jk}^{(0)})^2] = \mathrm{Var}[\tilde r_{jk}^{(0)}]
+ (\mathbb{E}[\tilde r_{jk}^{(0)}])^2
= r_{jk}^2 + \frac{1+r_{jk}^2}{d} + O\!\left(\frac{1}{Nd}\right).$$

**Substituting $r^2 = 0$ (unlinked SNPs)**: bias is $1/d$. At $d = 50$,
that's 0.02 per pair. Over a window of 1000 pairs, this adds ~20 to the
LD score from pure noise — comparable to the LD score itself.

**Bias correction is essential, not optional, at low $d$.**

## 4. The renormalized cosine estimator (what ldsc-rs actually computes)

After projection, each sketched column is renormalized in
`CountSketchProj::project_f64` (lines 240-249):

```rust
let scale = (n / nrm_sq).sqrt();
for v in dst.iter_mut() {
    *v *= scale;
}
```

so that $\|\tilde X_j'\|^2 = N$ exactly after renorm. The squared cosine
$$\tilde r'^2_{jk} = \frac{(\tilde X_j^\top \tilde X_k)^2}{\|\tilde X_j\|^2 \|\tilde X_k\|^2}$$

is the actual estimator that flows into the hot-loop r² extraction.

### Bias of $\tilde r'^2$: Taylor analysis

Let $A = \tilde X_j^\top \tilde X_k$, $B = \|\tilde X_j\|^2$,
$C = \|\tilde X_k\|^2$. Means: $\bar A = Nr$, $\bar B = \bar C = N$.

Variances and covariances (leading order in $1/N$, $1/d$):
- $\mathrm{Var}(A) \approx N^2(1+r^2)/d$
- $\mathrm{Var}(B) = \mathrm{Var}(C) \approx 2N^2/d$ (Lemma 2 with $j = k$, so $r_{jj}=1$)
- $\mathrm{Cov}(A, B) \approx 2N^2 r/d$
- $\mathrm{Cov}(A, C) \approx 2N^2 r/d$ (by symmetry)
- $\mathrm{Cov}(B, C) \approx 2N^2 r^2/d$

(Derivations for the covariances follow the same sign-pairing argument as
Lemma 2; subleading kurtosis terms dropped.)

Taylor-expand $f(A,B,C) = A^2/(BC)$ around $(Nr, N, N)$. The first-order
terms vanish in expectation. Partial derivatives at the mean:
- $f_{AA} = 2/(BC) \to 2/N^2$
- $f_{BB} = 2A^2/(B^3C) \to 2r^2/N^2$
- $f_{CC} = 2r^2/N^2$
- $f_{AB} = -2A/(B^2C) \to -2r/N^2$
- $f_{AC} = -2r/N^2$
- $f_{BC} = A^2/(B^2C^2) \to r^2/N^2$

Combining:

$$\mathbb{E}[\tilde r'^2] - r^2 \approx \frac{1}{d}\Big[
\underbrace{(1+r^2)}_{\mathrm{Var}(A)} + \underbrace{4r^2}_{\mathrm{Var}(B)+\mathrm{Var}(C)}
- \underbrace{8r^2}_{\mathrm{Cov}(A,B)+\mathrm{Cov}(A,C)}
+ \underbrace{2r^4}_{\mathrm{Cov}(B,C)}\Big]
= \frac{1 - 3r^2 + 2r^4}{d}$$

Factoring:
$$\boxed{\mathbb{E}[\tilde r'^2] - r^2 \approx \frac{(1 - r^2)(1 - 2r^2)}{d}}$$

The bias is exactly zero at $r^2 \in \{0.5, 1\}$, positive for $r^2 < 0.5$,
and negative for $r^2 > 0.5$.

### Empirical validation

Monte Carlo (N=2000, 30 000 trials, varying $r$):

| $r^2$ | $d$ | empirical bias | predicted $(1-r^2)(1-2r^2)/d$ |
|---|---|---|---|
| 0.000 | 200 | +0.00504 | +0.00500 |
| 0.273 | 200 | +0.00369 | +0.00394 |
| 0.501 | 200 | +0.00213 | +0.00186 |
| 0.693 | 200 | −0.00044 | +0.00010 |
| 0.892 | 200 | −0.00042 | −0.00060 |
| 0.000 | 50 | +0.02009 | +0.02000 |
| 0.501 | 50 | +0.00840 | +0.00745 |
| 0.892 | 50 | −0.00216 | −0.00242 |

Agreement is excellent. The empirical bias matches the analytical formula
to within Monte Carlo noise at all $r^2$ and $d$ tested.

## 5. The linear correction in `compute.rs`

The code applies:
$$f(\tilde r^2) = \tilde r^2 \cdot \frac{d}{d - 2} - \frac{1}{d - 2}$$

implemented in the hot loop as
```rust
*r2u_val = *ab_val * *ab_val * active_n_inv_sq_r2u_a + active_r2u_b;
```
where `active_n_inv_sq_r2u_a = n_inv² · r2u_a · d/(d-2)` and
`active_r2u_b = r2u_b - r2u_a/(d-2)` (lines 799-818). This folds the
sketch correction into the same linear-in-$\tilde r^2$ form as the
LDSC unbiased r² correction.

The base unbiased correction $\hat r^2 - (1 - \hat r^2)/(N-2)$ ldsc-rs
applies is from Bulik-Sullivan 2015, citing **Yin & Fan 2001, "Estimating
$R^2$ Shrinkage in Multiple Regression," J Experimental Education
69:203-224** — a general-statistics regression-shrinkage estimator
(Olkin-Pratt-style), not an LD-population-genetics result. This is the
only such citation in the LDSC math; it is mathematically reasonable for
diploid biallelic data treated as Pearson correlation. A more rigorous
LD-population-genetics alternative is **Ragsdale & Gravel 2020,
"Unbiased Estimation of Linkage Disequilibrium from Unphased Data" (Mol
Biol Evol 37:923, DOI 10.1093/molbev/msz265)**, which accounts for
Hardy-Weinberg / unphased-genotype sampling structure and notes that
"widely used estimators for r² and D² exhibit large and variable upward
biases." We do not adopt Ragsdale-Gravel here because (a) consistency
with reference LDSC outputs is the higher priority for replication
studies, and (b) the difference is small in our regime; flagged for
future consideration.

The correction is exact iff the bias model is $(1 - 2r^2)/d$, but the
true bias is $(1 - r^2)(1 - 2r^2)/d$. So the corrected estimator has
residual bias:

$$\mathbb{E}[f(\tilde r^2)] - r^2 = \frac{r^2(2r^2 - 1)}{d - 2}$$

**Properties**:
- Exact at $r^2 \in \{0, 1/2\}$
- Worst case at $r^2 = 1$: residual $= 1/(d-2) \approx 0.005$ at $d = 200$
- Negative residual for $r^2 \in (0, 1/2)$, positive for $r^2 \in (1/2, 1]$

### Empirical validation at $d = 200$

| $r^2$ | predicted residual | empirical residual |
|---|---|---|
| 0.000 | 0.000000 | −0.00001 ✓ |
| 0.092 | −0.000380 | +0.00015 (MC noise) |
| 0.256 | −0.000631 | −0.00073 ✓ |
| 0.497 | −0.000014 | −0.00003 ✓ |
| 0.808 | +0.002517 | +0.00222 ✓ |
| 0.980 | +0.004759 | +0.00478 ✓ |

Bottom line: **the linear correction is exact at $r^2 = 0$** (where most LD
pairs sit) and slightly biased for the high-$r^2$ pairs that contribute
most to LD score. The per-pair residual is at most ~0.005 at $d = 200$.

### Exact (quadratic) correction

To remove the residual bias, invert $\tilde r^2 = r^2 + (1-r^2)(1-2r^2)/d$
as a quadratic in $r^2$:
$$2r^4 + (d - 3) r^2 - (d \tilde r^2 - 1) = 0$$
$$r^2 = \frac{-(d-3) + \sqrt{(d-3)^2 + 8(d \tilde r^2 - 1)}}{4}$$

(taking the + root, valid for $\tilde r^2 \geq 0$). Sanity checks:
- $\tilde r^2 = 0$: $r^2 = ({-(d-3) + \sqrt{(d-3)^2 - 8}})/{4} \approx 0$ for large $d$
- $\tilde r^2 = 1$: $r^2 = ({-(d-3) + \sqrt{(d-3)^2 + 8(d-1)}})/{4} = ({-(d-3) + (d+1)})/{4} = 1$ ✓

**Implementation cost**: 1 sqrt + ~5 fp ops per pair. Modern AVX2 `vsqrtps`
does 8 sqrts in parallel at ~10-15 cycles latency, so vectorized cost is
~2 ns per pair. Compared to the current 2-mul/1-add hot loop (~0.5 ns
per pair), this is ~4× slower per pair, but the per-pair r² extraction
is a small fraction of total wall-clock (GEMM dominates). Estimated
total impact: 5-10% slower.

**Whether to ship it**: the residual bias is ~0.005 per pair, ~0.3% at the
LD-score level (per simulation in §7). For LD-score regression that's
well below the noise floor. **Probably not worth the complexity** unless
a user needs strict mathematical unbiasedness for individual r² values.

## 6. Variance reduction from renormalization

The current code comment claims "reduces per-pair variance by ~2× compared
to the raw estimator." This is **incorrect and misleading**. The actual
variance ratio $\mathrm{Var}(\tilde r_0^2) / \mathrm{Var}(\tilde r'^2)$
depends sharply on $r^2$:

| $r^2$ | variance ratio (raw / renormalized) at $d = 200$ |
|---|---|
| 0.000 | 1.03 |
| 0.092 | 1.43 |
| 0.256 | 2.39 |
| 0.497 | 6.16 |
| 0.808 | 49.3 |
| 0.980 | 5233 |

At low $r^2$ the renorm does almost nothing. At high $r^2$ it produces
enormous reductions (because $r'^2 \in [0, 1]$ is bounded while $r_0^2$
can swing wildly with the variance of $\|\tilde X_j\|^2$).

**The real reason renorm matters**: it provides heavy-tail protection.
A few high-LD pairs in a window can have raw squared sketches with
variance hundreds to thousands of times larger than the typical pair.
Renorm bounds the contribution.

At the LD-score level (200 neighbors, exponential LD decay), the variance
ratio is ~14×.

## 7. LD-score-level behavior

Simulation: 1 target SNP + 200 neighbors with $r_k = \exp(-k/30)$
(realistic LD decay), N=2000, 500 sketches per $d$ value.

Population LD score $\sum_k r_k^2 = 14.51$. Empirical LD score = 14.58.

| $d$ | mean corrected $\hat \ell$ | bias | std($\hat \ell$) | CV |
|---|---|---|---|---|
| 50 | 14.53 | −0.04 (0.3%) | 1.58 | 10.9% |
| 100 | 14.58 | +0.01 (0.05%) | 1.17 | 8.0% |
| 200 | 14.53 | −0.05 (0.3%) | 0.79 | 5.5% |
| 500 | 14.59 | +0.01 (0.05%) | 0.50 | 3.4% |

**Findings**:
- Bias is ~0.3% even at $d=50$, decaying to 0.05% at $d=500$.
- Variance scales as $1/d$ (CV halves when $d$ quadruples).
- At $d = 500$, the CV is 3.4% — well below the noise floor of downstream
  h² regression.

Comparison with the uncorrected ($\tilde r'^2$) and raw ($\tilde r_0^2$)
estimators:

| $d$ | bias corrected | bias uncorrected (no $d/(d-2)$ fix) | bias raw (no renorm) |
|---|---|---|---|
| 50 | −0.04 | +3.38 | +4.05 |
| 200 | −0.05 | +0.81 | +0.88 |
| 500 | +0.01 | +0.35 | +0.45 |

The linear correction is doing real work — without it, LD scores would
be biased by ~5% at $d = 50$ and ~3% even at $d = 500$.

## 8. Three concrete defects in the current implementation

### Defect 1: bias formula comment is wrong

**Location**: `src/l2/compute.rs:806`.

**Current**: "The squared cosine r̃² = (val̃'/N)² has bias ≈ (2r⁴ - 7r² + 1)/d"

**Correct**: "(2r⁴ - 3r² + 1)/d" (factors as $(1-r^2)(1-2r^2)/d$)

The $-7r^2$ should be $-3r^2$. The correction code itself doesn't use
this formula directly (it uses the linear correction), so this is
comment-only. But it would mislead any future implementer trying to
verify or extend the math.

### Defect 2: variance-reduction claim is misleading

**Location**: `src/l2/compute.rs:805`.

**Current**: "This reduces per-pair variance by ~2× compared to the raw estimator."

**Correct**: The reduction ratio is sharply $r^2$-dependent, ranging
from ~1× at $r^2 = 0$ to >5000× at $r^2 = 0.98$. The "2×" figure is
correct only at $r^2 \approx 0.25$. At the LD-score level the effective
ratio is ~14× because the high-LD pairs dominate the variance budget.

### Defect 3: linear correction has small residual bias

**Location**: `src/l2/compute.rs:799-818` (correction constants).

**Issue**: Linear correction has residual bias $r^2(2r^2-1)/(d-2)$ per
pair. At $r^2 = 0.9$, $d = 200$: residual ≈ +0.0036 per pair.

**Severity**: small. At the LD-score level the residual is ~0.3%, well
below the LDSC regression noise floor. **Probably not worth fixing**
unless strict per-pair unbiasedness is needed; the quadratic correction
in §5 would address it at ~5% wall-clock cost.

## 9. Concerns and unverified claims

**Concern A (real)**: Variance bound assumes Gaussian-like 4th moments
($K_{jk}/N^2 \ll 1$). Genotype data at low MAF can have $K_{jj}/N$ tens of
times larger. The accuracy bound $\mathrm{Var}(\tilde r) \approx (1+r^2)/d$
is loose for rare SNPs. The preprint should caveat this.

**Concern B (real, structural)**: CountSketch with a single hash has
$O(1/d)$ variance but only **Chebyshev-style tail bounds**. Individual
estimates can be wildly off with non-negligible probability. For
LD-score regression this is mitigated by summing many pairs, but
catastrophic outliers can occur.

**Concern C (real, conceptual)**: All sketched correlations in a window
share the same hash $S$, so they are correlated. The LD-score variance
is **not** the sum of per-pair variances. Empirical observation: at
moderate $r^2$ the correlation is weakly destructive (LD-score variance
is close to the sum of per-pair variances), but for windows dominated by
a single high-LD region the correlation could amplify systematic biases.

## 10. Alternative algorithms worth considering

### Option A: Median-of-k / Mean-of-k CountSketch (**implemented + rejected, 2026-05-12**)

Originally proposed: replace one $d$-dim sketch with $k$ independent
$(d/k)$-dim sketches and take the median of the $k$ estimates. Promised:
total memory and scatter cost unchanged, concentration improves from
Chebyshev (polynomial tails) to Hoeffding (exponential tails).

**Outcome: implemented as `--sketch-k <K>` on the `experiment/median-of-k`
branch (commits `364897e` and `5e2c4dd`), with both median and arithmetic-
mean aggregators tested. Both rejected, but for different reasons.**

#### Why median-of-K fails: the $\tilde r^2$ distribution is right-skewed at low $r^2$

The Alon-Matias-Szegedy median trick assumes the per-estimate noise is
symmetric around the mean. CountSketch's $\tilde r^2$ has a strongly
right-skewed distribution at low $r^2$. Monte Carlo on N=2,490 standardized
vectors with controlled true $r^2$ (`docs/countsketch-mc-distribution.py`,
`docs/countsketch-mc-aggregator-bias.py`):

| true $r^2$ | distribution skewness | $\mathbb{E}[X] - \mathrm{median}(X)$ |
|-----------:|---------------------:|-------------------------------------:|
| 0.0   | **+2.54** | +0.0032 |
| 0.01  | +1.74     | +0.0043 |
| 0.1   | +0.32     | +0.0014 |
| 0.3   | +0.02     | -0.0001 |
| 0.7   | -0.37     | -0.0023 |

The skewness is most severe at $r^2 = 0$ — exactly the regime that
dominates LD-score sums under rapid LD decay (most pairs in a 1 Mb window
are at $r^2 \approx 0$).

Median-of-K convergence: for K iid samples from a right-skewed
distribution, $\mathbb{E}[\text{median of K}]$ converges to the
*population median*, which is strictly less than the *population mean*.
Since the quadratic correction inverts the bias structure of the *mean*
($r^2 + (1-r^2)(1-2r^2)/d$), applying it to a median that's below the
mean leaves a residual systematic downward bias.

Per-pair median-of-3 corrected bias at $d_\text{sub} = 67$:

| true $r^2$ | median-of-K bias | mean-of-K bias |
|-----------:|-----------------:|---------------:|
| 0.0   | **-0.0052** | +0.0002 |
| 0.01  | -0.0054     | +0.0006 |
| 0.1   | -0.0070     | +0.0023 |
| 0.3   | -0.0008     | +0.0031 |
| 0.7   | +0.0035     | +0.0022 |

Summed over a typical chr22 1 Mb LD window (~2,400 pairs, mostly at
$r^2 \approx 0$), the median-of-3 path produces a per-SNP bias of
$\approx -12$ on the LD score. This matches the empirical end-to-end
measurement on chr22 1000G EUR with `--exact-bitpacked` as the reference
(mean L2 = 16.08):

| Mode | Pearson r vs exact | mean L2 | total bias |
|------|-------------------:|--------:|-----------:|
| `--sketch 200` (K=1)                    | 0.979 | 14.97 | -1.10  |
| `--sketch 450 --sketch-k 3` (d_sub=150) | 0.985 | 10.67 | -5.41  |
| `--sketch 300 --sketch-k 3` (d_sub=100) | 0.976 |  8.20 | -7.88  |
| `--sketch 201 --sketch-k 3` (d_sub=67)  | 0.952 |  4.85 | -11.23 |

Per-pair bias × window size ≈ measured bias, to ~0.5 LD-score units.

#### Why mean-of-K *doesn't help either*: variance at low $r^2$

Replacing median with arithmetic mean (also implemented on the experiment
branch, behind `LDSC_SKETCH_AGG=mean|median` env var) eliminates the
skewness bias: mean of K iid unbiased corrections is unbiased. But mean-
of-K at the same *total* budget DIM = K × d_sub is **not statistically
equivalent to single-sketch d=DIM**, because the variance bound
$\mathrm{Var}[\tilde r^2] \approx (1+r^2)/d$ is asymptotic in $d$ and
loose at $r^2 = 0$ for small $d$.

Empirical variance at fixed (X_j, X_k), 1,500 sketches each
(`docs/countsketch-mc-variance-check.py`):

| true $r^2$ | single d=200 SD | mean-of-3 d_sub=67 SD | ratio | single-d=67/√3 prediction |
|-----------:|----------------:|----------------------:|------:|--------------------------:|
| 0.0  | 0.00894 | **0.01416** | 1.58× | 0.01374 ✓ |
| 0.1  | 0.04028 | 0.03936     | 0.98× | 0.03991 ✓ |
| 0.3  | 0.05378 | 0.05460     | 1.02× | 0.05365 ✓ |

At $r^2 \geq 0.1$ the equivalence holds to within 2%, but at $r^2 = 0$
mean-of-3 has 1.58× higher SD than single-sketch d=200. The variance of
mean-of-K matches the canonical iid prediction $\mathrm{SD}_\text{single-d=67} / \sqrt{K}$
almost exactly — the loss vs single-d=200 isn't a bug; it's because
single-sketch d=67 has more than 3× the variance of single-d=200 at
$r^2 = 0$, breaking the "same total budget" claim.

End-to-end empirics with mean-of-K aggregator:

| Mode | Pearson r vs exact | mean L2 |
|------|-------------------:|--------:|
| `--sketch 200` (K=1)                                  | 0.979 | 14.97 |
| `--sketch 450 --sketch-k 3` (d_sub=150) + mean        | 0.978 | 14.99 |
| `--sketch 300 --sketch-k 3` (d_sub=100) + mean        | 0.963 | 14.37 |
| `--sketch 201 --sketch-k 3` (d_sub=67)  + mean        | 0.926 | 13.85 |

Mean-of-K is unbiased (mean L2 stays at the single-sketch baseline
~15.0), but Pearson r degrades as d_sub shrinks because each sub-sketch
contributes more variance than the asymptotic formula predicts.

#### Conclusion

Neither aggregator beats single-sketch d=DIM at fixed total compute
budget for LDSC:

- **Median-of-K**: systematically downward-biased in the regime that
  dominates LD-score sums. Catastrophic on real data (-11 LD-score
  units at K=3 d_sub=67).
- **Mean-of-K**: unbiased, but variance at low $r^2$ is sqrt(K)× higher
  than single-sketch d=DIM at the same compute. Pearson r drops; the
  user's path to higher accuracy is "increase $d$", not "add K".

The "Hoeffding concentration" promise of median-of-K applies only to
symmetric noise. For the LDSC use case the noise is bimodal in
skewness sign (right at low $r^2$, left at high $r^2$) and dominated by
the right-skew at low $r^2$ where most pair-contributions live.

Both implementations preserved on the `experiment/median-of-k` branch
along with the MC scripts under `docs/countsketch-mc-*.py` (the
experiment branch is the only place those scripts live).

### Option B: Bit-packed exact computation (**implemented + rejected, 2026-05-11**)

Genotype data is bit-packed in BED format (2 bits per genotype, 4 per
byte). For pair $(j, k)$, the inner product $X_j^\top X_k$ can be
computed via SIMD popcount and shuffle operations on raw BED bytes —
no f32 expansion. Compute the joint 3×3 genotype contingency table by
counting 2-bit-pair matches, then compute r² from the 9 counts.

This is what PLINK 2 does. Original estimate (this section, pre-implementation):
~4–8× faster than dense f32 GEMM and ~2× faster than `--sketch 200`,
making CountSketch obsolete for the LDSC use case.

**Outcome: this estimate was wrong.** Implemented as `--exact-bitpacked`
(Phase 1 scalar + Phase 2 AVX2 + NEON histogram kernels + intra-chr
parallelism) and benchmarked on AWS Batch (c6a.4xlarge, EPYC 7R13, full
1000G phase3 N=2,490, M=1.66M, --ld-wind-kb 1000, 3 runs each):

| Mode                                  | Time (mean ± σ) | vs sketch-200 | vs exact-f64 |
|---------------------------------------|-----------------|---------------|--------------|
| `--exact-bitpacked --mmap`            | 105.4s ± 4.5s   | 25.4× slower  | 2.6× slower  |
| `--snp-level-masking` (exact-f64)     | 41.2s ± 2.0s    | 9.9× slower   | 1.0×         |
| `--snp-level-masking --fast-f32`      | 16.9s ± 0.6s    | 4.1× slower   | 2.4× faster  |
| `--sketch 200`                        | 4.15s ± 0.03s   | 1.0×          | 9.9× faster  |

The "~1.6 × 10¹³ byte ops → ~16s" estimate assumed AVX2 popcount at
~1 TB/s throughput; actual per-pair cost is dominated by fixed overhead
(joint-code construction + 16-bin histogram + scatter to `l2[j]`,
`l2[k]`) that doesn't amortize the way GEMM's c² output entries per
chunk do. Per-pair on paper: GEMM ≈ N/8 cycles (one AVX2 FMA per
individual, fully pipelined) ≈ 311 cycles at N=2,490; bit-packed ≈ 16
chunks × ~130 cycles + scatter ≈ 2,200 cycles. ~7× per-pair gap, ~2.6×
wall after parallelism amortizes some of it. See `docs/perf-log.md`
2026-05-11 entry for the full analysis.

**Implication for this doc**: CountSketch is **not** obsolete.
`--sketch 200` remains the fastest path that produces statistically
useful LD scores, and `--snp-level-masking --fast-f32` remains the
fastest path for bit-stable LD scores. The bit-packed code has been
removed from main (commit `287902b`) and preserved on the
`experiment/bitpacked-exact` branch.

### Option C: int8 / lower-precision GEMM

The r² formula is robust to ~1% per-pair noise. AVX-512 VNNI instructions
do int8 matmul at several TOPS. ldsc-rs already has `--fast-f32`; a
`--fast-i8` would be ~4× faster than f32 with negligible accuracy loss.
faer or candle probably already supports int8 GEMM.

Less innovative than bit-packing but easier to integrate. Worth doing
as a low-hanging fruit if bit-packing is too ambitious.

### Option D: Structural compression

The current f32 GEMM path materializes an N×c f32 chunk before the GEMM,
which is 40 MB at N=50K, c=200 — spills out of L2/L3. The bit-packed
approach in Option B avoids this entirely. Otherwise consider:
- Streaming GEMM that processes one column at a time
- Packed-half-float (f16) storage with f32 accumulation (Intel AMX has
  native bf16 ops)

## 11. Recommendations

**Tier 1 (done 2026-05-11)**: fix the two comment defects in
`src/l2/compute.rs:803-826`.

**Tier 2 (done 2026-05-11)**: implement the quadratic correction
as the bias-correction strategy for `--sketch`. The earlier linear
fused-form correction was removed in the same change — it had a
known O(1/d) residual that the quadratic eliminates, and AWS
profiling showed no measurable wall-clock cost. d ≤ 50 emits a
runtime warning since the Taylor expansion breaks down there.

**Tier 3 (tested + rejected 2026-05-11)**: bit-packed exact computation
(Option B) was implemented as `--exact-bitpacked`, benchmarked on AWS,
and found to be 25× slower than `--sketch 200` and 2.6× slower than the
existing exact-f64 GEMM path at 1000G scale. The predicted ~2× speedup
over CountSketch did not materialize; the original throughput estimate
assumed per-byte AVX2 popcount would dominate, but per-pair fixed
overhead (joint-code construction + 16-bin histogram + scatter) is the
actual bottleneck. **CountSketch is NOT obsolete.** Implementation
removed from main (commit `287902b`) and preserved on the
`experiment/bitpacked-exact` branch for future reference (validation
oracle on adversarial datasets, or as a starting point for a PLINK-2-
style bit-plane reformat).

**Tier 4 (tested + rejected 2026-05-12)**: median-of-k AND mean-of-k
CountSketch (Option A) were both implemented as `--sketch-k <K>` on the
`experiment/median-of-k` branch, with the aggregator chosen at runtime
via `LDSC_SKETCH_AGG=mean|median`. K=1 byte-identical to current
`--sketch`; K>1 wall-clock scales K× (confirming the predicted
per-pair cost). Neither aggregator beats single-sketch `--sketch d`:
- *Median* is severely downward-biased because the per-sketch
  $\tilde r^2$ distribution is right-skewed at low $r^2$ (where most
  LD-window pairs live; skewness +2.5 at $r^2=0$). End-to-end bias on
  chr22 EUR: -11 LD-score units at K=3 d_sub=67.
- *Mean* eliminates the bias but has $\sqrt{K}$ × higher variance than
  single-sketch d=DIM at low $r^2$, because the asymptotic
  $(1+r^2)/d$ variance bound is loose at small $d$ — Pearson r vs
  exact drops from 0.979 (K=1 d=200) to 0.926 (K=3 d_sub=67) even
  though the mean is unbiased.

The right way to get higher LD-score accuracy is to increase
`--sketch d`, not split into K sub-sketches. See §10 Option A for the
full empirical analysis and the MC scripts on the `experiment/median-of-k`
branch.

## 13. Empirical comparison: quadratic vs linear correction (chr22 EUR, historical)

Before removing the linear correction (2026-05-11), ran both modes on
chr22 1000G EUR (N=503, M=18 627 after MAF≥0.05) at four sketch
dimensions. The linear correction is no longer in the codebase; this
section records the empirical evidence that motivated removing it.

### Aggregate metrics

| d | mode | mean L2 | mean bias | std bias | mean \|bias\| | Pearson r |
|---:|---|---:|---:|---:|---:|---:|
| (exact) | — | 18.849 | — | — | — | 1.0 |
| 50 | linear | 18.990 | +0.141 | 7.26 | 5.57 | 0.9159 |
| 50 | quadratic | 19.009 | +0.160 | 7.33 | 5.63 | 0.9138 |
| 100 | linear | 18.439 | −0.410 | 4.14 | 3.18 | 0.9683 |
| 100 | quadratic | 18.474 | −0.375 | 4.16 | 3.19 | 0.9678 |
| 200 | linear | 18.565 | −0.284 | 2.68 | 2.05 | 0.9866 |
| 200 | quadratic | 18.587 | −0.261 | 2.69 | 2.05 | 0.9865 |
| 500 | linear | 19.198 | +0.350 | 1.79 | 1.34 | 0.9940 |
| 500 | quadratic | 19.209 | +0.360 | 1.80 | 1.35 | 0.9940 |

| d | Δ\|mean bias\| (lin − quad) | Δstd (quad − lin) |
|---:|---:|---:|
| 50 | −0.019 (linear better) | +0.07 (quadratic noisier) |
| 100 | +0.036 (quadratic better) | +0.02 |
| 200 | +0.022 (quadratic better) | +0.005 |
| 500 | −0.011 (linear better, in MC noise) | +0.002 |

**Findings**:
- At the recommended d=200, quadratic reduces mean bias by 0.022
  (~0.12% relative) at the cost of marginally higher variance. MSE is
  essentially identical.
- At d=50 and d=500, the means are within MC noise (one sketch seed; a
  multi-seed average would be needed for tighter comparison).
- Pearson r is virtually identical between modes at all d.

### Per-bin analysis at d=200 (the recommended setting)

The per-bin breakdown reveals the quadratic correction *is* doing
predicted, signal-bearing work — even though the aggregate effect is
small.

| L2 bin | n | mean L2_exact | mean err_lin | mean err_quad | mean (q − l) shift |
|---|---:|---:|---:|---:|---:|
| (0, 5] | 1547 | 3.7 | −0.244 | −0.223 | **+0.021** |
| (5, 10] | 4558 | 7.5 | −0.278 | −0.253 | **+0.025** |
| (10, 20] | 6615 | 14.4 | −0.425 | −0.398 | **+0.027** |
| (20, 40] | 4022 | 27.2 | −0.431 | −0.405 | **+0.026** |
| (40, 80] | 1771 | 53.5 | +0.129 | +0.127 | −0.002 |
| (80, 200] | 114 | 103.1 | +5.953 | +5.818 | **−0.135** |

The math predicts that the linear correction's residual bias is
$r^2(2r^2 - 1)/(d-2)$:
- For windows dominated by **small r²** (low-LD-score SNPs), residual is
  negative, so linear under-shoots and quadratic correctly shifts upward
  (+0.02 to +0.03).
- For windows dominated by **large r²** (high-LD-score SNPs), residual
  is positive, so linear over-shoots and quadratic correctly shifts
  downward (−0.135).

The 114 SNPs in the (80, 200] L2 bin are precisely the ones where the
linear correction's bias is mathematically expected to be largest. The
quadratic correction reduces their error by ~2.3% (from +5.95 to +5.82).
Top 10 by shift magnitude are all in a single high-LD region (rs12484694
neighborhood, exact L2 ≈ 110), all consistently shifted by ≈ −0.21.

### Verdict

The quadratic correction is **mathematically more correct in the
predicted direction**: per-bin shifts match the leading-order bias
formula. The aggregate MSE improvement is small (within sketch
sampling variance at practical d values), but the bias-shift signal at
high-LD-score SNPs — the ones that drive LDSC regression weight — is
real and matches the formula prediction.

After this empirical confirmation, the linear correction was removed
from the codebase entirely. AWS biobank-scale profiling showed the
quadratic correction adds no detectable wall-clock cost, so there is
no tradeoff: quadratic is unconditionally better on correctness and
empirically equivalent on speed.

**Note on d ≤ 50**: at very low sketch dimensions, the Taylor expansion
underlying the quadratic inversion breaks down (higher-order bias
terms in the $1/d^2$, $1/d^3$ expansion become non-negligible, and the
quadratic inverse amplifies sketch noise). At d=50 the corrected
estimator is empirically no better than uncorrected at the LD-score
level. ldsc-rs prints a runtime warning when `--sketch d` is given
with d ≤ 50; users should pass d ≥ 100 (cost is essentially flat in d
below the GEMM crossover, so there is no reason to go lower).

## 12. Reproducibility

The Monte Carlo simulations underlying every empirical number in this
doc are reproducible. Two Python scripts (using
`/tmp/ldsc_py3_venv/bin/python3` with numpy):

### Simulation 1: bias of sketched estimators

```python
import numpy as np
rng = np.random.default_rng(42)
N = 2000; n_trials = 30000

def make_data(N, r, rng):
    z1 = rng.standard_normal(N); z2 = rng.standard_normal(N)
    x = z1; y = r*z1 + np.sqrt(max(0, 1-r*r))*z2
    return (x - x.mean())/x.std(), (y - y.mean())/y.std()

def cs(X, d, rng):
    b = rng.integers(0, d, N); s = rng.choice([-1, 1], N).astype(float)
    out = np.zeros((d, X.shape[1])); np.add.at(out, b, s[:, None]*X); return out

for r in [0.0, 0.3, 0.5, 0.7, 0.9]:
    x, y = make_data(N, r, rng); r_emp = (x @ y)/N; r2_emp = r_emp**2
    XY = np.stack([x, y], axis=1)
    for d in [50, 100, 200, 500]:
        bias_renorm = 0.0
        for _ in range(n_trials):
            t = cs(XY, d, rng)
            A = (t[:, 0]*t[:, 1]).sum()
            B = (t[:, 0]**2).sum(); C = (t[:, 1]**2).sum()
            bias_renorm += A*A/(B*C) - r2_emp
        bias_renorm /= n_trials
        pred = (1 - 3*r2_emp + 2*r2_emp**2) / d
        print(r, d, bias_renorm, pred)
```

### Simulation 2: LD-score level

See `docs/countsketch-math-analysis.py` (committed alongside this doc)
for the full LD-score simulation with M=200 neighbors and exponential
LD decay.

## Sources

- ldsc-rs CountSketch implementation: `src/l2/compute.rs` lines 166-725
- CountSketch original: Charikar, Chen, Farach-Colton, "Finding Frequent
  Items in Data Streams", ICALP 2002
- Variance bounds for sparse projections: Li, Hastie, Church, "Very
  Sparse Random Projections", KDD 2006
- Median-of-k construction: Alon, Matias, Szegedy, "The Space Complexity
  of Approximating the Frequency Moments", STOC 1996
- PLINK 2 bit-packed inner products: Chang et al., "Second-generation
  PLINK: rising to the challenge of larger and richer datasets",
  GigaScience 2015
- Source of the unbiased $r^2$ correction $\hat r^2 - (1-\hat r^2)/(N-2)$
  used by LDSC: Yin & Fan 2001, "Estimating $R^2$ Shrinkage in Multiple
  Regression," J Experimental Education 69:203-224 (cited by
  Bulik-Sullivan 2015 supplementary note)
- LD-genetics-rigorous alternative unbiased $r^2$ estimator:
  Ragsdale & Gravel 2020, "Unbiased Estimation of Linkage Disequilibrium
  from Unphased Data," Mol Biol Evol 37:923,
  doi:10.1093/molbev/msz265
