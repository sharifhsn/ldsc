# Incorporating LDSC++ Extensions into the Rust LDSC Engine

## Abstract

Linkage Disequilibrium Score Regression (LDSC) is a foundational method for estimating common variant heritability ($h^2$) and genetic correlations ($r_g$) from genome-wide association study (GWAS) summary statistics. The LDSC++ methodology extends standard LDSC to improve parameter estimation, specifically focusing on multivariate analyses with varying genomic overlaps.

This document outlines the theoretical basis of the five LDSC++ extensions and provides a concrete roadmap for implementing them within the statically-typed, parallelized Rust architecture.

---

## 1. Variable Block-Counts & 2. cM-Based Block Definitions

### Methodological Theory

Standard LDSC relies on block jackknife resampling using a fixed number of contiguous genomic blocks (typically 200). LDSC++ proposes allowing the number of blocks to vary, utilizing definitions based either on a set number of variants or recombination distance in centiMorgans (cM). This better accommodates datasets with heterogeneous missingness patterns and varying numbers of shared genetic variants across trait pairs.

### Implementation Plan

The Rust implementation currently enforces equal-sized blocks via integer division.

* **`cli.rs`**: Introduce arguments to specify block definitions, such as `--jackknife-cm` and `--jackknife-snps`.
* **`regressions.rs`**: The current `drop_ld_meta` function strips out metadata:
```rust
let meta = ["CHR", "BP", "CM"];

```


Modify this to conditionally retain the `CM` column when the recombination distance block definition is active, allowing downstream functions to group variants by genetic distance.
* **`jackknife.rs`**: Refactor the `jackknife` function signature. Instead of accepting `n_blocks: usize`, accept a slice of pre-calculated block boundary indices (`&[(usize, usize)]`).

## 3. Variable Block-Count Sampling (Segmented Regression)

### Methodological Theory

Rather than relying exclusively on block jackknife resampling to estimate the standard error, LDSC++ samples across genome blocks using a segmented (or piecewise) regression approach. The genome-wide covariance parameters are constructed as aggregates of these heterogeneous local parameters, which provides robust estimates without assuming a homoskedastic relationship across the entire genome.

For a partitioned model, the local expected chi-squared statistic for variant $i$ is modeled as:

$$E[\chi^2_i] \approx N_i \sum_k \tau_k \frac{\ell_{ik}}{M_k} + Na + 1$$

In segmented sampling, this equation is evaluated independently per block.

### Implementation Plan

Leveraging the existing zero-allocation Iteratively Re-Weighted Least Squares (IRWLS) implementation and Rayon thread pool, this can be executed efficiently.

* **`jackknife.rs`**: Introduce a new public function `piecewise_sampling()`.
* Use `rayon::prelude::into_par_iter` to map over the block boundaries defined in Extension 2.
* For each block, slice the design matrix `X` and response `y`, and call `irwls::irwls()`.
* Aggregate the resulting `IrwlsResult.est` vectors to compute the genome-wide parameter average and the empirical variance of the local estimates.

## 4. Imputation Quality (INFO) Weighting

### Methodological Theory

Standard LDSC protocols often involve applying a hard filter to remove variants with low imputation quality (e.g., INFO < 0.9). LDSC++ introduces an extended weighting scheme that controls for imputation quality smoothly, preventing the outright loss of information from filtered variants while appropriately down-weighting noisy signals.

### Implementation Plan

The Polars-based streaming pipeline is well-suited for adding this continuous weight.

* **`munge.rs`**: When parsing the summary statistics using `LazyFrame`, identify the imputation quality column (commonly `INFO`). Retain this column in the `.sumstats.gz` output instead of simply using it as a strict filter predicate.
* **`regressions.rs`**: Load the `INFO` column using the `extract_f64(&merged, "INFO")` helper. When calculating the initial regression weights, multiply the baseline heteroskedasticity weight by the INFO score to penalize variants with lower confidence.

## 5. Adjusted Correlation Correction Weighting

### Methodological Theory

LDSC weights correct for the heteroskedasticity of the regression and for correlations between neighboring test statistics. LDSC++ adjusts this weighting scheme to explicitly correct for the correlation between the LD score and the association statistic. This ensures that extreme LD values (close to zero, often due to noise) do not exert disproportionate leverage on the regression slope.

The standard weight computation incorporates the reference LD score ($\ell_{i}$) and the sample size ($N_i$). The variance parameter $v$ is typically approximated as:

$$v = \max\left(\frac{w_{\ell 2, i}}{M} + \frac{1}{N_i}, 10^{-9}\right)$$

The final weight $w_i = \frac{1}{2 v^2}$. LDSC++ applies a transformation to this term.

### Implementation Plan

This is an isolated algorithmic swap that can be gated behind a command-line flag.

* **`cli.rs`**: Add a boolean flag, e.g., `--ldscplus-weights`.
* **`regressions.rs`**: Inside the regression drivers (`run_h2_scalar`, `fit_h2_partitioned`, `run_rg`), locate the weight calculation loops:
```rust
for i in 0..n_obs {
    let v = (w_l2[i] / m_snps + 1.0 / n_vec[i]).max(1e-9);
    weights[i] = 1.0 / (2.0 * v * v);
}

```


Introduce a conditional branch here. If the LDSC++ weighting flag is active, calculate $w_i$ using the adjusted correlation correction formula provided in the LDSC++ methodology.

---

*End of Document*