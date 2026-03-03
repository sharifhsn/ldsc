# LD scores and matrices

## Overview[#](https://pan.ukbb.broadinstitute.org/docs/ld#overview "Direct link to heading")

We computed in-sample dosage-based LD matrices and scores for each of six ancestry group in UKBB. [LD matrices](https://pan.ukbb.broadinstitute.org/docs/ld#ld-matrices) are available in Hail's [BlockMatrix](https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html) format on Amazon AWS (see details [here](http://pan.ukbb.broadinstitute.org/docs/hail-format)). [LD scores](https://pan.ukbb.broadinstitute.org/docs/ld#ld-scores) are available in [LDSC](https://github.com/bulik/ldsc)\-compatible flat files (`.l2.ldscore.gz` and `.M_5_50`) [here](https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.ALL.ldscore.tar.gz). For large-scale analysis, you can also find a full LD score [Hail Table](https://hail.is/docs/0.2/hail.Table.html) (not restricted to the HapMap3 variants) on Amazon AWS (see details [here](http://pan.ukbb.broadinstitute.org/docs/hail-format))

For LD computation, please find technical details below. All the code is also publicly available [here](https://github.com/atgu/ukbb_pan_ancestry/blob/master/compute_ld_matrix.py). Detailed instruction for how to run LD score regression is available on [LDSC's website](https://github.com/bulik/ldsc/wiki).

## LD matrices[#](https://pan.ukbb.broadinstitute.org/docs/ld#ld-matrices "Direct link to heading")

-   The dosage-based genotype matrix XXX was column-wise mean-centered and normalized.
-   We applied the same variant QC filter used for the Pan-UKB GWAS (INFO > 0.8, MAC > 20 in each population; see details [here](https://pan.ukbb.broadinstitute.org/docs/qc#variant-qc))
-   For covariate correction, the residuals from the regression of genotype∼covariatesgenotype \\sim covariatesgenotype∼covariates were obtained via Xadj\=McXX\_{adj} = M\_cXXadj​\=Mc​X where Mc\=I−C(CTC)−1CTM\_c = I - C(C^TC)^{-1}C^TMc​\=I−C(CTC)−1CT, the residual-maker matrix, and CCC is the matrix of covariates.
-   We used the same covariates used for the Pan-UKB GWAS, namely ageageage, sexsexsex, age∗sexage\*sexage∗sex, age2age^2age2, age2∗sexage^2\*sexage2∗sex, and the first 10 PCs of the genotype matrix (see details [here](https://pan.ukbb.broadinstitute.org/docs/qc#gwas-model)).
-   We then computed LD matrix RRR via R\=XadjTXadjnR = \\frac{X\_{adj}^TX\_{adj}}{n}R\=nXadjT​Xadj​​ with a radius of **10 Mb**. Each element r^jk\\hat{r}\_{jk}r^jk​ of RRR represents the Pearson correlation coefficient of genotypes between variant jjj and kkk.
-   For X-chromosome, we computed a LD matrix jointly using both males and females where male genotypes are coded 0/1 and female genotypes are coded 0/1/2.

## LD scores[#](https://pan.ukbb.broadinstitute.org/docs/ld#ld-scores "Direct link to heading")

-   To account for an upward bias of the standard estimator of the Pearson correlation coefficient, we applied a bias adjustment for r^jk2\\hat{r}^2\_{jk}r^jk2​ using r~jk2\=n−1n−2r^jk2−1n−2\\tilde{r}^2\_{jk} = \\frac{n-1}{n-2}\\hat{r}^2\_{jk} - \\frac{1}{n-2}r~jk2​\=n−2n−1​r^jk2​−n−21​.
-   LD scores for variant jjj were subsequently computed via lj\=∑kr~jk2l\_j = \\sum\_k \\tilde{r}^2\_{jk}lj​\=∑k​r~jk2​ with a radius of **1 MB**.
-   For LDSC-compatible flat files, we only exported LD scores of high-quality HapMap 3 variants that are 1) in autosomes, 2) not in the MHC region, 3) biallelic SNPs, 4) with INFO > 0.9, and 5) MAF > 1% in UKB and gnomAD genome/exome (if available).
-   We note that, since we applied covariate adjustment above, these LD scores are equivalent to the covariate-adjusted LD scores as described in [Luo, Y. & Li, X. et al., 2020](https://www.biorxiv.org/content/10.1101/503144v4)