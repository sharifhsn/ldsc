*Nat. Genet.* **47, 291–295 (2015)**

# LD Score regression distinguishes confounding from polygenicity in genome-wide association studies

**Brendan K Bulik-Sullivan, Po-Ru Loh, Hilary K Finucane, Stephan Ripke, Jian Yang, Schizophrenia Working Group of the Psychiatric Genomics Consortium, Nick Patterson, Mark J Daly, Alkes L Price & Benjamin M Neale**

In the version of this supplementary file originally posted online on 2 February 2015, the institutional affiliation of Benedicto Crespo-Facorro was incorrectly given as University Hospital Marqués de Valdecilla, Instituto de Formación e Investigación Marqués de Valdecilla, University of Cantabria, E!39008 Santander, Spain. The affiliation should have read Department of Psychiatry, University Hospital Marqués de Valdecilla, School of Medicine, University of Cantabria–IDIVAL-CIBERSAM, Santander, Spain. The error has been corrected in this file as of 27 July 2015.

## SUPPLEMENTARY NOTE

## 1. LD Score in an Unstructured Sample

1.1. Model. We model phenotypes as generated from the equation

$$\phi = X\beta + \epsilon$$

where  $\phi$  is an  $N \times 1$  vector of (quantitative) phenotypes, X is an  $N \times M$  matrix of genotypes normalized to mean zero and variance one<sup>1</sup> (we ignore the distinction between normalizing and centering in our sample and in the population since the error so introduced has expectation zero and  $\mathcal{O}(1/N)$  variance),  $\beta$  is an  $M \times 1$  vector of per-normalized-genotype effect sizes and  $\epsilon$  is an  $N \times 1$  vector of environmental effects. We describe a model where all three variables on the right side of equation 1.1 are random. In this model,  $\mathbb{E}[\epsilon] = 0$ ,  $\mathrm{Var}[\epsilon] = (1 - h_g^2)I$ ,  $\mathbb{E}[\beta] = 0$  and  $\mathrm{Var}[\beta] = (h_g^2/M)I$ . To model genotypes, we assume that the genotype at variant j for individual i is independent of other individuals' genotypes, but we do incorporate linkage disequilibrium into the model: define  $r_{jk} := \mathbb{E}[X_{ij}X_{ik}]$ , which does not depend on i. Finally, we assume that X,  $\beta$  and  $\epsilon$  are mutually independent. We will relax the assumption that environmental effects are independent of genotype when we model population stratification in §2.2.

1.2. Relationship between LD and  $\chi^2$ -Statistics. For each variant j = 1, ..., M, we compute least-squares estimates of effect size  $\hat{\beta}_j := X_j^{\mathsf{T}} \phi/N$  (where  $X_j$  denotes the  $N \times 1$  vector of genotypes at variant j) and  $\chi^2$ -statistics  $\chi_j^2 := N \hat{\beta}_j^2$ . In this section, we compute  $\mathbb{E}[\chi_j^2]$  with the expectation taken over random X,  $\beta$ , and  $\epsilon$ .

**Proposition 1.** Define the LD Score of variant j as

(1.2) 
$$\ell_j := \sum_{k=1}^M r_{jk}^2.$$

Under the model described in §1.1, the expected  $\chi^2$ -statistic of variant j is

(1.3) 
$$\mathbb{E}[\chi_j^2] \approx \frac{Nh_g^2}{M} \ell_j + 1.$$

1

<sup>&</sup>lt;sup>1</sup>Note that the normalization to variance hides an implicit assumption that rare SNPs have larger effect sizes. We show via simulation in the main text that LD Score regression, as an inference procedure, is not particularly sensitive to this assumption.

<sup>&</sup>lt;sup>2</sup>This is almost the same assumption made in [1] (the difference is that they assume *i.i.d.* per-normalized genotype effect sizes for *genotyped* SNPs, and we make this assumption for all SNPs). If one wishes to specify a different variance structure for the per-normalized-genotype effect sizes, *e.g.*,  $\operatorname{Var}[\beta_j] = f_j$ , then all results presented herein hold with normalized genotypes  $(G_{ij} - 2p_j)/\sqrt{f_j}$  replacing the usual  $(G_{ij} - 2p_j)/\sqrt{2p_j(1-p_j)}$ , where  $G_{ij}$  denotes additively coded (0,1,2) genotypes.

*Proof.* Since  $\mathbb{E}[\hat{\beta}_j] = 0$ , observe that  $\mathbb{E}[\chi_j^2] = N \cdot \text{Var}[\hat{\beta}_j]$ . We will obtain the variance of  $\hat{\beta}_j$  via the law of total variance:

(1.4) 
$$\operatorname{Var}[\hat{\beta}_{j}] = \mathbb{E}[\operatorname{Var}[\hat{\beta}_{j} \mid X]] + \operatorname{Var}[\mathbb{E}[\hat{\beta}_{j} \mid X]]$$
$$= \mathbb{E}[\operatorname{Var}[\hat{\beta}_{j} \mid X]],$$

where the second line follows from the fact that  $\mathbb{E}[\hat{\beta}_j \mid X] = 0$ , irrespective of X. First,

(1.5) 
$$\operatorname{Var}[\hat{\beta}_{j} \mid X] = \frac{1}{N^{2}} \operatorname{Var}[X_{j}^{\mathsf{T}} \phi \mid X]$$
$$= \frac{1}{N^{2}} X_{j}^{\mathsf{T}} \operatorname{Var}[\phi \mid X] X_{j}$$
$$= \frac{1}{N^{2}} \left( \frac{h_{g}^{2}}{M} X_{j}^{\mathsf{T}} X X^{\mathsf{T}} X_{j} + N(1 - h_{g}^{2}) \right).$$

We can write the term on the left in terms of more familiar quantities as

(1.6) 
$$\frac{1}{N^2} X_j^{\mathsf{T}} X X^{\mathsf{T}} X_j = \sum_{k=1}^M \tilde{r}_{jk}^2,$$

where  $\tilde{r}_{jk} := \frac{1}{N} \sum_{i=1}^{N} X_{ij} X_{ik}$  denotes the sample correlation between additively-coded genotypes at variants j and k. Since

(1.7) 
$$\mathbb{E}[\hat{r}_{jk}^2] \approx r_{jk}^2 + (1 - r_{jk}^2)/N,$$

(where the approximation sign hides terms of order  $\mathcal{O}(1/N^2)$  and smaller; one can obtain this approximation via e.g., the  $\delta$ -method),

(1.8) 
$$\mathbb{E}\left[\sum_{k=1}^{M} \tilde{r}_{jk}^{2}\right] \approx \ell_{j} + \frac{M - \ell_{j}}{N}.$$

Thus,

(1.9) 
$$\mathbb{E}[\chi_j^2] \approx \frac{N(1-1/N)h_g^2}{M}\ell_j + 1$$
$$\approx \frac{Nh_g^2}{M}\ell_j + 1,$$

Values of N (study sample size) considered in the main text generally fall between  $10^4$  and  $10^5$ , so the approximation  $1 - 1/N \approx 1$  is appropriate.

# 2. LD Score with Population Stratification

2.1. Model of Population Structure. We model population structure induced by genetic drift in a mixture of two populations in equal proportions as follows: we draw a matrix of normalized genotypes X consisting of N/2 samples from population 1 and N/2 samples from population 2 (we will use the notation  $i \in P_m$  for  $m \in \{1,2\}$  to denote that individual i is a member of population m), subject to the following constraints:  $\operatorname{Var}[X_{ij}] = 1$ ,  $\mathbb{E}[X_{ij} \mid i \in P_1] = f_j$  and  $\mathbb{E}[X_{ij} \mid i \in P_2] = -f_j$ .

We model the drift term f as  $f \sim N(0, F_{ST}V)$ , where V is a correlation matrix<sup>3</sup> and  $F_{ST}$  is Wright's  $F_{ST}$  [2]. We postpone discussion of the off-diagonal entries of V (which might depend on LD in the ancestral population or recombination rates) until §2.2. Finally, if  $\ell_{j,m}$  denotes the LD Score of variant j in population m, we assume that  $\ell_{j,1} \approx \ell_{j,2} =: \ell_j$ . The last assumption warrants a brief explanation. Assuming approximately equal LD Scores in both populations is certainly not reasonable for very large values of  $F_{ST}$  (e.g., if population 1 and population 2 are from different continents) or in scenarios where one population has passed through a more severe bottleneck than the other (e.g., if population 1 is from Finland and population 2 is from West Africa). However, we are interested in modeling the population stratification that may remain after principal components analysis in GWAS that sample from non-admixed populations, and for this purpose the assumption that  $\ell_{j,1} \approx \ell_{j,2}$  seems reasonable, and is supported by the large values of  $R^2(\ell_{j,m},\ell_{j,n})$  that we observe for all pairs (m,n) of 1000 Genomes European subpopulations.

For reference, typical values of  $F_{ST}$  for human populations are  $\approx 0.1$  for populations from different continents,  $F_{ST} \approx 0.01$  for populations on the same continent, and  $F_{ST} < 0.01$  for subpopulations within the same country.

2.2. **LD** in a Mixture of Populations. Suppose j and k are unlinked variants such that  $r_{jk,1} = r_{jk,2} = 0$  and  $f_j$  is independent of  $f_k$ . In a mixture of populations, it will often hold that j and k will be in LD in the whole population even if they are in equilibrium in both component populations. Let  $r_{mix,jk}$  denote the correlation between SNPs j and k in such a mixture of populations. Conditional on f,

(2.1) 
$$\mathbb{E}[r_{mix,jk} \mid f] = \mathbb{E}[X_{ij}X_{ik} \mid f]$$
$$= \frac{1}{2} (\mathbb{E}[X_{ij}X_{ik} \mid f, i \in P_1] + \mathbb{E}[X_{ij}X_{ik} \mid f, i \in P_2])$$
$$= f_j f_k.$$

If we take the expectation over random  $f_j$  and  $f_j$ , then  $\mathbb{E}[r_{mix,jk}] = 0$ , because  $f_j$  and  $f_k$  are independent with expectation zero. We can use equation 2.1 to compute the variance,

(2.2) 
$$\operatorname{Var}[r_{mix,jk}] = \operatorname{Var}[\mathbb{E}[r_{mix,k} \mid f]] + \mathbb{E}[\operatorname{Var}[r_{mix,jk} \mid f]]$$
$$= \mathbb{E}[f_j^2 f_k^2] + 0$$
$$= \mathbb{E}[f_j^2] \mathbb{E}[f_k^2]$$
$$= F_{ST}^2.$$

Observe that since  $\mathbb{E}[r_{mix,jk}] = 0$ ,  $\text{Var}[r_{mix,jk}] = \mathbb{E}[r_{mix,jk}^2]$ . By equation 1.7, in a finite sample,

(2.3) 
$$\mathbb{E}[\hat{r}_{mix,jk}^2] \approx F_{ST}^2 + (1 - F_{ST}^2)/N.$$

 $<sup>^3</sup>$ In particular, we assume that the diagonal entries of V are all equal, or at least uncorrelated with LD Score. This assumption is unlikely to hold exactly: some parts of the genome drift faster than others, and the rate of drift may be correlated with LD Score (e.g., as a result of linked selection). Nevertheless, our simulations with real population stratification show that this is not likely to be a severe confounder in LD Score regression.

Thus, the sample LD Score is approximately

(2.4) 
$$\mathbb{E}[\tilde{\ell}_j] \approx \ell_j + MF_{ST}^2 + \frac{M(1 - F_{ST}^2)}{N}$$
$$\approx \ell_j + MF_{ST}^2 + \frac{M}{N}.$$

Note that we have ignored the case where j and k are linked and  $V_{jk} \neq 0$ . In this case,  $\mathbb{E}[f_j^2f_k^2] = F_{ST}^2 + 2F_{ST}^2V_{jk}^2$  (from the formula for the double second moments of a multivariate normal distribution). Even if for some variants j, the number of variants k such that  $V_{jk} > 0$  is  $\approx 10^3$ , this will make a negligible difference in  $\mathbb{E}[\tilde{\ell}_j]$ , because  $\sum_{k:V_{jk}>0} 2F_{ST}^2V_{jk}^2 < 2000F_{ST}^2 \ll MF_{ST}^2$  when  $M \approx 10^7$ .

2.3. **Model of Stratified Phenotype.** To model population stratification, we model phenotypes as generated by the equation

$$\phi = X\beta + S + \epsilon,$$

where X is as described in §2.1,  $\beta$  is as described in §1.1 and where S is an environmental stratification<sup>4</sup> term defined by

$$(2.6) S_i := \begin{cases} \sigma_s/2, & i \in P_1 \\ -\sigma_s/2, & i \in P_2. \end{cases}$$

Finally,  $\epsilon$  is as described in §1.1, except  $\operatorname{Var}[\epsilon] = (1 - h_g^2 - \sigma_s^2)$ , which assures that the variance of  $\phi$  in the population is 1<sup>5</sup>. We compute  $\chi^2$ -statistics as defined in §1.1. In this section, we compute  $\mathbb{E}[\chi_j^2]$  with the expectation taken over random X,  $\beta$ ,  $\epsilon$ , f but with S fixed to ensure population stratification.

# 2.4. Relationship between LD and Stratified $\chi^2$ -Statistics.

**Proposition 2.** Under the model described in §2.3, the expected  $\chi^2$ -statistic of variant j is

(2.7) 
$$\mathbb{E}[\chi_j^2] = \frac{Nh_g^2}{M}\ell_j + 1 + aNF_{ST},$$

where a is the expectation of squared difference in mean phenotypes between population1 and population 2.

*Proof.* Since  $\mathbb{E}[\hat{\beta}_j] = 0$ , observe that  $\mathbb{E}[\chi_j^2] = N \cdot \text{Var}[\hat{\beta}_j]$ . We will obtain the variance of  $\hat{\beta}_j$  via the law of total variance:

(2.8) 
$$\operatorname{Var}[\hat{\beta}_j] = \mathbb{E}[\operatorname{Var}[\hat{\beta}_j \mid X]] + \operatorname{Var}[\mathbb{E}[\hat{\beta}_j \mid X]].$$

Note that one can calculate f from X, so by conditioning on X we also implicitly condition on f. Unlike in equation 1.4,  $\mathbb{E}[\hat{\beta}_j \mid X] \neq 0$ , because of confounding from population stratification. The inner portion of the first term on the right side of equation 2.8 is the same as in equation 1.5,

<sup>&</sup>lt;sup>4</sup>Environmental population stratification occurs when environmental effects are correlated with ancestry. Genetic population stratification occurs when the alelle frequency of trait-increasing alleles is correlated with ancestry. Our model includes environmental stratification and a small amount of genetic stratification from drift. Stronger genetic stratification requires the action of natural selection on the phenotype in question (or a related phenotype). For a more thorough discussion of genetic stratification, see [3].

<sup>&</sup>lt;sup>5</sup>Note that we implicitly require  $1 - h_g^2 - \sigma_s^2 \ge 0$ .

(2.9) 
$$\operatorname{Var}[\hat{\beta}_{j} \mid X] = \frac{1}{N^{2}} \left( \frac{h_{g}^{2}}{M} X_{j}^{\mathsf{T}} X X^{\mathsf{T}} X_{j} + N(1 - h_{g}^{2}) \right).$$

We can take the expectation over random X (and therefore over random f) using the result from equation 2.4. Thus,

(2.10) 
$$\mathbb{E}[\operatorname{Var}[\hat{\beta}_j \mid X]] = \frac{1}{N^2} \left( \frac{h_g^2}{M} \mathbb{E}[X_j^\mathsf{T} X X^\mathsf{T} X_j] + N(1 - h_g^2) \right)$$
$$\approx \frac{h_g^2}{M} \ell_j + h_g^2 F_{ST}^2 + \frac{1}{N}.$$

Next, the inner portion of the second term on the right side of equation 2.8 is

(2.11) 
$$\mathbb{E}[\hat{\beta}_{j} \mid X] = \frac{1}{N} \mathbb{E}[X_{j}^{\mathsf{T}} X \beta + X_{j}^{\mathsf{T}} S + X_{j}^{\mathsf{T}} \epsilon]$$
$$= \frac{1}{N} X_{j}^{\mathsf{T}} S$$
$$= f \sigma_{s}.$$

Since f has variance  $F_{ST}$ ,  $Var[f\sigma_s] = \sigma_s^2 F_{ST}$ . Thus,

(2.12) 
$$\mathbb{E}[\chi_j^2] = N \cdot \operatorname{Var}[\hat{\beta}_j]$$
$$\frac{Nh_g^2}{M} \ell_j + 1 + NF_{ST}(\sigma_s^2 + h_g^2 F_{ST}).$$

We can interpret the final term,  $NF_{ST}(\sigma_S^2 + h_g^2 F_{ST})$ , as  $NF_{ST}$  times the expected squared mean difference in phenotype between populations, which has environmental component  $\sigma_s^2$  and genetic component  $h_g^2 F_{ST}$  (if we model X,  $\beta$  and f as random, there is zero genetic stratification on expectation, but with some small variance about zero). Precisely, if we let  $\bar{\phi}_m$  denote the mean phenotype in population  $m \in \{1, 2\}$ , then

$$(2\mathbb{E}[3\bar{\phi}_1 - \bar{\phi}_2)] = \sigma_s^2 + \sum_{j=1}^M \left[ \mathbb{E}[\beta_j^2] \left( \sum_{i \in P_1} \mathbb{E}[X_{ij}^2 \mid i \in P_1] + \sum_{i \in P_2} \mathbb{E}[X_{ij}^2 \mid i \in P_2] \right) \right]$$

$$= \sigma_s^2 + h_q^2 F_{ST}.$$

Set  $a := \mathbb{E}[(\bar{\phi}_1 - \bar{\phi}_2)^2]$ . Then we have

(2.14) 
$$\mathbb{E}[\chi_j^2] = \frac{Nh_g^2}{M}\ell_j + 1 + aNF_{ST},$$

as desired.  $\Box$ 

# 3. Variance

The results in §2.4 suggest a method for estimating the confounding term  $aNF_{ST}$  from summary statistics: if we regress  $\chi^2$  against LD Score, then the intercept minus one is an estimate of  $aNF_{ST}$ . Because the variance of  $\chi^2$  increases with LD Score, we can improve the efficiency of this estimator by weighting the regression by the reciprocal of the conditional variance function  $\text{Var}[\chi_j^2 | \ell_j]$ . We have derived the conditional expectation  $\mathbb{E}[\chi_j^2 | \ell_j]$  without making distributional assumptions on  $\beta$ 

or  $\epsilon$ ; however, we need stronger assumptions<sup>6</sup> in order to derive the conditional variance: we assume that N is large and  $\beta \sim N(0, h_g^2 I)$ , and  $\epsilon \sim N(0, (1 - h_g^2)I)^7$ . Then,

(3.1) 
$$\hat{\beta}_{j} = \frac{1}{N} (X_{j}^{\mathsf{T}} X \beta + X_{j}^{\mathsf{T}} \epsilon)$$

$$\sim N(0, h_{g}^{2} \ell_{j} / M + h_{g}^{2} / N) + N(0, (1 - h_{g}^{2}) / N)$$

$$\sim N(0, h_{g}^{2} \ell_{j} / M + 1 / N),$$

where the second line follows from a central limit theorem argument and §.1.2. Thus,  $\chi_j^2 = N \hat{\beta}_j^2$  follows a scaled  $\chi^2$  distribution with scale factor  $N h_g^2 \ell_j / M$ , and the conditional variance function is

(3.2) 
$$\operatorname{Var}[\chi_j^2 \mid \ell_j] = \left(1 + \frac{Nh_g^2}{M}\ell_j\right)^2.$$

Note that this is the correct conditional variance function for GWAS with no confounding bias. Since most published GWAS have taken steps to control for population stratification, the most likely use case will be a GWAS with at most a small amount of population stratification.

## 4. Meta-Analysis

Consider a GWAS for a quantitative trait consisting of t sub-studies (all of which sample from the same population) with sample sizes  $N_1, \ldots, N_k$  and total sample size N. For a SNP j, we compute z-scores  $z_{j1}, \ldots, z_{jt}$  (via linear regression:  $z_{j,s} := X_{j,s}^{\mathsf{T}} \phi_s / \sqrt{N_s}$ , where  $X_{j,s}$  is a vector of genotypes for SNP j in study s and  $\phi_s$  is a vector of phenotypes for study s), then perform sample size weighted meta-analysis with single genomic control to obtain test statistics

$$(4.1) z_{j,meta} := \sum_{s=1}^{t} z_{js} \sqrt{\frac{N_s}{\lambda_s N}},$$

and

$$\chi_{j,meta}^2 = z_{j,meta}^2.$$

For meta-analyses without genomic control, set  $\lambda_s = 1$  for all s.

<sup>&</sup>lt;sup>6</sup>Note that the regression weights do not affect the expectation of the parameter estimates, only the standard error. Therefore, if the distributional assumptions that we make in order to derive the conditional variance are violated, it will only increase the standard error. Concretely, if there are very few causal SNPs, or if the distribution of effect sizes is particularly leptokurtotic, then  $\operatorname{Var}[\chi_j^2 \mid \ell_j]$  will increase with  $\ell_j$  faster than the function that we derive in this section, and our estimates will be inefficient.

<sup>&</sup>lt;sup>7</sup>Normality is a stronger assumption than is necessary. We only need that  $\hat{\beta}_j$  be normally distributed, which can hold even if  $\beta$  is not normal, so long as N is large and there are sufficiently many causal SNPs so that  $\frac{1}{N}X_j^{\mathsf{T}}X\beta$  is approximately normal (which would follow from a CLT argument).

**Proposition 3.** Under the model of meta-analysis with genomic control described above, the expected meta-analysis  $\chi^2$ -statistic of variant j is

$$\frac{h_g^2}{MN} \left( \sum_{r=1}^t \sum_{s=1}^t \frac{N_r N_s}{\sqrt{\lambda_r \lambda_s}} \right) \ell_j + \frac{1}{N} \sum_{s=1}^t \frac{N_s}{\lambda_s}.$$

Proof. First, we need to compute  $\mathbb{E}[z_{jr}z_{js}]$ . Under a model where genotypes, environmental effects and SNP-effects are random,  $\mathbb{E}[z_{js}] = 0$  for all s, so  $\mathbb{E}[z_{jr}z_{js}]$  is equal to  $\text{Cov}[z_{jr}, z_{js}]$ . Let X denote the matrix of normalized genotypes in study r and Y the matrix of normalized genotypes in study s. Let  $\delta$  denote the vector of environmental effects in study r and  $\epsilon$  the vector of environmental effects in study s. Suppose that there is no sample overlap, and more generally no cryptic relatedness within or between studies. Further suppose that there is no population stratification. Finally, assume that (to a reasonable approximation), the LD structure in all population from which samples were drawn is approximately equal (though this allows for sampling variance around the population parameter in each study). Then,

$$\operatorname{Cov}[z_{jr}, z_{js}] = \frac{1}{N} \operatorname{Cov}[(X_{j}^{\mathsf{T}} X \beta + X_{j}^{\mathsf{T}} \delta), (Y_{j}^{\mathsf{T}} Y \beta + Y_{j}^{\mathsf{T}} \epsilon)] 
= \frac{1}{N} \left( \operatorname{Cov}[X_{j}^{\mathsf{T}} X \beta, Y_{j}^{\mathsf{T}} Y \beta] + \operatorname{Cov}[X_{j}^{\mathsf{T}} \delta, Y_{j}^{\mathsf{T}} \epsilon] \right).$$

We can evaluate these covariances with the law of total covariance. The term on the right is zero if we assume no sample overlap (because then  $\delta$  and  $\epsilon$  are independent). The term on the left is

(4.5) 
$$\operatorname{Cov}[X_{j}^{\mathsf{T}}X\beta, Y_{j}^{\mathsf{T}}Y\beta \,|\, X, Y] = X_{j}^{\mathsf{T}}X\operatorname{Cov}[\beta, \beta]Y^{\mathsf{T}}Y_{j}$$
$$= \frac{h_{g}^{2}}{M}X_{j}^{\mathsf{T}}XY^{\mathsf{T}}Y_{j}.$$

Removing the conditioning on X and Y,

$$\begin{aligned} \operatorname{Cov}[X_{j}^{\mathsf{T}}X\beta,Y_{j}^{\mathsf{T}}Y\beta] &= & \mathbb{E}[\operatorname{Cov}[X_{j}^{\mathsf{T}}X\beta,Y_{j}^{\mathsf{T}}Y\beta \mid X,Y]] + \operatorname{Cov}[\mathbb{E}[X_{j}^{\mathsf{T}}X\beta,\mid X],\mathbb{E}[Y_{j}^{\mathsf{T}}Y\beta \mid Y]] \\ &= & \mathbb{E}[\operatorname{Cov}[X_{j}^{\mathsf{T}}X\beta,Y_{j}^{\mathsf{T}}Y\beta \mid X,Y]] \\ &= & \frac{h_{g}^{2}}{M}\mathbb{E}[X_{j}^{\mathsf{T}}XY^{\mathsf{T}}Y_{j}] \\ &= & \frac{h_{g}^{2}}{M}\sum_{k=1}^{M}\sqrt{N_{r}N_{s}}\mathbb{E}[\hat{r}_{jk,r}\hat{r}_{jk,s}] \end{aligned}$$

$$(4.6) \qquad = & \frac{\sqrt{N_{r}N_{s}}h_{g}^{2}}{M}\ell_{j},$$

where  $\hat{r}_{jk,r}$  denotes the sample correlation between SNPs j and k in study r. Since we have assumed that the samples in studies r and s are independent,  $\hat{r}_{jk,r}$  and  $\hat{r}_{jk,s}$  are independent estimates of the parameter  $r_{jk}$ , so  $\mathbb{E}[\hat{r}_{jk,r}\hat{r}_{jk,s}] = \mathbb{E}[\hat{r}_{jk,r}]\mathbb{E}[\hat{r}_{jk,s}] =$ 

<sup>&</sup>lt;sup>8</sup>Note that the intercept term will generally be less than one. On the other hand, if application of GC correction in each study was warranted (*i.e.*, if, for each s, all inflation in  $\lambda_s$  reflects confounding), and there is no between-study cryptic relatedness (*e.g.*, sample overlap) then the intercept should be 1.

$$r_{jk}^2$$
. Then,

$$\mathbb{E}[\chi_{j,meta}^{2}] = \mathbb{E}\left[\left(\sum_{s=1}^{t} z_{js} \sqrt{\frac{N_{s}}{\lambda_{s}N}}\right)^{2}\right]$$

$$= \frac{1}{N} \sum_{r=1}^{t} \sum_{s=1}^{t} \sqrt{\frac{N_{r}N_{s}}{\lambda_{r}\lambda_{s}}} \mathbb{E}[z_{js}z_{jr}]$$

$$= \frac{1}{N} \sum_{r=1}^{t} \sum_{s=1}^{t} \frac{N_{r}N_{s}}{\sqrt{\lambda_{r}\lambda_{s}}} \left(\frac{h_{g}^{2}}{M}\ell_{j}\right) + \frac{1}{N} \sum_{s=1}^{t} \frac{N_{s}}{\lambda_{s}}$$

$$= \frac{h_{g}^{2}}{MN} \left(\sum_{r=1}^{t} \sum_{s=1}^{t} \frac{N_{r}N_{s}}{\sqrt{\lambda_{r}\lambda_{s}}}\right) \ell_{j} + \frac{1}{N} \sum_{s=1}^{t} \frac{N_{s}}{\lambda_{s}}.$$

$$(4.7)$$