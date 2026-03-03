# Downloads

The GWAS results are available for download in two main formats:

-   **Per-phenotype flat files**: for most analyses of **one or a few phenotypes**, we suggest using the per-phenotype flat files, available freely on Amazon AWS. More information on the file formats is available in the [Technical Details](https://pan.ukbb.broadinstitute.org/docs/per-phenotype-files).
    -   The phenotype manifest (browse on [Google Sheets](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=903887429) or download on [Amazon AWS](https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz)) contains the location and detailed information of all per-phenotype files for those phenotypes for which GWAS was run.
    -   The variant manifest contains detailed information on each variant (download on [Amazon AWS](https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz), [tbi](https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz.tbi)).
    -   **Please note that p-values are now stored as natural log p-values to avoid underflow (i.e., ln P, not -ln P or -log10 P).**
-   **Hail format**: For large-scale analyses of many phenotypes, we provide the full dataset in [Hail MatrixTable format](https://pan.ukbb.broadinstitute.org/docs/hail-format) on Google Cloud.
-   Please note that the previous iteration of release files have been archived at `s3://pan-ukb-us-east-1/archive_20200615/`.
-   All data are now additionally hosted on the [Genomics Data Lake - Azure Open Datasets](https://learn.microsoft.com/en-us/azure/open-datasets/dataset-panancestry-uk-bio-bank).

In addition, the LD matrices and scores are available in the following formats:

-   **LDSC-compatible flat files**: for running LD score regression, we suggest using the LD score flat files available on Amazon AWS (download the tarball file [here](https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.ALL.ldscore.tar.gz)). More information on the file formats is available on [the LDSC website](https://github.com/bulik/ldsc/wiki).
-   **Hail format**: For large-scale analyses, we provide the full LD matrices and scores in [Hail format](https://pan.ukbb.broadinstitute.org/docs/hail-format) on Amazon AWS.

All heritability estimates (see [here](https://pan.ukbb.broadinstitute.org/docs/heritability) for more information on our approach) are available for download in the following formats:

-   **Flat files**: the manifest flat file is available on AWS ([tarball here](https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/h2_manifest.tsv.bgz)) or on [Google Sheets](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1797288938). Our topline results can be found as part of the main phenotype manifest ([Amazon AWS](https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz) or [Google Sheets](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=903887429))
-   **Hail format**: for large-scale analyses and integration with our other datasets we provide heritability data in [Hail format](https://pan.ukbb.broadinstitute.org/docs/hail-format) on Google Cloud Platform.

The phenotype correlation matrix, used in the construction of the maximally independent set of phenotypes passing QC (see [here](https://pan.ukbb.broadinstitute.org/blog/2022/04/11/h2-qc-updated-sumstats) for details), is available on Amazon AWS (download the tarball file [here](https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/Pheno_pairwise_correlations_update.txt.bgz)).

The ancestry assignments (as well as corresponding principal components and covariates used in our analyses) are available for download through the UK Biobank portal as [Return 2442](https://biobank.ctsu.ox.ac.uk/showcase/dset.cgi?id=2442). These are available to researchers registered with the UK Biobank: refer to instructions within the AMS portal to download these results.