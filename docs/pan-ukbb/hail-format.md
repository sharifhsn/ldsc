# Hail Format

## Release files[#](https://pan.ukbb.broadinstitute.org/docs/hail-format#release-files "Direct link to heading")

The results of this analysis are released in two main files on Google Cloud Storage (file format compatible with Hail >= 0.2.42):

-   Summary statistics MatrixTable: `gs://ukb-diverse-pops-public/sumstats_release/results_full.mt` (12.78 T)
-   Meta-analysis MatrixTables (see the detailed descriptions [here](https://pan.ukbb.broadinstitute.org/docs/hail-format#meta-analysis-files))
    -   "High-quality" meta-analyses: `gs://ukb-diverse-pops-public/sumstats_release/meta_analysis.h2_qc.mt` (1.5 T)
    -   All ancestries (no QC) meta-analyses: `gs://ukb-diverse-pops-public/sumstats_release/meta_analysis.raw.mt` (12.5 T)

These are also available on Amazon S3:

-   Summary statistics MatrixTable: `s3://pan-ukb-us-east-1/sumstats_release/results_full.mt` (12.78 T)
-   Meta-analysis MatrixTables (see the detailed descriptions [here](https://pan.ukbb.broadinstitute.org/docs/hail-format#meta-analysis-files))
    -   "High-quality" meta-analyses: `s3://pan-ukb-us-east-1/sumstats_release/meta_analysis.h2_qc.mt` (1.5 T)
    -   All ancestries (no QC) meta-analyses: `s3://pan-ukb-us-east-1/sumstats_release/meta_analysis.raw.mt` (12.5 T)

In addition, in-sample full LD matrices and scores are available on Amazon S3:

-   LD BlockMatrix `s3://pan-ukb-us-east-1/ld_release/UKBB.{pop}.ldadj.bm` (43.3 T in total)
    -   Size by population: AFR: 12.0 T, AMR: 3.3 T, CSA: 6.4T, EAS: 2.6T, EUR: 14.1T, MID: 4.9T
-   Variant index Hail Table `s3://pan-ukb-us-east-1/ld_release/UKBB.{pop}.ldadj.variant.ht` (1.7 G in total)
-   LD score Hail Table `s3://pan-ukb-us-east-1/ld_release/UKBB.{pop}.ldscore.ht` (4.0 G in total)

where `{pop}` represents one of the population abbreviations (i.e., AFR, AMR, CSA, EAS, EUR, or MID).

### Requester pays[#](https://pan.ukbb.broadinstitute.org/docs/hail-format#requester-pays "Direct link to heading")

Note that the files in the Google Cloud Storage bucket are "requester pays." In order to compute over these files or download them, you will need to specify a project which may be billed for access and download costs. The data are stored in a US multi-region bucket: thus, access to the dataset is free for use for Compute Engine instances started within US regions, as well as for full downloads within the US and Canada. When performing large analyses on the dataset, we suggest "bringing the compute to the data" and starting a VM or Dataproc cluster in a US region. You can browse the directory structure in a requester pays bucket with the `-u` flag (and note the `hl.init` call below to access the data using Hail):

gsutil -u your\_project\_id ls gs://ukb-diverse-pops-public/sumstats\_release

Copy

## Using the libraries and files[#](https://pan.ukbb.broadinstitute.org/docs/hail-format#using-the-libraries-and-files "Direct link to heading")

##### note

**Please note that p-values are now stored as -log10 p-values to avoid underflow.**

The files on Google Cloud Platform can be accessed by cloning the [ukbb\_pan\_ancestry](https://github.com/atgu/ukbb_pan_ancestry) and the [ukb\_common](https://github.com/Nealelab/ukb_common) repos and accessing them programmatically. We recommend using these functions, as they allow for automatic application of our QC metrics as well as inclusion of all [QC flags](https://pan.ukbb.broadinstitute.org/docs/qc#quality-control-of-summary-statistics) and convenience metrics such as lambda GC. By default, when loading using `load_final_sumstats_mt`, the best practice QC parameters are used, which removes traits with a lambda GC < 0.5 or > 5 as well as applying all [QC filters](https://pan.ukbb.broadinstitute.org/docs/qc#quality-control-of-summary-statistics). This results in importing summary statistics for 527 traits; if it is preferable to load all traits with exported summary statistics (e.g., only applying the lambda GC < 0.5 or > 5 filter), use `load_final_sumstats_mt(filter_pheno_h2_qc=False)`, resulting in 7,228 traits. If any filtering is undesirable, use `load_final_sumstats_mt(filter_pheno_h2_qc=False, filter_phenos=False)`, which will import all 7,271 traits.

%%bash

git clone https://github.com/atgu/ukbb\_pan\_ancestry

git clone https://github.com/Nealelab/ukb\_common

Copy

from ukbb\_pan\_ancestry import \*

hl.init(spark\_conf={'spark.hadoop.fs.gs.requester.pays.mode': 'AUTO',

'spark.hadoop.fs.gs.requester.pays.project.id': 'your\_project\_id'})

\# loads all results for which sumstats were exported (lambda GC < 0.5 or > 5)

mt = load\_final\_sumstats\_mt(filter\_pheno\_h2\_qc=False)

\# use filter\_pheno\_h2\_qc=True to filter to just ancestry-trait pairs passing all QC

\# mt = load\_final\_sumstats\_mt(filter\_pheno\_h2\_qc=True)

mt.describe()

Copy

Additional arguments can be used to adjust the format of the p-values. Enabling `exponentiate_p` will produce absolute-scale p-values, while `legacy_exp_p_values` will produce legacy p-values that are in the format `ln p`.

## Results schema[#](https://pan.ukbb.broadinstitute.org/docs/hail-format#results-schema "Direct link to heading")

##### note

**Please note that p-values are now stored as -log10 p-values to avoid underflow.**

The basic summary statistics have the following schema:

\----------------------------------------

Column fields:

'trait\_type': str

'phenocode': str

'pheno\_sex': str

'coding': str

'modifier': str

'pheno\_data': struct {

n\_cases: int32,

n\_controls: int32,

heritability: struct {

estimates: struct {

ldsc: struct {

h2\_liability: float64,

h2\_liability\_se: float64,

h2\_z: float64,

h2\_observed: float64,

h2\_observed\_se: float64,

intercept: float64,

intercept\_se: float64,

ratio: float64,

ratio\_se: float64

},

sldsc\_25bin: struct {

h2\_liability: float64,

h2\_liability\_se: float64,

h2\_z: float64,

h2\_observed: float64,

h2\_observed\_se: float64,

intercept: float64,

intercept\_se: float64,

ratio: float64,

ratio\_se: float64

},

rhemc\_25bin: struct {

h2\_liability: float64,

h2\_liability\_se: float64,

h2\_z: float64,

h2\_observed: float64,

h2\_observed\_se: float64

},

rhemc\_8bin: struct {

h2\_liability: float64,

h2\_liability\_se: float64,

h2\_observed: float64,

h2\_observed\_se: float64,

h2\_z: float64

},

rhemc\_25bin\_50rv: struct {

h2\_observed: float64,

h2\_observed\_se: float64,

h2\_liability: float64,

h2\_liability\_se: float64,

h2\_z: float64

},

final: struct {

h2\_observed: float64,

h2\_observed\_se: float64,

h2\_liability: float64,

h2\_liability\_se: float64,

h2\_z: float64

}

},

qcflags: struct {

GWAS\_run: bool,

ancestry\_reasonable\_n: bool,

defined\_h2: bool,

significant\_z: bool,

in\_bounds\_h2: bool,

normal\_lambda: bool,

normal\_ratio: bool,

EUR\_plus\_1: bool,

pass\_all: bool

},

N\_ancestry\_QC\_pass: int32

},

saige\_version: str,

inv\_normalized: bool,

pop: str,

lambda\_gc: float64,

n\_variants: int64,

n\_sig\_variants: int64,

saige\_heritability: float64

}

'description': str

'description\_more': str

'coding\_description': str

'category': str

'n\_cases\_full\_cohort\_both\_sexes': int64

'n\_cases\_full\_cohort\_females': int64

'n\_cases\_full\_cohort\_males': int64

'pop\_index': int32

\----------------------------------------

Row fields:

'locus': locus<GRCh37>

'alleles': array<str>

'rsid': str

'varid': str

'vep': struct {

...

}

'freq': array<struct {

pop: str,

ac: float64,

af: float64,

an: int64,

gnomad\_exomes\_ac: int32,

gnomad\_exomes\_af: float64,

gnomad\_exomes\_an: int32,

gnomad\_genomes\_ac: int32,

gnomad\_genomes\_af: float64,

gnomad\_genomes\_an: int32

}>

'pass\_gnomad\_exomes': bool

'pass\_gnomad\_genomes': bool

'n\_passing\_populations': int32

'high\_quality': bool

'nearest\_genes': array<struct {

gene\_id: str,

gene\_name: str,

within\_gene: bool

}>

'info': float64

\----------------------------------------

Entry fields:

'summary\_stats': struct {

AF\_Allele2: float64,

imputationInfo: float64,

BETA: float64,

SE: float64,

\`p.value.NA\`: float64,

\`AF.Cases\`: float64,

\`AF.Controls\`: float64,

Pvalue: float64,

low\_confidence: bool

}

\----------------------------------------

Column key: \['trait\_type', 'phenocode', 'pheno\_sex', 'coding', 'modifier'\]

Row key: \['locus', 'alleles'\]

\----------------------------------------

Copy

### Columns (phenotypes)[#](https://pan.ukbb.broadinstitute.org/docs/hail-format#columns-phenotypes "Direct link to heading")

The columns are indexed by phenotype using a composite key of trait type, phenocode, pheno\_sex, coding, and modifier. Trait types have one of the values below. `phenocode` typically corresponds to the Field from UK Biobank, or the specific ICD code or phecode, or a custom moniker. `pheno_sex` designates which sexes were run, and is marked as `both_sexes` for most traits, though some phecodes were restricted to `females` or `males`. The `coding` field is primarily used for categorical variables, to indicate which one-hot encoding was used (e.g. [coding 2](http://biobank.ctsu.ox.ac.uk/showcase/coding.cgi?id=100434) for [field 1747](http://biobank.ctsu.ox.ac.uk/showcase/field.cgi?id=1747)). Finally, `modifier` refers to any downstream modifications of the phenotype (e.g. `irnt` for inverse-rank normal transformation).

By default, the MatrixTable loaded by `load_final_sumstats_mt` returns one column per phenotype-population pair. We can see the number of unique phenotypes for each `trait_type` by:

phenotype\_ht = mt.cols().collect\_by\_key() # Converting into one element per phenotype

phenotype\_ht.group\_by('trait\_type').aggregate(n\_phenos=hl.agg.count()).show()

\# results for all exported sumstats

+-----------------+----------+

| trait\_type | n\_phenos |

+-----------------+----------+

| str | int64 |

+-----------------+----------+

| "biomarkers" | 30 |

| "categorical" | 3686 |

| "continuous" | 820 |

| "icd10" | 921 |

| "phecode" | 1326 |

| "prescriptions" | 445 |

+-----------------+----------+

\# results for full QC pass-only sumstats

+-----------------+----------+

| trait\_type | n\_phenos |

+-----------------+----------+

| str | int64 |

+-----------------+----------+

| "biomarkers" | 23 |

| "categorical" | 179 |

| "continuous" | 206 |

| "icd10" | 34 |

| "phecode" | 64 |

| "prescriptions" | 21 |

+-----------------+----------+

Copy

You can explore the population-level data in more detail using (several fields removed for brevity):

phenotype\_ht = mt.cols()

phenotype\_ht.show(truncate=40, width=105)

+--------------+-----------+--------------+--------+----------+--------------------+

| trait\_type | phenocode | pheno\_sex | coding | modifier | pheno\_data.n\_cases |

+--------------+-----------+--------------+--------+----------+--------------------+

| str | str | str | str | str | int32 |

+--------------+-----------+--------------+--------+----------+--------------------+

| "biomarkers" | "30600" | "both\_sexes" | "" | "irnt" | 7694 |

| "biomarkers" | "30600" | "both\_sexes" | "" | "irnt" | 367192 |

| "biomarkers" | "30610" | "both\_sexes" | "" | "irnt" | 8422 |

| "biomarkers" | "30610" | "both\_sexes" | "" | "irnt" | 400988 |

| "biomarkers" | "30620" | "both\_sexes" | "" | "irnt" | 6214 |

| "biomarkers" | "30620" | "both\_sexes" | "" | "irnt" | 8407 |

| "biomarkers" | "30620" | "both\_sexes" | "" | "irnt" | 400822 |

| "biomarkers" | "30620" | "both\_sexes" | "" | "irnt" | 1499 |

| "biomarkers" | "30630" | "both\_sexes" | "" | "irnt" | 7679 |

| "biomarkers" | "30630" | "both\_sexes" | "" | "irnt" | 364987 |

+--------------+-----------+--------------+--------+----------+--------------------+

+-----------------------+------------------------------------------+

| pheno\_data.n\_controls | pheno\_data.heritability.estimates.lds... |

+-----------------------+------------------------------------------+

| int32 | float64 |

+-----------------------+------------------------------------------+

| NA | 1.62e-01 |

| NA | 1.18e-01 |

| NA | 1.98e-01 |

| NA | 2.18e-01 |

| NA | 1.27e-01 |

| NA | 1.48e-02 |

| NA | 1.14e-01 |

| NA | -2.67e-01 |

| NA | 1.30e-01 |

| NA | 1.89e-01 |

+-----------------------+------------------------------------------+

+------------------------------------------+------------------------------------------+

| pheno\_data.heritability.qcflags.norma... | pheno\_data.heritability.qcflags.norma... |

+------------------------------------------+------------------------------------------+

| bool | bool |

+------------------------------------------+------------------------------------------+

| True | True |

| True | True |

| True | True |

| True | True |

| True | True |

| True | True |

| True | True |

| True | True |

| True | True |

| True | True |

+------------------------------------------+------------------------------------------+

+------------------------------------------+------------------------------------------+

| pheno\_data.heritability.qcflags.EUR\_p... | pheno\_data.heritability.qcflags.pass\_all |

+------------------------------------------+------------------------------------------+

| bool | bool |

+------------------------------------------+------------------------------------------+

| True | True |

| True | True |

| True | True |

| True | True |

| True | True |

| True | True |

| True | True |

| True | True |

| True | True |

| True | True |

+------------------------------------------+------------------------------------------+

+------------------------------------------+--------------------------+---------------------------+

| pheno\_data.heritability.N\_ancestry\_QC... | pheno\_data.saige\_version | pheno\_data.inv\_normalized |

+------------------------------------------+--------------------------+---------------------------+

| int32 | str | bool |

+------------------------------------------+--------------------------+---------------------------+

| 2 | "SAIGE\_0.36.4" | False |

| 2 | "SAIGE\_0.36.4" | False |

| 2 | "SAIGE\_0.36.4" | False |

| 2 | "SAIGE\_0.44.5" | False |

| 4 | "SAIGE\_0.36.4" | False |

| 4 | "SAIGE\_0.36.4" | False |

| 4 | "SAIGE\_0.44.5" | False |

| 4 | "SAIGE\_0.36.4" | False |

| 3 | "SAIGE\_0.36.4" | False |

| 3 | "SAIGE\_0.44.5" | False |

+------------------------------------------+--------------------------+---------------------------+

+----------------+----------------------+-----------------------+---------------------------+

| pheno\_data.pop | pheno\_data.lambda\_gc | pheno\_data.n\_variants | pheno\_data.n\_sig\_variants |

+----------------+----------------------+-----------------------+---------------------------+

| str | float64 | int64 | int64 |

+----------------+----------------------+-----------------------+---------------------------+

| "CSA" | 1.03e+00 | 12200078 | 6 |

| "EUR" | 1.37e+00 | 20561726 | 37450 |

| "CSA" | 1.02e+00 | 12364741 | 772 |

| "EUR" | 1.67e+00 | 20739978 | 89683 |

| "AFR" | 1.02e+00 | 18630599 | 77 |

| "CSA" | 1.02e+00 | 12362444 | 0 |

| "EUR" | 1.42e+00 | 20739238 | 38220 |

| "MID" | 9.89e-01 | 12328418 | 0 |

| "CSA" | 1.01e+00 | 12195348 | 378 |

| "EUR" | 1.63e+00 | 20547047 | 62484 |

+----------------+----------------------+-----------------------+---------------------------+

+-------------------------------+----------------------------+------------------+--------------------+

| pheno\_data.saige\_heritability | description | description\_more | coding\_description |

+-------------------------------+----------------------------+------------------+--------------------+

| float64 | str | str | str |

+-------------------------------+----------------------------+------------------+--------------------+

| 2.41e-01 | "Albumin" | NA | NA |

| 6.45e-02 | "Albumin" | NA | NA |

| 3.92e-01 | "Alkaline phosphatase" | NA | NA |

| 1.31e-01 | "Alkaline phosphatase" | NA | NA |

| 3.94e-01 | "Alanine aminotransferase" | NA | NA |

| 2.19e-01 | "Alanine aminotransferase" | NA | NA |

| 6.28e-02 | "Alanine aminotransferase" | NA | NA |

| 2.05e-01 | "Alanine aminotransferase" | NA | NA |

| 3.58e-01 | "Apolipoprotein A" | NA | NA |

| 1.22e-01 | "Apolipoprotein A" | NA | NA |

+-------------------------------+----------------------------+------------------+--------------------+

+------------------------------------------+--------------------------------+

| category | n\_cases\_full\_cohort\_both\_sexes |

+------------------------------------------+--------------------------------+

| str | int64 |

+------------------------------------------+--------------------------------+

| "Biological samples > Assay results >... | 422605 |

| "Biological samples > Assay results >... | 422605 |

| "Biological samples > Assay results >... | 461525 |

| "Biological samples > Assay results >... | 461525 |

| "Biological samples > Assay results >... | 461326 |

| "Biological samples > Assay results >... | 461326 |

| "Biological samples > Assay results >... | 461326 |

| "Biological samples > Assay results >... | 461326 |

| "Biological samples > Assay results >... | 420088 |

| "Biological samples > Assay results >... | 420088 |

+------------------------------------------+--------------------------------+

+-----------------------------+---------------------------+-----------+

| n\_cases\_full\_cohort\_females | n\_cases\_full\_cohort\_males | pop\_index |

+-----------------------------+---------------------------+-----------+

| int64 | int64 | int32 |

+-----------------------------+---------------------------+-----------+

| 208336 | 179998 | 0 |

| 208336 | 179998 | 1 |

| 229156 | 194910 | 0 |

| 229156 | 194910 | 1 |

| 229118 | 194764 | 0 |

| 229118 | 194764 | 1 |

| 229118 | 194764 | 2 |

| 229118 | 194764 | 3 |

| 206413 | 179623 | 0 |

| 206413 | 179623 | 1 |

+-----------------------------+---------------------------+-----------+

showing top 10 rows

Copy

More information about the GWAS run is found in the `pheno_data` struct. This struct also includes all heritability information and QC flags. More details on heritability estimation methods are forthcoming, but can be previewed [here](https://pan.ukbb.broadinstitute.org/blog/2022/04/11/h2-qc-updated-sumstats). Descriptions of the heritability fields can be found [below](https://pan.ukbb.broadinstitute.org/docs/hail-format#heritability) and more information on QC flags can be found [here.](https://pan.ukbb.broadinstitute.org/docs/qc#quality-control-of-summary-statistics)

### Rows (variants)[#](https://pan.ukbb.broadinstitute.org/docs/hail-format#rows-variants "Direct link to heading")

The rows are indexed by locus and alleles. Direct annotations can be found in the `vep` schema, but we also provide a `nearest_genes` annotation for ease of analysis. Additionally, variant QC annotations are provided in the `high_quality` field (which is filtered to by default using `load_final_sumstats_mt` and can be switched off in the `filter_variants` parameter in that function).

### Entries (association tests)[#](https://pan.ukbb.broadinstitute.org/docs/hail-format#entries-association-tests "Direct link to heading")

##### note

Please note that p-values are now stored as log p-values to avoid underflow.

The entry fields house the summary statistics themselves. Note that there is a `low_confidence` annotation that indicates a possible low-quality association test (allele count in cases or controls <= 3, or overall minor allele count < 20).

The resulting dataset can be filtered and annotated as a standard Hail MatrixTable:

mt = mt.filter\_cols((mt.trait\_type == 'phecode') & (mt.lambda\_gc > 0.9) & (mt.lambda\_gc < 1.1))

Copy

## Meta-analysis files[#](https://pan.ukbb.broadinstitute.org/docs/hail-format#meta-analysis-files "Direct link to heading")

##### note

Please note that p-values are now stored as -log10 p-values to avoid underflow.

The meta-analysis results are in a similarly structured file which can be obtained as such:

meta\_mt = load\_meta\_analysis\_results()

Copy

By default, this function imports both the meta-analyses for all phenotype-ancestry pairs with completed GWAS as well as the "high-quality" meta-analyses, which are meta-analyses of only ancestry groups passing all [QC filters](https://pan.ukbb.broadinstitute.org/docs/qc#quality-control-of-summary-statistics). Naturally the high-quality meta-analyses are only available for those phenotypes for which at least two ancestries passed all QC (since a requirement for QC PASS is that a trait passes QC in EUR and at least one other ancestry). The meta-analyses for all pairs are found under the `meta_analysis` struct, while the high-quality meta-analyses are found in the `meta_analysis_hq` struct.

Entry fields:

'meta\_analysis': array<struct {

BETA: float64,

SE: float64,

Pvalue: float64,

Q: float64,

Pvalue\_het: float64,

N: int32,

N\_pops: int32,

AF\_Allele2: float64,

AF\_Cases: float64,

AF\_Controls: float64

}>

'meta\_analysis\_hq': array<struct {

BETA: float64,

SE: float64,

Pvalue: float64,

Q: float64,

Pvalue\_het: float64,

N: int32,

N\_pops: int32,

AF\_Allele2: float64,

AF\_Cases: float64,

AF\_Controls: float64

}>

Copy

Meta-analyses using only the high-quality pairs or using all pairs can be imported seperately by using `load_meta_analysis_results(h2_filter='pass')` or `load_meta_analysis_results(h2_filter='none')` respectively. These arguments always produce tables with the following entry schema:

Entry fields:

'meta\_analysis': array<struct {

BETA: float64,

SE: float64,

Pvalue: float64,

Q: float64,

Pvalue\_het: float64,

N: int32,

N\_pops: int32,

AF\_Allele2: float64,

AF\_Cases: float64,

AF\_Controls: float64

}>

Copy

All meta-analysis tables have results provided in an array, which includes the all-available-population meta-analysis in the 0th element `meta_mt.meta_analysis[0]` and leave-one-out meta-analyses in the remainder of the array.

As for the per-population summary statistics, additional arguments can be used to adjust the format of the p-values. Enabling `exponentiate_p` will produce absolute-scale p-values, while `legacy_exp_p_values` will produce legacy p-values that are in the format `ln p`.

### Combining the datasets[#](https://pan.ukbb.broadinstitute.org/docs/hail-format#combining-the-datasets "Direct link to heading")

We also provide a function to annotate the overall sumstats MatrixTable with the largest meta-analysis for that phenotype. Specification of `h2_filter` will determine if the annotated meta analysis is for all included phenotypes or for just those passing QC.

mt = load\_final\_sumstats\_mt()

mt = annotate\_mt\_with\_largest\_meta\_analysis(mt, h2\_filter='none')

Copy

If your analysis requires the simultaneous analysis of summary statistics from multiple populations (and not the meta-analysis), you can load the data with a similar structure to the meta-analysis MatrixTable (one column per phenotype, with population information packed into an array of entries and columns) using `load_final_sumstats_mt(separate_columns_by_pop=False)`.

## LD matrices[#](https://pan.ukbb.broadinstitute.org/docs/hail-format#ld-matrices "Direct link to heading")

The LD matrices are in [`BlockMatrix`](https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html) format. Please refer to [Hail's documentation](https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html) for available operations on `BlockMatrix`.

from hail.linalg import BlockMatrix

bm = BlockMatrix.read(get\_ld\_matrix\_path(pop='AFR'))

Copy

We note that the LD matrices were sparsified to a upper triangle (all elements of the lower triangle were zeroed out using [`BlockMatrix.sparsify_triangle`](https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html#hail.linalg.BlockMatrix.sparsify_triangle)).

### Variant indices[#](https://pan.ukbb.broadinstitute.org/docs/hail-format#variant-indices "Direct link to heading")

To determine which row/column corresponds to which variant, we provide variant indices for `BlockMatrix` in Hail `Table` format.

ht\_idx = hl.read\_table(get\_ld\_variant\_index\_path(pop='AFR'))

Copy

The variant indices table has the following schema and `idx` corresponds to a row/column index in `BlockMatrix`.

\----------------------------------------

Global fields:

'n\_samples': int32

'pop': str

\----------------------------------------

Row fields:

'locus': locus<GRCh37>

'alleles': array<str>

'idx': int64

\----------------------------------------

Key: \['locus', 'alleles'\]

\----------------------------------------

Copy

### Extracting a subset of LD matrix[#](https://pan.ukbb.broadinstitute.org/docs/hail-format#extracting-a-subset-of-ld-matrix "Direct link to heading")

To extract a subset of LD matrix, you first need to identify indices of your variants of interest. Here, we provide two examples:

\# filter by interval

interval = hl.parse\_locus\_interval('1:51572000-52857000')

ht\_idx = ht\_idx.filter(interval.contains(ht\_idx.locus))

\# or filter by a list of variant IDs (e.g., 1:51572412:A:G)

ht = hl.import\_table('/path/to/your/list')

ht = ht.transmute(\*\*hl.parse\_variant(ht.variant)).key\_by('locus', 'alleles')

ht\_idx = ht\_idx.join(ht, 'inner')

Copy

Then, you can filter the LD matrix into a subset using [`BlockMatrix.filter`](https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html#hail.linalg.BlockMatrix.filter):

idx = ht\_idx.idx.collect()

bm = bm.filter(idx, idx)

Copy

### Exporting a LD matrix to a flat file[#](https://pan.ukbb.broadinstitute.org/docs/hail-format#exporting-a-ld-matrix-to-a-flat-file "Direct link to heading")

Finally, to export a LD matrix to a flat file (txt file), you can use [`BlockMatrix.export`](https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html#hail.linalg.BlockMatrix.export):

\# Note: when you apply any operation on BlockMatrix,

\# you need to write it first to storage before export

bm.write('/path/to/tmp/bm', force\_row\_major=True)

BlockMatrix.export(

'/path/to/tmp/bm',

'/path/to/flat\_file.bgz',

delimiter=' '

)

Copy

If your matrix is small enough to fit on memory, you can also directly export it to numpy via [`BlockMatrix.to_numpy`](https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html#hail.linalg.BlockMatrix.to_numpy).

np\_mat = bm.to\_numpy()

Copy

## LD scores[#](https://pan.ukbb.broadinstitute.org/docs/hail-format#ld-scores "Direct link to heading")

The LD scores are in Hail [Table](https://hail.is/docs/0.2/hail.Table.html) format. For LDSC-compatible flat files, you can find them [here](https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.ALL.ldscore.tar.gz).

ht = hl.read\_table(get\_ld\_score\_ht\_path(pop='AFR'))

Copy

The LD score table has the following schema.

\----------------------------------------

Global fields:

None

\----------------------------------------

Row fields:

'locus': locus<GRCh37>

'alleles': array<str>

'rsid': str

'varid': str

'AF': float64

'ld\_score': float64

\----------------------------------------

Key: \['locus', 'alleles'\]

\----------------------------------------

Copy

## Heritability estimates[#](https://pan.ukbb.broadinstitute.org/docs/hail-format#heritability-estimates "Direct link to heading")

The heritability estimates can be found as a [flat file manifest](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1797288938) or as a Hail [Table](https://hail.is/docs/0.2/hail.Table.html).

ht = hl.read\_table(get\_h2\_ht())

Copy

This table has the following schema:

\----------------------------------------

Global fields:

None

\----------------------------------------

Row fields:

'trait\_type': str

'phenocode': str

'pheno\_sex': str

'coding': str

'modifier': str

'heritability': array<struct {

pop: str,

estimates: struct {

ldsc: struct {

h2\_liability: float64,

h2\_liability\_se: float64,

h2\_z: float64,

h2\_observed: float64,

h2\_observed\_se: float64,

intercept: float64,

intercept\_se: float64,

ratio: float64,

ratio\_se: float64

},

sldsc\_25bin: struct {

h2\_liability: float64,

h2\_liability\_se: float64,

h2\_z: float64,

h2\_observed: float64,

h2\_observed\_se: float64,

intercept: float64,

intercept\_se: float64,

ratio: float64,

ratio\_se: float64

},

rhemc\_25bin: struct {

h2\_liability: float64,

h2\_liability\_se: float64,

h2\_z: float64,

h2\_observed: float64,

h2\_observed\_se: float64

},

rhemc\_8bin: struct {

h2\_liability: float64,

h2\_liability\_se: float64,

h2\_observed: float64,

h2\_observed\_se: float64,

h2\_z: float64

},

rhemc\_25bin\_50rv: struct {

h2\_observed: float64,

h2\_observed\_se: float64,

h2\_liability: float64,

h2\_liability\_se: float64,

h2\_z: float64

},

final: struct {

h2\_observed: float64,

h2\_observed\_se: float64,

h2\_liability: float64,

h2\_liability\_se: float64,

h2\_z: float64

}

},

qcflags: struct {

GWAS\_run: bool,

ancestry\_reasonable\_n: bool,

defined\_h2: bool,

significant\_z: bool,

in\_bounds\_h2: bool,

normal\_lambda: bool,

normal\_ratio: bool,

EUR\_plus\_1: bool,

pass\_all: bool

},

N\_ancestry\_QC\_pass: int32

}>

\----------------------------------------

Key: \['trait\_type', 'phenocode', 'pheno\_sex', 'coding', 'modifier'\]

\----------------------------------------

Copy

Note that this is very similar to the heritability struct in the [results schema](https://pan.ukbb.broadinstitute.org/docs/hail-format#results-schema) -- the `load_final_sumstats_mt()` automatically uses `get_h2_ht()` to import the heritability table and annotate it into the column schema of the summary statistics results table.

Here the heritability field is an `array` of structs, where the elements of the array represent the ancestries for which heritability estimates are available. The corresponding ancestry is found in the `heritability.pop` field. To obtain one row for each ancestry-trait pair, use the following commmand:

ht = ht.explode('heritability')

Copy

The `heritability.estimates` struct contains point estimates and significance test results for:

-   Univariate LD score regression (`ldsc`), run using [LD score flat files](https://pan.ukbb.broadinstitute.org/docs/ld#ld-scores) using high-quality HapMap3 SNPs with MAF ≥\\geq≥ 0.01 with summary statistics exported from the [results table](https://pan.ukbb.broadinstitute.org/docs/hail-format#results-schema).
-   Stratified LD score regression (`sldsc_25bin`), run using the same summary statistcs as `ldsc` with LD scores generated from SNPs in 5 MAF and 5 LD score bins.
-   Randomized Haseman-Elston (`rhemc_8bin`) using genotype data with 2 MAF and 4 LD score bins using default settings. This was predominantly used to analyze non-EUR ancestry groups but includes a small set of estimates for traits in EUR.
-   Randomized Haseman-Elston (`rhemc_25bin`) using genotype data with 5 MAF and 5 LD score bins using default settings. Only run for non-EUR ancestry groups.
-   Randomized Haseman-Elston (`rhemc_25bin_50rv`) using genotype data with 5 MAF and 5 LD score bins using 50 random variables to reduce run-to-run variability at a slightly higher computational cost. Only run for non-EUR ancestry groups.

Within each of the above methods, we produce observed- and liability-scale estimates as well as standard errors and the z-score for the test of h2\>0h^2 > 0h2\>0. We also produce LDSC intercept and ratio estimates when relevant.

The `heritability.qcflags` struct contains the results from our sequential QC filtering scheme -- see [quality control](https://pan.ukbb.broadinstitute.org/docs/qc#quality-control-of-summary-statistics) for more details. `N_ancestry_QC_pass` refers to the number of ancestries passing all QC for a given trait.

More information can be found [here](https://pan.ukbb.broadinstitute.org/docs/heritability) on the heritability estimation approach, and important caveats can be found [here](https://pan.ukbb.broadinstitute.org/docs/qc#heritability).