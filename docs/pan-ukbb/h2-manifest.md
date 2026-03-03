This README is a description of the heritability manifest corresponding to the Pan-UK Biobank Project															
For a description of the project and details of the analysis, please see https://pan.ukbb.broadinstitute.org/															
The code used to generate the GWAS summary statistics is publicly available here: https://github.com/atgu/ukbb_pan_ancestry															
For information about summary statistics files, see: https://pan.ukbb.broadinstitute.org/docs/per-phenotype-files#per-phenotype-files															
Note: this manifest was created on 23-03-01.															
															
Notes about fields in the manifest:															
heritability is an array of structs which holds one element per heritability.pop.															
heritability.estimates contain all heritability point estimates.															
heritability.qcflags contains QC data as described below.															
Heritability methods: this field type refers to fields within the heritability.estimates struct which organize estimates from each tested heritability method.															
Heritability estimates: this field type refers to fields present in some (or all) of each heritability method struct. The intercept and ratio are only reported for LD score regression, while h2_liability, h2_observed, h2_z are present for all methods except Saige.															
Heritability estimates: saige only reports a pseudoheritability estimate which is reported in heritbability.estimates.saige.															
															
Field type	Field	Descriptor													
Phenotype ID	trait_type	One of the following: continuous, biomarkers, prescriptions, icd10, phecode, categorical													
Phenotype ID	phenocode	The code for the phenotype (for continuous, biomarkers, and categorical traits, this corresponds to the field ID as described by UKB, e.g. 21001 for BMI)													
Phenotype ID	pheno_sex	Indicating whether the phenotype was run for both sexes (pheno_sex="both_sexes") or in just females (pheno_sex="females") or males (pheno_sex="males"). In 0.1, this is only differentiated for phecodes													
Phenotype ID	coding	For categorical variables, this corresponds to the coding that was used (e.g. coding 2 for field 1747). For all other trait_types, this field is blank													
Phenotype ID	modifier	Refers to any miscellaneous downstream modifications of the phenotype (e.g. irnt for inverse-rank normal transformation). If the phenotype is updated, this field can be used to denote the update (e.g. the particular wave of COVID-19 data used).													
Phenotype ID	heritability.pop	Ancestry													
Heritability methods	heritability.estimates.ldsc.*	Univariate LD score regression													
Heritability methods	heritability.estimates.sldsc_25bin.*	Stratified LD score regression, 5 LD score bins x 5 MAF bins													
Heritability methods	heritability.estimates.rhemc_25bin.*	RHEmc (HE regression), 5 LD score bins x 5 MAF bins													
Heritability methods	heritability.estimates.rhemc_8bin.*	RHEmc (HE regression), 4 LD score bins x 2 MAF bins													
Heritability methods	heritability.estimates.rhemc_25bin_50rv.*	RHEmc (HE regression), 5 LD score bins x 5 MAF bins; 50 random variables for improved power													
Heritability methods	heritability.estimates.saige	Saige pseudo-h2 estimate													
Heritability methods	heritability.estimates.final.*	Final estimates; 25 bin SLDSC for EUR and 25 bin, 50 RV RHEmc for non-EUR (these are also present in the full manifest)													
Heritability estimates	heritability.estimates.*.h2_liability	Liability scale results, using in-sample prevalence													
Heritability estimates	heritability.estimates.*.h2_observed	Observed scale results													
Heritability estimates	heritability.estimates.*.h2_z	Heritability Z-scores (for test of h2 > 0)													
Heritability estimates	heritability.estimates.*.intercept	LDSC intercept													
Heritability estimates	heritability.estimates.*.ratio	LDSC ratio													
QC flags	heritability.N_ancestry _QC_pass	Number of ancestries passing QC per trait													
QC flags	heritability.qcflags.GWAS_run	Whether or not GWAS was run													
QC flags	heritability.qcflags.ancestry_reasonable_n	If the ancestry has reasonable sample size; removes AMR													
QC flags	heritability.qcflags.defined_h2	Whether h2 was non-missing													
QC flags	heritability.qcflags.significant_z	If h2 Z-score > 4													
QC flags	heritability.qcflags.in_bounds_h2	If h2 was within 0 and 1													
QC flags	heritability.qcflags.normal_lambda	If lambda was > 0.9													
QC flags	heritability.qcflags.normal_ratio	If ratio was < 0.3 or Z < 4 in EUR/CSA/AFR													
QC flags	heritability.qcflags.EUR_plus_1	If phenotype passed QC in EUR + 1 other ancestry													
QC flags	heritability.qcflags.pass_all	If QC was passed for this ancestry-trait pair (all of the above)													