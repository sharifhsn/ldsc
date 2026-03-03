This README is a description of the phenotype manifest corresponding to the Pan-UK Biobank Project																					
For a description of the project and details of the analysis, please see https://pan.ukbb.broadinstitute.org/																					
The code used to generate the GWAS summary statistics is publicly available here: https://github.com/atgu/ukbb_pan_ancestry																					
For information about summary statistics files, see: https://pan.ukbb.broadinstitute.org/docs/per-phenotype-files#per-phenotype-files																					
Note: this manifest was updated on 23-03-01. The old manifest is archived at OLD_20220411_phenotype_manifest.																					
Please note: the p-values in the summary statistic flat files are negative log10 p-values. P-value columns have been renamed with a "neglog10" suffix to indicate this.																					
To cite the Pan-UKB data:	Karczewski, Gupta, Kanai et al., 2024																				
																					
Notes about fields in the manifest:																					
Phenotype ID: The first 5 fields are guaranteed to be unique																					
Cohort case counts: If a trait is quantitative (trait_type is "continuous" or "biomarkers"), all samples are considered to be "cases". Thus, the number of cases is equivalent to the number of samples.																					
Cohort case counts: A  meta analysis was performed using ancestries for which GWAS we performed per trait -- see pops. The corresponding case counts are in n_cases_full_cohort*. All meta analyses use only non-low-confidence SNPs.																					
Cohort case counts: A seperate "high quality" meta analysis was performed using only ancestries passing QC per trait -- see pops_pass_qc. The corresponding case counts for these passing ancestries per-trait are in n_cases_hq_cohort*.																					
Population-specific: The variable pop is a placeholder for a 3-letter ancestry code. For example, n_cases_AFR is the number of cases with AFR ancestry.																					
Population-specific: If a trait is quantitative (trait_type is "continuous" or "biomarkers"), all samples are considered to be "cases". Thus, the number of cases is equivalent to the number of samples.																					
Population-specific: For heritability fields (e.g. h2_observed, h2_liability, h2_z), we used a genotype-based Haseman-Elston regression approach (RHEmc) to maxmize power for non-EUR ancestries, while using stratified LD-score regression (SLDSC) for EUR ancestry. This is indicated by the field name. See h2 manifest for more details.																					
Population-specific: For heritability fields (e.g. h2_observed, h2_liability, h2_z), the corresponding results can all be found in the estimates.final.* columns in the h2 manifest, or in the estimates.rhemc_25bin_50rv.* and estimates.sldsc_25bin.* columns for non-EUR and EUR respectively.																					
File information: For each field in this section there also exists a field with the suffix _tabix, which contains the equivalent information for the tabix file. For instance, filename_tabix contains the name of the tabix file.																					
																					
Field type	Field	Descriptor																			
Phenotype ID	trait_type	One of the following: continuous, biomarkers, prescriptions, icd10, phecode, categorical																			
Phenotype ID	phenocode	The code for the phenotype (for continuous, biomarkers, and categorical traits, this corresponds to the field ID as described by UKB, e.g. 21001 for BMI)																			
Phenotype ID	pheno_sex	Indicating whether the phenotype was run for both sexes (pheno_sex="both_sexes") or in just females (pheno_sex="females") or males (pheno_sex="males"). In 0.1, this is only differentiated for phecodes																			
Phenotype ID	coding	For categorical variables, this corresponds to the coding that was used (e.g. coding 2 for field 1747). For all other trait_types, this field is blank																			
Phenotype ID	modifier	Refers to any miscellaneous downstream modifications of the phenotype (e.g. irnt for inverse-rank normal transformation). If the phenotype is updated, this field can be used to denote the update (e.g. the particular wave of COVID-19 data used).																			
Phenotype ID	description	A shorter description of the phenotype (for continuous, biomarkers, and categorical variables, corresponds to the Description on the showcase). For phecodes, this is the "description" column in the phecodes definition file.																			
Phenotype ID	description_more	A longer description of the phenotype (for continuous and categorical variables, corresponds to the Notes page on the showcase).																			
Phenotype ID	coding_description	For categorical variables, a description of the particular coding that was used (the Meaning column on the showcase page for that coding).																			
Phenotype ID	category	A categorization of the phenotype. For continuous, biomarkers, and categorical traits, this corresponds to the Category at the top of the showcase page. For ICD codes, this corresponds to the Chapter of the ICD code; for phecodes, this is the "group" column in the phecodes definition file; for prescriptions, this corresponds to a semi-manual categorization of prescription drugs.																			
Phenotype ID	in_max_independent_set	If the phenotype is in the maximially independent set, a set of traits passing all QC and showing low correlations with one another.																			
Cohort case counts	n_cases_full_cohort_both_sexes	Number of cases (or individuals phenotyped for quantitative traits) across all ancestry groups, females and males combined. Should be similar to the sum of per-ancestry n_cases, but may include some ancestry outliers and samples that failed QC.																			
Cohort case counts	n_cases_full_cohort_females	Number of female cases (or individuals phenotyped for quantitative traits) across all ancestry groups. May include some ancestry outliers and samples that failed QC.																			
Cohort case counts	n_cases_full_cohort_males	Number of male cases (or individuals phenotyped for quantitative traits) across all ancestry groups. May include some ancestry outliers and samples that failed QC.																			
Cohort case counts	n_cases_hq_cohort_both_sexes	Number of cases (or individuals phenotyped for quantitative traits) across ancestry groups passing stringent phenotype QC (see pops_pass_qc), females and males combined. Should be similar to the sum of per-ancestry n_cases for relevant ancestries, but may include some ancestry outliers and samples that failed QC.																			
Cohort case counts	n_cases_hq_cohort_females	Number of female cases (or individuals phenotyped for quantitative traits) across ancestry groups passing stringent phenotype QC. May include some ancestry outliers and samples that failed QC.																			
Cohort case counts	n_cases_hq_cohort_males	Number of male cases (or individuals phenotyped for quantitative traits) across ancestry groups passing stringent phenotype QC. May include some ancestry outliers and samples that failed QC.																			
Ancestry metadata	pops	Comma-delimited list of ancestry codes for which this phenotypes was GWASed.																			
Ancestry metadata	num_pops	Number of ancestry groups for which this phenotype was GWASed.																			
Ancestry metadata	pops_pass_qc	Comma-delimited list of ancestry codes for which this phenotype passes QC (see h2 manifest and phenotype_qc_{pop} field).																			
Ancestry metadata	num_pops_pass_qc	Number of ancestry groups for which this phenotype passes QC.																			
Population-specific	n_cases_{pop}	Number of cases (or individuals phenotyped for quantitative traits) in the GWAS analysis. Excludes ancestry outliers and samples that failed QC.																			
Population-specific	n_controls_{pop}	Number of controls in the GWAS analysis. Excludes ancestry outliers and samples that failed QC.																			
Population-specific	{rhemc_25bin_50rv/sldsc_25bin}_h2_observed_{pop}	Observed scale heritability estimates using 25 MAF/LD bins with RHEmc (non-EUR) or SLDSC (EUR).																			
Population-specific	{rhemc_25bin_50rv/sldsc_25bin}_h2_observed_se_{pop}	Observed scale heritability standard error using 25 MAF/LD bins with RHEmc (non-EUR) or SLDSC (EUR).																			
Population-specific	{rhemc_25bin_50rv/sldsc_25bin}_h2_liability_{pop}	Libaility scale heritability estimates using 25 MAF/LD bins with RHEmc (non-EUR) or SLDSC (EUR), transformed using per-ancestry in-sample prevalence.																			
Population-specific	{rhemc_25bin_50rv/sldsc_25bin}_h2_liability_se_{pop}	Liability scale heritability standard error using 25 MAF/LD bins with RHEmc (non-EUR) or SLDSC (EUR), transformed using per-ancestry in-sample prevalence.																			
Population-specific	{rhemc_25bin_50rv/sldsc_25bin}_h2_z_{pop}	Heritability Z-scores (for test of h2 > 0) using per-ancestry-trait pair h2 estimates and standard errors.																			
Population-specific	lambda_gc_{pop}	The genomic control (lambda GC) calculated from the summary statistics for pop with low-confidence statistics removed and only considering high-quality variants.																			
Population-specific	phenotype_qc_{pop}	Phenotype QC outcome for each ancestry-trait pair. Filters are described in the h2 manifest in more detail.  Filters are applied sequentially; this field specifies either PASS or the reason for failure.																			
File information	filename	Name of summary statistics file.																			
File information	filename_tabix	Name of summary statistics tabix index.																			
File information	aws_link	(online manifest only) Link to download summary statistics file from Amazon AWS.																			
File information	aws_link_tabix	(online manifest only) Link to download summary statistics file tabix index from Amazon AWS.																			
File information	aws_path	(flat file manifest only) AWS s3 path to summary statistics file.																			
File information	aws_path_tabix	(flat file manifest only) AWS s3 path to summary statistics file tabix index.																			
File information	wget	wget command to download summary statistics file.																			
File information	wget_tabix	wget command to download summary statistics tabix file.																			
File information	size_in_bytes	Size of summary statistics file in bytes.																			
File information	size_in_bytes_tabix	Size of summary statistics file in bytes for tabix file.																			
File information	md5_hex	MD5 hexadecimal hash.																			
File information	md5_hex_tabix	MD5 hexadecimal hash for tabix file.																			
																					
Notes:																					
ICD codes were defined as any of UK Biobank fields 41202, 41204, 41201, 40001																					