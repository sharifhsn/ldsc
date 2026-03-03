# Pan-UKBB Docs Summary (Local Snapshot)

This summary distills the key items from `docs/pan-ukbb/*.md` that matter for LDSC workflows and performance testing.

## Data Releases and Where They Live
- Per-phenotype flat files: one `tsv.bgz` per phenotype, indexed, for single-phenotype analysis. The phenotype manifest lists file paths.
- Hail format: full MatrixTable and meta-analysis datasets on GCS/S3 for large-scale analysis.
- LD scores: LDSC-compatible flat files (`.l2.ldscore.gz` and `.M_5_50`) available as `UKBB.ALL.ldscore.tar.gz` on AWS.
- LD matrices: Hail BlockMatrix format on AWS for large-scale LD analysis.
- Heritability results: available as flat files and in Hail; phenotype and heritability manifests enumerate results.
- Phenotype correlation matrix: available as a flat file on AWS.

## Per-Phenotype File Format (Important for LDSC)
- Format: `tsv.bgz` (bgzip). Use `zcat`/`bgzip` or pandas/gzip to read.
- Variant fields: `chr`, `pos`, `ref`, `alt` (GRCh37). `alt` is the effect allele.
- Meta-analysis fields include `beta_meta`, `se_meta`, and p-value fields.
- Population-specific fields also include per-ancestry effect size and p-values.

## P-Value Encoding (Critical)
The docs contain a discrepancy:
- `downloads.md` says p-values are stored as **natural log p-values** (`ln P`).
- `per-phenotype-files.md` says p-values are stored as **-log10 p-values**, and lists fields like `neglog10_pval_meta` and `neglog10_pval_{pop}`.

Action: always **inspect the actual column names** in your downloaded file to decide how to transform to raw P for LDSC.

## LD Scores and Matrices (for LDSC)
- LD scores were computed in-sample per ancestry group.
- LD matrix window: 10 Mb; LD score window: 1 Mb.
- Variant filters for LDSC flat files: high-quality HapMap3 variants, autosomes, not in MHC, biallelic SNPs, INFO > 0.9, MAF > 1% in UKB and gnomAD.
- Covariate-adjusted LD scores (age, sex, age*sex, age^2, age^2*sex, first 10 PCs).

## Quality Control Highlights
- Variant QC uses INFO > 0.8 and MAC > 20 per population in LD computation.
- Variants are compared to gnomAD frequencies; discrepant variants are flagged in QC metrics.
- QC flags exist in the phenotype/heritability manifests to filter unreliable traits or ancestry-trait pairs.

## Practical Consequences for LDSC
- Use `UKBB.<POP>.rsid.l2.ldscore.gz` if your sumstats use rsids; otherwise use the non-rsid LD scores with `chr:pos:ref:alt` IDs.
- For LDMS partitioned scores, use `UKBB.<POP>.8LDMS.*` or `25LDMS.*` as `--ref-ld` and keep `--w-ld` on the base (single-column) LD score.
- Convert the log p-values from Pan-UKBB to raw P before running LDSC.

