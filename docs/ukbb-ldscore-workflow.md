# UKBB LD Score Bundle (data/UKBB.ALL.ldscore) — What It Is and How to Use It

## What This Data Is
The archive `data/UKBB.ALL.ldscore.tar.gz` unpacks to `data/UKBB.ALL.ldscore/` and contains **precomputed LD scores** for UK Biobank reference panels across multiple populations:

- Populations: `AFR`, `AMR`, `CSA`, `EAS`, `EUR`, `MID`
- Variants keyed either by:
  - `chr:pos:ref:alt` in `UKBB.<POP>.l2.ldscore.gz`, or
  - `rsid` in `UKBB.<POP>.rsid.l2.ldscore.gz`
- Partitioned LD scores by LD/MAF bins:
  - `UKBB.<POP>.8LDMS.l2.ldscore.gz`
  - `UKBB.<POP>.25LDMS.l2.ldscore.gz`
- Sidecar files for partitioned LD scores:
  - `.l2.M` and `.l2.M_5_50` (annotation SNP counts)
  - `.l2.annot.gz` (annotation matrix)
  - `rsid`-keyed variants of the above

These LD scores are **inputs** to `h2`/`rg` (and stratified LD score regression), not raw genotypes.

## Choosing the Right Files
Pick **population** and **variant ID scheme** to match your sumstats:

- If sumstats SNP IDs are `rsid`, use `UKBB.<POP>.rsid.l2.ldscore.gz` (and rsid variants of LDMS files).
- If sumstats SNP IDs are `chr:pos:ref:alt`, use `UKBB.<POP>.l2.ldscore.gz`.

## Recommended Workflows

### 1) Standard h2 / rg (single-column LD scores)
Use the base LD score file as both reference and weights unless you have a separate weights panel.

Python 3 (ldsc_py3):
```
uv run --project /Users/sharif/Code/ldsc/ldsc_py3 \
  python /Users/sharif/Code/ldsc/ldsc_py3/ldsc.py \
  --h2 /path/to/trait.sumstats.gz \
  --ref-ld /Users/sharif/Code/ldsc/data/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.ldscore.gz \
  --w-ld /Users/sharif/Code/ldsc/data/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.ldscore.gz \
  --out /path/to/out/trait_h2
```

Rust (release binary):
```
/Users/sharif/Code/ldsc/target/release/ldsc h2 \
  --h2 /path/to/trait.sumstats.gz \
  --ref-ld /Users/sharif/Code/ldsc/data/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.ldscore.gz \
  --w-ld /Users/sharif/Code/ldsc/data/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.ldscore.gz \
  --out /path/to/out/trait_h2
```

### 2) Stratified/LDMS h2 (multi-column LD scores)
Use the LDMS files as the reference LD scores; keep the **weights** as the single-column base LD scores.

Example (8 LDMS bins, rsid):
```
--ref-ld /Users/sharif/Code/ldsc/data/UKBB.ALL.ldscore/UKBB.EUR.8LDMS.rsid.l2.ldscore.gz \
--w-ld   /Users/sharif/Code/ldsc/data/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.ldscore.gz
```

Notes:
- `.l2.M_5_50` is discovered automatically alongside the LD score file. Use `--not-M-5-50` if you want `.l2.M` instead.
- Use `--overlap-annot` only if your annotations overlap. LDMS bins are typically disjoint.

## Practical Checklist
- Ensure sumstats SNP IDs match the LD score file (`rsid` vs `chr:pos:ref:alt`).
- Use `--ref-ld` (non-chr-split) with these files.
- For stratified LD scores (8LDMS/25LDMS), keep weights on the base LD score file.

