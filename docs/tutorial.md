# ldsc-rs tutorial

This walk-through takes you from a fresh shell to a published-comparable
heritability estimate on real GWAS data, using the same tooling and
canonical reference panel as the original LDSC pipeline. Expected wall
time end-to-end on a 2024 Apple Silicon laptop: **about 4 minutes**
(plus a one-time LD score computation that's ~1 minute on the European
1000 Genomes panel).

If you'd rather not install anything, every step below is also doable
in your browser at <https://sharifhsn.github.io/ldsc/> — but the CLI
gives you the same numerical results and integrates with downstream
pipelines.

## 1. Install

### Option A — Cargo (recommended for Rust users)

```sh
cargo install ldsc                       # from crates.io once published
# or, from source:
git clone https://github.com/sharifhsn/ldsc.git
cd ldsc
cargo build --release --features mimalloc
# binary lands at ./target/release/ldsc
```

### Option B — Pre-built binary

Grab the appropriate archive from the
[GitHub Releases page](https://github.com/sharifhsn/ldsc/releases),
extract `ldsc` (or `ldsc.exe` on Windows), and place it on your `$PATH`.
Binaries are built for Linux x86_64 (musl, statically linked), macOS
arm64, and Windows x86_64.

### Option C — Docker

```sh
docker pull ghcr.io/sharifhsn/ldsc:latest
docker run --rm -v "$PWD:/data" ghcr.io/sharifhsn/ldsc:latest ldsc --help
```

### Sanity check

```sh
ldsc --version          # → ldsc 0.x.y
ldsc l2 --help          # → flag listing for the LD score subcommand
```

## 2. Get the 1000 Genomes European reference panel

LDSC's standard reference panel is the European subset of 1000 Genomes
Phase 3, filtered to common SNPs. If you already have a PLINK BED/BIM/FAM
trio, skip this step.

```sh
# Filtered European-ancestry panel from the original LDSC project:
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz
tar xzf 1000G_Phase3_plinkfiles.tgz
# Or use the trio bundled with the ldsc-rs repository at data/1000G_eur.{bed,bim,fam}
```

## 3. Compute LD scores

```sh
ldsc l2 \
  --bfile data/1000G_eur \
  --ld-wind-kb 1000 \
  --out reference/eur_ld
# wall: ~1 minute on Apple M-series; produces:
#   reference/eur_ld{1..22}.l2.ldscore.gz
#   reference/eur_ld{1..22}.l2.M
#   reference/eur_ld{1..22}.l2.M_5_50
```

If you want exact per-SNP windows (the LDSC paper's mathematical
definition; differs from Python LDSC's chunked approximation by ~12%
on mean L2):

```sh
ldsc l2 --bfile data/1000G_eur --ld-wind-kb 1000 \
        --snp-level-masking --out reference/eur_ld_strict
```

If you want bit-identical Python LDSC output:

```sh
ldsc l2 --bfile data/1000G_eur --ld-wind-kb 1000 \
        --python-compat --out reference/eur_ld_pycompat
```

## 4. Munge a GWAS summary statistics file

Most GWAS results are distributed in tab-separated form with columns
roughly like `SNP A1 A2 BETA SE P N` (with variations). `ldsc
munge-sumstats` standardizes the format and filters to HapMap3 SNPs.

Download a real example — the Yengo et al. 2018 BMI meta-analysis:

```sh
curl -L -o data/BMI_Yengo2018.txt.gz \
  https://portals.broadinstitute.org/collaboration/giant/images/1/14/Bmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz
```

Then munge:

```sh
ldsc munge-sumstats \
  --sumstats data/BMI_Yengo2018.txt.gz \
  --merge-alleles reference/w_hm3.snplist \
  --out sumstats/BMI_Yengo2018
# wall: ~30 seconds; produces sumstats/BMI_Yengo2018.sumstats.gz
```

(The `--merge-alleles` flag filters to HapMap3 SNPs with valid alleles.
Download `w_hm3.snplist` from the LDSC project's reference data bundle
if you don't have it.)

## 5. Estimate heritability

```sh
ldsc h2 \
  --h2 sumstats/BMI_Yengo2018.sumstats.gz \
  --ref-ld-chr reference/eur_ld@ \
  --w-ld-chr reference/eur_ld@ \
  --out h2/BMI_Yengo2018
# wall: ~1-2 seconds
```

Read the output:

```
ldsc h2 --h2 sumstats/BMI_Yengo2018.sumstats.gz ...

Heritability of phenotype 1
---------------------------
Total Observed scale h2:    0.21XX (0.00XX)
Lambda GC:                  1.XX
Mean Chi^2:                 1.XX
Intercept:                  1.0XX (0.0XX)
Ratio:                      0.0XX (0.0XX)
```

You should see `h^2 ≈ 0.21`, consistent with Yengo's published BMI
heritability estimate.

## 6. Genetic correlation (optional)

If you have a second GWAS, compute genetic correlation:

```sh
ldsc rg \
  --rg sumstats/BMI_Yengo2018.sumstats.gz,sumstats/Trait2.sumstats.gz \
  --ref-ld-chr reference/eur_ld@ \
  --w-ld-chr  reference/eur_ld@ \
  --out rg/BMI_vs_trait2
```

Output reports `rg`, SE, P-value, and `gcov_int`.

## 7. Speed knobs for biobank-scale data

If your reference panel has $N >$ 10,000 individuals, the default
exact mode can take an hour. Use CountSketch with masking.

**What `d` means.** CountSketch compresses your N individuals into
`d` "synthetic individuals" by random ±1-weighted hashing: each real
individual is randomly thrown into one of `d` buckets with a random
sign, then the buckets are summed. The compressed data preserves the
*pairwise patterns between SNPs* (which is what LD scores measure)
with per-pair noise ∝ 1/d, while making all the downstream linear
algebra `d`-dimensional instead of `N`-dimensional. So `d` controls
accuracy, and `N` controls the one-time scatter-add cost — they
decouple, which is why the recommended `d` is the same regardless
of how big your panel is.

The universal sweet spot is `--sketch 1000 --snp-level-masking`, which
the N×d sweep (`preprint/data/dn_sweep_full.csv`) verified gives h²
within ~0.003 of per-SNP exact across N = 503 to 100,000:

```sh
ldsc l2 --bfile data/biobank_50k --ld-wind-kb 1000 \
        --sketch 1000 --snp-level-masking \
        --out reference/biobank_ld
# wall: ~10 s on Apple M5 Pro for 1.66M SNPs at N=50K;
#       ~14 s at N=100K
```

Going higher buys diminishing accuracy returns: `--sketch 5000
--snp-level-masking` reaches Pearson r ≥ 0.999 vs exact at ~30 s wall,
appropriate when per-SNP accuracy matters (partitioned heritability,
fine-mapping). Going lower (`--sketch 200`) is ~2 s faster but drifts
h² by ~0.01 from per-SNP exact, so it's only appropriate for
exploratory QC / visualization, not downstream regression.

## 8. No-install alternative

The same pipeline runs in your browser at
<https://sharifhsn.github.io/ldsc/>. Drag-and-drop the BED/BIM/FAM trio
and the `.sumstats.gz` file. All compute happens client-side via
WebAssembly + SharedArrayBuffer worker threads — no upload of genotype
data to any server. Suitable for sensitive data where local install is
inconvenient or prohibited.

## Next steps

- For partitioned heritability across functional annotations, see
  `ldsc h2 --overlap-annot --frqfile-chr ... --print-coefficients`.
- For cell-type-specific heritability (CTS): `ldsc h2-cts ...`.
- For the full set of flags: `ldsc <subcommand> --help`. Every Python
  LDSC flag is accepted with the same semantics.
- Bug reports and questions: <https://github.com/sharifhsn/ldsc/issues>.
