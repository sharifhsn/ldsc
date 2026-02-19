#!/usr/bin/env bash
# test_rust.sh — Rust-only smoke test for all ldsc subcommands.
#
# Exercises all 5 subcommands and their key flags without Python.
#
# Usage:
#   bash test_rust.sh [--build]
#
# For ldscore/h2/rg tests, expects 1000G reference data at:
#   ../data/1000G_phase3_common_norel.{bed,bim,fam}

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="$SCRIPT_DIR/target/release/ldsc"
BFILE="$SCRIPT_DIR/../data/1000G_phase3_common_norel"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

PASS=0
FAIL=0
declare -a ERRORS=()

pass()  { echo "  PASS: $1"; PASS=$((PASS+1)); }
fail()  { echo "  FAIL: $1"; FAIL=$((FAIL+1)); ERRORS+=("$1"); }
run_test() {
    local desc="$1"; shift
    if "$@" >/dev/null 2>&1; then pass "$desc"; else fail "$desc"; fi
}
assert_file() {
    if [[ -f "$2" ]]; then pass "$1"; else fail "$1 (missing: $2)"; fi
}
section() { echo; echo "── $* ──"; }

# ─── Build ────────────────────────────────────────────────────────────────────
if [[ "${1:-}" == "--build" ]]; then
    echo "Building release binary..."
    (cd "$SCRIPT_DIR" && cargo build --release 2>&1) || { echo "Build FAILED"; exit 1; }
fi

[[ -x "$BIN" ]] || {
    echo "Binary not found: $BIN"
    echo "Run: cargo build --release  (or pass --build)"
    exit 1
}

echo "Testing: $("$BIN" --version 2>&1)"

# ─── Help (all subcommands) ───────────────────────────────────────────────────
section "Help"
run_test "munge-sumstats --help"  "$BIN" munge-sumstats --help
run_test "ldscore --help"         "$BIN" ldscore --help
run_test "h2 --help"              "$BIN" h2 --help
run_test "rg --help"              "$BIN" rg --help
run_test "make-annot --help"      "$BIN" make-annot --help

# ─── munge-sumstats ───────────────────────────────────────────────────────────
section "munge-sumstats"

# BETA + SE → Z (most common case)
{
    printf 'SNP\tA1\tA2\tBETA\tSE\tN\n'
    printf 'rs1\tA\tG\t0.05\t0.01\t50000\n'
    printf 'rs2\tC\tT\t-0.10\t0.02\t50000\n'
    printf 'rs3\tA\tG\t0.20\t0.03\t48000\n'
    printf 'rs4\tC\tT\t0.03\t0.01\t52000\n'
} > "$TMP/ss.txt"
run_test "BETA+SE → Z" \
    "$BIN" munge-sumstats --sumstats "$TMP/ss.txt" --out "$TMP/munge_basic"
assert_file "BETA+SE: output .sumstats.gz" "$TMP/munge_basic.sumstats.gz"

# P + OR → Z with --signed-sumstats
{
    printf 'SNP\tA1\tA2\tOR\tP\tN\n'
    printf 'rs1\tA\tG\t1.10\t0.001\t50000\n'
    printf 'rs2\tC\tT\t0.92\t0.005\t50000\n'
    printf 'rs3\tA\tG\t1.25\t0.0001\t48000\n'
} > "$TMP/ss_or.txt"
run_test "--signed-sumstats OR,1" \
    "$BIN" munge-sumstats --sumstats "$TMP/ss_or.txt" \
        --signed-sumstats "OR,1" --out "$TMP/munge_or"

# P-only with --a1-inc (A1 always increases risk → Z always positive)
{
    printf 'SNP\tA1\tA2\tP\tN\n'
    printf 'rs1\tA\tG\t0.001\t50000\n'
    printf 'rs2\tC\tT\t0.005\t50000\n'
    printf 'rs3\tA\tG\t0.05\t48000\n'
} > "$TMP/ss_p.txt"
run_test "--a1-inc" \
    "$BIN" munge-sumstats --sumstats "$TMP/ss_p.txt" \
        --a1-inc --out "$TMP/munge_a1inc"

# --no-alleles (skip allele columns; no strand-ambiguity filter)
{
    printf 'SNP\tZ\tN\n'
    printf 'rs1\t2.5\t50000\n'
    printf 'rs2\t-1.8\t50000\n'
    printf 'rs3\t3.2\t48000\n'
} > "$TMP/ss_noalleles.txt"
run_test "--no-alleles" \
    "$BIN" munge-sumstats --sumstats "$TMP/ss_noalleles.txt" \
        --no-alleles --out "$TMP/munge_noalleles"

# --keep-maf (retain MAF column in output)
{
    printf 'SNP\tA1\tA2\tBETA\tSE\tN\tMAF\n'
    printf 'rs1\tA\tG\t0.05\t0.01\t50000\t0.30\n'
    printf 'rs2\tC\tT\t-0.10\t0.02\t50000\t0.45\n'
    printf 'rs3\tA\tG\t0.20\t0.03\t50000\t0.05\n'
} > "$TMP/ss_frq.txt"
run_test "--keep-maf" \
    "$BIN" munge-sumstats --sumstats "$TMP/ss_frq.txt" \
        --keep-maf --out "$TMP/munge_frq"

# --n (override sample size constant) + --ignore (drop a redundant column)
{
    printf 'SNP\tA1\tA2\tBETA\tSE\tN\tEXTRA_COL\n'
    printf 'rs1\tA\tG\t0.05\t0.01\t50000\tjunk\n'
    printf 'rs2\tC\tT\t-0.10\t0.02\t50000\tjunk\n'
} > "$TMP/ss_extra.txt"
run_test "--n and --ignore" \
    "$BIN" munge-sumstats --sumstats "$TMP/ss_extra.txt" \
        --n 60000 --ignore "EXTRA_COL" --out "$TMP/munge_n_ignore"

# --info-list (mean of multiple INFO columns for imputation filter)
{
    printf 'SNP\tA1\tA2\tBETA\tSE\tN\tINFO_EUR\tINFO_EAS\n'
    printf 'rs1\tA\tG\t0.05\t0.01\t50000\t0.99\t0.95\n'
    printf 'rs2\tC\tT\t-0.10\t0.02\t50000\t0.85\t0.88\n'
    printf 'rs3\tA\tG\t0.20\t0.03\t50000\t0.70\t0.65\n'
} > "$TMP/ss_info.txt"
run_test "--info-list" \
    "$BIN" munge-sumstats --sumstats "$TMP/ss_info.txt" \
        --info-list "INFO_EUR,INFO_EAS" --out "$TMP/munge_infolist"

# --n-cas-col / --n-con-col (case/control counts summed to N)
{
    printf 'SNP\tA1\tA2\tBETA\tSE\tNCAS\tNCON\n'
    printf 'rs1\tA\tG\t0.05\t0.01\t5000\t45000\n'
    printf 'rs2\tC\tT\t-0.10\t0.02\t5000\t45000\n'
} > "$TMP/ss_cascon.txt"
run_test "--n-cas-col / --n-con-col" \
    "$BIN" munge-sumstats --sumstats "$TMP/ss_cascon.txt" \
        --n-cas-col "NCAS" --n-con-col "NCON" --out "$TMP/munge_cascon"

# --nstudy / --nstudy-min (filter SNPs by number of contributing studies)
{
    printf 'SNP\tA1\tA2\tBETA\tSE\tN\tNSTUDY\n'
    printf 'rs1\tA\tG\t0.05\t0.01\t50000\t5\n'
    printf 'rs2\tC\tT\t-0.10\t0.02\t50000\t2\n'
    printf 'rs3\tA\tG\t0.20\t0.03\t50000\t7\n'
} > "$TMP/ss_nstudy.txt"
run_test "--nstudy / --nstudy-min" \
    "$BIN" munge-sumstats --sumstats "$TMP/ss_nstudy.txt" \
        --nstudy "NSTUDY" --nstudy-min 3 --out "$TMP/munge_nstudy"

# --merge-alleles (harmonise against a reference SNP list)
run_test "--merge-alleles" \
    "$BIN" munge-sumstats --sumstats "$TMP/ss.txt" \
        --merge-alleles "$TMP/ss.txt" --out "$TMP/munge_merge"

# column name overrides: --snp-col, --a1-col, --a2-col, --n-col
{
    printf 'RSID\tEFFECT_ALLELE\tOTHER_ALLELE\tBETA\tSE\tSAMPLE_N\n'
    printf 'rs1\tA\tG\t0.05\t0.01\t50000\n'
    printf 'rs2\tC\tT\t-0.10\t0.02\t50000\n'
} > "$TMP/ss_nonstandard.txt"
run_test "--snp-col / --a1-col / --a2-col / --n-col" \
    "$BIN" munge-sumstats --sumstats "$TMP/ss_nonstandard.txt" \
        --snp-col "RSID" --a1-col "EFFECT_ALLELE" \
        --a2-col "OTHER_ALLELE" --n-col "SAMPLE_N" \
        --out "$TMP/munge_colnames"

# ─── make-annot ───────────────────────────────────────────────────────────────
section "make-annot"

# Minimal BIM file (chr1 and chr2)
{
    printf '1\trs1\t0.0\t100000\tA\tG\n'
    printf '1\trs2\t0.0\t200000\tC\tT\n'
    printf '1\trs3\t0.0\t300000\tA\tG\n'
    printf '1\trs4\t0.0\t400000\tC\tT\n'
    printf '2\trs5\t0.0\t100000\tA\tG\n'
    printf '2\trs6\t0.0\t200000\tC\tT\n'
} > "$TMP/mini.bim"

# BED file (chr1 and chr2 regions)
{
    printf 'chr1\t150000\t250001\n'
    printf 'chr2\t50000\t150001\n'
} > "$TMP/regions.bed"

run_test "BED → annot.gz" \
    "$BIN" make-annot \
        --bimfile "$TMP/mini.bim" \
        --bed-file "$TMP/regions.bed" \
        --annot-file "$TMP/mini.annot.gz"
assert_file "annot.gz output exists" "$TMP/mini.annot.gz"

run_test "--windowsize 50000" \
    "$BIN" make-annot \
        --bimfile "$TMP/mini.bim" \
        --bed-file "$TMP/regions.bed" \
        --windowsize 50000 \
        --annot-file "$TMP/mini_windowed.annot.gz"

run_test "--nomerge" \
    "$BIN" make-annot \
        --bimfile "$TMP/mini.bim" \
        --bed-file "$TMP/regions.bed" \
        --nomerge \
        --annot-file "$TMP/mini_nomerge.annot.gz"

# Gene set + coordinate file (subset of target genes annotated by genomic region)
{
    printf 'GENE1\n'
    printf 'GENE2\n'
} > "$TMP/geneset.txt"
{
    printf 'GENE1\t1\t90000\t210000\n'
    printf 'GENE2\t2\t90000\t210000\n'
    printf 'GENE3\t1\t290000\t310000\n'
} > "$TMP/genecoords.txt"
run_test "--gene-set-file + --gene-coord-file" \
    "$BIN" make-annot \
        --bimfile "$TMP/mini.bim" \
        --gene-set-file "$TMP/geneset.txt" \
        --gene-coord-file "$TMP/genecoords.txt" \
        --annot-file "$TMP/geneset.annot.gz"

# ─── ldscore / h2 / rg ───────────────────────────────────────────────────────
if [[ ! -f "${BFILE}.bim" ]]; then
    echo
    echo "  Skipping ldscore/h2/rg: data not found at ${BFILE}.{bed,bim,fam}"
    echo "  Place 1000G_phase3_common_norel.{bed,bim,fam} in ../data/ to enable."
else
    # Extract chr22 SNP list and BIM (no Python — pure awk)
    awk '$1 == 22 {print $2}' "${BFILE}.bim" > "$TMP/chr22_snps.txt"
    awk '$1 == 22'           "${BFILE}.bim" > "$TMP/chr22.bim"
    N22=$(wc -l < "$TMP/chr22_snps.txt" | tr -d ' ')
    echo
    echo "  Found data: using $N22 chr22 SNPs"

    mkdir -p "$TMP/ld"
    LD="$TMP/ld"

    section "ldscore"

    # --ld-wind-cm (default window; generates per-chr .l2.ldscore.gz/.M/.M_5_50)
    run_test "--ld-wind-cm 1" \
        "$BIN" ldscore \
            --bfile "$BFILE" --extract "$TMP/chr22_snps.txt" \
            --ld-wind-cm 1 --out "$LD/cm1"
    assert_file "chr22 .l2.ldscore.gz" "$LD/cm1.22.l2.ldscore.gz"
    assert_file "chr22 .l2.M"          "$LD/cm1.22.l2.M"
    assert_file "chr22 .l2.M_5_50"     "$LD/cm1.22.l2.M_5_50"

    # --ld-wind-kb (physical-distance window)
    run_test "--ld-wind-kb 1000" \
        "$BIN" ldscore \
            --bfile "$BFILE" --extract "$TMP/chr22_snps.txt" \
            --ld-wind-kb 1000 --out "$LD/kb1000"

    # --ld-wind-snp (fixed flanking-SNP window)
    run_test "--ld-wind-snp 200" \
        "$BIN" ldscore \
            --bfile "$BFILE" --extract "$TMP/chr22_snps.txt" \
            --ld-wind-snp 200 --out "$LD/snp200"

    # --maf (exclude low-frequency SNPs from output)
    run_test "--maf 0.05" \
        "$BIN" ldscore \
            --bfile "$BFILE" --extract "$TMP/chr22_snps.txt" \
            --ld-wind-cm 1 --maf 0.05 --out "$LD/maf05"

    # --print-snps (output subset; all SNPs still contribute to LD windows)
    head -5000 "$TMP/chr22_snps.txt" > "$TMP/chr22_print.txt"
    run_test "--print-snps" \
        "$BIN" ldscore \
            --bfile "$BFILE" --extract "$TMP/chr22_snps.txt" \
            --print-snps "$TMP/chr22_print.txt" \
            --ld-wind-cm 1 --out "$LD/printsnps"

    # --per-allele (weight r² by 2p(1-p) of target SNP)
    run_test "--per-allele" \
        "$BIN" ldscore \
            --bfile "$BFILE" --extract "$TMP/chr22_snps.txt" \
            --ld-wind-cm 1 --per-allele --out "$LD/perallele"

    # --chunk-size (BLAS tile size; larger values change window approximation)
    run_test "--chunk-size 100" \
        "$BIN" ldscore \
            --bfile "$BFILE" --extract "$TMP/chr22_snps.txt" \
            --ld-wind-cm 1 --chunk-size 100 --out "$LD/chunk100"

    # --blas-threads (OpenBLAS thread count)
    run_test "--blas-threads 2" \
        "$BIN" ldscore \
            --bfile "$BFILE" --extract "$TMP/chr22_snps.txt" \
            --ld-wind-cm 1 --blas-threads 2 --out "$LD/blas2"

    # --keep (restrict to a subset of individuals)
    awk 'NR<=100' "${BFILE}.fam" > "$TMP/keep100.txt"
    run_test "--keep 100 individuals" \
        "$BIN" ldscore \
            --bfile "$BFILE" --extract "$TMP/chr22_snps.txt" \
            --keep "$TMP/keep100.txt" \
            --ld-wind-cm 1 --out "$LD/keep100"

    # --annot: partitioned LD scores (one annotation column from a BED region)
    printf 'chr22\t16050000\t17050000\n' > "$TMP/chr22_region.bed"
    if "$BIN" make-annot \
            --bimfile "$TMP/chr22.bim" \
            --bed-file "$TMP/chr22_region.bed" \
            --annot-file "$TMP/chr22annot.22.annot.gz" \
            >/dev/null 2>&1; then
        run_test "--annot (partitioned LD scores)" \
            "$BIN" ldscore \
                --bfile "$BFILE" --extract "$TMP/chr22_snps.txt" \
                --annot "$TMP/chr22annot." \
                --ld-wind-cm 1 --out "$LD/annot"

        # --thin-annot: same data but without CHR/BP/SNP/CM header columns
        zcat "$TMP/chr22annot.22.annot.gz" | cut -f5- \
            > "$TMP/chr22thin.22.annot"
        run_test "--thin-annot" \
            "$BIN" ldscore \
                --bfile "$BFILE" --extract "$TMP/chr22_snps.txt" \
                --annot "$TMP/chr22thin." \
                --thin-annot \
                --ld-wind-cm 1 --out "$LD/thinannot"
    fi

    section "h2 and rg"

    # Synthetic sumstats using chr22 SNP IDs (pure awk — no Python required).
    # Z values cycle [2.1, -1.3, 0.5] — not realistic but exercises all code paths.
    {
        printf 'SNP\tA1\tA2\tZ\tN\n'
        awk 'NR<=5000 {
            z = (NR % 3 == 1) ? 2.1 : (NR % 3 == 2) ? -1.3 : 0.5
            print $1 "\tA\tG\t" z "\t50000"
        }' "$TMP/chr22_snps.txt"
    } | gzip > "$TMP/sim_ss1.sumstats.gz"

    {
        printf 'SNP\tA1\tA2\tZ\tN\n'
        awk 'NR<=5000 {
            z = (NR % 5 == 1) ? 1.8 : (NR % 5 == 2) ? -0.9 : \
                (NR % 5 == 3) ? 1.4 : (NR % 5 == 4) ? -2.0 : 0.3
            print $1 "\tA\tG\t" z "\t60000"
        }' "$TMP/chr22_snps.txt"
    } | gzip > "$TMP/sim_ss2.sumstats.gz"

    # h2: basic
    run_test "h2: basic" \
        "$BIN" h2 \
            --h2 "$TMP/sim_ss1.sumstats.gz" \
            --ref-ld-chr "$LD/cm1." \
            --w-ld-chr   "$LD/cm1." \
            --out "$TMP/h2_basic"

    # h2: --no-intercept (constrain intercept = 1)
    run_test "h2: --no-intercept" \
        "$BIN" h2 \
            --h2 "$TMP/sim_ss1.sumstats.gz" \
            --ref-ld-chr "$LD/cm1." \
            --w-ld-chr   "$LD/cm1." \
            --no-intercept --out "$TMP/h2_noint"

    # h2: --intercept-h2 (fix intercept to specific value)
    run_test "h2: --intercept-h2 1.02" \
        "$BIN" h2 \
            --h2 "$TMP/sim_ss1.sumstats.gz" \
            --ref-ld-chr "$LD/cm1." \
            --w-ld-chr   "$LD/cm1." \
            --intercept-h2 1.02 --out "$TMP/h2_fixedint"

    # h2: --two-step (step 1 estimates intercept on chi2 ≤ 30 SNPs)
    run_test "h2: --two-step 30" \
        "$BIN" h2 \
            --h2 "$TMP/sim_ss1.sumstats.gz" \
            --ref-ld-chr "$LD/cm1." \
            --w-ld-chr   "$LD/cm1." \
            --two-step 30 --out "$TMP/h2_twostep"

    # h2: --chisq-max (filter extreme chi-squared statistics)
    run_test "h2: --chisq-max 80" \
        "$BIN" h2 \
            --h2 "$TMP/sim_ss1.sumstats.gz" \
            --ref-ld-chr "$LD/cm1." \
            --w-ld-chr   "$LD/cm1." \
            --chisq-max 80 --out "$TMP/h2_chisqmax"

    # h2: --n-blocks (jackknife block count)
    run_test "h2: --n-blocks 100" \
        "$BIN" h2 \
            --h2 "$TMP/sim_ss1.sumstats.gz" \
            --ref-ld-chr "$LD/cm1." \
            --w-ld-chr   "$LD/cm1." \
            --n-blocks 100 --out "$TMP/h2_nblocks"

    # h2: --print-cov + --print-delete-vals (diagnostic output)
    run_test "h2: --print-cov + --print-delete-vals" \
        "$BIN" h2 \
            --h2 "$TMP/sim_ss1.sumstats.gz" \
            --ref-ld-chr "$LD/cm1." \
            --w-ld-chr   "$LD/cm1." \
            --print-cov --print-delete-vals --out "$TMP/h2_diag"

    # h2: --not-M-5-50 (use .l2.M instead of .l2.M_5_50 for denominator)
    run_test "h2: --not-M-5-50" \
        "$BIN" h2 \
            --h2 "$TMP/sim_ss1.sumstats.gz" \
            --ref-ld-chr "$LD/cm1." \
            --w-ld-chr   "$LD/cm1." \
            --not-m-5-50 --out "$TMP/h2_notm550"

    # h2: liability-scale conversion (--samp-prev + --pop-prev)
    run_test "h2: --samp-prev + --pop-prev" \
        "$BIN" h2 \
            --h2 "$TMP/sim_ss1.sumstats.gz" \
            --ref-ld-chr "$LD/cm1." \
            --w-ld-chr   "$LD/cm1." \
            --samp-prev 0.5 --pop-prev 0.01 --out "$TMP/h2_liability"

    # rg: basic pair
    run_test "rg: basic pair" \
        "$BIN" rg \
            --rg "$TMP/sim_ss1.sumstats.gz,$TMP/sim_ss2.sumstats.gz" \
            --ref-ld-chr "$LD/cm1." \
            --w-ld-chr   "$LD/cm1." \
            --out "$TMP/rg_basic"

    # rg: --no-intercept (fix gencov intercept = 0)
    run_test "rg: --no-intercept" \
        "$BIN" rg \
            --rg "$TMP/sim_ss1.sumstats.gz,$TMP/sim_ss2.sumstats.gz" \
            --ref-ld-chr "$LD/cm1." \
            --w-ld-chr   "$LD/cm1." \
            --no-intercept --out "$TMP/rg_noint"

    # rg: --intercept-gencov (fix per-pair gencov intercept)
    run_test "rg: --intercept-gencov 0.0" \
        "$BIN" rg \
            --rg "$TMP/sim_ss1.sumstats.gz,$TMP/sim_ss2.sumstats.gz" \
            --ref-ld-chr "$LD/cm1." \
            --w-ld-chr   "$LD/cm1." \
            --intercept-gencov 0.0 --out "$TMP/rg_intgencov"

    # rg: --two-step (two-step estimator for gencov)
    run_test "rg: --two-step 30" \
        "$BIN" rg \
            --rg "$TMP/sim_ss1.sumstats.gz,$TMP/sim_ss2.sumstats.gz" \
            --ref-ld-chr "$LD/cm1." \
            --w-ld-chr   "$LD/cm1." \
            --two-step 30 --out "$TMP/rg_twostep"

    # rg: --samp-prev + --pop-prev (liability-scale rg)
    run_test "rg: --samp-prev + --pop-prev" \
        "$BIN" rg \
            --rg "$TMP/sim_ss1.sumstats.gz,$TMP/sim_ss2.sumstats.gz" \
            --ref-ld-chr "$LD/cm1." \
            --w-ld-chr   "$LD/cm1." \
            --samp-prev 0.5,0.5 --pop-prev 0.01,0.01 --out "$TMP/rg_liability"

    # rg: --print-cov + --print-delete-vals
    run_test "rg: --print-cov + --print-delete-vals" \
        "$BIN" rg \
            --rg "$TMP/sim_ss1.sumstats.gz,$TMP/sim_ss2.sumstats.gz" \
            --ref-ld-chr "$LD/cm1." \
            --w-ld-chr   "$LD/cm1." \
            --print-cov --print-delete-vals --out "$TMP/rg_diag"

    # h2 partitioned (if --annot LD scores were computed)
    if [[ -f "$LD/annot.22.l2.ldscore.gz" ]]; then
        run_test "h2: partitioned (--annot LD scores)" \
            "$BIN" h2 \
                --h2 "$TMP/sim_ss1.sumstats.gz" \
                --ref-ld-chr "$LD/annot." \
                --w-ld-chr   "$LD/cm1." \
                --print-coefficients \
                --out "$TMP/h2_partitioned"
    fi
fi

# ─── Summary ─────────────────────────────────────────────────────────────────
echo
echo "══════════════════════════════"
echo "  Results: $PASS passed, $FAIL failed"
if [[ ${#ERRORS[@]} -gt 0 ]]; then
    echo "  Failed tests:"
    for e in "${ERRORS[@]}"; do echo "    - $e"; done
fi
echo "══════════════════════════════"

[[ $FAIL -eq 0 ]]
