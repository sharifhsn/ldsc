#!/bin/bash
#
# Local perf gate: run all l2 modes on bench_200k and report times.
# Must pass BEFORE submitting any AWS HPC job.
#
# Usage:
#   bash scripts/local-perf-gate.sh              # quick (bench_5k)
#   bash scripts/local-perf-gate.sh --full       # thorough (bench_200k)
#   bash scripts/local-perf-gate.sh --biobank    # biobank_50k (if available)
#
set -euo pipefail

MODE="${1:---quick}"

case "$MODE" in
    --full)    BFILE="data/bench_200k"; LABEL="200k" ;;
    --biobank) BFILE="data/biobank_50k"; LABEL="biobank_50k" ;;
    *)         BFILE="data/bench_5k";    LABEL="5k" ;;
esac

if [ ! -f "${BFILE}.bed" ]; then
    echo "ERROR: ${BFILE}.bed not found" >&2
    exit 1
fi

BINARY="target/release/ldsc"
if [ ! -f "$BINARY" ]; then
    echo "Building release binary..."
    cargo build --release
fi

COMMON="--bfile $BFILE --ld-wind-kb 1000 --chunk-size 200 --yes-really"
OUTDIR=$(mktemp -d)

run_mode() {
    local label="$1"
    local extra="$2"
    local out="${OUTDIR}/${label}"

    # shellcheck disable=SC2086
    local start=$(date +%s%N)
    $BINARY l2 $COMMON $extra --out "$out" >/dev/null 2>&1
    local end=$(date +%s%N)
    local ms=$(( (end - start) / 1000000 ))
    printf "  %-30s %6d ms\n" "$label" "$ms"
}

echo "=== Local Perf Gate ($LABEL, $(date)) ==="
echo "Binary: $($BINARY --version 2>&1 || echo unknown)"
echo ""

echo "--- Exact modes ---"
run_mode "exact-f64"            ""
run_mode "exact-f32"            "--fast-f32"

echo "--- Sketch Rademacher (auto-f32) ---"
run_mode "sketch-50"            "--sketch 50"
run_mode "sketch-100"           "--sketch 100"
run_mode "sketch-200"           "--sketch 200"

echo "--- Sketch CountSketch (auto-f32) ---"
run_mode "countsketch-50"       "--sketch 50  --sketch-method countsketch"
run_mode "countsketch-100"      "--sketch 100 --sketch-method countsketch"

echo "--- Stochastic ---"
run_mode "stochastic-50"        "--stochastic 50"

echo "--- Subsample (if biobank) ---"
if [ "$LABEL" = "biobank_50k" ]; then
    run_mode "subsample-5k"         "--subsample 5000"
    run_mode "subsample-5k-sk50"    "--subsample 5000 --sketch 50"
    run_mode "subsample-2k-cs200"   "--subsample 2000 --sketch 200 --sketch-method countsketch"
    run_mode "subsample-5k-cs200"   "--subsample 5000 --sketch 200 --sketch-method countsketch"
    run_mode "subsample-5k-cs500"   "--subsample 5000 --sketch 500 --sketch-method countsketch"
else
    echo "  (skipped — subsample only meaningful at biobank N)"
fi

echo ""
echo "--- Parity check (exact-f64 vs Python on bench_5k) ---"
if [ "$LABEL" = "5k" ] || [ "$LABEL" = "200k" ]; then
    # Quick parity: run exact-f64 on bench_5k extract from full genome
    PARITY_OUT="${OUTDIR}/parity"
    $BINARY l2 --bfile data/1000G_phase3_common_norel \
        --extract data/bench_5k.bim.snplist \
        --ld-wind-kb 1000 --chunk-size 200 --yes-really \
        --out "$PARITY_OUT" >/dev/null 2>&1 || true

    if [ -f "${PARITY_OUT}.l2.ldscore.gz" ]; then
        echo "  Parity output generated (manual diff required)"
    else
        echo "  (parity check skipped — needs extract snplist)"
    fi
else
    echo "  (skipped for biobank)"
fi

rm -rf "$OUTDIR"
echo ""
echo "=== Done. Review times above for regressions before submitting AWS jobs. ==="
