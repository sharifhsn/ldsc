#!/bin/bash
#
# Biobank-scale benchmark: exact vs sketch × {f64, f32} at various d values.
# Designed for AWS Batch with the biobank_50k dataset (1.66M SNPs × 50K individuals).
#
set -euo pipefail

BFILE="/data/biobank_50k"
RUNS="${BENCH_RUNS:-3}"
WARMUP="${BENCH_WARMUP:-1}"

echo "=== Biobank 50K Benchmark ==="
echo "Dataset: $BFILE"
echo "Runs: $RUNS, Warmup: $WARMUP"
echo ""

COMMON="--bfile $BFILE --ld-wind-kb 1000 --chunk-size 200 --yes-really"

run_bench() {
    local label="$1"; local extra="$2"; local out_json="/tmp/${label}.json"
    echo "--- $label ---"
    # shellcheck disable=SC2086
    hyperfine --warmup "$WARMUP" --runs "$RUNS" --export-json "$out_json" \
        "ldsc l2 $COMMON $extra --out /tmp/${label}_out"
    echo ""
}

run_bench "exact-f64"       ""
run_bench "exact-f32"       "--fast-f32"
run_bench "sketch-50-f64"   "--sketch 50"
run_bench "sketch-100-f64"  "--sketch 100"
run_bench "sketch-200-f64"  "--sketch 200"
run_bench "sketch-50-f32"   "--sketch 50  --fast-f32"
run_bench "sketch-100-f32"  "--sketch 100 --fast-f32"
run_bench "sketch-200-f32"  "--sketch 200 --fast-f32"

echo ""
echo "=== RAW JSON ==="
for label in exact-f64 exact-f32 sketch-50-f64 sketch-100-f64 sketch-200-f64 \
             sketch-50-f32 sketch-100-f32 sketch-200-f32; do
    f="/tmp/${label}.json"
    echo "--- $label ---"
    cat "$f" 2>/dev/null || echo "(not available)"
    echo ""
done
