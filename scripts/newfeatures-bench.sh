#!/bin/bash
#
# Focused new-feature benchmark: CountSketch and --subsample.
# Designed for AWS Batch with the biobank_50k dataset.
# Use ON_DEMAND to avoid SPOT preemption during multi-hour runs.
#
set -euo pipefail

BFILE="/data/biobank_50k"
RUNS="${BENCH_RUNS:-3}"
WARMUP="${BENCH_WARMUP:-1}"

echo "=== New Features Benchmark (CountSketch + Subsample) ==="
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

# Rademacher sketch baseline (non-fused, restored to correct performance)
run_bench "sketch-50-f32-baseline"   "--sketch 50  --fast-f32"
run_bench "sketch-100-f32-baseline"  "--sketch 100 --fast-f32"

echo "=== CountSketch ==="
echo ""
run_bench "countsketch-50-f32"   "--sketch 50  --fast-f32 --sketch-method countsketch"
run_bench "countsketch-100-f32"  "--sketch 100 --fast-f32 --sketch-method countsketch"
run_bench "countsketch-200-f32"  "--sketch 200 --fast-f32 --sketch-method countsketch"
run_bench "countsketch-50-f64"   "--sketch 50  --sketch-method countsketch"
run_bench "countsketch-100-f64"  "--sketch 100 --sketch-method countsketch"

echo "=== Subsample ==="
echo ""
run_bench "subsample-5k-exact-f32"         "--fast-f32 --subsample 5000"
run_bench "subsample-10k-exact-f32"        "--fast-f32 --subsample 10000"
run_bench "subsample-5k-sketch50-f32"      "--fast-f32 --subsample 5000 --sketch 50"
run_bench "subsample-5k-exact-f64"         "--subsample 5000"
run_bench "subsample-5k-sketch50-f64"      "--subsample 5000 --sketch 50"

echo ""
echo "=== RAW JSON ==="
for label in sketch-50-f32-baseline sketch-100-f32-baseline \
             countsketch-50-f32 countsketch-100-f32 countsketch-200-f32 \
             countsketch-50-f64 countsketch-100-f64 \
             subsample-5k-exact-f32 subsample-10k-exact-f32 \
             subsample-5k-sketch50-f32 subsample-5k-exact-f64 subsample-5k-sketch50-f64; do
    f="/tmp/${label}.json"
    echo "--- $label ---"
    cat "$f" 2>/dev/null || echo "(not available)"
    echo ""
done
