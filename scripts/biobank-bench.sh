#!/bin/bash
#
# Biobank-scale benchmark: exact and CountSketch modes.
# --sketch auto-enables f32 (bit-identical to f64, ~1.3× faster), so no separate
# f64 sketch runs are needed.
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

# Exact baselines (--global-pass required: per-chr at N=50K needs ~70GB, exceeds 28GB container)
run_bench "exact-f64"       "--global-pass"
run_bench "exact-f32"       "--global-pass --fast-f32"

# CountSketch (auto-f32, fused decode+normalize+scatter)
# Per-chr parallel is safe here: d=10000 × 22 chr = ~7GB, well within 28GB.
run_bench "countsketch-50"    "--sketch 50"
run_bench "countsketch-100"   "--sketch 100"
run_bench "countsketch-200"   "--sketch 200"
run_bench "countsketch-500"   "--sketch 500"
run_bench "countsketch-1000"  "--sketch 1000"
run_bench "countsketch-2000"  "--sketch 2000"
run_bench "countsketch-5000"  "--sketch 5000"
run_bench "countsketch-10000" "--sketch 10000"

echo ""
echo "=== RAW JSON ==="
for label in exact-f64 exact-f32 \
             countsketch-50 countsketch-100 countsketch-200 countsketch-500 \
             countsketch-1000 countsketch-2000 countsketch-5000 countsketch-10000; do
    f="/tmp/${label}.json"
    echo "--- $label ---"
    cat "$f" 2>/dev/null || echo "(not available)"
    echo ""
done
