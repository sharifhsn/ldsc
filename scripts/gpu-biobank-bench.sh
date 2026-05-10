#!/bin/bash
#
# GPU biobank benchmark: compare GPU vs CPU on exact path at N=50K.
# GPU only accelerates the exact (non-sketch) GEMM path.
# Designed for AWS Batch on g5.2xlarge (A10G 24GB VRAM, 32GB RAM).
#
set -euo pipefail

BFILE="/data/biobank_50k"
RUNS="${BENCH_RUNS:-3}"
WARMUP="${BENCH_WARMUP:-1}"

echo "=== GPU Biobank 50K Benchmark ==="
echo "Dataset: $BFILE"
echo "Runs: $RUNS, Warmup: $WARMUP"
echo ""

# Check GPU availability
if nvidia-smi >/dev/null 2>&1; then
    echo "GPU detected:"
    nvidia-smi --query-gpu=name,memory.total,compute_cap --format=csv,noheader
else
    echo "WARNING: No GPU detected — GPU benchmarks will fall back to CPU"
fi
echo ""

COMMON="--bfile $BFILE --ld-wind-kb 1000 --chunk-size 200 --yes-really --global-pass"

run_bench() {
    local label="$1"; local extra="$2"; local out_json="/tmp/${label}.json"
    echo "--- $label ---"
    # shellcheck disable=SC2086
    hyperfine --warmup "$WARMUP" --runs "$RUNS" --export-json "$out_json" \
        "ldsc l2 $COMMON $extra --out /tmp/${label}_out"
    echo ""
}

# === CPU baselines (exact path, same --global-pass for fair comparison) ===
run_bench "cpu-exact-f64"       ""
run_bench "cpu-exact-f32"       "--fast-f32"

# === GPU variants ===
run_bench "gpu-exact-f32"       "--gpu"
run_bench "gpu-exact-flex32"    "--gpu --gpu-flex32"

# === CountSketch baselines (no --global-pass; per-chr is their natural mode) ===
# These use the fused BED-decode-normalize-scatter path (no GEMM, no GPU).
COMMON_CS="--bfile $BFILE --ld-wind-kb 1000 --chunk-size 200 --yes-really"

run_bench_cs() {
    local label="$1"; local extra="$2"; local out_json="/tmp/${label}.json"
    echo "--- $label ---"
    # shellcheck disable=SC2086
    hyperfine --warmup "$WARMUP" --runs "$RUNS" --export-json "$out_json" \
        "ldsc l2 $COMMON_CS $extra --out /tmp/${label}_out"
    echo ""
}

run_bench_cs "countsketch-50"      "--sketch 50"
run_bench_cs "countsketch-200"     "--sketch 200"
run_bench_cs "countsketch-1000"    "--sketch 1000"
run_bench_cs "countsketch-5000"    "--sketch 5000"

echo ""
echo "=== RAW JSON ==="
for label in cpu-exact-f64 cpu-exact-f32 \
             gpu-exact-f32 gpu-exact-flex32 \
             countsketch-50 countsketch-200 countsketch-1000 countsketch-5000; do
    f="/tmp/${label}.json"
    echo "--- $label ---"
    cat "$f" 2>/dev/null || echo "(not available)"
    echo ""
done
