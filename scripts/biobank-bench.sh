#!/bin/bash
#
# Biobank-scale benchmark: exact vs sketch at various d values.
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

COMMON="--bfile $BFILE --ld-wind-kb 1000 --chunk-size 200 --yes-really --verbose-timing"

# Exact baseline
echo "--- Exact ---"
hyperfine --warmup "$WARMUP" --runs "$RUNS" --export-json /tmp/exact.json \
    "ldsc l2 $COMMON --out /tmp/exact_out"

# Sketch sweep
for d in 50 100 200; do
    echo "--- Sketch d=$d ---"
    hyperfine --warmup "$WARMUP" --runs "$RUNS" --export-json "/tmp/s${d}.json" \
        "ldsc l2 $COMMON --sketch $d --out /tmp/s${d}_out"
done

echo ""
echo "=== RAW JSON ==="
for f in /tmp/exact.json /tmp/s50.json /tmp/s100.json /tmp/s200.json; do
    echo "--- $(basename "$f") ---"
    cat "$f" 2>/dev/null || echo "(not available)"
    echo ""
done
