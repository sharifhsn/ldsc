#!/bin/bash
set -euo pipefail

BFILE="/data/1000G_phase3_common_norel"
S3_BUCKET="${S3_DATA_BUCKET:-ldsc-bench-data-270497617191}"
RUNS="${BENCH_RUNS:-10}"
WARMUP="${BENCH_WARMUP:-2}"

# ── Download benchmark data from S3 ──────────────────────────────────────────
echo "=== Downloading benchmark data from s3://${S3_BUCKET}/ ==="
mkdir -p /data
aws s3 cp "s3://${S3_BUCKET}/1000G_phase3_common_norel.bed" /data/ --quiet
aws s3 cp "s3://${S3_BUCKET}/1000G_phase3_common_norel.bim" /data/ --quiet
aws s3 cp "s3://${S3_BUCKET}/1000G_phase3_common_norel.fam" /data/ --quiet
aws s3 cp "s3://${S3_BUCKET}/bench_5k.bed" /data/ --quiet
aws s3 cp "s3://${S3_BUCKET}/bench_5k.bim" /data/ --quiet
aws s3 cp "s3://${S3_BUCKET}/bench_5k.fam" /data/ --quiet
echo "Download complete."
echo ""

# ── System info ───────────────────────────────────────────────────────────────
echo "=== LDSC Benchmark ==="
echo "Binary: $(ldsc --version 2>&1 || echo 'unknown')"
echo "CPU: $(nproc) vCPUs"
echo "CPU Model: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs)"
echo "Memory: $(free -h 2>/dev/null | awk '/Mem:/{print $2}' || echo 'unknown')"
echo "Date: $(date -u)"
echo ""

# ── Run benchmark ─────────────────────────────────────────────────────────────
if [ $# -eq 0 ]; then
    # Default: full 1000G benchmark with verbose timing
    echo "=== Running default benchmark (${RUNS} runs, ${WARMUP} warmup) ==="
    hyperfine \
        --warmup "$WARMUP" \
        --runs "$RUNS" \
        --export-json /tmp/results.json \
        "ldsc l2 --bfile $BFILE --ld-wind-kb 1000 --chunk-size 200 --out /tmp/bench --yes-really --verbose-timing"

    echo ""
    echo "=== RESULTS JSON ==="
    cat /tmp/results.json
else
    # Custom command passed as arguments
    echo "=== Running custom command ==="
    exec "$@"
fi
