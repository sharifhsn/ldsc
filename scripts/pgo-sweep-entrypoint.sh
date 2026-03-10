#!/bin/bash
set -euo pipefail

BFILE="/data/1000G_phase3_common_norel"
S3_BUCKET="${S3_DATA_BUCKET:-ldsc-bench-data-270497617191}"
RUNS="${BENCH_RUNS:-10}"
WARMUP="${BENCH_WARMUP:-2}"
CHUNK_SIZES="${CHUNK_SIZES:-100,200,400,800}"
TARGET="x86_64-unknown-linux-musl"

# ── Download benchmark data from S3 ──────────────────────────────────────────
echo "=== Downloading benchmark data from s3://${S3_BUCKET}/ ==="
mkdir -p /data
aws s3 cp "s3://${S3_BUCKET}/1000G_phase3_common_norel.bed" /data/ --quiet
aws s3 cp "s3://${S3_BUCKET}/1000G_phase3_common_norel.bim" /data/ --quiet
aws s3 cp "s3://${S3_BUCKET}/1000G_phase3_common_norel.fam" /data/ --quiet
echo "Download complete."
echo ""

# ── System info ───────────────────────────────────────────────────────────────
echo "=== PGO Sweep Benchmark ==="
echo "CPU: $(nproc) vCPUs"
echo "CPU Model: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs)"
echo "Memory: $(free -h 2>/dev/null | awk '/Mem:/{print $2}' || echo 'unknown')"
echo "Date: $(date -u)"
echo "Chunk sizes: ${CHUNK_SIZES}"
echo ""

# ── Step 1: Build standard (non-PGO) binary ──────────────────────────────────
echo "=== Step 1: Building standard binary ==="
cd /build
rm -f target/${TARGET}/release/ldsc target/${TARGET}/release/deps/ldsc-*
cargo build --release --features mimalloc --target ${TARGET} 2>&1 | tail -3
strip target/${TARGET}/release/ldsc
cp target/${TARGET}/release/ldsc /usr/local/bin/ldsc-std
echo "Standard binary built."
echo ""

# ── Step 2: Build instrumented binary for PGO ────────────────────────────────
echo "=== Step 2: Building PGO-instrumented binary ==="
PGO_DIR="/tmp/pgo-profiles"
mkdir -p "$PGO_DIR"
rm -f target/${TARGET}/release/ldsc target/${TARGET}/release/deps/ldsc-*
RUSTFLAGS="-Cprofile-generate=$PGO_DIR" \
    cargo build --release --features mimalloc --target ${TARGET} 2>&1 | tail -3
echo "Instrumented binary built."
echo ""

# ── Step 3: PGO training on full 1.66M genome ────────────────────────────────
echo "=== Step 3: PGO training (full 1000G, 1.66M SNPs, --ld-wind-kb 1000) ==="
echo "This takes ~60-90s..."
time target/${TARGET}/release/ldsc l2 \
    --bfile "$BFILE" \
    --ld-wind-kb 1000 \
    --chunk-size 200 \
    --out /tmp/pgo_train \
    --yes-really \
    --verbose-timing 2>&1
echo ""

PROFRAW_COUNT=$(ls "$PGO_DIR"/*.profraw 2>/dev/null | wc -l)
echo "Profile files generated: $PROFRAW_COUNT"

# ── Step 4: Merge profiles and rebuild ────────────────────────────────────────
echo "=== Step 4: Merging profiles and rebuilding with PGO ==="
llvm-profdata merge -output="$PGO_DIR/merged.profdata" "$PGO_DIR"/*.profraw
echo "Profiles merged."

rm -f target/${TARGET}/release/ldsc target/${TARGET}/release/deps/ldsc-*
RUSTFLAGS="-Cprofile-use=$PGO_DIR/merged.profdata -Cllvm-args=-pgo-warn-missing-function" \
    cargo build --release --features mimalloc --target ${TARGET} 2>&1 | tail -5
strip target/${TARGET}/release/ldsc
cp target/${TARGET}/release/ldsc /usr/local/bin/ldsc-pgo
echo "PGO binary built."
echo ""

# ── Step 5: Run benchmarks ───────────────────────────────────────────────────
echo "=== Step 5: Benchmark sweep ==="
echo ""

IFS=',' read -ra CHUNKS <<< "$CHUNK_SIZES"

# Standard binary benchmarks
for chunk in "${CHUNKS[@]}"; do
    echo "--- standard, chunk_size=${chunk} ---"
    cp /usr/local/bin/ldsc-std /usr/local/bin/ldsc
    hyperfine \
        --warmup "$WARMUP" \
        --runs "$RUNS" \
        --export-json "/tmp/results-std-c${chunk}.json" \
        "ldsc l2 --bfile $BFILE --ld-wind-kb 1000 --chunk-size ${chunk} --out /tmp/bench --yes-really --verbose-timing"
    echo ""
done

# PGO binary benchmarks
for chunk in "${CHUNKS[@]}"; do
    echo "--- PGO, chunk_size=${chunk} ---"
    cp /usr/local/bin/ldsc-pgo /usr/local/bin/ldsc
    hyperfine \
        --warmup "$WARMUP" \
        --runs "$RUNS" \
        --export-json "/tmp/results-pgo-c${chunk}.json" \
        "ldsc l2 --bfile $BFILE --ld-wind-kb 1000 --chunk-size ${chunk} --out /tmp/bench --yes-really --verbose-timing"
    echo ""
done

# ── Summary table ─────────────────────────────────────────────────────────────
echo ""
echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║              COMPREHENSIVE PGO + CHUNK SIZE SWEEP                  ║"
echo "║  AWS c7a.4xlarge (EPYC 7R13, 16 vCPU)                             ║"
echo "║  Full 1000G (1.66M SNPs), --ld-wind-kb 1000, ${RUNS} runs         ║"
echo "╠══════════════════════════════════════════════════════════════════════╣"
printf "║  %-12s │ %-22s │ %-22s ║\n" "chunk_size" "standard (mean±σ)" "PGO (mean±σ)"
echo "╠══════════════════════════════════════════════════════════════════════╣"

for chunk in "${CHUNKS[@]}"; do
    std_result=$(python3 -c "
import json
d = json.load(open('/tmp/results-std-c${chunk}.json'))
r = d['results'][0]
print(f\"{r['mean']:.3f}s ± {r['stddev']:.3f}s\")
" 2>/dev/null || echo "N/A")

    pgo_result=$(python3 -c "
import json
d = json.load(open('/tmp/results-pgo-c${chunk}.json'))
r = d['results'][0]
print(f\"{r['mean']:.3f}s ± {r['stddev']:.3f}s\")
" 2>/dev/null || echo "N/A")

    printf "║  %-12s │ %-22s │ %-22s ║\n" "$chunk" "$std_result" "$pgo_result"
done

echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Baseline (non-PGO, chunk=200, prev run): 43.6s ± 0.11s"
echo "Python reference: 1548.5s"

# ── Dump all JSON results ─────────────────────────────────────────────────────
echo ""
echo "=== RAW RESULTS JSON ==="
for chunk in "${CHUNKS[@]}"; do
    echo "--- std c${chunk} ---"
    cat "/tmp/results-std-c${chunk}.json" 2>/dev/null || echo "{}"
    echo ""
    echo "--- pgo c${chunk} ---"
    cat "/tmp/results-pgo-c${chunk}.json" 2>/dev/null || echo "{}"
    echo ""
done
