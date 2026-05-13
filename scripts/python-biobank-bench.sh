#!/bin/bash
#
# Python LDSC biobank-scale benchmark.
# Runs `python ldsc.py --l2` on biobank_50k (N=50K, 1.66M SNPs) to get a real
# wall-clock measurement of Python LDSC at biobank scale, since the
# previously-cited "1548s extrapolated" number was actually 1000G-scale.
# Designed for long-running on-demand AWS Batch (~9-12 hours expected).
#
set -euo pipefail

S3_BUCKET="${S3_DATA_BUCKET:-ldsc-bench-data-270497617191}"
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
S3_OUT="s3://${S3_BUCKET}/python_biobank_${RUN_ID}"

echo "=== Python LDSC biobank benchmark ==="
echo "RUN_ID: $RUN_ID"
echo "S3 output: $S3_OUT"
echo ""

echo "=== System info ==="
echo "CPU: $(nproc) vCPUs"
echo "CPU model: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs)"
echo "Memory: $(free -h 2>/dev/null | awk '/Mem:/{print $2}' || echo 'unknown')"
echo "Date: $(date -u)"
echo ""

# ── Download biobank data ──────────────────────────────────────────────────
mkdir -p /data
echo "=== Downloading biobank_50k from s3://${S3_BUCKET}/ ==="
for f in biobank_50k.bed biobank_50k.bim biobank_50k.fam; do
    echo "  $f ..."
    aws s3 cp "s3://${S3_BUCKET}/$f" /data/ --quiet
done
echo "Download complete."
ls -lh /data/biobank_50k*
echo ""

# ── Run Python LDSC l2 ──────────────────────────────────────────────────────
# Match the Rust biobank measurement parameters from docs/perf-log.md:
#   --ld-wind-kb 1000 --chunk-size 200 --yes-really
echo "=== Running Python LDSC l2 ==="
echo "Started at $(date -u)"
START=$(date +%s)

python /opt/ldsc_py3/ldsc.py \
    --l2 \
    --bfile /data/biobank_50k \
    --ld-wind-kb 1000 \
    --chunk-size 200 \
    --yes-really \
    --out /tmp/biobank_py3 \
    2>&1 | tail -60

END=$(date +%s)
ELAPSED=$((END - START))
echo ""
echo "=== Finished at $(date -u) ==="
echo "Wall-clock elapsed: ${ELAPSED}s ($((ELAPSED / 60))m $((ELAPSED % 60))s)"
echo "                    $(echo "scale=2; $ELAPSED / 3600" | bc) hours"
echo ""

# ── Upload outputs ──────────────────────────────────────────────────────────
echo "=== Uploading outputs ==="
ls -lh /tmp/biobank_py3* 2>/dev/null | head -30
tar czf /tmp/python_biobank.tar.gz -C /tmp \
    $(cd /tmp && ls biobank_py3*) 2>/dev/null || true
aws s3 cp /tmp/python_biobank.tar.gz "${S3_OUT}/python_biobank.tar.gz" --quiet || true

# Also save a small summary file
cat > /tmp/python_biobank_summary.txt <<EOF
Python LDSC biobank benchmark
==============================
RUN_ID:           ${RUN_ID}
Dataset:          biobank_50k (N=50,000, 1.66M SNPs)
Command:          python ldsc.py --l2 --bfile /data/biobank_50k --ld-wind-kb 1000 --chunk-size 200 --out /tmp/biobank_py3
Hardware:         $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs)
                  $(nproc) vCPUs, $(free -h 2>/dev/null | awk '/Mem:/{print $2}')
Wall-clock:       ${ELAPSED}s ($((ELAPSED / 60))m $((ELAPSED % 60))s)
                  $(echo "scale=2; $ELAPSED / 3600" | bc) hours
EOF
cat /tmp/python_biobank_summary.txt
aws s3 cp /tmp/python_biobank_summary.txt "${S3_OUT}/summary.txt" --quiet

echo ""
echo "=== Done. Outputs at: $S3_OUT ==="
