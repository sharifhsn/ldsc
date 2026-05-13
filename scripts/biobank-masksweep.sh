#!/bin/bash
#
# Biobank sketch-with-masking sweep: tests `--sketch d --snp-level-masking`
# combinations to see if sketch + per-SNP windowing reaches the per-SNP-exact
# truth at sketch speed. Designed for AWS Batch with biobank_50k.
#
# RUN_ID prefix must start with "biobank_dsweep_" to satisfy the existing
# instance-role S3 PutObject policy scope.
#
set -euo pipefail

BFILE="/data/biobank_50k"
S3_BUCKET="${S3_DATA_BUCKET:-ldsc-bench-data-270497617191}"
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)_masked}"
S3_OUT="s3://${S3_BUCKET}/biobank_dsweep_${RUN_ID}"

echo "=== Biobank mask-sweep ==="
echo "Dataset: $BFILE"
echo "Output:  $S3_OUT"
echo ""

COMMON="--bfile $BFILE --ld-wind-kb 1000 --chunk-size 200 --yes-really --mmap"

run_and_upload() {
    local label="$1"; local flags="$2"; shift 2
    local prefix="/tmp/${label}"
    echo "--- $label ($flags) ---"
    local t0=$SECONDS

    # shellcheck disable=SC2086
    ldsc l2 $COMMON $flags --out "$prefix"
    local t_ldsc=$((SECONDS - t0))
    echo "elapsed_ldsc ${label}: ${t_ldsc}s"

    local tar_path="/tmp/${label}.tar.gz"
    (cd /tmp && tar czf "${label}.tar.gz" \
        "${label}.l2.ldscore.gz" "${label}.l2.M" "${label}.l2.M_5_50" \
        "${label}"[1-9].l2.ldscore.gz "${label}"[1-9].l2.M "${label}"[1-9].l2.M_5_50 \
        "${label}"[1-9][0-9].l2.ldscore.gz "${label}"[1-9][0-9].l2.M "${label}"[1-9][0-9].l2.M_5_50 \
        2>/dev/null || true)
    local tar_size; tar_size=$(stat -c %s "$tar_path" 2>/dev/null || stat -f %z "$tar_path" 2>/dev/null || echo 0)
    echo "tarball ${label}.tar.gz: ${tar_size} bytes"

    aws s3 cp "$tar_path" "${S3_OUT}/${label}.tar.gz" --quiet
    local t_total=$((SECONDS - t0))
    echo "elapsed_total ${label}: ${t_total}s (ldsc ${t_ldsc}s, upload $((t_total - t_ldsc))s)"

    rm -f "${prefix}".l2.ldscore.gz "${prefix}".l2.M "${prefix}".l2.M_5_50
    rm -f "${prefix}"[1-9].l2.ldscore.gz "${prefix}"[1-9].l2.M "${prefix}"[1-9].l2.M_5_50
    rm -f "${prefix}"[1-9][0-9].l2.ldscore.gz "${prefix}"[1-9][0-9].l2.M "${prefix}"[1-9][0-9].l2.M_5_50
    rm -f "$tar_path"
    echo ""
}

# Sketch + per-SNP masking sweep over d
run_and_upload "v1m_d200"  "--sketch 200  --snp-level-masking"
run_and_upload "v1m_d500"  "--sketch 500  --snp-level-masking"
run_and_upload "v1m_d1000" "--sketch 1000 --snp-level-masking"
run_and_upload "v1m_d1600" "--sketch 1600 --snp-level-masking"

echo ""
echo "=== Done. All outputs at: $S3_OUT ==="
