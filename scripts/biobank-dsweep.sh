#!/bin/bash
#
# Biobank d-sweep: produce LD scores at multiple sketch d values for V1
# (sketch only) and V3 (sketch-hybrid-fused), plus exact references.
# Each variant's outputs are bundled into a single tarball and uploaded to
# S3 (instead of 23 separate per-chr files) — fast enough to survive short
# spot sessions. Designed for AWS Batch with biobank_50k.
#
set -euo pipefail

BFILE="/data/biobank_50k"
S3_BUCKET="${S3_DATA_BUCKET:-ldsc-bench-data-270497617191}"
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
S3_OUT="s3://${S3_BUCKET}/biobank_dsweep_${RUN_ID}"

echo "=== Biobank d-sweep ==="
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

    # Pack all output files for this variant into one tarball — one S3 upload
    # instead of 67 per-file uploads. Fast enough to fit in short spot sessions.
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

    # Free disk for the next variant (biobank ldsc outputs ≈ 150MB per variant)
    rm -f "${prefix}".l2.ldscore.gz "${prefix}".l2.M "${prefix}".l2.M_5_50
    rm -f "${prefix}"[1-9].l2.ldscore.gz "${prefix}"[1-9].l2.M "${prefix}"[1-9].l2.M_5_50
    rm -f "${prefix}"[1-9][0-9].l2.ldscore.gz "${prefix}"[1-9][0-9].l2.M "${prefix}"[1-9][0-9].l2.M_5_50
    rm -f "$tar_path"
    echo ""
}

# V1: sketch-only at varying d
run_and_upload "v1_d200"  "--sketch 200"
run_and_upload "v1_d500"  "--sketch 500"
run_and_upload "v1_d1000" "--sketch 1000"
run_and_upload "v1_d1600" "--sketch 1600"

# V3: hybrid-fused (exact B×B, sketched A×B) at varying d
run_and_upload "v3_d200"  "--sketch 200 --sketch-hybrid --sketch-hybrid-fused"
run_and_upload "v3_d500"  "--sketch 500 --sketch-hybrid --sketch-hybrid-fused"
run_and_upload "v3_d1000" "--sketch 1000 --sketch-hybrid --sketch-hybrid-fused"

# Truth references at N=50K. --global-pass required for exact paths
# (per-chr parallel exact OOMs at biobank scale per biobank-bench.sh).
run_and_upload "chunk_exact" "--fast-f32 --global-pass"
run_and_upload "snp_exact"   "--snp-level-masking --fast-f32 --global-pass"

echo ""
echo "=== Done. All outputs at: $S3_OUT ==="
echo "Download with:"
echo "  aws s3 cp $S3_OUT . --recursive --profile AdministratorAccess-270497617191"
