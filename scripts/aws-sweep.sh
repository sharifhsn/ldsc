#!/bin/bash
#
# Comprehensive PGO + chunk-size sweep on AWS Batch.
#
# Builds a fat Docker image with Rust toolchain + llvm-profdata, submits
# a single AWS Batch job that:
#   1. Builds standard binary
#   2. Builds PGO-instrumented binary → trains on full 1.66M genome → rebuilds
#   3. Runs chunk-size sweep for both standard and PGO binaries
#
# Usage:
#   ./scripts/aws-sweep.sh                    # full sweep (chunk 100,200,400,800)
#   ./scripts/aws-sweep.sh --skip-build       # reuse existing image
#   ./scripts/aws-sweep.sh --chunks 200,400   # custom chunk sizes
#   ./scripts/aws-sweep.sh --runs 5           # fewer runs per config
#
set -euo pipefail

PROFILE="AdministratorAccess-270497617191"
REGION="us-east-1"
ACCOUNT_ID="270497617191"
PREFIX="ldsc-bench"
ECR_REPO="${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com/${PREFIX}"
S3_BUCKET="${PREFIX}-data-${ACCOUNT_ID}"
DNS_SERVER="10.64.0.1"

AWS="aws --profile $PROFILE --region $REGION --no-cli-pager"

SKIP_BUILD=false
CHUNK_SIZES="100,200,400,800"
RUNS=10
WARMUP=2
VCPUS=16
MEMORY=16384

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-build) SKIP_BUILD=true; shift ;;
        --chunks)     CHUNK_SIZES="$2"; shift 2 ;;
        --runs)       RUNS="$2"; shift 2 ;;
        --warmup)     WARMUP="$2"; shift 2 ;;
        --vcpus)      VCPUS="$2"; shift 2 ;;
        --memory)     MEMORY="$2"; shift 2 ;;
        *)            echo "Unknown option: $1"; exit 1 ;;
    esac
done

TIMESTAMP=$(date +%Y%m%d-%H%M%S)
IMAGE_TAG="pgo-sweep-${TIMESTAMP}"
DEPS_IMAGE="${PREFIX}-deps:latest"

# ── ECR Login ─────────────────────────────────────────────────────────────────
$AWS ecr get-login-password | docker login --username AWS --password-stdin \
    "${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com" 2>/dev/null

# ── Build ─────────────────────────────────────────────────────────────────────
if [ "$SKIP_BUILD" = false ]; then
    echo "=== Building PGO sweep image (fat: Rust + llvm + hyperfine + awscli) ==="

    # Check deps image exists
    if ! docker image inspect "$DEPS_IMAGE" >/dev/null 2>&1; then
        echo "ERROR: Deps image '$DEPS_IMAGE' not found. Run scripts/aws-build-deps.sh first." >&2
        exit 1
    fi

    # Step 1: Add llvm-profdata to deps image + copy source + compile
    echo "  [1/3] Compiling src/ with llvm-profdata..."
    docker rm -f ldsc-pgo-build 2>/dev/null || true
    # Start container from deps image (has cached target/)
    docker create --dns "$DNS_SERVER" --name ldsc-pgo-build \
        -w /build \
        -e CARGO_INCREMENTAL=0 \
        "$DEPS_IMAGE" \
        sh -c '
            apk add --no-cache llvm bash &&
            rm -f target/x86_64-unknown-linux-musl/release/ldsc target/x86_64-unknown-linux-musl/release/deps/ldsc-* &&
            echo "=== Deps already cached, compiling src/ only ===" &&
            cargo build --release --features mimalloc --target x86_64-unknown-linux-musl 2>&1 | tail -5 &&
            echo "=== Build complete ==="
        '
    # Copy source files INTO the container (not bind mount — so they persist in commit)
    # docker cp copies src/ INTO existing /build/src/ as /build/src/src/ — avoid this
    # by copying contents via tar
    tar -cf - -C "$(pwd)" src Cargo.toml Cargo.lock | docker cp - ldsc-pgo-build:/build/
    # Now start the container to compile
    docker start -a ldsc-pgo-build
    # Commit the container (has compiled target/ + llvm-profdata + source)
    docker commit ldsc-pgo-build "${PREFIX}-pgo-builder:latest"
    docker rm -f ldsc-pgo-build

    # Step 2: From pgo-builder, add hyperfine + awscli + entrypoint
    echo "  [2/3] Adding hyperfine + awscli + entrypoint..."
    docker rm -f ldsc-pgo-runtime 2>/dev/null || true
    docker run --dns "$DNS_SERVER" --name ldsc-pgo-runtime \
        -v "$(pwd)/scripts/pgo-sweep-entrypoint.sh:/tmp/entrypoint.sh:ro" \
        "${PREFIX}-pgo-builder:latest" \
        sh -c '
            apk add --no-cache curl python3 py3-pip &&
            curl -sL "https://github.com/sharkdp/hyperfine/releases/download/v1.19.0/hyperfine-v1.19.0-x86_64-unknown-linux-musl.tar.gz" \
                | tar xz -C /usr/local/bin --strip-components=1 hyperfine-v1.19.0-x86_64-unknown-linux-musl/hyperfine &&
            pip install --break-system-packages awscli &&
            cp /tmp/entrypoint.sh /entrypoint.sh &&
            chmod +x /entrypoint.sh
        '
    docker commit \
        --change 'ENTRYPOINT ["/entrypoint.sh"]' \
        --change 'WORKDIR /build' \
        ldsc-pgo-runtime "${PREFIX}:${IMAGE_TAG}"
    docker rm -f ldsc-pgo-runtime

    echo "  [3/3] Pushing to ECR..."
    docker tag "${PREFIX}:${IMAGE_TAG}" "${ECR_REPO}:${IMAGE_TAG}"
    docker push "${ECR_REPO}:${IMAGE_TAG}" 2>&1 | tail -5
    echo "Image pushed: ${ECR_REPO}:${IMAGE_TAG}"
else
    echo "=== Skipping build (--skip-build) ==="
    # Find the most recent pgo-sweep tag
    IMAGE_TAG=$($AWS ecr describe-images --repository-name "${PREFIX}" \
        --query "sort_by(imageDetails, &imagePushedAt)[-1].imageTags[0]" \
        --output text 2>/dev/null || echo "latest")
    echo "Using image tag: $IMAGE_TAG"
fi

# ── Register job definition ──────────────────────────────────────────────────
echo ""
echo "=== Registering job definition ==="
$AWS batch register-job-definition \
    --job-definition-name "${PREFIX}-job" \
    --type container \
    --container-properties "{
        \"image\": \"${ECR_REPO}:${IMAGE_TAG}\",
        \"resourceRequirements\": [
            {\"type\": \"VCPU\", \"value\": \"${VCPUS}\"},
            {\"type\": \"MEMORY\", \"value\": \"${MEMORY}\"}
        ]
    }" --query "jobDefinitionArn" --output text

# ── Submit job ────────────────────────────────────────────────────────────────
JOB_NAME="pgo-sweep-${TIMESTAMP}"
echo ""
echo "=== Submitting PGO sweep job: $JOB_NAME ==="
echo "Chunk sizes: $CHUNK_SIZES"
echo "Runs per config: $RUNS (warmup: $WARMUP)"
echo ""

OVERRIDES=$(cat << EOF
{
  "resourceRequirements": [
    {"type": "VCPU", "value": "${VCPUS}"},
    {"type": "MEMORY", "value": "${MEMORY}"}
  ],
  "environment": [
    {"name": "BENCH_RUNS", "value": "${RUNS}"},
    {"name": "BENCH_WARMUP", "value": "${WARMUP}"},
    {"name": "CHUNK_SIZES", "value": "${CHUNK_SIZES}"},
    {"name": "S3_DATA_BUCKET", "value": "${S3_BUCKET}"}
  ]
}
EOF
)

JOB_ID=$($AWS batch submit-job \
    --job-name "$JOB_NAME" \
    --job-queue "${PREFIX}-queue" \
    --job-definition "${PREFIX}-job" \
    --container-overrides "$OVERRIDES" \
    --query "jobId" --output text)

echo "Job ID: $JOB_ID"
echo "Console: https://${REGION}.console.aws.amazon.com/batch/home?region=${REGION}#jobs/detail/${JOB_ID}"

# ── Stream logs ──────────────────────────────────────────────────────────────
echo ""
echo "=== Waiting for job to start ==="
LOG_STREAM=""
PREV_STATUS=""

while true; do
    JOB_JSON=$($AWS batch describe-jobs --jobs "$JOB_ID")
    STATUS=$(echo "$JOB_JSON" | python3 -c "import json,sys; print(json.load(sys.stdin)['jobs'][0]['status'])")

    if [ "$STATUS" != "$PREV_STATUS" ]; then
        echo "$(date +%H:%M:%S) Status: $STATUS"
        PREV_STATUS="$STATUS"
    fi

    case "$STATUS" in
        SUCCEEDED|FAILED)
            LOG_STREAM=$(echo "$JOB_JSON" | python3 -c "
import json, sys
j = json.load(sys.stdin)['jobs'][0]
print(j.get('container', {}).get('logStreamName', ''))" 2>/dev/null || true)
            break
            ;;
        RUNNING)
            LOG_STREAM=$(echo "$JOB_JSON" | python3 -c "
import json, sys
j = json.load(sys.stdin)['jobs'][0]
print(j.get('container', {}).get('logStreamName', ''))" 2>/dev/null || true)
            if [ -n "$LOG_STREAM" ]; then
                break
            fi
            ;;
    esac
    sleep 5
done

if [ -n "$LOG_STREAM" ]; then
    echo ""
    echo "=== Streaming logs (${LOG_STREAM}) ==="
    echo ""

    NEXT_TOKEN=""
    while true; do
        CUR_STATUS=$($AWS batch describe-jobs --jobs "$JOB_ID" \
            --query "jobs[0].status" --output text)

        if [ -z "$NEXT_TOKEN" ]; then
            LOG_OUTPUT=$($AWS logs get-log-events \
                --log-group-name "/aws/batch/${PREFIX}" \
                --log-stream-name "$LOG_STREAM" \
                --start-from-head \
                --output json 2>/dev/null || echo '{"events":[],"nextForwardToken":""}')
        else
            LOG_OUTPUT=$($AWS logs get-log-events \
                --log-group-name "/aws/batch/${PREFIX}" \
                --log-stream-name "$LOG_STREAM" \
                --next-token "$NEXT_TOKEN" \
                --start-from-head \
                --output json 2>/dev/null || echo '{"events":[],"nextForwardToken":""}')
        fi

        echo "$LOG_OUTPUT" | python3 -c "
import json, sys
for e in json.load(sys.stdin).get('events', []):
    print(e['message'])" 2>/dev/null || true

        NEW_TOKEN=$(echo "$LOG_OUTPUT" | python3 -c "
import json, sys
print(json.load(sys.stdin).get('nextForwardToken', ''))" 2>/dev/null || true)

        if [ "$NEW_TOKEN" = "$NEXT_TOKEN" ] && { [ "$CUR_STATUS" = "SUCCEEDED" ] || [ "$CUR_STATUS" = "FAILED" ]; }; then
            break
        fi
        NEXT_TOKEN="$NEW_TOKEN"
        sleep 5
    done
fi

echo ""
FINAL_STATUS=$($AWS batch describe-jobs --jobs "$JOB_ID" \
    --query "jobs[0].status" --output text)
echo "=== Job $FINAL_STATUS ==="

if [ "$FINAL_STATUS" = "FAILED" ]; then
    REASON=$($AWS batch describe-jobs --jobs "$JOB_ID" \
        --query "jobs[0].container.reason" --output text 2>/dev/null || echo "unknown")
    echo "Failure reason: $REASON"
    exit 1
fi
