#!/bin/bash
#
# Build, push, and run ldsc benchmark on AWS Batch.
# Streams CloudWatch logs back to terminal in real time.
#
# Usage:
#   ./scripts/aws-bench.sh                          # default benchmark
#   ./scripts/aws-bench.sh --runs 20                # override run count
#   ./scripts/aws-bench.sh --skip-build             # reuse existing image
#   ./scripts/aws-bench.sh --command "hyperfine ..." # custom command
#   ./scripts/aws-bench.sh --vcpus 8 --memory 4096  # override resources
#   ./scripts/aws-bench.sh --dataset biobank_50k    # use 50K-individual synthetic data
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
LOG_GROUPS=("/aws/batch/${PREFIX}" "/aws/batch/job")  # try custom first, then default

# ── Parse arguments ───────────────────────────────────────────────────────────
SKIP_BUILD=false
BENCH_RUNS=10
BENCH_WARMUP=2
CUSTOM_CMD=""
VCPUS=16
MEMORY=8192
BENCH_DATASET="1000G"
JOB_NAME="bench-$(date +%Y%m%d-%H%M%S)"

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-build)  SKIP_BUILD=true; shift ;;
        --runs)        BENCH_RUNS="$2"; shift 2 ;;
        --warmup)      BENCH_WARMUP="$2"; shift 2 ;;
        --command)     CUSTOM_CMD="$2"; shift 2 ;;
        --vcpus)       VCPUS="$2"; shift 2 ;;
        --memory)      MEMORY="$2"; shift 2 ;;
        --dataset)     BENCH_DATASET="$2"; shift 2 ;;
        --name)        JOB_NAME="$2"; shift 2 ;;
        *)             echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Auto-increase memory for biobank dataset (working memory ~2GB, plus page cache headroom)
if [ "$BENCH_DATASET" = "biobank_50k" ] && [ "$MEMORY" -lt 16384 ]; then
    MEMORY=16384
    echo "NOTE: Auto-increased memory to ${MEMORY}MB for biobank_50k dataset"
fi

# ── Build & Push ──────────────────────────────────────────────────────────────
#
# Uses pre-built deps image for fast incremental builds (~30s for src/ only).
# Run `scripts/aws-build-deps.sh` first (one-time, or when deps change).
#
# Docker daemon DNS is broken (8.8.8.8 unreachable); use --dns workaround.
#
DEPS_IMAGE="${PREFIX}-deps:latest"
RUNTIME_IMAGE="${PREFIX}-runtime:latest"

if [ "$SKIP_BUILD" = false ]; then
    # Check that deps image exists
    if ! docker image inspect "$DEPS_IMAGE" >/dev/null 2>&1; then
        echo "ERROR: Deps image '$DEPS_IMAGE' not found." >&2
        echo "Run 'scripts/aws-build-deps.sh' first (one-time setup)." >&2
        exit 1
    fi
    if ! docker image inspect "$RUNTIME_IMAGE" >/dev/null 2>&1; then
        echo "ERROR: Runtime image '$RUNTIME_IMAGE' not found." >&2
        echo "Run 'scripts/aws-build-deps.sh' first (one-time setup)." >&2
        exit 1
    fi

    echo "=== Compiling src/ only (deps cached) ==="
    docker rm -f ldsc-cargo-build 2>/dev/null || true
    docker run --dns "$DNS_SERVER" --name ldsc-cargo-build \
        -v "$(pwd)/src:/build/src:ro" \
        -v "$(pwd)/Cargo.toml:/build/Cargo.toml:ro" \
        -v "$(pwd)/Cargo.lock:/build/Cargo.lock:ro" \
        -v "$(pwd)/.cargo/config.toml:/build/.cargo/config.toml:ro" \
        -w /build \
        -e CARGO_INCREMENTAL=0 \
        "$DEPS_IMAGE" \
        sh -c "rm -f target/x86_64-unknown-linux-musl/release/ldsc target/x86_64-unknown-linux-musl/release/deps/ldsc-* && \
               cargo build --release --features mimalloc --target x86_64-unknown-linux-musl && \
               strip target/x86_64-unknown-linux-musl/release/ldsc"
    docker cp ldsc-cargo-build:/build/target/x86_64-unknown-linux-musl/release/ldsc target/release/ldsc.musl
    docker rm -f ldsc-cargo-build

    echo ""
    echo "=== Assembling runtime image ==="
    chmod +x target/release/ldsc.musl
    docker rm -f ldsc-bench-assemble 2>/dev/null || true
    docker run --entrypoint true --name ldsc-bench-assemble "$RUNTIME_IMAGE"
    docker cp target/release/ldsc.musl ldsc-bench-assemble:/usr/local/bin/ldsc
    docker cp scripts/bench-entrypoint.sh ldsc-bench-assemble:/entrypoint.sh
    docker cp scripts/biobank-bench.sh ldsc-bench-assemble:/usr/local/bin/biobank-bench
    docker commit \
        --change 'ENTRYPOINT ["/entrypoint.sh"]' \
        ldsc-bench-assemble "${PREFIX}:latest"
    docker rm -f ldsc-bench-assemble

    echo "Image built: ${PREFIX}:latest"

    echo ""
    echo "=== Pushing to ECR ==="
    $AWS ecr get-login-password | docker login --username AWS --password-stdin \
        "${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com" 2>/dev/null

    docker tag "${PREFIX}:latest" "${ECR_REPO}:latest"
    docker push "${ECR_REPO}:latest" 2>&1 | tail -5
    echo "Push complete."
else
    echo "=== Skipping build (--skip-build) ==="
fi

# ── Submit Batch Job ──────────────────────────────────────────────────────────
echo ""
echo "=== Submitting Batch job: $JOB_NAME ==="

# Build container overrides
OVERRIDES=$(cat << EOF
{
  "resourceRequirements": [
    {"type": "VCPU", "value": "${VCPUS}"},
    {"type": "MEMORY", "value": "${MEMORY}"}
  ],
  "environment": [
    {"name": "BENCH_RUNS", "value": "${BENCH_RUNS}"},
    {"name": "BENCH_WARMUP", "value": "${BENCH_WARMUP}"},
    {"name": "BENCH_DATASET", "value": "${BENCH_DATASET}"},
    {"name": "S3_DATA_BUCKET", "value": "${S3_BUCKET}"}
  ]
}
EOF
)

# Add custom command override if specified
if [ -n "$CUSTOM_CMD" ]; then
    OVERRIDES=$(echo "$OVERRIDES" | CUSTOM_CMD="$CUSTOM_CMD" python3 -c "
import json, sys, shlex, os
o = json.load(sys.stdin)
o['command'] = shlex.split(os.environ['CUSTOM_CMD'])
json.dump(o, sys.stdout)
")
fi

# Use biobank queue (50GB EBS) for biobank dataset, default queue otherwise
if [ "$BENCH_DATASET" = "biobank_50k" ]; then
    JOB_QUEUE="${PREFIX}-biobank-queue"
else
    JOB_QUEUE="${PREFIX}-queue"
fi

JOB_ID=$($AWS batch submit-job \
    --job-name "$JOB_NAME" \
    --job-queue "$JOB_QUEUE" \
    --job-definition "${PREFIX}-job" \
    --container-overrides "$OVERRIDES" \
    --query "jobId" --output text)

echo "Job ID: $JOB_ID"
echo "Console: https://${REGION}.console.aws.amazon.com/batch/home?region=${REGION}#jobs/detail/${JOB_ID}"

# ── Poll job status until RUNNING ────────────────────────────────────────────
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

# ── Stream CloudWatch logs ────────────────────────────────────────────────────
if [ -n "$LOG_STREAM" ]; then
    echo ""
    echo "=== Streaming logs (${LOG_STREAM}) ==="
    echo ""

    # Auto-detect which log group contains this stream
    LOG_GROUP=""
    for g in "${LOG_GROUPS[@]}"; do
        if $AWS logs get-log-events \
            --log-group-name "$g" \
            --log-stream-name "$LOG_STREAM" \
            --start-from-head --limit 1 \
            --output json >/dev/null 2>&1; then
            LOG_GROUP="$g"
            break
        fi
    done
    if [ -z "$LOG_GROUP" ]; then
        echo "WARNING: Could not find log stream in any log group (tried: ${LOG_GROUPS[*]})"
        LOG_GROUP="/aws/batch/${PREFIX}"  # fall back
    else
        echo "Using log group: $LOG_GROUP"
    fi

    NEXT_TOKEN=""
    while true; do
        CUR_STATUS=$($AWS batch describe-jobs --jobs "$JOB_ID" \
            --query "jobs[0].status" --output text)

        if [ -z "$NEXT_TOKEN" ]; then
            LOG_OUTPUT=$($AWS logs get-log-events \
                --log-group-name "$LOG_GROUP" \
                --log-stream-name "$LOG_STREAM" \
                --start-from-head \
                --output json 2>/dev/null || echo '{"events":[],"nextForwardToken":""}')
        else
            LOG_OUTPUT=$($AWS logs get-log-events \
                --log-group-name "$LOG_GROUP" \
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
        sleep 3
    done

    echo ""
    echo "=== Job $STATUS ==="
else
    echo "WARNING: Could not find log stream for job $JOB_ID"
    echo "Check the AWS Batch console for details."
fi

# ── Print final status ────────────────────────────────────────────────────────
FINAL_STATUS=$($AWS batch describe-jobs --jobs "$JOB_ID" \
    --query "jobs[0].status" --output text)
echo "Final status: $FINAL_STATUS"

if [ "$FINAL_STATUS" = "FAILED" ]; then
    REASON=$($AWS batch describe-jobs --jobs "$JOB_ID" \
        --query "jobs[0].container.reason" --output text 2>/dev/null || echo "unknown")
    echo "Failure reason: $REASON"
    exit 1
fi
