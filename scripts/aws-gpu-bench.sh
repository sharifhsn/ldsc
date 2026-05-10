#!/bin/bash
#
# Build, push, and run ldsc GPU benchmark on AWS Batch (g5.2xlarge A10G).
# Streams CloudWatch logs back to terminal in real time.
#
# Usage:
#   ./scripts/aws-gpu-bench.sh                    # default: gpu-biobank-bench sweep
#   ./scripts/aws-gpu-bench.sh --skip-build       # reuse existing image
#   ./scripts/aws-gpu-bench.sh --command "..."     # custom command
#   ./scripts/aws-gpu-bench.sh --runs 5            # override hyperfine runs
#
# Prerequisites:
#   1. Run scripts/aws-gpu-setup.sh once (creates GPU compute env + queue + job def)
#   2. Run scripts/aws-gpu-build-deps.sh once (builds deps + runtime base images)
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
LOG_GROUPS=("/aws/batch/${PREFIX}" "/aws/batch/job")

# ── Parse arguments ───────────────────────────────────────────────────────────
SKIP_BUILD=false
BENCH_RUNS=3
BENCH_WARMUP=1
CUSTOM_CMD=""
VCPUS=8
MEMORY=28672
CUDA_VERSION=12040  # CUDA 12.4 — safe default for AWS GPU AMIs
JOB_NAME="gpu-bench-$(date +%Y%m%d-%H%M%S)"

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-build)  SKIP_BUILD=true; shift ;;
        --runs)        BENCH_RUNS="$2"; shift 2 ;;
        --warmup)      BENCH_WARMUP="$2"; shift 2 ;;
        --command)     CUSTOM_CMD="$2"; shift 2 ;;
        --vcpus)       VCPUS="$2"; shift 2 ;;
        --memory)      MEMORY="$2"; shift 2 ;;
        --cuda-version) CUDA_VERSION="$2"; shift 2 ;;
        --name)        JOB_NAME="$2"; shift 2 ;;
        *)             echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ── Build & Push ──────────────────────────────────────────────────────────────
GPU_DEPS_IMAGE="${PREFIX}-gpu-deps:latest"
GPU_RUNTIME_IMAGE="${PREFIX}-gpu-runtime:latest"

if [ "$SKIP_BUILD" = false ]; then
    # Check that GPU deps image exists
    if ! docker image inspect "$GPU_DEPS_IMAGE" >/dev/null 2>&1; then
        echo "ERROR: GPU deps image '$GPU_DEPS_IMAGE' not found." >&2
        echo "Run 'scripts/aws-gpu-build-deps.sh' first (one-time setup)." >&2
        exit 1
    fi
    if ! docker image inspect "$GPU_RUNTIME_IMAGE" >/dev/null 2>&1; then
        echo "ERROR: GPU runtime image '$GPU_RUNTIME_IMAGE' not found." >&2
        echo "Run 'scripts/aws-gpu-build-deps.sh' first (one-time setup)." >&2
        exit 1
    fi

    echo "=== Compiling src/ only with GPU feature (CUDA ${CUDA_VERSION}, deps cached) ==="
    docker rm -f ldsc-gpu-build 2>/dev/null || true
    docker run --dns "$DNS_SERVER" --name ldsc-gpu-build \
        -v "$(pwd)/src:/build/src:ro" \
        -v "$(pwd)/Cargo.toml:/build/Cargo.toml:ro" \
        -v "$(pwd)/Cargo.lock:/build/Cargo.lock:ro" \
        -v "$(pwd)/build.rs:/build/build.rs:ro" \
        -v "$(pwd)/.cargo:/build/.cargo:ro" \
        -w /build \
        -e CARGO_INCREMENTAL=0 \
        -e CUDARC_CUDA_VERSION="$CUDA_VERSION" \
        "$GPU_DEPS_IMAGE" \
        sh -c "rm -f target/release/ldsc target/release/deps/ldsc-* && \
               cargo build --release --features gpu,mimalloc && \
               strip target/release/ldsc"
    docker cp ldsc-gpu-build:/build/target/release/ldsc target/release/ldsc.gpu
    docker rm -f ldsc-gpu-build

    echo ""
    echo "=== Assembling GPU runtime image ==="
    chmod +x target/release/ldsc.gpu
    docker rm -f ldsc-gpu-assemble 2>/dev/null || true
    docker run --entrypoint true --name ldsc-gpu-assemble "$GPU_RUNTIME_IMAGE"
    docker cp target/release/ldsc.gpu ldsc-gpu-assemble:/usr/local/bin/ldsc
    docker cp scripts/bench-entrypoint.sh ldsc-gpu-assemble:/entrypoint.sh
    docker cp scripts/gpu-biobank-bench.sh ldsc-gpu-assemble:/usr/local/bin/gpu-biobank-bench
    docker cp scripts/biobank-bench.sh ldsc-gpu-assemble:/usr/local/bin/biobank-bench
    docker commit \
        --change 'ENTRYPOINT ["/entrypoint.sh"]' \
        ldsc-gpu-assemble "${PREFIX}:gpu-latest"
    docker rm -f ldsc-gpu-assemble

    echo "Image built: ${PREFIX}:gpu-latest"

    echo ""
    echo "=== Pushing to ECR ==="
    $AWS ecr get-login-password | docker login --username AWS --password-stdin \
        "${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com" 2>/dev/null

    docker tag "${PREFIX}:gpu-latest" "${ECR_REPO}:gpu-latest"
    docker push "${ECR_REPO}:gpu-latest" 2>&1 | tail -5
    echo "Push complete."
else
    echo "=== Skipping build (--skip-build) ==="
fi

# ── Submit Batch Job ──────────────────────────────────────────────────────────
echo ""
echo "=== Submitting GPU Batch job: $JOB_NAME ==="

OVERRIDES=$(cat << EOF
{
  "resourceRequirements": [
    {"type": "VCPU", "value": "${VCPUS}"},
    {"type": "MEMORY", "value": "${MEMORY}"},
    {"type": "GPU", "value": "1"}
  ],
  "environment": [
    {"name": "BENCH_RUNS", "value": "${BENCH_RUNS}"},
    {"name": "BENCH_WARMUP", "value": "${BENCH_WARMUP}"},
    {"name": "BENCH_DATASET", "value": "biobank_50k"},
    {"name": "S3_DATA_BUCKET", "value": "${S3_BUCKET}"},
    {"name": "GPU_BENCH", "value": "1"}
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

JOB_ID=$($AWS batch submit-job \
    --job-name "$JOB_NAME" \
    --job-queue "${PREFIX}-gpu-queue-v2" \
    --job-definition "${PREFIX}-gpu-job" \
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
        LOG_GROUP="/aws/batch/${PREFIX}"
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

        # Print new events
        echo "$LOG_OUTPUT" | python3 -c "
import json, sys
data = json.load(sys.stdin)
for evt in data.get('events', []):
    print(evt['message'])
" 2>/dev/null

        NEW_TOKEN=$(echo "$LOG_OUTPUT" | python3 -c "
import json, sys
print(json.load(sys.stdin).get('nextForwardToken', ''))" 2>/dev/null || true)
        if [ -n "$NEW_TOKEN" ] && [ "$NEW_TOKEN" != "$NEXT_TOKEN" ]; then
            NEXT_TOKEN="$NEW_TOKEN"
        fi

        if [ "$CUR_STATUS" = "SUCCEEDED" ] || [ "$CUR_STATUS" = "FAILED" ]; then
            # Final flush
            sleep 2
            if [ -n "$NEXT_TOKEN" ]; then
                $AWS logs get-log-events \
                    --log-group-name "$LOG_GROUP" \
                    --log-stream-name "$LOG_STREAM" \
                    --next-token "$NEXT_TOKEN" \
                    --start-from-head \
                    --output json 2>/dev/null | python3 -c "
import json, sys
data = json.load(sys.stdin)
for evt in data.get('events', []):
    print(evt['message'])
" 2>/dev/null || true
            fi
            echo ""
            echo "=== Job ${CUR_STATUS} ==="
            break
        fi
        sleep 3
    done
fi
