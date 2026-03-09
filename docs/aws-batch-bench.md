# AWS Batch Benchmarking Pipeline

A self-contained guide to running reproducible benchmarks on AWS Batch with EC2 Spot instances. Supports both Rust and Python (uv-managed) applications.

Drop this doc, the referenced scripts, and Dockerfiles into any repo to get a thermally stable, isolated benchmarking environment for ~$0.02/run.

---

## 1. Why AWS Batch?

Local workstations hit thermal throttling under sustained load, corrupting timing distributions. HPC clusters have queue wait times and shared-filesystem contention. AWS Batch on Spot gives:

- **Consistent compute**: Dedicated vCPUs on bare-metal-class instances (no neighbor noise).
- **Reproducibility**: Same AMI, same instance type, same region every time.
- **Cost**: A 10-run hyperfine benchmark on `c7a.4xlarge` (16 vCPU) takes ~15 minutes wall time = **~$0.02** at Spot rates.
- **No idle costs**: Instances provision on job submit and terminate on completion.

## 2. Architecture

```
┌─────────────┐    ┌──────────┐    ┌───────────────────┐    ┌────────────┐
│ Local machine│───▶│ ECR      │───▶│ AWS Batch         │───▶│ CloudWatch │
│ (build+push) │    │ (images) │    │ (EC2 Spot c7a.4xl)│    │ (logs)     │
└─────────────┘    └──────────┘    └───────────────────┘    └────────────┘
       │                                    │                       │
       │              S3 bucket ◀───────────┘                       │
       │           (benchmark data)                                 │
       └──────────── streams logs in real time ◀────────────────────┘
```

### Components

| Resource | Purpose |
|---|---|
| **ECR Repository** | Stores benchmark Docker images |
| **S3 Bucket** | Holds input data files (downloaded by container at startup) |
| **Batch Compute Environment** | EC2 Spot pool (`SPOT_CAPACITY_OPTIMIZED`), `c7a` family |
| **Batch Job Queue** | Routes jobs to the compute environment |
| **Batch Job Definition** | Container config (image, default vCPUs/memory) |
| **CloudWatch Log Group** | Captures stdout/stderr from containers |
| **IAM Roles** | Batch execution role (ECR pull + S3 read + CW write) |

## 3. One-Time AWS Setup

All commands use the AWS CLI. Set these variables first:

```bash
REGION="us-east-1"
ACCOUNT_ID="$(aws sts get-caller-identity --query Account --output text)"
PREFIX="my-bench"  # change per project
ECR_REPO="${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com/${PREFIX}"
S3_BUCKET="${PREFIX}-data-${ACCOUNT_ID}"
```

### 3a. ECR Repository

```bash
aws ecr create-repository --repository-name "$PREFIX" --region "$REGION"
```

### 3b. S3 Data Bucket

```bash
aws s3 mb "s3://${S3_BUCKET}" --region "$REGION"
# Upload your benchmark data:
aws s3 cp my_data/ "s3://${S3_BUCKET}/" --recursive
```

### 3c. IAM Roles

Create two roles:

**Batch Service Role** (AWS-managed, usually auto-created):
```bash
# Usually exists already as AWSBatchServiceRole
```

**ECS Task Execution Role** (for ECR pull + S3 read + CloudWatch):
```bash
cat > /tmp/trust.json << 'EOF'
{
  "Version": "2012-10-17",
  "Statement": [{
    "Effect": "Allow",
    "Principal": {"Service": "ecs-tasks.amazonaws.com"},
    "Action": "sts:AssumeRole"
  }]
}
EOF

aws iam create-role \
    --role-name "${PREFIX}-task-role" \
    --assume-role-policy-document file:///tmp/trust.json

# Attach policies
aws iam attach-role-policy --role-name "${PREFIX}-task-role" \
    --policy-arn arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy
aws iam attach-role-policy --role-name "${PREFIX}-task-role" \
    --policy-arn arn:aws:iam::aws:policy/AmazonS3ReadOnlyAccess
```

### 3d. Batch Compute Environment

```bash
cat > /tmp/ce.json << EOF
{
  "type": "MANAGED",
  "allocationStrategy": "SPOT_CAPACITY_OPTIMIZED",
  "minvCpus": 0,
  "maxvCpus": 64,
  "desiredvCpus": 0,
  "instanceTypes": ["c7a.4xlarge", "c7a.2xlarge", "c6a.4xlarge"],
  "subnets": ["<YOUR-SUBNET-ID>"],
  "securityGroupIds": ["<YOUR-SG-ID>"],
  "instanceRole": "ecsInstanceRole",
  "type": "SPOT",
  "spotIamFleetRole": "arn:aws:iam::${ACCOUNT_ID}:role/AmazonEC2SpotFleetRole"
}
EOF

aws batch create-compute-environment \
    --compute-environment-name "${PREFIX}-ce" \
    --type MANAGED \
    --compute-resources file:///tmp/ce.json \
    --region "$REGION"
```

> **Subnet/SG**: Use your default VPC's public subnet and default security group. The container needs outbound internet for S3 access.

### 3e. Job Queue + Job Definition

```bash
# Job queue
aws batch create-job-queue \
    --job-queue-name "${PREFIX}-queue" \
    --priority 1 \
    --compute-environment-order "order=1,computeEnvironment=${PREFIX}-ce" \
    --region "$REGION"

# Job definition
aws batch register-job-definition \
    --job-definition-name "${PREFIX}-job" \
    --type container \
    --container-properties '{
        "image": "'${ECR_REPO}':latest",
        "resourceRequirements": [
            {"type": "VCPU", "value": "16"},
            {"type": "MEMORY", "value": "8192"}
        ],
        "jobRoleArn": "arn:aws:iam::'${ACCOUNT_ID}':role/'${PREFIX}'-task-role",
        "executionRoleArn": "arn:aws:iam::'${ACCOUNT_ID}':role/'${PREFIX}'-task-role",
        "logConfiguration": {
            "logDriver": "awslogs",
            "options": {
                "awslogs-group": "/aws/batch/'${PREFIX}'",
                "awslogs-region": "'${REGION}'",
                "awslogs-stream-prefix": "bench"
            }
        }
    }' \
    --region "$REGION"

# Create the log group
aws logs create-log-group \
    --log-group-name "/aws/batch/${PREFIX}" \
    --region "$REGION"
```

### 3f. CloudWatch Log Group

```bash
aws logs create-log-group \
    --log-group-name "/aws/batch/${PREFIX}" \
    --region "$REGION"
```

## 4. Dockerfiles

### 4a. Rust Application

Multi-stage build. Pin the Debian version in the build stage to avoid glibc mismatch with the runtime stage.

```dockerfile
# ── build ────────────────────────────────────────────────────────────────
FROM rust:1-bookworm AS builder
WORKDIR /app
COPY Cargo.toml Cargo.lock ./
COPY src/ src/
RUN cargo build --release && strip target/release/my-app

# ── runtime ──────────────────────────────────────────────────────────────
FROM debian:bookworm-slim AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl unzip && \
    rm -rf /var/lib/apt/lists/*

# hyperfine
ARG HYPERFINE_VERSION=1.19.0
RUN curl -sL "https://github.com/sharkdp/hyperfine/releases/download/v${HYPERFINE_VERSION}/hyperfine_${HYPERFINE_VERSION}_amd64.deb" \
    -o /tmp/hf.deb && dpkg -i /tmp/hf.deb && rm /tmp/hf.deb

# AWS CLI (for S3 data download at container start)
RUN curl -sL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o /tmp/awscli.zip && \
    cd /tmp && unzip -q awscli.zip && ./aws/install && \
    rm -rf /tmp/aws /tmp/awscli.zip

COPY --from=builder /app/target/release/my-app /usr/local/bin/my-app
COPY scripts/bench-entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
```

> **CRITICAL: Pin `rust:1-bookworm`, not `rust:1`.** The untagged `rust:1` image silently
> tracks the latest Debian release. When Debian trixie became stable, `rust:1` jumped from
> bookworm (glibc 2.36) to trixie (glibc 2.39). Binaries built on trixie reference
> `GLIBC_2.39` symbols (e.g. tokio's `pidfd_getpid`) that don't exist on bookworm-slim
> runtime containers, causing `version 'GLIBC_2.39' not found` crashes. Always pin the
> Debian codename in both build and runtime stages.

#### CentOS 7 Compatibility (glibc 2.17)

If your target HPC runs CentOS 7 or RHEL 7, build inside a CentOS 7 container so the binary links against glibc 2.17:

```dockerfile
FROM centos:7 AS builder
# CentOS 7 is EOL — point repos at vault
RUN sed -i 's|mirrorlist=|#mirrorlist=|g' /etc/yum.repos.d/CentOS-*.repo && \
    sed -i 's|#baseurl=http://mirror.centos.org|baseurl=http://vault.centos.org|g' /etc/yum.repos.d/CentOS-*.repo
RUN yum install -y gcc-c++ make perl-IPC-Cmd && yum clean all -y
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable
ENV PATH="/root/.cargo/bin:$PATH"
WORKDIR /app
COPY Cargo.toml Cargo.lock ./
COPY src/ src/
RUN cargo build --release && strip target/release/my-app
```

#### musl Static Build (runs on any Linux)

For a fully static binary with zero glibc dependency:

```dockerfile
FROM rust:1-alpine AS builder
RUN apk add --no-cache musl-dev gcc g++ make
WORKDIR /app
COPY Cargo.toml Cargo.lock ./
COPY src/ src/
RUN cargo build --release --target x86_64-unknown-linux-musl && \
    strip target/x86_64-unknown-linux-musl/release/my-app
```

> **Performance note**: musl's allocator is slower than glibc's for multi-threaded
> workloads. If you see a regression, add `jemallocator` as a dependency and set it as
> the global allocator. Benchmark both before deciding.

### 4b. Python Application (uv-managed)

```dockerfile
# ── build ────────────────────────────────────────────────────────────────
FROM python:3.12-slim-bookworm AS builder
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

WORKDIR /app

# Install dependencies first (cached unless lock/toml change)
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --locked --no-install-project --no-editable

# Then install the project itself
COPY . /app
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-editable

# ── runtime ──────────────────────────────────────────────────────────────
FROM python:3.12-slim-bookworm AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl unzip && \
    rm -rf /var/lib/apt/lists/*

# hyperfine
ARG HYPERFINE_VERSION=1.19.0
RUN curl -sL "https://github.com/sharkdp/hyperfine/releases/download/v${HYPERFINE_VERSION}/hyperfine_${HYPERFINE_VERSION}_amd64.deb" \
    -o /tmp/hf.deb && dpkg -i /tmp/hf.deb && rm /tmp/hf.deb

# AWS CLI
RUN curl -sL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o /tmp/awscli.zip && \
    cd /tmp && unzip -q awscli.zip && ./aws/install && \
    rm -rf /tmp/aws /tmp/awscli.zip

# Copy only the venv from the builder (no uv needed at runtime)
COPY --from=builder /app/.venv /app/.venv
ENV PATH="/app/.venv/bin:$PATH"

COPY scripts/bench-entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
```

> Add `.venv` and `__pycache__` to `.dockerignore` to keep the build context clean.

## 5. Entrypoint Script

Generic entrypoint that downloads data from S3, prints system info, and runs hyperfine. Customize the benchmark command for your application.

```bash
#!/bin/bash
set -euo pipefail

S3_BUCKET="${S3_DATA_BUCKET:-my-bench-data-000000000000}"
RUNS="${BENCH_RUNS:-10}"
WARMUP="${BENCH_WARMUP:-2}"

# ── Download data ────────────────────────────────────────────────────────
echo "=== Downloading data from s3://${S3_BUCKET}/ ==="
mkdir -p /data
aws s3 cp "s3://${S3_BUCKET}/" /data/ --recursive --quiet
echo "Download complete."
echo ""

# ── System info ──────────────────────────────────────────────────────────
echo "=== Benchmark Environment ==="
echo "CPU: $(nproc) vCPUs"
echo "CPU Model: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs)"
echo "Memory: $(free -h 2>/dev/null | awk '/Mem:/{print $2}' || echo 'unknown')"
echo "Kernel: $(uname -r)"
echo "Date: $(date -u)"
echo ""

# ── Benchmark ────────────────────────────────────────────────────────────
if [ $# -eq 0 ]; then
    echo "=== Running benchmark (${RUNS} runs, ${WARMUP} warmup) ==="
    hyperfine \
        --warmup "$WARMUP" \
        --runs "$RUNS" \
        --export-json /tmp/results.json \
        "my-app --input /data/input.dat --output /tmp/out"

    echo ""
    echo "=== RESULTS JSON ==="
    cat /tmp/results.json
else
    echo "=== Running custom command ==="
    exec "$@"
fi
```

## 6. Deployment Script

The `aws-bench.sh` script handles the full cycle: build, push, submit, and stream logs.

```bash
#!/bin/bash
set -euo pipefail

REGION="us-east-1"
ACCOUNT_ID="$(aws sts get-caller-identity --query Account --output text)"
PREFIX="my-bench"
ECR_REPO="${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com/${PREFIX}"
S3_BUCKET="${PREFIX}-data-${ACCOUNT_ID}"

# ── Parse arguments ──────────────────────────────────────────────────────
SKIP_BUILD=false
BENCH_RUNS=10
BENCH_WARMUP=2
CUSTOM_CMD=""
VCPUS=16
MEMORY=8192
JOB_NAME="bench-$(date +%Y%m%d-%H%M%S)"

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-build)  SKIP_BUILD=true; shift ;;
        --runs)        BENCH_RUNS="$2"; shift 2 ;;
        --warmup)      BENCH_WARMUP="$2"; shift 2 ;;
        --command)     CUSTOM_CMD="$2"; shift 2 ;;
        --vcpus)       VCPUS="$2"; shift 2 ;;
        --memory)      MEMORY="$2"; shift 2 ;;
        --name)        JOB_NAME="$2"; shift 2 ;;
        *)             echo "Unknown: $1"; exit 1 ;;
    esac
done

# ── Build & Push ─────────────────────────────────────────────────────────
if [ "$SKIP_BUILD" = false ]; then
    echo "=== Building Docker image ==="
    docker buildx build --network=host -f Dockerfile.bench -t "${PREFIX}:latest" .

    echo "=== Pushing to ECR ==="
    aws ecr get-login-password --region "$REGION" | \
        docker login --username AWS --password-stdin \
        "${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com"
    docker tag "${PREFIX}:latest" "${ECR_REPO}:latest"
    docker push "${ECR_REPO}:latest" 2>&1 | tail -5
fi

# ── Submit Job ───────────────────────────────────────────────────────────
echo "=== Submitting: $JOB_NAME ==="

OVERRIDES=$(cat << EOF
{
  "resourceRequirements": [
    {"type": "VCPU", "value": "${VCPUS}"},
    {"type": "MEMORY", "value": "${MEMORY}"}
  ],
  "environment": [
    {"name": "BENCH_RUNS", "value": "${BENCH_RUNS}"},
    {"name": "BENCH_WARMUP", "value": "${BENCH_WARMUP}"},
    {"name": "S3_DATA_BUCKET", "value": "${S3_BUCKET}"}
  ]
}
EOF
)

if [ -n "$CUSTOM_CMD" ]; then
    OVERRIDES=$(echo "$OVERRIDES" | python3 -c "
import json, sys, shlex
o = json.load(sys.stdin)
o['command'] = shlex.split('''$CUSTOM_CMD''')
json.dump(o, sys.stdout)")
fi

JOB_ID=$(aws batch submit-job \
    --job-name "$JOB_NAME" \
    --job-queue "${PREFIX}-queue" \
    --job-definition "${PREFIX}-job" \
    --container-overrides "$OVERRIDES" \
    --region "$REGION" \
    --query "jobId" --output text)

echo "Job ID: $JOB_ID"
echo "Console: https://${REGION}.console.aws.amazon.com/batch/home?region=${REGION}#jobs/detail/${JOB_ID}"

# ── Poll until running ───────────────────────────────────────────────────
echo ""
echo "=== Waiting for job ==="
LOG_STREAM=""
PREV_STATUS=""

while true; do
    JOB_JSON=$(aws batch describe-jobs --jobs "$JOB_ID" --region "$REGION")
    STATUS=$(echo "$JOB_JSON" | python3 -c "import json,sys; print(json.load(sys.stdin)['jobs'][0]['status'])")

    if [ "$STATUS" != "$PREV_STATUS" ]; then
        echo "$(date +%H:%M:%S) $STATUS"
        PREV_STATUS="$STATUS"
    fi

    case "$STATUS" in
        SUCCEEDED|FAILED)
            LOG_STREAM=$(echo "$JOB_JSON" | python3 -c "
import json, sys
print(json.load(sys.stdin)['jobs'][0].get('container',{}).get('logStreamName',''))" 2>/dev/null || true)
            break ;;
        RUNNING)
            LOG_STREAM=$(echo "$JOB_JSON" | python3 -c "
import json, sys
print(json.load(sys.stdin)['jobs'][0].get('container',{}).get('logStreamName',''))" 2>/dev/null || true)
            [ -n "$LOG_STREAM" ] && break ;;
    esac
    sleep 5
done

# ── Stream CloudWatch logs ───────────────────────────────────────────────
if [ -n "$LOG_STREAM" ]; then
    echo ""
    echo "=== Logs ==="
    NEXT_TOKEN=""
    while true; do
        CUR_STATUS=$(aws batch describe-jobs --jobs "$JOB_ID" --region "$REGION" \
            --query "jobs[0].status" --output text)

        if [ -z "$NEXT_TOKEN" ]; then
            LOG_OUTPUT=$(aws logs get-log-events \
                --log-group-name "/aws/batch/${PREFIX}" \
                --log-stream-name "$LOG_STREAM" \
                --start-from-head --region "$REGION" \
                --output json 2>/dev/null || echo '{"events":[],"nextForwardToken":""}')
        else
            LOG_OUTPUT=$(aws logs get-log-events \
                --log-group-name "/aws/batch/${PREFIX}" \
                --log-stream-name "$LOG_STREAM" \
                --next-token "$NEXT_TOKEN" \
                --start-from-head --region "$REGION" \
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
fi

FINAL=$(aws batch describe-jobs --jobs "$JOB_ID" --region "$REGION" \
    --query "jobs[0].status" --output text)
echo "Final: $FINAL"
[ "$FINAL" = "FAILED" ] && exit 1
```

## 7. Usage

```bash
# Full build + benchmark (10 runs, 2 warmup)
./scripts/aws-bench.sh

# Skip rebuild, reuse existing image
./scripts/aws-bench.sh --skip-build

# More runs for tighter confidence intervals
./scripts/aws-bench.sh --runs 20 --warmup 5

# Custom command (e.g. A/B test two binaries)
./scripts/aws-bench.sh --command "hyperfine --warmup 2 --runs 10 'my-app-v1 ...' 'my-app-v2 ...'"

# Override compute resources
./scripts/aws-bench.sh --vcpus 8 --memory 4096
```

## 8. Variance Mitigation

Cloud VMs introduce microarchitectural noise that bare-metal doesn't have. To get reliable numbers:

| Technique | How |
|---|---|
| **Warmup runs** | `--warmup 2` (minimum). First runs page in data and warm caches. |
| **High run count** | `--runs 10` gives σ/μ < 1% on `c7a` Spot. Use 20 for sub-0.5%. |
| **Dedicated instances** | `c7a` family has consistent EPYC cores. Avoid burstable (`t3`). |
| **Spot stability** | `SPOT_CAPACITY_OPTIMIZED` avoids interruption-prone pools. |
| **Pin vCPUs** | Request exact vCPUs matching the instance (16 for c7a.4xlarge) to avoid co-scheduling. |

## 9. Cost Estimate

| Scenario | Instance | Spot $/hr | Typical Duration | Cost |
|---|---|---|---|---|
| Quick bench (10 runs, small data) | c7a.4xlarge | ~$0.08 | 15 min | **$0.02** |
| Full bench (20 runs, large data) | c7a.4xlarge | ~$0.08 | 30 min | **$0.04** |
| Heavy compute (96 vCPU) | c7a.24xlarge | ~$1.40 | 30 min | **$0.70** |

Spot instances auto-terminate when the job finishes. **No idle costs.** S3 and ECR storage are negligible (pennies/month).

## 10. Docker DNS Workaround

Docker BuildKit has a known DNS bug ([moby/buildkit#5009](https://github.com/moby/buildkit/issues/5009)) where containers can't resolve hostnames. Two workarounds:

**Option A**: Use `--network=host` with BuildKit:
```bash
docker buildx build --network=host -f Dockerfile.bench -t my-bench .
```

**Option B**: Use `docker run --dns` for manual builds:
```bash
docker run --dns 10.64.0.1 --name build-step debian:bookworm-slim \
    bash -c 'apt-get update && apt-get install -y curl ...'
docker cp my-binary build-step:/usr/local/bin/
docker commit build-step my-bench:latest
```

Replace `10.64.0.1` with your host's DNS server (`cat /etc/resolv.conf`).

## 11. Troubleshooting

| Problem | Cause | Fix |
|---|---|---|
| `GLIBC_2.xx not found` | Build image has newer glibc than runtime | Pin both stages to same Debian codename (e.g. `bookworm`) |
| Job stuck in `RUNNABLE` | No Spot capacity or subnet misconfigured | Check Batch CE status; add more instance types to the pool |
| S3 download fails | Missing IAM permissions or wrong bucket | Verify task role has `s3:GetObject` on the bucket |
| High variance (>2% CV) | Burstable instance or too few runs | Use `c7a`/`c6a` family; increase `--runs` |
| Docker DNS failure | BuildKit DNS bug | Use `--network=host` or `--dns` flag (see section 10) |
| Container OOM killed | Memory limit too low | Increase `--memory` in job overrides |
