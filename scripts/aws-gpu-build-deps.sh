#!/bin/bash
#
# One-time: build and push the GPU deps-only Docker image + GPU runtime base image.
# Run this when Cargo.toml dependencies change (rare).
#
# After this, `aws-gpu-bench.sh` only recompiles src/ (~30-60s) instead of all deps.
#
set -euo pipefail

PROFILE="AdministratorAccess-270497617191"
REGION="us-east-1"
ACCOUNT_ID="270497617191"
PREFIX="ldsc-bench"
ECR_REPO="${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com/${PREFIX}"
DNS_SERVER="10.64.0.1"

AWS="aws --profile $PROFILE --region $REGION --no-cli-pager"

# CUDA version for cudarc symbol loading (must be <= target driver's CUDA version).
# 12040 = CUDA 12.4 (AWS default), 11080 = CUDA 11.8 (older HPC)
CUDA_TARGET_VERSION="${CUDA_TARGET_VERSION:-12040}"

echo "=== Building GPU deps image (glibc + gpu + mimalloc, CUDA ${CUDA_TARGET_VERSION}) ==="
echo "This compiles all dependencies and caches the target/ directory."
echo "Only needs to run when Cargo.toml/Cargo.lock change."
echo ""

# Step 1: Compile deps only (dummy main.rs), cache target/ in image
docker rm -f ldsc-gpu-deps-build 2>/dev/null || true
docker run --dns "$DNS_SERVER" --name ldsc-gpu-deps-build \
    -v "$(pwd)/Cargo.toml:/build/Cargo.toml:ro" \
    -v "$(pwd)/Cargo.lock:/build/Cargo.lock:ro" \
    -v "$(pwd)/.cargo:/build/.cargo:ro" \
    -w /build \
    -e CARGO_INCREMENTAL=0 \
    -e CUDARC_CUDA_VERSION="$CUDA_TARGET_VERSION" \
    rust:1-bookworm \
    sh -c '
        mkdir -p src &&
        echo "fn main() {}" > src/main.rs &&
        echo "fn main() {}" > build.rs &&
        cargo build --release --features gpu,mimalloc 2>&1 &&
        echo "=== GPU deps compiled ==="
    '
docker commit ldsc-gpu-deps-build "${PREFIX}-gpu-deps:latest"
docker rm -f ldsc-gpu-deps-build
echo "GPU deps image built: ${PREFIX}-gpu-deps:latest"

# Step 2: Build GPU runtime base image (CUDA runtime + hyperfine + awscli)
echo ""
echo "=== Building GPU runtime base image ==="
docker rm -f ldsc-gpu-runtime-build 2>/dev/null || true
docker run --dns "$DNS_SERVER" --name ldsc-gpu-runtime-build \
    ubuntu:22.04 \
    bash -c '
        apt-get update &&
        apt-get install -y --no-install-recommends ca-certificates curl unzip gnupg &&
        curl -fsSL https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/3bf863cc.pub \
            | gpg --dearmor -o /usr/share/keyrings/cuda-archive-keyring.gpg &&
        echo "deb [signed-by=/usr/share/keyrings/cuda-archive-keyring.gpg] https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/ /" \
            > /etc/apt/sources.list.d/cuda.list &&
        apt-get update &&
        apt-get install -y --no-install-recommends cuda-nvrtc-12-6 &&
        rm -rf /var/lib/apt/lists/* &&
        curl -sL "https://github.com/sharkdp/hyperfine/releases/download/v1.19.0/hyperfine_1.19.0_amd64.deb" \
            -o /tmp/hf.deb && dpkg -i /tmp/hf.deb && rm /tmp/hf.deb &&
        curl -sL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o /tmp/awscli.zip &&
        cd /tmp && unzip -q awscli.zip && ./aws/install &&
        rm -rf /tmp/aws /tmp/awscli.zip
    '
docker cp scripts/bench-entrypoint.sh ldsc-gpu-runtime-build:/entrypoint.sh
docker commit ldsc-gpu-runtime-build "${PREFIX}-gpu-runtime:latest"
docker rm -f ldsc-gpu-runtime-build

# Make entrypoint executable
docker run --name ldsc-gpu-runtime-chmod "${PREFIX}-gpu-runtime:latest" chmod +x /entrypoint.sh
docker commit --change 'ENTRYPOINT ["/entrypoint.sh"]' ldsc-gpu-runtime-chmod "${PREFIX}-gpu-runtime:latest"
docker rm -f ldsc-gpu-runtime-chmod

echo "GPU runtime base image built: ${PREFIX}-gpu-runtime:latest"

# Step 3: Push both to ECR
echo ""
echo "=== Pushing to ECR ==="
$AWS ecr get-login-password | docker login --username AWS --password-stdin \
    "${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com" 2>/dev/null

# Ensure repos exist (reuse main ldsc-bench repo with gpu- prefix tags)
for repo in "${PREFIX}-gpu-deps" "${PREFIX}-gpu-runtime"; do
    $AWS ecr describe-repositories --repository-names "$repo" >/dev/null 2>&1 || \
        $AWS ecr create-repository --repository-name "$repo" >/dev/null
done

docker tag "${PREFIX}-gpu-deps:latest" "${ECR_REPO}-gpu-deps:latest"
docker push "${ECR_REPO}-gpu-deps:latest" 2>&1 | tail -3

docker tag "${PREFIX}-gpu-runtime:latest" "${ECR_REPO}-gpu-runtime:latest"
docker push "${ECR_REPO}-gpu-runtime:latest" 2>&1 | tail -3

echo ""
echo "Done! Now use 'scripts/aws-gpu-bench.sh' for fast iterations."
