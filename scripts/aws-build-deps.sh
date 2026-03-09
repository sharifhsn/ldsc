#!/bin/bash
#
# One-time: build and push the deps-only Docker image + runtime base image.
# Run this when Cargo.toml dependencies change (rare).
#
# After this, `aws-bench.sh` only recompiles src/ (~30s) instead of all deps (~11min).
#
set -euo pipefail

PROFILE="AdministratorAccess-270497617191"
REGION="us-east-1"
ACCOUNT_ID="270497617191"
PREFIX="ldsc-bench"
ECR_REPO="${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com/${PREFIX}"
DNS_SERVER="10.64.0.1"

AWS="aws --profile $PROFILE --region $REGION --no-cli-pager"

echo "=== Building deps image (musl + mimalloc) ==="
echo "This compiles all dependencies and caches the target/ directory."
echo "Only needs to run when Cargo.toml/Cargo.lock change."
echo ""

# Step 1: Compile deps only (dummy main.rs), cache target/ in image
docker rm -f ldsc-deps-build 2>/dev/null || true
docker run --dns "$DNS_SERVER" --name ldsc-deps-build \
    -v "$(pwd)/Cargo.toml:/build/Cargo.toml:ro" \
    -v "$(pwd)/Cargo.lock:/build/Cargo.lock:ro" \
    -w /build \
    -e CARGO_INCREMENTAL=0 \
    rust:1-alpine \
    sh -c '
        apk add --no-cache musl-dev &&
        mkdir -p src &&
        echo "fn main() {}" > src/main.rs &&
        cargo build --release --features mimalloc --target x86_64-unknown-linux-musl 2>&1 &&
        echo "=== Deps compiled ==="
    '
docker commit ldsc-deps-build "${PREFIX}-deps:latest"
docker rm -f ldsc-deps-build
echo "Deps image built: ${PREFIX}-deps:latest"

# Step 2: Build runtime base image (hyperfine + awscli, no ldsc binary yet)
echo ""
echo "=== Building runtime base image ==="
docker rm -f ldsc-runtime-build 2>/dev/null || true
docker run --dns "$DNS_SERVER" --name ldsc-runtime-build \
    debian:bookworm-slim \
    bash -c '
        apt-get update &&
        apt-get install -y --no-install-recommends ca-certificates curl unzip &&
        rm -rf /var/lib/apt/lists/* &&
        curl -sL "https://github.com/sharkdp/hyperfine/releases/download/v1.19.0/hyperfine_1.19.0_amd64.deb" \
            -o /tmp/hf.deb && dpkg -i /tmp/hf.deb && rm /tmp/hf.deb &&
        curl -sL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o /tmp/awscli.zip &&
        cd /tmp && unzip -q awscli.zip && ./aws/install &&
        rm -rf /tmp/aws /tmp/awscli.zip
    '
# Add entrypoint script
docker cp scripts/bench-entrypoint.sh ldsc-runtime-build:/entrypoint.sh
docker commit ldsc-runtime-build "${PREFIX}-runtime:latest"
docker rm -f ldsc-runtime-build

# Make entrypoint executable
docker run --name ldsc-runtime-chmod "${PREFIX}-runtime:latest" chmod +x /entrypoint.sh
docker commit --change 'ENTRYPOINT ["/entrypoint.sh"]' ldsc-runtime-chmod "${PREFIX}-runtime:latest"
docker rm -f ldsc-runtime-chmod

echo "Runtime base image built: ${PREFIX}-runtime:latest"

# Step 3: Push both to ECR
echo ""
echo "=== Pushing to ECR ==="
$AWS ecr get-login-password | docker login --username AWS --password-stdin \
    "${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com" 2>/dev/null

# Ensure repos exist
for repo in "${PREFIX}-deps" "${PREFIX}-runtime"; do
    $AWS ecr describe-repositories --repository-names "$repo" >/dev/null 2>&1 || \
        $AWS ecr create-repository --repository-name "$repo" >/dev/null
done

docker tag "${PREFIX}-deps:latest" "${ECR_REPO}-deps:latest"
docker push "${ECR_REPO}-deps:latest" 2>&1 | tail -3

docker tag "${PREFIX}-runtime:latest" "${ECR_REPO}-runtime:latest"
docker push "${ECR_REPO}-runtime:latest" 2>&1 | tail -3

echo ""
echo "Done! Now use 'scripts/aws-bench.sh' for fast iterations (only recompiles src/)."
