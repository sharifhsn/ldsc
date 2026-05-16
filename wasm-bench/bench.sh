#!/usr/bin/env bash
# Build + run the GEMM microbench under wasmtime.
# Usage: ./bench.sh [extra wasmtime args]
set -euo pipefail

cd "$(dirname "$0")"

# Pick a rustup toolchain that has wasm32-wasip1.
RUSTUP_BIN="${RUSTUP_BIN:-/Users/sharif/.rustup/toolchains/stable-aarch64-apple-darwin/bin}"
export PATH="$RUSTUP_BIN:$PATH"

cargo build --release --target wasm32-wasip1 -q

wasmtime run \
    -W simd=y \
    -W relaxed-simd=y \
    -W bulk-memory=y \
    -W max-memory-size=4294967296 \
    target/wasm32-wasip1/release/ldsc-wasm-bench.wasm "$@"
