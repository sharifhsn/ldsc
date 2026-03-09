#!/usr/bin/env bash
# build_pgo.sh — Profile-Guided Optimization build for ldsc.
#
# LLVM PGO uses runtime profiles to optimize branch prediction, function
# layout, and inlining decisions. Requires llvm-profdata (from the llvm
# package). On Arch: `sudo pacman -S llvm`.
#
# Steps:
#   1. Build instrumented binary (profile-generate)
#   2. Run benchmark on representative data (generates .profraw files)
#   3. Merge profiles with llvm-profdata
#   4. Rebuild with merged profile (profile-use)
#
# The PGO binary is written to target/release/ldsc (same path, replaces normal build).
# Save a copy of the non-PGO binary first if you want to compare.
#
# Usage:
#   bash scripts/build_pgo.sh
#   BFILE=/path/to/custom bash scripts/build_pgo.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BFILE="${BFILE:-$ROOT/data/1000G_phase3_common_norel}"
PGO_DIR="${PGO_DIR:-/tmp/ldsc-pgo}"
MERGED="$PGO_DIR/merged.profdata"

if [[ ! -f "${BFILE}.bim" ]]; then
    echo "ERROR: BIM file not found: ${BFILE}.bim"
    echo "Set BFILE= to a path with .bed/.bim/.fam files."
    exit 1
fi

if ! command -v llvm-profdata &>/dev/null; then
    echo "ERROR: llvm-profdata not found."
    echo "On Arch Linux: sudo pacman -S llvm"
    echo "On Ubuntu/Debian: sudo apt install llvm"
    exit 1
fi

mkdir -p "$PGO_DIR"
rm -f "$PGO_DIR"/*.profraw "$PGO_DIR"/*.profdata

echo "=== Step 1: Build instrumented binary ==="
RUSTFLAGS="-Cprofile-generate=$PGO_DIR" \
    cargo build --release 2>&1

echo ""
echo "=== Step 2: Run profile workload (full genome, --ld-wind-kb 1000) ==="
echo "This takes ~80s — same as a normal benchmark run."
"$ROOT/target/release/ldsc" l2 \
    --bfile "$BFILE" \
    --ld-wind-kb 1000 \
    --out "$PGO_DIR/pgo_train" \
    2>&1

echo ""
echo "=== Step 3: Merge profile data ==="
llvm-profdata merge -output="$MERGED" "$PGO_DIR"/*.profraw
echo "Merged profile: $MERGED"

echo ""
echo "=== Step 4: Build optimized binary ==="
RUSTFLAGS="-Cprofile-use=$MERGED -Cllvm-args=-pgo-warn-missing-function" \
    cargo build --release 2>&1

echo ""
echo "=== PGO build complete ==="
echo "Binary: $ROOT/target/release/ldsc"
echo ""
echo "To benchmark:"
echo "  hyperfine --warmup 1 --runs 3 \\"
echo "    '$ROOT/target/release/ldsc l2 --bfile $BFILE --ld-wind-kb 1000 --out /tmp/bench_pgo'"
