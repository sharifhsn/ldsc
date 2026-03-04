#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="$ROOT/perf/l2_threads"
RUST_BIN="$ROOT/target/release/ldsc"
BFILE="${BFILE:-$ROOT/data/1000G.EUR.QC}"
EXTRACT="${EXTRACT:-}"
EXTRACT_N="${EXTRACT_N:-50000}"
THREADS_LIST="${THREADS_LIST:-1 2 4 8 16 32}"
LD_WIND_CM="${LD_WIND_CM:-1}"

mkdir -p "$OUT_DIR"

if [[ ! -x "$RUST_BIN" ]]; then
  echo "Missing release binary: $RUST_BIN" >&2
  echo "Build it first with: cargo build --release" >&2
  exit 1
fi

if [[ -z "$EXTRACT" && "$EXTRACT_N" != "0" ]]; then
  EXTRACT_FILE="$OUT_DIR/extract_${EXTRACT_N}.snps"
  awk -v n="$EXTRACT_N" 'NR<=n {print $2} NR==n {exit}' "${BFILE}.bim" > "$EXTRACT_FILE"
  EXTRACT="$EXTRACT_FILE"
fi

TIME_BIN=""
if command -v /usr/bin/time >/dev/null 2>&1; then
  TIME_BIN="/usr/bin/time -p"
fi

run_timed() {
  local label="$1"; shift
  echo "==> $label"
  if [[ -n "$TIME_BIN" ]]; then
    # shellcheck disable=SC2086
    $TIME_BIN "$@" 2>"$OUT_DIR/${label}.time"
  else
    time "$@" 2>"$OUT_DIR/${label}.time"
  fi
}

for t in $THREADS_LIST; do
  OUT_PREFIX="$OUT_DIR/rust_l2_t${t}"
  CMD=("$RUST_BIN" --rayon-threads "$t" l2 --bfile "$BFILE" --out "$OUT_PREFIX" --ld-wind-cm "$LD_WIND_CM")
  if [[ -n "$EXTRACT" ]]; then
    CMD+=(--extract "$EXTRACT")
  fi
  run_timed "rust_l2_t${t}" "${CMD[@]}"
done

echo "Done. Timing files are in $OUT_DIR/*.time"
