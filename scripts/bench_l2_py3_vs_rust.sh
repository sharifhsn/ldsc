#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="$ROOT/perf/l2"
PY3_PROJECT="$ROOT/ldsc_py3"
RUST_BIN="$ROOT/target/release/ldsc"
PY3_PYTHON="${PY3_PYTHON:-3.9}"
BFILE="${BFILE:-$ROOT/data/1000G.EUR.QC}"
EXTRACT="${EXTRACT:-}"
EXTRACT_N="${EXTRACT_N:-500000}"
LOG_LEVEL="${LOG_LEVEL:-}"
LD_WIND_CM="${LD_WIND_CM:-1}"

mkdir -p "$OUT_DIR"

if ! command -v uv >/dev/null 2>&1; then
  echo "uv is required but not found in PATH." >&2
  exit 1
fi

if [[ ! -x "$RUST_BIN" ]]; then
  echo "Missing release binary: $RUST_BIN" >&2
  echo "Build it first with: cargo build --release" >&2
  exit 1
fi

PY_WRAPPER="$ROOT/scripts/_ldsc_wrapper.py"
if [[ ! -f "$PY_WRAPPER" ]]; then
  echo "Missing Python wrapper: $PY_WRAPPER" >&2
  exit 1
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

if [[ -z "$EXTRACT" && "$EXTRACT_N" != "0" ]]; then
  EXTRACT_FILE="$OUT_DIR/extract_${EXTRACT_N}.snps"
  awk -v n="$EXTRACT_N" 'NR<=n {print $2} NR==n {exit}' "${BFILE}.bim" > "$EXTRACT_FILE"
  EXTRACT="$EXTRACT_FILE"
fi

PY_CMD=(env LDSC_SCRIPT="$PY3_PROJECT/ldsc.py" uv run --python "$PY3_PYTHON" --script "$PY_WRAPPER" --l2 --bfile "$BFILE" --out "$OUT_DIR/py3_l2" --ld-wind-cm "$LD_WIND_CM")
RS_CMD=("$RUST_BIN" l2 --bfile "$BFILE" --out "$OUT_DIR/rust_l2" --ld-wind-cm "$LD_WIND_CM")
if [[ -n "$EXTRACT" ]]; then
  PY_CMD+=(--extract "$EXTRACT")
  RS_CMD+=(--extract "$EXTRACT")
fi
if [[ -n "$LOG_LEVEL" ]]; then
  PY_CMD+=(--log-level "$LOG_LEVEL")
  RS_CMD+=(--log-level "$LOG_LEVEL")
fi

run_timed "py3_l2" "${PY_CMD[@]}"
if [[ -n "$LOG_LEVEL" ]]; then
  run_timed "rust_l2" env RUST_LOG="ldsc=$LOG_LEVEL" "${RS_CMD[@]}"
else
  run_timed "rust_l2" "${RS_CMD[@]}"
fi

echo "Done. Timing files are in $OUT_DIR/*.time"
