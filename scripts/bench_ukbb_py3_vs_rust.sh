#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="$ROOT/data/UKBB.ALL.ldscore"
OUT_DIR="$ROOT/perf/ukbb"
PY3_PROJECT="$ROOT/ldsc_py3"
RUST_BIN="$ROOT/target/release/ldsc"

POP="${POP:-EUR}"
LDMS="${LDMS:-none}"   # none | 8 | 25
ID_SCHEME="${ID_SCHEME:-rsid}"  # rsid | chrpos
SUMSTATS="${SUMSTATS:-}"
SUMSTATS2="${SUMSTATS2:-}"
OVERLAP_ANNOT="${OVERLAP_ANNOT:-0}"
LOG_LEVEL="${LOG_LEVEL:-}"

if [[ -z "$SUMSTATS" ]]; then
  echo "SUMSTATS is required (path to .sumstats.gz)." >&2
  exit 1
fi

if [[ ! -d "$DATA_DIR" ]]; then
  if [[ -f "$ROOT/data/UKBB.ALL.ldscore.tar.gz" ]]; then
    tar -xzf "$ROOT/data/UKBB.ALL.ldscore.tar.gz" -C "$ROOT/data"
  else
    echo "Missing $ROOT/data/UKBB.ALL.ldscore.tar.gz" >&2
    exit 1
  fi
fi

suffix=""
if [[ "$ID_SCHEME" == "rsid" ]]; then
  suffix="rsid."
fi

case "$LDMS" in
  none)
    REF_LD="$DATA_DIR/UKBB.${POP}.${suffix}l2.ldscore.gz"
    ;;
  8)
    REF_LD="$DATA_DIR/UKBB.${POP}.8LDMS.${suffix}l2.ldscore.gz"
    ;;
  25)
    REF_LD="$DATA_DIR/UKBB.${POP}.25LDMS.${suffix}l2.ldscore.gz"
    ;;
  *)
    echo "LDMS must be one of: none | 8 | 25" >&2
    exit 1
    ;;
 esac

W_LD="$DATA_DIR/UKBB.${POP}.${suffix}l2.ldscore.gz"

# Python LDSC expects --ref-ld/--w-ld as a prefix (it appends .l2.ldscore.gz)
REF_LD_PY="$REF_LD"
W_LD_PY="$W_LD"
if [[ "$REF_LD_PY" == *.l2.ldscore.gz ]]; then
  REF_LD_PY="${REF_LD_PY%.l2.ldscore.gz}"
fi
if [[ "$W_LD_PY" == *.l2.ldscore.gz ]]; then
  W_LD_PY="${W_LD_PY%.l2.ldscore.gz}"
fi

# If M sidecar exists for single-file LD scores, pass --M to Rust to match Python.
M_FLAG=()
M_SUFFIX=".l2.M_5_50"
if [[ "${NOT_M_5_50:-0}" == "1" ]]; then
  M_SUFFIX=".l2.M"
fi
M_PATH="${REF_LD}"
if [[ "$M_PATH" == *.l2.ldscore.gz ]]; then
  M_PATH="${M_PATH%.l2.ldscore.gz}${M_SUFFIX}"
else
  M_PATH="${M_PATH}${M_SUFFIX}"
fi
if [[ -f "$M_PATH" ]]; then
  M_LINE="$(head -n 1 "$M_PATH")"
  if [[ "$M_LINE" != *$'\t'* && "$M_LINE" != *' '* ]]; then
    M_VAL="$(echo "$M_LINE" | tr -d '\t ')"
    if [[ -n "$M_VAL" ]]; then
      M_FLAG=(--M "$M_VAL")
    fi
  fi
fi

for f in "$REF_LD" "$W_LD"; do
  if [[ ! -f "$f" ]]; then
    echo "Missing required LD score file: $f" >&2
    exit 1
  fi
done

mkdir -p "$OUT_DIR"

if ! command -v uv >/dev/null 2>&1; then
  echo "uv is required but not found in PATH." >&2
  exit 1
fi

# Python runner (uv PEP 723 script)
PY3_PYTHON="${PY3_PYTHON:-3.9}"
PY_WRAPPER="$ROOT/scripts/_ldsc_wrapper.py"
if [[ ! -f "$PY_WRAPPER" ]]; then
  echo "Missing Python wrapper: $PY_WRAPPER" >&2
  exit 1
fi

# Require prebuilt Rust release binary.
if [[ ! -x "$RUST_BIN" ]]; then
  echo "Missing release binary: $RUST_BIN" >&2
  echo "Build it first with: cargo build --release" >&2
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

PY3_CMD=(env LDSC_SCRIPT="$PY3_PROJECT/ldsc.py" uv run --python "$PY3_PYTHON" --script "$PY_WRAPPER")
RUST_CMD=("$RUST_BIN")

EXTRA_FLAGS=()
if [[ "$OVERLAP_ANNOT" == "1" ]]; then
  EXTRA_FLAGS+=(--overlap-annot)
fi
if [[ -n "${TWO_STEP:-}" ]]; then
  EXTRA_FLAGS+=(--two-step "$TWO_STEP")
fi
if [[ -n "${CHISQ_MAX:-}" ]]; then
  EXTRA_FLAGS+=(--chisq-max "$CHISQ_MAX")
fi

PY3_H2_CMD=("${PY3_CMD[@]}" --h2 "$SUMSTATS" --ref-ld "$REF_LD_PY" --w-ld "$W_LD_PY" --out "$OUT_DIR/py3_h2")
RUST_H2_CMD=("${RUST_CMD[@]}" h2 --h2 "$SUMSTATS" --ref-ld "$REF_LD" --w-ld "$W_LD" --out "$OUT_DIR/rust_h2")
if [[ -n "$LOG_LEVEL" ]]; then
  PY3_H2_CMD+=(--log-level "$LOG_LEVEL")
  RUST_H2_CMD+=(--log-level "$LOG_LEVEL")
fi
if [[ "${#M_FLAG[@]}" -gt 0 ]]; then
  RUST_H2_CMD+=("${M_FLAG[@]}")
fi
if [[ "${#EXTRA_FLAGS[@]}" -gt 0 ]]; then
  PY3_H2_CMD+=("${EXTRA_FLAGS[@]}")
  RUST_H2_CMD+=("${EXTRA_FLAGS[@]}")
fi

run_timed "py3_h2" "${PY3_H2_CMD[@]}"
if [[ -n "$LOG_LEVEL" ]]; then
  run_timed "rust_h2" env RUST_LOG="ldsc=$LOG_LEVEL" "${RUST_H2_CMD[@]}"
else
  run_timed "rust_h2" "${RUST_H2_CMD[@]}"
fi

if [[ -n "$SUMSTATS2" ]]; then
  PY3_RG_CMD=("${PY3_CMD[@]}" --rg "$SUMSTATS,$SUMSTATS2" --ref-ld "$REF_LD_PY" --w-ld "$W_LD_PY" --out "$OUT_DIR/py3_rg")
  RUST_RG_CMD=("${RUST_CMD[@]}" rg --rg "$SUMSTATS,$SUMSTATS2" --ref-ld "$REF_LD" --w-ld "$W_LD" --out "$OUT_DIR/rust_rg")
  if [[ -n "$LOG_LEVEL" ]]; then
    PY3_RG_CMD+=(--log-level "$LOG_LEVEL")
    RUST_RG_CMD+=(--log-level "$LOG_LEVEL")
  fi
  if [[ "${#M_FLAG[@]}" -gt 0 ]]; then
    RUST_RG_CMD+=("${M_FLAG[@]}")
  fi
  if [[ "${#EXTRA_FLAGS[@]}" -gt 0 ]]; then
    PY3_RG_CMD+=("${EXTRA_FLAGS[@]}")
    RUST_RG_CMD+=("${EXTRA_FLAGS[@]}")
  fi
  run_timed "py3_rg" "${PY3_RG_CMD[@]}"
  if [[ -n "$LOG_LEVEL" ]]; then
    run_timed "rust_rg" env RUST_LOG="ldsc=$LOG_LEVEL" "${RUST_RG_CMD[@]}"
  else
    run_timed "rust_rg" "${RUST_RG_CMD[@]}"
  fi
fi

echo "Done. Timing files are in $OUT_DIR/*.time"
