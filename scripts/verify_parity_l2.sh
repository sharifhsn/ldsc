#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="$ROOT/perf/parity_l2"
RUST_BIN="$ROOT/target/release/ldsc"
BFILE="${BFILE:-$ROOT/data/1000G_phase3_common_norel}"
EXTRACT_N="${EXTRACT_N:-50000}"
LD_WIND_KB="${LD_WIND_KB:-1000}"
CHUNK_SIZE="${CHUNK_SIZE:-200}"

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

EXTRACT_FILE=""
if [[ "$EXTRACT_N" != "0" ]]; then
  EXTRACT_FILE="$OUT_DIR/extract.snps"
  awk -v n="$EXTRACT_N" 'NR<=n {print $2} NR==n {exit}' "${BFILE}.bim" > "$EXTRACT_FILE"
fi

PY_OUT="$OUT_DIR/py3_l2"
RS_OUT="$OUT_DIR/rust_l2"
rm -f "${PY_OUT}.l2.ldscore.gz" "${PY_OUT}.l2.M" "${PY_OUT}.l2.M_5_50"
rm -f "${RS_OUT}"*.l2.ldscore.gz "${RS_OUT}"*.l2.M "${RS_OUT}"*.l2.M_5_50

PY_CMD=(uv run --project "$ROOT/ldsc_py" --with "bitarray==2" python "$ROOT/ldsc_py/ldsc.py" --l2 --bfile "$BFILE" --out "$PY_OUT" --ld-wind-kb "$LD_WIND_KB" --chunk-size "$CHUNK_SIZE")
RS_CMD=("$RUST_BIN" l2 --bfile "$BFILE" --out "$RS_OUT" --ld-wind-kb "$LD_WIND_KB" --chunk-size "$CHUNK_SIZE")
if [[ -n "$EXTRACT_FILE" ]]; then
  PY_CMD+=(--extract "$EXTRACT_FILE")
  RS_CMD+=(--extract "$EXTRACT_FILE")
fi

echo "Running Python LDSC..."
"${PY_CMD[@]}" >"$OUT_DIR/py3_l2.stdout" 2>&1
echo "Running Rust LDSC..."
"${RS_CMD[@]}" >"$OUT_DIR/rust_l2.stdout" 2>&1

CMP_SCRIPT="$OUT_DIR/_compare_l2.py"
cat >"$CMP_SCRIPT" <<'PY'
# /// script
# requires-python = ">=3.9"
# dependencies = []
# ///
import gzip
import sys
from pathlib import Path

py_path = Path(sys.argv[1])
rust_paths = [Path(p) for p in sys.argv[2:]]

def read_ldscore(path):
    rows = {}
    with gzip.open(path, 'rt') as f:
        header = f.readline().strip().split('\t')
        l2_idx = header.index('L2')
        snp_idx = header.index('SNP')
        for line in f:
            parts = line.rstrip('\n').split('\t')
            snp = parts[snp_idx]
            l2 = float(parts[l2_idx])
            rows[snp] = l2
    return rows

py = read_ldscore(py_path)
rs = {}
for p in rust_paths:
    rs.update(read_ldscore(p))

if len(py) != len(rs):
    raise SystemExit(f"row count mismatch: py={len(py)} rust={len(rs)}")

max_diff = 0.0
worst = None
for snp, py_v in py.items():
    rs_v = rs.get(snp)
    if rs_v is None:
        raise SystemExit(f"missing snp in rust: {snp}")
    diff = abs(py_v - rs_v)
    if diff > max_diff:
        max_diff = diff
        worst = (snp, py_v, rs_v)

print(f"max_abs_diff={max_diff:.6g}")
if worst is not None:
    print(f"worst={worst[0]} py={worst[1]:.6g} rust={worst[2]:.6g}")

if max_diff > 1e-3:
    raise SystemExit("diff exceeds tolerance 1e-3")
print("OK")
PY

RUST_FILES=()
if [[ -f "${RS_OUT}.l2.ldscore.gz" ]]; then
  RUST_FILES=("${RS_OUT}.l2.ldscore.gz")
else
  RUST_FILES=("${RS_OUT}"*.l2.ldscore.gz)
  if [[ "${#RUST_FILES[@]}" -eq 1 && "${RUST_FILES[0]}" == "${RS_OUT}*.l2.ldscore.gz" ]]; then
    echo "No Rust l2 output files found at ${RS_OUT}*.l2.ldscore.gz" >&2
    exit 1
  fi
fi

uv run --script "$CMP_SCRIPT" "${PY_OUT}.l2.ldscore.gz" "${RUST_FILES[@]}"

echo "L2 parity check passed."
