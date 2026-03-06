#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${OUT_DIR:-$ROOT/perf/tiny_l2_check}"
PY_PROJECT="${PY_PROJECT:-$ROOT/ldsc_py}"
RUST_BIN="${RUST_BIN:-$ROOT/target/release/ldsc}"
BFILE="${BFILE:-$ROOT/data/1000G_phase3_common_norel}"
EXTRACT_N="${EXTRACT_N:-2000}"
LD_WIND_CM="${LD_WIND_CM:-1}"
CHUNK_SIZE="${CHUNK_SIZE:-50}"
TOL="${TOL:-1e-3}"
PY_WITH="${PY_WITH:-bitarray==2}"

mkdir -p "$OUT_DIR"

if ! command -v uv >/dev/null 2>&1; then
  echo "uv is required but not found in PATH." >&2
  exit 1
fi

if [[ ! -x "$RUST_BIN" ]]; then
  echo "Missing Rust binary: $RUST_BIN" >&2
  echo "Build it first (example): cargo build --release" >&2
  exit 1
fi

if [[ ! -f "$PY_PROJECT/ldsc.py" ]]; then
  echo "Missing Python LDSC entrypoint: $PY_PROJECT/ldsc.py" >&2
  exit 1
fi

for ext in bed bim fam; do
  if [[ ! -f "${BFILE}.${ext}" ]]; then
    echo "Missing input file: ${BFILE}.${ext}" >&2
    exit 1
  fi
done

EXTRACT_FILE="$OUT_DIR/extract_${EXTRACT_N}.snps"
awk -v n="$EXTRACT_N" 'NR<=n {print $2} NR==n {exit}' "${BFILE}.bim" > "$EXTRACT_FILE"
if [[ "$(wc -l < "$EXTRACT_FILE")" -ne "$EXTRACT_N" ]]; then
  echo "Failed to create extract list with $EXTRACT_N SNPs." >&2
  exit 1
fi

PY_OUT="$OUT_DIR/py_l2"
RS_OUT="$OUT_DIR/rust_l2"
rm -f "${PY_OUT}"* "${RS_OUT}"*

PY_CMD=(uv run --project "$PY_PROJECT")
if [[ -n "$PY_WITH" ]]; then
  PY_CMD+=(--with "$PY_WITH")
fi
PY_CMD+=(python "$PY_PROJECT/ldsc.py" --l2 --bfile "$BFILE" --out "$PY_OUT" --ld-wind-cm "$LD_WIND_CM" --chunk-size "$CHUNK_SIZE" --extract "$EXTRACT_FILE" --yes-really)
RS_CMD=("$RUST_BIN" l2 --bfile "$BFILE" --out "$RS_OUT" --ld-wind-cm "$LD_WIND_CM" --chunk-size "$CHUNK_SIZE" --extract "$EXTRACT_FILE" --yes-really)

run_timed() {
  local label="$1"
  shift
  local t0 t1 elapsed
  t0="$(date +%s.%N)"
  "$@" >"$OUT_DIR/${label}.stdout" 2>"$OUT_DIR/${label}.stderr"
  t1="$(date +%s.%N)"
  elapsed="$(awk -v a="$t0" -v b="$t1" 'BEGIN { printf "%.6f", (b - a) }')"
  printf "%s\n" "$elapsed" >"$OUT_DIR/${label}.time"
}

run_timed py "${PY_CMD[@]}"
run_timed rust "${RS_CMD[@]}"

RUST_FILES=()
if [[ -f "${RS_OUT}.l2.ldscore.gz" ]]; then
  RUST_FILES=("${RS_OUT}.l2.ldscore.gz")
else
  RUST_FILES=("${RS_OUT}"*.l2.ldscore.gz)
  if [[ "${#RUST_FILES[@]}" -eq 1 && "${RUST_FILES[0]}" == "${RS_OUT}*.l2.ldscore.gz" ]]; then
    echo "No Rust l2 outputs found." >&2
    exit 1
  fi
fi

python3 - "$TOL" "${PY_OUT}.l2.ldscore.gz" "${RUST_FILES[@]}" <<'PY'
import gzip
import sys
from pathlib import Path

tol = float(sys.argv[1])
py_path = Path(sys.argv[2])
rust_paths = [Path(p) for p in sys.argv[3:]]

def read_ldscore(path: Path):
    rows = {}
    with gzip.open(path, "rt") as f:
        header = f.readline().strip().split("\t")
        snp_idx = header.index("SNP")
        l2_idx = header.index("L2")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            rows[parts[snp_idx]] = float(parts[l2_idx])
    return rows

py_rows = read_ldscore(py_path)
rs_rows = {}
for p in rust_paths:
    rs_rows.update(read_ldscore(p))

if len(py_rows) != len(rs_rows):
    raise SystemExit(f"row count mismatch: py={len(py_rows)} rust={len(rs_rows)}")

max_diff = 0.0
worst = None
for snp, py_v in py_rows.items():
    rs_v = rs_rows.get(snp)
    if rs_v is None:
        raise SystemExit(f"missing snp in rust output: {snp}")
    d = abs(py_v - rs_v)
    if d > max_diff:
        max_diff = d
        worst = (snp, py_v, rs_v)

print(f"max_abs_diff={max_diff:.6g}")
if worst is not None:
    print(f"worst={worst[0]} py={worst[1]:.6g} rust={worst[2]:.6g}")

if max_diff > tol:
    raise SystemExit(f"parity FAIL: max_abs_diff={max_diff:.6g} > tol={tol}")

print("parity OK")
PY

PY_REAL="$(cat "$OUT_DIR/py.time" 2>/dev/null || true)"
RS_REAL="$(cat "$OUT_DIR/rust.time" 2>/dev/null || true)"

echo ""
echo "Tiny l2 check complete"
echo "  extract_n: $EXTRACT_N"
echo "  py_real_s: ${PY_REAL:-n/a}"
echo "  rust_real_s: ${RS_REAL:-n/a}"
echo "  output_dir: $OUT_DIR"
