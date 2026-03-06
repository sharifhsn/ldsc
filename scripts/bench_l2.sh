#!/usr/bin/env bash
# bench_l2.sh — L2 parity + performance benchmark using a pre-built BED file.
#
# Unlike check_l2_tiny_py_vs_rust.sh this script does NOT use --extract.
# Both Python and Rust read the BED file directly, so timing reflects real
# compute cost rather than whole-file I/O overhead.
#
# Uses --ld-wind-kb (not --ld-wind-cm) because the 1000G BED files have the
# CM column zeroed out; --ld-wind-cm 1 with all-zero CM would make every SNP
# on every chromosome fall in every other's window, giving nonsensical scores.
# --ld-wind-kb 1000 uses BP coordinates and naturally respects chr boundaries.
#
# Default dataset: data/bench_5k  (5k SNPs, fast smoke benchmark)
#   plink1.9 --bfile data/bench_200k --thin-count 5000 --seed 42 \
#            --make-bed --out data/bench_5k
# For a heavier run:
#   BFILE=$ROOT/data/bench_200k bash scripts/bench_l2.sh
#   BFILE=$ROOT/data/1000G_phase3_common_norel bash scripts/bench_l2.sh
#
# Usage:
#   bash scripts/bench_l2.sh                         # default 5k dataset
#   BFILE=$ROOT/data/bench_200k bash scripts/bench_l2.sh   # 200k dataset
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BFILE="${BFILE:-$ROOT/data/bench_5k}"
PY_PROJECT="${PY_PROJECT:-$ROOT/ldsc_py}"
RUST_BIN="${RUST_BIN:-$ROOT/target/release/ldsc}"
OUT_DIR="${OUT_DIR:-$ROOT/perf/bench_l2}"
LD_WIND_KB="${LD_WIND_KB:-1000}"
CHUNK_SIZE="${CHUNK_SIZE:-50}"
TOL="${TOL:-1e-3}"
PY_WITH="${PY_WITH:-bitarray==2}"

mkdir -p "$OUT_DIR"

# ---------------------------------------------------------------------------
# Preflight checks
# ---------------------------------------------------------------------------
if ! command -v uv >/dev/null 2>&1; then
  echo "uv is required but not found in PATH." >&2
  exit 1
fi
if [[ ! -x "$RUST_BIN" ]]; then
  echo "Missing release binary: $RUST_BIN" >&2
  echo "Build it first: cargo build --release" >&2
  exit 1
fi
if [[ ! -f "$PY_PROJECT/ldsc.py" ]]; then
  echo "Missing Python LDSC: $PY_PROJECT/ldsc.py" >&2
  exit 1
fi
for ext in bed bim fam; do
  if [[ ! -f "${BFILE}.${ext}" ]]; then
    echo "Missing input file: ${BFILE}.${ext}" >&2
    echo "Create it with: plink1.9 --bfile data/1000G_phase3_common_norel --thin 0.12 --seed 42 --make-bed --out data/bench_200k" >&2
    exit 1
  fi
done

N_SNPS="$(wc -l < "${BFILE}.bim")"
echo "Benchmark dataset: ${BFILE} (${N_SNPS} SNPs)"
echo "Window: --ld-wind-kb ${LD_WIND_KB}"
echo ""

PY_OUT="$OUT_DIR/py_l2"
RS_OUT="$OUT_DIR/rust_l2"
rm -f "${PY_OUT}"*.l2.ldscore.gz "${PY_OUT}"*.l2.M* \
      "${RS_OUT}"*.l2.ldscore.gz "${RS_OUT}"*.l2.M*

# ---------------------------------------------------------------------------
# Timed runs — no --extract
# ---------------------------------------------------------------------------
run_timed() {
  local label="$1"; shift
  local t0 t1 elapsed
  t0="$(date +%s.%N)"
  "$@" >"$OUT_DIR/${label}.stdout" 2>"$OUT_DIR/${label}.stderr"
  t1="$(date +%s.%N)"
  elapsed="$(awk -v a="$t0" -v b="$t1" 'BEGIN { printf "%.3f", (b - a) }')"
  printf "%s\n" "$elapsed" >"$OUT_DIR/${label}.time"
  echo "${label}: ${elapsed}s"
}

PY_CMD=(uv run --project "$PY_PROJECT")
[[ -n "${PY_WITH:-}" ]] && PY_CMD+=(--with "$PY_WITH")
PY_CMD+=(python "$PY_PROJECT/ldsc.py" --l2
         --bfile "$BFILE" --out "$PY_OUT"
         --ld-wind-kb "$LD_WIND_KB" --chunk-size "$CHUNK_SIZE"
         --yes-really)

RS_CMD=("$RUST_BIN" l2
        --bfile "$BFILE" --out "$RS_OUT"
        --ld-wind-kb "$LD_WIND_KB" --chunk-size "$CHUNK_SIZE"
        --yes-really)

echo "Running Python..."
run_timed py "${PY_CMD[@]}"

echo "Running Rust..."
run_timed rust "${RS_CMD[@]}"

# ---------------------------------------------------------------------------
# Parity check
# ---------------------------------------------------------------------------
RUST_FILES=()
if [[ -f "${RS_OUT}.l2.ldscore.gz" ]]; then
  RUST_FILES=("${RS_OUT}.l2.ldscore.gz")
else
  mapfile -t RUST_FILES < <(ls "${RS_OUT}"*.l2.ldscore.gz 2>/dev/null)
fi
if [[ ${#RUST_FILES[@]} -eq 0 ]]; then
  echo "No Rust l2 outputs found — check $OUT_DIR/rust.stderr" >&2
  exit 1
fi

python3 - "$TOL" "${PY_OUT}.l2.ldscore.gz" "${RUST_FILES[@]}" <<'PY'
import gzip, sys
from pathlib import Path

tol = float(sys.argv[1])
py_path = Path(sys.argv[2])
rust_paths = [Path(p) for p in sys.argv[3:]]

def read_ldscore(path):
    rows = {}
    with gzip.open(path, "rt") as f:
        header = f.readline().strip().split("\t")
        snp_idx = header.index("SNP")
        l2_idx  = header.index("L2")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            rows[parts[snp_idx]] = float(parts[l2_idx])
    return rows

py_rows = read_ldscore(py_path)
rs_rows = {}
for p in rust_paths:
    rs_rows.update(read_ldscore(p))

if len(py_rows) != len(rs_rows):
    sys.exit(f"row count mismatch: py={len(py_rows)} rust={len(rs_rows)}")

max_diff = 0.0
worst = None
for snp, py_v in py_rows.items():
    rs_v = rs_rows.get(snp)
    if rs_v is None:
        sys.exit(f"missing SNP in Rust output: {snp}")
    d = abs(py_v - rs_v)
    if d > max_diff:
        max_diff = d
        worst = (snp, py_v, rs_v)

print(f"max_abs_diff={max_diff:.6g}")
if worst:
    print(f"worst_snp={worst[0]}  py={worst[1]:.6g}  rust={worst[2]:.6g}")
if max_diff > tol:
    sys.exit(f"parity FAIL: max_abs_diff={max_diff:.6g} > tol={tol}")
print("parity OK")
PY

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
PY_S="$(cat "$OUT_DIR/py.time")"
RS_S="$(cat "$OUT_DIR/rust.time")"
SPEEDUP="$(awk -v p="$PY_S" -v r="$RS_S" 'BEGIN { printf "%.1f", p/r }')"

echo ""
echo "=== bench_l2 summary ==="
echo "  bfile:       ${BFILE}"
echo "  snps:        ${N_SNPS}"
echo "  python:      ${PY_S}s"
echo "  rust:        ${RS_S}s"
echo "  speedup:     ${SPEEDUP}×"
echo "  output_dir:  ${OUT_DIR}"
