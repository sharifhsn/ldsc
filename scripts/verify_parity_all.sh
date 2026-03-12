#!/usr/bin/env bash
#
# Comprehensive parity check: all mode × precision combinations vs Python LDSC.
#
# Modes tested (each run against Python's exact f64 output as ground truth):
#   exact-f64       --                         (must match Python, max_abs_diff < 1e-3)
#   exact-f32       --fast-f32                 (must match Python, max_abs_diff < 2e-2)
#   sketch-50-f64   --sketch 50                (reports Pearson r, median % error)
#   sketch-100-f64  --sketch 100
#   sketch-200-f64  --sketch 200
#   sketch-50-f32   --sketch 50  --fast-f32
#   sketch-100-f32  --sketch 100 --fast-f32
#   sketch-200-f32  --sketch 200 --fast-f32
#
# Uses chr22-only data by default (single chromosome, clean output naming).
#
# Usage:
#   ./scripts/verify_parity_all.sh [--bfile path] [--ld-wind-kb 1000]
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="$ROOT/perf/parity_all"
RUST_BIN="$ROOT/target/release/ldsc"
BFILE="${BFILE:-$ROOT/data/bench_5k_chr22}"
LD_WIND_KB="${LD_WIND_KB:-1000}"
CHUNK_SIZE="${CHUNK_SIZE:-200}"

while [[ $# -gt 0 ]]; do
    case $1 in
        --bfile) BFILE="$2"; shift 2 ;;
        --ld-wind-kb) LD_WIND_KB="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

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

COMMON_ARGS="--ld-wind-kb $LD_WIND_KB --chunk-size $CHUNK_SIZE"

# ── Step 1: Python reference ──────────────────────────────────────────────────
PY_DIR="$OUT_DIR/py"; mkdir -p "$PY_DIR"
PY_OUT="$PY_DIR/l2"
echo "=== Running Python LDSC (reference) ==="
uv run --project "$ROOT/ldsc_py" --with "bitarray==2" python "$ROOT/ldsc_py/ldsc.py" \
    --l2 --bfile "$BFILE" --out "$PY_OUT" $COMMON_ARGS \
    >"$PY_DIR/stdout.txt" 2>&1
# Python outputs a single combined file
PY_FILE="$PY_OUT.l2.ldscore.gz"
echo "Python done: $(zcat "$PY_FILE" | wc -l) lines (incl. header)"
echo ""

# ── Step 2: Define test variants ──────────────────────────────────────────────
declare -a LABELS=()
declare -a EXTRA_ARGS=()
# tolerance: "exact64" | "exact32" | "approx"
declare -a MODE=()

LABELS+=("exact-f64");    EXTRA_ARGS+=("");                           MODE+=("exact64")
LABELS+=("exact-f32");    EXTRA_ARGS+=("--fast-f32");                 MODE+=("exact32")
LABELS+=("sketch-50-f64");  EXTRA_ARGS+=("--sketch 50");              MODE+=("approx")
LABELS+=("sketch-100-f64"); EXTRA_ARGS+=("--sketch 100");             MODE+=("approx")
LABELS+=("sketch-200-f64"); EXTRA_ARGS+=("--sketch 200");             MODE+=("approx")
LABELS+=("sketch-50-f32");  EXTRA_ARGS+=("--sketch 50 --fast-f32");   MODE+=("approx")
LABELS+=("sketch-100-f32"); EXTRA_ARGS+=("--sketch 100 --fast-f32");  MODE+=("approx")
LABELS+=("sketch-200-f32"); EXTRA_ARGS+=("--sketch 200 --fast-f32");  MODE+=("approx")

# ── Step 3: Run each Rust variant ─────────────────────────────────────────────
echo "=== Running Rust variants ==="
for i in "${!LABELS[@]}"; do
    label="${LABELS[$i]}"
    extra="${EXTRA_ARGS[$i]}"
    rs_dir="$OUT_DIR/$label"; mkdir -p "$rs_dir"
    rs_out="$rs_dir/l2"
    printf "  %-22s ... " "$label"
    # shellcheck disable=SC2086
    "$RUST_BIN" l2 --bfile "$BFILE" --out "$rs_out" $COMMON_ARGS $extra \
        >"$rs_dir/stdout.txt" 2>&1
    echo "done"
done
echo ""

# ── Step 4: Comparison script ─────────────────────────────────────────────────
CMP_SCRIPT="$OUT_DIR/_compare.py"
cat >"$CMP_SCRIPT" <<'PY'
# /// script
# requires-python = ">=3.9"
# dependencies = []
# ///
import gzip, sys, math
from pathlib import Path
from glob import glob

py_path = Path(sys.argv[1])
rs_glob = sys.argv[2]   # directory prefix, e.g. /tmp/exact-f64/l2
mode    = sys.argv[3]   # "exact64" | "exact32" | "approx"

def read_ldscore(path):
    rows = {}
    with gzip.open(path, 'rt') as f:
        header = f.readline().strip().split('\t')
        l2_idx = header.index('L2')
        snp_idx = header.index('SNP')
        for line in f:
            parts = line.rstrip('\n').split('\t')
            rows[parts[snp_idx]] = float(parts[l2_idx])
    return rows

py = read_ldscore(py_path)

rs = {}
for p in sorted(glob(rs_glob + "*.l2.ldscore.gz")):
    rs.update(read_ldscore(p))

if len(py) != len(rs):
    print(f"  ROW MISMATCH: py={len(py)} rust={len(rs)}")
    sys.exit(1)

py_vals, rs_vals = [], []
for snp, pv in py.items():
    rv = rs.get(snp)
    if rv is None:
        print(f"  MISSING SNP: {snp}")
        sys.exit(1)
    py_vals.append(pv)
    rs_vals.append(rv)

n = len(py_vals)
max_diff = max(abs(p - r) for p, r in zip(py_vals, rs_vals))
mean_py  = sum(py_vals) / n
mean_rs  = sum(rs_vals) / n
rmse     = math.sqrt(sum((p - r)**2 for p, r in zip(py_vals, rs_vals)) / n)

num   = sum((p - mean_py) * (r - mean_rs) for p, r in zip(py_vals, rs_vals))
denom = math.sqrt(sum((p - mean_py)**2 for p in py_vals) *
                  sum((r - mean_rs)**2  for r in rs_vals))
pearson_r = num / denom if denom > 0 else 0.0

pct_errs = sorted(100 * abs(p - r) / abs(p) for p, r in zip(py_vals, rs_vals) if abs(p) > 0.1)
median_pct = pct_errs[len(pct_errs) // 2] if pct_errs else 0.0

if mode == "exact64":
    tol = 1e-3
    ok  = max_diff < tol
    status = "OK" if ok else f"FAIL (exceeds {tol:.0e})"
    print(f"  max_abs_diff={max_diff:.3g}  RMSE={rmse:.3g}  r={pearson_r:.6f}  {status}")
    if not ok: sys.exit(1)
elif mode == "exact32":
    tol = 2e-2
    ok  = max_diff < tol
    status = "OK" if ok else f"FAIL (exceeds {tol:.0e})"
    print(f"  max_abs_diff={max_diff:.3g}  RMSE={rmse:.3g}  r={pearson_r:.6f}  {status}")
    if not ok: sys.exit(1)
else:
    print(f"  max_abs_diff={max_diff:.3g}  RMSE={rmse:.3g}  r={pearson_r:.6f}  median_pct_err={median_pct:.1f}%")
PY

# ── Step 5: Compare each variant ──────────────────────────────────────────────
echo "=== Parity Results (vs Python exact f64) ==="
printf "%-22s  %s\n" "Mode" "Stats"
printf "%-22s  %s\n" "----" "-----"

FAILED=()
for i in "${!LABELS[@]}"; do
    label="${LABELS[$i]}"
    mode="${MODE[$i]}"
    rs_prefix="$OUT_DIR/$label/l2"
    printf "%-22s " "$label"
    if result=$(uv run --script "$CMP_SCRIPT" "$PY_FILE" "$rs_prefix" "$mode" 2>&1); then
        echo "$result"
    else
        echo "FAILED: $result"
        FAILED+=("$label")
    fi
done

echo ""
if [[ ${#FAILED[@]} -eq 0 ]]; then
    echo "All parity checks passed."
else
    echo "FAILED variants: ${FAILED[*]}"
    exit 1
fi
