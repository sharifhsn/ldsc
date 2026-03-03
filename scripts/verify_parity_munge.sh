#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="$ROOT/perf/parity_munge"
PY3_PROJECT="$ROOT/ldsc_py3"
RUST_BIN="$ROOT/target/release/ldsc"
UV_PYTHON="${UV_PYTHON:-3.9}"

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

INPUT_SUMSTATS="${INPUT_SUMSTATS:-}"
LINES="${LINES:-200000}"

if [[ -z "$INPUT_SUMSTATS" ]]; then
  if [[ -f "$ROOT/data/biomarkers-30600-both_sexes-irnt.sample.tsv" ]]; then
    INPUT_SUMSTATS="$OUT_DIR/sample_input.tsv"
    head -n $((LINES + 1)) "$ROOT/data/biomarkers-30600-both_sexes-irnt.sample.tsv" > "$INPUT_SUMSTATS"
  else
    INPUT_SUMSTATS="$OUT_DIR/toy_input.tsv"
    TOY_SCRIPT="$OUT_DIR/_toy_input.py"
    cat >"$TOY_SCRIPT" <<'PY'
# /// script
# requires-python = ">=3.9"
# dependencies = []
# ///
from pathlib import Path
import random
out = Path('/Users/sharif/Code/ldsc/perf/parity_munge/toy_input.tsv')
random.seed(0)
with out.open('w') as f:
    f.write('SNP\tA1\tA2\tBETA\tSE\tP\tN\tINFO\tFRQ\n')
    for i in range(1, 1001):
        snp = f'rs{i}'
        a1, a2 = ('A','C') if i % 2 == 0 else ('G','T')
        beta = (i % 10 - 5) * 0.01
        se = 0.05 + (i % 7) * 0.005
        p = max(1e-6, min(1.0, 0.5 + (i % 9 - 4) * 0.05))
        n = 10000 + (i % 3) * 100
        info = 0.95 if i % 13 else 0.5
        frq = 0.05 + (i % 20) * 0.01
        f.write(f'{snp}\t{a1}\t{a2}\t{beta:.6f}\t{se:.6f}\t{p:.6g}\t{n}\t{info:.3f}\t{frq:.3f}\n')
print('Wrote', out)
PY
    uv run --python "$UV_PYTHON" --script "$TOY_SCRIPT" >/dev/null
  fi
fi

PY_OUT="$OUT_DIR/py3_munge"
RS_OUT="$OUT_DIR/rust_munge"

rm -f "${PY_OUT}.sumstats.gz" "${RS_OUT}.sumstats.gz"

MUNGE_WRAPPER="$OUT_DIR/_munge_wrapper.py"
cat >"$MUNGE_WRAPPER" <<'PY'
# /// script
# requires-python = ">=3.9"
# dependencies = [
#   "numpy==1.21.5",
#   "pandas==1.3.3",
#   "scipy==1.7.3",
#   "bitarray==2",
#   "python-dateutil==2.8.2",
#   "pytz==2022.4",
#   "six==1.16.0",
# ]
# ///
import os
import runpy
import sys

script = os.environ["MUNGE_SCRIPT"]
project_root = os.path.dirname(script)
if project_root not in sys.path:
    sys.path.insert(0, project_root)

sys.argv = [script] + sys.argv[1:]
runpy.run_path(script, run_name="__main__")
PY

MUNGE_SCRIPT="$PY3_PROJECT/munge_sumstats.py" \
  uv run --python "$UV_PYTHON" --script "$MUNGE_WRAPPER" --sumstats "$INPUT_SUMSTATS" --out "$PY_OUT" >/dev/null

"$RUST_BIN" munge-sumstats --sumstats "$INPUT_SUMSTATS" --out "$RS_OUT" >/dev/null

CMP_SCRIPT="$OUT_DIR/_compare.py"
cat >"$CMP_SCRIPT" <<'PY'
# /// script
# requires-python = ">=3.9"
# dependencies = []
# ///
import csv
import gzip
from pathlib import Path

py_path = Path('/Users/sharif/Code/ldsc/perf/parity_munge/py3_munge.sumstats.gz')
rs_path = Path('/Users/sharif/Code/ldsc/perf/parity_munge/rust_munge.sumstats.gz')

def read(path):
    with gzip.open(path, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = {}
        cols = reader.fieldnames
        for row in reader:
            rows[row['SNP']] = row
    return rows, cols

py_rows, py_cols = read(py_path)
rs_rows, rs_cols = read(rs_path)

print(f'Python rows: {len(py_rows)}  Rust rows: {len(rs_rows)}')

if set(py_cols) != set(rs_cols):
    print('Column mismatch:', sorted(py_cols), sorted(rs_cols))

set_py = set(py_rows)
set_rs = set(rs_rows)
if set_py != set_rs:
    print('SNP set mismatch: py-only', len(set_py - set_rs), 'rs-only', len(set_rs - set_py))
    raise SystemExit(1)

num_cols = [c for c in ['Z', 'N', 'FRQ'] if c in py_cols]
max_diff = {c: 0.0 for c in num_cols}

for snp in sorted(set_py):
    py = py_rows[snp]
    rs = rs_rows[snp]
    for col in ['A1', 'A2']:
        if col in py_cols and py[col] != rs[col]:
            print(f'{col} mismatch at {snp}: {py[col]} vs {rs[col]}')
            raise SystemExit(1)
    for col in num_cols:
        a = float(py[col])
        b = float(rs[col])
        diff = abs(a - b)
        if diff > max_diff[col]:
            max_diff[col] = diff

for col in num_cols:
    print(f'{col}: max_abs_diff={max_diff[col]:.6g}')
    if max_diff[col] > 1e-6:
        raise SystemExit(1)

print('Munge parity OK')
PY

uv run --python "$UV_PYTHON" --script "$CMP_SCRIPT"
