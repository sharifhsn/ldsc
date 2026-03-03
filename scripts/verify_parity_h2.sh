#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="$ROOT/perf/parity"
PY3_PROJECT="$ROOT/ldsc_py3"
RUST_BIN="$ROOT/target/release/ldsc"

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

PY3_PYTHON="${PY3_PYTHON:-3.9}"
PY_WRAPPER="$ROOT/scripts/_ldsc_wrapper.py"
if [[ ! -f "$PY_WRAPPER" ]]; then
  echo "Missing Python wrapper: $PY_WRAPPER" >&2
  exit 1
fi
PY3_CMD=(env LDSC_SCRIPT="$PY3_PROJECT/ldsc.py" uv run --python "$PY3_PYTHON" --script "$PY_WRAPPER")

GEN_SCRIPT="$OUT_DIR/_gen_toy_h2_cases.py"
cat >"$GEN_SCRIPT" <<'PY'
# /// script
# requires-python = ">=3.9"
# dependencies = []
# ///
import gzip, random
from pathlib import Path

out_dir = Path('/Users/sharif/Code/ldsc/perf/parity')
out_dir.mkdir(parents=True, exist_ok=True)

# Case 1: unsorted LD order, constant N
snps = [f'rs{i}' for i in range(1, 21)]
random.seed(0)

sumstats_order = snps[:]
random.shuffle(sumstats_order)

sumstats_path = out_dir / 'toy.sumstats.gz'
with gzip.open(sumstats_path, 'wt') as f:
    f.write('SNP\tZ\tN\n')
    for i, snp in enumerate(sumstats_order):
        z = (i - 5) * 0.1
        f.write(f'{snp}\t{z:.4f}\t1000\n')

ld_order = snps[:]
random.shuffle(ld_order)

ld_path = out_dir / 'ref.l2.ldscore.gz'
with gzip.open(ld_path, 'wt') as f:
    f.write('CHR\tSNP\tBP\tL2\n')
    for i, snp in enumerate(ld_order):
        chr_ = 1 + (i % 2)
        bp = 100 + i * 10
        l2 = 10 + i
        f.write(f'{chr_}\t{snp}\t{bp}\t{l2:.3f}\n')

w_path = out_dir / 'w.l2.ldscore.gz'
with gzip.open(w_path, 'wt') as f:
    f.write('CHR\tSNP\tBP\tL2\n')
    for i, snp in enumerate(ld_order):
        chr_ = 1 + (i % 2)
        bp = 100 + i * 10
        l2 = 5 + i
        f.write(f'{chr_}\t{snp}\t{bp}\t{l2:.3f}\n')

# Case 2: varying N + L2 < 1 to test clamp + median even case
snps = [f'rs{i}' for i in range(1, 31)]
random.seed(1)

sumstats_order = snps[:]
random.shuffle(sumstats_order)

sumstats_path = out_dir / 'toy_varn.sumstats.gz'
with gzip.open(sumstats_path, 'wt') as f:
    f.write('SNP\tZ\tN\n')
    for i, snp in enumerate(sumstats_order):
        z = (i - 10) * 0.05
        n = 800 + (i % 7) * 50
        f.write(f'{snp}\t{z:.4f}\t{n}\n')

ld_order = snps[:]
random.shuffle(ld_order)

ld_path = out_dir / 'ref_varn.l2.ldscore.gz'
with gzip.open(ld_path, 'wt') as f:
    f.write('CHR\tSNP\tBP\tL2\n')
    for i, snp in enumerate(ld_order):
        chr_ = 1 + (i % 2)
        bp = 100 + i * 7
        l2 = 0.5 + (i % 5)  # includes <1.0 to test clamp
        f.write(f'{chr_}\t{snp}\t{bp}\t{l2:.3f}\n')

w_path = out_dir / 'w_varn.l2.ldscore.gz'
with gzip.open(w_path, 'wt') as f:
    f.write('CHR\tSNP\tBP\tL2\n')
    for i, snp in enumerate(ld_order):
        chr_ = 1 + (i % 2)
        bp = 100 + i * 7
        l2 = 0.3 + (i % 4)  # includes <1.0
        f.write(f'{chr_}\t{snp}\t{bp}\t{l2:.3f}\n')

# Case 3: chisq-max boundary
snps = [f'rs{i}' for i in range(1, 11)]

sumstats_path = out_dir / 'toy_chisq.sumstats.gz'
with gzip.open(sumstats_path, 'wt') as f:
    f.write('SNP\tZ\tN\n')
    zs = [5.0, 5.4772256, 5.5, 3.0, 2.0, 1.0, 0.5, 4.0, 5.4772256, 6.0]  # includes sqrt(30)
    for snp, z in zip(snps, zs):
        f.write(f'{snp}\t{z:.7f}\t1000\n')

ld_path = out_dir / 'ref_chisq.l2.ldscore.gz'
with gzip.open(ld_path, 'wt') as f:
    f.write('CHR\tSNP\tBP\tL2\n')
    for i, snp in enumerate(snps):
        f.write(f'1\t{snp}\t{100+i}\t{10+i}\n')

w_path = out_dir / 'w_chisq.l2.ldscore.gz'
with gzip.open(w_path, 'wt') as f:
    f.write('CHR\tSNP\tBP\tL2\n')
    for i, snp in enumerate(snps):
        f.write(f'1\t{snp}\t{100+i}\t{5+i}\n')

print('Wrote toy datasets to', out_dir)
PY
uv run --python "$PY3_PYTHON" --script "$GEN_SCRIPT" >/dev/null

CMP_SCRIPT="$OUT_DIR/_compare_h2_stdout.py"
cat >"$CMP_SCRIPT" <<'PY'
# /// script
# requires-python = ">=3.9"
# dependencies = []
# ///
import re
import sys
from pathlib import Path

def parse(path):
    text = Path(path).read_text()
    out = {}
    m = re.search(r'Total Observed scale h2: ([-0-9.]+) \\(([-0-9.]+)\\)', text)
    if m:
        out['h2'] = float(m.group(1))
        out['h2_se'] = float(m.group(2))
    m = re.search(r'Lambda GC: ([-0-9.]+)', text)
    if m:
        out['lambda_gc'] = float(m.group(1))
    m = re.search(r'Mean Chi\\^2: ([-0-9.]+)', text)
    if m:
        out['mean_chi2'] = float(m.group(1))
    m = re.search(r'Intercept: ([-0-9.]+) \\(([-0-9.]+)\\)', text)
    if m:
        out['intercept'] = float(m.group(1))
        out['intercept_se'] = float(m.group(2))
    elif 'Intercept: constrained to' in text:
        m = re.search(r'Intercept: constrained to ([-0-9.]+)', text)
        out['intercept'] = float(m.group(1))
    if 'Ratio < 0' in text:
        out['ratio_tag'] = 'lt0'
    else:
        m = re.search(r'Ratio: ([-0-9.]+) \\(([-0-9.]+)\\)', text)
        if m:
            out['ratio'] = float(m.group(1))
            out['ratio_se'] = float(m.group(2))
        elif 'Ratio: NA' in text:
            out['ratio_tag'] = 'na'
    return out

def close(a, b, tol=5e-4):
    return abs(a - b) <= tol

py = parse(sys.argv[1])
rs = parse(sys.argv[2])

keys = ['h2','h2_se','lambda_gc','mean_chi2','intercept','intercept_se']
for k in keys:
    if k in py or k in rs:
        if k not in py or k not in rs:
            print('Missing', k, 'py', k in py, 'rs', k in rs)
            sys.exit(1)
        if not close(py[k], rs[k]):
            print('Mismatch', k, py[k], rs[k])
            sys.exit(1)

print('OK')
PY

run_case() {
  local name="$1"
  local sumstats="$2"
  local ref_prefix="$3"
  local w_prefix="$4"
  local m_val="$5"
  local n_blocks="${6:-}"
  local extra_py=()
  local extra_rs=()

  shift 6
  while [[ "$#" -gt 0 ]]; do
    extra_py+=("$1")
    extra_rs+=("$1")
    shift
  done

  local py_out="$OUT_DIR/py3_${name}"
  local rs_out="$OUT_DIR/rust_${name}"

  echo "==> ${name}: python"
  local py_cmd=("${PY3_CMD[@]}" --h2 "$sumstats" --ref-ld "$ref_prefix" --w-ld "$w_prefix" \
    --out "$py_out" --M "$m_val" --two-step 30)
  if [[ -n "$n_blocks" ]]; then
    py_cmd+=(--n-blocks "$n_blocks")
  fi
  if (( ${#extra_py[@]} )); then
    py_cmd+=("${extra_py[@]}")
  fi
  "${py_cmd[@]}" >"$py_out.stdout" 2>&1

  echo "==> ${name}: rust"
  local rs_cmd=("$RUST_BIN" h2 --h2 "$sumstats" --ref-ld "${ref_prefix}.l2.ldscore.gz" \
    --w-ld "${w_prefix}.l2.ldscore.gz" --out "$rs_out" --M "$m_val" --two-step 30)
  if [[ -n "$n_blocks" ]]; then
    rs_cmd+=(--n-blocks "$n_blocks")
  fi
  if (( ${#extra_rs[@]} )); then
    rs_cmd+=("${extra_rs[@]}")
  fi
  "${rs_cmd[@]}" >"$rs_out.stdout" 2>&1

  uv run --python "$PY3_PYTHON" --script "$CMP_SCRIPT" "$py_out.stdout" "$rs_out.stdout"
}

run_case "toy" "$OUT_DIR/toy.sumstats.gz" "$OUT_DIR/ref" "$OUT_DIR/w" 1000 20
run_case "toy_varn" "$OUT_DIR/toy_varn.sumstats.gz" "$OUT_DIR/ref_varn" "$OUT_DIR/w_varn" 1000 30
run_case "toy_chisq" "$OUT_DIR/toy_chisq.sumstats.gz" "$OUT_DIR/ref_chisq" "$OUT_DIR/w_chisq" 1000 5 --chisq-max 30

echo "All parity checks passed."
