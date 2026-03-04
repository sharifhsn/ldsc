#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${OUT_DIR:-$ROOT/perf/l2_rust_hyperfine}"
RUST_BIN="${RUST_BIN:-$ROOT/target/release/ldsc}"
BFILE="${BFILE:-$ROOT/data/1000G.EUR.QC}"
LD_WIND_CM="${LD_WIND_CM:-1}"
EXTRACT_N="${EXTRACT_N:-10000}"
EXTRACT_FILE="${EXTRACT_FILE:-$OUT_DIR/extract_${EXTRACT_N}.snps}"
RAYON_THREADS="${RAYON_THREADS:-4}"
RUNS="${RUNS:-3}"
WARMUP="${WARMUP:-1}"
VERIFY_EXTRACT="${VERIFY_EXTRACT:-1}"

mkdir -p "$OUT_DIR"

if ! command -v hyperfine >/dev/null 2>&1; then
  echo "hyperfine is required but not found in PATH." >&2
  exit 1
fi

if [[ ! -x "$RUST_BIN" ]]; then
  echo "Missing release binary: $RUST_BIN" >&2
  echo "Build it first with: cargo build --release" >&2
  exit 1
fi

for ext in bed bim fam; do
  if [[ ! -f "${BFILE}.${ext}" ]]; then
    echo "Missing input file: ${BFILE}.${ext}" >&2
    exit 1
  fi
done

if [[ ! -s "$EXTRACT_FILE" ]] || [[ "$(wc -l < "$EXTRACT_FILE")" -ne "$EXTRACT_N" ]]; then
  awk -v n="$EXTRACT_N" 'NR<=n {print $2} NR==n {exit}' "${BFILE}.bim" > "$EXTRACT_FILE"
fi

if [[ "$(wc -l < "$EXTRACT_FILE")" -ne "$EXTRACT_N" ]]; then
  echo "Extract file '$EXTRACT_FILE' does not contain $EXTRACT_N SNP IDs." >&2
  exit 1
fi

OUT_PREFIX="$OUT_DIR/rust_l2_extract_${EXTRACT_N}"
HF_JSON="$OUT_DIR/hyperfine_rust_l2_extract_${EXTRACT_N}.json"
HF_MD="$OUT_DIR/hyperfine_rust_l2_extract_${EXTRACT_N}.md"
HF_TXT="$OUT_DIR/hyperfine_rust_l2_extract_${EXTRACT_N}.txt"

CMD="$RUST_BIN --rayon-threads $RAYON_THREADS l2 --bfile $BFILE --out $OUT_PREFIX --ld-wind-cm $LD_WIND_CM --extract $EXTRACT_FILE"

if [[ "$VERIFY_EXTRACT" == "1" ]]; then
  VERIFY_LOG="$OUT_DIR/preflight_extract_${EXTRACT_N}.log"
  rm -f "${OUT_PREFIX}"* "$VERIFY_LOG"
  # Preflight once so we can prove --extract reduced the active SNP set.
  $CMD >"$VERIFY_LOG" 2>&1
  EXTRACTED_ACTUAL="$(sed -n 's/^  After --extract: \([0-9]\+\) SNPs$/\1/p' "$VERIFY_LOG" | tail -n 1)"
  if [[ -z "$EXTRACTED_ACTUAL" ]]; then
    echo "Could not verify extracted SNP count from $VERIFY_LOG" >&2
    exit 1
  fi
  if [[ "$EXTRACTED_ACTUAL" != "$EXTRACT_N" ]]; then
    echo "Extract mismatch: requested $EXTRACT_N SNPs, ldsc reported $EXTRACTED_ACTUAL" >&2
    exit 1
  fi
  echo "Verified extract shrink: requested=$EXTRACT_N, ldsc_after_extract=$EXTRACTED_ACTUAL"
  rm -f "${OUT_PREFIX}"*
fi

hyperfine \
  --warmup "$WARMUP" \
  --runs "$RUNS" \
  --prepare "rm -f ${OUT_PREFIX}*" \
  --cleanup "rm -f ${OUT_PREFIX}*" \
  --export-json "$HF_JSON" \
  --export-markdown "$HF_MD" \
  --command-name "rust_l2_extract_${EXTRACT_N}" \
  "$CMD" | tee "$HF_TXT"

echo "Saved benchmark outputs:" 

echo "  $HF_TXT"
echo "  $HF_JSON"
echo "  $HF_MD"
