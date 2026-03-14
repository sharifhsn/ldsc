#!/bin/bash
#
# Benchmark --sketch at various d values on the full 1.66M SNP dataset.
# Designed to run on AWS Batch (c6a.4xlarge, 16 vCPU EPYC 7R13).
#
# Produces:
#   1. hyperfine timing for exact + sketch d=25,50,100,200,500,1000
#   2. LD score accuracy comparison (mean L2, bias, median relative error, Pearson r)
#
set -euo pipefail

DATA_DIR="${DATA_DIR:-/data}"
LDSC="${LDSC_BIN:-ldsc}"
BFILE="${DATA_DIR}/1000G_phase3_common_norel"
OUT_DIR="/tmp/sketch_bench"
RUNS="${BENCH_RUNS:-5}"
WARMUP="${BENCH_WARMUP:-1}"

mkdir -p "$OUT_DIR"

echo "=== Sketch Benchmark ==="
echo "Binary: $LDSC"
echo "Data: $BFILE"
echo "Runs: $RUNS, Warmup: $WARMUP"
echo ""

# 1. Exact baseline (also produces reference LD scores)
echo "=== Step 1: Exact baseline ==="
$LDSC l2 --bfile "$BFILE" --out "${OUT_DIR}/exact" --ld-wind-kb 1000 --verbose-timing 2>&1 | tail -3
echo ""

# 2. Sketch at various d values
for d in 25 50 100 200 500 1000; do
    echo "=== Step 2: Sketch d=$d ==="
    $LDSC l2 --bfile "$BFILE" --out "${OUT_DIR}/sketch_${d}" --ld-wind-kb 1000 --sketch "$d" --verbose-timing 2>&1 | tail -3
    echo ""
done

# 3. Hyperfine timing
echo "=== Step 3: Hyperfine timing ==="
hyperfine --warmup "$WARMUP" --runs "$RUNS" --export-json "${OUT_DIR}/timing.json" \
    --command-name "exact" \
    "$LDSC l2 --bfile $BFILE --out ${OUT_DIR}/hf_exact --ld-wind-kb 1000" \
    --command-name "sketch_d25" \
    "$LDSC l2 --bfile $BFILE --out ${OUT_DIR}/hf_s25 --ld-wind-kb 1000 --sketch 25" \
    --command-name "sketch_d50" \
    "$LDSC l2 --bfile $BFILE --out ${OUT_DIR}/hf_s50 --ld-wind-kb 1000 --sketch 50" \
    --command-name "sketch_d100" \
    "$LDSC l2 --bfile $BFILE --out ${OUT_DIR}/hf_s100 --ld-wind-kb 1000 --sketch 100" \
    --command-name "sketch_d200" \
    "$LDSC l2 --bfile $BFILE --out ${OUT_DIR}/hf_s200 --ld-wind-kb 1000 --sketch 200" \
    --command-name "sketch_d500" \
    "$LDSC l2 --bfile $BFILE --out ${OUT_DIR}/hf_s500 --ld-wind-kb 1000 --sketch 500" \
    --command-name "sketch_d1000" \
    "$LDSC l2 --bfile $BFILE --out ${OUT_DIR}/hf_s1000 --ld-wind-kb 1000 --sketch 1000"
echo ""

# 4. Accuracy comparison
echo "=== Step 4: LD score accuracy ==="
python3 -c "
import gzip, math, os

def read_all_ldscore(prefix):
    data = {}
    for c in range(1, 23):
        path = f'{prefix}{c}.l2.ldscore.gz'
        if not os.path.exists(path):
            continue
        with gzip.open(path, 'rt') as f:
            lines = f.readlines()
        header = lines[0].strip().split()
        l2_idx = header.index('L2')
        snp_idx = header.index('SNP')
        for line in lines[1:]:
            parts = line.strip().split()
            data[parts[snp_idx]] = float(parts[l2_idx])
    return data

def pearson_r(xs, ys):
    n = len(xs)
    mx = sum(xs)/n; my = sum(ys)/n
    cov = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    sx = math.sqrt(sum((x-mx)**2 for x in xs))
    sy = math.sqrt(sum((y-my)**2 for y in ys))
    return cov / (sx * sy) if sx > 0 and sy > 0 else 0

out = '${OUT_DIR}'
exact = read_all_ldscore(f'{out}/exact')
print(f'Exact: {len(exact)} SNPs, mean L2={sum(exact.values())/len(exact):.3f}')
print()
print(f'{\"d\":>5s} {\"mean_L2\":>8s} {\"bias\":>8s} {\"bias%\":>7s} {\"med_rel%\":>8s} {\"p90_rel%\":>8s} {\"pearson\":>8s} {\"spearman\":>8s}')

for d in [25, 50, 100, 200, 500, 1000]:
    sketch = read_all_ldscore(f'{out}/sketch_{d}')
    common = sorted(set(exact) & set(sketch))
    if not common:
        continue
    xs = [exact[s] for s in common]
    ys = [sketch[s] for s in common]
    me = sum(xs)/len(xs); ms = sum(ys)/len(ys)
    rel = sorted([abs(y-x)/max(abs(x),0.01) for x,y in zip(xs,ys)])
    med_rel = rel[len(rel)//2]
    p90_rel = rel[int(0.9*len(rel))]
    r = pearson_r(xs, ys)
    # Spearman
    n = len(xs)
    rx = sorted(range(n), key=lambda i: xs[i])
    ry = sorted(range(n), key=lambda i: ys[i])
    rank_x = [0]*n; rank_y = [0]*n
    for i,v in enumerate(rx): rank_x[v] = i
    for i,v in enumerate(ry): rank_y[v] = i
    rho = pearson_r(rank_x, rank_y)
    print(f'{d:5d} {ms:8.3f} {ms-me:+8.3f} {100*(ms-me)/me:+7.1f}% {100*med_rel:8.1f} {100*p90_rel:8.1f} {r:8.4f} {rho:8.4f}')
"

echo ""
echo "=== Benchmark complete ==="
