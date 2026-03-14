# Implementation Plan: Fused BED-Decode, Subsampling, CountSketch

Three independent optimizations to `compute_ldscore_global` in `src/l2.rs`.
Each is a separate commit. None touches the existing exact (non-sketch) path.

---

## Dependencies Between Optimizations

```
[1] Fused BED-Decode-Normalize-Project ──(independent)──> [2] Subsample
                                        ──(independent)──> [3] CountSketch
```

All three are independent. Implement in order 1 → 2 → 3 because:
- (1) is the highest-value optimization (32x memory traffic reduction) and is **exact-parity** with the current sketch path, making it the safest to validate.
- (2) is a CLI convenience that leverages existing `--keep` infrastructure.
- (3) introduces a new sketch method that replaces the projection GEMM.

CountSketch (3) and the Fused kernel (1) can compose: CountSketch's hash-scatter can be tiled the same way as the fused Rademacher projection. But the initial implementation of (3) should be standalone (using the existing `b_full` buffer for decode+normalize, then a custom scatter-add). A combined fused-CountSketch can follow as a second pass.

---

## Optimization 1: Fused BED-Decode-Normalize-Project

### 1.1 Verification: Mathematical Correctness

**Current path** (sketch mode):
```
BED bytes → decode → b_full[N × c] (f32/f64)     // ~40MB at N=50K, c=200
for j in 0..c:
    normalize b_full[:,j] in-place (mean=0, var=1)
b_mat[d × c] = P[d × N] × b_full[N × c]          // faer GEMM
for j in 0..c:
    rescale b_mat[:,j] so ||b_mat[:,j]||² = N      // ratio estimator
```

**Proposed fused path** (same output, different memory access pattern):
```
Pass 1 — Statistics (read ~2.5MB packed BED bytes):
    for j in 0..c:
        scan packed bytes for SNP j → compute sum, count, sum_sq
        derive mean_j, inv_std_j (identical to normalize_col_{f32,f64}_with_stats)

Pass 2 — Tiled Decode+Normalize+Project:
    accum[d × c] = zeros   // output accumulator, ~200×200×8 = 320KB (fits L2)
    for tile_start in (0..N).step_by(TILE):
        tile_end = min(tile_start + TILE, N)
        tile_n = tile_end - tile_start
        // Decode tile: b_tile[tile_n × c] from packed bytes
        for j in 0..c:
            for i in 0..tile_n:
                geno = lut[byte >> shift & 0x3]
                b_tile[i][j] = (geno - mean_j) * inv_std_j   // or 0 if NaN
        // Project tile: accum += P[:,tile_start..tile_end]^T × b_tile is wrong shape
        // Actually: accum[d × c] += P[d × tile_n] × b_tile[tile_n × c]
        // P_tile = P[:, tile_start..tile_end]  (d × tile_n, a submatrix view)
        matmul(accum, Accum::Add, P_tile, b_tile, 1.0, Par::rayon(0))

    // Ratio estimator re-normalization (unchanged)
    for j in 0..c:
        nrm_sq = ||accum[:,j]||²
        accum[:,j] *= sqrt(N / nrm_sq)
    b_mat[:, 0..c] = accum
```

**Why this is identical**: The projection is linear: `P × B = P × [B_tile1; B_tile2; ...] = P_tile1 × B_tile1 + P_tile2 × B_tile2 + ...`. Matrix multiplication distributes over block-row decomposition. The normalization (center + scale) is applied per-SNP with the same mean/std, so each tile element is identically transformed. QED.

**Memory traffic analysis**:
- Pass 1: Read 2.5MB packed BED bytes (bytes_per_snp × c = 12,500 × 200 at N=50K).
- Pass 2: Read 2.5MB packed BED bytes again + read P_tile (d × TILE × 8B).
  With TILE=512, d=100: P_tile = 100×512×8 = 400KB per tile. b_tile = 512×200×4 = 400KB.
  Total per tile: ~800KB working set (fits L2 cache).
  Total across all tiles: 2×2.5MB (BED) + d×N×8B (P read once) ≈ 5MB + 4MB = 9MB.
- Current path: 40MB (b_full) + 4MB (P) + 40MB (GEMM reads b_full) ≈ 84MB.
- **Reduction: ~9x memory traffic** (not 32x as initially estimated — P must still be read).

**Parallelism improvement**: The current sketch projection `P × b_full` is a single GEMM with d=100, c=200, N=50K → 100×200×50K = 1B FLOPs. This is large enough for rayon. But the tiled approach naturally parallelizes across tiles (each tile's GEMM is independent into a thread-local accumulator), giving TILE-level parallelism. With N=50K, TILE=512: 98 tiles × c=200 columns each. We can `par_iter` over tiles, with each thread owning a local `accum[d × c]` that gets summed at the end.

However, there is a subtlety: **faer's `Par::rayon(0)` already parallelizes the single large GEMM internally**. Tiling and parallelizing across tiles with thread-local accumulators may not beat a single optimized GEMM — it depends on whether the working-set reduction (L2 vs L3) outweighs the reduced GEMM efficiency on smaller tiles. This needs benchmarking.

**Recommended approach**: Implement the tiled path with a single-threaded accumulation first (sequential tile loop, Par::Seq for each small tile GEMM). Benchmark against current. If slower, try Par::rayon(0) on a single large GEMM but with the b_tile approach (avoid materializing b_full). If that's also slower, fall back to the simplest version: just avoid allocating b_full by decoding directly into smaller tiles but still using a single GEMM.

### 1.2 Implementation Steps

**File: `src/bed.rs`** — Add raw byte access for statistics pass.

```rust
/// Read raw packed BED bytes for a contiguous block of SNPs into a pre-allocated buffer.
/// Returns a slice of `count * bytes_per_snp` bytes.
impl ChunkReader<T> {
    pub fn raw_bytes(&self) -> &[u8] {
        &self.block_buf
    }

    /// Read raw bytes for `count` SNPs sequentially (no decode).
    /// After this call, `raw_bytes()[..count * bytes_per_snp]` contains the packed data.
    pub fn read_raw_next(&mut self, bed: &mut Bed, count: usize) -> Result<()> {
        let total_bytes = count * self.bytes_per_snp;
        bed.read_next_block(count, &mut self.block_buf[..total_bytes])?;
        Ok(())
    }
}
```

Expose `bytes_per_snp` and `n_iid()` on ChunkReader (currently private fields).

**File: `src/l2.rs`** — New function `sketch_fused_bed_project`.

Step-by-step:

1. **Add a new function** (~150 lines):

```rust
/// Fused BED-decode-normalize-project for sketch mode.
/// Reads raw BED bytes, computes per-SNP statistics, then tiles over individuals
/// to decode+normalize+project without materializing the full N×c intermediate matrix.
///
/// Returns the projected matrix (d × c) with ratio estimator normalization applied.
fn sketch_fused_bed_project(
    raw_bytes: &[u8],          // c * bytes_per_snp packed BED bytes
    bytes_per_snp: usize,
    n_indiv: usize,
    c: usize,                  // number of SNPs in this chunk
    d: usize,                  // sketch dimension
    proj: &MatF,               // P matrix (d × n_indiv), f64
    n: f64,                    // n_indiv as f64
    maf_out: &mut [f64],       // output: MAF per SNP (length c)
    out: &mut MatF,            // output: projected+renormalized (d × c)
    iid_positions: &[IidPos],  // precomputed byte/shift positions
    all_iids: bool,            // fast path flag
    lut: &[f64; 4],            // BED decode LUT (for f64)
) {
    // --- Pass 1: compute per-SNP statistics from packed bytes ---
    // Use a byte-level LUT: for each possible byte value, precompute
    // (sum_of_4_genotypes, count_of_non_missing_4, sum_sq_of_4_genotypes).
    // This is identical to compute_snp_stats but returns mean/inv_std.

    struct SnpStats {
        mean: f64,
        inv_std: f64,
        maf: f64,
        count: usize,
    }
    let stats: Vec<SnpStats> = (0..c).map(|j| {
        let snp_bytes = &raw_bytes[j * bytes_per_snp .. (j+1) * bytes_per_snp];
        // ... byte-level accumulation with LUT ...
        // Return mean, inv_std, maf (same math as normalize_col_f64_with_stats)
    }).collect();

    for (j, s) in stats.iter().enumerate() {
        maf_out[j] = s.maf;
    }

    // --- Pass 2: tiled decode + normalize + project ---
    const TILE: usize = 512;
    // Zero the output accumulator
    for j in 0..c {
        let col = out.col_mut(j).try_as_col_major_mut().unwrap().as_slice_mut();
        for v in &mut col[..d] { *v = 0.0; }
    }

    // Temporary tile buffer
    let mut b_tile = MatF::zeros(TILE, c);  // NOTE: allocated once, reused

    for tile_start in (0..n_indiv).step_by(TILE) {
        let tile_end = (tile_start + TILE).min(n_indiv);
        let tile_n = tile_end - tile_start;

        // Decode + normalize into b_tile[tile_n × c]
        if all_iids {
            for j in 0..c {
                let snp_bytes = &raw_bytes[j * bytes_per_snp .. (j+1) * bytes_per_snp];
                let mean = stats[j].mean;
                let inv_std = stats[j].inv_std;
                let col = b_tile.col_mut(j).try_as_col_major_mut().unwrap().as_slice_mut();
                for i in 0..tile_n {
                    let iid = tile_start + i;
                    let byte = snp_bytes[iid / 4];
                    let bits = (byte >> ((iid % 4) * 2)) & 0x3;
                    let geno = lut[bits as usize];
                    col[i] = if geno.is_nan() { 0.0 } else { (geno - mean) * inv_std };
                }
            }
        } else {
            // Subset path: use iid_positions[tile_start..tile_end]
            // But iid_positions indices are NOT contiguous w.r.t. individuals.
            // We need to iterate over iid_positions[tile_start..tile_end].
            for j in 0..c {
                let snp_bytes = &raw_bytes[j * bytes_per_snp .. (j+1) * bytes_per_snp];
                let mean = stats[j].mean;
                let inv_std = stats[j].inv_std;
                let col = b_tile.col_mut(j).try_as_col_major_mut().unwrap().as_slice_mut();
                for (i, pos) in iid_positions[tile_start..tile_end].iter().enumerate() {
                    let byte = snp_bytes[pos.byte_idx];
                    let bits = (byte >> pos.shift) & 0x3;
                    let geno = lut[bits as usize];
                    col[i] = if geno.is_nan() { 0.0 } else { (geno - mean) * inv_std };
                }
            }
        }

        // Project: out[d × c] += P[:, tile_start..tile_end] × b_tile[tile_n × c]
        let p_tile = mat_slice(proj.as_ref(), 0..d, tile_start..tile_end);
        let b_sl = mat_slice(b_tile.as_ref(), 0..tile_n, 0..c);
        matmul_to(
            mat_slice_mut(out.as_mut(), 0..d, 0..c),
            p_tile,
            b_sl,
            1.0,
            Accum::Add,  // accumulate across tiles
            Par::Seq,    // small GEMM per tile; outer parallelism TBD
        );
    }

    // Ratio estimator re-normalization
    for j in 0..c {
        let col = out.col(j).try_as_col_major().unwrap().as_slice();
        let mut nrm_sq = 0.0f64;
        for &v in &col[..d] { nrm_sq += v * v; }
        if nrm_sq > 0.0 {
            let scale = (n / nrm_sq).sqrt();
            let col_mut = out.col_mut(j).try_as_col_major_mut().unwrap().as_slice_mut();
            for v in &mut col_mut[..d] { *v *= scale; }
        }
    }
}
```

2. **Wire into the main loop**: Replace the sketch-mode block (lines 960-1130) with a conditional:

```rust
if sketching && use_fused_sketch {
    // Read raw bytes (already done by chunk_reader.read_raw_next or prefetch)
    // Call sketch_fused_bed_project(...)
    // Write results into bufs.b_mat
} else if sketching {
    // Existing path: decode → b_full → normalize → P×b_full → renorm
} else {
    // Existing exact path (UNTOUCHED)
}
```

3. **CLI flag**: No new flag needed — the fused path is activated automatically when `--sketch` is set. It produces identical output. Could add `--no-fused-sketch` escape hatch for debugging.

4. **f32 variant**: Duplicate the function for f32 (`sketch_fused_bed_project_f32`) or make it generic over a trait. Given the existing code's pattern of explicit f32/f64 branches, an explicit f32 variant is more consistent with the codebase style.

5. **Byte-level statistics LUT**: Pre-compute a 256-entry table mapping each possible BED byte to (sum, count, sum_sq) of the 4 genotypes it encodes. This is identical to what `compute_snp_stats` already does (lines 2140-2146). Factor out the shared LUT.

### 1.3 Adversarial Checklist

| Risk | Mitigation |
|------|------------|
| **Tile boundary misalignment**: last tile has tile_n < TILE | `tile_end = min(...)` and use `tile_n` consistently |
| **NaN handling**: missing genotypes must become 0.0 after normalization | Check `lut[bits]` for NaN/missing (BED code 0b01 = missing); set to 0.0 |
| **Statistics precision**: byte-LUT accumulation matches f32-based sum_sumsq_f32 | Both use exact integer arithmetic (genotypes are {0,1,2}); LUT sums are u32, converted to f64. Must match exactly. **Key difference**: current code calls `sum_sumsq_f32` on the decoded f32 column, which accumulates in f32 then converts to f64. The byte-LUT approach accumulates in u32/u64 → f64. These give **identical results** because BED values are {0,1,2,NaN}, all exact in f32. BUT: the NaN-count detection differs. Current: `sum.is_nan()`. Byte LUT: count < n_indiv. Must handle identically. |
| **b_tile allocation inside loop**: TILE=512, c=200, f64 → 800KB per tile, allocated once | Allocate `b_tile` once before the chunk loop, reuse. 800KB fits comfortably. |
| **iid_positions tiling with --keep**: iid_positions may not correspond to contiguous BED individuals | The tile iterates over `iid_positions[tile_start..tile_end]`, not over BED individuals. Each IidPos stores the byte_idx/shift for the actual BED position. This is correct. |
| **Accum::Add vs Accum::Replace**: first tile must be Add into zeros | Zero the output before the tile loop. All tiles use Accum::Add. |
| **P matrix column ordering**: P[:, i] corresponds to individual i in BED | When `--keep` is active, individuals are subsetted. P is d × n_indiv where n_indiv is the filtered count. iid_positions[i] maps to BED position of the i-th kept individual. P[:, i] is the projection weight for that same individual. The tile over `i = 0..n_indiv` correctly pairs P[:, i] with iid_positions[i]. |
| **Prefetch-bed path**: prefetch returns decoded f32 Mat, not raw bytes | Need to add a raw-bytes prefetch path, OR only use fused kernel with sequential reader. Simplest: fused kernel only activates with sequential reader (non-prefetch). With prefetch, fall back to current path. Document this. |
| **Off-by-one in bytes_per_snp**: `(n_iid + 3) / 4` must match | Use `bed.bytes_per_snp()` directly. |

### 1.4 Testing Plan

1. **Bit-exact parity with current sketch path**:
   ```bash
   # Current sketch path
   cargo build --release
   target/release/ldsc l2 --bfile data/bench_5k --out /tmp/sketch_current \
       --ld-wind-kb 1000 --sketch 50

   # Fused sketch path (same binary, fused is automatic)
   target/release/ldsc l2 --bfile data/bench_5k --out /tmp/sketch_fused \
       --ld-wind-kb 1000 --sketch 50

   # Compare: max_abs_diff must be 0.0
   diff <(zcat /tmp/sketch_current.*.l2.ldscore.gz) <(zcat /tmp/sketch_fused.*.l2.ldscore.gz)
   ```

2. **Test with --keep** (individual subset):
   ```bash
   head -100 data/bench_5k.fam > /tmp/keep100.txt
   target/release/ldsc l2 --bfile data/bench_5k --out /tmp/sketch_keep \
       --ld-wind-kb 1000 --sketch 50 --keep /tmp/keep100.txt
   ```

3. **Test with --extract** (SNP subset):
   ```bash
   head -1000 data/bench_5k.bim | awk '{print $2}' > /tmp/extract1k.txt
   target/release/ldsc l2 --bfile data/bench_5k --out /tmp/sketch_extract \
       --ld-wind-kb 1000 --sketch 50 --extract /tmp/extract1k.txt
   ```

4. **Test with --annot** (partitioned):
   - Sketch + partitioned LD scores with annotation file.

5. **Test with --fast-f32**:
   - Verify f32 fused path matches f32 current path.

6. **Edge cases**: chunk where c < chunk_c (last chunk), monomorphic SNPs (inv_std=0), all-missing SNPs.

### 1.5 Benchmarking Plan

**Local quick check**:
```bash
hyperfine --warmup 1 \
    'target/release/ldsc l2 --bfile data/bench_200k --out /tmp/bench --ld-wind-kb 1000 --sketch 100' \
    --min-runs 5
```

**AWS HPC** (via `scripts/aws-bench.sh`):
- Full 1.66M SNPs, N=2490, `--sketch 100`: compare fused vs current.
- Biobank 50K, `--sketch 100`: this is where the memory traffic reduction matters most.
- Report: wall time, user time (parallelism ratio), and `--verbose-timing` breakdown.

**Key metric**: At N=50K, sketch d=100, the current path's bottleneck is the P×b_full GEMM. The fused path should show:
- Reduced `t_sketch` (no separate GEMM)
- Increased `t_norm` (fused into decode)
- Net wall time reduction proportional to memory traffic saved

---

## Optimization 2: Individual Subsampling (`--subsample N'`)

### 2.1 Verification: Statistical Correctness

LD scores measure pairwise correlation (r²) between SNPs across individuals. For common variants (MAF > 5%), the sample r² converges as O(1/N) in variance. At N=2490 (1000G), LD scores are already well-estimated. At N=50K (biobank), there is massive redundancy — N'=5000-10000 individuals give nearly identical LD scores for the typical MAF spectrum.

**Mathematical argument**: For SNPs i,j with population correlation ρ:
- `E[r²] = ρ² + (1-ρ²)/N + O(1/N²)` (the (1-ρ²)/N term is the r²_unbiased correction)
- `Var[r²] ≈ 4ρ²(1-ρ²)²/N + O(1/N²)` for the dominant term
- At N=50K: `Var[r²] ≈ 8×10⁻⁵ × ρ²(1-ρ²)²`
- At N'=5K: `Var[r²] ≈ 8×10⁻⁴ × ρ²(1-ρ²)²` — 10x higher but still small
- LD scores sum ~1000 such r² terms, so SE scales as `sqrt(1000) × sqrt(Var) ≈ 0.9%` relative error at N'=5K.

**Performance impact**: GEMM is O(N × c²) per chunk. Reducing N from 50K to 5K gives 10x speedup in GEMM. BED I/O also reduces ~10x (via --keep subsetting). Total speedup: ~10x for exact computation at N'=5K. This makes exact-at-N'=5K competitive with sketch-d=100-at-N=50K, but with better accuracy.

### 2.2 Implementation Steps

**File: `src/cli.rs`** — Add flag:

```rust
/// Randomly subsample N' individuals from the reference panel for LD computation.
/// Reduces both I/O and GEMM cost proportionally. For biobank-scale data (N>10K),
/// N'=5000-10000 gives accurate LD scores for common variants (MAF>5%).
/// Overridden by --keep (which provides explicit control over which individuals).
#[arg(long, value_name = "N_PRIME", conflicts_with = "keep")]
pub subsample: Option<usize>,
```

**File: `src/l2.rs`** — In the `run_l2` function (or wherever `iid_indices` is constructed):

```rust
// If --subsample is set and no --keep, generate random subset
let iid_indices: Option<Vec<isize>> = if let Some(n_prime) = args.subsample {
    let n_total = /* iid_count from FAM */;
    if n_prime >= n_total {
        eprintln!("WARNING: --subsample {} >= n_indiv {}; using all individuals", n_prime, n_total);
        None
    } else {
        // Fisher-Yates partial shuffle to select n_prime indices
        let mut rng = fastrand::Rng::with_seed(42);  // deterministic
        let mut indices: Vec<usize> = (0..n_total).collect();
        for i in 0..n_prime {
            let j = rng.usize(i..n_total);
            indices.swap(i, j);
        }
        indices.truncate(n_prime);
        indices.sort_unstable();  // sort for sequential BED access
        let iid_isize: Vec<isize> = indices.iter().map(|&i| i as isize).collect();
        Some(iid_isize)
    }
} else if let Some(ref keep_file) = args.keep {
    // Existing --keep logic
    ...
} else {
    None
};
```

**Key detail**: Sort the subsampled indices so BED reads remain sequential-friendly. The `ChunkReader` with `iid_positions` already handles arbitrary subsets.

**File: `src/l2.rs`** — No changes to `compute_ldscore_global` needed. It already accepts `iid_indices: Option<&[isize]>`.

### 2.3 Adversarial Checklist

| Risk | Mitigation |
|------|------------|
| **Seed determinism**: results must be reproducible | Use `fastrand::Rng::with_seed(42)`, same as sketch projection |
| **n_prime > n_indiv**: must not panic | Guard with >= check, warn and proceed with all individuals |
| **n_prime = 0**: degenerate case | Reject: `anyhow::ensure!(n_prime > 0, "--subsample must be > 0")` |
| **Interaction with --keep**: both specify individuals | `conflicts_with = "keep"` in clap |
| **r2_unbiased correction uses wrong N**: must use n_prime, not n_total | `n_indiv` in `compute_ldscore_global` is already derived from the actual number of individuals used (after --keep filtering). The subsample indices are passed as `iid_indices`, so `n_indiv` will correctly be `n_prime`. Verify this in the caller. |
| **M/M_5_50 files**: should reflect the reference panel, not the subsample | M files count SNPs, not individuals. No change needed. |
| **MAF computation**: uses subsampled individuals | This is correct — MAF should be computed from the same individuals used for LD. The subsampled MAF will be noisier but unbiased. |

### 2.4 Testing Plan

1. **Sanity check**: `--subsample N` where N = n_indiv should give identical output to no flag.
2. **Parity with --keep**: generate a keep file with the same random individuals, verify identical output.
   ```bash
   # Generate the same random subset as --subsample would
   # (need a small script that uses the same Fisher-Yates with seed 42)
   target/release/ldsc l2 --bfile data/bench_5k --out /tmp/sub500 \
       --ld-wind-kb 1000 --subsample 500
   # Compare with --keep using same individuals
   target/release/ldsc l2 --bfile data/bench_5k --out /tmp/keep500 \
       --ld-wind-kb 1000 --keep /tmp/keep500.txt
   diff <(zcat /tmp/sub500.*.l2.ldscore.gz) <(zcat /tmp/keep500.*.l2.ldscore.gz)
   ```
3. **Accuracy at N'=5K on biobank_50k**: compare LD scores at N'=5K vs exact at N=50K. Report Pearson r and median absolute error.

### 2.5 Benchmarking Plan

**AWS HPC**, biobank_50k dataset:
```bash
hyperfine \
    'ldsc l2 --bfile biobank_50k --out /tmp/exact --ld-wind-kb 1000' \
    'ldsc l2 --bfile biobank_50k --out /tmp/sub5k --ld-wind-kb 1000 --subsample 5000' \
    'ldsc l2 --bfile biobank_50k --out /tmp/sub10k --ld-wind-kb 1000 --subsample 10000' \
    'ldsc l2 --bfile biobank_50k --out /tmp/sketch100 --ld-wind-kb 1000 --sketch 100' \
    --min-runs 3
```

Expected: `--subsample 5000` should be ~10x faster than exact, comparable to `--sketch 100`, but with much better accuracy.

---

## Optimization 3: CountSketch Projection (`--sketch-method countsketch`)

### 3.1 Verification: Mathematical Correctness

**CountSketch**: Each of N individuals is assigned to one of d buckets via hash function h(i) ∈ {0,...,d-1}, and a random sign σ(i) ∈ {+1,-1}. The sketch of column x is:

```
x̃[k] = Σ_{i: h(i)=k} σ(i) × x[i]    for k = 0, ..., d-1
```

**JL property**: For any two columns x, y:
- `E[<x̃, ỹ>] = <x, y>` (unbiased)
- `Var[<x̃, ỹ>] = (1/d) × (||x||²||y||² + <x,y>² - Σ_i x_i² y_i²)`

For normalized columns (||x|| = ||y|| = sqrt(N)):
- Leading variance term: `(N² + <x,y>²) / d`
- For Rademacher projection: `(N² + <x,y>²) / d` (identical leading term!)
- But CountSketch has a higher-order correction: `- Σ_i x_i² y_i² / d` which is negative (reduces variance slightly for non-uniform entries).

**Key difference from Rademacher**: CountSketch is strictly sparser — each column of P has exactly one nonzero. This means:
- No P matrix in memory (just hash function + signs, O(N) storage vs O(d×N) for dense P).
- Application is O(N×c) instead of O(d×N×c) for dense GEMM. For d=100, N=50K: 5M vs 1B operations.
- **But**: the scatter-add pattern is memory-unfriendly (random writes to d buckets). At d=100-200 the output fits L1, so this is fine.

**Bias correction**: Same ratio estimator applies. After projection, rescale each column so ||x̃'|| = sqrt(N). Then apply the same d/(d-2) correction.

**Accuracy**: At same d, CountSketch has slightly higher variance than Rademacher (the `Σ_i x_i² y_i²` correction favors Rademacher for non-uniform entries). Recommendation: use d=100-200 with CountSketch for equivalent quality to d=50 Rademacher.

### 3.2 Implementation Steps

**File: `src/cli.rs`** — Add flag:

```rust
/// Sketch projection method: "rademacher" (default, dense random ±1/√d) or
/// "countsketch" (hash-based, O(N) per column instead of O(d×N)).
/// CountSketch is faster but slightly noisier at same d; use d=100-200 for
/// equivalent quality to Rademacher d=50.
#[arg(long, value_name = "METHOD", default_value = "rademacher")]
pub sketch_method: String,
```

**File: `src/l2.rs`** — New data structure and projection function.

```rust
/// CountSketch projection: hash each individual to a bucket with a random sign.
struct CountSketchProj {
    bucket: Vec<u32>,    // h(i) for each individual, length n_indiv
    sign: Vec<f64>,      // σ(i) ∈ {+1, -1} for each individual, length n_indiv
    d: usize,
}

impl CountSketchProj {
    fn new(n_indiv: usize, d: usize, seed: u64) -> Self {
        let mut rng = fastrand::Rng::with_seed(seed);
        let bucket: Vec<u32> = (0..n_indiv).map(|_| rng.u32(..d as u32)).collect();
        let sign: Vec<f64> = (0..n_indiv).map(|_| {
            if rng.bool() { 1.0 } else { -1.0 }
        }).collect();
        Self { bucket, sign, d }
    }

    /// Project a column x[0..n] into sketch[0..d] via scatter-add.
    /// sketch must be zeroed before calling.
    #[inline]
    fn project_col(&self, x: &[f64], sketch: &mut [f64]) {
        debug_assert!(sketch.len() >= self.d);
        for (i, &val) in x.iter().enumerate() {
            let b = self.bucket[i] as usize;
            sketch[b] += self.sign[i] * val;
        }
    }

    /// Project multiple columns: out[d × c] from x[n × c].
    fn project_matrix(&self, x: &MatF, n: usize, c: usize, out: &mut MatF) {
        for j in 0..c {
            let x_col = x.col(j).try_as_col_major().unwrap().as_slice();
            let out_col = out.col_mut(j).try_as_col_major_mut().unwrap().as_slice_mut();
            // Zero
            for v in &mut out_col[..self.d] { *v = 0.0; }
            self.project_col(&x_col[..n], out_col);
        }
    }
}
```

**Integration into the sketch setup block (lines 524-611)**:

```rust
enum SketchMethod {
    Rademacher { proj_f64: Option<MatF>, proj_f32: Option<MatF32> },
    CountSketch { cs: CountSketchProj },
}

// Replace proj_f64/proj_f32 with SketchMethod enum.
let sketch_method = if let Some(d) = sketch_dim {
    if args.sketch_method == "countsketch" {
        SketchMethod::CountSketch { cs: CountSketchProj::new(n_indiv, d, 42) }
    } else {
        // Existing Rademacher setup...
        SketchMethod::Rademacher { proj_f64, proj_f32 }
    }
} else {
    // No sketch
    ...
};
```

**Integration into the projection block (lines 1060-1130)**:

```rust
if sketching {
    match &sketch_method {
        SketchMethod::Rademacher { .. } => {
            // Existing P × b_full GEMM path
        }
        SketchMethod::CountSketch { cs } => {
            // b_full is already normalized (same as current path)
            match bufs {
                GemmBufs::F64 { ref mut b_mat, .. } => {
                    cs.project_matrix(&b_full_f64, n_indiv, c, b_mat);
                    // Ratio estimator re-normalization (same as Rademacher)
                    for j in 0..c {
                        // ... same renorm code ...
                    }
                }
                // f32 variant
            }
        }
    }
}
```

**Memory savings**: No d×N projection matrix. At d=100, N=50K: saves 100×50K×8 = 40MB (f64) or 20MB (f32). The CountSketch projection itself is O(N×c) = O(50K×200) = 10M operations per chunk, compared to 1B for the dense GEMM. **~100x fewer FLOPs for the projection step**.

**However**, the projection is memory-bound (sequential reads of x, random writes to sketch). At N=50K, c=200: we read 80MB (b_full f64) and write 160KB (d×c output). The bottleneck shifts entirely to the BED read + normalization, which is exactly what Optimization 1 addresses.

### 3.3 Adversarial Checklist

| Risk | Mitigation |
|------|------------|
| **Hash quality**: fastrand may produce correlated bucket assignments | fastrand's WyRand is high-quality for this use case. 2-universal hashing would be theoretically cleaner but unnecessary in practice. |
| **f32 variant**: need CountSketch with f32 sign/accumulation | Add `project_col_f32` with `sign: Vec<f32>`. Accumulate in f32 (genotypes are small integers, d=100-200 buckets, no overflow risk). |
| **Interaction with fused kernel (Opt 1)**: CountSketch can be fused with BED decode | Initial implementation: use b_full buffer (same as Rademacher). Future: fuse CountSketch with tile decode (even simpler than fusing Rademacher, since CountSketch is just a scatter-add per element). |
| **Bias correction constants**: same d/(d-2) formula? | Yes. The ratio estimator corrects for the sketch dimension d regardless of the projection method. The leading bias of CountSketch is the same as Rademacher. |
| **Deterministic seeding**: must be reproducible | Use seed=42, same as Rademacher. But use a DIFFERENT seed from the Rademacher P matrix to avoid any correlation if both are ever compared. Use seed=43 for CountSketch. Actually, just use 42 — they're never used simultaneously. |
| **Bucket collisions**: at d=100, N=50K, expect 500 individuals per bucket | This is fine. CountSketch's JL guarantee holds regardless of collision rate (it's in the variance formula). |
| **--keep interaction**: CountSketch indices must match subsetted individuals | `n_indiv` passed to `CountSketchProj::new` is already the filtered count. project_col iterates over x[0..n_indiv], which corresponds to the subsetted individuals. Correct. |
| **Zero-variance SNPs**: projection of all-zeros column | Results in all-zero sketch. nrm_sq=0, no rescaling. Correct — matches Rademacher behavior. |

### 3.4 Testing Plan

1. **Statistical correctness**: Compare CountSketch LD scores vs exact at various d.
   ```bash
   for d in 50 100 200; do
       target/release/ldsc l2 --bfile data/bench_200k --out /tmp/cs_d${d} \
           --ld-wind-kb 1000 --sketch $d --sketch-method countsketch
       target/release/ldsc l2 --bfile data/bench_200k --out /tmp/rad_d${d} \
           --ld-wind-kb 1000 --sketch $d --sketch-method rademacher
   done
   # Compare both against exact; report Pearson r, median abs error
   ```

2. **Verify unbiasedness**: Average over multiple seeds (if we add a --sketch-seed flag).

3. **Edge cases**: d=3 (minimum), d=n_indiv-1 (maximum), monomorphic SNPs, all-missing SNPs.

### 3.5 Benchmarking Plan

**AWS HPC**, full 1.66M and biobank_50k:
```bash
# 1000G (N=2490)
hyperfine \
    'ldsc l2 --bfile 1000G --out /tmp/rad100 --ld-wind-kb 1000 --sketch 100' \
    'ldsc l2 --bfile 1000G --out /tmp/cs100 --ld-wind-kb 1000 --sketch 100 --sketch-method countsketch' \
    --min-runs 5

# Biobank (N=50K) — this is where CountSketch shines
hyperfine \
    'ldsc l2 --bfile biobank_50k --out /tmp/exact --ld-wind-kb 1000' \
    'ldsc l2 --bfile biobank_50k --out /tmp/rad100 --ld-wind-kb 1000 --sketch 100' \
    'ldsc l2 --bfile biobank_50k --out /tmp/cs100 --ld-wind-kb 1000 --sketch 100 --sketch-method countsketch' \
    --min-runs 3
```

**Expected results**:
- At N=2490: CountSketch slightly faster than Rademacher (no d×N GEMM), but the GEMM is small (100×2490), so marginal improvement.
- At N=50K: CountSketch projection is ~100x faster than Rademacher GEMM. If BED read + norm dominates, total speedup is modest (~10-20%). Combined with Optimization 1 (fused), the entire sketch path becomes BED-read-bound.

---

## Summary: Implementation Order and Commit Plan

### Commit 1: feat: fused BED-decode-normalize-project for sketch mode
- Files changed: `src/bed.rs` (expose raw bytes), `src/l2.rs` (new fused function + wiring)
- Parity: bit-exact with current sketch path (max_abs_diff=0)
- Risk: medium (new code path, but mathematically identical)

### Commit 2: feat: `--subsample N'` for individual subsampling
- Files changed: `src/cli.rs` (new flag), `src/l2.rs` (subsample logic in run_l2)
- Parity: not applicable (different individuals = different results, by design)
- Risk: low (thin wrapper around existing --keep infrastructure)

### Commit 3: feat: `--sketch-method countsketch` hash-based projection
- Files changed: `src/cli.rs` (new flag), `src/l2.rs` (CountSketchProj struct + integration)
- Parity: not applicable (different projection = different results, by design)
- Risk: low-medium (new projection method, but same bias correction framework)

### Optional Follow-up: feat: fused CountSketch (combines Opt 1 + Opt 3)
- Fuse CountSketch scatter-add directly into the BED decode tile loop.
- Eliminates b_full buffer entirely for CountSketch mode.
- Reads 2.5MB packed BED + O(N) hash/sign arrays. Total: ~2.9MB vs 84MB current.
- This is the ultimate memory-traffic-optimal path.

---

## Open Questions

1. **TILE size**: 512 is a guess. Need to benchmark 256, 512, 1024 on AWS. Too small = GEMM overhead. Too large = exceeds L2. Target: b_tile + P_tile < 1MB (L2 per core).

2. **Par::Seq vs Par::rayon(0) for tile GEMMs**: Each tile GEMM is d×TILE×c ≈ 100×512×200 = 10M FLOPs. This is borderline for rayon overhead. Benchmark both.

3. **Parallel tile accumulation**: Can we `par_iter` over tiles with thread-local accumulators? Adds `d×c×8B×n_threads` memory (320KB×16 = 5MB). Need to benchmark whether the parallelism gain outweighs the final reduction cost.

4. **CountSketch + f32**: Should the scatter-add accumulate in f32 or f64? At d=100, ~500 values per bucket, max value ~2 (genotype), accumulated sum ~1000. f32 can represent integers up to 2^23 exactly, so f32 accumulation is safe. But after normalization, values are fractional — f32 accumulation may lose precision. Benchmark both.

5. **Seed strategy**: CountSketch and Rademacher use the same seed (42). If we ever want to compare them in the same run, they should use different seeds. For now, only one is active at a time, so this is fine.
