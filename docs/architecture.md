# Architecture & Source Map

## Source Map

```
src/
├── main.rs          Clap dispatch — parses CLI, calls into subcommand modules.
│
├── cli.rs           All argument structs (MungeArgs, L2Args, H2Args, RgArgs,
│                    MakeAnnotArgs, CtsAnnotArgs). No logic — pure clap derive macros.
│
├── gpu.rs           [cfg(feature = "gpu")] CUDA matmul via CubeCL/cubek-matmul.
│                    · GpuContext — holds ComputeClient + Strategy
│                    · matmul_tn — A^T × B single-shot GPU GEMM
│                    · matmul_tn_tiled — tiled variant for VRAM-limited windows
│
├── parse.rs         File I/O helpers:
│                    · scan_sumstats / scan_ldscore  → Polars LazyFrame
│                    · concat_chrs(prefix, suffix)   → concat per-chr files
│                    · read_m_total / read_m_vec      → .l2.M files
│                    · read_annot(prefix, thin)       → MatF + col names
│
├── munge.rs         munge-sumstats pipeline (Polars LazyFrame, no data loaded until
│                    collect). Internal functions:
│                    · apply_ignore, apply_col_overrides, normalize_columns
│                    · apply_info_list, apply_n_override
│                    · derive_z (BETA/SE → Z; P + sign → Z; --a1-inc)
│                    · filter_snps, apply_nstudy_filter
│                    · write_sumstats_gz (gzip TSV output)
│
├── l2/              LD score computation (split into submodules):
│   ├── mod.rs       · run — orchestrates BIM read, --extract / --annot / --keep,
│   │                  calls compute_ldscore_global, writes per-chr .l2.ldscore.gz
│   │                  and .l2.M / .l2.M_5_50 files
│   ├── compute.rs   · GemmBufs — enum holding f32 or f64 scratch buffers (--fast-f32)
│   │                · compute_ldscore_global — ring-buffer GEMM loop (sequential,
│   │                  scalar and partitioned + sketch + stochastic paths)
│   │                · r2_unbiased — r² − (1−r²)/(n−2)
│   ├── window.rs    · WindowMode — Cm / Kb / Snp enum
│   │                · get_block_lefts_f64, get_block_lefts_by_chr — window boundaries
│   ├── normalize.rs · normalize_col_f{32,64}_with_stats — impute NaN → mean,
│   │                  centre, unit-variance; AVX2+FMA sum_sumsq_f32
│   ├── snp_stats.rs · compute_snp_stats — fast BED scan for MAF + het/missing
│   └── io.rs        · parse_bim, count_fam, parse_fam — PLINK file parsers
│                    · load_individual_indices — --keep FID/IID → isize indices
│                    · write_ldscore_refs — gzip TSV output
│                    · load_snp_set — HashSet<String> from --extract / --print-snps
│
├── irwls.rs         Iteratively Re-Weighted Least Squares.
│                    · IrwlsResult — est + optional jackknife fields
│                    · irwls(x, y, weights, n_iter) — pre-alloc xw/yw, SVD solve,
│                      reweight on fitted values; zero-alloc inner loop
│
├── jackknife.rs     Block jackknife variance estimation.
│                    · jackknife(x, y, weights, n_blocks, n_iter) →
│                      full-data IRWLS + n_blocks parallel leave-one-out refits
│                      (rayon par_iter) → pseudo-values → SE + covariance matrix
│
└── regressions.rs   h2 and rg regression drivers.
                     · run_h2 — loads sumstats + LD scores, inner-joins on SNP,
                       detects K annotation columns, dispatches to scalar (K=1)
                       or partitioned (K>1) path; supports --two-step, --no-intercept,
                       --intercept-h2, --print-coefficients, liability-scale output
                     · run_h2_partitioned — K-column design matrix, per-annotation
                       enrichment, resolves M vector from per-annotation M files
                     · run_h2_scalar — shared by standalone h2 and rg univariate sub-fits
                     · run_rg — iterates trait pairs; gencov regression + univariate h2
                       per trait; --two-step, --intercept-gencov, --no-intercept,
                       liability-scale rg; prints summary table
                     · load_ld_ref / load_ld — LazyFrame readers for ref and weight LD
                     · resolve_m / resolve_m_vec — reads .l2.M_5_50 or falls back to n_obs
                     · liability_scale_h2 — observed → liability scale conversion
                     · print_jackknife_diagnostics — --print-cov / --print-delete-vals

make_annot.rs        BED → 0/1 annotation generator.
                     · annotate_from_bed — loads BED intervals per chromosome,
                       sorts and merges, binary-search annotation per SNP
                     · annotate_from_gene_set — gene symbols → coordinate lookup →
                       same interval merge/annotate pipeline
                     · write_annot_file — CHR BP SNP CM ANNOT TSV (optional .gz)
```

## Key data-flow invariants

- `l2 --annot prefix` reads `{prefix}{chr}.annot[.gz]` for every chromosome found in
  the BIM (not a single `prefix.annot.gz` file).
- `--extract` filters the BIM *before* window computation; `--print-snps` filters only the output.
- `bed_idx` (original BIM row index) differs from `pos` (index in the filtered `all_snps` slice)
  when `--extract` is active; `bed_idx_to_pos` in `run()` maps between them.
- `--keep` passes `iid_indices: Option<&[isize]>` to the internal BED reader; `n_indiv_actual` (not
  the FAM total) is used for normalization and the r²-unbiased correction.
- The ring buffer `ring_size = max_window + chunk_c` guarantees no live window slot is overwritten
  before it has been consumed in the A×B product.
- `rg --no-intercept` propagates `fixed_intercept = Some(1.0)` to both the gencov regression and
  each univariate `run_h2_scalar` call, matching Python's behaviour.

## Dependency rationale

| Crate | Version | Role |
|-------|---------|------|
| internal `bed` module | - | minimal PLINK .bed reader tailored to LDSC |
| `polars` | 0.53 | lazy CSV streaming (munge + LD score file loading) |
| `faer` | 0.24 | dense matrix algebra; SVD for IRWLS |
| `rayon` | 1 | data-parallel jackknife blocks |
| `statrs` | 0.18 | Normal CDF/quantile for P→Z conversion |
| `clap` | 4 | derive-macro CLI argument parsing |
| `anyhow` | 1 | error propagation |
| `flate2` | 1 | gzip output for .sumstats.gz and .ldscore.gz |
| `cubecl` | 0.10 | (optional, `gpu` feature) multi-backend GPU compute |
| `cubek-matmul` | 0.2 | (optional, `gpu` feature) autotuned GPU matmul |
| `fastrand` | 2 | Rademacher random generation for `--stochastic` and `--sketch` modes |
| `mimalloc` | 0.1 | (optional, `mimalloc` feature) fast allocator for musl builds |
