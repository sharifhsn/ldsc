# Performance Log

This log captures **post-change** performance measurements and parity checks.
Each entry should list the dataset, command, and key timings so we can track
how changes affect runtime.

## 2026-03-03
- Change: low-risk l2 optimizations (buffer reuse, contiguous BED indexing, `general_mat_mul`).
- Dataset: 1000G.EUR.QC, extract 500k SNPs (`/Users/sharif/Code/ldsc/perf/l2/extract_500000.snps`).
- Command: `ldsc l2 --bfile /Users/sharif/Code/ldsc/data/1000G.EUR.QC --out /Users/sharif/Code/ldsc/perf/l2/rust_l2_after --ld-wind-cm 1 --extract /Users/sharif/Code/ldsc/perf/l2/extract_500000.snps`
- Rust timing: `real 22.39s` (from `/Users/sharif/Code/ldsc/perf/l2/rust_l2_after.time`).
- Previous baseline (pre-opt): `real 27.19s` (from `/Users/sharif/Code/ldsc/perf/l2/rust_l2.time`).
- Speedup: ~17.7%.
- Parity check (Python vs Rust, 500k extract): **FAILED**.
  - Script: `EXTRACT_N=500000 /Users/sharif/Code/ldsc/scripts/verify_parity_l2.sh`
  - Result: `max_abs_diff=0.001` at `rs4275489` (tolerance 1e-3 exceeded).

## 2026-03-03 (later)
- Change: reverted low-risk l2 optimizations (removed buffer reuse, contiguous BED index fast path, and `general_mat_mul` usage).
- Parity check (Python vs Rust, 500k extract): **FAILED**.
  - Script: `EXTRACT_N=500000 /Users/sharif/Code/ldsc/scripts/verify_parity_l2.sh`
  - Result: `max_abs_diff=0.001` at `rs4275489` (tolerance 1e-3 exceeded).
  - Implication: parity drift persists without those optimizations.

## 2026-03-03 (later, option D)
- Change: added f64 default normalization with opt-in `--fast-f32`.
- Dataset: 1000G.EUR.QC (full dataset, no extract).
- Command (f64 default): `ldsc l2 --bfile /Users/sharif/Code/ldsc/data/1000G.EUR.QC --out /Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64 --ld-wind-cm 1`
- Rust timing (f64): `real 298.17s` (from `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64.time`).
- Command (fast f32): `ldsc l2 --bfile /Users/sharif/Code/ldsc/data/1000G.EUR.QC --out /Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f32 --ld-wind-cm 1 --fast-f32`
- Rust timing (fast f32): `real 329.21s` (from `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f32.time`).
- Result: `--fast-f32` is ~10.4% slower on this run.

## 2026-03-03 (baseline: f64-only)
- Change: removed `--fast-f32`; L2 normalization is f64-only.
- Dataset: 1000G.EUR.QC (full dataset, no extract).
- Command: `ldsc l2 --bfile /Users/sharif/Code/ldsc/data/1000G.EUR.QC --out /Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace --ld-wind-cm 1`
- Trace log: `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.stdout`
- Timing (from trace): `maf_prefilter=52.005s`, `compute_ldscore=215.833s` with
  `bed_read=54.602s`, `norm=11.225s`, `bb_dot=9.473s`, `ab_dot=93.586s`, `r2u=17.356s`,
  and `write_outputs=6.830s`.
- Wall time: `real 276.58s` (from `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.stderr`).
- Output fingerprints: `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.sha256` (69 files).

## 2026-03-03 (ab_dot optimizations)
- Change: preallocated `ab` GEMM buffer, contiguous ring-buffer fast path for `A` windows, and precomputed pq weights.
- Dataset: 1000G.EUR.QC (full dataset, no extract).
- Command: `ldsc l2 --bfile /Users/sharif/Code/ldsc/data/1000G.EUR.QC --out /Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace --ld-wind-cm 1`
- Trace log: `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.stdout`
- Timing (from trace): `maf_prefilter=49.605s`, `compute_ldscore=207.446s` with
  `bed_read=54.576s`, `norm=11.607s`, `bb_dot=9.268s`, `ab_dot=90.946s`, `r2u=34.483s`,
  and `write_outputs=7.102s`.
- Wall time: `real 267.37s` (from `/Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.stderr`).
- Output parity: `shasum -a 256 -c /Users/sharif/Code/ldsc/perf/l2/rust_l2_full_f64_trace.sha256` → OK.

## 2026-03-04 (baseline for kb-window perf work)
- Change: none (baseline before ring-buffer slack tweak).
- Dataset: 1000G_phase3_common_norel (full dataset, no extract).
- Command: `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_full_kb1000 --ld-wind-kb 1000 --chunk-size 50`
- Trace log: `/home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_trace_kb1000.log`
- Timing (from trace): `compute_ldscore=63.264s` with
  `bed_read=5.117s`, `norm=9.261s`, `bb_dot=17.772s`, `ab_dot=37.089s`, `r2u_ab=3.788s`,
  `ab_wrap_copy=9.507s`, `ab_contig_matmul=28.384s`, `ab_wrap_matmul=8.706s`.

## 2026-03-04 (ring-buffer slack: max_window + 2*chunk)
- Change: increased ring buffer size to reduce wrap copies.
- Dataset: 1000G_phase3_common_norel (extract 200k SNPs).
- Command: `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_small_kb1000 --ld-wind-kb 1000 --chunk-size 50 --extract /home/sharif/Code/ldsc/perf/small_l2_trace/extract.snps`
- Trace logs: `/home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_trace_kb1000.run{1,2,3}.log` (run0 warmup).
- Timing (avg of runs 1-3): `compute_ldscore=11.740s` with
  `bed_read=0.617s`, `norm=1.105s`, `bb_dot=2.111s`, `ab_dot=4.251s`, `r2u_ab=0.438s`,
  `ab_wrap_copy=2.368s`, `ab_contig_matmul=1.983s`, `ab_wrap_matmul=2.268s`.

## 2026-03-04 (Par::Seq heuristic for small matmuls)
- Change: use `Par::Seq` for small matmuls based on output size (rolled back after regression).
- Dataset: 1000G_phase3_common_norel (extract 200k SNPs).
- Command: `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_small_kb1000 --ld-wind-kb 1000 --chunk-size 50 --extract /home/sharif/Code/ldsc/perf/small_l2_trace/extract.snps`
- Trace logs: `/home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_trace_kb1000_parheur.run{1,2,3}.log` (run0 warmup).
- Timing (avg of runs 1-3): `compute_ldscore=18.045s` with
  `bed_read=0.567s`, `norm=0.977s`, `bb_dot=0.944s`, `ab_dot=12.171s`, `r2u_ab=0.391s`,
  `ab_wrap_copy=2.219s`, `ab_contig_matmul=5.404s`, `ab_wrap_matmul=6.767s`.
- Notes: regression (~54% slower compute_ldscore); heuristic removed.

## 2026-03-04 (annot path r2u on-the-fly, single run)
- Change: compute r2u_ab then accumulate annot contributions via explicit loops (no matmul on r2u_ab).
- Dataset: 1000G_phase3_common_norel (extract 200k SNPs).
- Command: `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_small_kb1000 --ld-wind-kb 1000 --chunk-size 50 --extract /home/sharif/Code/ldsc/perf/small_l2_trace/extract.snps`
- Trace log: `/home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_trace_kb1000_r2u_onthefly.log`
- Timing (single run): `compute_ldscore=12.278s` with
  `bed_read=0.636s`, `norm=1.150s`, `bb_dot=2.165s`, `ab_dot=4.534s`, `r2u_ab=0.444s`,
  `ab_wrap_copy=2.494s`, `ab_contig_matmul=2.167s`, `ab_wrap_matmul=2.367s`.

## 2026-03-04 (full 1000G kb-window, single run after revert)
- Change: reverted annot r2u on-the-fly; baseline ring-buffer slack retained.
- Dataset: 1000G_phase3_common_norel (full dataset, no extract).
- Command: `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_full_kb1000 --ld-wind-kb 1000 --chunk-size 50`
- Trace log: `/home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_trace_kb1000.log`
- Timing (single run): `compute_ldscore=89.029s` with
  `bed_read=5.062s`, `norm=9.274s`, `bb_dot=17.341s`, `ab_dot=37.292s`, `r2u_ab=3.790s`,
  `ab_wrap_copy=9.207s`, `ab_contig_matmul=28.459s`, `ab_wrap_matmul=8.833s`,
  `write_outputs=2.599s`.
- Aggregate ab_chunk: `n=33297`, `contig=77.4%`, `avg_w=694.99`, `avg_c=50.00`,
  `avg_ab_ms=1.120`, `avg_r2u_ms=0.114`.

## 2026-03-04 (full 1000G ring slack sweep, single runs)
- Dataset: 1000G_phase3_common_norel (full dataset, no extract), `--ld-wind-kb 1000`, `--chunk-size 50`.
- ring +8c:
  - Trace log: `/home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_trace_kb1000_ring8.log`
  - `compute_ldscore=89.151s`
  - `ab_wrap_copy=8.412s`, `ab_contig_matmul=28.936s`, `ab_wrap_matmul=8.061s`
  - `ab_contig=26563`, `ab_wrap=6734`
- ring +16c:
  - Trace log: `/home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_trace_kb1000_ring16.log`
  - `compute_ldscore=89.813s`
  - `ab_wrap_copy=7.700s`, `ab_contig_matmul=30.468s`, `ab_wrap_matmul=7.551s`
  - `ab_contig=27266`, `ab_wrap=6031`
- ring +32c:
  - Trace log: `/home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_trace_kb1000_ring32.log`
  - `compute_ldscore=89.215s`
  - `ab_wrap_copy=6.362s`, `ab_contig_matmul=32.299s`, `ab_wrap_matmul=5.996s`
  - `ab_contig=28402`, `ab_wrap=4895`
- ring +4c:
  - Trace log: `/home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_trace_kb1000_ring4.log`
  - `compute_ldscore=89.653s`
  - `ab_wrap_copy=8.982s`, `ab_contig_matmul=29.132s`, `ab_wrap_matmul=8.362s`
  - `ab_contig=26035`, `ab_wrap=7262`

## 2026-03-04 (ring-buffer slack: max_window + 8*chunk, single run)
- Change: increased ring buffer size to reduce wrap copies further.
- Dataset: 1000G_phase3_common_norel (extract 200k SNPs).
- Command: `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_small_kb1000 --ld-wind-kb 1000 --chunk-size 50 --extract /home/sharif/Code/ldsc/perf/small_l2_trace/extract.snps`
- Trace log: `/home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_trace_kb1000_ring8.log`
- Timing (single run): `compute_ldscore=11.293s` with
  `bed_read=0.629s`, `norm=1.112s`, `bb_dot=2.146s`, `ab_dot=4.260s`, `r2u_ab=0.441s`,
  `ab_wrap_copy=1.855s`, `ab_contig_matmul=2.527s`, `ab_wrap_matmul=1.733s`.

## 2026-03-04 (ring-buffer slack: max_window + 16*chunk, single run)
- Change: increased ring buffer size again to reduce wrap copies.
- Dataset: 1000G_phase3_common_norel (extract 200k SNPs).
- Command: `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_small_kb1000 --ld-wind-kb 1000 --chunk-size 50 --extract /home/sharif/Code/ldsc/perf/small_l2_trace/extract.snps`
- Trace log: `/home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_trace_kb1000_ring16.log`
- Timing (single run): `compute_ldscore=10.568s` with
  `bed_read=0.601s`, `norm=1.099s`, `bb_dot=1.976s`, `ab_dot=4.152s`, `r2u_ab=0.433s`,
  `ab_wrap_copy=1.465s`, `ab_contig_matmul=2.781s`, `ab_wrap_matmul=1.372s`.

## 2026-03-04 (ring-buffer slack: max_window + 32*chunk, single run)
- Change: increased ring buffer size again to reduce wrap copies.
- Dataset: 1000G_phase3_common_norel (extract 200k SNPs).
- Command: `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_small_kb1000 --ld-wind-kb 1000 --chunk-size 50 --extract /home/sharif/Code/ldsc/perf/small_l2_trace/extract.snps`
- Trace log: `/home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_trace_kb1000_ring32.log`
- Timing (single run): `compute_ldscore=11.107s` with
  `bed_read=0.608s`, `norm=1.113s`, `bb_dot=2.326s`, `ab_dot=4.718s`, `r2u_ab=0.443s`,
  `ab_wrap_copy=1.048s`, `ab_contig_matmul=3.646s`, `ab_wrap_matmul=1.072s`.

## 2026-03-04 (remove r2u_ab zeroing, single run)
- Change: removed explicit zeroing of r2u_ab before overwrite.
- Dataset: 1000G_phase3_common_norel (extract 200k SNPs), ring +8c.
- Command: `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_small_kb1000 --ld-wind-kb 1000 --chunk-size 50 --extract /home/sharif/Code/ldsc/perf/small_l2_trace/extract.snps`
- Trace log: `/home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_trace_kb1000_nozero.log`
- Timing (single run): `compute_ldscore=10.594s` with
  `bed_read=0.603s`, `norm=1.097s`, `bb_dot=1.835s`, `ab_dot=3.945s`, `r2u_ab=0.436s`,
  `ab_wrap_copy=1.835s`, `ab_contig_matmul=2.335s`, `ab_wrap_matmul=1.610s`.

## 2026-03-04 (column-major loop order, single run)
- Change: loop j outer, wi inner in scalar accumulation to match column-major access.
- Dataset: 1000G_phase3_common_norel (extract 200k SNPs), ring +8c.
- Command: `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_small_kb1000 --ld-wind-kb 1000 --chunk-size 50 --extract /home/sharif/Code/ldsc/perf/small_l2_trace/extract.snps`
- Trace log: `/home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_trace_kb1000_colmajor.log`
- Timing (single run): `compute_ldscore=11.289s` with
  `bed_read=0.628s`, `norm=1.119s`, `bb_dot=2.044s`, `ab_dot=4.499s`, `r2u_ab=0.448s`,
  `ab_wrap_copy=1.707s`, `ab_contig_matmul=2.664s`, `ab_wrap_matmul=1.835s`.

## 2026-03-04 (chunk size sweep, 200k extract)
- Dataset: 1000G_phase3_common_norel (extract 200k SNPs), ring +8c, `--ld-wind-kb 1000`.
- Command template: `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_small_kb1000_c{C} --ld-wind-kb 1000 --chunk-size {C} --extract /home/sharif/Code/ldsc/perf/small_l2_trace/extract.snps`
- Trace logs: `/home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_trace_kb1000_chunk{C}.log`
- Results (single runs):
  - `chunk=25`: `compute=15.952s` `bb_dot=2.811s` `ab_dot=5.831s` `r2u_ab=0.444s` `wrap_copy=4.221s`
  - `chunk=50`: `compute=10.551s` `bb_dot=1.830s` `ab_dot=3.882s` `r2u_ab=0.438s` `wrap_copy=1.846s`
  - `chunk=100`: `compute=8.854s` `bb_dot=1.396s` `ab_dot=3.657s` `r2u_ab=0.446s` `wrap_copy=0.791s`
  - `chunk=200`: `compute=7.441s` `bb_dot=1.153s` `ab_dot=2.896s` `r2u_ab=0.489s` `wrap_copy=0.316s`
  - `chunk=400`: `compute=8.639s` `bb_dot=1.866s` `ab_dot=3.411s` `r2u_ab=0.586s` `wrap_copy=0.105s`
  - Best on this machine: `chunk=200`.

## 2026-03-04 (chunk size sweep, full 1000G)
- Dataset: 1000G_phase3_common_norel (full dataset), ring +8c, `--ld-wind-kb 1000`.
- Command template: `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_full_kb1000_c{C} --ld-wind-kb 1000 --chunk-size {C}`
- Trace logs: `/home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_trace_kb1000_chunk{C}.log`
- Results (single runs):
  - `chunk=50`:  `compute=88.816s` `bb_dot=17.486s` `ab_dot=37.617s` `r2u_ab=3.797s` `wrap_copy=8.429s`
  - `chunk=100`: `compute=70.450s` `bb_dot=10.618s` `ab_dot=30.346s` `r2u_ab=3.971s` `wrap_copy=3.863s`
  - `chunk=200`: `compute=66.656s` `bb_dot=9.993s` `ab_dot=27.191s` `r2u_ab=4.428s` `wrap_copy=1.753s`
  - `chunk=400`: `compute=76.829s` `bb_dot=14.644s` `ab_dot=30.248s` `r2u_ab=5.399s` `wrap_copy=0.671s`
  - Best on this machine: `chunk=200`.

## 2026-03-04 (ring alignment strategy, chunk=200)
- Change: ring_size aligned to chunk size (`ring_size % chunk == 0`).
- Dataset: 1000G_phase3_common_norel, `--ld-wind-kb 1000`, `--chunk-size 200`.
- 200k extract:
  - Trace log: `/home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_trace_kb1000_chunk200_align200k.log`
  - `compute_ldscore=7.586s`
  - `ab_wrap_copy=0.324s`, `ab_contig_matmul=2.669s`, `ab_wrap_matmul=0.773s`
- Full dataset:
  - Trace log: `/home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_trace_kb1000_chunk200_align_full.log`
  - `compute_ldscore=60.678s`
  - `ab_wrap_copy=1.711s`, `ab_contig_matmul=22.961s`, `ab_wrap_matmul=3.998s`

## 2026-03-04 (remove pq_window build, 200k extract)
- Change: compute `pq_k` on the fly instead of building `pq_window`.
- Dataset: 1000G_phase3_common_norel (extract 200k SNPs), ring aligned, `--chunk-size 200`.
- Command: `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_small_kb1000_c200_align --ld-wind-kb 1000 --chunk-size 200 --extract /home/sharif/Code/ldsc/perf/small_l2_trace/extract.snps`
- Trace log: `/home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_trace_kb1000_no_pq_window.log`
- Timing (single run): `compute_ldscore=6.884s` with
  `bed_read=0.631s`, `norm=1.102s`, `bb_dot=1.093s`, `ab_dot=2.916s`, `r2u_ab=0.490s`,
  `ab_wrap_copy=0.317s`, `ab_contig_matmul=2.271s`, `ab_wrap_matmul=0.645s`.

## 2026-03-04 (wrap column copy, 200k extract)
- Change: copy wrap window via column slices (`mat_copy_from`) instead of per-element loop.
- Dataset: 1000G_phase3_common_norel (extract 200k SNPs), ring aligned, `--chunk-size 200`.
- Command: `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_small_kb1000_c200_align --ld-wind-kb 1000 --chunk-size 200 --extract /home/sharif/Code/ldsc/perf/small_l2_trace/extract.snps`
- Trace log: `/home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_trace_kb1000_colcopy.log`
- Timing (single run): `compute_ldscore=7.654s` with
  `bed_read=0.630s`, `norm=1.141s`, `bb_dot=1.291s`, `ab_dot=3.147s`, `r2u_ab=0.491s`,
  `ab_wrap_copy=0.570s`, `ab_contig_matmul=2.391s`, `ab_wrap_matmul=0.756s`.

## 2026-03-04 (full l2 workflow checkpoint, best config)
- Config: ring aligned, `--chunk-size 200`, `--ld-wind-kb 1000`, no `pq_window` build.
- Trace log: `/home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_trace_kb1000_c200_align_best.log`
- Workflow timings:
  - `parse_bim=0.336s`, `count_fam=0.000093s`, `read_annot=0.000001s`
  - `maf_prefilter=7.247s`
  - `compute_ldscore=59.853s` with `bed_read=5.377s`, `norm=9.308s`,
    `bb_dot=9.533s`, `ab_dot=26.749s`, `r2u_bb=0.600s`, `r2u_ab=4.441s`,
    `ab_wrap_copy=1.714s`, `ab_contig_matmul=22.777s`, `ab_wrap_matmul=3.971s`
  - `write_outputs=2.480s`

## 2026-03-04 (ld score computation structure summary)
- Pipeline shape:
  - Parse BIM/FAM → optional SNP/individual filters → optional annotation load
  - Optional MAF prefilter pass over BED (computes per-SNP MAF + het/missing)
  - Main LD-score pass:
    1. Read BED chunk into `b_mat`, normalize columns
    2. Compute `BᵀB` (bb_dot) for within-chunk LD
    3. Maintain ring buffer of prior chunks to compute `AᵀB` (ab_dot)
    4. Apply unbiased r² transform and accumulate L2 into per-SNP output
  - Write outputs (`.ldscore.gz`, `.M`, `.M_5_50`, per-chr splits)
- Hot section: `compute_ldscore` dominated by `ab_dot` then `bb_dot`, with smaller `bed_read` + `norm` + `r2u` overheads.

## 2026-03-04 (noise check, 200k extract)
- Dataset: 1000G_phase3_common_norel (extract 200k SNPs), ring aligned, `--chunk-size 200`.
- Trace logs: `/home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_trace_kb1000_noise.run{0,1,2}.log`
- `compute_ldscore`: mean `6.9855s`, stdev `0.0643s`, min `6.9028s`, max `7.0597s`.

## 2026-03-04 (r2u constants, warmup + 3 runs)
- Change: precomputed r2u denominator constant (remove per-call branch/division).
- Dataset: 1000G_phase3_common_norel (extract 200k SNPs), ring aligned, `--chunk-size 200`.
- Trace logs: `/home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_trace_kb1000_r2uconst.run{0,1,2,3}.log` (run0 warmup).
- `compute_ldscore`: mean `7.0291s`, stdev `0.0432s`, min `6.9684s`, max `7.0651s`.
- Mean timings: `bed_read=0.639s`, `norm=1.105s`, `bb_dot=1.161s`, `ab_dot=2.954s`,
  `r2u_ab=0.482s`, `ab_wrap_copy=0.352s`, `ab_contig_matmul=2.270s`, `ab_wrap_matmul=0.684s`.

## 2026-03-04 (trace-branch hoist, warmup + 3 runs)
- Change: avoid `trace_idx` checks inside inner loop when tracing is disabled.
- Dataset: 1000G_phase3_common_norel (extract 200k SNPs), ring aligned, `--chunk-size 200`.
- Trace logs: `/home/sharif/Code/ldsc/perf/small_l2_trace/rust_l2_trace_kb1000_notrace.run{0,1,2,3}.log` (run0 warmup).
- `compute_ldscore`: mean `7.2729s`, stdev `0.1472s`, min `7.0697s`, max `7.4136s`.

## 2026-03-04 (maf prefilter fast path, full 1000G)
- Change: compute MAF/het stats directly from BED bytes (no intermediate matrix); contiguous BED block reads when possible.
- Dataset: 1000G_phase3_common_norel (full dataset), ring aligned, `--chunk-size 200`, `--ld-wind-kb 1000`.
- Command: `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_full_kb1000_chunk200_align_mafpre --ld-wind-kb 1000 --chunk-size 200`
- Trace log: `/home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_trace_kb1000_chunk200_align_mafpre.log`
- Timing (from trace): `maf_prefilter=0.630s`, `compute_ldscore=58.737s` with
  `bed_read=5.266s`, `norm=9.519s`, `bb_dot=9.625s`, `ab_dot=25.689s`, `r2u_ab=4.368s`,
  `ab_wrap_copy=1.755s`, `ab_contig_matmul=21.899s`, `ab_wrap_matmul=3.790s`.
- Wall time: `real 62.550s` (from `/home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_trace_kb1000_chunk200_align_mafpre.time`).

## 2026-03-04 (parity + perf checkpoint, 200k extract)
- Dataset: 1000G_phase3_common_norel (extract 200k SNPs), `--ld-wind-kb 1000`, `--chunk-size 200`.
- Python command: `uv run --project /home/sharif/Code/ldsc/ldsc_py --with bitarray==2 python /home/sharif/Code/ldsc/ldsc_py/ldsc.py --l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/parity_l2_200k/py_l2 --ld-wind-kb 1000 --chunk-size 200 --extract /home/sharif/Code/ldsc/perf/parity_l2_200k/extract.snps --yes-really --log-level warn`
- Rust command: `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/parity_l2_200k/rust_l2 --ld-wind-kb 1000 --chunk-size 200 --extract /home/sharif/Code/ldsc/perf/parity_l2_200k/extract.snps --yes-really`
- Parity: `max_abs_diff=0`, OK.
- Timing: Python `real 50.067s`, Rust `real 7.913s` → ~6.33x faster.

## 2026-03-04 (ld-wind-kb sweep, 200k extract, wiki-based)
- Rationale: LD Score Estimation Tutorial recommends `--ld-wind-cm 1`, but our BIM CM column is zeroed; we sweep `--ld-wind-kb` values around ~1 Mb and smaller windows consistent with wiki guidance that LD decays within ~100kb in high-recombination regions.
- Dataset: 1000G_phase3_common_norel (extract 200k SNPs), `--chunk-size 200`.
- Extract: `/home/sharif/Code/ldsc/perf/parity_l2_200k/extract.snps`
- Python command template: `uv run --project /home/sharif/Code/ldsc/ldsc_py --with bitarray==2 python /home/sharif/Code/ldsc/ldsc_py/ldsc.py --l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/parity_l2_200k_kb/kb{KB}/py_l2 --ld-wind-kb {KB} --chunk-size 200 --extract /home/sharif/Code/ldsc/perf/parity_l2_200k/extract.snps --yes-really --log-level warn`
- Rust command template: `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/parity_l2_200k_kb/kb{KB}/rust_l2 --ld-wind-kb {KB} --chunk-size 200 --extract /home/sharif/Code/ldsc/perf/parity_l2_200k/extract.snps --yes-really`
- Parity: `max_abs_diff=0` for all KB values.
- Results:
  - `kb=100`: Python `44.247s`, Rust `5.113s` → `8.65x`
  - `kb=500`: Python `48.429s`, Rust `6.213s` → `7.79x`
  - `kb=1000`: Python `53.719s`, Rust `8.614s` → `6.24x`
  - `kb=2000`: Python `61.834s`, Rust `12.884s` → `4.80x`
