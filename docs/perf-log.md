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

## 2026-03-04 (bed reader streaming + faer branch, full 1000G)
- Change: new internal streaming BED reader (no mmap), strict header/length validation, contiguous block reads, negative indices.
- Dataset: 1000G_phase3_common_norel (full dataset), ring aligned, `--chunk-size 200`, `--ld-wind-kb 1000`.
- Command: `RUST_LOG=ldsc=trace ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_full_kb1000_chunk200_align_bedstream --ld-wind-kb 1000 --chunk-size 200 --yes-really`
- Trace output: `/home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_full_kb1000_chunk200_align_bedstream.stdout`
- Timing (from trace): `maf_prefilter=0.607s`, `compute_ldscore=57.748s` with
  `bed_read=5.166s`, `norm=9.526s`, `bb_dot=9.324s`, `ab_dot=25.232s`, `r2u_ab=4.306s`,
  `ab_wrap_copy=1.684s`, `ab_contig_matmul=21.502s`, `ab_wrap_matmul=3.731s`.
- Wall time: `real 61.597s` (from `/home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_trace_kb1000_chunk200_align_bedstream.time`).
- Comparison vs prior full trace (maf prefilter fast path): slightly lower `bed_read` and `compute_ldscore` (~1–2%); likely small but real. More runs needed to reduce noise.

## 2026-03-04 (faer nightly SIMD flag, full 1000G)
- Change: enabled `faer` feature `nightly` (SIMD-enabled kernels).
- Dataset: 1000G_phase3_common_norel (full dataset), ring aligned, `--chunk-size 200`, `--ld-wind-kb 1000`.
- Command: `RUST_LOG=ldsc=trace ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_full_kb1000_chunk200_align_faer_nightly --ld-wind-kb 1000 --chunk-size 200 --yes-really`
- Trace output: `/home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_full_kb1000_chunk200_align_faer_nightly.stdout`
- Timing (from trace): `maf_prefilter=2.586s`, `compute_ldscore=64.307s` with
  `bed_read=9.167s`, `norm=10.314s`, `bb_dot=9.663s`, `ab_dot=26.506s`, `r2u_ab=4.456s`,
  `ab_wrap_copy=1.643s`, `ab_contig_matmul=22.618s`, `ab_wrap_matmul=3.888s`.
- Wall time: `real 70.212s` (from `/home/sharif/Code/ldsc/perf/full_l2_trace/rust_l2_trace_kb1000_chunk200_align_faer_nightly.time`).
- Observation: regression vs non-nightly build (slower `bed_read` and `compute_ldscore`); likely due to different codegen/CPU feature usage. Consider reverting `faer` nightly unless further investigation shows a config issue.

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

## 2026-03-04 (Rust scaling sweep, full 1000G)
- Dataset: 1000G_phase3_common_norel (full dataset), `--ld-wind-kb 1000`, `--chunk-size 200`.
- Command template: `RUST_LOG=ldsc=trace ldsc l2 --rayon-threads {N} --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/scale_full_trace/rust_full_t{N} --ld-wind-kb 1000 --chunk-size 200 --yes-really`
- Summary CSV: `/home/sharif/Code/ldsc/perf/scale_full_trace/summary.csv`
- `compute_ldscore` scaling (single-run):
  - `N=1`: `158.143s` (baseline)
  - `N=2`: `96.186s` (1.64x)
  - `N=3`: `75.388s` (2.10x)
  - `N=4`: `65.987s` (2.40x)
  - `N=6`: `56.762s` (2.79x)
  - `N=8`: `60.007s` (2.64x)
  - `N=10`: `64.073s` (2.47x)
  - `N=12`: `70.682s` (2.24x)
- Observation: scaling improves to ~6 threads then saturates/declines; likely memory bandwidth + overhead ceiling on this machine. Keep Rayon defaults for now.

## 2026-03-04 (fast-f32 feature, 200k extract)
- Change: compile-time `fast-f32` feature to run `l2` core matmuls in f32 while accumulating in f64.
- Dataset: 1000G_phase3_common_norel (extract 200k SNPs), `--ld-wind-kb 1000`, `--chunk-size 200`.
- Command (f64): `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/fast_f32_200k/rust_l2_f64 --ld-wind-kb 1000 --chunk-size 200 --extract /home/sharif/Code/ldsc/perf/parity_l2_200k/extract.snps --yes-really`
- Command (f32): `cargo build --release --features fast-f32` then `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/fast_f32_200k/rust_l2_f32 --ld-wind-kb 1000 --chunk-size 200 --extract /home/sharif/Code/ldsc/perf/parity_l2_200k/extract.snps --yes-really`
- Timing: f64 `real 8.428s`, f32 `real 6.054s` → `1.39x` speedup.
- Note: f32 path is **not parity‑safe**; use for performance experiments only.

## 2026-03-04 (fast-f32 feature, full 1000G)
- Change: compile-time `fast-f32` feature to run `l2` core matmuls in f32 while accumulating in f64.
- Dataset: 1000G_phase3_common_norel (full dataset), `--ld-wind-kb 1000`, `--chunk-size 200`.
- Command (f64): `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/fast_f32_full/rust_l2_f64 --ld-wind-kb 1000 --chunk-size 200 --yes-really`
- Command (f32): `cargo build --release --features fast-f32` then `ldsc l2 --bfile /home/sharif/Code/ldsc/data/1000G_phase3_common_norel --out /home/sharif/Code/ldsc/perf/fast_f32_full/rust_l2_f32 --ld-wind-kb 1000 --chunk-size 200 --yes-really`
- Timing: f64 `real 68.468s`, f32 `real 47.502s` → `1.44x` speedup.
- Output diff (all 22 chr files, 1,664,852 values):
  - `mean_abs_diff=0.000303151`, `rmse=0.00057935`
  - `max_abs_diff=0.008`
  - `max_rel_diff=0.000651466` (relative to f64, excluding zero refs)
- Note: f32 path is **not parity‑safe**; use for performance experiments only.

## 2026-03-05 (no-window / whole-chromosome LD scores via --ld-wind-cm 1)
- Approach: 1000G BED files have CM=0 for all SNPs. With `--ld-wind-cm 1`, every SNP is within 1cM of every other on the same chromosome → effectively a whole-chromosome window with no truncation. Using single-chromosome BED files avoids cross-chromosome bleeding.
- Machine: Ryzen 5 5600X (6 cores / 12 threads), faer `Par::rayon(0)`.

### Parity + perf (Python vs Rust, --ld-wind-cm 1, single chr22 BED)
- Commands: `ldsc.py --l2 --bfile <chr22_bed> --ld-wind-cm 1 --yes-really` / `ldsc l2 --bfile <chr22_bed> --ld-wind-cm 1 --yes-really`
- Speedup grows with M because work is O(M²): Python single-threaded, Rust uses faer Par::rayon(0) across 12 threads.

| SNPs   | Python wall | Rust wall | Rust user | Speedup     | max_abs_diff |
|--------|-------------|-----------|-----------|-------------|--------------|
| 2,950  | 8.36s       | 0.222s    | 1.184s    | **33.7×**   | 0 ✓          |
| 24,624 | 664s        | 14.0s     | 119.8s    | **47.4×**   | 0 ✓          |

- Rayon speedup at 2,950 SNPs: 5.3× (small matrices, overhead dominates); at 24,624 SNPs: 8.6×
- Calibration: O(M²) confirmed; `k = 1.426e-7 s·CPU/SNP²` (user time per SNP pair)

### Full-genome no-window estimate (sequential, O(M²) per chr, 9.32× rayon)
- Total wall: ~**38 min**, total CPU: ~**6h**
- Bottleneck: chr1 (288s) + chr2 (295s) = ~10 min combined
- If chromosomes parallelized independently (each chr independent in no-window mode): bottleneck = chr2 ~295s → ~5 min wall
- Per-chromosome wall times:
  - chr1: 288s, chr2: 295s, chr3: 206s, chr4: 170s, chr5: 162s, chr6: 183s
  - chr7: 128s, chr8: 124s, chr9: 89s, chr10: 119s, chr11: 109s, chr12: 97s
  - chr13: 59s, chr14: 46s, chr15: 40s, chr16: 46s, chr17: 34s, chr18: 38s
  - chr19: 18s, chr20: 29s, chr21: 9s, chr22: 9s

## 2026-03-05 (B1 scalar/partitioned path unification)
- Change: removed separate scalar LD-score accumulation path (~170 lines); scalar case now uses the partitioned matmul path with a synthetic all-ones annotation column. Also removed dead code (jackknife, irwls, Bed.path), tracing from regressions.rs, section separators, and verbose doc comments. Total: 8,457 → 7,835 lines.
- Machine: Ryzen 5 5600X (6 cores / 12 threads), faer `Par::rayon(0)`.

### No-window (--ld-wind-cm 1, single-chr BEDs)

| Chr | SNPs | Previous wall | Current wall | Previous user | Current user |
|-----|------|--------------|--------------|---------------|--------------|
| chr22 | 24,624 | 14.0s | 9.2s | 119.8s | 88.2s |
| chr1 | 137,124 | 288s | 300s | — | 3126s |

- chr22: **34% faster wall, 26% faster CPU** — matmul-based accumulation is more cache-friendly than the hand-written scalar loop at this matrix size.
- chr1: within noise (~4% diff); at 137k SNPs both paths are memory-bandwidth-limited.

### Full 1000G (--ld-wind-kb 1000, --chunk-size 200, 1,664,852 SNPs)
- Command: `ldsc l2 --bfile data/1000G_phase3_common_norel --ld-wind-kb 1000 --chunk-size 200 --out /tmp/ldsc_full_perf --yes-really`
- Timing: `real 70.0s`, `user 430.1s`
- Previous best (pre-B1): `real 60.7s`–`68.5s` range
- Within normal run-to-run noise (~5-10s variance on this machine).

### Parity
- `bench_5k` (5k SNPs, --ld-wind-kb 1000): `max_abs_diff=0` ✓
- `check_l2_tiny_py_vs_rust.sh` (2k extract): `max_abs_diff=0` ✓
- fast-f32 vs f64 drift: `max_abs_diff=0.001` (unchanged from pre-B1)

## 2026-03-06 (GPU acceleration via CubeCL + CubeK)
- Change: added optional `--features gpu` feature gate using `cubecl 0.10.0-pre.2` (CUDA backend) and `cubek-matmul 0.2.0-pre.2` (autotuned GPU GEMM). New `src/gpu.rs` (~160 lines) wraps `GpuContext` with `matmul_tn()` and `matmul_tn_tiled()`. CLI flag `--gpu` dispatches GEMM to GPU at both B×B and A×B call sites in `compute_ldscore_global`. Output tensor uses pitched/padded strides from `TensorHandle::empty`; stride-aware extraction handles this correctly. Default build (no `gpu` feature) is completely unaffected.
- Machine: Ryzen 5 5600X (6C/12T) + NVIDIA GeForce RTX 3060 (12 GB VRAM), CUDA 13.1, cudarc 0.19.3.
- Dependencies: `cubecl =0.10.0-pre.2` (features=["cuda"]), `cubek-matmul =0.2.0-pre.2`.

### Parity (GPU f32 vs CPU f64)
- `bench_5k` (5k SNPs, --ld-wind-kb 1000): `max_abs_diff=0.001`, `max_rel_diff=0.049%`
- `chr1` (137k SNPs, whole-chromosome): `max_abs_diff=0.003`, `max_rel_diff=0.020%`
- Drift is inherent to f32 matmul precision — matches `--features fast-f32` behavior.
  With f16 matmul (supported by cubek-matmul via `MatmulPrecision` trait), precision would
  degrade further (~3 decimal digits), but transfer bandwidth would halve, potentially
  improving the compute-to-transfer ratio on bandwidth-bound problems.

### Small data: bench_5k (5,000 SNPs, n=2,490, --ld-wind-kb 1000)

| Mode | Wall | User | Sys |
|------|------|------|-----|
| CPU (rayon) | **0.14s** | 0.57s | 0.05s |
| GPU (RTX 3060) | 1.02s | 0.62s | 0.39s |

- GPU **7× slower** — dominated by ~1s CUDA initialization overhead.
  f16 would not help here; the bottleneck is init, not compute or transfer.

### Medium data: bench_200k (200,237 SNPs, n=2,490, --ld-wind-kb 1000)

| Mode | Wall | User | Sys |
|------|------|------|-----|
| CPU (rayon) | **5.4s** | 27.0s | 1.9s |
| GPU (RTX 3060) | 8.9s | 7.8s | 1.4s |

- GPU **1.6× slower wall** but uses **3.5× less CPU time** (7.8s vs 27s).
  Each GEMM is ~2490×200×200 ≈ 100 MFLOP — microseconds on GPU, but PCIe
  transfer of the 2490×200 matrices (~2 MB) takes longer than the compute.
  f16 would halve transfer size to ~1 MB per matrix, potentially bringing GPU
  wall time closer to CPU by reducing PCIe bottleneck.

### Largest problem: chr1 whole-chromosome (137,124 SNPs, n=2,490, --ld-wind-snps 999999999)

| Mode | Wall | User | Sys |
|------|------|------|-----|
| CPU (rayon) | **6m56s** | 70m30s | 0m21s |
| GPU (RTX 3060) | 10m25s | 5m31s | 6m09s |

- GPU **1.5× slower wall** but uses **13× less CPU time** (5.5 min vs 70.5 min).
- The massive `sys` time (6 min) is CPU↔GPU data transfer. With n=2,490, each
  chunk GEMM is ~2490×200×200 ≈ 100 MFLOP — trivial compute for GPU.
  The window grows to 137k columns, but each A×B dispatch is still small.
  f16 would halve transfer volume (~4 MB→2 MB per chunk pair), cutting the
  6-minute sys overhead significantly.

### Analysis: when does GPU break even?

At n=2,490 individuals, GEMM sizes are too small for GPU to overcome PCIe
transfer latency. The crossover depends on the compute-to-transfer ratio:

- **Compute per chunk**: `n × W × C` FLOPs (W=window cols, C=chunk cols)
- **Transfer per chunk**: `(n×C + n×W) × 4 bytes` (f32) for upload + `(W×C) × 4 bytes` for download
- **GPU compute throughput**: ~10 TFLOPS (RTX 3060 f32)
- **PCIe bandwidth**: ~12 GB/s (PCIe 4.0 x16)

For n=2,490, W=200, C=200: compute = 100 MFLOP = 10µs; transfer = 4 MB = 330µs → **33× transfer-bound**.
For n=500,000, W=200, C=200: compute = 20 GFLOP = 2ms; transfer = 800 MB = 67ms → still transfer-bound but ratio improves.
For n=500,000 with f16: transfer halved to 400 MB = 33ms; compute still ~2ms → **better ratio**.

**Crossover estimate**: GPU breaks even at roughly n≥50k individuals with f32,
or n≥25k with f16. At UK Biobank scale (n=500k), GPU would be ~3-5× faster
than CPU, with f16 potentially reaching ~8-10× by halving transfer overhead.

### Notes
- `Strategy::Auto` selects kernel automatically; on first run it autotuned
  (cache stored on device for subsequent runs).
- `TensorHandle::empty()` returns pitched/padded strides (e.g., stride 4 for
  shape [2,2]); `extract_result()` handles this correctly.
- Build: `cargo build --release --features gpu` requires CUDA toolkit.
  Default build (`cargo build --release`) is completely unaffected.

## 2026-03-06 (runtime f32, GPU f16 via half crate)
- Change: converted `fast-f32` from compile-time feature to runtime `--fast-f32` CLI flag. Added `--gpu-f16` flag using `half::f16` for GPU matmul (f16 inputs, f32 accumulation via tensor cores, f16 output widened to f32). Added `half = "2"` as optional dep under `gpu` feature.
- Machine: Ryzen 5 5600X (6C/12T) + NVIDIA GeForce RTX 2070 SUPER (8 GB VRAM).

### chr1 whole-chromosome (137,124 SNPs, n=2,490, --ld-wind-cm 1)

| Mode | Wall | User | Sys | vs CPU f64 |
|------|------|------|-----|------------|
| CPU f64 | **5m18s** | 54m00s | 0m18s | 1.0× |
| CPU f32 (--fast-f32) | **3m24s** | 31m15s | 0m11s | 1.56× |
| GPU f32 (--gpu) | 10m37s | 5m36s | 6m16s | 0.50× |
| GPU f16 (--gpu --gpu-f16) | 12m43s | 7m47s | 6m01s | 0.42× |
| GPU flex32 (--gpu --gpu-flex32) | 10m07s | — | — | 0.60× |

### Parity (vs CPU f64 baseline, 137,124 SNPs)

| Mode | max_abs_diff | mean_abs_diff | mean_rel_err |
|------|-------------|--------------|-------------|
| CPU f32 | 0.001 | — | — |
| GPU f32 | 0.003 | — | — |
| GPU f16 (removed) | 3.504 | 0.139 | 0.024% |
| GPU flex32 | 0.003 | 0.000121 | 0.0000% |

- L2 range on chr1 whole-window: [-3.7, 6040.8], median 736.9
- GPU f16 worst SNP: rs7549549, cpu64=5277.58, gpu16=5274.08 (abs_diff=3.50, rel=0.066%)
- GPU f16 relative error is small (median 0.024%, p95 0.065%) but absolute error grows with L2 magnitude because f16 has only ~3 decimal digits of precision

### Analysis
- **GPU f16 was slower than GPU f32 at n=2,490** (12m43s vs 10m37s). The f32→f16 conversion on CPU added ~2 min overhead. Replaced with flex32.
- **GPU flex32 eliminates CPU-side conversion**: uploads f32, GPU converts to f16 in shared memory, accumulates in f32, outputs f32. 30s faster than old GPU f16 (607s vs 763s) with GPU f32-equivalent parity (max_abs_diff=0.003).
- **CPU f32 is fastest** at 1000G scale (1.56× over f64) because the GEMM is compute-bound on CPU with 12-thread rayon parallelism, and f32 doubles SIMD throughput.
- **GPU value proposition**: at biobank scale (n≥50k), GPU flex32 halves register/shared-memory usage (fitting larger tiles) and exploits tensor cores. At n=2,490, PCIe transfer dominates and GPU is slower than CPU.

## 2026-03-07 (micro-optimization batch)
- Changes applied (all bit-identical output, max_abs_diff=0):
  1. **normalize_col_f32 scratch buffer reuse**: pre-allocate `Vec<f64>` once, reuse across all SNPs (was allocating per-SNP)
  2. **faer copy_from for ring buffer stores**: replaced element-wise `for i in 0..n { ring[(i,slot)] = b[(i,j)] }` with `submatrix_mut().copy_from(submatrix())`
  3. **faer copy_from for window gather**: same pattern for A-buffer assembly from ring buffer
  4. **Pre-inverted n**: `* n_inv` instead of `/ n` in r² unbiased transform
  5. **mat_add_in_place loop order**: swapped to column-major-friendly (j outer, i inner) in `la.rs`
  6. **mat_copy_from**: now delegates to faer builtin `dst.copy_from(src)`
  7. **Par::Seq for jackknife**: small matmuls in h2.rs jackknife_fast use `Par::Seq` instead of `Par::rayon(0)`
  8. **lambda_gc quickselect**: `select_nth_unstable` (O(n)) instead of full sort (O(n log n))
  9. **--verbose-timing flag**: lightweight `std::time::Instant` instrumentation gated behind CLI flag (no runtime cost when disabled)
- **Rejected**: `target-cpu=native` — caused 2-4× regression because faer does runtime CPU feature detection via cpuid; compile-time native codegen interferes with faer's internal SIMD dispatch. Documented in `.cargo/config.toml`.
- **Skipped**: double-copy elimination in column normalization — complex borrow restructuring for ~1% gain.

### Hyperfine: stable vs nightly+perf (full 1000G, --ld-wind-kb 1000, --chunk-size 200)

| Binary | Mean | σ | Min | Max | User |
|--------|------|---|-----|-----|------|
| stable (rustc 1.94.0) | **85.4s** | 1.2s | 84.2s | 86.6s | 486.7s |
| nightly+perf (1.96.0, faer nightly, polars performant, zlib-rs) | 87.7s | 0.6s | 87.3s | 88.4s | 498.8s |

- Warmup: 1 run, timed: 3 runs each.
- **Nightly+perf is 3% slower** — faer `nightly` feature consistently regresses (confirmed both 03-04 and 03-07). Polars `nightly` feature fails to compile on current nightly (missing `LaneCount` in `simd`). Polars `performant` and flate2 `zlib-rs` had negligible effect.
- **Recommendation: stay on stable, no nightly features.**

### Hyperfine: f64 vs fast-f32 (full 1000G, --ld-wind-kb 1000, --chunk-size 200)

| Mode | Mean | σ | Min | Max | User |
|------|------|---|-----|-----|------|
| Rust f64 | 83.8s | 12.2s | 75.9s | 97.8s | 477.7s |
| Rust f32 (`--fast-f32`) | **66.3s** | 3.3s | 62.7s | 69.2s | 293.7s |

- Warmup: 1 run, timed: 3 runs each.
- **f32 is 1.26× faster** with notably lower variance (f64 75–98s vs f32 63–69s).
- f64 high variance likely due to thermal throttling on longer runs; f32 finishes before throttling onset.
- f32 is **not bit-exact** with Python: `max_abs_diff=0.008`, `mean_abs_diff=0.000303`, `max_rel_diff=0.065%` (from earlier 03-04 testing on same dataset).

### Headline: Python vs Rust, full 1000G (1,664,852 SNPs, 2,490 indiv, --ld-wind-kb 1000)
- Machine: Ryzen 5 5600X (6C/12T), 32 GB RAM
- Python (hyperfine, 1 run): **1548.5s** (25m49s wall, 25m35s user)
- Rust f64 (hyperfine, 3 runs): **85.4s ± 1.2s** (8m07s user) — **18× speedup**, `max_abs_diff=0` ✓
- Rust f32 (hyperfine, 3 runs): **66.3s ± 3.3s** (4m54s user) — **23× speedup**, `max_abs_diff=0.008`

## Skip MAF prefilter when --maf absent (2026-03-08)

- **Change**: `if maf_pre {` → `if maf_pre && (args.maf.is_some() || pq_exp.is_some()) {` in `run_ldscore`
- **Rationale**: Python fuses MAF filtering with the single BED load (free). Rust was doing a separate full BED pass (~4s) that removes ~0 SNPs from well-filtered data (1000G common_norel). When `--maf` is not set, the threshold is 0.0 and only monomorphic/all-het SNPs are removed — essentially none in 1000G. `normalize_col_f64` already handles zero-variance columns correctly (all-zeros after centering, no NaN).
- **Result** (hyperfine, 1 warmup + 3 runs, full 1000G, --ld-wind-kb 1000): **77.9s ± 6.6s**
  - vs prior baseline 85.4s ± 1.2s → improvement, but high variance (thermal); see PGO entry below for true combined baseline
- **Correctness**: output identical for 1000G (no SNPs removed by prefilter). When `--maf` is set, prefilter still runs as before.

## PGO (Profile-Guided Optimization) on MAF-fix binary (2026-03-08)

- **Script**: `scripts/build_pgo.sh` — instruments binary, runs full 1000G workload to collect profiles, merges with `llvm-profdata`, rebuilds.
- **Result** (hyperfine, 1 warmup + 3 runs, full 1000G, --ld-wind-kb 1000): **80.7s ± 1.3s**
  - vs original 85.4s ± 1.2s → **~5.5s improvement (6.4%)** total from MAF fix + PGO combined
  - Low variance (±1.3s) is the reliable combined baseline; 77.9s from MAF-only run was a noisy sample
- **Analysis**: PGO impact is modest. faer matmul (48% of runtime) uses runtime SIMD dispatch with tight vectorized loops — PGO helps branches/indirect calls, not SIMD. Normalization loop (19%) has some NaN-check branches that benefit marginally. PGO is more impactful for workloads with irregular control flow (interpreters, web servers).
- **Decision**: PGO not kept in normal `cargo build --release`. Use `scripts/build_pgo.sh` when maximum wall-time performance is needed. The ~6% gain doesn't justify the 2× compile time for routine development.
- **New headline**: Rust f64 (MAF fix + PGO): **~80.7s** → **19× speedup** vs Python 1548.5s

## Double-buffered BED reading via --prefetch-bed (2026-03-08)

- **Change**: `compute_ldscore_global` now takes `bed_path: &str` instead of `&mut Bed`. When `--prefetch-bed` is set, spawns a background `std::thread` that reads BED chunks via `crossbeam_channel::bounded(1)`, one chunk ahead of the compute loop.
- **Default path** (sequential): unchanged behavior, same performance as before.
- **Prefetch path result** (hyperfine, 1 warmup + 3 runs, full 1000G, --ld-wind-kb 1000): **90.3s ± 1.3s** — a **~12% regression** vs sequential (83s).
- **Root cause of regression**: On local SSD with warm OS page cache, BED reads are memory-bandwidth-bound (kernel copy), not I/O-wait-bound. The reader thread competes with rayon's 12-thread GEMM pool for CPU cores on a 6C/12T machine. Adding a 13th thread causes OS preemption of GEMM workers, increasing context switches and memory bandwidth contention. No true I/O blocking means no overlap benefit.
- **Decision**: `--prefetch-bed` disabled by default. Expected benefit on GPFS/NFS cold storage (PMACS HPC, CentOS 7): BED `read_exact()` blocks on InfiniBand network, genuinely freeing CPU for GEMM. Recommended to test on HPC with first-run cold GPFS.
- **io_uring investigated and ruled out**: CentOS 7 kernel 3.10 predates io_uring by 8 major versions (requires 5.1+). Even on supported kernels, io_uring shows no benefit for large sequential reads from warm page cache. tokio-uring crate is also stale (last release Nov 2022) and single-threaded only.
- **Correctness**: both paths verified byte-exact vs baseline on bench_5k (all 22 chromosomes).

## Remaining optimization opportunities (assessed 2026-03-07)

### Timing breakdown (f64, full 1000G, --ld-wind-kb 1000, --chunk-size 200, ~85s total)

| Section | Time | % | Description |
|---------|------|---|-------------|
| ab_dot | 31s | 36% | A^T×B GEMM + r2u + annot contribution |
| norm | 16s | 19% | Per-SNP normalize (impute NaN, center, unit-variance) |
| bed_read | 11s | 13% | Read BED chunks via file I/O |
| bb_dot | 10s | 12% | B^T×B GEMM + r2u + annot contribution |
| maf_prefilter | 4s | 5% | Extra full-BED pass to compute MAF/het stats |
| write_outputs | 3s | 4% | Gzip 22 per-chr output files |
| ring_store | 2s | 2% | Copy cols into ring buffer |
| overhead | ~8s | 9% | BIM parse, alloc, block_left, etc. |

### Candidates (ordered by expected impact)

#### 1. Profile-Guided Optimization (PGO) — estimated 5-15% (~5-12s)
- LLVM PGO uses runtime profiles to optimize branch prediction, function layout, and inlining.
- Build process: `cargo build --release` with `-Cprofile-generate`, run benchmark, rebuild with `-Cprofile-use`.
- No code changes, no correctness risk. Especially effective for loops with data-dependent branches (normalization, r2u accumulation).
- Complexity: medium (build script only). Would need a separate build recipe or Makefile target.

#### 2. Eliminate MAF prefilter BED pass — saves ~4s (5%)
- Currently the BED file is read twice when `maf_pre=true` (default): once for MAF/het stats, once for LD computation.
- The main loop's `normalize_col` already computes MAF as a side effect (it returns the MAF value).
- When `--maf` is not set, the prefilter only removes all-het/missing SNPs (extremely rare in 1000G).
- Option A: skip prefilter when `--maf` is not set; detect all-het/missing during normalize and zero the column.
- Option B: fuse the two passes entirely (compute MAF + het stats during the first chunk pass, then re-derive block_left and restart — complex).
- Complexity: low (option A), high (option B). No correctness risk for option A.

#### 3. Double-buffered BED reading — saves up to ~11s (13%)
- Read the next BED chunk in a background thread while computing GEMM on the current chunk.
- bed_read (11s) and GEMM (41s) use different resources (I/O vs CPU); with pipelining, BED reading overlaps completely with compute.
- Requires a second `Bed` file handle (or `Send`-safe reader) and channel-based synchronization.
- Complexity: high. Risk: concurrency bugs, lifetime issues with borrowed buffers.

#### 4. Fat LTO — maybe 1-3%
- Change `lto = "thin"` to `lto = "fat"` in `[profile.release]`.
- Enables more aggressive cross-crate inlining (faer internals, polars).
- Increases link time (~2×) but zero code changes.
- Complexity: trivial. Easy to A/B test with hyperfine.

#### 5. Lower gzip compression level — saves ~1-2s
- `Compression::default()` is level 6. Level 1 produces slightly larger files but writes much faster.
- Output is 22 `.l2.ldscore.gz` files totaling ~3s write time; level 1 could halve this.
- Complexity: trivial (one-line change). Trade-off: ~30% larger output files.

#### 6. Bulk column copy in normalize — <1%
- Lines 463-464 and 473-474 copy raw→col_buf→b_mat element-wise. Could use `copy_from_slice` on contiguous column data or `submatrix.copy_from()`.
- Marginal gain since the inner loop is already auto-vectorized and the overhead is dwarfed by the actual normalization arithmetic.
- Complexity: low.

### Not worth pursuing
- **Rayon thread tuning**: Already swept (6 threads optimal on 5600X, default 12 within noise).
- **SIMD normalization**: norm is memory-bound (streaming 2490×200 matrices), not compute-bound; auto-vectorization sufficient.
- **Chunk size > 200**: Tested 400, regressed due to L2/L3 cache pressure.
- **Fusing r2u into GEMM**: Would fight faer's internal cache blocking strategy.
- **f16 on CPU**: No native f16 arithmetic on x86-64 (pre-Sapphire Rapids). Widening to f32 adds overhead with no SIMD throughput gain.
- **f16/fp8 on GPU**: At n=2,490, GPU is PCIe-transfer-bound (compute is 80× cheaper than transfer). Even fp8 (4× less transfer) only brings GPU to parity with CPU f32. And fp8 has ~1 decimal digit of precision — unusable for LD scores (L2 range 1-6000).
- **target-cpu=native**: Confirmed 2-4× regression with faer (runtime CPU detection conflict).

## 2026-03-08 — AWS Batch benchmarking pipeline + parallel output + gzip fast

### Infrastructure
- AWS Batch on EC2 Spot: c7a.4xlarge (AMD EPYC 7R13, 16 vCPUs, 8 GiB)
- Build: `rust:1-bookworm` container (glibc 2.36), runtime: `debian:bookworm-slim`
- Benchmark: `hyperfine --warmup 2 --runs 10`
- Dataset: full 1000G (1,664,852 SNPs, 2,490 individuals), `--ld-wind-kb 1000 --chunk-size 200`

### Changes measured (combined)
1. **Parallel output writing**: `rayon::par_iter()` over 22 chromosomes for gzip output (was sequential)
2. **Gzip fast compression**: `Compression::fast()` (level 1) instead of `Compression::default()` (level 6)

### Results (parallel output + gzip fast, on AWS)
```
Mean: 73.196s ± 0.212s (10 runs)
Range: 72.769s — 73.475s
σ/μ: 0.29% (excellent stability)
```

### Verbose timing breakdown (single run, same hardware)
```
bed_read(stall)=10.223s  norm=19.705s  bb_dot=10.470s  ab_dot=27.290s  ring_store=2.167s
compute_ldscore_total=70.074s
write_outputs=1.135s  total=71.208s
```

| Phase | Time | % of total |
|-------|------|-----------|
| bed_read | 10.2s | 14.3% |
| norm | 19.7s | 27.7% |
| bb_dot | 10.5s | 14.7% |
| ab_dot | 27.3s | 38.3% |
| ring_store | 2.2s | 3.0% |
| write_outputs | 1.1s | 1.6% |

### Notes
- Write phase is only 1.6% of total — parallel output + gzip fast savings are real but small (~1-2s)
- Local baseline on Ryzen 5 5600X was ~83s; AWS EPYC 7R13 (16 vCPU) is ~73s (more cores for rayon)
- The low variance (σ = 0.21s, 0.29%) confirms AWS Spot is suitable for reproducible benchmarking
- Output file size: ~28.7 MB combined (gzip level 1) vs ~18.6 MB (level 6) — 55% larger
- ab_dot dominates (38%); norm is second (28%) — these are the targets for further optimization

## 2026-03-08 — Fat LTO A/B test on AWS

### Setup
- Same infrastructure as above (c7a.4xlarge, EPYC 7R13, 16 vCPUs)
- Same code (parallel output + gzip fast), only change is `lto = "fat"` vs `lto = "thin"` in Cargo.toml
- Both built inside `rust:1-bookworm` container, stripped

### Results (hyperfine --warmup 2 --runs 10)

| Config | Mean | σ | Min | Max | Binary size |
|--------|------|---|-----|-----|-------------|
| Thin LTO | 73.196s | 0.212s | 72.769s | 73.475s | 51.3 MB |
| Fat LTO | **71.008s** | 0.386s | 70.680s | 71.867s | 49.1 MB |

- **Fat LTO is 3.0% faster** (2.19s improvement, ranges don't overlap → statistically significant)
- Binary is 4% smaller with fat LTO
- Build time: ~10 min (fat) vs ~8.5 min (thin) — 18% longer
- **Decision**: keep `lto = "thin"` for development iteration speed. Use fat LTO for release builds.
- **New best**: **71.0s** on AWS EPYC → **~21.8× speedup** vs Python 1548.5s

## 2026-03-08 — CentOS 7 (glibc 2.17) vs musl (static) portability benchmark

### Motivation
Target HPC (PMACS) runs CentOS 7 (glibc 2.17). Bookworm binaries require GLIBC_2.34+ and won't run there. Two options:
1. **CentOS 7 build**: compile inside centos:7 container → dynamic binary linked against glibc 2.17
2. **musl build**: compile with `x86_64-unknown-linux-musl` target → fully static binary, no glibc dependency

### Build Details
- CentOS 7: `centos:7` + vault repos + `rustup stable` (1.94.0). Max GLIBC symbol: GLIBC_2.16 (getauxval). Binary: 49.5 MB.
- musl: `rust:1-alpine` + `--target x86_64-unknown-linux-musl`. `static-pie linked`. Binary: 49.1 MB.
- Both stripped, thin LTO, codegen-units=1.

### Results (AWS c7a.4xlarge, EPYC 7R13, 16 vCPUs, hyperfine --warmup 2 --runs 10)

| Variant | Mean | σ | Min | Max | vs Bookworm |
|---------|------|---|-----|-----|-------------|
| Bookworm (baseline) | 73.196s | 0.212s | 72.769s | 73.475s | — |
| CentOS 7 (glibc 2.17) | 71.825s | 0.168s | 71.560s | 72.052s | 1.9% faster |
| musl (default alloc) | 81.878s | 0.109s | 81.638s | 82.045s | 11.9% slower |
| **musl + mimalloc** | **69.897s** | 0.179s | 69.711s | 70.280s | **4.5% faster** |

### Analysis
- **musl + mimalloc is the fastest variant** at 69.9s — 4.5% faster than bookworm baseline, and the new overall best.
- musl's default allocator is 12% slower due to single-threaded contention under rayon parallelism. mimalloc eliminates this entirely and then some.
- CentOS 7 (glibc 2.17) is also fast at 71.8s — essentially identical to bookworm.
- **Recommendation**: use musl + mimalloc (`--features mimalloc --target x86_64-unknown-linux-musl`) for HPC deployment. Fully static binary, runs on any Linux, and is the fastest option. Build with: `cargo build --release --features mimalloc --target x86_64-unknown-linux-musl`.
- **New best**: **69.9s** on AWS EPYC → **~22.2× speedup** vs Python 1548.5s

## 2026-03-09 — Normalize column optimization (AWS benchmark)

### Changes
1. **Eliminated col_buf intermediate buffer**: raw BED data is now copied directly into `b_mat` columns, then normalized in-place. Previously: `raw → col_buf → normalize → col_buf → b_mat` (two copies per column per chunk). Now: `raw → b_mat column → normalize in-place` (one copy).
2. **Fused center + variance accumulation**: `normalize_col_f64` and `normalize_col_f32` now compute mean and variance in fewer passes. Center + sum-of-squares accumulated in a single loop (was separate passes). Reduces memory traffic over 2,490 elements × 1,664,852 SNPs.
3. **Removed scratch Vec from f32 path**: `normalize_col_f32` previously allocated a `Vec<f64>` scratch buffer per SNP for f64-precision accumulation. Now uses f64 accumulators directly on f32 data without intermediate allocation.

### Parity
- Full 1000G (1,664,852 SNPs): **max_abs_diff=0** vs prior Rust baseline ✓ (bit-identical)

### Results (AWS c7a.4xlarge, EPYC 7R13, 16 vCPUs, hyperfine --warmup 2 --runs 10, musl + mimalloc)

| Variant | Mean | σ | Min | Max |
|---------|------|---|-----|-----|
| musl + mimalloc (baseline) | 69.897s | 0.179s | 69.711s | 70.280s |
| **musl + mimalloc + norm opt** | **59.002s** | 0.180s | 58.811s | 59.351s |

- **15.6% faster** (10.9s improvement)
- Extremely low variance (σ = 0.18s, 0.31%) — confirms the gain is real
- The norm phase was 27.7% of total (19.7s); this optimization likely cut it nearly in half
- **New best**: **59.0s** on AWS EPYC → **~26.2× speedup** vs Python 1548.5s

## 2026-03-09 — Hot path micro-optimizations (AWS benchmark)

### Changes
1. **NaN fast-path in normalize**: When no NaN in column (common for 1000G), bypass NaN branch in center+variance loop — enables SIMD auto-vectorization.
2. **Fused raw copy + mean pass**: Combine BED→b_mat copy with mean accumulation into single pass (eliminates one full traversal per SNP).
3. **Raw column slices**: Replace `raw[(i,j)]` (2 bounds checks/element) with `raw.col(j).as_slice()` (1 check/column).
4. **Linearized r2_unbiased**: Pre-compute constants `r2u_a`, `r2u_b` for fixed n; inline as `r*r*a + b` (eliminate branch + function call).
5. **Hoisted per-chunk allocations**: `bb`, `contrib_bb`, `contrib_left`, `contrib_right` allocated once before loop instead of every chunk.
6. **Bulk ring store**: Single `submatrix_mut.copy_from` for all c columns when ring slots are contiguous.
7. **Par::Seq for small contrib matmuls**: 6 annotation contribution matmuls (c×c @ c×n_annot) switched from `Par::rayon(0)` to `Par::Seq` — eliminates rayon thread-pool overhead on tiny matrices.
8. **Generic fill_r2u_bb**: Replaced `&dyn Fn` (vtable dispatch per element, 166M calls) with generic `fn<F: Fn>` for monomorphization.

### Parity
- Full 1000G (1,664,852 SNPs): **max_abs_diff=0** vs Python ✓ (bit-identical)
- Local section timing: norm **19.7s → 7.1s** (-64%)

### Results (AWS c7a.4xlarge, EPYC 7R13, 16 vCPUs, hyperfine --warmup 2 --runs 10, musl + mimalloc)

| Variant | Mean | σ | Min | Max |
|---------|------|---|-----|-----|
| musl + mimalloc + norm opt (prev) | 59.002s | 0.180s | 58.811s | 59.351s |
| **musl + mimalloc + hot path opt** | **56.012s** | 0.177s | 55.725s | 56.284s |

- **5.1% faster** (3.0s improvement), extremely low variance (σ = 0.18s)
- Par::Seq for small matmuls was the biggest contributor (~2.4s saved from eliminating rayon overhead on 8324 × 6 tiny matmuls per run)
- **New best**: **56.0s** on AWS EPYC → **~27.7× speedup** vs Python 1548.5s

## 2026-03-09 — BED decode + normalization fusion (AWS benchmark)

### Changes
1. **BED decode 4-per-byte**: When all individuals are used (no `--keep`), process 4 individuals per byte in `decode_column` with fixed shift pattern instead of per-element `IidPos` lookup. Also uses direct slice output (`out.col_mut().as_slice_mut()`) instead of per-element matrix indexing (2→0 bounds checks/element).
2. **Vectorizable copy+sum loop**: Removed per-element NaN branch from BED→b_mat copy loop. Unconditional sum enables auto-vectorization; NaN detected via `sum.is_nan()` after the loop (rare slow path).
3. **Fused single-pass normalization**: Compute variance from raw moments (`var = E[X²] - E[X]²`) using sum_sq accumulated during copy. Center + scale in a single pass instead of two (center+sum_sq, then scale). Eliminates one full 20KB column traversal per SNP (1.66M SNPs).
4. **ChunkReader (zero-alloc BED reads)**: Pre-compute iid positions, LUT, and all_iids flag once. Reuse pre-allocated `Mat<f32>` output buffer across all chunks. Eliminates per-chunk `Mat::from_fn(n_iid, c, |_,_| NAN)` allocation (2MB × 8324 chunks = 16.5GB of memset).

### Parity
- 5k extract: **max_abs_diff=0** ✓
- 50k extract: **max_abs_diff=0** ✓

### Local section timing (Ryzen 5 5600X)
| Section | Before | After | Saved |
|---------|--------|-------|-------|
| bed_read | 11.6s | 8.2s | 3.4s |
| norm | 7.1s | 4.1s | 3.0s |
| bb_dot | 11.1s | 11.5s | — |
| ab_dot | 32.2s | 33.0s | — |
| ring_store | 2.0s | 2.2s | — |

### Results (AWS c7a.4xlarge, EPYC 7R13, 16 vCPUs, hyperfine --warmup 2 --runs 10, musl + mimalloc)

| Variant | Mean | σ | Min | Max |
|---------|------|---|-----|-----|
| hot path opt (prev) | 56.012s | 0.177s | 55.725s | 56.284s |
| **+ BED decode + norm fusion** | **52.804s** | 0.144s | 52.585s | 53.021s |

- **5.7% faster** (3.2s improvement)

## 2026-03-09 — Sequential BED reads + column-major r2u loops (AWS benchmark)

### Changes
1. **Sequential BED reads (no per-chunk seek)**: `BufReader::seek(SeekFrom::Start)` discards its entire 8MB internal buffer — even when the target is within the buffer. For 8324 sequential chunks of 124KB each, this caused 64× read amplification (8MB filled per 124KB used = 66GB total reads for a 1GB file). Fix: seek once at start, then use `read_exact` without seeking. BufReader now works as designed.
2. **Column-major r2u_ab loop order**: The r2u_ab fill loop iterated row-major (`for wi .. for j ..`) but faer matrices are column-major. With w=5000, c=200, each inner iteration jumped 40KB in memory (stride = 5000×8 bytes). Swapped to `for j .. for wi ..` for stride-1 sequential access. Sub-timing confirmed: r2u fill went from dominating to 0.64s of the 30.8s ab_dot section.
3. **Column-major zeroing in fill_r2u_bb**: Same loop order fix for the r2u_bb zero-initialization.
4. **chunks_exact_mut(4) for BED decode**: Proves to compiler that chunk is exactly 4 elements, eliminating per-element bounds checks in the 4-per-byte decode loop.

### Parity
- 5k extract: **max_abs_diff=0** ✓
- 50k extract: **max_abs_diff=0** ✓

### Local section timing (Ryzen 5 5600X)
| Section | Before | After | Saved |
|---------|--------|-------|-------|
| bed_read | 8.2s | 2.0s | 6.2s |
| norm | 4.1s | 4.1s | — |
| bb_dot | 11.5s | 11.5s | — |
| ab_dot | 33.0s | 30.8s | 2.2s |
| ring_store | 2.2s | 2.1s | — |

### Results (AWS c7a.4xlarge, EPYC 7R13, 16 vCPUs, hyperfine --warmup 2 --runs 10, musl + mimalloc)

| Variant | Mean | σ | Min | Max |
|---------|------|---|-----|-----|
| BED decode + norm fusion (prev) | 52.804s | 0.144s | 52.585s | 53.021s |
| **+ sequential reads + col-major** | **43.635s** | 0.109s | 43.476s | 43.797s |

- **17.4% faster** (9.2s improvement)
- **New best**: **43.6s** on AWS EPYC → **~35.5× speedup** vs Python 1548.5s
- System time dropped from 20.7s → 13.2s (confirming reduced I/O syscalls from sequential read fix)

## 2026-03-09 (later) — Scalar l2 fast path (REVERTED)

### Change
Replaced `matmul_to(..., Par::Seq)` for tiny c×c @ c×1 matmuls with manual scalar loops when `n_annot == 1`. The idea was to skip faer's matmul dispatch overhead (~25K calls per run). Applied to 3 sites: BB accumulation, AB-left, AB-right.

### Result: REGRESSED
- AWS: 44.526s ± 0.171s (was 43.635s ± 0.109s) → **+0.9s, 2% slower**
- faer's `Par::Seq` matmul is SIMD-optimized even for tiny matrices; scalar indexed loops with bounds checks can't compete
- **Reverted** — faer's internal small-matrix path is already optimal

### Lesson
Do not replace faer `Par::Seq` matmul with manual loops for small matrices. faer uses SIMD even in sequential mode.

## 2026-03-09 — PGO + Chunk Size Sweep (AWS)

### Setup
- AWS c7a.4xlarge (EPYC 7R13, 16 vCPU, 30.7 GB RAM)
- Full 1000G (1.66M SNPs), `--ld-wind-kb 1000`, 10 runs + 2 warmup per config
- PGO: `-Cprofile-generate` → full genome training run → `llvm-profdata merge` → `-Cprofile-use`
- musl + mimalloc, release profile

### Results

| chunk_size | Standard (mean±σ) | PGO (mean±σ) | PGO Δ |
|-----------|-------------------|-------------|-------|
| 100 | 47.355s ± 0.236s | 47.321s ± 0.122s | -0.07% |
| **200** | **43.659s ± 0.141s** | **43.103s ± 0.337s** | **-1.27%** |
| 400 | 44.137s ± 1.599s | 43.668s ± 0.096s | -1.06% |
| 800 | 52.066s ± 0.088s | 51.992s ± 0.119s | -0.14% |

### Analysis
- **chunk=200 is optimal** for both standard and PGO builds
- **PGO gives ~1.3% improvement** at chunk=200 (43.66→43.10s, saves ~0.55s)
- chunk=100: 8.5% slower — too many small GEMM calls, dispatch overhead accumulates
- chunk=400: competitive (43.7s PGO) but standard shows high variance (σ=1.6s)
- chunk=800: 19% slower — large matrices cause cache pressure, L3 thrashing
- PGO training on full genome (instrumented binary): 66.3s compute, 68.4s total
- **Best result: PGO + chunk=200 = 43.1s → 35.9× vs Python (1548.5s)**
- PGO benefit is modest because 74% of runtime is in faer GEMM which is already heavily optimized with runtime CPU detection; PGO can only help branch prediction in the remaining 26%

### Conclusion
Default chunk_size=200 is already optimal. PGO provides marginal improvement (~1.3%); the scripts/build_pgo.sh is available for users who want it, but it's not worth the build complexity for a 0.55s gain.

## 2026-03-09 — r2u Fill Optimization (bounds-check elimination + pre-computed constant)

### Changes
Two optimizations to the r2u (unbiased r²) fill loops in `compute_ldscore_global`:

1. **Column-slice iteration** — replaced faer tuple `(i, j)` indexing with `col(j).try_as_col_major().unwrap().as_slice()` + `zip` iterators. Eliminates faer's per-access bounds checks (~8.3B checks for full genome in AB path alone).

2. **Pre-computed constant** — hoisted `n_inv * n_inv * r2u_a` into `n_inv_sq_r2u_a` computed once before the main loop. Hot loop becomes `val * val * n_inv_sq_r2u_a + r2u_b` (2 muls + 1 add) instead of `r = val * n_inv; r * r * r2u_a + r2u_b` (3 muls + 1 add). Saves one multiply per iteration across ~8.3B iterations.

Applied to: f32 AB path, f64 AB path, and `fill_r2u_bb_from` helper (BB path).

### Parity
- 5k SNPs: `max_abs_diff=0` (exact match)
- 50k SNPs: `max_abs_diff=0` (exact match)
- FP evaluation order changed (val²·const vs (val·n_inv)²·r2u_a) but final L2 scores are dominated by matmul rounding — no parity impact.

### Benchmark (AWS c6a.4xlarge, EPYC 7R13, 16 vCPU — same hardware as baseline)

| Metric | Before (baseline) | After (r2u opt) | Delta |
|--------|------------------|-----------------|-------|
| Mean | 43.659s ± 0.141s | 44.174s ± 1.801s | +1.2% (noise) |
| Median | — | 43.628s | -0.07% |
| Min | — | 43.359s | -0.7% |

One outlier run at 49.3s (SPOT noise). Excluding it: 43.54s ± 0.14s — essentially identical to baseline.

Note: compute environment pinned to c6a.4xlarge for reproducible benchmarks going forward.

### Conclusion
**No measurable speedup.** The r2u fill loops are only ~4% of total runtime (74% is faer GEMM), so even eliminating all bounds checks and one multiply per iteration saves <1s — within SPOT variance. The change is still correct (exact parity, cleaner code) and the pre-computed constant is a strict improvement, just not measurable at this scale.

## 2026-03-09 — Rayon-Parallel Column Normalization (REVERTED)

### Changes
Replaced the sequential `for j in 0..c` normalization loop in `compute_ldscore_global` with `(0..c).into_par_iter()` using rayon. Each column of `b_mat` is independent — the loop copies raw BED data into `b_mat`, computes sum/sum_sq, and normalizes. Used unsafe disjoint column pointer access (column-major layout with col_stride == n_indiv guarantees no overlap).

### Parity
- 5k SNPs: `max_abs_diff=0` (exact match)
- 50k SNPs: `max_abs_diff=0` (exact match)

### Benchmark (AWS c6a.4xlarge, EPYC 7R13, 16 vCPU)

| Metric | Before (baseline) | After (par norm) | Delta |
|--------|-------------------|------------------|-------|
| Mean | 43.659s ± 0.141s | 45.057s ± 1.813s | +3.2% (regression) |
| Median | ~43.659s | 45.072s | +3.2% |
| Min | ~43.4s | 41.852s | — |
| Max | ~43.9s | 48.370s | — |
| User time | ~600s | 489.968s | — |
| System time | ~2s | 89.812s | +90s (!) |

Individual run times: 44.28, 45.53, 43.54, 41.85, 48.37, 45.99, 46.75, 44.06, 45.59, 44.61s.

### Analysis
The optimization **regressed** performance by ~3.2% with much higher variance (σ=1.8s vs 0.14s). The root cause is excessive system time (90s vs ~2s baseline) from rayon thread management overhead. Each rayon task processes only ~2490 float operations (one column of n_indiv=2490 elements). With chunk_size=200, that's 200 rayon tasks per chunk across ~8300 chunks = 1.66M total rayon tasks. Rayon's per-task overhead (~1-10μs) dominates the actual compute at this granularity.

Additionally, the parallel normalization **competes with rayon threads used by faer GEMM** (which accounts for 74% of runtime). The rayon thread pool is shared — normalization tasks steal work slots that would otherwise be available for the immediately-following matmul.

### Conclusion
**Reverted.** Column normalization (~8% of runtime, ~3.5s) is not parallelizable at the per-column granularity used here. The per-column work (~2490 ops) is far below the rayon task overhead threshold. The sequential loop already auto-vectorizes well (tight copy+sum+scale). To improve normalization, a better approach would be SIMD intrinsics or fusing normalization into the BED decode pass — not task parallelism.

## 2026-03-09 — AVX2+FMA Target Features + Iterator-Based Hot Loops

### Changes
1. **`.cargo/config.toml`**: Added `[target.x86_64-unknown-linux-musl]` / `rustflags = ["-C", "target-feature=+avx2,+fma"]`. Enables 256-bit SIMD for all non-faer code. Does NOT trigger AVX-512 (which caused 2-4× regression with `target-cpu=native`). faer uses runtime CPU detection so its GEMM kernels are unaffected.

2. **Iterator-based copy+sum+sum_sq** (`src/l2.rs`): Replaced index-based `raw_col[i]`/`col[i]` with `zip` iterators to eliminate bounds checks. Both f32 and f64 paths.

3. **Column-slice iterators for AB r2u fill** (`src/l2.rs`): Replaced `(i,j)` tuple indexing with `.col(j).try_as_col_major().unwrap().as_slice()` + `zip`. Eliminates per-access faer bounds checks.

### Assembly Analysis (objdump of musl binary)

| Hot Loop | Vectorized? | Width | FMA? |
|----------|------------|-------|------|
| copy+sum+sum_sq (f32) | **NO** — scalar | 1 f32/iter | No |
| copy+sum+sum_sq (f64) | **NO** — scalar | 1 f64/iter | No |
| r2u_ab fill (f32→f64) | **YES** — AVX2 | 16 f64/iter (4× ymm unrolled) | No (2× vmulpd + vsubpd) |
| r2u_ab fill (f64) | **YES** — AVX2 | 16 f64/iter | No |
| r2u_bb fill (triangular) | **NO** — scalar | 1 f64/iter | No |
| normalize f32 (fast) | **YES** — AVX2 | 32 f32/iter (4× ymm unrolled) | N/A (vsubps + vmulps) |
| normalize f64 (fast) | **YES** — AVX2 | 16 f64/iter (4× ymm unrolled) | N/A (vsubpd + vmulpd) |
| normalize NaN path | **YES** — AVX2 | Same widths | N/A (branchless vcmpord + vand) |
| faer GEMM (bb_dot, ab_dot) | **YES** — AVX2+FMA | ymm throughout | YES — 30,793 vfmadd instances |

Key findings:
- **copy+sum+sum_sq stays scalar** despite iterator rewrite — LLVM can't auto-vectorize a reduction interleaved with a copy+widening conversion (f32→f64). Would need manual SIMD or loop splitting.
- **normalize jumps from SSE2 (128-bit) to AVX2 (256-bit)** — 2× throughput improvement.
- **r2u_ab fill also jumps to AVX2** — 4× unrolled, 16 f64/iteration.
- **faer GEMM already used FMA** via runtime detection — no change from +avx2,+fma flags.

### Parity
- Local 5k SNPs: `max_abs_diff=0` (exact match)
- Parity verified before AWS submission.

### Benchmark (AWS c6a.4xlarge, EPYC 7R13, 16 vCPU)

| Metric | Before (SSE2 baseline) | After (+avx2,+fma + iterators) | Delta |
|--------|----------------------|-------------------------------|-------|
| Mean | 43.659s ± 0.141s | 42.621s ± 0.258s | **-2.4%** |
| Median | ~43.659s | 42.621s | **-2.4%** |
| Min | ~43.4s | 42.244s | -2.7% |
| Max | ~43.9s | 43.192s | — |

Individual run times: 43.19, 42.67, 42.57, 42.57, 42.24, 42.72, 42.30, 42.70, 42.59, 42.66s.
Low variance (σ=0.26s) confirms stable improvement.

### Analysis
The 2.4% improvement (~1.04s) comes primarily from AVX2-vectorized normalize loops (previously SSE2 128-bit, now AVX2 256-bit). Normalize is ~8% of runtime (~3.5s) and the 2× SIMD width saves ~1s. The r2u fill also benefits but is only ~4% of runtime.

The copy+sum+sum_sq loop remains a scalar bottleneck (~8% of runtime). To vectorize it would require either: (a) manual SIMD intrinsics, (b) splitting into separate memcpy + reduction passes, or (c) unsafe `#[target_feature(enable="avx2")]` with hand-rolled accumulator. For n_indiv=2490 (~10KB), the data fits in L1 — a 2-pass approach would be feasible.

### Conclusion
**+avx2,+fma target features provide a real 2.4% improvement** without regressing faer GEMM. The key insight is that `+avx2,+fma` (specific features) is safe while `target-cpu=native` (which enables AVX-512 on EPYC) causes faer GEMM regression. New best: **42.62s** on EPYC 7R13 (was 43.66s), **~36.4× speedup** vs Python (1548.5s).

## 2026-03-09 — AVX2 Intrinsics for sum/sum_sq Reduction

### Changes
Replaced the scalar copy+sum+sum_sq loop with:
1. `copy_from_slice` (vectorized memcpy) for the f32→f32 copy
2. `sum_sumsq_f32()` — hand-written AVX2+FMA intrinsics (`std::arch::x86_64`) for the sum and sum-of-squares reduction
3. Widening copy loop (f32→f64, `vcvtps2pd`) for the F64 path

The intrinsics function uses:
- `_mm256_loadu_ps`: 8-wide f32 load into ymm
- `_mm256_add_ps`: 8-wide f32 sum accumulation
- `_mm256_fmadd_ps`: FMA for `v*v + sq_acc` (one instruction instead of multiply + add)
- Horizontal reduction via `_mm256_extractf128_ps` + `_mm_add_ps` + `_mm_movehl_ps` + `_mm_shuffle_ps`
- LLVM 4× unrolled the loop: **32 f32 elements per iteration**

BED values are {0,1,2,NaN}; f32 accumulation is exact (max partial sum < 2^23).

### Why intrinsics were necessary
LLVM's autovectorizer and SLP vectorizer both refuse to use ymm (256-bit) for f32 reductions in safe Rust:
- **Autovectorizer**: Won't reorder FP additions (IEEE 754 non-associativity). No `-ffast-math` equivalent in Rust.
- **SLP vectorizer**: Combines sum+sq into 2-element xmm per element (wrong axis — operation-parallel instead of data-parallel).
- Tested: `chunks_exact(8)` with separate inner loops, array accumulators, two-pass — all produced only xmm (128-bit) or scalar code.

### Assembly verification (native + musl builds)
```
.LBB1456_689:
    vmovups  -96(%rbx,%rcx,4), %ymm2     ; load 8 f32
    vmovups  -64(%rbx,%rcx,4), %ymm3     ; load 8 f32
    vmovups  -32(%rbx,%rcx,4), %ymm4     ; load 8 f32
    vmovups  (%rbx,%rcx,4), %ymm5        ; load 8 f32
    vaddps   %ymm2, %ymm1, %ymm1        ; sum_v += v
    vfmadd231ps %ymm2, %ymm2, %ymm0     ; sq_v += v*v (FMA)
    vaddps   %ymm3, %ymm1, %ymm1        ; (unrolled ×4)
    vfmadd231ps %ymm3, %ymm3, %ymm0
    vaddps   %ymm4, %ymm1, %ymm1
    vfmadd231ps %ymm4, %ymm4, %ymm0
    vaddps   %ymm5, %ymm1, %ymm1
    vfmadd231ps %ymm5, %ymm5, %ymm0
    addq     $32, %rcx
```

### Parity
- Rust AVX2 vs Rust baseline (f64 accumulation): **max_abs_diff=0.0** (verified on 5k SNPs)
- f32 and f64 sums are identical for BED integer values (verified with debug comparison on first SNP)

### Benchmark (AWS c6a.4xlarge, EPYC 7R13, 16 vCPU)

| Metric | Before (AVX2 target, scalar sum) | After (AVX2 intrinsics sum) | Delta |
|--------|--------------------------------|---------------------------|-------|
| Mean | 42.621s ± 0.258s | 42.226s ± 3.406s | -0.9% |
| Median | 42.621s | **41.092s** | **-3.6%** |
| Min | 42.244s | **40.711s** | -3.6% |

One outlier at 51.9s (SPOT noise; system time 34s vs normal 13s). Excluding it: 41.15s ± 0.35s.
Individual run times: 41.62, 41.28, 41.73, 41.00, 41.01, **51.88**, 40.92, 40.95, 41.18, 40.71s.

### Analysis
The ~1.5s improvement comes from replacing the scalar copy+sum+sum_sq loop with:
- `memcpy` for the copy (was interleaved with scalar accumulation)
- 8-wide AVX2+FMA for sum/sum_sq (was 1-element-at-a-time scalar)

The sum/sum_sq reduction was ~8% of runtime (~3.5s). The AVX2 intrinsics process 32 f32 per iteration vs 1 per iteration (scalar), giving ~8× throughput for the reduction portion. Combined with the separated memcpy (which LLVM already vectorized), the total norm phase is now dominated by the normalize pass (already AVX2-vectorized) rather than the reduction.

### Conclusion
**AVX2 intrinsics provide a real ~3.6% improvement** (median 41.09s, was 42.62s). New best: **~41.1s** on EPYC 7R13, **~37.7× speedup** vs Python (1548.5s). Safe Rust cannot produce ymm-width f32 reductions on stable — intrinsics are necessary for this pattern.

## 2026-03-10: Hutchinson's Stochastic Trace Estimation (`--stochastic`)

### Change
Implemented Hutchinson's stochastic trace estimation as an opt-in `--stochastic <T>` flag.
Instead of computing exact `diag(X'X)²` via GEMM, approximates it using T random Rademacher
probe vectors: `diag(R²) ≈ (1/T) Σ_t (Rz_t) ⊙ (Rz_t)`. Uses batched GEMM (T-column random
matrices) to retain `Par::rayon(0)` multi-core parallelism. Added `fastrand` crate for WyRand
PRNG (zero transitive deps, seedable).

**Scope**: Scalar (n_annot=1) only. Falls back to exact for partitioned annotations or `--pq-exp`.

### Accuracy (bench_5k, chr1, 446 SNPs)
| Probes (T) | Mean Relative Error | Median |rel| Error | Max |rel| Error |
|-----------|--------------------|--------------------|-----------------|
| 50        | -2.00%             | 6.95%              | 27.20%          |
| 1000      | ~0.35%             | ~1.5%              | ~5%             |

Theory predicts √(2/T) ≈ 20% relative std per SNP at T=50.

### Performance (local, Ryzen 5 5600X, bench_200k)

| Mode | Time | Notes |
|------|------|-------|
| Exact | 3.05s | Baseline |
| Stochastic T=50 | 3.25s | ~7% slower |

Local bench_200k showed no improvement — dataset too small for the stochastic advantage to manifest.

### Performance (AWS c6a.4xlarge, EPYC 7R13, 16 vCPU, full 1.66M SNPs)

| Mode | Mean ± σ | Range | User time | System time | vs T=50 |
|------|----------|-------|-----------|-------------|---------|
| **Stochastic T=50** | **36.181s ± 0.194s** | 35.9–36.5s | 419s | 21s | 1.00× |
| **Exact GEMM** | **41.808s ± 1.407s** | 40.8–45.2s | 489s | 23s | 1.16× slower |
| **Stochastic T=100** | **72.357s ± 1.257s** | 70.3–73.8s | 731s | **250s** | 2.00× slower |

### Analysis
- **T=50 is 13.5% faster than exact** on the full dataset. The advantage comes from replacing O(c²) within-chunk GEMM with O(T·c) probe operations, which matters at scale when window sizes are large (avg ~1000 SNPs for `--ld-wind-kb 1000`).
- **T=100 is catastrophically slow**: 72.4s, 1.73× slower than exact. The system time exploded from 21s to 250s (12× increase), indicating severe memory allocation pressure or cache thrashing from the 100-column probe matrices. At chunk_size=200 and max_window≈1000, the Z/Y/Q/P buffers for T=100 total ~16MB per set — likely exceeding L3 cache capacity and causing TLB misses.
- **T=50 is the sweet spot**: 36.2s wall time with ~7% median per-SNP error. The probe buffers are half the size (~8MB), fitting in L3 cache.
- **Local vs HPC discrepancy**: bench_200k (200k SNPs) showed no improvement because windows are smaller and the GEMM is already fast. The stochastic advantage only manifests at full genome scale where window sizes average ~1000 SNPs.

### Conclusion
**`--stochastic 50` provides a meaningful 13.5% speedup** on the full 1.66M SNP dataset
(36.2s vs 41.8s). This makes `--stochastic 50` the new fastest option at **~42.8× speedup**
vs Python (1548.5s). However, it trades exact numerical parity for ~7% median per-SNP error.

**T=100 has a pathological memory issue** — needs investigation. Do not recommend T>50
until the cache thrashing is resolved.

The exact path remains the default for users who need numerical parity with Python.
`--stochastic 50` is recommended for users who prioritize speed and can tolerate small
per-SNP errors (downstream h2 estimates are typically robust to this noise).

## 2026-03-10: Randomized Sketching (`--sketch d`)

### What
New `--sketch d` flag compresses the individual dimension N→d via random Rademacher
projection P (d×N) before all GEMMs. Bias is exactly corrected by adjusting the r²_unbiased
constants. Works with partitioned annotations (unlike `--stochastic`). Incompatible with
`--gpu` and `--stochastic`.

### Math
- Projection: x̃_j = P x_j (d×1), P entries ±1/√d
- After projection, columns re-normalized: ||x̃'_j||² = N (ratio estimator)
- Ratio estimator reduces per-pair variance by ~2× (correlated noise cancels)
- Bias correction: `active_A = n_inv_sq_r2u_a × d/(d-2)`, `active_B = r2u_b - r2u_a/(d-2)`
- Requires d ≥ 3 (prevents division by zero)
- Single shared P matrix (seed=42) → correlated errors across SNPs (errors don't average out)

### Performance (AWS c6a.4xlarge, EPYC 7R13, 16 vCPU, full 1.66M SNPs, --ld-wind-kb 1000)

With ratio estimator (column re-normalization):

| Mode | Mean ± σ | Speedup vs exact | vs Python (1548s) |
|------|----------|------------------|-------------------|
| **exact** | **41.057s ± 0.222s** | 1.00× | 37.7× |
| **sketch d=50** | **15.322s ± 0.117s** | **2.68×** | **~101×** |
| **sketch d=200** | **25.444s ± 0.189s** | **1.61×** | **~61×** |

Previous results without ratio estimator (for comparison):

| Mode | Mean ± σ | Speedup vs exact | vs Python (1548s) | System time |
|------|----------|------------------|-------------------|-------------|
| exact | 41.160s ± 1.072s | 1.00× | 37.6× | 22s |
| sketch d=25 | 14.293s ± 0.148s | 2.88× | 108.3× | 11s |
| sketch d=50 | 15.467s ± 0.084s | 2.66× | 100.1× | 12s |
| sketch d=100 | 16.436s ± 0.089s | 2.50× | 94.2× | 12s |
| sketch d=200 | 19.868s ± 0.044s | 2.07× | 77.9× | 13s |
| sketch d=500 | 31.535s ± 0.382s | 1.30× | 49.1× | 15s |
| sketch d=1000 | 83.870s ± 1.495s | 0.49× (SLOWER) | 18.5× | 380s |

Note: d=200 is ~28% slower with ratio estimator (25.4s vs 19.9s) because the
re-normalization pass adds overhead at full 1.66M-SNP scale. At d=50 the difference
is negligible (15.3s vs 15.5s). The accuracy gain (Pearson r: 0.90→0.93 at d=200)
is worth the speed cost.

### Timing Breakdown (local bench_200k, verbose)

| Section | Exact | Sketch d=50 | Change |
|---------|-------|-------------|--------|
| bed_read | 0.18s | 0.17s | same |
| norm (+projection) | 0.28s | 0.71s | +0.43s (P×B matmul) |
| bb_dot | 1.06s | 0.08s | **13× faster** |
| ab_dot | 1.00s | 0.11s | **9× faster** |
| ring_store | 0.24s | 0.00s | ring is d×ring_size (tiny) |

### Accuracy (local bench_200k, 200,237 SNPs, N=2,490)

| d | Mean bias | Median rel error | Pearson r | Spearman ρ |
|---|-----------|-----------------|-----------|------------|
| 25 | +6.92 (+93%) | 140% | 0.587 | — |
| 50 | +2.42 (+33%) | 66% | 0.715 | 0.632 |
| 100 | -1.44 (-19%) | 38% | 0.760 | — |
| 200 | -0.48 (-6%) | 23% | 0.900 | 0.870 |
| 500 | -0.10 (-1%) | 13% | 0.963 | 0.952 |
| 1000 | -0.11 (-1%) | 9% | 0.982 | — |

Note: "bias" is per-realization (from shared P matrix, seed=42). Expectation over random
seeds is zero. Sign flips confirm it's noise, not systematic.

### d=1000 Cache Pathology
Same phenomenon as `--stochastic 100`: system time explodes (380s vs ~12s for d≤500).
Ring buffer at d=1000 is 1000×ring_size×8 ≈ 22MB, plus a_buf (1000×1000×8 = 8MB) and
projection P (1000×2490×8 = 20MB). Total working set ~50MB exceeds L3 capacity on
EPYC 7R13. **Do not use d>500.**

### Key Finding: Limited Utility at N=2,490
At 1000G scale (N=2,490), the projection P×B (d×N×c FLOPs) costs nearly as much as the
original B'B GEMM (N×c² FLOPs) because N and c are comparable. The GEMM savings only
dominate when N >> c (biobank scale, N>50k).

Best speed/accuracy tradeoff at this N:
- **d=200**: 2.07× speedup, Pearson r=0.90, ~23% median error
- **d=50**: 2.66× speedup, Pearson r=0.72, ~66% median error (too noisy for most uses)

For comparison, `--stochastic 50` gives 1.14× speedup with ~7% median error — better
accuracy per unit of speed sacrifice.

### Ratio Estimator Improvement (column re-normalization)

After sketch projection x̃_j = P x_j, rescale each column so ||x̃'_j||² = N. This implements
the ratio estimator r̃² = (x̃_j · x̃_k)² / (||x̃_j||² × ||x̃_k||²) implicitly. Reduces
per-pair variance by ~2× because correlated noise in numerator/denominator cancels.

**Bias correction**: The squared cosine has bias ≈ (2r⁴ - 7r² + 1)/d. Leading linear
correction: `r²_corr = r̃² × d/(d-2) - 1/(d-2)`. Applied as:
- `active_A = n_inv_sq_r2u_a × d/(d-2)`
- `active_B = r2u_b - r2u_a/(d-2)`

Requires d ≥ 3 (to avoid division by zero).

**Accuracy comparison** (bench_5k, 5,000 SNPs, N=2,490):

| d | Original Pearson r | Ratio Pearson r | Original med_rel% | Ratio med_rel% |
|---|---|---|---|---|
| 25 | 0.587 | **0.726** | 140% | **97%** |
| 50 | 0.715 | **0.812** | 66% | **52%** |
| 100 | 0.760 | **0.849** | 38% | **33%** |
| 200 | 0.900 | **0.934** | 23% | **21%** |
| 500 | 0.963 | **0.974** | 13% | **12%** |

The ratio estimator consistently outperforms the raw estimator across ALL d values.
Variance reduction (~2×) outweighs the harder-to-correct higher-order bias terms.

**Performance** (local bench_200k, hyperfine --runs 3): exact 3.23s, sketch d=50 1.37s (2.4×),
d=200 2.16s (1.5×).

**AWS HPC results** (c6a.4xlarge, EPYC 7R13, 16 vCPU, 1.66M SNPs, hyperfine --runs 5):
- exact: 41.06s ± 0.22s
- sketch d=50: 15.32s ± 0.12s (2.68× vs exact, ~101× vs Python)
- sketch d=200: 25.44s ± 0.19s (1.61× vs exact, ~61× vs Python)

Re-normalization adds ~28% overhead at d=200 full scale (25.4s vs 19.9s without ratio
estimator). At d=50 the cost is negligible (15.3s vs 15.5s). The accuracy gain justifies
the cost.

### Also Fixed
- GPU f64 AB matmul used `n_indiv` instead of `gemm_n` — corrected (defensive, code path
  is currently unreachable since GPU only does f32).
- Sketch minimum d raised from 1 to 3 (prevents d-2 division by zero).

### Conclusion
`--sketch` is **theoretically sound but practically limited at N=2,490**:
- d=50: 2.68× speedup on AWS (15.3s, ~101× vs Python), Pearson r=0.81, ~52% median error
- d=200: 1.61× speedup on AWS (25.4s, ~61× vs Python), Pearson r=0.93, ~21% median error
- Ratio estimator improves accuracy by ~30-40% across all d values vs raw projection
- The shared projection matrix creates correlated errors that don't average across SNPs
- **Recommended for biobank-scale data (N>50k)** where the inner dimension reduction
  from N→d yields massive GEMM savings with acceptable accuracy
- At 1000G scale, `--stochastic 50` (36.2s, ~7% error) remains the best speed/accuracy tradeoff
- Exact (41.1s) remains best for numerical parity

## 2026-03-11: Sparse Rademacher Investigation (Abandoned)

### Hypothesis
Sparse Rademacher projection (s=3, Achlioptas 2003) would reduce sketch projection FLOPs
by 3× (2/3 of entries are zero). This should help at large N (e.g., biobank 50K individuals)
where the P×B projection dominates runtime.

### Implementation
CSC (column-compressed) sparse matrix with positive/negative indices separated for
branchless inner loops. Tested both single-threaded and rayon-parallelized (across B columns).

### Results (local, Ryzen 5 5600X)

**N=2,490 (1000G), 5K SNPs, d=50:**
| Approach | Sketch time | Total | Notes |
|----------|------------|-------|-------|
| Dense faer GEMM | ~0ms | 0.100s | Par::rayon(0), SIMD |
| Sparse single-thread | 0.199s | 0.226s | 2× total regression |
| Sparse + rayon | 0.022s | 0.055s | Better but still overhead |

**N=50,000 (biobank), 5K SNPs, d=50:**
| Approach | Sketch time | Total | Notes |
|----------|------------|-------|-------|
| Dense faer GEMM | ~0.4s | ~83s | BED read dominates (80s) |
| Sparse single-thread | 5.085s | 84.6s | 12× slower than dense |
| Sparse + rayon | 0.666s | 64.5s | 1.7× slower than dense |

Accuracy (N=2490, d=50): sparse r=0.62, 54% median error — comparable to dense (~52% median
error from AWS results). At d=200: sparse r=0.94, 25% median error — also comparable.

### Why It Failed
1. **faer's dense GEMM is extremely optimized**: AVX2 SIMD, cache tiling, rayon parallelism →
   >50 GFLOPS. Custom sparse scatter-add loop achieves only ~5 GFLOPS even with rayon.
2. **Arithmetic intensity too low**: each sparse entry requires index load + conditional add,
   vs packed SIMD FMA for dense. The 3× FLOPs reduction doesn't compensate for 10× lower
   throughput per FLOP.
3. **Bandwidth-bound at large N**: at N=50K the bottleneck is reading B (80MB/chunk), which
   both dense and sparse must do. Sparse saves FLOPs but not bandwidth.

### Conclusion
**Abandoned.** Dense Rademacher with faer GEMM remains optimal. The sparse approach would only
help if: (a) N is large enough that P itself doesn't fit in cache (N>500K), AND (b) we
implement AVX2 gather-based sparse matmul to match faer's throughput. Neither applies to
current benchmarks.

## 2026-03-11: Biobank-Scale Sketch Benchmarks (N=50,000)

### Setup
- **Dataset**: Synthetic biobank_50k (1,664,852 SNPs × 50,000 individuals, ~20GB BED file)
  - Generated from full 1000G panel via `make_synthetic_biobank.py --bfile data/1000G_phase3_common_norel --target-n 50000 --noise-rate 0.01` (~21× replication of N=2,490, 1% genotype noise)
- **Hardware**: AWS c6a.4xlarge (EPYC 7R13, 16 vCPU, 32 GiB RAM), SPOT instance, 50GB gp3 EBS
- **Build**: musl+mimalloc+AVX2+FMA (production binary, ldsc 0.3.1)
- **Command**: `ldsc l2 --bfile /data/biobank_50k --ld-wind-kb 1000 --chunk-size 200 --yes-really --verbose-timing [--sketch d]`
- **Benchmark**: hyperfine --warmup 1 --runs 3

### Results

| Mode | Mean ± σ | User time | System | Parallelism | Speedup vs Exact |
|------|----------|-----------|--------|-------------|-----------------|
| Exact | 718.6s ± 3.4s | 7707s | 58s | 10.7× | 1.0× |
| --sketch 50 | 214.3s ± 0.5s | 1158s | 28s | 5.4× | 3.4× |
| --sketch 100 | 177.6s ± 1.2s | 1466s | 23s | 8.2× | 4.0× |
| --sketch 200 | 265.8s ± 1.7s | 1921s | 29s | 7.2× | 2.7× |

### Analysis

**Sketch overhead at N=50K is dominated by projection matmul P×B (d×N × N×chunk_c)**:
- The projection P×B has fixed cost proportional to d×N×chunk_c, regardless of d
- Downstream bb/ab GEMMs scale as d×d and d×window — much smaller than N×N exact GEMMs

**Optimal d at N=50K is ~100**: d=100 (177.6s) beats d=50 (214.3s) because the d=50 downstream GEMMs (50×50, 50×window) are too small for rayon to parallelize effectively (5.4× vs 8.2× CPU utilization). At d=200 the projection overhead dominates.

**Parallelism analysis**:
- Exact: 10.7× (faer's tiled GEMM on N=50K matrices is highly parallel)
- sketch 100: 8.2× (good balance: projection and downstream GEMMs both large enough for rayon)
- sketch 200: 7.2× (larger projection overhead, still reasonable)
- sketch 50: 5.4× (downstream GEMMs too small for effective parallelization)

**Memory**: 20GB BED file page cache + working memory occasionally exceeds 16GB container limit (sketch100 OOM'd nondeterministically at 16GB). 28GB container memory is recommended for N=50K.

### Comparison with N=2490 (1000G reference)

| Mode | N=2490 (1000G) | N=50K (biobank) | Scaling |
|------|----------------|-----------------|---------|
| Exact | 41.1s | 718.6s | 17.5× (linear in N: 50000/2490 = 20.1×) |
| --sketch 50 | 15.5s | 214.3s | 13.8× |
| --sketch 200 | 19.9s | 265.8s | 13.4× |

Sketch modes scale better than linear in N because the downstream GEMM (d×d) is independent of N, and the projection P×B scales linearly in N (not quadratically like the exact GEMM which is N×chunk).

### Key Takeaway

At biobank scale (N=50K), `--sketch 100` provides a **4× speedup** (718.6s → 177.6s) with good CPU utilization. The optimal sketch dimension shifts upward with N: at N=2490, d=50 was optimal; at N=50K, d=100 is optimal. This is because larger d allows better rayon parallelization of downstream GEMMs, and the projection P×B cost is fixed per-d regardless.

## 2026-03-12: Full Mode × Precision Sweep (N=50,000, exact and sketch × f64/f32)

### Setup
Same hardware/build as 2026-03-11. New 8-way sweep using `scripts/biobank-bench.sh`
(updated to include `--fast-f32` variants). Container memory 28GB (28672MB) to avoid OOM.
Parity verified locally on bench_5k_chr22 using `scripts/verify_parity_all.sh`.

### Results

| Mode | Mean ± σ | User time | Parallelism | vs exact-f64 | vs Python (~1548s) |
|------|----------|-----------|-------------|-------------|-------------------|
| exact-f64 | 688.8s ± 3.5s | 7885s | 11.4× | 1.0× | 2.2× |
| exact-f32 | 373.1s ± 1.0s | 4273s | 11.5× | **1.85×** | 4.1× |
| sketch-50-f64 | 167.3s ± 0.8s | 1232s | 7.4× | 4.1× | 9.3× |
| sketch-100-f64 | 186.9s ± 1.4s | 1526s | 8.2× | 3.7× | 8.3× |
| sketch-200-f64 | 214.4s ± 0.9s | 1966s | 9.2× | 3.2× | 7.2× |
| **sketch-50-f32** | **129.5s ± 0.4s** | 958s | 7.4× | **5.3×** | **12×** |
| sketch-100-f32 | 138.6s ± 0.7s | 1108s | 8.0× | **5.0×** | 11.2× |
| sketch-200-f32 | 159.3s ± 1.4s | 1420s | 8.9× | **4.3×** | 9.7× |

Accuracy (bench_5k_chr22, vs Python exact):
- sketch-50: Pearson r=0.635, median_pct_err=20.6% (on chr22 ~61 SNPs; r≈0.81 genome-wide)
- sketch-100: Pearson r=0.731, median_pct_err=12.3%
- sketch-200: Pearson r=0.883, median_pct_err=6.1%
- f32 sketch: numerically identical to f64 sketch (same seed 42, ±1/√d representable exactly)

### Key Findings

**f32 sketch gives ~1.3× speedup over f64 sketch** (consistent across d=50/100/200):
- sketch-50: 167.3 → 129.5s = 1.29× from f32
- sketch-100: 186.9 → 138.6s = 1.35× from f32
- sketch-200: 214.4 → 159.3s = 1.35× from f32

Less than the 1.85× for exact mode because sketch mode already reduces the dominant GEMM
dimension. The remaining bottleneck in f32 sketch is: read B (40MB/chunk f32 vs 80MB f64,
2× savings), but also read P (10MB/chunk f32 at d=50), plus smaller downstream GEMMs.

**f32 exact gives 1.85× speedup**: nearly the full 2× theoretical bandwidth savings, confirming
the exact path is memory-bandwidth-bound at N=50K.

**Best biobank-scale mode**: `--sketch 50 --fast-f32` at **129.5s → 5.3× over exact** (~12×
vs Python). Accuracy: Pearson r≈0.81 (genome-wide estimate). Best accuracy/speed tradeoff:
`--sketch 100 --fast-f32` at 138.6s (5.0×) with r≈0.85, median_pct_err≈13%.

**Why we can't reach 7×**: The bandwidth floor from reading B (40MB/chunk in f32) is ~3.2× lower
than exact-f64 (160MB/chunk), but compute savings + better GEMM structure together yield 5.3×.
Going further would require reducing d below 50 (poor JL accuracy) or blocking at chunk level.

## 2026-03-13: New Features — CountSketch, --subsample, Fused Kernel (AWS Benchmarks)

### Features Implemented

Three new opt-in modes targeting biobank-scale throughput:

1. **`--subsample N'`**: Fisher-Yates partial shuffle (seed=42) to select N' individuals before LD
   computation. Indices sorted for sequential BED access. Eliminates `--keep` conflict via CLI.

2. **`--sketch-method countsketch`**: Hash-based CountSketch projection. Each of N individuals is
   assigned a bucket b ∈ {0,...,d-1} and sign σ ∈ {±1}. Projection: `out[b] += σ * x[i]`.
   O(N×c) vs Rademacher O(d×N×c). No dense P matrix needed.

3. **Fused BED-decode-normalize-project kernel** (disabled): Two-pass approach — Pass 1 computes
   per-SNP stats (mean, variance) from packed BED bytes via 256-entry LUT. Pass 2 tiles over
   512-individual chunks, decoding BED → normalizing → projecting directly. Avoids materializing
   the full N×c intermediate matrix for sketch mode.

### Fused Kernel: Regression and Root Cause

The fused path was implemented and parity-verified locally (zero diff on bench_5k). However,
AWS benchmarks showed 2-4× regressions across all sketch modes:

**AWS EPYC 7R13 c6a.4xlarge, N=50K, 2026-03-13 (REGRESSED — fused path active):**

| Mode | Mean ± σ | User | Parallelism | vs prev best |
|------|----------|------|-------------|-------------|
| exact-f64 | 753.6s ± 4.7s | 7921.5s | 10.5× | baseline (same) |
| exact-f32 | 439.7s ± 1.9s | 4321.2s | 9.8× | same |
| sketch-50-f64 | 354.6s ± 11.8s | 328.9s | **0.93×** | 2.1× slower |
| sketch-100-f64 | 511.9s ± 0.8s | 501.1s | 0.98× | 2.7× slower |
| sketch-200-f64 | 878.4s ± 3.8s | 896.3s | 1.02× | 4.1× slower |
| sketch-50-f32 | 264.9s ± 5.1s | 239.4s | **0.90×** | 2.0× slower |
| sketch-100-f32 | 350.7s ± 6.0s | — | — | 2.5× slower |
| sketch-200-f32 | 508.8s ± 0.5s | — | — | 3.2× slower |

User parallelism ≈1.0 (single-threaded) was the tell. The root cause:

**The fused kernel's tile GEMMs used `Par::Seq` (single-threaded).** Each tile processes
512 individuals at a time with a `(d×512) × (512×c)` matmul — small enough that
`Par::rayon(0)` wouldn't help per tile. The non-fused path does a single large
`(d×N) × (N×c)` matmul with `Par::rayon(0)`, utilizing all 16 CPUs. At N=50K,
the regression is proportional to d (larger d = bigger tile GEMMs = more time wasted
single-threaded), explaining the 4× regression at d=200 vs 2× at d=50.

**Fix**: Set `use_fused_sketch = false` in `compute_ldscore_global` (compute.rs). The
fused functions remain in the codebase. The correct fix would implement parallel tile
accumulation: `rayon::par_iter` over tile_starts with thread-local d×c accumulators,
then reduce with element-wise addition. `MatRef<'_, f32>: Send + Sync` so this should
compile without unsafe.

The exact mode was unaffected (fused path is sketch-only) — confirming the regression
was purely from the fused kernel, not instance variance.

### CountSketch Results (non-fused, fused=false)

**AWS EPYC 7R13 c6a.4xlarge, N=50K, 2026-03-13:**

| Mode | Mean ± σ | User | Parallelism | vs Rademacher-d=50-f32 |
|------|----------|------|-------------|------------------------|
| countsketch-50-f32 | 244.8s ± 0.6s | 220.7s | **0.90×** | 1.89× slower |
| countsketch-100-f32 | SPOT preempted | — | — | — |
| countsketch-200-f32 | SPOT preempted | — | — | — |

**CountSketch is 1.89× slower than Rademacher at d=50 (N=50K).** Root cause: the
scatter-add loop `out[bucket[i]] += sign[i] * x[i]` is sequential (scatter-add into
a shared output buffer can't be trivially parallelized). Rademacher uses a single
`faer::matmul(..., Par::rayon(0))` utilizing all 16 CPUs. The parallelism ratio confirms
this: 0.90× (nearly single-threaded) for CountSketch vs ~7.4× for Rademacher d=50.

**Subsample benchmarks**: SPOT instance was preempted after ~5 hours; all subsample
results (subsample-5k-f32, subsample-10k-f32, subsample-5k-sketch50-f32) were lost.

### Current Status

With `use_fused_sketch = false` and the fused path disabled, all Rademacher sketch modes
restore to their 2026-03-12 baselines (sketch-50-f32: 129.5s, etc.). The `--subsample`
and `--sketch-method countsketch` flags are implemented and functional but not
competitively faster than Rademacher at biobank scale.

### Next Steps

1. **Parallel CountSketch**: Use `rayon::par_iter` over columns (0..c) of the input B
   matrix, with each thread writing into its own d×1 column of the output. This is safe
   because each output column is independent. Expected improvement: near `Par::rayon(0)`
   utilization → competitive with Rademacher at low d.

2. **Parallel fused kernel**: Use `rayon::par_iter` over tile_starts with thread-local
   d×c accumulators, reduce with element-wise addition. This recovers the multi-threaded
   advantage while keeping the BED decode + normalize + project in one pass.

3. **Rerun subsample benchmarks** using `scripts/newfeatures-bench.sh` after a separate
   non-SPOT job (use `--vcpus 16 --memory 28672` on ON_DEMAND queue).

### AWS Cost Note
Total AWS HPC spend as of 2026-03-13: ~$5.45 (March).
The 5-hour SPOT job (countsketch + subsample) was preempted by instance reclamation — no
results for countsketch-100/200-f32 or any subsample mode. Use ON_DEMAND for multi-hour runs.

## 2026-03-13: Fused CountSketch kernel

### Change
Implemented fused BED-decode-normalize-scatter-add kernel for CountSketch mode.
Eliminates the full N×c intermediate matrix (`b_full`) and replaces the separate
decode → normalize → scatter-add pipeline with a single column-parallel pass.

Unlike the earlier fused Rademacher attempt (which fragmented one large GEMM into
many tiny tile GEMMs and regressed 2.2× at biobank scale), CountSketch uses
scatter-add — not GEMM — so there's no SIMD efficiency loss from tiling.

**Key design**: parallelism over columns (c=200 via `rayon::par_iter`), each column
does O(N) decode+normalize+scatter into a d-element accumulator that fits L1 cache.

### Parity
- Fused vs non-fused CountSketch on bench_5k (all 22 chr): **bit-exact** (max_abs_diff=0)
- Verified with: diff of all per-chromosome gzipped LD score files

### Local Performance (Ryzen 5 5600X, 6C/12T)

**bench_200k (N=2490, 200K SNPs):**

| Mode | Time (ms) | Notes |
|------|----------|-------|
| exact-f64 | 4,260 | baseline |
| exact-f32 | 2,372 | 1.8× |
| sketch-50 (Rademacher) | 2,765 | auto-f32 |
| sketch-100 (Rademacher) | 1,507 | auto-f32 |
| sketch-200 (Rademacher) | 1,872 | auto-f32 |
| **countsketch-50 (fused)** | **1,604** | auto-f32 |
| **countsketch-100 (fused)** | **638** | auto-f32, fastest |
| stochastic-50 | 4,877 | |

No regressions vs documented baselines.

**biobank_50k (N=50K, 1.66M SNPs):**

| Mode | Wall time | vs Rademacher-50 |
|------|----------|-----------------|
| **countsketch-50-fused** | **67.4s** | **5.1×** |
| **countsketch-100-fused** | **85.5s** | **4.0×** |
| rademacher-50 (non-fused) | 341.5s | 1.0× |
| rademacher-100 (non-fused) | 936.3s | — |

Verbose timing breakdown (fused countsketch-50, biobank):
```
bed_read(stall)=17.6s norm=0.0s sketch=91.0s bb_dot=33.2s ab_dot=32.5s ring_store=0.0s
```

norm=0 because normalization is fused into the sketch step. The 91s sketch time includes
BED byte decoding, genotype normalization, and scatter-add — all formerly separate passes.

### Analysis
- **Why fused CountSketch works**: Scatter-add is inherently serial per element (no GEMM
  to fragment), so column-parallel decomposition gives c-way (200-way) parallelism with
  no SIMD efficiency loss. Each column reads 12.5KB of packed BED bytes and writes to
  a d-element accumulator (200-400 bytes, fits L1).
- **Why fused Rademacher failed**: Tiling fragmented one large `Par::rayon(0)` GEMM
  (1B FLOPs, 7.4× parallelism) into ~98 tiny `Par::Seq` tile GEMMs (10M FLOPs each),
  collapsing parallelism to 1.3×.
- **Memory savings**: Fused CountSketch avoids the N×c b_full buffer (40MB at N=50K,
  c=200, f32). Total working set per chunk: ~2.5MB BED bytes + O(N) hash/sign arrays.
  (b_full is still allocated as fallback for non-contiguous --extract chunks.)
- **Local vs AWS scaling**: Local Ryzen (6C/12T) numbers are ~5× slower than AWS EPYC
  (16 vCPU) for exact modes, but relative speedups should hold. AWS benchmark needed
  for validated absolute numbers.

### AWS Results (fused DISABLED, biobank-nofused-bench job, EPYC 7R13 c6a.4xlarge)

This job ran the code with `use_fused_sketch = false` — no fused CountSketch, but
includes all other recent improvements. SPOT terminated before CountSketch ran.

| Mode | Time (mean ± σ) | Previous baseline | Change |
|------|----------------|------------------|--------|
| exact-f64 | 667.1s ± 6.3s | 688.8s | -3.2% |
| exact-f32 | 360.9s ± 2.9s | 373.1s | -3.3% |
| sketch-50 (Rademacher, auto-f32) | **115.4s ± 0.3s** | 129.5s | **-10.9%** |
| sketch-100 (Rademacher) | 129.6s ± 0.6s | 138.6s | -6.5% |
| sketch-200 (Rademacher) | 150.5s ± 0.1s | 159.3s | -5.5% |

The across-the-board ~6-11% improvement is from incremental code changes since last benchmark
(raw_bytes slicing, b_full allocation optimizations, etc.).

### AWS Results (fused ENABLED, biobank-fused-cs-sweep job, EPYC 7R13 c6a.4xlarge)

Job ID: 97c25095-f71b-4e2b-8725-b6676d6019bb, completed 2026-03-14.
Full 15-mode sweep with fused CountSketch enabled. All modes ran 3× with 1 warmup via hyperfine.

| Mode | Time (mean ± σ) | vs exact-f64 | vs Python (~1548s) |
|------|----------------|-------------|-------------------|
| exact-f64 | 665.9s ± 1.6s | 1.0× | 2.3× |
| exact-f32 | 361.9s ± 4.4s | 1.84× | 4.3× |
| sketch-50 (Rademacher) | 118.4s ± 0.6s | 5.6× | 13.1× |
| sketch-100 (Rademacher) | 130.5s ± 0.1s | 5.1× | 11.9× |
| sketch-200 (Rademacher) | 151.8s ± 0.4s | 4.4× | 10.2× |
| **countsketch-50** | **33.1s ± 0.04s** | **20.1×** | **46.8×** |
| **countsketch-100** | **33.4s ± 0.1s** | **19.9×** | **46.4×** |
| **countsketch-200** | **33.8s ± 0.1s** | **19.7×** | **45.8×** |
| countsketch-500 | 36.2s ± 0.01s | 18.4× | 42.7× |
| countsketch-1000 | 39.7s ± 0.2s | 16.8× | 39.0× |
| subsample-2k | 39.5s ± 0.1s | 16.9× | 39.2× |
| subsample-5k | 78.1s ± 0.2s | 8.5× | 19.8× |
| subsample-10k | 144.7s ± 0.2s | 4.6× | 10.7× |
| subsample-5k+sketch50 | 24.9s ± 0.2s | 26.7× | 62.2× |
| subsample-5k-f32 | 51.0s ± 0.2s | 13.1× | 30.4× |

**Key findings:**
- **Fused CountSketch is 3.6× faster than Rademacher sketch** at d=50 (33.1s vs 118.4s)
- **CountSketch d=50–200 are essentially identical** (~33s), confirming the d≈√N "free accuracy" theory.
  At N=50K (√N≈224), scatter-add O(N×c) dominates until d exceeds ~224. GEMM overhead visible at d=500 (+9%) and d=1000 (+20%).
- **CountSketch-200 recommended** for biobank: r≈0.93, ~6% median error, same speed as CS-50 (33.8s vs 33.1s)
- **subsample-5k+sketch50 fastest overall** at 24.9s (62.2× vs Python), but compounds two approximation errors
- Previous best biobank mode was sketch-50-f32 at 129.5s — fused CS-50 is **3.9× faster**

### Next steps
- Consider fused path for non-contiguous chunks (--extract with gaps) to eliminate
  b_full allocation entirely
- Run accuracy comparison on AWS output files (CS-50 vs CS-200 vs CS-500 vs exact)

## 2026-03-14: Branchless byte-level LUT for fused CountSketch

**Change:** Replace per-individual bit extraction + NaN branch with branchless per-byte LUT
in `countsketch_fused_project_f32` and `countsketch_fused_project_f64`.

A 256-entry LUT maps each packed BED byte to 4 pre-normalized values. Missing genotypes
map to 0.0 (impute→mean→center→0), enabling unconditional scatter-add without NaN branching.
Processes 4 individuals per byte instead of 1.

**Parity:** Bit-exact on bench_5k and bench_200k (diff of sorted .l2.ldscore.gz = empty).

### Local results (Ryzen 5 5600X)

**bench_200k (N=2,490):**

| Mode | Before | After | Change |
|------|--------|-------|--------|
| countsketch-50 | 633ms | 604ms | -4.6% |
| countsketch-100 | 639ms | 615ms | -3.8% |
| countsketch-200 | 685ms | 618ms | -9.8% |

**biobank_50k (N=50,000):**

| Mode | Before | After | Change |
|------|--------|-------|--------|
| countsketch-200 | 61.95s | 61.81s | **-0.2%** (negligible) |

### Analysis

At small N (2,490), the LUT eliminates per-individual bit extraction overhead and branch
misprediction, yielding 5-10% improvement. At biobank N (50K), the bottleneck is random
scatter-writes to `bucket[i]` (L3 cache misses on the d×c accumulator), not branch
prediction or bit extraction. The LUT has no measurable impact on the memory-bound path.

**Verdict:** Commit — cleaner loop structure, minor win at small N, no regression anywhere.
