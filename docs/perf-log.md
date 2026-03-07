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
