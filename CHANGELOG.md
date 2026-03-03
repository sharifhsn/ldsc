# Changelog

All notable changes to this project will be documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.3] — 2026-02-23

### Added
- `cts-annot` continuous-annotation binning (Python `--cts-bin`) and `.annot` output.
- Cell-type-specific `h2-cts` support with `.ldcts` inputs and `--print-all-cts`.
- Overlap-aware partitioned h2 via `--overlap-annot` with frequency inputs.
- Windows (MSVC) system-BLAS support via vcpkg (OpenBLAS + Clapack).
- GitHub Release artifacts for Linux/macOS/Windows with checksums.

### Changed
- Release packaging uses system OpenBLAS on Linux/macOS.
- README reorganized for a faster install/usage path.

## [0.1.0] — 2026-02-19

### Added
- `munge-sumstats` — Polars LazyFrame streaming pipeline; full Python API parity
  (column-name overrides, `--signed-sumstats`, `--info-list`, `--nstudy`, `--a1-inc`, etc.)
- `l2` — ring-buffer DGEMM loop; global sequential pass matching Python's
  cross-chromosome LD window behaviour; 7× speedup over Python on 1000G data
  - `--annot` (partitioned LD scores, per-chromosome annot files)
  - `--extract`, `--print-snps`, `--maf`, `--keep`, `--per-allele`
- `h2` — scalar and partitioned heritability; `--two-step`, `--intercept-h2`,
  `--no-intercept`, `--print-coefficients`, liability-scale output
- `rg` — genetic correlation for all trait pairs; `--two-step`, `--intercept-gencov`,
  `--no-intercept`, `--samp-prev`/`--pop-prev`
- `make-annot` — BED-file and gene-set annotation generators
- Statically linked OpenBLAS (no runtime library needed)
- 40 unit tests; integration smoke test for l2

[Unreleased]: https://github.com/sharifhsn/ldsc/compare/v0.1.3...HEAD
[0.1.3]: https://github.com/sharifhsn/ldsc/releases/tag/v0.1.3
[0.1.0]: https://github.com/sharifhsn/ldsc/releases/tag/v0.1.0
