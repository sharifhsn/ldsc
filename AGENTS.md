# Repository Guidelines

## Project Structure & Module Organization
- `src/` holds the Rust implementation. Entry points and subcommands are in `main.rs` and `cli.rs`; core pipelines live in modules like `munge.rs`, `ldscore.rs`, `regressions.rs`, `jackknife.rs`, and `irwls.rs`.
- `tests/` contains integration tests (`tests/integration.rs`) plus fixtures under `tests/fixtures/`.
- `data/` is for local test data (e.g., 1000G reference files used by the smoke test).
- `test_rust.sh` is a repo-level CLI smoke test script.
- `target/` and `web/target/` are build artifacts.

## Build, Test, and Development Commands
- `cargo build --release` builds the production binary at `target/release/ldsc`.
- `cargo test` runs unit tests embedded in `src/*.rs` and standard integration tests.
- `cargo test --test integration` runs only the end-to-end integration suite.
- `bash test_rust.sh` runs a full CLI smoke test; add `--build` to build first. Requires `data/1000G_phase3_common_norel.{bed,bim,fam}`.
- `docker build -t ldsc .` builds the container image (see `Dockerfile`).
- Run `cargo fmt` before committing changes.
- Always run `cargo clippy --release` before committing changes and use it as the feedback loop.

## Key Behavior & Interfaces (from README)
- Subcommands: `munge-sumstats`, `ldscore`, `h2`, `rg`, `make-annot`.
- Inputs accept plain, `.gz`, and `.bz2` for sumstats/ldscores/annots where noted in CLI help.
- `--ref-ld-chr` / `--w-ld-chr` follow the `@` placeholder convention (e.g., `chr@` → `chr22`).

## Coding Style & Naming Conventions
- Rust edition: 2024 (`Cargo.toml`). Use standard `rustfmt` output (4-space indentation).
- Prefer idiomatic Rust naming: `snake_case` for functions/modules, `CamelCase` for types/traits, `SCREAMING_SNAKE_CASE` for constants.
- If adding complex logic, include brief comments for non-obvious invariants (see `src/ldscore.rs` and `src/regressions.rs`).

## Testing Guidelines
- Unit tests live alongside modules under `#[cfg(test)]` (e.g., `src/parse.rs`).
- Integration tests target the compiled binary in `tests/integration.rs`.
- Naming: `test_<behavior>` for unit tests; integration tests should describe the scenario (see `tests/integration.rs`).
- For larger end-to-end checks, prefer `bash test_rust.sh` and document any required data files in the PR.

## Commit & Pull Request Guidelines
- Commit messages follow Conventional Commits style seen in history: `feat:`, `fix:`, `docs:`, `chore:`, `ci:` with optional scopes (e.g., `fix(ldscore): ...`).
- PRs should include:
  - A short summary of changes and rationale.
  - Tests run and their commands (or note if not run).
  - Any data prerequisites (e.g., 1000G files in `data/`).
  - For CLI behavior changes, include example invocations and expected output files.

## Configuration Notes
- Default BLAS: static OpenBLAS build (`blas-openblas-static`).
- CI uses system OpenBLAS: `--no-default-features --features blas-openblas-system` (see `.github/workflows/ci.yml`).
- Runtime tuning flags (lower priority; not heavily battle-tested): `--blas-threads`, `--rayon-threads`, `--polars-threads`.
- `ldsc --version` prints the compiled BLAS backend.
- Building requires Rust >= 1.85 and a C/Fortran toolchain for static OpenBLAS; system BLAS needs `libopenblas-dev` + `pkg-config`.
- Large datasets are not stored in Git; keep local-only files under `data/`.

## CI Notes
- CI workflow: `.github/workflows/ci.yml`. On failure it dumps OpenBLAS build logs.
