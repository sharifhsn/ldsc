# syntax=docker/dockerfile:1
#
# Multi-stage build for ldsc.
#
# Stages
# ------
#   base     — Rust toolchain + build tools + sccache + cargo-chef
#   planner  — generate cargo-chef recipe (dependency fingerprint)
#   builder  — compile deps then source; both stages use BuildKit cache mounts
#   runtime  — minimal Debian image with only the stripped binary
#
# BuildKit cache mounts keep the cargo registry and sccache artifacts on the
# build host between runs, so individual crate recompilations are skipped when
# only application source changes.

# ── base: toolchain + tools ───────────────────────────────────────────────────
FROM rust:1.85 AS base

# OpenBLAS is compiled from source (ndarray-linalg openblas-static feature)
# and statically linked into the final binary.
RUN apt-get update && apt-get install -y --no-install-recommends \
        cmake \
        gfortran \
    && rm -rf /var/lib/apt/lists/*

# Compile sccache and cargo-chef once; cached as a Docker layer.
RUN cargo install sccache --locked && cargo install cargo-chef --locked

ENV RUSTC_WRAPPER=sccache \
    SCCACHE_DIR=/sccache \
    CARGO_INCREMENTAL=0

# ── planner: generate dependency recipe ──────────────────────────────────────
FROM base AS planner
WORKDIR /app
COPY . .
RUN --mount=type=cache,target=/usr/local/cargo/registry \
    --mount=type=cache,target=/usr/local/cargo/git \
    --mount=type=cache,target=$SCCACHE_DIR,sharing=locked \
    cargo chef prepare --recipe-path recipe.json

# ── builder: compile deps then source ────────────────────────────────────────
FROM base AS builder
WORKDIR /app

# Cook dependencies — this layer is cache-hit whenever recipe.json is unchanged,
# i.e. whenever Cargo.toml / Cargo.lock are unmodified.
COPY --from=planner /app/recipe.json recipe.json
RUN --mount=type=cache,target=/usr/local/cargo/registry \
    --mount=type=cache,target=/usr/local/cargo/git \
    --mount=type=cache,target=$SCCACHE_DIR,sharing=locked \
    cargo chef cook --release --recipe-path recipe.json

# Compile the application (only application source, deps already baked in).
COPY . .
RUN --mount=type=cache,target=/usr/local/cargo/registry \
    --mount=type=cache,target=/usr/local/cargo/git \
    --mount=type=cache,target=$SCCACHE_DIR,sharing=locked \
    cargo build --release

# ── runtime: minimal final image ─────────────────────────────────────────────
FROM debian:bookworm-slim AS runtime

# libgfortran5: Fortran runtime linked by the OpenBLAS build.
# ca-certificates: needed for HTTPS downloads (LD score files, summary stats).
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        libgfortran5 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/target/release/ldsc /usr/local/bin/ldsc

ENTRYPOINT ["ldsc"]
