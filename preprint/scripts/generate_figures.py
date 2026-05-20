#!/usr/bin/env python3
"""Generate all figures for the ldsc-rs preprint."""

import gzip
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FIGDIR = Path(__file__).parent.parent / "figures"
FIGDIR.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "Helvetica", "Arial"],
    "font.size": 9,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "figure.dpi": 300,
})


def fig7_optimal_d():
    """3-panel: wall(d|N), accuracy(d|N), h2(d|N) — optimal-d characterization."""
    import csv
    csv_path = Path(__file__).parent.parent / "data" / "dn_sweep_full.csv"
    if not csv_path.exists():
        print("  SKIP fig7: dn_sweep_full.csv missing — run /tmp/dn_analyze.py first")
        return
    rows = []
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            rows.append(r)

    # Parse mode → d (mask suffix); exact-mask treated as d=∞ baseline
    def mode_to_d(m):
        if m == "exact-mask":
            return None
        # mode like "sk1600-mask" → 1600
        s = m.removeprefix("sk").removesuffix("-mask")
        return int(s)

    by_n = {}
    for r in rows:
        N = int(r["N"]); d = mode_to_d(r["mode"])
        by_n.setdefault(N, {})[d] = {
            "wall": float(r["wall_s"]),
            "rss": int(r["peak_rss_mb"]),
            "r": float(r["pearson_r"]) if r["pearson_r"] not in ("", "nan") else float("nan"),
            "h2": float(r["h2"]) if r["h2"] else float("nan"),
            "h2_se": float(r["h2_se"]) if r["h2_se"] else float("nan"),
        }

    Ns = sorted(by_n.keys())
    cmap = plt.get_cmap("viridis")
    n_colors = {N: cmap(i / max(1, len(Ns) - 1)) for i, N in enumerate(Ns)}

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4))

    # Panel A: wall vs d
    ax = axes[0]
    for N in Ns:
        cells = [(d, v["wall"]) for d, v in by_n[N].items()
                 if d is not None]
        if not cells:
            continue
        cells.sort()
        ds, ws = zip(*cells)
        ax.plot(ds, ws, marker="o", markersize=4, color=n_colors[N],
                linewidth=1.0, label=f"N = {N:,}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Sketch dimension d")
    ax.set_ylabel("Wall-clock time (s)")
    ax.set_title("Time vs d", fontsize=9)
    ax.legend(fontsize=6.5, loc="upper left", framealpha=0.9)
    ax.grid(True, which="major", alpha=0.3, linewidth=0.4)
    ax.grid(True, which="minor", alpha=0.15, linewidth=0.3)

    # Panel B: Pearson r vs d (sketch vs exact-mask baseline at same N)
    ax = axes[1]
    for N in Ns:
        cells = [(d, v["r"]) for d, v in by_n[N].items()
                 if d is not None and not np.isnan(v["r"])]
        if not cells:
            continue
        cells.sort()
        ds, rs = zip(*cells)
        ax.plot(ds, rs, marker="o", markersize=4, color=n_colors[N],
                linewidth=1.0, label=f"N = {N:,}")
    ax.set_xscale("log")
    ax.axhline(0.999, color="gray", linestyle=":", linewidth=0.5, alpha=0.7)
    ax.text(60, 0.9985, "r = 0.999", fontsize=6, color="gray")
    ax.set_xlabel("Sketch dimension d")
    ax.set_ylabel("Pearson r (sketch vs exact f32 + mask)")
    ax.set_title("Accuracy vs d", fontsize=9)
    ax.set_ylim(0.7, 1.005)
    ax.grid(True, which="major", alpha=0.3, linewidth=0.4)

    # Panel C: h2 vs d, normalized to exact-mask baseline at same N
    ax = axes[2]
    for N in Ns:
        if None not in by_n[N]:
            continue
        truth_h2 = by_n[N][None]["h2"]
        if np.isnan(truth_h2):
            continue
        cells = [(d, v["h2"] - truth_h2) for d, v in by_n[N].items()
                 if d is not None and not np.isnan(v["h2"])]
        if not cells:
            continue
        cells.sort()
        ds, dh = zip(*cells)
        ax.plot(ds, dh, marker="o", markersize=4, color=n_colors[N],
                linewidth=1.0, label=f"N = {N:,}")
    ax.set_xscale("log")
    ax.axhline(0, color="gray", linestyle="-", linewidth=0.6, alpha=0.6)
    ax.axhspan(-0.001, 0.001, color="green", alpha=0.1)
    ax.text(60, 0.0012, "truth cluster (±0.001)",
            fontsize=6, color="green")
    ax.set_xlabel("Sketch dimension d")
    ax.set_ylabel(r"$\hat{h}^2_{sketch} - \hat{h}^2_{exact}$ (BMI Yengo 2018)")
    ax.set_title("Heritability drift vs d", fontsize=9)
    ax.grid(True, which="major", alpha=0.3, linewidth=0.4)

    for a in axes:
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(FIGDIR / "fig7_optimal_d.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  fig7_optimal_d.png")


def fig6_maf_facet():
    """Per-SNP L2 error distribution faceted by MAF bin.

    Inputs (set via env or hardcoded /tmp paths):
      /tmp/gcta_1000G_eur.score.ld    (full-genome GCTA: SNP MAF + L2)
      /tmp/ldsc_1000G_eur*.l2.ldscore.gz (full-genome ldsc-rs exact, per-chr)
      /tmp/ldsc_1000G_eur_sketch200*.l2.ldscore.gz (optional sketch=200 per-chr)
    """
    import csv
    gcta_path = Path("/tmp/gcta_1000G_eur.score.ld")
    if not gcta_path.exists():
        print("  SKIP fig6: missing GCTA output")
        return

    # Read GCTA output: SNP chr bp MAF mean_rsq snp_num max_rsq ldscore
    gcta = {}
    maf = {}
    with open(gcta_path) as f:
        next(f)  # header
        for line in f:
            parts = line.split()
            if len(parts) < 8:
                continue
            snp = parts[0]
            gcta[snp] = float(parts[7])
            maf[snp] = float(parts[3])

    # Read all ldsc-rs per-chr exact outputs
    rust_exact = {}
    for chr_idx in range(1, 23):
        p = Path(f"/tmp/ldsc_1000G_eur{chr_idx}.l2.ldscore.gz")
        if not p.exists():
            continue
        with gzip.open(p, "rt") as f:
            next(f)
            for line in f:
                parts = line.strip().split("\t")
                rust_exact[parts[1]] = float(parts[3])

    if not rust_exact:
        print("  SKIP fig6: missing ldsc-rs full-genome L2 outputs")
        return

    # Optional sketch=200 outputs
    rust_sketch = {}
    for chr_idx in range(1, 23):
        p = Path(f"/tmp/ldsc_1000G_eur_sketch200_{chr_idx}.l2.ldscore.gz")
        if not p.exists():
            p = Path(f"/tmp/ldsc_1000G_eur_sketch200{chr_idx}.l2.ldscore.gz")
            if not p.exists():
                continue
        with gzip.open(p, "rt") as f:
            next(f)
            for line in f:
                parts = line.strip().split("\t")
                rust_sketch[parts[1]] = float(parts[3])

    snps_ge = sorted(set(gcta) & set(rust_exact))
    if not snps_ge:
        print("  SKIP fig6: no SNP overlap between GCTA + ldsc-rs")
        return

    diff_gcta_rust = np.array(
        [gcta[s] - rust_exact[s] for s in snps_ge]
    )
    maf_arr = np.array([maf[s] for s in snps_ge])

    bins = [(0.05, 0.10), (0.10, 0.20), (0.20, 0.35), (0.35, 0.50)]
    bin_labels = [f"[{a:.2f}, {b:.2f})" for a, b in bins]

    fig, axes = plt.subplots(1, 2 if rust_sketch else 1,
                             figsize=(7.0, 3.2), sharey=False)
    if not rust_sketch:
        axes = [axes]

    # Panel A: GCTA - ldsc-rs exact
    bin_data_a = []
    for lo, hi in bins:
        mask = (maf_arr >= lo) & (maf_arr < hi)
        bin_data_a.append(diff_gcta_rust[mask])
    bp = axes[0].boxplot(bin_data_a, labels=bin_labels, showfliers=False,
                         patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("#98df8a")
    axes[0].axhline(0, color="gray", linestyle=":", linewidth=0.6)
    axes[0].set_xlabel("MAF bin")
    axes[0].set_ylabel("L2 difference (GCTA − ldsc-rs exact)")
    axes[0].set_title(
        f"(a) GCTA vs ldsc-rs exact (n = {len(snps_ge):,} SNPs)",
        fontsize=8,
    )

    # Panel B: Rust exact - Rust sketch=200 (if available)
    if rust_sketch:
        snps_es = sorted(set(rust_exact) & set(rust_sketch) & set(maf))
        diff_es = np.array(
            [rust_exact[s] - rust_sketch[s] for s in snps_es]
        )
        maf_es = np.array([maf[s] for s in snps_es])
        bin_data_b = []
        for lo, hi in bins:
            mask = (maf_es >= lo) & (maf_es < hi)
            bin_data_b.append(diff_es[mask])
        bp = axes[1].boxplot(bin_data_b, labels=bin_labels, showfliers=False,
                             patch_artist=True)
        for patch in bp["boxes"]:
            patch.set_facecolor("#1f77b4")
        axes[1].axhline(0, color="gray", linestyle=":", linewidth=0.6)
        axes[1].set_xlabel("MAF bin")
        axes[1].set_ylabel("L2 difference (exact − sketch d=200)")
        axes[1].set_title(
            f"(b) Sketch (d=200) error vs MAF (n = {len(snps_es):,})",
            fontsize=8,
        )

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, axis="y", alpha=0.25, linewidth=0.4)

    fig.tight_layout()
    fig.savefig(FIGDIR / "fig6_maf_facet.png", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  fig6_maf_facet.png (n_GCTA_join={len(snps_ge):,})")


def fig5_scaling_curve():
    """T vs N (log-log) for 4 modes, with cost-model fit."""
    import csv
    rows = []
    csv_path = Path(__file__).parent.parent / "data" / "scaling_bench.csv"
    if not csv_path.exists():
        print("  SKIP fig5: missing scaling_bench.csv")
        return
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            rows.append((int(r["N"]), r["mode"], float(r["wall_s"]),
                         float(r["peak_rss_mb"])))

    modes = ["sketch-1600-mask", "sketch-200", "sketch-1000",
             "exact-f32", "exact-f64"]
    labels = {"sketch-1600-mask": "sketch-1600 + mask (truth cluster)",
              "sketch-200": "sketch-200",
              "sketch-1000": "sketch-1000",
              "exact-f32": "exact-f32", "exact-f64": "exact-f64"}
    colors = {"sketch-1600-mask": "#d62728",
              "sketch-200": "#aec7e8", "sketch-1000": "#9edae5",
              "exact-f32": "#98df8a", "exact-f64": "#2ca02c"}
    markers = {"sketch-1600-mask": "D",
               "sketch-200": "o", "sketch-1000": "s",
               "exact-f32": "^", "exact-f64": "v"}

    fig, (ax_t, ax_m) = plt.subplots(1, 2, figsize=(7.2, 3.4))

    # Time panel (left)
    for mode in modes:
        pts = [(n, w) for (n, m, w, _r) in rows if m == mode]
        if not pts:
            continue
        ns, ws = zip(*sorted(pts))
        ax_t.plot(ns, ws, color=colors[mode], marker=markers[mode],
                  markersize=5, linewidth=1.1, label=labels[mode])
    ax_t.set_xscale("log"); ax_t.set_yscale("log")
    ax_t.set_xlabel("Number of individuals (N)")
    ax_t.set_ylabel("Wall-clock time (seconds)")
    ax_t.set_title("Time vs N", fontsize=9)
    ax_t.legend(loc="upper left", fontsize=7, framealpha=0.9)
    ax_t.grid(True, which="major", alpha=0.3, linewidth=0.5)
    ax_t.grid(True, which="minor", alpha=0.15, linewidth=0.4)

    # Memory panel (right)
    for mode in modes:
        pts = [(n, r) for (n, m, _w, r) in rows if m == mode]
        if not pts:
            continue
        ns, rss = zip(*sorted(pts))
        ax_m.plot(ns, rss, color=colors[mode], marker=markers[mode],
                  markersize=5, linewidth=1.1, label=labels[mode])
    ax_m.set_xscale("log"); ax_m.set_yscale("log")
    ax_m.set_xlabel("Number of individuals (N)")
    ax_m.set_ylabel("Peak resident set size (MB)")
    ax_m.set_title("Memory vs N", fontsize=9)
    ax_m.grid(True, which="major", alpha=0.3, linewidth=0.5)
    ax_m.grid(True, which="minor", alpha=0.15, linewidth=0.4)

    for ax in (ax_t, ax_m):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(FIGDIR / "fig5_scaling_curve.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  fig5_scaling_curve.png")


def fig0_architecture():
    """Workflow / architecture diagram: data flow through ldsc-rs."""
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")

    def box(x, y, w, h, text, color="#e6f0ff", edge="#1f77b4", fontsize=7.5,
            bold=False):
        ax.add_patch(FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.04,rounding_size=0.12",
            linewidth=0.9, facecolor=color, edgecolor=edge,
        ))
        weight = "bold" if bold else "normal"
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=fontsize, fontweight=weight)

    def arrow(x0, y0, x1, y1, style="-|>", color="#333", lw=0.9, ls="-"):
        ax.add_patch(FancyArrowPatch(
            (x0, y0), (x1, y1),
            arrowstyle=style, mutation_scale=10, lw=lw, color=color,
            linestyle=ls, shrinkA=2, shrinkB=2,
        ))

    # Main vertical path (left column)
    box(0.3, 8.7, 4.0, 0.9,
        "PLINK BED / BIM / FAM\n(genotype + metadata)",
        color="#fff4e6", edge="#ff7f0e", bold=True)
    box(0.3, 7.4, 4.0, 0.9,
        "BED decode + normalize\n(256-entry LUT, AVX2 / FMA)",
        color="#e6f0ff", edge="#1f77b4")
    box(0.3, 6.1, 4.0, 0.9,
        "Ring-buffer chunked GEMM\n($B^\\top B$, $A^\\top B$ via faer)",
        color="#e6f0ff", edge="#1f77b4")
    box(0.3, 4.8, 4.0, 0.7,
        "Per-chromosome parallel pool (rayon)",
        color="#e6f0ff", edge="#1f77b4")
    box(0.3, 3.6, 4.0, 0.8,
        "L2 / M / M$_{5,50}$ per-chr outputs",
        color="#e6f7e6", edge="#2ca02c", bold=True)
    box(0.3, 1.9, 4.0, 1.2,
        "h² / r$_g$ / partitioned h²\n(IRWLS + 200-block jackknife)",
        color="#e6f7e6", edge="#2ca02c")

    # Main-path arrows
    for y0, y1 in [(8.7, 8.3), (7.4, 7.0), (6.1, 5.5), (4.8, 4.4), (3.6, 3.1)]:
        arrow(2.3, y0, 2.3, y1)

    # Right column: alt mode branches
    box(5.5, 7.4, 4.2, 0.9,
        "Alt: --sketch d\nFused CountSketch O(Nc), bypasses GEMM",
        color="#fff4f0", edge="#d62728")
    box(5.5, 6.1, 4.2, 0.9,
        "Alt: --snp-level-masking\nPost-GEMM r² mask, <1% overhead",
        color="#fff4f0", edge="#d62728")
    box(5.5, 4.8, 4.2, 0.7,
        "Alt: --python-compat (bit-identical)",
        color="#fff4f0", edge="#d62728")

    # Lateral arrows from main GEMM box to alt modes
    arrow(4.3, 7.85, 5.5, 7.85, style="-|>", color="#888", lw=0.7, ls="--")
    arrow(4.3, 6.55, 5.5, 6.55, style="-|>", color="#888", lw=0.7, ls="--")
    arrow(4.3, 5.15, 5.5, 5.15, style="-|>", color="#888", lw=0.7, ls="--")

    # Deployment row (bottom)
    box(0.3, 0.3, 4.0, 1.0,
        "Native CLI binary\n(musl, AVX2 / FMA, mimalloc)",
        color="#f0f0f0", edge="#444")
    box(5.5, 0.3, 4.2, 1.0,
        "ldsc-web (WASM)\nsharifhsn.github.io/ldsc/\nSAB workers, hand-rolled SIMD",
        color="#f0f0f0", edge="#444")
    arrow(2.3, 1.9, 2.3, 1.3)
    arrow(2.3, 1.9, 7.6, 1.3, style="-|>", color="#888", lw=0.7, ls=":")

    fig.tight_layout()
    fig.savefig(FIGDIR / "fig0_architecture.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  fig0_architecture.png")


def fig1_performance_1000g():
    """Bar chart: 1000G performance comparison."""
    labels = [
        "Python\nLDSC", "ldsc-rs\nexact f64", "ldsc-rs\nexact f32",
        "sketch\nd=50", "sketch\nd=200", "sketch\nd=500", "sketch\nd=1000",
    ]
    # AWS exact times; local sketch times (Ryzen 5 5600X)
    times = [1548.5, 41.1, 30.4, 3.7, 5.2, 8.4, 14.2]
    speedups = ["1.0x", "37.7x", "~51x", "~418x", "~298x", "~184x", "~109x"]
    colors = ["#d62728", "#2ca02c", "#98df8a", "#1f77b4", "#1f77b4", "#1f77b4", "#1f77b4"]

    fig, ax = plt.subplots(figsize=(6.5, 3.2))
    bars = ax.bar(range(len(labels)), times, color=colors, edgecolor="black", linewidth=0.6)

    ax.set_yscale("log")
    ax.set_ylabel("Wall-clock time (seconds)")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylim(1, 3000)

    for i, (t, label) in enumerate(zip(times, speedups)):
        ax.text(i, t * 1.2, label, ha="center", va="bottom", fontsize=6, fontweight="bold")

    ax.set_title("1000 Genomes Phase 3 (m=1.66M SNPs, N=2,490)", fontsize=9, pad=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig1_performance_1000g.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  fig1_performance_1000g.png")


def fig2_biobank_scaling():
    """Dual-axis plot: biobank time + accuracy vs sketch dimension, with cost model."""
    # Measured data (AWS EPYC 7R13, 16 vCPU, 26 data points)
    d_measured = np.array([25,50,75,100,150,200,300,400,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3500,4000,4500,5000,6000,7500,10000])
    t_measured = np.array([16,16,16,16,16,16,17,16,18,17,20,21,21,22,22,24,25,27,28,30,32,35,36,40,49,59])
    # Exact baselines
    t_exact_f64 = 656.0
    t_exact_f32 = 424.0

    # Cost model: T(d) = 14.41 + 4.45 * d/1000 (fit to d>=1000, R²=0.9955)
    T_SCATTER = 14.41
    T_GEMM_PER_1K = 4.45
    d_model = np.linspace(50, 12000, 200)
    t_model = T_SCATTER + T_GEMM_PER_1K * d_model / 1000

    # Accuracy model: r(d) ≈ 1 - 8/d (approximate, valid for d ≤ 1000)
    C_ACC = 8.0
    # Measured on full biobank (N=50K, 1.66M SNPs) — 26 data points
    d_acc = np.array([25,50,75,100,150,200,300,400,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3500,4000,4500,5000,6000,7500,10000])
    r_acc = np.array([0.724058,0.842493,0.899158,0.922367,0.945361,0.959743,0.972237,0.977193,0.981231,0.987793,0.990241,0.991484,0.992600,0.993300,0.993785,0.993932,0.994661,0.994977,0.995019,0.995506,0.995731,0.995890,0.996104,0.996431,0.996649,0.996979])
    r_model = 1 - C_ACC / d_model

    # Crossover point
    d_star = T_SCATTER / (T_GEMM_PER_1K / 1000)  # ~3238

    fig, ax1 = plt.subplots(figsize=(6, 3.8))
    color_time = "#2ca02c"
    color_acc = "#1f77b4"
    color_model = "#d62728"

    # Time axis (left)
    ax1.plot(d_model, t_model, "--", color=color_model, linewidth=1.2, alpha=0.7,
             label=f"Model: T = {T_SCATTER:.1f} + {T_GEMM_PER_1K:.1f}d/1000")
    ax1.plot(d_measured, t_measured, "o", color=color_time, markersize=6, zorder=4,
             label="Measured time")
    ax1.axhline(y=T_SCATTER, color="gray", linestyle=":", linewidth=0.6, alpha=0.5)
    ax1.text(11500, T_SCATTER + 2, f"scatter-add floor ({T_SCATTER:.1f}s)",
             fontsize=6, color="gray", ha="right")
    ax1.axvline(x=d_star, color="gray", linestyle="--", linewidth=0.7, alpha=0.4)
    ax1.text(d_star + 100, 58, f"$d^* \\approx {d_star:.0f}$\n(crossover)",
             fontsize=6, color="gray")
    ax1.set_xlabel("Sketch dimension (d)")
    ax1.set_ylabel("Wall-clock time (seconds)", color=color_time)
    ax1.tick_params(axis="y", labelcolor=color_time)
    ax1.set_ylim(0, 75)
    ax1.set_xlim(0, 12000)

    # Accuracy axis (right)
    ax2 = ax1.twinx()
    ax2.plot(d_acc, r_acc, "s", color=color_acc, markersize=5, zorder=3,
             label="Measured accuracy")
    ax2.plot(d_model, r_model, "-", color=color_acc, linewidth=1, alpha=0.5,
             label=f"Model: r = 1 - {C_ACC:.0f}/d")
    ax2.set_ylabel("Pearson r with exact LD scores", color=color_acc)
    ax2.tick_params(axis="y", labelcolor=color_acc)
    ax2.set_ylim(0.75, 1.005)

    ax1.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="center right", fontsize=6.5,
               framealpha=0.9)

    ax1.set_title("Biobank scale (m=1.66M SNPs, N=50,000): cost model", fontsize=9, pad=10)
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig2_biobank_scaling.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  fig2_biobank_scaling.png")


def _load_ldscore(path):
    """Load LD scores from a gzipped file, deduplicating by SNP ID."""
    scores = {}
    with gzip.open(path, "rt") as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            snp = parts[1]
            if snp not in scores:  # keep first occurrence (dedup)
                scores[snp] = float(parts[3])
    return scores


def fig3_parity_scatter(python_path, rust_exact_path, sketch_path):
    """Scatter plot: (a) Python vs Rust exact, (b) exact vs sketch."""
    print("  Loading LD scores for scatter plot...")
    py_scores = _load_ldscore(python_path)
    rust_scores = _load_ldscore(rust_exact_path)
    sketch_scores = _load_ldscore(sketch_path)

    # Panel A: Python vs Rust exact
    common_a = sorted(set(py_scores) & set(rust_scores))
    n_a = len(common_a)
    rng = np.random.RandomState(42)
    idx_a = rng.choice(n_a, min(50000, n_a), replace=False)
    sub_a = [common_a[i] for i in sorted(idx_a)]
    py_vals = np.array([py_scores[s] for s in sub_a])
    rust_vals = np.array([rust_scores[s] for s in sub_a])
    mad_a = np.max(np.abs(py_vals - rust_vals))
    r_a = np.corrcoef(py_vals, rust_vals)[0, 1]

    # Panel B: Rust exact vs sketch
    common_b = sorted(set(rust_scores) & set(sketch_scores))
    n_b = len(common_b)
    idx_b = rng.choice(n_b, min(50000, n_b), replace=False)
    sub_b = [common_b[i] for i in sorted(idx_b)]
    exact_vals = np.array([rust_scores[s] for s in sub_b])
    sketch_vals = np.array([sketch_scores[s] for s in sub_b])
    r_b = np.corrcoef(exact_vals, sketch_vals)[0, 1]

    fig, axes = plt.subplots(1, 2, figsize=(7, 3.2))

    # Panel A: real Python vs Rust exact
    ax = axes[0]
    ax.scatter(py_vals, rust_vals, s=0.5, alpha=0.3, color="#2ca02c", rasterized=True)
    lims = [min(py_vals.min(), rust_vals.min()), min(max(py_vals.max(), rust_vals.max()), 400)]
    ax.plot(lims, lims, "k--", linewidth=0.5, alpha=0.5)
    ax.set_xlabel("Python LDSC LD scores")
    ax.set_ylabel("ldsc-rs exact LD scores")
    ax.set_title(f"(a) Exact mode: max_abs_diff = {mad_a:.1f}", fontsize=8)
    ax.text(0.05, 0.92, f"n = {n_a:,} SNPs", transform=ax.transAxes, fontsize=7)
    ax.text(0.05, 0.84, f"r = {r_a:.6f}", transform=ax.transAxes, fontsize=7)
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    # Panel B: exact vs sketch
    ax = axes[1]
    ax.scatter(exact_vals, sketch_vals, s=0.5, alpha=0.3, color="#1f77b4", rasterized=True)
    lims = [min(exact_vals.min(), sketch_vals.min()), min(max(exact_vals.max(), sketch_vals.max()), 400)]
    ax.plot(lims, lims, "k--", linewidth=0.5, alpha=0.5)
    ax.set_xlabel("Exact LD scores")
    ax.set_ylabel("CountSketch (d=200) LD scores")
    ax.set_title(f"(b) Sketch mode: r = {r_b:.3f}", fontsize=8)
    ax.text(0.05, 0.92, f"n = {len(sub_b):,} SNPs (subsampled)", transform=ax.transAxes, fontsize=7)
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(FIGDIR / "fig3_parity_scatter.png", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  fig3_parity_scatter.png (panel a: mad={mad_a:.4f}, r={r_a:.6f}, n={n_a}; panel b: r={r_b:.4f}, n={n_b})")


def fig4_sketch_accuracy():
    """Accuracy vs sketch dimension — both datasets + universal 1-C/d model."""
    # 1000G measured (N=2490) — 16 data points
    d_1k = np.array([25,50,75,100,150,200,300,400,500,750,1000,1250,1500,1750,2000,2250])
    r_1k = np.array([0.674082,0.822052,0.894936,0.911202,0.943675,0.945535,0.966553,0.972597,0.981291,0.984958,0.992505,0.993045,0.994060,0.995594,0.996160,0.996723])

    # Biobank (N=50K) — measured on full 1.66M SNPs, 26 data points
    d_bio = np.array([25,50,75,100,150,200,300,400,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3500,4000,4500,5000,6000,7500,10000])
    r_bio = np.array([0.724058,0.842493,0.899158,0.922367,0.945361,0.959743,0.972237,0.977193,0.981231,0.987793,0.990241,0.991484,0.992600,0.993300,0.993785,0.993932,0.994661,0.994977,0.995019,0.995506,0.995731,0.995890,0.996104,0.996431,0.996649,0.996979])

    # Approximate accuracy model: r ≈ 1 - C/d, C ≈ 8 (valid for d ≤ 1000)
    C = 8.0
    d_model = np.linspace(30, 12000, 300)
    r_model = 1 - C / d_model

    # Cost model for biobank time (right axis)
    T_SCATTER = 12.79
    T_GEMM_PER_1K = 5.026
    t_model = T_SCATTER + T_GEMM_PER_1K * d_model / 1000

    fig, ax1 = plt.subplots(figsize=(6, 3.8))

    color_1k = "#1f77b4"
    color_bio = "#ff7f0e"
    color_model = "gray"
    color_time = "#d62728"

    # Accuracy data
    ax1.plot(d_1k, r_1k, "o", color=color_1k, markersize=5, zorder=3,
             label="1000G (N=2,490)")
    ax1.plot(d_bio, r_bio, "s", color=color_bio, markersize=5, zorder=3,
             label="Biobank (N=50,000)")
    # Universal model
    ax1.plot(d_model, r_model, "-", color=color_model, linewidth=1.2, alpha=0.6,
             label=f"Model: r = 1 − {C:.0f}/d", zorder=2)

    ax1.set_xlabel("Sketch dimension (d)")
    ax1.set_ylabel("Pearson r with exact LD scores")
    ax1.set_ylim(0.75, 1.005)
    ax1.set_xlim(0, 11000)
    ax1.axhline(y=1.0, color="gray", linestyle=":", linewidth=0.5, alpha=0.3)

    # Annotate key accuracy thresholds
    for r_target, d_opt, label in [(0.99, 900, "r=0.99\nd=900"), (0.995, 1800, "r=0.995\nd=1800")]:
        ax1.annotate(label, xy=(d_opt, r_target), xytext=(d_opt + 800, r_target - 0.04),
                     arrowprops=dict(arrowstyle="->", color="gray", lw=0.6),
                     fontsize=6, color="gray")

    # Crossover annotation
    d_star = 3238
    ax1.axvline(x=d_star, color="gray", linestyle="--", linewidth=0.6, alpha=0.3)
    ax1.text(d_star + 100, 0.77, f"Biobank $d^*$≈{d_star}", fontsize=6, color="gray")

    # Time on secondary axis
    ax2 = ax1.twinx()
    ax2.plot(d_model, t_model, "--", color=color_time, linewidth=1, alpha=0.5,
             label="Biobank time (model)")
    d_meas = np.array([25,50,75,100,150,200,300,400,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3500,4000,4500,5000,6000,7500,10000])
    t_meas = np.array([16,16,16,16,16,16,17,16,18,17,20,21,21,22,22,24,25,27,28,30,32,35,36,40,49,59])
    ax2.plot(d_meas, t_meas, "^", color=color_time, markersize=5, alpha=0.7,
             label="Biobank time (measured)")
    ax2.set_ylabel("Wall-clock time, N=50K (seconds)", color=color_time, fontsize=8)
    ax2.tick_params(axis="y", labelcolor=color_time)
    ax2.set_ylim(0, 75)

    ax1.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="center right", fontsize=6.5,
               framealpha=0.9)

    ax1.set_title("CountSketch: accuracy and cost vs. sketch dimension (42 points)", fontsize=9, pad=10)
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig4_sketch_accuracy.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  fig4_sketch_accuracy.png")


if __name__ == "__main__":
    print("Generating figures...")
    fig0_architecture()
    fig1_performance_1000g()
    fig2_biobank_scaling()
    fig4_sketch_accuracy()
    fig5_scaling_curve()
    fig6_maf_facet()
    fig7_optimal_d()

    # Fig 3 needs Python + Rust exact + sketch LD score files
    python_path = Path("/tmp/py_ld.l2.ldscore.gz")
    rust_exact_path = Path("/tmp/rust_ld.l2.ldscore.gz")
    sketch_path = Path("/tmp/rust_ld_sketch200.l2.ldscore.gz")
    if python_path.exists() and rust_exact_path.exists() and sketch_path.exists():
        fig3_parity_scatter(str(python_path), str(rust_exact_path), str(sketch_path))
    else:
        missing = [str(p) for p in [python_path, rust_exact_path, sketch_path] if not p.exists()]
        print(f"  SKIP fig3: missing {', '.join(missing)}")

    print("Done. Figures in", FIGDIR)
