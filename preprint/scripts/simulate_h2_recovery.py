#!/usr/bin/env python3
"""
Simulation study: does --snp-level-masking reduce h² bias?

Recovers true h² from synthetic phenotypes using LD scores computed both with
and without --snp-level-masking, on chromosome 22 of 1000G.

If masked h² is closer to truth across replicates, the masking flag is a
systematic correction, not just a difference. Closes preprint review_notes
item #5.

Usage:
    python preprint/scripts/simulate_h2_recovery.py prep
    python preprint/scripts/simulate_h2_recovery.py simulate
    python preprint/scripts/simulate_h2_recovery.py aggregate
    python preprint/scripts/simulate_h2_recovery.py all  # runs all three
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parents[2]
BFILE = ROOT / "data" / "1000G_eur"   # 503 EUR samples, matches paper's homogeneous cohort
LDSC_BIN = ROOT / "target" / "release" / "ldsc"
PLINK2_BIN = Path("/tmp/plink2/plink2")  # canonical GWAS tool

OUT_DIR = ROOT / "preprint" / "data"
LDSCORE_DIR = OUT_DIR / "sim_ldscores"
SUMSTATS_DIR = OUT_DIR / "sim_sumstats"
PLINK_DIR = OUT_DIR / "sim_plink"
RESULTS_CSV = OUT_DIR / "h2_simulation_results.csv"
SNPLIST = OUT_DIR / "chr22_maf05.snplist"

CHR_TARGET = "22"
MAF_MIN = 0.05
LD_WIND_KB = 1000          # canonical LDSC window (≈ 1 cM)
LD_WIND_KB_TINY = 1        # control condition for sanity check
H2_TRUE_VALUES = [0.2, 0.5]
N_REPLICATES = 50
H2_RECOVERY_TOL_SE = 2.0   # mean(h2_masked) within this many SE of h2_true

# h2 stdout regexes (format strings live in src/regressions.rs:745,747,749,751)
RE_H2 = re.compile(r"Total Observed scale h2: ([-0-9.]+) \(([-0-9.]+)\)")
RE_INTERCEPT_FREE = re.compile(r"Intercept: ([-0-9.]+) \(([-0-9.]+)\)")
RE_INTERCEPT_FIXED = re.compile(r"Intercept: constrained to ([-0-9.]+)")
RE_MEAN_CHI2 = re.compile(r"Mean Chi\^2: ([-0-9.]+)")


# ---------------------------------------------------------------------------
# Helpers — BED I/O, simulation, regex parsers
# ---------------------------------------------------------------------------

def read_bim(prefix: Path) -> list[tuple[str, str, int, str, str, str]]:
    """Return list of (chr, snp, cm, bp, a1, a2)."""
    rows = []
    with open(f"{prefix}.bim") as f:
        for line in f:
            parts = line.rstrip("\n").split()
            rows.append((parts[0], parts[1], parts[2], int(parts[3]), parts[4], parts[5]))
    return rows


def count_fam(prefix: Path) -> int:
    with open(f"{prefix}.fam") as f:
        return sum(1 for _ in f)


def select_chr_snp_indices(bim: list, chrom: str) -> list[int]:
    return [i for i, row in enumerate(bim) if row[0] == chrom]


def load_bed_subset_standardized(
    prefix: Path,
    snp_indices: list[int],
    n_indiv: int,
    maf_min: float,
) -> tuple[list[int], np.ndarray]:
    """
    Read selected SNPs from a SNP-major BED, decode 2-bit genotypes,
    apply MAF filter, standardize columns to mean 0 / var 1, NaN→0.

    Returns (kept_snp_indices, G[n_indiv, n_kept] float32).

    Lifted from scripts/make_synthetic_biobank.py:81-117 BED unpack pattern.
    """
    bytes_per_snp = (n_indiv + 3) // 4

    with open(f"{prefix}.bed", "rb") as f:
        magic = f.read(3)
        assert magic == b"\x6c\x1b\x01", "Invalid BED magic"
        # Read only the SNPs we need by seeking, in SNP-index order
        # snp_indices is already chromosome-sorted (BIM order)
        chunks = []
        for snp_i in snp_indices:
            f.seek(3 + snp_i * bytes_per_snp)
            chunks.append(f.read(bytes_per_snp))
        raw = np.frombuffer(b"".join(chunks), dtype=np.uint8).reshape(
            len(snp_indices), bytes_per_snp
        )

    # Unpack each byte to 4 2-bit codes: 00,01,10,11 → genotype slot
    expanded = np.zeros((len(snp_indices), bytes_per_snp * 4), dtype=np.uint8)
    for bp in range(4):
        expanded[:, bp::4] = (raw >> (bp * 2)) & 0x03
    raw_codes = expanded[:, :n_indiv]  # (M, N)

    # PLINK 1.9 BED encoding:
    #   00 → hom A1 (dosage 2)
    #   01 → missing (NaN)
    #   10 → het (dosage 1)
    #   11 → hom A2 (dosage 0)
    G = np.full(raw_codes.shape, np.nan, dtype=np.float32)
    G[raw_codes == 0x00] = 2.0
    G[raw_codes == 0x02] = 1.0
    G[raw_codes == 0x03] = 0.0
    # 0x01 stays NaN

    # MAF filter: based on observed alt-allele freq (treat NaN as missing)
    # MAF = min(p, 1-p) where p = mean(G)/2 over observed
    n_obs = np.sum(~np.isnan(G), axis=1)
    sums = np.nansum(G, axis=1)
    p = np.where(n_obs > 0, sums / (2.0 * n_obs), 0.0)
    maf = np.minimum(p, 1.0 - p)
    keep = maf >= maf_min

    G = G[keep]                                      # (M_kept, N)
    kept_snp_indices = [snp_indices[i] for i in range(len(snp_indices)) if keep[i]]

    # Standardize each row (SNP) to mean 0, variance 1, impute NaN→0 (= mean)
    # Note: returns G as (N, M) for downstream (individuals × SNPs)
    means = np.nanmean(G, axis=1, keepdims=True)
    G = np.where(np.isnan(G), means, G)
    centered = G - means
    var = np.mean(centered * centered, axis=1, keepdims=True)
    std = np.sqrt(np.where(var > 0, var, 1.0)).astype(np.float32)
    G_std = (centered / std).astype(np.float32)

    return kept_snp_indices, G_std.T  # (N, M_kept)


def simulate_phenotype(G: np.ndarray, h2: float, seed: int) -> np.ndarray:
    """Infinitesimal model: y = G·β + ε with var(y) ≈ 1, true heritability = h².

    G is standardized (N, M). β ~ N(0, h²/M). ε ~ N(0, 1-h²).
    Returns y of length N with var(y) ≈ 1.
    """
    n, m = G.shape
    rng = np.random.default_rng(seed)
    beta = rng.standard_normal(m).astype(np.float32) * np.sqrt(h2 / m)
    eps = rng.standard_normal(n).astype(np.float32) * np.sqrt(1.0 - h2)
    y = G @ beta + eps
    return y


def read_fam_iids(prefix: Path) -> list[tuple[str, str]]:
    """Return list of (FID, IID) tuples from .fam file."""
    rows = []
    with open(f"{prefix}.fam") as f:
        for line in f:
            parts = line.rstrip("\n").split()
            rows.append((parts[0], parts[1]))
    return rows


def write_pheno_batch(
    fid_iid: list[tuple[str, str]],
    pheno_names: list[str],
    pheno_matrix: np.ndarray,  # shape (n_samples, n_phenos)
    path: Path,
) -> None:
    """Write a plink2 batched .pheno file: #FID\tIID\tPHENO1\tPHENO2\t...

    All phenotype columns are written in one file so plink2 --glm processes them
    in a single invocation (much faster than 100 separate calls).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("#FID\tIID\t" + "\t".join(pheno_names) + "\n")
        for i, (fid, iid) in enumerate(fid_iid):
            vals = "\t".join(f"{v:.6f}" for v in pheno_matrix[i])
            f.write(f"{fid}\t{iid}\t{vals}\n")


def run_plink_glm_batch(
    bfile: Path,
    pheno_path: Path,
    snplist: Path,
    maf: float,
    out_prefix: Path,
) -> None:
    """Run plink2 --glm --linear on all phenotypes in pheno_path. Writes
    {out_prefix}.{pheno_name}.glm.linear for each phenotype column."""
    args = [
        str(PLINK2_BIN),
        "--bfile", str(bfile),
        "--extract", str(snplist),
        "--maf", str(maf),
        "--pheno", str(pheno_path),
        "--glm", "allow-no-covars", "hide-covar",
        "--out", str(out_prefix),
    ]
    proc = subprocess.run(args, check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"plink2 --glm failed (exit {proc.returncode})\n"
            f"args: {args}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )


def plink_glm_to_munge_input(plink_path: Path, tsv_path: Path) -> int:
    """Convert plink2 .glm.linear → munge-sumstats input TSV.

    plink2 cols: #CHROM POS ID REF ALT PROVISIONAL_REF? A1 OMITTED A1_FREQ TEST OBS_CT BETA SE T_STAT P [ERRCODE]
    munge cols:  SNP A1 A2 N BETA P
       (A2 = OMITTED, the non-tested allele; munge derives Z from P with sign from BETA)

    Returns the number of valid (non-NA) rows written.
    """
    n_written = 0
    with open(plink_path) as fin, open(tsv_path, "w") as fout:
        header = fin.readline().rstrip("\n").lstrip("#").split("\t")
        ix = {name: i for i, name in enumerate(header)}
        fout.write("SNP\tA1\tA2\tN\tBETA\tP\n")
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if parts[ix["P"]] == "NA" or parts[ix["BETA"]] == "NA":
                continue
            fout.write(
                f"{parts[ix['ID']]}\t{parts[ix['A1']]}\t{parts[ix['OMITTED']]}\t"
                f"{parts[ix['OBS_CT']]}\t{parts[ix['BETA']]}\t{parts[ix['P']]}\n"
            )
            n_written += 1
    return n_written


def run_ldsc_munge(input_tsv: Path, n: int, out_prefix: Path) -> None:
    """Run `ldsc munge-sumstats` to canonically convert P + sign(BETA) → Z via
    statrs::erfc_inv (src/munge.rs:481). Output is {out_prefix}.sumstats.gz."""
    args = [
        str(LDSC_BIN), "munge-sumstats",
        "--sumstats", str(input_tsv),
        "--N", str(n),
        "--signed-sumstats", "BETA,0",
        "--no-alleles",   # don't enforce A1/A2 strand checks (we trust plink output)
        "--out", str(out_prefix),
    ]
    proc = subprocess.run(args, check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"ldsc munge-sumstats failed (exit {proc.returncode})\n"
            f"args: {args}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )


def parse_h2_stdout(text: str) -> dict | None:
    """Parse `ldsc h2` stdout. Returns None on parse failure."""
    out: dict = {}
    m = RE_H2.search(text)
    if not m:
        return None
    out["h2"] = float(m.group(1))
    out["h2_se"] = float(m.group(2))
    m = RE_INTERCEPT_FREE.search(text)
    if m:
        out["intercept"] = float(m.group(1))
        out["intercept_se"] = float(m.group(2))
    else:
        m = RE_INTERCEPT_FIXED.search(text)
        if m:
            out["intercept"] = float(m.group(1))
            out["intercept_se"] = float("nan")
    m = RE_MEAN_CHI2.search(text)
    if m:
        out["mean_chi2"] = float(m.group(1))
    return out


def run_ldsc(args: list[str]) -> str:
    """Run an ldsc subprocess, return stdout. Raises on non-zero exit."""
    proc = subprocess.run(
        [str(LDSC_BIN), *args],
        check=False,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"ldsc failed (exit {proc.returncode})\nargs: {args}\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    return proc.stdout


# ---------------------------------------------------------------------------
# Subcommand: prep — subset chr22 + compute LD scores both ways
# ---------------------------------------------------------------------------

def cmd_prep(args) -> None:
    if not LDSC_BIN.exists():
        sys.exit(f"Missing {LDSC_BIN}. Run: cargo build --release")
    if not (BFILE.with_suffix(".bed")).exists():
        sys.exit(f"Missing {BFILE}.bed. Stage data from S3 first.")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LDSCORE_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Read BIM, select chr22
    print(f"[prep] Reading BIM from {BFILE}.bim ...")
    t0 = time.time()
    bim = read_bim(BFILE)
    n_indiv = count_fam(BFILE)
    chr22_idx = select_chr_snp_indices(bim, CHR_TARGET)
    print(f"[prep]   total SNPs={len(bim):,}  chr22 SNPs={len(chr22_idx):,}  N={n_indiv}  ({time.time()-t0:.1f}s)")

    # 2. Decode chr22 BED + apply MAF filter (Python side, used by simulate)
    print(f"[prep] Decoding chr22 BED + applying MAF≥{MAF_MIN} ...")
    t0 = time.time()
    kept_idx, _ = load_bed_subset_standardized(BFILE, chr22_idx, n_indiv, MAF_MIN)
    kept_snp_ids = [bim[i][1] for i in kept_idx]
    print(f"[prep]   kept SNPs after MAF filter: {len(kept_snp_ids):,}  ({time.time()-t0:.1f}s)")

    # 3. Write SNP list for ldsc l2 --extract
    SNPLIST.write_text("\n".join(kept_snp_ids) + "\n")
    print(f"[prep] Wrote {SNPLIST}")

    # 4. Run ldsc l2 twice (chunked + masked) at canonical 1000kb window
    for variant in ("chunked", "masked"):
        out_prefix = LDSCORE_DIR / variant
        print(f"[prep] Computing LD scores: {variant} → {out_prefix}.l2.ldscore.gz")
        cli = [
            "l2",
            "--bfile", str(BFILE),
            "--extract", str(SNPLIST),
            "--ld-wind-kb", str(LD_WIND_KB),
            "--maf", str(MAF_MIN),
            "--out", str(out_prefix),
            "--yes-really",
        ]
        if variant == "masked":
            cli.append("--snp-level-masking")
        t0 = time.time()
        run_ldsc(cli)
        print(f"[prep]   done ({time.time()-t0:.1f}s)")

    # 5. Also compute the tiny-window control LD scores for sanity check
    for variant in ("chunked", "masked"):
        out_prefix = LDSCORE_DIR / f"control_{variant}"
        print(f"[prep] Computing control LD scores (window={LD_WIND_KB_TINY}kb): {variant}")
        cli = [
            "l2",
            "--bfile", str(BFILE),
            "--extract", str(SNPLIST),
            "--ld-wind-kb", str(LD_WIND_KB_TINY),
            "--maf", str(MAF_MIN),
            "--out", str(out_prefix),
            "--yes-really",
        ]
        if variant == "masked":
            cli.append("--snp-level-masking")
        t0 = time.time()
        run_ldsc(cli)
        print(f"[prep]   done ({time.time()-t0:.1f}s)")

    print("[prep] Complete.")


# ---------------------------------------------------------------------------
# Subcommand: simulate — generate sumstats and run h² regressions
# ---------------------------------------------------------------------------

def cmd_simulate(args) -> None:
    if not LDSC_BIN.exists():
        sys.exit(f"Missing {LDSC_BIN}. Run: cargo build --release")
    if not PLINK2_BIN.exists():
        sys.exit(f"Missing {PLINK2_BIN}. Install plink2 alpha7 macOS arm64.")
    if not SNPLIST.exists():
        sys.exit("Missing SNP list. Run `prep` first.")

    SUMSTATS_DIR.mkdir(parents=True, exist_ok=True)
    PLINK_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Decode chr22 + standardize G ONCE (only used for phenotype simulation)
    print(f"[simulate] Loading + standardizing chr22 genotypes ...")
    t0 = time.time()
    bim = read_bim(BFILE)
    chr22_idx = select_chr_snp_indices(bim, CHR_TARGET)
    n_indiv = count_fam(BFILE)
    _, G = load_bed_subset_standardized(BFILE, chr22_idx, n_indiv, MAF_MIN)
    n, m = G.shape
    print(f"[simulate]   G.shape=({n}, {m})  RAM={G.nbytes/1e6:.0f} MB  ({time.time()-t0:.1f}s)")

    # 2. Generate all phenotypes (simulate y = Gβ + ε), one column per (h2, rep)
    print(f"[simulate] Generating {len(H2_TRUE_VALUES) * N_REPLICATES} phenotypes ...")
    t0 = time.time()
    pheno_names: list[str] = []
    pheno_meta: list[tuple[int, float, float, str]] = []  # (rep, h2_true, var_y, name)
    pheno_cols: list[np.ndarray] = []
    for h2_true in H2_TRUE_VALUES:
        for rep in range(N_REPLICATES):
            seed = 1000 * int(round(h2_true * 100)) + rep
            y = simulate_phenotype(G, h2_true, seed=seed)
            name = f"sim_h2_{int(round(h2_true*100)):03d}_rep{rep:03d}"
            pheno_names.append(name)
            pheno_cols.append(y)
            pheno_meta.append((rep, h2_true, float(y.var()), name))
    pheno_matrix = np.stack(pheno_cols, axis=1)  # (n_indiv, n_phenos)
    del pheno_cols
    print(f"[simulate]   done in {time.time()-t0:.1f}s, shape={pheno_matrix.shape}")

    # Sanity check #1: phenotype variance must be in [0.5, 1.5] for every rep.
    # SD of var(y) on chr22 EUR is ~0.07 due to LD-induced quadratic-form variance.
    bad_var = [(r, h, v) for (r, h, v, _) in pheno_meta if not (0.5 <= v <= 1.5)]
    if bad_var:
        sys.exit(f"[FAIL sanity #1: phenotype variance] {len(bad_var)} reps out of range. "
                 f"First: rep={bad_var[0][0]} h2={bad_var[0][1]} var={bad_var[0][2]:.3f}")
    var_arr = np.array([v for (_, _, v, _) in pheno_meta])
    print(f"[simulate] sanity #1 OK: var(y) range [{var_arr.min():.3f}, {var_arr.max():.3f}], "
          f"mean={var_arr.mean():.3f}, sd={var_arr.std():.3f}")

    # 3. Write batched .pheno file and run plink2 --glm ONCE for all phenotypes
    fid_iid = read_fam_iids(BFILE)
    pheno_path = PLINK_DIR / "all_phenos.pheno"
    write_pheno_batch(fid_iid, pheno_names, pheno_matrix, pheno_path)
    plink_out_prefix = PLINK_DIR / "glm"
    print(f"[simulate] Running plink2 --glm on {len(pheno_names)} phenotypes (single invocation) ...")
    t0 = time.time()
    run_plink_glm_batch(BFILE, pheno_path, SNPLIST, MAF_MIN, plink_out_prefix)
    print(f"[simulate]   plink2 done in {time.time()-t0:.1f}s")

    # 4. Convert each plink output → munge input TSV → ldsc munge-sumstats → .sumstats.gz
    print(f"[simulate] Converting plink outputs through ldsc munge-sumstats (statrs P→Z) ...")
    t0 = time.time()
    sumstats_paths: dict[str, str] = {}  # pheno_name → sumstats.gz path
    for i, name in enumerate(pheno_names):
        plink_out = PLINK_DIR / f"glm.{name}.glm.linear"
        if not plink_out.exists():
            sys.exit(f"[FAIL] missing plink output for {name}: {plink_out}")
        munge_input = PLINK_DIR / f"{name}.tsv"
        n_snps_w = plink_glm_to_munge_input(plink_out, munge_input)
        if n_snps_w == 0:
            sys.exit(f"[FAIL] no valid SNPs in plink output for {name}")
        munge_out_prefix = SUMSTATS_DIR / name
        run_ldsc_munge(munge_input, n_indiv, munge_out_prefix)
        sumstats_paths[name] = f"{munge_out_prefix}.sumstats.gz"
        if (i + 1) % 20 == 0:
            print(f"[simulate]   munged {i+1}/{len(pheno_names)} ({time.time()-t0:.0f}s)")
    print(f"[simulate]   munge complete in {time.time()-t0:.0f}s")

    # 5. Run ldsc h² for each (sumstats × LD score variant), sequentially.
    rows = []
    chunked_ref = LDSCORE_DIR / "chunked"
    masked_ref = LDSCORE_DIR / "masked"
    M_total = read_M_5_50(chunked_ref)

    n_runs = len(pheno_meta) * 2
    print(f"[simulate] Running {n_runs} h² regressions ...")
    t0 = time.time()
    for i, (rep, h2_true, var_y, name) in enumerate(pheno_meta):
        sumstats_path = sumstats_paths[name]
        for variant, ref in (("chunked", chunked_ref), ("masked", masked_ref)):
            stdout = run_ldsc([
                "h2",
                "--h2", sumstats_path,
                "--ref-ld", str(ref) + ".l2.ldscore.gz",
                "--w-ld", str(ref) + ".l2.ldscore.gz",
                "--M", str(M_total),
                "--out", str(SUMSTATS_DIR / f".tmp_{variant}_{name}"),
            ])
            parsed = parse_h2_stdout(stdout)
            if parsed is None:
                sys.exit(f"[FAIL] could not parse h² stdout for {name} {variant}\n{stdout}")
            rows.append({
                "rep": rep,
                "h2_true": h2_true,
                "masking": variant,
                "h2_est": parsed["h2"],
                "h2_se": parsed["h2_se"],
                "intercept": parsed.get("intercept", float("nan")),
                "intercept_se": parsed.get("intercept_se", float("nan")),
                "mean_chi2": parsed.get("mean_chi2", float("nan")),
                "var_y": var_y,
            })
        if (i + 1) % 10 == 0:
            elapsed = time.time() - t0
            done = (i + 1) * 2
            eta = elapsed * (n_runs - done) / done
            print(f"[simulate]   {done}/{n_runs} runs  elapsed={elapsed:.0f}s  ETA={eta:.0f}s")
    print(f"[simulate]   {n_runs} h² runs done in {time.time()-t0:.0f}s")

    # Sanity check #2: mean χ² ≈ 1 + N·h²·ℓ̄/M (relaxed: heavy-tailed across reps,
    # use median which is more robust)
    print("[simulate] Sanity check #2: median χ² vs predicted (LDSC theory) ...")
    ldbar_chunked = mean_ld_score(chunked_ref)
    for h2_true in H2_TRUE_VALUES:
        chi2_vals = [r["mean_chi2"] for r in rows
                     if r["h2_true"] == h2_true and r["masking"] == "chunked"]
        observed_mean = float(np.mean(chi2_vals))
        observed_median = float(np.median(chi2_vals))
        predicted = 1.0 + n_indiv * h2_true * ldbar_chunked / M_total
        # Use median (heavy tails) — relaxed check, just looking for "in the ballpark"
        rel_err = abs(observed_median - predicted) / predicted
        status = "OK " if rel_err < 0.30 else "WARN"
        print(f"[simulate]   h2={h2_true}  median χ²={observed_median:.3f}  mean χ²={observed_mean:.3f}  "
              f"predicted={predicted:.3f}  median rel_err={rel_err*100:.1f}%  [{status}]")

    # 6. Write CSV
    RESULTS_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[simulate] Wrote {RESULTS_CSV} ({len(rows)} rows)")

    # Clean up temp h² output files
    for f in SUMSTATS_DIR.glob(".tmp_*"):
        f.unlink()


def read_M_5_50(prefix: Path) -> int:
    """Read .l2.M_5_50 (single integer)."""
    with open(f"{prefix}.l2.M_5_50") as f:
        return int(f.read().strip().split()[0])


def mean_ld_score(prefix: Path) -> float:
    """Read .l2.ldscore.gz and return mean L2."""
    vals = []
    with gzip.open(f"{prefix}.l2.ldscore.gz", "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        l2_col = header.index("L2")
        for line in f:
            vals.append(float(line.rstrip("\n").split("\t")[l2_col]))
    return float(np.mean(vals))


# ---------------------------------------------------------------------------
# Subcommand: aggregate — print results table + run remaining sanity checks
# ---------------------------------------------------------------------------

def cmd_aggregate(args) -> None:
    if not RESULTS_CSV.exists():
        sys.exit(f"Missing {RESULTS_CSV}. Run `simulate` first.")

    rows = []
    with open(RESULTS_CSV) as f:
        for r in csv.DictReader(f):
            rows.append({**r,
                         "rep": int(r["rep"]),
                         "h2_true": float(r["h2_true"]),
                         "h2_est": float(r["h2_est"]),
                         "h2_se": float(r["h2_se"]),
                         "intercept": float(r["intercept"]),
                         })

    print(f"\n=== Aggregate over {len(rows)} h² runs ===\n")

    summary = []  # (h2_true, est_chunked±se, est_masked±se, paired_diff, p_value)
    for h2_true in sorted(set(r["h2_true"] for r in rows)):
        per_rep = {}  # rep → {chunked, masked}
        for r in rows:
            if r["h2_true"] != h2_true:
                continue
            per_rep.setdefault(r["rep"], {})[r["masking"]] = r["h2_est"]
        # paired diffs (masked - chunked)
        diffs = [d["masked"] - d["chunked"] for d in per_rep.values()
                 if "masked" in d and "chunked" in d]
        chunked = np.array([d["chunked"] for d in per_rep.values() if "chunked" in d])
        masked = np.array([d["masked"] for d in per_rep.values() if "masked" in d])
        n_pairs = len(diffs)
        mean_chunked, se_chunked = chunked.mean(), chunked.std(ddof=1) / math.sqrt(len(chunked))
        mean_masked, se_masked = masked.mean(), masked.std(ddof=1) / math.sqrt(len(masked))
        bias_chunked = mean_chunked - h2_true
        bias_masked = mean_masked - h2_true
        diff_arr = np.array(diffs)
        diff_mean = diff_arr.mean()
        diff_se = diff_arr.std(ddof=1) / math.sqrt(n_pairs)
        ci95 = (diff_mean - 1.96 * diff_se, diff_mean + 1.96 * diff_se)

        summary.append({
            "h2_true": h2_true,
            "n": n_pairs,
            "mean_chunked": mean_chunked,
            "se_chunked": se_chunked,
            "bias_chunked": bias_chunked,
            "rel_bias_chunked": bias_chunked / h2_true,
            "mean_masked": mean_masked,
            "se_masked": se_masked,
            "bias_masked": bias_masked,
            "rel_bias_masked": bias_masked / h2_true,
            "diff_mean": diff_mean,
            "diff_ci_lo": ci95[0],
            "diff_ci_hi": ci95[1],
        })

    # Sanity check #4: mean(masked) within H2_RECOVERY_TOL_SE × SE of true h²
    print("Sanity check #4: masked recovery (mean(h2_masked) close to true h²)")
    for s in summary:
        delta_se = abs(s["bias_masked"]) / s["se_masked"]
        ok = delta_se <= H2_RECOVERY_TOL_SE
        flag = "OK " if ok else "WARN"
        print(f"  h2_true={s['h2_true']:.2f}  mean_masked={s['mean_masked']:.4f}  "
              f"bias={s['bias_masked']:+.4f}  ({delta_se:.2f} SE)  [{flag}]")
    print()

    # Markdown table
    print("## Markdown table\n")
    print("| true h² | n |  ĥ²_chunked (SE)   |  ĥ²_masked (SE)    | bias_chunked | bias_masked | paired Δ (masked−chunked) [95% CI] |")
    print("|---------|---|--------------------|--------------------|--------------|-------------|------------------------------------|")
    for s in summary:
        print(f"| {s['h2_true']:.2f}   | {s['n']:2d} | {s['mean_chunked']:.4f} ({s['se_chunked']:.4f}) "
              f"| {s['mean_masked']:.4f} ({s['se_masked']:.4f}) "
              f"| {s['bias_chunked']:+.4f} ({100*s['rel_bias_chunked']:+.1f}%) "
              f"| {s['bias_masked']:+.4f} ({100*s['rel_bias_masked']:+.1f}%) "
              f"| {s['diff_mean']:+.4f} [{s['diff_ci_lo']:+.4f}, {s['diff_ci_hi']:+.4f}] |")
    print()

    # Typst fragment
    print("## Typst fragment for preprint/main.typ\n")
    print("```typst")
    print("#figure(")
    print("  table(")
    print("    columns: 5,")
    print("    align: (right, right, right, right, right),")
    print("    table.header(")
    print('      [*True $h^2$*], [*$hat(h)^2$ chunked*], [*$hat(h)^2$ masked*],')
    print('      [*Bias chunked*], [*Bias masked*],')
    print("    ),")
    print("    table.hline(stroke: 0.5pt),")
    for s in summary:
        print(f"    [{s['h2_true']:.2f}], "
              f"[{s['mean_chunked']:.4f} ± {s['se_chunked']:.4f}], "
              f"[{s['mean_masked']:.4f} ± {s['se_masked']:.4f}], "
              f"[{100*s['rel_bias_chunked']:+.1f}%], "
              f"[{100*s['rel_bias_masked']:+.1f}%],")
    print("  ),")
    print("  caption: [")
    print(f"    Recovery of true $h^2$ from {N_REPLICATES} simulated phenotypes per condition.")
    print(f"    Chromosome 22 of 1000 Genomes EUR-only cohort ($N = 503$ unrelated Europeans,")
    print(f"    ~18.6K SNPs after MAF $>= {MAF_MIN}$), `--ld-wind-kb {LD_WIND_KB}` (canonical LDSC")
    print(f"    window, $\\approx$1 cM). Infinitesimal genetic architecture. Z-scores generated by")
    print("    PLINK2 `--glm --linear`, converted to sumstats via `ldsc munge-sumstats`")
    print("    (using `statrs::erfc_inv` for canonical $P \\to Z$). Both estimators show")
    print("    upward bias at this sample size; the paired difference between masked and")
    print("    chunked LD scores is small ($\\sim$1% of $h^2$).")
    print("  ],")
    print(") <tbl:masking-bias>")
    print("```")
    print()

    # Sanity check #3 was based on a flawed assumption (that 1kb window would be
    # too small for chunked vs masked to differ). In fact at narrow windows, the
    # chunked-vs-masked gap is LARGEST because the chunk size c=200 dwarfs the
    # number of true in-window SNPs, so chunked over-counts heavily.
    # The check has been removed; the meaningful sanity check is #2 (median χ²
    # matches LDSC theory), which validates the simulation harness end-to-end.
    if (LDSCORE_DIR / "control_chunked.l2.ldscore.gz").exists():
        c = mean_ld_score(LDSCORE_DIR / "control_chunked")
        m = mean_ld_score(LDSCORE_DIR / "control_masked")
        print(f"Reference: 1kb-window LD scores (NOT a sanity check, just informational): "
              f"mean(L2_chunked)={c:.3f}, mean(L2_masked)={m:.3f}, "
              f"chunked/masked={c/m:.2f}x  (large ratio expected at narrow windows)")
        print()


# ---------------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------------

def main() -> None:
    desc = (__doc__ or "").strip().splitlines()[0] if __doc__ else "h² masking simulation"
    p = argparse.ArgumentParser(description=desc)
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("prep", help="subset chr22 + compute LD scores (chunked + masked)")
    sub.add_parser("simulate", help="generate sumstats + run h² regressions")
    sub.add_parser("aggregate", help="parse CSV + print results table")
    sub.add_parser("all", help="run prep, simulate, aggregate end-to-end")
    args = p.parse_args()

    if args.cmd == "all":
        cmd_prep(args)
        cmd_simulate(args)
        cmd_aggregate(args)
    elif args.cmd == "prep":
        cmd_prep(args)
    elif args.cmd == "simulate":
        cmd_simulate(args)
    elif args.cmd == "aggregate":
        cmd_aggregate(args)


if __name__ == "__main__":
    main()
