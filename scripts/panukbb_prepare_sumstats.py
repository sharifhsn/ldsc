#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["polars"]
# ///
"""
Prepare Pan-UKBB per-phenotype summary stats for LDSC munge_sumstats.

Outputs a gzipped TSV with columns: SNP, A1, A2, BETA, SE, P, N
"""
from __future__ import annotations

import argparse
import gzip
import sys
from typing import Optional

import polars as pl

POP_CHOICES = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID"]
MODE_CHOICES = ["pop", "meta", "meta_hq"]
PVAL_ENCODING = ["auto", "neglog10", "ln", "raw"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert Pan-UKBB per-phenotype files to LDSC-ready TSV.gz"
    )
    parser.add_argument("--input", required=True, help="Path to .tsv.bgz input")
    parser.add_argument("--output", required=True, help="Path to output .tsv.gz")
    parser.add_argument("--mode", choices=MODE_CHOICES, default="pop")
    parser.add_argument("--pop", choices=POP_CHOICES, default="EUR")
    parser.add_argument(
        "--pval-encoding",
        choices=PVAL_ENCODING,
        default="auto",
        help="How p-values are encoded in the file",
    )
    parser.add_argument(
        "--pval-col",
        default=None,
        help="Override p-value column name (e.g., neglog10_pval_EUR)",
    )
    parser.add_argument(
        "--beta-col",
        default=None,
        help="Override beta column name (e.g., beta_EUR)",
    )
    parser.add_argument(
        "--se-col",
        default=None,
        help="Override SE column name (e.g., se_EUR)",
    )
    parser.add_argument(
        "--n",
        type=float,
        default=None,
        help="Sample size to use when no N column is present",
    )
    parser.add_argument(
        "--n-col",
        default=None,
        help="Override N column name if present",
    )
    parser.add_argument(
        "--drop-low-confidence",
        action="store_true",
        default=False,
        help="Drop rows flagged as low_confidence for the chosen population",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=1_000_000,
        help="Row batch size for streaming read",
    )
    parser.add_argument(
        "--log-every",
        type=int,
        default=2_000_000,
        help="Log progress every N rows (0 to disable)",
    )
    return parser.parse_args()


def infer_cols(mode: str, pop: str) -> tuple[str, str, str, Optional[str]]:
    if mode == "pop":
        beta_col = f"beta_{pop}"
        se_col = f"se_{pop}"
        pval_col = f"neglog10_pval_{pop}"
    elif mode == "meta":
        beta_col = "beta_meta"
        se_col = "se_meta"
        pval_col = "neglog10_pval_meta"
    else:
        beta_col = "beta_meta_hq"
        se_col = "se_meta_hq"
        pval_col = "neglog10_pval_meta_hq"
    return beta_col, se_col, pval_col, None


def infer_pval_encoding(pval_col: str, encoding: str) -> str:
    if encoding != "auto":
        return encoding
    lower = pval_col.lower()
    if "neglog10" in lower:
        return "neglog10"
    if "ln" in lower or "log_p" in lower or "logp" in lower:
        return "ln"
    return "raw"


def pval_expr(pval_col: str, encoding: str) -> pl.Expr:
    col = pl.col(pval_col)
    if encoding == "neglog10":
        return pl.lit(10.0).pow(-col)
    if encoding == "ln":
        return col.exp()
    return col


def main() -> int:
    args = parse_args()

    beta_col, se_col, pval_col, _ = infer_cols(args.mode, args.pop)
    beta_col = args.beta_col or beta_col
    se_col = args.se_col or se_col
    pval_col = args.pval_col or pval_col

    n_col = args.n_col
    low_conf_col = f"low_confidence_{args.pop}"

    required = ["chr", "pos", "ref", "alt", beta_col, se_col, pval_col]
    if n_col:
        required.append(n_col)
    if args.drop_low_confidence:
        required.append(low_conf_col)

    encoding = infer_pval_encoding(pval_col, args.pval_encoding)

    reader = pl.read_csv_batched(
        args.input,
        separator="\t",
        has_header=True,
        columns=required,
        batch_size=args.batch_size,
        schema_overrides={"chr": pl.Utf8},
    )

    total_rows = 0
    if args.output.endswith((".gz", ".bgz")):
        out_handle = gzip.open(args.output, "wb")
    else:
        out_handle = open(args.output, "wb")

    with out_handle as out:
        first = True
        while True:
            batches = reader.next_batches(1)
            if not batches:
                break
            df = batches[0].with_columns(
                [
                    pl.col("chr").cast(pl.Utf8),
                    pl.col("pos").cast(pl.Utf8),
                    pl.col("ref").cast(pl.Utf8),
                    pl.col("alt").cast(pl.Utf8),
                    pl.col(beta_col).cast(pl.Float64, strict=False),
                    pl.col(se_col).cast(pl.Float64, strict=False),
                    pl.col(pval_col).cast(pl.Float64, strict=False),
                ]
            )
            df = df.drop_nulls([beta_col, se_col, pval_col])

            if args.drop_low_confidence:
                df = df.filter(
                    pl.col(low_conf_col)
                    .cast(pl.Utf8, strict=False)
                    .str.to_lowercase()
                    .ne("true")
                )

            n_expr: pl.Expr
            if n_col:
                n_expr = pl.col(n_col)
            elif args.n is not None:
                n_expr = pl.lit(float(args.n))
            else:
                raise SystemExit("No N column found. Provide --n or --n-col.")

            df = df.with_columns(
                [
                    pl.concat_str(
                        [pl.col("chr"), pl.col("pos"), pl.col("ref"), pl.col("alt")],
                        separator=":",
                    ).alias("SNP"),
                    pl.col("alt").alias("A1"),
                    pl.col("ref").alias("A2"),
                    pl.col(beta_col).alias("BETA"),
                    pl.col(se_col).alias("SE"),
                    pval_expr(pval_col, encoding).clip(1e-300, 1.0).alias("P"),
                    n_expr.alias("N"),
                ]
            ).select(["SNP", "A1", "A2", "BETA", "SE", "P", "N"])

            df.write_csv(out, separator="\t", include_header=first)
            first = False

            total_rows += df.height
            if args.log_every and total_rows % args.log_every < df.height:
                print(f"processed {total_rows} rows", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
