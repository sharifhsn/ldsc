#!/usr/bin/env python3
"""
Generate synthetic biobank-scale PLINK BED/BIM/FAM from an existing dataset.

Replicates individuals K times with ~1% genotype noise to create a larger N
while preserving realistic LD structure. BED format is SNP-major.

Uses numpy for vectorized byte manipulation — generates ~20GB in ~10 minutes.

Usage:
    python scripts/make_synthetic_biobank.py \
        --bfile data/bench_5k \
        --target-n 50000 \
        --noise-rate 0.01 \
        --seed 42 \
        --out data/bench_5k_50k
"""

import argparse
import os
import shutil
import sys
import time

import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Generate synthetic biobank PLINK data")
    parser.add_argument("--bfile", required=True, help="Input PLINK prefix")
    parser.add_argument("--target-n", type=int, required=True, help="Target number of individuals")
    parser.add_argument("--noise-rate", type=float, default=0.01, help="Fraction of genotypes to flip (default: 0.01)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--out", required=True, help="Output PLINK prefix")
    args = parser.parse_args()

    rng = np.random.RandomState(args.seed)

    # Read FAM to get original N
    with open(f"{args.bfile}.fam") as f:
        fam_lines = f.readlines()
    n_orig = len(fam_lines)

    # Read BIM to get M
    with open(f"{args.bfile}.bim") as f:
        bim_lines = f.readlines()
    n_snps = len(bim_lines)

    # Compute replication factor
    k = (args.target_n + n_orig - 1) // n_orig
    n_new = min(k * n_orig, args.target_n)

    bytes_per_snp_orig = (n_orig + 3) // 4
    bytes_per_snp_new = (n_new + 3) // 4
    bed_size_gb = (3 + n_snps * bytes_per_snp_new) / (1024**3)

    print(f"Source: {n_orig} individuals, {n_snps} SNPs")
    print(f"Replication factor: {k} → {n_new} individuals")
    print(f"Noise rate: {args.noise_rate}")
    print(f"Output BED: {bed_size_gb:.2f} GB ({bytes_per_snp_new} bytes/SNP)")

    # Write FAM
    print("Writing FAM...")
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(f"{args.out}.fam", "w") as fout:
        count = 0
        for rep in range(k):
            for line in fam_lines:
                if count >= n_new:
                    break
                parts = line.strip().split()
                parts[0] = f"FAM{rep}"
                parts[1] = f"IID{rep}_{parts[1]}"
                fout.write("\t".join(parts) + "\n")
                count += 1

    # Copy BIM
    print("Writing BIM...")
    shutil.copy2(f"{args.bfile}.bim", f"{args.out}.bim")

    # Read entire source BED into memory
    print("Reading source BED...")
    with open(f"{args.bfile}.bed", "rb") as f:
        magic = f.read(3)
        assert magic == b'\x6c\x1b\x01', "Invalid BED file"
        bed_data = np.frombuffer(f.read(), dtype=np.uint8)
    bed_data = bed_data.reshape(n_snps, bytes_per_snp_orig)

    # Process in batches of SNPs to limit memory
    batch_size = 1000  # SNPs per batch
    t0 = time.time()

    print(f"Writing BED ({n_snps} SNPs)...")
    with open(f"{args.out}.bed", "wb") as fout:
        fout.write(b'\x6c\x1b\x01')

        for batch_start in range(0, n_snps, batch_size):
            batch_end = min(batch_start + batch_size, n_snps)
            batch_m = batch_end - batch_start

            # Extract 2-bit genotypes for this batch: (batch_m, n_orig)
            raw = bed_data[batch_start:batch_end]  # (batch_m, bytes_per_snp_orig)

            # Unpack bytes to 2-bit genotypes
            # Each byte has 4 genotypes: bits [1:0], [3:2], [5:4], [7:6]
            geno = np.zeros((batch_m, n_orig), dtype=np.uint8)
            for bit_pair in range(4):
                start_idx = bit_pair * (n_orig // 4 + (1 if bit_pair < n_orig % 4 else 0))
                # Simpler: unpack all 4 positions from each byte
                pass

            # Actually, let's unpack properly using bit shifts
            # Expand each byte to 4 genotypes
            expanded = np.zeros((batch_m, bytes_per_snp_orig * 4), dtype=np.uint8)
            for bp in range(4):
                expanded[:, bp::4] = (raw >> (bp * 2)) & 0x03
            geno = expanded[:, :n_orig]  # trim to actual n_orig

            # Replicate K times: (batch_m, n_new)
            geno_rep = np.tile(geno, (1, k))[:, :n_new]  # (batch_m, n_new)

            # Add noise to copies (not the first copy)
            if args.noise_rate > 0:
                # Noise mask: which entries to flip (only in replicated copies, not first)
                noise_mask = np.zeros((batch_m, n_new), dtype=bool)
                noise_mask[:, n_orig:] = rng.random((batch_m, n_new - n_orig)) < args.noise_rate

                # Don't flip missing values (01 = 0x01)
                is_missing = geno_rep == 0x01
                noise_mask &= ~is_missing

                # Noise: 00->10, 10->{00 or 11}, 11->10
                # 00 (hom ref) -> 10 (het)
                flip_00 = noise_mask & (geno_rep == 0x00)
                geno_rep[flip_00] = 0x02  # het

                # 10 (het) -> 00 or 11
                flip_10 = noise_mask & (geno_rep == 0x02)
                het_choice = rng.random(np.sum(flip_10)) < 0.5
                vals = np.where(het_choice, 0x00, 0x03)
                geno_rep[flip_10] = vals

                # 11 (hom alt) -> 10 (het)
                flip_11 = noise_mask & (geno_rep == 0x03)
                geno_rep[flip_11] = 0x02  # het

            # Pack back to bytes: 4 genotypes per byte
            # Pad to multiple of 4
            padded_n = bytes_per_snp_new * 4
            if padded_n > n_new:
                geno_padded = np.zeros((batch_m, padded_n), dtype=np.uint8)
                geno_padded[:, :n_new] = geno_rep
            else:
                geno_padded = geno_rep

            # Pack: combine 4 consecutive 2-bit values into 1 byte
            packed = (
                geno_padded[:, 0::4] |
                (geno_padded[:, 1::4] << 2) |
                (geno_padded[:, 2::4] << 4) |
                (geno_padded[:, 3::4] << 6)
            ).astype(np.uint8)

            fout.write(packed.tobytes())

            elapsed = time.time() - t0
            snps_done = batch_end
            rate = snps_done / elapsed if elapsed > 0 else 0
            eta = (n_snps - snps_done) / rate if rate > 0 else 0
            print(f"  {snps_done}/{n_snps} SNPs ({snps_done*100//n_snps}%), "
                  f"{elapsed:.1f}s elapsed, ETA {eta:.0f}s", end="\r")

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.1f}s. Output: {args.out}.{{bed,bim,fam}}")
    print(f"  {n_new} individuals, {n_snps} SNPs")
    print(f"  BED size: {bed_size_gb:.2f} GB")


if __name__ == "__main__":
    main()
