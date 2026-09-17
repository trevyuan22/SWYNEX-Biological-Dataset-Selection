#!/usr/bin/env python3
"""
Compute basic assembly statistics for Leptodirus hochenwartii.
Usage: python assembly_stats.py <fasta_file> <label>
"""

import sys
import os
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import statistics


def n50_l50(lengths):
    """Return N50 and L50 for a list of sequence lengths."""
    s = sorted(lengths, reverse=True)
    half = sum(s) / 2
    cum = 0
    for i, l in enumerate(s, start=1):
        cum += l
        if cum >= half:
            return l, i
    return 0, 0


def nx(lengths, x):
    """Return Nx for percentile x (e.g., nx(lengths, 75) = N75)."""
    s = sorted(lengths, reverse=True)
    target = sum(s) * (x / 100)
    cum = 0
    for l in s:
        cum += l
        if cum >= target:
            return l
    return 0


def aun(lengths):
    """auN: threshold-free contiguity metric."""
    total = sum(lengths)
    return sum(l * l for l in lengths) / total if total else 0


def analyze(fasta, label):
    print("=" * 70)
    print(f"  {label}")
    print(f"  File: {os.path.basename(fasta)}")
    print("=" * 70)

    lengths = []
    counts = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0, "other": 0}

    for rec in SeqIO.parse(fasta, "fasta"):
        L = len(rec.seq)
        lengths.append(L)
        for base in rec.seq:
            b = base.upper()
            if b in counts:
                counts[b] += 1
            else:
                counts["other"] += 1

    total = sum(lengths)
    n = len(lengths)

    print("\n--- Contiguity ---")
    print(f"  Number of sequences           : {n:,}")
    print(f"  Total assembly length (bp)    : {total:,}")
    print(f"  Longest sequence (bp)         : {max(lengths):,}")
    print(f"  Shortest sequence (bp)        : {min(lengths):,}")
    print(f"  Mean sequence length (bp)     : {statistics.mean(lengths):,.0f}")
    print(f"  Median sequence length (bp)   : {statistics.median(lengths):,.0f}")
    print(f"  N50 (bp)                      : {n50_l50(lengths)[0]:,}")
    print(f"  L50 (count)                   : {n50_l50(lengths)[1]:,}")
    print(f"  N75 (bp)                      : {nx(lengths, 75):,}")
    print(f"  N90 (bp)                      : {nx(lengths, 90):,}")
    print(f"  auN (bp)                      : {aun(lengths):,.0f}")

    gc = counts["G"] + counts["C"]
    gc_pct = (gc / total * 100) if total else 0

    print("\n--- Composition ---")
    print(f"  GC content (%)                : {gc_pct:.2f}")
    print(f"  Number of N bases             : {counts['N']:,} "
          f"({counts['N']/total*100:.2f}%)")

    print("\n--- Nucleotide Composition ---")
    for b in ["A", "T", "G", "C", "N", "other"]:
        pct = counts[b] / total * 100 if total else 0
        print(f"  {b:>5}: {counts[b]:>15,}  ({pct:6.2f}%)")

    return {
        "label": label, "n": n, "total": total,
        "longest": max(lengths), "shortest": min(lengths),
        "mean": statistics.mean(lengths), "median": statistics.median(lengths),
        "n50": n50_l50(lengths)[0], "l50": n50_l50(lengths)[1],
        "n75": nx(lengths, 75), "n90": nx(lengths, 90),
        "aun": aun(lengths), "gc": gc_pct, "counts": counts,
    }


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python assembly_stats.py <fasta_file> <label>")
        sys.exit(1)
    analyze(sys.argv[1], sys.argv[2])