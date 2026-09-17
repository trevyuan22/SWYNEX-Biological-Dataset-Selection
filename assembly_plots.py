import sys
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt
import numpy as np

def plot(fasta, label):
    lengths, gcs = [], []
    for rec in SeqIO.parse(fasta, "fasta"):
        lengths.append(len(rec.seq))
        gcs.append(gc_fraction(rec.seq) * 100)

    lengths = np.array(lengths)
    gcs = np.array(gcs)

    fig, ax = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Assembly Diagnostics — {label}", fontsize=14, fontweight="bold")

    # 1. Length distribution
    ax[0, 0].hist(lengths, bins=50, color="steelblue", edgecolor="black")
    ax[0, 0].set_xscale("log")
    ax[0, 0].set_title("Length distribution (log)")
    ax[0, 0].set_xlabel("Sequence length (bp)")
    ax[0, 0].set_ylabel("Count")

    # 2. Nx curve
    s = np.sort(lengths)[::-1]
    cum = np.cumsum(s) / s.sum() * 100
    ax[0, 1].plot(np.arange(1, len(s)+1), cum, color="darkorange")
    ax[0, 1].axhline(50, color="red", linestyle="--", label="N50 threshold")
    ax[0, 1].set_xscale("log")
    ax[0, 1].set_title("Nx curve")
    ax[0, 1].set_xlabel("Sequence rank")
    ax[0, 1].set_ylabel("Cumulative length (%)")
    ax[0, 1].legend()

    # 3. GC content distribution
    ax[1, 0].hist(gcs, bins=50, color="seagreen", edgecolor="black")
    ax[1, 0].axvline(gcs.mean(), color="red", linestyle="--", label=f"Mean = {gcs.mean():.2f}%")
    ax[1, 0].set_title("GC content distribution")
    ax[1, 0].set_xlabel("GC content (%)")
    ax[1, 0].set_ylabel("Count")
    ax[1, 0].legend()

    # 4. Length vs GC scatter
    ax[1, 1].scatter(lengths, gcs, alpha=0.4, s=10, color="purple")
    ax[1, 1].set_xscale("log")
    ax[1, 1].set_title("Length vs GC (contamination check)")
    ax[1, 1].set_xlabel("Sequence length (bp)")
    ax[1, 1].set_ylabel("GC content (%)")

    plt.tight_layout()
    out = label.replace(" ", "_") + "_plots.png"
    plt.savefig(out, dpi=200)
    print(f"Saved: {out}")

if __name__ == "__main__":
    plot(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "Assembly")