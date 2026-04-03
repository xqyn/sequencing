'''
project: sequencing
XQ - Leiden UMC
UMI analyu
update:
    - april 2: init
    
    
Usage
-----
python3 umi_qc_analysis.py \
    --fastq        reads_with_umi.fastq.gz \
    --input_raw    raw.fastq.gz \
    --outdir       qc_plots/ \
    --umi_len      12 \
    --umi_pos      5 \
    --subsample    200000 \
    --threads      4
'''

import argparse
import itertools
import logging
import math
import os
import random
import sys
from collections import Counter, defaultdict
from typing import Literal
 
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
 
# custom FASTQ / sequence helpers (project-local)
from fastq import SeqQuery
from sequence import ascii_to_phred, hamming


# ── logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="  %(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ── styling ───────────────────────────────────────────────────────────────────
PALETTE = {
    "bg":      "#FFFFFF",
    "panel":   "#FFFFFF",
    "accent1": "#4fc3f7",
    "accent2": "#f06292",
    "accent3": "#a5d6a7",
    "accent4": "#ffb74d",
    "text":    "#000000",
    "grid":    "#cccccc",
}
NCOLORS = {"A": "#4fc3f7", 
           "C": "#f06292", 
           "G": "#a5d6a7", 
           "T": "#ffb74d", 
           "N": "#9e9e9e"}
BASE_CMAP = LinearSegmentedColormap.from_list(
    "umi", ["#e0f0ff", "#4fc3f7", "#f06292"], N=256
)
 
 
def apply_style():
    plt.rcParams.update({
        "figure.facecolor":  PALETTE["bg"],
        "axes.facecolor":    PALETTE["panel"],
        "axes.edgecolor":    PALETTE["grid"],
        "axes.labelcolor":   PALETTE["text"],
        "xtick.color":       PALETTE["text"],
        "ytick.color":       PALETTE["text"],
        "text.color":        PALETTE["text"],
        "grid.color":        PALETTE["grid"],
        "grid.linestyle":    "--",
        "grid.alpha":        0.5,
        "font.family":       "monospace",
        "axes.titlesize":    11,
        "axes.labelsize":    9,
        "xtick.labelsize":   8,
        "ytick.labelsize":   8,
        "legend.fontsize":   8,
        "savefig.dpi":       150,
        "savefig.bbox":      "tight",
        "savefig.facecolor": PALETTE["bg"],
    })
 
 
apply_style()

# ── plot helpers ──────────────────────────────────────────────────────────────
def _ax_spine(ax):
    for spine in ax.spines.values():
        spine.set_edgecolor(PALETTE["grid"])
 
 
def _save(fig, path):
    fig.savefig(path)
    plt.close(fig)
    log.info("  Saved: %s", path)
    

# =============================================================================
# Data collection
# =============================================================================

def collect_umi_data(fastq_path: str,
                     umi_len: int = 8,
                     umi_end: Literal[3, 5] = 5,
                     subsample: int | None = None,
                     min_qual: int = 0,
                     seed: int = 920) -> dict:
    """
    Stream through a FASTQ file once and collect raw UMI quality-control data.
 
    Applies a per-read quality gate, performs Algorithm-L reservoir sampling,
    and returns raw arrays for downstream vectorised analysis.
 
    Args:
        fastq_path  (str)            : path to the FASTQ file (.fastq or .fastq.gz).
        umi_len     (int)            : expected UMI length in bases (default: 8).
        umi_end     (Literal[3, 5]) : UMI position; 5 for 5'-end, 3 for 3'-end (default: 5).
        subsample   (int | None)     : reservoir size for Algorithm-L sampling; None disables.
        min_qual    (int)            : minimum mean Phred score to pass the quality gate (default: 0).
        seed        (int)            : random seed for reproducibility (default: 920).
 
    Returns:
        dict with keys: total_reads, umi_seqs, umi_phreds, umi_lengths,
        per_base_nuc, per_base_qual, gc, n_content, umi_mean_qual,
        read_mean_qual, homopolymer_lens, reservoir.
 
    Raises:
        ValueError: If umi_end is not 3 or 5.
    """
    random.seed(seed)
    np.random.seed(seed)
 
    if umi_end not in (3, 5):
        raise ValueError(f"umi_end must be 3 or 5, got {umi_end}")
 
    # ── raw stores ────────────────────────────────────────────────────────────
    umi_seqs          = []
    umi_phreds        = []
    umi_lengths       = []
    per_base_nuc      = defaultdict(Counter)
    per_base_qual_sum = defaultdict(float)
    per_base_qual_cnt = defaultdict(int)
    gc_perc_list      = []
    n_perc_list       = []
    umi_mean_qual     = []
    read_mean_qual    = []      # FIX 6: was missing entirely
    homopolymer_lens  = []
    total_reads       = 0
 
    # ── Algorithm L state ─────────────────────────────────────────────────────
    reservoir  = []             # holds (umi_seq, umi_qual, q_mean)
    _W         = None
    _next_skip = None
    _sampling  = subsample is not None
    n_passing  = 0              # FIX 9: track post-filter count for Algorithm L
 
    log.info("Streaming FASTQ: %s", fastq_path)
 
    for ind, query in enumerate(SeqQuery(fastq_filename=fastq_path)):
        total_reads += 1
 
        if umi_end == 5:
            umi_seq  = query.seq[:umi_len]
            umi_qual = query.qual[:umi_len]
        else:
            umi_seq  = query.seq[-umi_len:]
            umi_qual = query.qual[-umi_len:]
 
        umi_phred = ascii_to_phred(umi_qual)
        q_mean    = float(np.mean(umi_phred))
 
        # quality gate
        if q_mean < min_qual:
            continue
 
        n_passing += 1          # FIX 9: increment after gate
 
        # ── store raw data ────────────────────────────────────────────────────
        umi_seqs.append(umi_seq)
        umi_phreds.append(umi_phred)
        umi_lengths.append(len(umi_seq))
 
        # ── per-base composition & quality ────────────────────────────────────
        for pos, (nuc, qval) in enumerate(zip(umi_seq, umi_phred)):
            per_base_nuc[pos][nuc]  += 1
            per_base_qual_sum[pos]  += qval
            per_base_qual_cnt[pos]  += 1
 
        # ── summary metrics ───────────────────────────────────────────────────
        gc_perc = (umi_seq.count("G") + umi_seq.count("C")) / max(len(umi_seq), 1)
        gc_perc_list.append(gc_perc)
 
        n_perc = umi_seq.count("N") / max(len(umi_seq), 1)
        n_perc_list.append(n_perc)
 
        umi_mean_qual.append(q_mean)
        read_mean_qual.append(float(np.mean(ascii_to_phred(query.qual))))   # FIX 6
 
        homopolymer = max(
            (sum(1 for _ in grp) for _, grp in itertools.groupby(umi_seq)),
            default=0,
        )
        homopolymer_lens.append(homopolymer)
 
        # ── Algorithm L reservoir sampling ────────────────────────────────────
        if _sampling:
            record = (umi_seq, umi_qual, q_mean)
 
            if len(reservoir) < subsample:
                reservoir.append(record)
 
                if len(reservoir) == subsample:
                    # FIX 9: use n_passing (post-filter index) not ind
                    _W         = math.exp(math.log(random.random()) / subsample)
                    _next_skip = n_passing + math.floor(
                        math.log(random.random()) / math.log(1 - _W)
                    ) + 1
 
            elif n_passing == _next_skip:       # FIX 9: compare against n_passing
                reservoir[random.randint(0, subsample - 1)] = record
                _W         *= math.exp(math.log(random.random()) / subsample)
                _next_skip += math.floor(
                    math.log(random.random()) / math.log(1 - _W)
                ) + 1
 
        if total_reads % 100_000 == 0:
            log.info("  Processed %d reads …", total_reads)
 
    log.info("Total reads processed: %d  |  passing filter: %d", total_reads, n_passing)
 
    return {
        "total_reads":      total_reads,
        "umi_seqs":         umi_seqs,
        "umi_phreds":       umi_phreds,
        "umi_lengths":      umi_lengths,
        "per_base_nuc":     dict(per_base_nuc),
        "per_base_qual":    {p: per_base_qual_sum[p] / per_base_qual_cnt[p]
                             for p in per_base_qual_sum},
        "gc":               gc_perc_list,
        "n_content":        n_perc_list,
        "umi_mean_qual":    umi_mean_qual,
        "read_mean_qual":   read_mean_qual,     # FIX 6
        "homopolymer_lens": homopolymer_lens,
        "reservoir":        reservoir,
    }
 
    
# =============================================================================
# Derived metrics
# =============================================================================
 
def saturation_curve(umi_seqs: list, n_points: int = 30):
    """Downsample to increasing fractions, count unique UMIs at each."""
    n     = len(umi_seqs)
    sizes = np.unique(np.linspace(1, n, n_points, dtype=int))
    arr   = np.array(umi_seqs)
    unique = []
    for s in sizes:
        idx = np.random.choice(n, size=s, replace=False)
        unique.append(len(set(arr[idx])))
    return sizes, unique
 
 
def sample_hamming(umi_list: list, n_pairs: int = 50_000, seed: int = 42) -> list:
    """Sample random pairs and compute Hamming distances."""
    rng  = np.random.default_rng(seed)
    arr  = np.array(umi_list)
    n    = len(arr)
    idxA = rng.integers(0, n, n_pairs)
    idxB = rng.integers(0, n, n_pairs)
    return [hamming(arr[i], arr[j]) for i, j in zip(idxA, idxB) if i != j]
 
 
def positional_entropy(per_base_nuc: dict) -> dict:
    """Shannon entropy (bits) at each UMI position."""
    entropy = {}
    for pos, cnt in per_base_nuc.items():
        total = sum(cnt.values())
        H = -sum((v / total) * math.log2(v / total) for v in cnt.values() if v > 0)
        entropy[pos] = H
    return entropy


# =============================================================================
# Plotting functions
# =============================================================================
 
# ── 01. UMI length distribution ───────────────────────────────────────────────
def plot_length_dist(lengths: list, umi_len: int, outdir: str):
    """FIX 5: function was called in main() but missing from file."""
    cnt  = Counter(lengths)
    x, y = zip(*sorted(cnt.items()))
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(x, y, color=PALETTE["accent1"], width=0.6)
    ax.axvline(umi_len, color=PALETTE["accent2"], lw=1.5, ls="--",
               label=f"expected ({umi_len} nt)")
    ax.set_xlabel("UMI length (nt)")
    ax.set_ylabel("Read count")
    ax.set_title("UMI Length Distribution")
    ax.legend()
    ax.grid(axis="y")
    _ax_spine(ax)
    _save(fig, f"{outdir}/01_umi_length_distribution.png")


# ── 02. Per-base nucleotide composition ───────────────────────────────────────
def plot_per_base_composition(per_base_nuc: dict, outdir: str):
    positions = sorted(per_base_nuc.keys())
    data = {b: [] for b in "ACGTN"}
    for pos in positions:
        total = sum(per_base_nuc[pos].values())
        for b in "ACGTN":
            data[b].append(per_base_nuc[pos].get(b, 0) / total * 100)
 
    fig, ax = plt.subplots(figsize=(max(8, len(positions) * 0.45), 4))
    bottom = np.zeros(len(positions))
    for b in "ACGTN":
        vals = np.array(data[b])
        ax.bar(positions, vals, bottom=bottom, color=NCOLORS[b], label=b, width=0.8)
        bottom += vals
 
    ax.axhline(25, color=PALETTE["grid"], lw=0.7, ls=":", alpha=0.8, label="25% (uniform)")
    ax.set_xlim(-0.5, len(positions) - 0.5)
    ax.set_ylim(0, 100)
    ax.set_xlabel("UMI position")
    ax.set_ylabel("Nucleotide frequency (%)")
    ax.set_title("Per-base Nucleotide Composition")
    ax.legend(loc="upper right", ncol=5)
    ax.grid(axis="y")
    _ax_spine(ax)
    _save(fig, f"{outdir}/02_per_base_composition.png")


# ── 3. Per-base mean Phred quality ────────────────────────────────────────────
def plot_per_base_qual(per_base_qual: dict, outdir: str):
    positions = sorted(per_base_qual.keys())
    quals     = [per_base_qual[p] for p in positions]
    colors    = [PALETTE["accent3"] if q >= 30 else
                 PALETTE["accent4"] if q >= 20 else
                 PALETTE["accent2"] for q in quals]

    fig, ax = plt.subplots(figsize=(max(8, len(positions) * 0.45), 4))
    ax.bar(positions, quals, color=colors, width=0.8)
    ax.axhline(30, color=PALETTE["accent3"], lw=1.2, ls="--", alpha=0.7, label="Q30")
    ax.axhline(20, color=PALETTE["accent4"], lw=1.2, ls="--", alpha=0.7, label="Q20")
    ax.set_ylim(0, 42)
    ax.set_xlabel("UMI position")
    ax.set_ylabel("Mean Phred score")
    ax.set_title("Per-base UMI Quality Score")
    ax.legend()
    ax.grid(axis="y")
    _ax_spine(ax)
    _save(fig, f"{outdir}/03_per_base_quality.png")


# ── 03. Per-base mean Phred quality ───────────────────────────────────────────
def plot_per_base_qual(per_base_qual: dict, outdir: str):
    positions = sorted(per_base_qual.keys())
    quals     = [per_base_qual[p] for p in positions]
    colors    = [PALETTE["accent3"] if q >= 30 else
                 PALETTE["accent4"] if q >= 20 else
                 PALETTE["accent2"] for q in quals]
 
    fig, ax = plt.subplots(figsize=(max(8, len(positions) * 0.45), 4))
    ax.bar(positions, quals, color=colors, width=0.8)
    ax.axhline(30, color=PALETTE["accent3"], lw=1.2, ls="--", alpha=0.7, label="Q30")
    ax.axhline(20, color=PALETTE["accent4"], lw=1.2, ls="--", alpha=0.7, label="Q20")
    ax.set_ylim(0, 42)
    ax.set_xlabel("UMI position")
    ax.set_ylabel("Mean Phred score")
    ax.set_title("Per-base UMI Quality Score")
    ax.legend()
    ax.grid(axis="y")
    _ax_spine(ax)
    _save(fig, f"{outdir}/03_per_base_quality.png")   # FIX 3: missing ) was here
 
 
# ── 04. GC content distribution ───────────────────────────────────────────────
def plot_gc_content(gc_list: list, outdir: str):
    gc_arr  = np.array(gc_list) * 100
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(gc_arr, bins=50, color=PALETTE["accent1"], edgecolor="none", density=True)
    ax.axvline(50, color=PALETTE["accent2"], lw=1.5, ls="--", label="50% GC")
    mean_gc = gc_arr.mean()
    ax.axvline(mean_gc, color=PALETTE["accent4"], lw=1.5, ls="-",
               label=f"mean={mean_gc:.1f}%")
    ax.set_xlabel("GC content (%)")
    ax.set_ylabel("Density")
    ax.set_title("UMI GC-Content Distribution")
    ax.legend()
    ax.grid(axis="y")
    ax.set_xticks(range(0, 101, 5))
    ax.set_xlim(0, 100)
    _ax_spine(ax)
    _save(fig, f"{outdir}/04_gc_content.png")
 
 
# ── 05. Homopolymer distribution ──────────────────────────────────────────────
def plot_homopolymer(hp_lens: list, umi_len: int, outdir: str):
    cnt  = Counter(hp_lens)
    x, y = zip(*sorted(cnt.items()))
    total = sum(y)
 
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(x, [v / total * 100 for v in y],
           color=[PALETTE["accent2"] if xi >= umi_len * 0.5
                  else PALETTE["accent1"] for xi in x],
           width=0.7)
    ax.axvline(umi_len * 0.5, color=PALETTE["accent2"], lw=1.2, ls="--", alpha=0.6,
               label=f"≥50% homopolymer (>{umi_len // 2} nt)")
    ax.set_xlabel("Longest homopolymer run in UMI (nt)")
    ax.set_ylabel("Reads (%)")
    ax.set_title("Homopolymer Run Distribution")
    ax.legend()
    ax.grid(axis="y")
    _ax_spine(ax)
    _save(fig, f"{outdir}/05_homopolymer_distribution.png")
 
 
# ── 06. UMI mean quality histogram ────────────────────────────────────────────
def plot_umi_qual_hist(umi_qual: list, outdir: str):
    arr     = np.array(umi_qual)
    q30_pct = (arr >= 30).mean() * 100
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(arr, bins=80, color=PALETTE["accent3"], edgecolor="none", density=True)
    ax.axvline(30, color=PALETTE["accent2"], lw=1.5, ls="--", label="Q30")
    ax.set_xlabel("Mean Phred score (UMI region)")
    ax.set_ylabel("Density")
    ax.set_title(f"UMI Mean Quality Distribution  (Q30≥: {q30_pct:.1f}%)")
    ax.legend()
    ax.grid(axis="y")
    _ax_spine(ax)
    _save(fig, f"{outdir}/06_umi_mean_quality_hist.png")
 
 
# ── 07. Saturation curve ──────────────────────────────────────────────────────
def plot_saturation(umi_seqs: list, outdir: str):
    sizes, unique = saturation_curve(umi_seqs)
    sizes_m  = sizes / 1e6
    unique_m = np.array(unique) / 1e6
 
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(sizes_m, unique_m, color=PALETTE["accent1"], lw=2)
    ax.fill_between(sizes_m, unique_m, alpha=0.15, color=PALETTE["accent1"])
    ax.set_xlabel("Reads sampled (M)")
    ax.set_ylabel("Unique UMIs (M)")
    ax.set_title("Library Saturation Curve")
    ax.grid()
    _ax_spine(ax)
 
    final_rate = unique_m[-1] / sizes_m[-1] if sizes_m[-1] > 0 else 0
    ax.annotate(
        f"slope ≈ {final_rate:.2f}",
        xy=(sizes_m[-1], unique_m[-1]),
        xytext=(sizes_m[-1] * 0.6, unique_m[-1] * 0.9),
        color=PALETTE["accent4"],
        fontsize=8,
        arrowprops=dict(arrowstyle="->", color=PALETTE["accent4"]),
    )
    _save(fig, f"{outdir}/07_saturation_curve.png")
 
 
# ── 08. Hamming distance distribution ─────────────────────────────────────────
def plot_hamming(dists: list, umi_len: int, outdir: str):
    from scipy.stats import binom
    cnt  = Counter(dists)
    x, y = zip(*sorted(cnt.items()))
    total = sum(y)
 
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(x, [v / total * 100 for v in y], color=PALETTE["accent1"],
           width=0.7, label="observed")
    expected = [binom.pmf(d, umi_len, 0.75) * 100 for d in x]
    ax.plot(x, expected, color=PALETTE["accent2"], lw=1.5, ls="--",
            label="expected (random)")
    ax.set_xlabel("Hamming distance between UMI pairs")
    ax.set_ylabel("Pair fraction (%)")
    ax.set_title("Pairwise Hamming Distance Distribution (sub-sampled)")
    ax.legend()
    ax.grid(axis="y")
    _ax_spine(ax)
    _save(fig, f"{outdir}/08_hamming_distance.png")
 
 
# ── 09. Positional entropy ────────────────────────────────────────────────────
def plot_entropy(entropy: dict, outdir: str):
    positions = sorted(entropy.keys())
    H         = [entropy[p] for p in positions]
    max_H     = math.log2(4)
 
    fig, ax = plt.subplots(figsize=(max(8, len(positions) * 0.45), 4))
    ax.bar(positions, H, color=[
        PALETTE["accent3"] if h >= max_H * 0.9 else
        PALETTE["accent4"] if h >= max_H * 0.7 else
        PALETTE["accent2"] for h in H
    ], width=0.8)
    ax.axhline(max_H, color=PALETTE["grid"], lw=1, ls="--", alpha=0.8,
               label="max entropy (2 bits)")
    ax.set_ylim(0, max_H * 1.1)
    ax.set_xlabel("UMI position")
    ax.set_ylabel("Shannon entropy (bits)")
    ax.set_title("Per-position Entropy  (higher = more diverse / less biased)")
    ax.legend()
    ax.grid(axis="y")
    _ax_spine(ax)
    _save(fig, f"{outdir}/09_positional_entropy.png")
 
 
# ── 10. Read quality vs UMI quality scatter ───────────────────────────────────
def plot_qual_scatter(read_qual: list, umi_qual: list, outdir: str,
                      max_pts: int = 50_000):
    rq = np.array(read_qual)
    uq = np.array(umi_qual)
    if len(rq) > max_pts:
        idx = np.random.choice(len(rq), max_pts, replace=False)
        rq, uq = rq[idx], uq[idx]
 
    fig, ax = plt.subplots(figsize=(6, 6))
    h = ax.hexbin(rq, uq, gridsize=60, cmap=BASE_CMAP, mincnt=1, linewidths=0)
    plt.colorbar(h, ax=ax, label="Read count")
    ax.axhline(30, color=PALETTE["accent2"], lw=0.8, ls="--", alpha=0.7)
    ax.axvline(30, color=PALETTE["accent2"], lw=0.8, ls="--", alpha=0.7)
    ax.set_xlabel("Mean read quality (Phred)")
    ax.set_ylabel("Mean UMI quality (Phred)")
    ax.set_title("Read Quality vs UMI Quality")
    corr = np.corrcoef(rq, uq)[0, 1]
    ax.text(0.05, 0.95, f"r = {corr:.3f}",
            transform=ax.transAxes, color=PALETTE["text"], fontsize=9)
    _ax_spine(ax)
    _save(fig, f"{outdir}/10_read_vs_umi_quality_scatter.png")
 
 
# ── 11. N-content distribution ────────────────────────────────────────────────
def plot_n_content(n_list: list, outdir: str):
    arr         = np.array(n_list) * 100
    cnt_nonzero = (arr > 0).sum()
    total       = len(arr)
 
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(arr, bins=np.arange(0, 101, 2), color=PALETTE["accent2"],
            edgecolor="none", density=True)
    ax.set_xlabel("N bases in UMI (%)")
    ax.set_ylabel("Density")
    ax.set_title(f"UMI N-Content  ({cnt_nonzero / total * 100:.2f}% of UMIs contain ≥1 N)")
    ax.grid(axis="y")
    _ax_spine(ax)
    _save(fig, f"{outdir}/11_n_content.png")
 
 
# ── 12. Top over-represented UMIs ─────────────────────────────────────────────
def plot_top_umis(umi_seqs: list, outdir: str, top_n: int = 30):
    cnt   = Counter(umi_seqs)
    total = len(umi_seqs)
    top   = cnt.most_common(top_n)
    labels, counts = zip(*top)
    freqs = np.array(counts) / total * 100
 
    fig, ax = plt.subplots(figsize=(10, max(4, top_n * 0.28)))
    ax.barh(range(len(labels)), freqs[::-1], color=PALETTE["accent1"])
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels[::-1], fontsize=7)
    ax.set_xlabel("Frequency (%)")
    ax.set_title(f"Top {top_n} Most Frequent UMIs")
    ax.grid(axis="x")
    _ax_spine(ax)
    _save(fig, f"{outdir}/12_top_umis.png")
 
 
# ── 13. PCR duplication summary ───────────────────────────────────────────────
def plot_duplication_summary(umi_seqs: list, outdir: str):
    """FIX 4: function was called in main() but missing from file."""
    cnt         = Counter(umi_seqs)
    total       = len(umi_seqs)
    unique_umis = len(cnt)
    dup_reads   = sum(v - 1 for v in cnt.values())
    depth_cnt   = Counter(cnt.values())
 
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
 
    # pie: unique vs duplicate reads
    ax = axes[0]
    labels = ["Unique reads", "PCR duplicates"]
    sizes  = [total - dup_reads, dup_reads]
    colors = [PALETTE["accent3"], PALETTE["accent2"]]
    ax.pie(sizes, labels=labels, colors=colors, autopct="%1.1f%%", startangle=90,
           textprops={"color": PALETTE["text"]})
    ax.set_title(f"PCR Duplication Estimate\n"
                 f"Total: {total:,}  |  Unique UMIs: {unique_umis:,}")
 
    # depth distribution
    ax2 = axes[1]
    depths, freq = zip(*sorted(depth_cnt.items()))
    ax2.bar(depths, freq, color=PALETTE["accent1"], width=0.8)
    ax2.set_yscale("log")
    ax2.set_xlim(0, np.percentile(list(cnt.values()), 99))
    ax2.set_xlabel("Reads per unique UMI (depth)")
    ax2.set_ylabel("Count (log scale)")
    ax2.set_title("UMI Depth Distribution")
    ax2.grid(axis="y")
    _ax_spine(ax2)
 
    fig.tight_layout()
    _save(fig, f"{outdir}/13_duplication_summary.png")
 
 
# ── 00. Summary dashboard ─────────────────────────────────────────────────────
def plot_dashboard(metrics: dict, outdir: str):
    """
    Extended single-page QC summary dashboard.
 
    Metric boxes (row 1): total reads, reads passing, unique UMIs, duplication rate
    Metric boxes (row 2): median UMI Q, % ≥ Q30, mean GC%, % UMIs with N
    Metric boxes (row 3): median entropy, min entropy, % homopolymer, top UMI %
    Flags section       : auto-generated QC warnings below the boxes
    """
    m   = metrics
    fig = plt.figure(figsize=(14, 8))
    fig.patch.set_facecolor(PALETTE["bg"])
    ax  = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_facecolor(PALETTE["bg"])
 
    # ── title ────────────────────────────────────────────────────────────────
    ax.text(0.5, 0.97, "UMI QC Summary Dashboard",
            ha="center", va="top", fontsize=17, fontweight="bold",
            color=PALETTE["accent1"], transform=ax.transAxes)
    ax.text(0.5, 0.925, m.get("fastq_path", ""),
            ha="center", va="top", fontsize=8,
            color=PALETTE["text"], transform=ax.transAxes)
 
    # ── metric boxes (3 rows × 4 cols) ───────────────────────────────────────
    boxes = [
        # row 1 — read counts
        ("Total reads",        f"{m['total_reads']:,}",                 PALETTE["accent1"]),
        ("Reads passing",      f"{m['reads_passing']:,}",               PALETTE["accent1"]),
        ("% reads passing",    f"{m['pct_passing']:.1f}%",             PALETTE["accent3"]),
        ("Unique UMIs",        f"{m['unique_umis']:,}",                 PALETTE["accent3"]),
        # row 2 — quality
        ("Duplication rate",   f"{m['dup_rate'] * 100:.1f}%",          PALETTE["accent2"]),
        ("Median UMI Q",       f"Q{m['median_umi_q']:.1f}",            PALETTE["accent3"]),
        ("% UMI ≥ Q30",        f"{m['pct_q30']:.1f}%",                 PALETTE["accent3"]),
        ("Mean GC%",           f"{m['mean_gc']:.1f}%",                  PALETTE["accent4"]),
        # row 3 — diversity & bias
        ("% UMIs with N",      f"{m['pct_n']:.2f}%",                   PALETTE["accent2"]),
        ("Median entropy",     f"{m['median_entropy']:.2f} bits",       PALETTE["accent1"]),
        ("Min entropy",        f"{m['min_entropy']:.2f} bits",          PALETTE["accent4"]),
        ("% homopolymer UMIs", f"{m['pct_homopolymer']:.1f}%",         PALETTE["accent2"]),
    ]
 
    cols = 4
    w, h = 0.215, 0.12
    x0, y0 = 0.025, 0.83
 
    for i, (label, value, color) in enumerate(boxes):
        col = i % cols
        row = i // cols
        x   = x0 + col * (w + 0.018)
        y   = y0 - row * (h + 0.018)
        rect = plt.Rectangle((x, y), w, h, transform=ax.transAxes,
                              facecolor="#f5f5f5",
                              edgecolor=color, linewidth=1.5)
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h * 0.63, value,
                ha="center", va="center", fontsize=13, fontweight="bold",
                color=color, transform=ax.transAxes)
        ax.text(x + w / 2, y + h * 0.22, label,
                ha="center", va="center", fontsize=7.5,
                color="#444444", transform=ax.transAxes)
 
    # ── saturation slope bar ──────────────────────────────────────────────────
    slope      = m.get("saturation_slope", None)
    bar_y      = y0 - 3 * (h + 0.018) - 0.01
    bar_label_x = x0
    if slope is not None:
        slope_color = (PALETTE["accent3"] if slope < 0.2 else
                       PALETTE["accent4"] if slope < 0.5 else
                       PALETTE["accent2"])
        bar_w = 4 * w + 3 * 0.018          # full width of 4-col grid
        # background track
        bg = plt.Rectangle((bar_label_x, bar_y), bar_w, 0.035,
                            transform=ax.transAxes,
                            facecolor="#eeeeee", edgecolor=PALETTE["grid"], linewidth=0.8)
        ax.add_patch(bg)
        # filled bar proportional to slope
        filled = plt.Rectangle((bar_label_x, bar_y), bar_w * min(slope, 1.0), 0.035,
                                transform=ax.transAxes,
                                facecolor=slope_color, alpha=0.7)
        ax.add_patch(filled)
        slope_interp = ("over-sequenced" if slope < 0.1 else
                        "well saturated" if slope < 0.3 else
                        "moderate" if slope < 0.6 else
                        "under-sequenced")
        ax.text(bar_label_x + bar_w / 2, bar_y + 0.018,
                f"Saturation slope = {slope:.2f}  ({slope_interp})",
                ha="center", va="center", fontsize=8.5,
                color="#111111", transform=ax.transAxes)
 
    # ── QC flags ──────────────────────────────────────────────────────────────
    flags = []
    if m["pct_passing"] < 90:
        flags.append(f"⚠  Only {m['pct_passing']:.1f}% reads passed quality filter")
    if m["dup_rate"] > 0.5:
        flags.append(f"⚠  High duplication rate ({m['dup_rate'] * 100:.1f}%) — over-sequenced or low complexity")
    if m["pct_q30"] < 80:
        flags.append(f"⚠  Only {m['pct_q30']:.1f}% UMIs pass Q30 — check sequencing quality")
    if m["pct_n"] > 1.0:
        flags.append(f"⚠  {m['pct_n']:.2f}% UMIs contain N bases — check base calling")
    if m["min_entropy"] < 1.0:
        flags.append(f"⚠  Lowest positional entropy = {m['min_entropy']:.2f} bits — severe bias at that position")
    elif m["median_entropy"] < 1.5:
        flags.append(f"⚠  Median entropy {m['median_entropy']:.2f} bits — possible UMI synthesis bias")
    if m["mean_gc"] < 30 or m["mean_gc"] > 70:
        flags.append(f"⚠  Extreme GC content ({m['mean_gc']:.1f}%) — check synthesis protocol")
    if m["pct_homopolymer"] > 5:
        flags.append(f"⚠  {m['pct_homopolymer']:.1f}% UMIs have homopolymer run ≥ 50% of length")
    if slope is not None and slope > 0.7:
        flags.append(f"⚠  Saturation slope {slope:.2f} — library is under-sequenced")
    if not flags:
        flags = ["✓  No major QC flags detected"]
 
    flag_y = bar_y - 0.045
    ax.text(x0, flag_y + 0.03, "QC Flags",
            fontsize=10, fontweight="bold",
            color=PALETTE["accent4"], transform=ax.transAxes)
    for fl in flags:
        color = PALETTE["accent2"] if fl.startswith("⚠") else PALETTE["accent3"]
        ax.text(x0, flag_y, fl, fontsize=8.5, color=color, transform=ax.transAxes)
        flag_y -= 0.045
 
    _save(fig, f"{outdir}/00_summary_dashboard.png")
 
 
# =============================================================================
# TSV export
# =============================================================================
 
def export_tsv(data: dict, metrics: dict, outdir: str):
    # Per-base stats
    positions = sorted(data["per_base_nuc"].keys())
    rows = []
    for pos in positions:
        cnt   = data["per_base_nuc"][pos]
        total = sum(cnt.values())
        row   = {"position": pos}
        for b in "ACGTN":
            row[f"freq_{b}"] = cnt.get(b, 0) / total
        row["mean_qual"] = data["per_base_qual"].get(pos, np.nan)
        rows.append(row)
    pd.DataFrame(rows).to_csv(f"{outdir}/per_base_stats.tsv", index=False, sep="\t")
 
    # Scalar summary
    pd.DataFrame([metrics]).to_csv(f"{outdir}/umi_qc_summary.tsv", index=False, sep="\t")
 
    # UMI count table (top 10 000)
    cnt = Counter(data["umi_seqs"])
    df  = pd.DataFrame(cnt.most_common(10_000), columns=["umi", "count"])
    df["fraction"] = df["count"] / len(data["umi_seqs"])
    df.to_csv(f"{outdir}/umi_counts_top10k.tsv", index=False, sep="\t")
 
    log.info("  TSV tables written to %s", outdir)
 
 
# =============================================================================
# CLI
# =============================================================================
 
def parse_args():
    p = argparse.ArgumentParser(
        description="Comprehensive UMI quality control for gRNA FASTQ files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--fastq",      required=True,  help="Input FASTQ (.fastq or .fastq.gz)")
    p.add_argument("--outdir",     required=True,  help="Output directory for plots + TSVs")
    p.add_argument("--umi_len",    type=int, required=True, help="UMI length (nt)")
    p.add_argument("--umi_pos",    type=int, default=5, choices=[5, 3],
                   help="UMI at 5' (5) or 3' (3) end")
    p.add_argument("--subsample",  type=int, default=200_000,
                   help="Reservoir size for Hamming / saturation analyses")
    p.add_argument("--min_qual",   type=float, default=0,
                   help="Minimum mean UMI Phred to include a read")
    p.add_argument("--threads",    type=int, default=4, help="(reserved)")
    p.add_argument("--seed",       type=int, default=920)
    return p.parse_args()
 
 
# =============================================================================
# Main
# =============================================================================
 
def main():
    args   = parse_args()
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
 
    log.info("UMI QC  |  umi_len=%d  umi_end=%d'  min_qual=%s",
             args.umi_len, args.umi_pos, args.min_qual)
 
    # ── collect ───────────────────────────────────────────────────────────────
    data = collect_umi_data(          # FIX 2: was collect_umi_data → now defined as such
        fastq_path = args.fastq,
        umi_len    = args.umi_len,
        umi_end    = args.umi_pos,
        subsample  = args.subsample,
        min_qual   = args.min_qual,
        seed       = args.seed,
    )
 
    umi_seqs = data["umi_seqs"]
    if not umi_seqs:
        log.error("No UMIs extracted — check --umi_len / --umi_pos / --fastq")
        sys.exit(1)
 
    # ── derived metrics ───────────────────────────────────────────────────────
    log.info("Computing derived metrics …")
    entropy     = positional_entropy(data["per_base_nuc"])
    umi_arr     = np.array(data["umi_mean_qual"])
    cnt         = Counter(umi_seqs)
    total       = len(umi_seqs)
    unique_umis = len(cnt)
    dup_reads   = sum(v - 1 for v in cnt.values())
    dup_rate    = dup_reads / total if total else 0
 
    # saturation slope (use reservoir to stay memory-safe)
    reservoir_seqs = [r[0] for r in data["reservoir"]]   # FIX 7: extract str from tuple
    sat_sizes, sat_unique = saturation_curve(reservoir_seqs)
    saturation_slope = (sat_unique[-1] / sat_sizes[-1]
                        if sat_sizes[-1] > 0 else 0.0)
 
    # Hamming distances
    ham_dists = None
    try:
        from scipy.stats import binom as _binom  # noqa – confirm scipy available
        log.info("Sampling Hamming distances …")
        ham_dists = sample_hamming(reservoir_seqs, n_pairs=50_000, seed=args.seed)
    except ImportError:
        log.warning("scipy not installed — skipping Hamming plot")
 
    # homopolymer rate
    hp_threshold   = args.umi_len * 0.5
    pct_homopolymer = sum(1 for h in data["homopolymer_lens"] if h >= hp_threshold) / total * 100
 
    # top UMI frequency
    top_umi_freq = cnt.most_common(1)[0][1] / total * 100 if cnt else 0.0
 
    # scalar metrics dict — extended for dashboard
    metrics = {
        "fastq_path":        args.fastq,
        "total_reads":       data["total_reads"],
        "reads_passing":     total,
        "pct_passing":       round(total / data["total_reads"] * 100, 2) if data["total_reads"] else 0,
        "unique_umis":       unique_umis,
        "dup_rate":          round(dup_rate, 4),
        "median_umi_q":      round(float(np.median(umi_arr)), 2),
        "pct_q30":           round(float((umi_arr >= 30).mean() * 100), 2),
        "mean_gc":           round(float(np.mean(data["gc"]) * 100), 2),
        "pct_n":             round(float(np.mean([x > 0 for x in data["n_content"]]) * 100), 2),
        "median_entropy":    round(float(np.median(list(entropy.values()))), 4),
        "min_entropy":       round(float(min(entropy.values())), 4),
        "pct_homopolymer":   round(pct_homopolymer, 2),
        "top_umi_pct":       round(top_umi_freq, 4),
        "saturation_slope":  round(saturation_slope, 4),
    }
 
    # ── plots ─────────────────────────────────────────────────────────────────
    log.info("Generating plots → %s", outdir)
 
    plot_dashboard(metrics, outdir)
    plot_length_dist(data["umi_lengths"], args.umi_len, outdir)
    plot_per_base_composition(data["per_base_nuc"], outdir)
    plot_per_base_qual(data["per_base_qual"], outdir)
    plot_gc_content(data["gc"], outdir)
    plot_homopolymer(data["homopolymer_lens"], args.umi_len, outdir)
    plot_umi_qual_hist(data["umi_mean_qual"], outdir)
    plot_saturation(reservoir_seqs, outdir)
 
    if ham_dists:
        plot_hamming(ham_dists, args.umi_len, outdir)
 
    plot_entropy(entropy, outdir)
    plot_qual_scatter(data["read_mean_qual"], data["umi_mean_qual"], outdir)  # FIX 6
    plot_n_content(data["n_content"], outdir)
    plot_top_umis(umi_seqs, outdir)
    plot_duplication_summary(umi_seqs, outdir)   # FIX 4
 
    # ── TSV export ────────────────────────────────────────────────────────────
    log.info("Exporting TSV tables …")
    export_tsv(data, metrics, outdir)
 
    # ── console summary ───────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("  UMI QC Summary")
    print("=" * 60)
    for k, v in metrics.items():
        if k == "fastq_path":
            continue
        print(f"  {k:<24} {v}")
    print("=" * 60)
    print(f"  Plots saved to: {outdir}/")
    print("=" * 60 + "\n")
 
 
if __name__ == "__main__":
    main()