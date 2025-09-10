#!/usr/bin/env python3
import argparse
import gzip
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")   # headless
import matplotlib.pyplot as plt
from collections import Counter, defaultdict

def smart_open(p): 
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")

# ---------- Loaders ----------
def load_bins(counts_path, chrom_filter=None):
    xs = []  # (chrom, start, end, x1_ESC, x2_NPC)
    with smart_open(counts_path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("#") or ln.startswith("-"):
                continue
            if ln.lower().startswith("chrom"):
                continue
            parts = ln.split("\t")
            if len(parts) < 5:
                continue
            chrom, s, e = parts[0], int(parts[1]), int(parts[2])
            x1, x2 = int(parts[3]), int(parts[4])
            if chrom_filter and chrom != chrom_filter:
                continue
            xs.append((chrom, s, e, x1, x2))
    return xs

def load_dhms_bins(path):
    # BED6-ish with state in col6: A1 (ESC-enriched) or A2 (NPC-enriched)
    by_chr = defaultdict(list)
    a1 = a2 = 0
    with smart_open(path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom, s, e, _, _, state = parts[:6]
            s, e = int(s), int(e)
            by_chr[chrom].append((s, e, state))
            if state == "A1":
                a1 += 1
            elif state == "A2":
                a2 += 1
    for c in by_chr:
        by_chr[c].sort(key=lambda t: t[0])
    return by_chr, a1, a2

def load_dhms_regions(path):
    lens = []
    per_chr = Counter()
    with smart_open(path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, s, e = parts[0], int(parts[1]), int(parts[2])
            L = max(0, e - s)
            lens.append(L)
            per_chr[chrom] += 1
    return lens, per_chr

# ---------- Helpers ----------
def region_series(bins, chrom, start, end):
    pos, xesc, xnpc = [], [], []
    for c, s, e, x1, x2 in bins:
        if c != chrom:
            continue
        if e <= start:
            continue
        if s >= end:
            break
        pos.append((s + e) // 2)
        xesc.append(x1)
        xnpc.append(x2)
    return np.array(pos), np.array(xesc), np.array(xnpc)

def safe_log2_ratio(num, den, eps=1.0):
    return np.log2((num + eps) / (den + eps))

# ---------- Figures ----------
def fig_tracks(counts_path, region_chr, region_start, region_end, out_png):
    allbins = load_bins(counts_path, chrom_filter=region_chr)
    pos, xesc, xnpc = region_series(allbins, region_chr, region_start, region_end)
    if len(pos) == 0:
        print(f"[warn] no bins in region {region_chr}:{region_start}-{region_end}")
        return
    # Flip to ESC/NPC as requested
    lr = safe_log2_ratio(xesc, xnpc)

    plt.figure(figsize=(10, 6))
    ax1 = plt.subplot(3,1,1)
    ax1.fill_between(pos, xesc, step="mid", alpha=0.9)
    ax1.set_ylabel("ESC\ncounts")
    ax1.set_xlim(region_start, region_end)

    ax2 = plt.subplot(3,1,2, sharex=ax1)
    ax2.fill_between(pos, xnpc, step="mid", alpha=0.9)
    ax2.set_ylabel("NPC\ncounts")

    ax3 = plt.subplot(3,1,3, sharex=ax1)
    ax3.plot(pos, lr, linewidth=1)
    ax3.axhline(0, color="black", linewidth=0.8)
    ax3.set_ylabel("log2(ESC/NPC)")
    ax3.set_xlabel(f"{region_chr}:{region_start}-{region_end}")

    plt.tight_layout()
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    plt.savefig(out_png, dpi=150)
    plt.close()
    print(f"[ok] region tracks saved -> {out_png}")

def fig_tracks_with_calls(counts_path, dhms_bins_path, region_chr, region_start, region_end, out_png):
    allbins = load_bins(counts_path, chrom_filter=region_chr)
    pos, xesc, xnpc = region_series(allbins, region_chr, region_start, region_end)
    if len(pos) == 0:
        print(f"[warn] no bins in region {region_chr}:{region_start}-{region_end}")
        return
    lr = safe_log2_ratio(xesc, xnpc)  # ESC/NPC

    calls_by_chr, _, _ = load_dhms_bins(dhms_bins_path)
    calls = [t for t in calls_by_chr.get(region_chr, []) if not (t[1] <= region_start or t[0] >= region_end)]

    plt.figure(figsize=(10, 7))
    ax1 = plt.subplot(4,1,1)
    ax1.fill_between(pos, xesc, step="mid", alpha=0.9)
    ax1.set_ylabel("ESC")
    ax1.set_xlim(region_start, region_end)

    ax2 = plt.subplot(4,1,2, sharex=ax1)
    ax2.fill_between(pos, xnpc, step="mid", alpha=0.9)
    ax2.set_ylabel("NPC")

    ax3 = plt.subplot(4,1,3, sharex=ax1)
    ax3.plot(pos, lr, linewidth=1)
    ax3.axhline(0, color="black", linewidth=0.8)
    ax3.set_ylabel("log2(ESC/NPC)")

    ax4 = plt.subplot(4,1,4, sharex=ax1)
    ax4.set_ylim(0, 1)
    ax4.set_yticks([])
    ax4.set_ylabel("DHMS calls")
    ax4.set_xlabel(f"{region_chr}:{region_start}-{region_end}")
    for s, e, state in calls:
        color = "#377eb8" if state == "A1" else "#e41a1c"  # blue ESC-enriched, red NPC-enriched
        ax4.axvspan(max(s, region_start), min(e, region_end), color=color, alpha=0.35)

    plt.tight_layout()
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    plt.savefig(out_png, dpi=150)
    plt.close()
    print(f"[ok] region w/ calls saved -> {out_png}")

def fig_ma(counts_path, ri_chr, out_png, sample_n=300000):
    xs = load_bins(counts_path, chrom_filter=None if ri_chr=="all" else ri_chr)
    if not xs:
        print(f"[warn] no bins for MA on {ri_chr}")
        return
    xesc = np.array([t[3] for t in xs], dtype=float)
    xnpc = np.array([t[4] for t in xs], dtype=float)
    A = np.log2(xesc + xnpc + 1.0)           # log-intensity
    M = np.log2((xesc + 1.0)/(xnpc + 1.0))   # ESC/NPC

    n = len(A)
    if n > sample_n:
        idx = np.random.RandomState(1).choice(n, size=sample_n, replace=False)
        A = A[idx]
        M = M[idx]

    plt.figure(figsize=(6,5))
    plt.scatter(A, M, s=2, alpha=0.5)
    plt.axhline(0, linewidth=0.8, color="black")
    plt.xlabel("log2 intensity")
    plt.ylabel("log2 ESC/NPC")
    plt.title(f"MA plot ({ri_chr})")
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[ok] MA plot saved -> {out_png}")

def fig_hist(counts_path, out_png):
    xs = load_bins(counts_path)
    if not xs:
        print("[warn] no bins for histogram")
        return
    xesc = np.array([t[3] for t in xs], dtype=float)
    xnpc = np.array([t[4] for t in xs], dtype=float)
    lr = np.log2((xesc+1.0)/(xnpc+1.0))  # ESC/NPC
    plt.figure(figsize=(6,4))
    plt.hist(lr, bins=101)
    plt.xlabel("log2 ESC/NPC per 1kb bin")
    plt.ylabel("bin count")
    plt.title("Genome-wide distribution (ESC vs NPC)")
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[ok] histogram saved -> {out_png}")

def fig_dhms_state_bar(dhms_bins_path, out_png):
    _, a1, a2 = load_dhms_bins(dhms_bins_path)
    if a1 + a2 == 0:
        print("[warn] no DHMS bins for state bar")
        return
    plt.figure(figsize=(5,4))
    bars = plt.bar(["ESC-enriched\n(A1)", "NPC-enriched\n(A2)"], [a1, a2])
    plt.ylabel("DHMS bins")
    plt.title("DHMS calls by state")
    for b, n in zip(bars, [a1, a2]):
        plt.text(b.get_x()+b.get_width()/2, n, str(n), ha="center", va="bottom")
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[ok] DHMS state bar saved -> {out_png}")

def fig_dhms_len_hist(dhms_regions_path, out_png):
    lens, _ = load_dhms_regions(dhms_regions_path)
    if not lens:
        print("[warn] no DHMS regions found for length histogram")
        return
    kb = np.array(lens, dtype=float) / 1000.0
    plt.figure(figsize=(6,4))
    plt.hist(kb, bins=50)
    plt.xlabel("DHMS region length (kb)")
    plt.ylabel("# regions")
    plt.title("DHMS region length distribution")
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[ok] DHMS length hist saved -> {out_png}")

def fig_dhms_per_chrom(dhms_regions_path, out_png):
    _, per_chr = load_dhms_regions(dhms_regions_path)
    if not per_chr:
        print("[warn] no DHMS regions for per-chrom chart")
        return
    items = sorted(per_chr.items(), key=lambda x: x[0])
    labels = [k for k,_ in items]
    counts = [v for _,v in items]
    plt.figure(figsize=(10,4))
    plt.bar(labels, counts)
    plt.ylabel("# DHMS regions")
    plt.title("DHMS regions per chromosome")
    plt.xticks(rotation=90)
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    plt.savefig(out_png, dpi=150)
    plt.close()
    print(f"[ok] DHMS per-chrom saved -> {out_png}")

# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(description="Quick visuals for ChIPDiff-faithful outputs (simple, result-focused).")
    ap.add_argument("--counts", required=True, help="results/bin_counts.tsv.gz")
    ap.add_argument("--dhms-bins", default="chipdiff-faithful/results/dhms_bins.bed.gz")
    ap.add_argument("--dhms-regions", default="chipdiff-faithful/results/dhms_regions.bed.gz")
    ap.add_argument("--region-chr", default="chr14")
    ap.add_argument("--region-start", type=int, default=117100000)
    ap.add_argument("--region-end",   type=int, default=117130000)
    ap.add_argument("--ri-chr", default="all", help="'all' or a chromosome (e.g. chr19)")
    ap.add_argument("--out-dir", default="chipdiff-faithful/results/plots")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Raw bin context
    fig_tracks(args.counts, args.region_chr, args.region_start, args.region_end,
               os.path.join(args.out_dir, "region_tracks_ESCoverNPC.png"))
    fig_tracks_with_calls(args.counts, args.dhms_bins, args.region_chr, args.region_start, args.region_end,
               os.path.join(args.out_dir, "region_tracks_with_calls.png"))
    fig_ma(args.counts, args.ri_chr,
           os.path.join(args.out_dir, f"ma_{args.ri_chr}_ESCoverNPC.png"))
    fig_hist(args.counts, os.path.join(args.out_dir, "logratio_hist_ESCoverNPC.png"))

    # HMM result-centric
    if os.path.exists(args.dhms_bins):
        fig_dhms_state_bar(args.dhms_bins, os.path.join(args.out_dir, "dhms_state_bar.png"))
    if os.path.exists(args.dhms_regions):
        fig_dhms_len_hist(args.dhms_regions, os.path.join(args.out_dir, "dhms_length_hist.png"))
        fig_dhms_per_chrom(args.dhms_regions, os.path.join(args.out_dir, "dhms_per_chrom.png"))

if __name__ == "__main__":
    main()
