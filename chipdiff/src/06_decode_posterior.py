#!/usr/bin/env python3
import sys
import os
import gzip
import argparse
import numpy as np
from collections import defaultdict

def smart_open(p, mode="rt"):
    return gzip.open(p, mode) if p.endswith(".gz") else open(p, mode)

def parse_header_counts(path):
    with smart_open(path, "rt") as f:
        head = f.readline()
        if not head.startswith("#"):
            raise SystemExit("Header manquant dans bin_counts.")
        toks = head.strip("# \n").replace("\t"," ").split()
        vals = {}
        for t in toks:
            if "=" in t:
                k,v = t.split("=",1)
                vals[k]=int(v)
        for k in ("n1","n2","m","bin"):
            if k not in vals:
                raise SystemExit(f"{k} manquant dans header.")
        # skip optional column header
        _ = f.readline()
    return vals["n1"], vals["n2"], vals["m"], vals["bin"]

def load_bin_counts(path):
    """Map (chrom,bin_idx) -> (start,end,x1,x2) ; et bins par chrom."""
    n1,n2,m,bin_size = parse_header_counts(path)
    counts = {}
    chrom_bins = defaultdict(int)
    with smart_open(path,"rt") as f:
        _ = f.readline()
        for line in f:
            if not line or line.startswith("#"):
                continue
            if line.lower().startswith("chrom"):
                continue
            chrom,s_s,e_s,x1_s,x2_s = line.rstrip("\n").split("\t")
            s,e = int(s_s), int(e_s)
            x1,x2 = int(x1_s), int(x2_s)
            bidx = s // bin_size
            counts[(chrom,bidx)] = (s,e,x1,x2)
            if bidx+1 > chrom_bins[chrom]:
                chrom_bins[chrom] = bidx+1
    return counts, dict(chrom_bins), n1,n2,m,bin_size

def load_regions(path):
    regs=[]
    with smart_open(path,"rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            regs.append((parts[0], int(parts[1]), int(parts[2])))
    return regs

def load_lut(path):
    d = np.load(path, allow_pickle=True)
    return d["lut"], list(d["states"])

def load_transitions(path):
    d = np.load(path, allow_pickle=True)
    A = d["A"]          # (3,3) in prob
    A_log = np.log(A + 1e-300)
    return A, A_log

def build_sequences(counts, bin_size, regions):
    """Retourne liste de séquences: [(chrom, bidx, start, end, x1, x2), ...] par région."""
    seqs = []
    for chrom, start, end in regions:
        if end <= start:
            continue
        b0 = start // bin_size
        b1 = (end - 1) // bin_size
        seq=[]
        for b in range(b0, b1+1):
            rec = counts.get((chrom,b))
            if rec is None:
                continue
            s,e,x1,x2 = rec
            seq.append((chrom,b,s,e,x1,x2))
        if seq:
            seqs.append(seq)
    return seqs

def emission_logprobs_for_seq(seq, lut):
    T = len(seq)
    B = np.empty((T,3), dtype=np.float64)
    max_x1 = lut.shape[0] - 1
    max_x2 = lut.shape[1] - 1
    for t,(_,_,_,_,x1,x2) in enumerate(seq):
        i = x1 if x1<=max_x1 else max_x1
        j = x2 if x2<=max_x2 else max_x2
        B[t,:] = lut[i,j,:]
    return B

def logsumexp(a, axis=None):
    amax = np.max(a, axis=axis, keepdims=True)
    s = np.log(np.sum(np.exp(a-amax), axis=axis, keepdims=True) + 1e-300) + amax
    return np.squeeze(s, axis=axis)

def forward_backward_log(B, A_log, pi_log):
    T,K = B.shape
    la = np.empty((T,K))
    la[0,:] = pi_log + B[0,:]
    for t in range(1,T):
        la[t,:] = B[t,:] + logsumexp(la[t-1,:][:,None] + A_log, axis=0)
    loglik = logsumexp(la[-1,:], axis=0)

    lb = np.empty((T,K))
    lb[-1,:] = 0.0
    for t in range(T-2, -1, -1):
        lb[t,:] = logsumexp(A_log + (B[t+1,:] + lb[t+1,:])[None,:], axis=1)

    log_gamma = la + lb - loglik  # (T,K)
    gamma = np.exp(log_gamma)     # en proba
    return gamma, float(loglik)

def merge_dhms_bins(dhms_bins, gap_bp=0):
    """
    Merge overlapping or nearby genomic bins (DHMs) into consolidated intervals.

    Parameters
    ----------
    dhms_bins : Iterable[Tuple[str, int, int]]
        Sequence of bins where each element is a tuple (chromosome, start, end).
        Coordinates are treated as integers; the function assumes start <= end but
        does not enforce it.
    gap_bp : int, optional
        Maximum gap in base pairs between adjacent intervals on the same chromosome
        that will be bridged when merging. Default is 0 (only overlapping or
        directly touching intervals are merged).

    Returns
    -------
    List[Tuple[str, int, int]]
        A list of merged bins as (chromosome, start, end) tuples. The output is
        sorted by chromosome and start coordinate. If the input is empty, an empty
        list is returned.

    Behavior and notes
    ------------------
    - Input is sorted by (chromosome, start, end) before merging, so original order
      is not preserved.
    - Two intervals on the same chromosome are merged if the next interval's start
      is <= current_end + gap_bp.
    - Intervals on different chromosomes are never merged across chromosomes.
    - Time complexity is O(n log n) due to the initial sort; the merge pass is O(n).
    """
    if not dhms_bins:
        return []
    # trier
    dhms_bins.sort(key=lambda r: (r[0], r[1], r[2]))
    merged=[]
    c, s, e = dhms_bins[0][0], dhms_bins[0][1], dhms_bins[0][2]
    for chrom, start, end in dhms_bins[1:]:
        if chrom==c and start <= e + gap_bp:
            if end > e:
                e = end
        else:
            merged.append((c,s,e))
            c,s,e = chrom, start, end
    merged.append((c,s,e))
    return merged

def main():
    ap = argparse.ArgumentParser(
        description=("Décodage forward–backward avec A appris ; appel DHMS bin si P(A1)>rho ou P(A2)>rho ; "
                     "fusion des bins DHMS adjacents (gap=0). Papier-fidèle.")
    )
    ap.add_argument("--obs", required=True, help="results/bin_counts.tsv.gz")
    ap.add_argument("--putative", required=True, help="results/putative_regions.bed.gz")
    ap.add_argument("--lut", required=True, help="results/emission_lut.npz")
    ap.add_argument("--trans", required=True, help="results/transitions.npz")
    ap.add_argument("--rho", type=float, default=0.95, help="seuil postérieur (default 0.95)")
    ap.add_argument("--out-bins", default="results/dhms_bins.bed.gz")
    ap.add_argument("--out-regions", default="results/dhms_regions.bed.gz")
    args = ap.parse_args()

    counts, chrom_bins, n1,n2,m,bin_size = load_bin_counts(args.obs)
    regions = load_regions(args.putative)
    lut, _states = load_lut(args.lut)
    A, A_log = load_transitions(args.trans)

    # π = [1,0,0] en log
    pi_log = np.array([0.0, -np.inf, -np.inf], dtype=np.float64)

    os.makedirs(os.path.dirname(args.out_bins) or ".", exist_ok=True)
    outb = smart_open(args.out_bins, "wt")

    dhms_bin_list = []  # for merging later
    total_bins = called_bins = 0
    total_ll = 0.0

    seqs = build_sequences(counts, bin_size, regions)
    for seq in seqs:
        if not seq:
            continue
        B = emission_logprobs_for_seq(seq, lut)
        gamma, ll = forward_backward_log(B, A_log, pi_log)
        total_ll += ll

        # for each bin: write postA1, postA2 and state if exceeds rho
        for t, (chrom,b,s,e,x1,x2) in enumerate(seq):
            total_bins += 1
            pA1 = gamma[t,1]
            pA2 = gamma[t,2]
            state = "NONE"
            if pA1 > args.rho and pA1 >= pA2:
                state = "A1"
                dhms_bin_list.append((chrom, s, e))
                called_bins += 1
            elif pA2 > args.rho and pA2 > pA1:
                state = "A2"
                dhms_bin_list.append((chrom, s, e))
                called_bins += 1
            outb.write(f"{chrom}\t{s}\t{e}\t{pA1:.6g}\t{pA2:.6g}\t{state}\n")

    outb.close()

    # Merging (gap=0, "no gap")
    regions_merged = merge_dhms_bins(dhms_bin_list, gap_bp=0)
    os.makedirs(os.path.dirname(args.out_regions) or ".", exist_ok=True)
    with smart_open(args.out_regions, "wt") as fr:
        for chrom, s, e in regions_merged:
            fr.write(f"{chrom}\t{s}\t{e}\n")

    sys.stderr.write(f"[decode] total bins traversed: {total_bins}, called DHMS bins: {called_bins}\n")
    sys.stderr.write(f"[decode] merged DHMS regions: {len(regions_merged)}\n")
    sys.stderr.write(f"[decode] total loglik across regions: {total_ll:.3f}\n")

if __name__ == "__main__":
    main()
