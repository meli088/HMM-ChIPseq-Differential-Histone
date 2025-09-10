#!/usr/bin/env python3
import argparse
import gzip
import os
import re
from collections import defaultdict, namedtuple

def smart_open(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")

def parse_gtf_for_promoters(gtf_path, pad=1000):
    """
    Parse a GTF file and extract promoter regions (TSS-centered windows) for each gene.

    The function scans a GTF file (opened via the provided gtf_path with smart_open)
    and uses rows with feature == "transcript" to identify transcription start site
    (TSS) candidates. For each gene it collects all TSS candidates (one per transcript)
    and then chooses a single representative TSS per gene:
    - for genes on the "+" strand: the minimum transcript start is used as TSS
    - for genes on the "-" strand: the maximum transcript end is used as TSS

    A promoter interval is produced by padding the selected TSS by `pad` bases on
    both sides, producing the interval [tss - pad, tss + pad). Promoter start
    coordinates are clipped to be >= 0.

    Attributes and name resolution
    - The attribute field is parsed with a simple regex to extract key/"value" pairs.
    - The function prefers the "gene_name" attribute if present; otherwise it uses
        "gene_id"; if neither is present it uses the string "NA" as the gene key.

    Filtering and conflict resolution
    - Only lines with feature == "transcript" are considered (to get robust TSS).
    - Only transcripts with strand in ("+", "-") are used.
    - If the same gene is encountered on multiple chromosomes or strands, the first
        encountered chrom/strand pair is kept and subsequent inconsistent entries are
        skipped.

    Return value
    - Returns a collections.defaultdict(list) mapping chromosome name (str) to a
        list of promoter tuples sorted by promoter start:
            (promoter_start, promoter_end, gene_name_or_id, strand)
        promoter_start and promoter_end are integers in 0-based coordinates (start is
        inclusive, end is exclusive as usual in Python conventions used here).

    Parameters
    - gtf_path (str or file-like): Path to the GTF file (opened using smart_open).
    - pad (int, optional): Number of bases to pad upstream and downstream of the
        selected TSS to form the promoter window. Default is 1000.

    Notes and assumptions
    - The GTF start coordinate is converted from 1-based inclusive to 0-based
        half-open by subtracting 1 from the reported start.
    - The attribute parser uses a simple regex r'(\w+)\s+"([^"]+)"' and may not
        capture all valid/complex attribute formats in every GTF variant.
    - The function requires that smart_open and re are available in the calling
        module's scope.
    - Promoters are returned per chromosome and sorted by start position.
    """
    # gene -> {chrom, strand, tss_candidates(list)}
    G = {}
    # Prefer 'gene_name' if present, else 'gene_id'
    def parse_attrs(s):
        d = {}
        for m in re.finditer(r'(\w+)\s+"([^"]+)"', s):
            d[m.group(1)] = m.group(2)
        return d

    with smart_open(gtf_path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature != "transcript":  # use transcript rows for robust TSS
                continue
            start = int(start) - 1  # GTF is 1-based inclusive; convert to 0-based half-open
            end   = int(end)
            a = parse_attrs(attrs)
            gene = a.get("gene_name") or a.get("gene_id") or "NA"
            if strand not in ("+", "-"):
                continue
            rec = G.setdefault(gene, {"chrom":chrom, "strand":strand, "tss": []})
            # If same gene appears on multiple chrom/strand (rare), keep first consistent set
            if rec["chrom"] != chrom or rec["strand"] != strand:
                # skip conflicting entries
                continue
            # TSS candidate
            tss = start if strand=="+" else end
            rec["tss"].append(tss)

    promoters_by_chrom = defaultdict(list)
    for gene, rec in G.items():
        chrom, strand, tss_list = rec["chrom"], rec["strand"], rec["tss"]
        if not tss_list:
            continue
        tss = min(tss_list) if strand=="+" else max(tss_list)
        ps = max(0, tss - pad)
        pe = tss + pad
        promoters_by_chrom[chrom].append((ps, pe, gene, strand))

    for c in promoters_by_chrom:
        promoters_by_chrom[c].sort(key=lambda x: x[0])
    return promoters_by_chrom

Bin = namedtuple("Bin", "start end x1 x2")

def load_bins(counts_path):
    """
    Load genomic bins and their counts from a tab-delimited counts file.

    Each non-empty, non-comment line is expected to contain at least five columns:
    chrom, start, end, x1, x2. The start, end, x1 and x2 columns are parsed as integers.
    Lines that are empty, begin with '-' or '#', begin with a header starting with
    'chrom' (case-insensitive), contain fewer than five columns, or contain values
    that cannot be converted to integers are silently skipped.

    Parameters
    ----------
    counts_path : str or os.PathLike
        Path to the counts file. The function uses `smart_open` to open the file,
        so compressed files and other openable resources supported by `smart_open`
        are accepted.

    Returns
    -------
    defaultdict(list)
        A mapping from chromosome name (str) to a list of Bin objects for that
        chromosome. For each chromosome the list is sorted in ascending order by
        the Bin.start attribute.

    Behavior / Notes
    ----------------
    - Each accepted input line is converted to a Bin(s, e, x1, x2) and appended
      to the list for the corresponding chromosome key.
    - The function performs best-effort parsing and does not raise on malformed
      lines; such lines are skipped.
    - The returned container is a collections.defaultdict(list), so missing keys
      will produce empty lists when accessed.

    Example
    -------
    >>> bins = load_bins('counts.tsv')
    >>> chrom_bins = bins['chr1']  # list of Bin objects sorted by start
    """
    bins = defaultdict(list)
    with smart_open(counts_path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            if ln.startswith("-") or ln.startswith("#"):
                continue
            # skip the tab header
            if ln.lower().startswith("chrom"):
                continue
            parts = ln.split("\t")
            if len(parts) < 5:
                continue
            try:
                chrom = parts[0]
                s = int(parts[1])
                e = int(parts[2])
                x1 = int(parts[3])
                x2 = int(parts[4])
            except ValueError:
                # any unexpected line, just skip
                continue
            bins[chrom].append(Bin(s, e, x1, x2))
    for c in bins:
        bins[c].sort(key=lambda b: b.start)
    return bins


def sum_counts_overlaps(promoters, bins):
    """
    Compute summed counts of overlapping bins for promoter regions.

    Parameters
    ----------
    promoters : dict
        Mapping from chromosome name (str) to a sequence of promoter tuples.
        Each promoter tuple must be (start, end, gene, strand) where:
          - start, end : int
              0-based (or 1-based) genomic coordinates defining the promoter interval.
          - gene : hashable
              Gene identifier used as a key for per-gene aggregation.
          - strand : object
              Strand information (kept but not used for aggregation).
        Promoter intervals for each chromosome are assumed to be sorted by start
        coordinate (non-decreasing). Overlapping or adjacent promoter intervals are
        permitted but affect per-gene counting according to input tuples.

    bins : dict
        Mapping from chromosome name (str) to a sequence (e.g. list) of bin
        objects. Each bin object must expose attributes:
          - start, end : int
              Genomic coordinates of the bin interval.
          - x1, x2 : numeric
              Two count/value fields to be summed when a bin overlaps a promoter.
        Bin intervals for each chromosome are assumed to be sorted by start
        coordinate (non-decreasing). The implementation advances by the earlier
        interval end, so correct ordering is required for linear-time behavior.

    Returns
    -------
    tuple
        A pair (global_sums, per_gene) where:
          - global_sums : tuple (sum_x1, sum_x2, overlapped_bins)
              sum_x1 : numeric
                  Total sum of x1 across all bins that overlap any promoter.
              sum_x2 : numeric
                  Total sum of x2 across all bins that overlap any promoter.
              overlapped_bins : int
                  Total number of bin intervals that overlapped at least one
                  promoter (each overlapping bin is counted once per overlap
                  detection in the two-list sweep).
          - per_gene : collections.defaultdict(list)
              Mapping gene -> [gene_x1_sum, gene_x2_sum, gene_overlapped_bin_count].
              The defaultdict is populated with a 3-element list for each
              encountered gene: cumulative x1, cumulative x2, and count of bins
              observed overlapping that gene.

    Behavior and notes
    ------------------
    - The function uses a two-pointer sweep per chromosome to detect overlaps
      between promoter intervals and bin intervals. This yields linear-time
      behavior O(P + B) per chromosome assuming inputs are sorted.
    - When a bin overlaps a promoter, the bin's x1 and x2 are added to both the
      global sums and to the corresponding gene's cumulative values, and the
      overlapped_bins and per-gene count are incremented.
    - If a bin overlaps multiple promoters on the same chromosome, it will be
      processed separately for each overlapping promoter as determined by the
      sweep (i.e., counts may be added multiple times if promoter intervals
      overlap), because aggregation keys are per promoter's gene identifier.
    - Coordinates and count types should be numeric; the function does not
      validate types and will raise on missing attributes or incompatible types.
    - If a chromosome appears in only one of the inputs, it is treated as having
      zero intervals in the other mapping (no overlaps will be counted for that
      chromosome).
    - The function has no side effects other than populating and returning the
      per_gene defaultdict; it does not modify the input mappings.
    """
    sum_x1 = sum_x2 = 0
    overlapped_bins = 0
    per_gene = defaultdict(lambda: [0,0,0])
    for chrom in promoters.keys() | bins.keys():
        P = promoters.get(chrom, [])
        B = bins.get(chrom, [])
        i = j = 0
        while i < len(P) and j < len(B):
            ps, pe, gene, strand = P[i]
            bs, be, x1, x2 = B[j].start, B[j].end, B[j].x1, B[j].x2
            if be <= ps:
                j += 1
            elif pe <= bs:
                i += 1
            else:
                # overlap
                sum_x1 += x1
                sum_x2 += x2
                overlapped_bins += 1
                g = per_gene[gene]
                g[0] += x1
                g[1] += x2
                g[2] += 1
                # advance the earlier end
                if be < pe:
                    j += 1
                else:
                    i += 1
    return (sum_x1, sum_x2, overlapped_bins), per_gene

def write_promoters_bed(promoters_by_chrom, out_path):
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as fo:
        for chrom in sorted(promoters_by_chrom):
            for ps, pe, gene, strand in promoters_by_chrom[chrom]:
                fo.write(f"{chrom}\t{ps}\t{pe}\t{gene}\t0\t{strand}\n")

def main():
    ap = argparse.ArgumentParser(
        description="Build promoters (±1 kb) from mm9.refGene.gtf(.gz) and test if NPC has less H3K27me3 at promoters than ESC using bin_counts."
    )
    ap.add_argument("--gtf", required=True, help="refs/mm9.refGene.gtf.gz")
    ap.add_argument("--counts", required=True, help="results/bin_counts.tsv.gz")
    ap.add_argument("--pad", type=int, default=1000, help="half window around TSS (default 1000)")
    ap.add_argument("--promoters-out", default="results/mm9_promoters_pm1kb.bed",
                    help="write promoters BED (default results/mm9_promoters_pm1kb.bed)")
    ap.add_argument("--per-gene-out", default="results/promoter_counts_per_gene.tsv",
                    help="per-gene summed x1/x2 counts across promoter-overlapping bins")
    args = ap.parse_args()

    promoters = parse_gtf_for_promoters(args.gtf, pad=args.pad)
    write_promoters_bed(promoters, args.promoters_out)

    bins = load_bins(args.counts)
    (sx1, sx2, nbin), per_gene = sum_counts_overlaps(promoters, bins)

    # Global summary: repression proxy = H3K27me3; expect NPC<ESC if repression is lost in NPC
    fc_global = (sx2 / sx1) if sx1 > 0 else float("nan")
    trend = "NPC < ESC (depletion in NPC)" if sx2 < sx1 else "NPC ≥ ESC"
    print("Promoter H3K27me3 summary (bin overlap aggregation):")
    print(f"  promoters written:     {sum(len(v) for v in promoters.values())}")
    print(f"  promoter-overlap bins: {nbin}")
    print(f"  sum counts ESC (x1):   {sx1}")
    print(f"  sum counts NPC (x2):   {sx2}")
    print(f"  NPC/ESC ratio:         {fc_global:.4f}  --> {trend}")

    # Per-gene table (optional inspection)
    os.makedirs(os.path.dirname(args.per_gene_out) or ".", exist_ok=True)
    with open(args.per_gene_out, "w") as fo:
        fo.write("gene\tn_bins\tsum_x1_ESC\tsum_x2_NPC\tratio_x2_over_x1\n")
        for g, (gx1, gx2, gn) in sorted(per_gene.items()):
            ratio = (gx2 / gx1) if gx1 > 0 else "NA"
            fo.write(f"{g}\t{gn}\t{gx1}\t{gx2}\t{ratio}\n")

if __name__ == "__main__":
    main()
