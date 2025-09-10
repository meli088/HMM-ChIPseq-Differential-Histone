#!/usr/bin/env python3
import sys
import os
import gzip
import argparse

def smart_open(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def parse_header_and_prepare(fp):
    """
    Parse the custom header at the start of a data file and position the file
    pointer ready to read data rows.

    Parameters
    ----------
    fp : file-like object
        Open text-mode file object (or file-like) that supports readline(),
        tell() and seek(). The function reads from and may reposition the
        file pointer.

    Returns
    -------
    tuple
        A tuple of four integers: (n1, n2, m, bin_size) parsed from the
        header line.

    Behavior
    --------
    - Expects the first non-empty line to be a header beginning with '#'
      and containing tokens of the form: n1=INT, n2=INT, m=INT, bin=INT.
      Tokens may be separated by spaces or tabs and may appear in any order.
      Example: "# n1=5 n2=4 m=3 bin=200"
    - After parsing that line the function reads the next line and, if it
      starts (case-insensitive) with "chrom", treats it as an optional
      column header and leaves the file pointer after it. If the second
      line does not start with "chrom", the function seeks the file
      pointer back to the start of that line so that iteration over the
      file will begin at that line (this is safe for gzip files because
      only a single readline() was advanced).
    - The file pointer is left positioned at the first data row (i.e., just
      after the header lines), so callers can continue iterating over fp.

    Errors
    ------
    Raises SystemExit with a descriptive message in these cases:
    - The first line is missing or does not start with '#'.
    - One or more of n1, n2, m, or bin tokens are missing from the header.
    - The file ends unexpectedly after the header line (no second line).
    """
    # 1) header line with n1/n2/m/bin
    line = fp.readline()
    if not line or not line.startswith("#"):
        raise SystemExit("Expected first line to start with '# n1=... n2=... m=... bin=...'")

    n1 = n2 = m = bin_size = None
    # tokens may be space or tab separated
    for tok in line.strip("# \n").replace("\t", " ").split():
        if tok.startswith("n1="):
            n1 = int(tok.split("=", 1)[1])
        elif tok.startswith("n2="):
            n2 = int(tok.split("=", 1)[1])
        elif tok.startswith("m="):
            m = int(tok.split("=", 1)[1])
        elif tok.startswith("bin="):
            bin_size = int(tok.split("=", 1)[1])

    if any(v is None for v in (n1, n2, m, bin_size)):
        raise SystemExit("Missing one of n1/n2/m/bin in header.")

    # 2) optional column header line
    pos = fp.tell()
    line2 = fp.readline()
    if not line2:
        raise SystemExit("File ended after header.")
    if not line2.lower().startswith("chrom"):
        # not the column header; put file pointer back at start of this line
        # (safe on gzip because we only moved forward one readline)
        fp.seek(pos)

    return n1, n2, m, bin_size

def find_putative_bins(bin_counts_path, eta, out_bins_path=None):
    """
    Find putative genomic bins from per-bin count file using a simple threshold on normalized counts.

    Parameters
    ----------
    bin_counts_path : str or os.PathLike
        Path to a tab-delimited input file containing a header (parsed by parse_header_and_prepare)
        followed by lines with columns: chrom, start, end, x1, x2. The file is opened with smart_open
        so it may be a local path or a remote/streaming path supported by smart_open.
    eta : float
        Tuning parameter used to compute the detection threshold. The function computes
        thresh = 2.0 / (m * eta), where m is obtained from the input file header by
        parse_header_and_prepare.
    out_bins_path : str or os.PathLike, optional
        If provided, the function writes the putative bins (chrom, start, end) to this path
        as a tab-delimited file. Any missing parent directories will be created. If None, no file
        is written.

    Returns
    -------
    tuple
        A pair (putative, thresh) where:
          - putative is a list of tuples (chrom: str, start: int, end: int) for which the
            statistic Fi = (x1 / n1) + (x2 / n2) exceeds the computed threshold.
            n1 and n2 are obtained from the input file header by parse_header_and_prepare.
          - thresh is the float threshold value used for selection.

    Behavior and notes
    ------------------
    - The function skips empty lines, comment lines (lines starting with '#'), and header
      lines that start with (case-insensitive) "chrom".
    - Each data line is expected to have exactly five tab-separated fields:
      chrom, start, end, x1, x2. start, end, x1 and x2 are parsed as integers.
    - The file header must be parsable by parse_header_and_prepare, which must return
      (n1, n2, m, bin_size). The bin_size value is not used in this function but is part
      of the expected header format.
    - If out_bins_path is provided, the output file will contain one tab-delimited row per
      putative bin with columns chrom, start, end.
    - Errors due to malformed input lines, missing files, or failures in smart_open or
      parse_header_and_prepare will propagate as exceptions (e.g., ValueError, FileNotFoundError).

    Example
    -------
    >>> putative, thresh = find_putative_bins("bins.tsv.gz", eta=0.5, out_bins_path="putative.bed")
    """
    with smart_open(bin_counts_path, "rt") as f:
        n1, n2, m, bin_size = parse_header_and_prepare(f)
        thresh = 2.0 / (m * eta)

        putative = []
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            if line.lower().startswith("chrom"):
                continue
            chrom, start_s, end_s, x1_s, x2_s = line.rstrip("\n").split("\t")
            start = int(start_s)
            end = int(end_s)
            x1 = int(x1_s)
            x2 = int(x2_s)
            Fi = (x1 / n1) + (x2 / n2)
            if Fi > thresh:
                putative.append((chrom, start, end))

    if out_bins_path:
        os.makedirs(os.path.dirname(out_bins_path) or ".", exist_ok=True)
        with smart_open(out_bins_path, "wt") as outb:
            for chrom, start, end in putative:
                outb.write(f"{chrom}\t{start}\t{end}\n")

    return putative, thresh

def merge_within_1kb(bins, gap_bp=1000):
    """
    Merge consecutive genomic bins that are within a specified gap.

    Parameters
    ----------
    bins : list of tuple
        Iterable of genomic bins, each represented as (chrom, start, end). Chromosome should be
        comparable (e.g. string), and start/end should be integers. The input may be unsorted.
    gap_bp : int, optional
        Maximum allowed gap in base pairs between adjacent bins on the same chromosome for them to be
        merged. Default is 1000.

    Returns
    -------
    list of tuple
        A list of merged regions as (chrom, start, end), sorted by chromosome and start. Consecutive
        bins on the same chromosome whose gap (start of the next bin minus current end) is less than
        or equal to gap_bp will be merged into a single region.

    Behavior and assumptions
    ------------------------
    - The function sorts the input by (chrom, start, end) before merging; the original order is not preserved.
    - Overlapping bins (where start <= current_end) are considered to have zero gap and will be merged.
    - Each returned region has the minimal start and maximal end covering the merged bins.
    - Assumes each bin has end >= start and that coordinates are integers.
    - If `bins` is empty, an empty list is returned.

    Examples
    --------
    >>> bins = [('chr1', 100, 200), ('chr1', 250, 300), ('chr1', 1300, 1400)]
    >>> merge_within_1kb(bins, gap_bp=1000)
    [('chr1', 100, 300), ('chr1', 1300, 1400)]
    """
    if not bins:
        return []

    # sort by chrom, then start
    bins.sort(key=lambda t: (t[0], t[1], t[2]))

    merged = []
    cur_chrom, cur_start, cur_end = bins[0]

    for chrom, start, end in bins[1:]:
        if chrom == cur_chrom and start <= cur_end + gap_bp:
            # extend region
            if end > cur_end:
                cur_end = end
        else:
            merged.append((cur_chrom, cur_start, cur_end))
            cur_chrom, cur_start, cur_end = chrom, start, end

    merged.append((cur_chrom, cur_start, cur_end))
    return merged

def main():
    """
    Main entry point for computing putative histone-modification (HM) bins and
    merged HM regions from binned count data.

    This function parses command-line arguments, identifies putative HM bins using
    a thresholding rule based on the statistic F(i) = x1/n1 + x2/n2, and merges
    adjacent putative bins into contiguous HM regions.

    Behavior:
    - Command-line arguments:
        --in           : Path to input bin counts file (expected format: bin_counts.tsv.gz)
        --eta          : Valid bin fraction (float). Default: 0.7 (recommended for mouse).
        --out-bins     : Output path for BED-like putative bins (chrom, start, end).
                         Default: results/putative_bins.bed.gz.
        --out-regions  : Output path for merged putative HM regions (BED). Default:
                         results/putative_regions.bed.gz.
    - Calls find_putative_bins(inp, eta, out_bins_path=...) which computes putative
      bins based on F(i) and selects bins with F(i) > 2 / (m * eta). That function
      returns the list of putative bins and the computed threshold value.
    - Calls merge_within_1kb(putative_bins, gap_bp=1000) to merge consecutive putative
      bins that are within 1 kb into larger putative HM regions (paper-faithful).
    - Ensures the output directory for regions exists and writes merged regions in
      BED format (chrom, start, end) to --out-regions using smart_open.
    - Emits summary information to stderr including the eta value, the numeric
      threshold, the number of bins passing the threshold, and the number of merged
      regions.

    Return:
    - None. The function performs file I/O and writes results to the specified output
      files. Exceptions from I/O or from the underlying helper functions propagate to
      the caller.

    Example (CLI):
        python 03_find_putative_sites.py --in bin_counts.tsv.gz --eta 0.7
    """
    ap = argparse.ArgumentParser(
        description=("Compute putative histone-modification bins via F(i)=x1/n1+x2/n2 "
                     "with threshold > 2/(m*eta), then merge consecutive putative bins "
                     "within 1 kb into putative HM regions (paper-faithful).")
    )
    ap.add_argument("--in", dest="inp", required=True, help="bin_counts.tsv.gz from Script 02")
    ap.add_argument("--eta", type=float, default=0.7, help="valid bin fraction (default 0.7 for mouse)")
    ap.add_argument("--out-bins", default="results/putative_bins.bed.gz", help="BED-like putative bins (chrom start end)")
    ap.add_argument("--out-regions", default="results/putative_regions.bed.gz", help="merged putative HM regions (BED)")
    args = ap.parse_args()

    putative_bins, thresh = find_putative_bins(args.inp, args.eta, out_bins_path=args.out_bins)
    regions = merge_within_1kb(putative_bins, gap_bp=1000)

    os.makedirs(os.path.dirname(args.out_regions) or ".", exist_ok=True)
    with smart_open(args.out_regions, "wt") as out:
        for chrom, start, end in regions:
            out.write(f"{chrom}\t{start}\t{end}\n")

    sys.stderr.write(f"[putative] eta={args.eta} -> threshold=2/(m*eta)={thresh:.12g}\n")
    sys.stderr.write(f"[putative] bins passing threshold: {len(putative_bins)}\n")
    sys.stderr.write(f"[putative] merged regions: {len(regions)}\n")

if __name__ == "__main__":
    main()
