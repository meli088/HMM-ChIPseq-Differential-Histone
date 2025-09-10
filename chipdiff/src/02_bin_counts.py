#!/usr/bin/env python3
import sys
import os
import gzip
import argparse
from collections import defaultdict

def smart_open(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def read_chrom_sizes(path):
    """
    Read chromosome sizes from a text file and return them as a list of tuples.

    Each non-empty, non-comment line in the input file is expected to contain at
    least two whitespace-separated fields: a chromosome name and its size. Lines
    that are empty, start with '#' (treated as comments), contain fewer than two
    fields, or have a non-integer size are silently skipped.

    Parameters
    ----------
    path : str or path-like
        Path or URI to the chromosome sizes file. The file is opened using
        smart_open, so local files and supported remote paths (e.g., S3) are allowed.

    Returns
    -------
    list of tuple
        A list of (chromosome, size) tuples where `chromosome` is a str and
        `size` is an int, in the same order they appear in the file.

    Notes
    -----
    - Only the first two fields on each valid line are used; extra fields are ignored.
    - Leading and trailing whitespace is trimmed from each line before parsing.
    - Invalid lines are skipped without raising exceptions.

    Example
    -------
    Given a file "example.sizes" with lines:
        chr1 248956422
        chr2 242193529

    The function returns:
        [('chr1', 248956422), ('chr2', 242193529)]
    """
    sizes = []
    with smart_open(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            chrom, size_s = parts[0], parts[1]
            try:
                size = int(size_s)
            except ValueError:
                continue
            sizes.append((chrom, size))
    return sizes

def load_centers(path):
    """
    Yield (chrom, pos) tuples parsed from a BED-like centers file.

    Parameters
    ----------
    path : str or path-like
        Path or URI to a tab-delimited BED-like file. Each non-comment line is
        expected to contain at least three columns: chrom, start, end. Additional
        columns (e.g. strand) are permitted but ignored.

    Yields
    ------
    tuple
        A 2-tuple (chrom, pos) where:
        - chrom (str): chromosome/contig name from the first column.
        - pos (int): integer position parsed from the 'start' field (second column).
          The 'start' value is treated as the point coordinate (consistent with
          BED-style 0-based starts).

    Behavior
    --------
    - Lines that are empty or start with '#' are skipped.
    - Lines with fewer than three tab-separated fields are skipped.
    - If the start field cannot be converted to an integer, that line is skipped
      (no exception is propagated).
    - The function does not validate the 'end' value or perform bounds checking.
    - The file is opened with smart_open, so 'path' may refer to local files or
      supported remote URIs.

    Examples
    --------
    >>> # centers.bed (tab-separated)
    >>> #chrom\tstart\tend
    >>> chr1\t1000\t1001
    >>> chr2\t2000\t2001
    >>> list(load_centers('centers.bed'))
    [('chr1', 1000), ('chr2', 2000)]
    """
    """Yield (chrom, pos) for each center (BED-like: chrom start end strand). We use 'start' as the point."""
    with smart_open(path, "rt") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                pos = int(parts[1])  # start == center position
            except ValueError:
                continue
            yield chrom, pos

def main():
    """
    Tile the genome into fixed-size bins, count point-centers from two samples per bin,
    and write a gzipped TSV of per-bin counts plus summary header.

    This function is intended to be used as a command-line entry point. It parses
    arguments (see argparse usage in the function) to obtain:
    - --esc: path to ESC centers (BED-like single-position records; .gz accepted)
    - --npc: path to NPC centers (BED-like single-position records; .gz accepted)
    - --chrom-sizes: chrom sizes file (tab-separated: chrom, length)
    - --bin: integer bin size in base pairs (default 1000)
    - --out: output path for a gzipped TSV

    Behavior:
    - Reads chromosome sizes via read_chrom_sizes(args.chrom_sizes). If no chromosomes
        are parsed, the function exits with SystemExit.
    - Iterates over centers from the two input files via load_centers(), counting
        how many centers fall into each (chromosome, bin_index) where
        bin_index = position // bin_size. Positions are treated as 0-based integers.
    - Centers with chromosomes not present in the chrom sizes file or with positions
        < 0 or >= chrom_length are ignored.
    - Computes the total number of bins m as the sum over chromosomes of
        ceil(chrom_length / bin_size). Bins cover intervals [start, end) where
        start = b * bin_size and end = min(start + bin_size, chrom_length). The
        final bin on a chromosome may be smaller than bin_size.
    - Ensures the output directory exists, opens the output via smart_open(..., "wt")
        so that .gz output is supported, and writes a header line of the form
        "# n1={n1}\tn2={n2}\tm={m}\tbin={bin_size}" followed by a column header
        "chrom\tstart\tend\tx1\tx2".
    - Writes one TSV row per bin, in the same order as chroms in the chrom sizes
        input, giving per-bin counts x1 (ESC) and x2 (NPC).
    - After writing the file, prints a short summary to stderr:
        "[bin] wrote {total_bins} bins; n1={n1}, n2={n2}"

    Notes and assumptions:
    - The implementation relies on helper functions read_chrom_sizes, load_centers,
        and smart_open to exist in the module or be imported; their expected behavior
        is as used above (read_chrom_sizes returns an iterable/list of (chrom, length),
        load_centers yields (chrom, pos) pairs, and smart_open supports writing to .gz).
    - The function returns None and performs its effects via file I/O and stderr.
    """
    ap = argparse.ArgumentParser(
        description=("Tile genome into fixed-size bins (default 1kb) using chrom.sizes; "
                     "count ESC/NPC centers per bin to produce x1,x2. "
                     "Outputs gzipped TSV with header carrying n1,n2,m,bin.")
    )
    ap.add_argument("--esc", required=True, help="ESC centers (BED-like points; .gz ok).")
    ap.add_argument("--npc", required=True, help="NPC centers (BED-like points; .gz ok).")
    ap.add_argument("--chrom-sizes", required=True, help="chrom sizes file (e.g., refs/mm9.chrom.sizes).")
    ap.add_argument("--bin", type=int, default=1000, help="bin size (default 1000).")
    ap.add_argument("--out", required=True, help="output .tsv.gz")
    args = ap.parse_args()

    bin_size = args.bin
    chrom_sizes = read_chrom_sizes(args.chrom_sizes)
    if not chrom_sizes:
        raise SystemExit(f"No chrom sizes parsed from {args.chrom_sizes}")

    # count centers per (chrom, bin_idx)
    c1 = defaultdict(int)
    c2 = defaultdict(int)

    # helper: increment counts only if position falls within known chrom size
    valid_chrom = {c for c, _ in chrom_sizes}
    chrom_len = {c: L for c, L in chrom_sizes}

    n1 = 0
    for chrom, pos in load_centers(args.esc):
        if chrom not in valid_chrom: 
            continue
        if pos < 0 or pos >= chrom_len[chrom]:
            continue
        bin_idx = pos // bin_size
        c1[(chrom, bin_idx)] += 1
        n1 += 1

    n2 = 0
    for chrom, pos in load_centers(args.npc):
        if chrom not in valid_chrom:
            continue
        if pos < 0 or pos >= chrom_len[chrom]:
            continue
        bin_idx = pos // bin_size
        c2[(chrom, bin_idx)] += 1
        n2 += 1

    # Write bins in genomic order
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with smart_open(args.out, "wt") as out:
        # number of bins m (sum over chrom of ceil(size/bin))
        m = 0
        for chrom, size in chrom_sizes:
            m += (size + bin_size - 1) // bin_size
        out.write(f"# n1={n1}\tn2={n2}\tm={m}\tbin={bin_size}\n")
        out.write("chrom\tstart\tend\tx1\tx2\n")

        total_bins = 0
        for chrom, size in chrom_sizes:
            nb = (size + bin_size - 1) // bin_size
            for b in range(nb):
                start = b * bin_size
                end = min(start + bin_size, size)
                x1 = c1.get((chrom, b), 0)
                x2 = c2.get((chrom, b), 0)
                out.write(f"{chrom}\t{start}\t{end}\t{x1}\t{x2}\n")
                total_bins += 1

    sys.stderr.write(f"[bin] wrote {total_bins} bins; n1={n1}, n2={n2}\n")

if __name__ == "__main__":
    main()
