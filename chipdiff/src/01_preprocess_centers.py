#!/usr/bin/env python3
import sys
import os
import gzip
import argparse

def smart_open(path, mode="rt"):
    """
    Open a file normally or with gzip depending on its extension.

    Parameters
    ----------
    path : str
        Path to the input file. If it ends with ".gz", it is opened with gzip.
    mode : str, optional
        File opening mode (default: "rt" = read text).

    Returns
    -------
    file object
        Opened file handle (either gzip or normal).
    """
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def parse_args():
    """Parse command-line arguments for preprocessing fragment centers.

    This function builds and returns an argparse.Namespace describing the
    command-line options used to preprocess aligned reads into approximate
    fragment centers, following the paper's procedure. The preprocessing
    expected workflow is: collapse PCR duplicates by (chrom, tag_site, strand),
    then shift reads by a fixed distance toward the strand to approximate
    fragment centers.

    Command-line options:
    - --in (dest="inp", required): Input TSV (can be gzipped). Expected to have
        at least four whitespace- or tab-separated columns: chrom, start, end, strand.
        Additional columns (e.g., read_id, copies, sequence) are ignored.
    - --out (dest="out", required): Output BED-like file for centers (can be gzipped).
    - --shift (int, default=100): Number of base pairs to shift each read toward
        its strand to approximate the fragment center (positive integer; default 100).
    - --min-pos (int, default=0): Minimum genomic coordinate to keep; centers with
        positions less than this value are dropped (default 0).

    Returns:
            argparse.Namespace containing the parsed arguments with attributes:
            - inp: input path (str)
            - out: output path (str)
            - shift: shift distance (int)
            - min_pos: minimum allowed position (int)

    Example:
            args = parse_args()
            # args.inp, args.out, args.shift, args.min_pos can now be used by the script.
    """
    ap = argparse.ArgumentParser(
        description=("Paper-faithful preprocessing: collapse PCR duplicates "
                     "by (chrom, tag_site, strand), then shift ±100 bp toward "
                     "the strand to approximate fragment centers. "
                     "Input must have at least 4 columns: chrom,start,end,strand; "
                     "extra columns (read_id,copies,sequence) are ignored.")
    )
    ap.add_argument("--in", dest="inp", required=True, help="Input TSV (.gz ok): chrom start end strand [extra cols]")
    ap.add_argument("--out", dest="out", required=True, help="Output BED-like centers (.gz ok)")
    ap.add_argument("--shift", type=int, default=100, help="Shift toward strand (default 100, per paper)")
    ap.add_argument("--min-pos", type=int, default=0, help="Drop centers < this genomic position (default 0)")
    return ap.parse_args()

def main():
    """
    Preprocess tag centers from a single-end BED-like TSV and write shifted 1-bp centers.

    Reads an input TSV file (a.inp) with at least four tab-separated columns per non-empty, non-comment line:
        chrom, start, end, strand
    - Lines starting with '#' or empty lines are ignored.
    - start and end must parse as integers.
    - strand must be '+' or '-'; otherwise a ValueError is raised.
    - Malformed lines raise ValueError with the offending line number.

    Behavior:
    - For each valid record the "tag site" is computed using the single-end convention:
        site = start if strand == '+' else end - 1
    - Duplicates are removed based on the tuple (chrom, site, strand). Deduplication uses a set and therefore
      does not preserve input order; the output order is arbitrary.
    - Each unique tag is shifted to produce a 1-bp "center":
        center = site + a.shift  (for '+')
        center = site - a.shift  (for '-')
    - Centers strictly less than a.min_pos are dropped.
    - The remaining centers are written to the output file (a.out) as 4-column TSV lines:
        chrom<TAB>center<TAB>center+1<TAB>strand

    Side effects and logging:
    - Ensures the output directory exists (os.makedirs(..., exist_ok=True)).
    - Writes per-run summary counts to stderr:
        - number of unique tags before shifting
        - number of centers kept and number dropped (center < a.min_pos)

    Expected attributes on the parsed arguments object `a`:
    - a.inp: path to input file (string)
    - a.out: path to output file (string)
    - a.shift: integer shift amount
    - a.min_pos: integer minimum allowed center position

    Returns:
    - None

    Exceptions:
    - Raises ValueError for malformed input lines (wrong column count, non-integer start/end, invalid strand).
    """
    a = parse_args()
    dedup = set()

    with smart_open(a.inp, "rt") as fin:
        for lineno, line in enumerate(fin, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                raise ValueError(f"Line {lineno}: need ≥4 columns (chrom,start,end,strand). Got: {line}")
            chrom, s_s, e_s, strand = parts[0], parts[1], parts[2], parts[3]
            if strand not in {"+", "-"}:
                raise ValueError(f"Line {lineno}: strand must be '+' or '-', got '{strand}'")
            try:
                start = int(s_s)
                end = int(e_s)
            except ValueError:
                raise ValueError(f"Line {lineno}: start/end not integers: {s_s}/{e_s}")

            # tag site BEFORE shifting (single-end convention)
            site = start if strand == "+" else end - 1
            dedup.add((chrom, site, strand))

    kept = 0
    dropped = 0  # noqa: E702
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with smart_open(a.out, "wt") as fout:
        for chrom, site, strand in dedup:
            center = site + a.shift if strand == "+" else site - a.shift
            if center < a.min_pos:
                dropped += 1
                continue
            fout.write(f"{chrom}\t{center}\t{center+1}\t{strand}\n")
            kept += 1

    sys.stderr.write(f"[preprocess] unique tags (pre-shift): {len(dedup)}\n")
    sys.stderr.write(f"[preprocess] centers kept: {kept}, dropped(<{a.min_pos}): {dropped}\n")

if __name__ == "__main__":
    main()
