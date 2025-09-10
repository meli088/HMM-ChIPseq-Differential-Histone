#!/usr/bin/env python3
import argparse
import gzip
import random

def smart_open(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def main():
    """
    Split a BED file of tag centers into two reproducible halves (replicate 1 and replicate 2).

    This command-line entrypoint reads an input BED file (plain or gzipped) containing tag
    center records, randomly shuffles the non-empty, non-comment lines, and writes the
    first half to one output and the second half to another output. The shuffle is
    deterministic when a fixed random seed is provided.

    Behavior
    - Input lines that are empty or start with '#' are ignored.
    - The remaining lines are shuffled using Python's random module seeded by --seed.
    - If the number of lines is odd, the first output (rep1) will contain floor(n/2)
        lines and the second output (rep2) will contain the remaining ceil(n/2) lines.
    - Lines are written verbatim to the output files; no BED validation or reformatting
        is performed.
    - The program prints a one-line summary indicating the number of input and output
        lines.

    Command-line arguments
    - --in (required): Path to the input BED[.gz] file.
    - --out1 (required): Path to the output file for replicate 1 (BED[.gz] allowed).
    - --out2 (required): Path to the output file for replicate 2 (BED[.gz] allowed).
    - --seed (optional, int): Random seed for reproducible splits. Defaults to 123.
        Using the same seed will produce the same split every time.
    """
    ap = argparse.ArgumentParser(
        description="Split a BED file of tag centers into two reproducible halves (rep1, rep2)."
    )
    ap.add_argument("--in", dest="infile", required=True, help="Input BED[.gz] file")
    ap.add_argument("--out1", required=True, help="Output BED[.gz] for replicate 1")
    ap.add_argument("--out2", required=True, help="Output BED[.gz] for replicate 2")
    ap.add_argument("--seed", type=int, default=123, help="Random seed (default 123)")
    args = ap.parse_args()

    random.seed(args.seed)

    with smart_open(args.infile, "rt") as fi:
        lines = [ln for ln in fi if ln.strip() and not ln.startswith("#")]

    random.shuffle(lines)

    half = len(lines) // 2
    rep1, rep2 = lines[:half], lines[half:]

    for out, subset in [(args.out1, rep1), (args.out2, rep2)]:
        with smart_open(out, "wt") as fo:
            for ln in subset:
                fo.write(ln)

    print(f"[split] input={len(lines)} lines → rep1={len(rep1)}, rep2={len(rep2)}")

if __name__ == "__main__":
    main()
