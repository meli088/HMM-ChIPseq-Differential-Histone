#!/usr/bin/env python3
import gzip
import argparse
import os
import numpy as np

def smart_open(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")

def file_exists(p): 
    return p and os.path.exists(p)

def _np_to_py(x):
    """
    Convert NumPy arrays and NumPy scalar types into native Python types.

    This function normalizes values that come from NumPy so they are more readily
    serializable and easier to work with in pure-Python contexts. Behaviors:

    - If x is a numpy.ndarray:
        - If x.shape == () (a 0-d array), it is converted to the corresponding
            Python scalar using .item().
        - Otherwise, the array is converted to a nested list, and the conversion is
            applied recursively to each element.
    - If x is an instance of a NumPy scalar type (np.generic), it is converted to
        the corresponding Python scalar via .item().
    - For all other types, the value is returned unchanged.

    Parameters
    ----------
    x : object
            The value to convert. Typically a numpy.ndarray, numpy scalar, or any other
            Python object.

    Returns
    -------
    object
            A Python-native representation of the input: lists for arrays, Python
            scalars for NumPy scalar values or 0-d arrays, or the original object for
            non-NumPy inputs.

    Notes
    -----
    - The conversion is recursive for array contents.
    - The function is defensive around isinstance checks of NumPy types to remain
        tolerant of unusual NumPy subclasses or runtime environments.
    """
    if isinstance(x, np.ndarray):
        if x.shape == ():           # 0-d array
            return _np_to_py(x.item())
        return [ _np_to_py(v) for v in x.tolist() ]
    try:
        if isinstance(x, (np.generic,)):
            return x.item()
    except Exception:
        pass
    return x

def main():
    ap = argparse.ArgumentParser(description="Summarize DHMS results (short terminal + full files).")
    ap.add_argument("--bins", required=True, help="results/dhms_bins.bed[.gz]")
    ap.add_argument("--regions", required=True, help="results/dhms_regions.bed[.gz]")
    ap.add_argument("--lut", default="results/emission_lut.npz", help="results/emission_lut.npz (metadata)")
    ap.add_argument("--trans", default="results/transitions.npz", help="results/transitions.npz (metadata)")
    ap.add_argument("--out-txt", default="results/summary.txt", help="output text summary")
    ap.add_argument("--out-json", default="results/summary.json", help="output JSON summary")
    args = ap.parse_args()

    # Count bins and states
    n_bins = a1 = a2 = 0
    with smart_open(args.bins) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            state = parts[5]
            if state == "A1":
                n_bins += 1
                a1 += 1
            elif state == "A2":
                n_bins += 1
                a2 += 1

    # Regions
    n_regions = sum(1 for line in smart_open(args.regions) if line.strip())

    esc_pct = (100.0 * a1 / (a1 + a2)) if (a1 + a2) > 0 else None

    # Metadata (save but don’t print)
    meta = {}
    if file_exists(args.lut):
        d = np.load(args.lut, allow_pickle=True)
        meta["lut"] = { k: _np_to_py(d[k]) for k in d.files }
    if file_exists(args.trans):
        t = np.load(args.trans, allow_pickle=True)
        meta["train"] = { k: _np_to_py(t[k]) for k in t.files }

    # Short terminal output
    print("Summary of DHMS results:")
    print(f"  DHMS bins:    {n_bins}")
    print(f"  DHMS regions: {n_regions}")
    if esc_pct is not None:
        print(f"  ESC-enriched: {a1} ({esc_pct:.1f}%)")
        print(f"  NPC-enriched: {a2} ({100.0 - esc_pct:.1f}%)")
    print(f"\n[info] Full metadata written to {args.out_txt} and {args.out_json}")

    # Save full TXT
    with open(args.out_txt, "w") as fo:
        fo.write(f"DHMS bins:    {n_bins}\n")
        fo.write(f"DHMS regions: {n_regions}\n")
        fo.write(f"ESC-enriched: {a1} ({esc_pct:.1f}%)\n")
        fo.write(f"NPC-enriched: {a2} ({100.0 - esc_pct:.1f}%)\n\n")


if __name__ == "__main__":
    main()
