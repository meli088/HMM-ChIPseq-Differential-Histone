#!/usr/bin/env python3
import argparse
import gzip
import re
import os
from collections import defaultdict, namedtuple

def smart_open(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")

# --- helpers to classify cell content -------------------------------------------------
CHR_RE   = re.compile(r"^chr[0-9XYM]+$", re.I)
STRAND_RE= re.compile(r"^[+-]$")
INT_RE   = re.compile(r"^\d+$")
NUM_RE   = re.compile(r"^[0-9\s.,Ee+-]+$")

def coerce_float(s):
    if s is None:
        return None
    x = str(s).strip()
    if x == "" or x.lower() in {"na","nan"}:
        return None
    x = x.replace(" ", "").replace(",", ".")  # European decimal
    try:
        return float(x)
    except Exception:
        return None

def norm_gene(s):
    return (s or "").strip().replace(" ", "").upper()

# --- load DHMS regions ---------------------------------------------------------------
Interval = namedtuple("Interval", "start end")
def load_regions(bed_path):
    """
    Load genomic intervals from a BED-like file and return them grouped and sorted by chromosome.

    Parameters
    ----------
    bed_path : str or os.PathLike
        Path to a BED-like file. The file may be plain text or compressed (smart_open is used).
        Lines that are empty or that start with "#", "track" or "browser" are ignored.

    Returns
    -------
    defaultdict(list)
        A dictionary mapping chromosome name (str) to a list of Interval objects sorted
        by their start coordinate. Each valid input line must contain at least three
        tab-separated fields: chrom, start, end. The start and end fields are parsed as int.

    Notes
    -----
    - The function does not merge or deduplicate overlapping intervals; all parsed Interval
      objects are preserved and then sorted by their start attribute.
    - Coordinates are used as provided in the file (no conversion between 0-based and
      1-based conventions is performed).
    - If a line has fewer than three columns it is skipped. If start/end cannot be parsed
      as integers a ValueError will propagate.
    - Requires an Interval class (with a .start attribute) and smart_open to be available
      in the calling module.
    """
    by_chrom = defaultdict(list)
    with smart_open(bed_path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith(("#","track","browser")):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, s, e = parts[0], int(parts[1]), int(parts[2])
            by_chrom[chrom].append(Interval(s, e))
    for c in by_chrom:
        by_chrom[c].sort(key=lambda iv: iv.start)
    return by_chrom

def overlap_any(ps, pe, regs):
    """
    Determine whether a query interval overlaps any region in a sorted list of regions.

    The function checks whether the half-open interval [ps, pe) overlaps at least one region
    in the provided sequence `regs`. Each element of `regs` is expected to have numeric
    `.start` and `.end` attributes (or properties). `regs` must be sorted in ascending order
    by `.start`.

    Parameters
    ----------
    ps : int | float
        Start coordinate of the query interval (inclusive).
    pe : int | float
        End coordinate of the query interval (exclusive).
    regs : Sequence
        Sequence (e.g. list) of region-like objects sorted by `.start`. Each region must
        expose `.start` and `.end` numeric attributes representing a half-open interval
        [start, end).

    Returns
    -------
    bool
        True if any region in `regs` overlaps the interval [ps, pe), False otherwise.

    Behavior / overlap definition
    -----------------------------
    Two intervals are considered to overlap when they share any interior point. Using
    half-open semantics, [a, b) and [c, d) overlap iff not (b <= c or d <= a).
    This function returns True as soon as such an overlap is found.

    Assumptions and complexity
    --------------------------
    - `regs` must be sorted by `.start`. If not sorted, results may be incorrect.
    - The implementation advances to the first region with `.end > ps` and then scans
      forward until `.start >= pe`, so it inspects only regions that could possibly
      overlap the query interval.
    - Time complexity is O(k) where k is the number of scanned regions (worst-case O(n)).
    - If `regs` is empty, the function returns False.
    """
    # regs sorted; two-pointer scan
    if not regs:
        return False
    _lo, hi = 0, len(regs)
    # advance while end <= ps
    i = 0
    while i < hi and regs[i].end <= ps:
        i += 1
    while i < hi and regs[i].start < pe:
        if not (regs[i].end <= ps or regs[i].start >= pe):
            return True
        i += 1
    return False

# --- XLS parser (robust per-row field detection) -------------------------------------
def parse_hcne_xls(xls_path, sheet):
    """
    Parse a single sheet from an Excel workbook that contains HCNE / gene-table data
    and return a list of normalized row dictionaries.

    This function is heuristic: it reads the sheet as raw strings, detects the header
    row by looking for common column labels, and then interprets every subsequent
    non-empty row by scanning tokens and classifying them as chromosome, strand,
    integers, numeric values, boolean-like flags, or free text. The goal is to
    extract a compact representation of each gene/region row with these fields:
        - symbol: gene symbol (string)
        - chrom: chromosome string (e.g. "chr1")
        - strand: strand, typically "+" or "-" (string)
        - txstart: transcript start coordinate (int)
        - txend: transcript end coordinate (int)
        - esc: ESC expression/score (float or None)
        - npc: NPC expression/score (float or None)
        - with_k27: True/False if a K27-like flag is present, otherwise None

    Parameters
    ----------
    xls_path : str or os.PathLike
            Path to the Excel (.xls/.xlsx) file to read.
    sheet : str or int
            Sheet name or sheet index to parse.

    Behavior / heuristics
    ---------------------
    - Uses pandas.read_excel(..., header=None, dtype=str) to preserve raw cell text.
    - Detects the header row by scanning each row for tokens like "RefSeq",
        "Symbol", "chrom", "strand", "tx", "expression", "K27", "DHMS" (case-insensitive).
        The first matching row is assumed to be the header; data rows start after it.
    - Empty cells are normalized to empty strings and skipped when building token lists.
    - Token classification:
            - CHR_RE is used to detect chromosome tokens (assign to `chrom`).
            - STRAND_RE is used to detect strand tokens (assign to `strand`).
            - "yes","y","true" (case-insensitive) set with_k27 = True;
                "no","n","false" set with_k27 = False.
            - INT_RE matches integer-looking tokens -> collected as ints.
            - NUM_RE matches numeric-looking tokens -> coerced with coerce_float and
                collected as floats.
            - Anything else is treated as free text and collected for symbol selection.
    - Gene symbol selection:
            - From the free-text tokens, selects an "alphabetic-ish" candidate: the token
                with the largest count of alphabetic characters (ties broken by length).
    - txstart/txend:
            - If integer tokens exist, txstart = min(ints), txend = max(ints).
            - Otherwise, uses floats > 1e3 converted to ints for coordinates.
    - ESC / NPC signals:
            - From the remaining floats (excluding ones that equal txstart/txend),
                picks the two largest values; assigns the smaller of those to `esc` and
                the larger to `npc` (heuristic).
    - Rows missing any of chrom, strand, symbol, txstart, or txend are skipped.
    - Returns a list of dicts for the successfully parsed rows.

    Return
    ------
    List[dict]
            Each dict contains the keys: "symbol" (str), "chrom" (str), "strand" (str),
            "txstart" (int), "txend" (int), "esc" (float or None), "npc" (float or None),
            "with_k27" (bool or None).

    Dependencies / notes
    --------------------
    - Expects the module-level names CHR_RE, STRAND_RE, INT_RE, NUM_RE and the helper
        function coerce_float to be defined and appropriate for the input data.
    - Assumes numeric tokens may contain commas or other locale formatting; coerce_float
        should handle those cases if present.
    - Will raise pandas-related errors if the file/sheet cannot be read.
    - Because of the many heuristics, inspect results on a representative sheet and
        adjust regexes or coerce_float as needed for your data.
    """
    import pandas as pd
    # read raw, keep everything as string
    df = pd.read_excel(xls_path, sheet_name=sheet, header=None, dtype=str)
    # find header row by looking for typical labels
    header_like = df.apply(
        lambda r: r.astype(str).str.contains(
            r"(RefSeq|Symbol|chrom|strand|tx|expression|K27|DHMS)",
            case=False, regex=True, na=False
        ).any(), axis=1)
    hdr = header_like[header_like].index[0] if header_like.any() else 0
    data = df.iloc[hdr+1:].fillna("")

    rows = []
    for _, r in data.iterrows():
        cells = [str(x).strip() for x in r.tolist() if str(x).strip() != ""]
        if not cells:
            continue

        chrom = strand = symbol = None
        txstart = txend = None
        esc = npc = None
        with_k27 = None

        ints = []
        floats = []
        texts = []

        # first pass: identify obvious tokens
        for c in cells:
            lc = c.lower()
            if CHR_RE.match(c):
                chrom = c
                continue
            if STRAND_RE.match(c):
                strand = c
                continue
            if lc in {"yes","y","true"}:
                with_k27 = True
                continue
            if lc in {"no","n","false"}:
                with_k27 = False
                continue
            if INT_RE.match(c):
                try:
                    ints.append(int(c))
                except Exception:
                    pass
                continue
            if NUM_RE.match(c):
                v = coerce_float(c)
                if v is not None:
                    floats.append(v)
                    continue
            texts.append(c)

        # pick gene symbol: prefer a short alpha token from texts
        # the table often has Symbol + RefSeq in same row; Symbol is nice short text
        if texts:
            # take the alphabetic-ish token with max letters
            cand = sorted(texts, key=lambda s: (-sum(ch.isalpha() for ch in s), len(s)))[0]
            symbol = cand

        # txstart/txend: take min/max of large ints; if none, try from floats (as ints)
        if ints:
            txstart, txend = min(ints), max(ints)
        elif floats:
            ints2 = [int(v) for v in floats if v and v > 1e3]
            if ints2:
                txstart, txend = min(ints2), max(ints2)

        # expressions ESC/NPC: take two remaining largest floats (not equal to tx coords)
        esc_npc_pool = [v for v in floats if not (txstart and abs(v-txstart)<1) and not (txend and abs(v-txend)<1)]
        esc_npc_pool.sort(reverse=True)
        if len(esc_npc_pool) >= 2:
            # WARNING: some sheets log-transform; your screenshot looks like rawish numbers-with-commas
            # We'll assume raw signals; order on sheet is ESC then NPC, but we can’t rely on it.
            # Take smaller as ESC (often lower), larger as NPC if fold>1 points to upreg in NPC.
            a, b = esc_npc_pool[1], esc_npc_pool[0]
            esc, npc = a, b

        # sanity checks
        if not(chrom and strand and symbol and txstart is not None and txend is not None):
            continue

        rows.append({
            "symbol": symbol,
            "chrom": chrom,
            "strand": strand,
            "txstart": int(txstart),
            "txend": int(txend),
            "esc": coerce_float(esc),
            "npc": coerce_float(npc),
            "with_k27": bool(with_k27) if with_k27 is not None else None
        })
    return rows

# --- main -----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Parse HCNE XLS (messy columns), build ±pad promoters, check NPC/ESC fold change and overlap with DHMS regions."
    )
    ap.add_argument("--xls", required=True, help="refs/bioinf-2008-0517-File004.xls")
    ap.add_argument("--sheet", default="HCNE_genes")
    ap.add_argument("--regions", required=True, help="results/dhms_regions.bed.gz")
    ap.add_argument("--pad", type=int, default=1000, help="promoter half-window around TSS")
    ap.add_argument("--fold", type=float, default=4.0, help="NPC/ESC ≥ fold counts as up-regulated")
    ap.add_argument("--out", default="results/hcne_overlap_summary.txt")
    ap.add_argument("--tsv", default="results/hcne_per_gene_overlap.tsv")
    args = ap.parse_args()

    rows = parse_hcne_xls(args.xls, args.sheet)
    if not rows:
        print("[error] Could not parse any rows from XLS. Check sheet name/contents.")
        return

    # Build promoters and call up-regulated
    promoters = defaultdict(list)  # chrom -> list of (start,end,gene)
    up_genes = set()
    all_genes = set()
    for r in rows:
        tss = r["txstart"] if r["strand"] == "+" else r["txend"]
        ps = max(0, tss - args.pad)
        pe = tss + args.pad
        g = norm_gene(r["symbol"])
        promoters[r["chrom"]].append((ps, pe, g))
        all_genes.add(g)
        if r["esc"] is not None and r["esc"] > 0 and r["npc"] is not None:
            if (r["npc"]/r["esc"]) >= args.fold:
                up_genes.add(g)

    for c in promoters:
        promoters[c].sort(key=lambda x: x[0])

    # Load DHMS
    regs = load_regions(args.regions)

    # Overlap
    overlap_genes = set()
    for chrom in set(list(promoters.keys())+list(regs.keys())):
        P = promoters.get(chrom, [])
        R = regs.get(chrom, [])
        if not P or not R:
            continue
        i = j = 0
        while i < len(P) and j < len(R):
            ps, pe, g = P[i]
            rs, re = R[j].start, R[j].end
            if re <= ps:
                j += 1
            elif pe <= rs:
                i += 1
            else:
                overlap_genes.add(g)
                if re < pe:
                    j += 1
                else:
                    i += 1

    others = all_genes - up_genes
    n_up   = len(up_genes)
    n_oth  = len(others)
    up_hits  = sum((g in overlap_genes) for g in up_genes)
    oth_hits = sum((g in overlap_genes) for g in others)

    up_pct  = (100.0*up_hits/n_up) if n_up else 0.0
    oth_pct = (100.0*oth_hits/n_oth) if n_oth else 0.0

    # Print summary
    print("HCNE promoter overlap with DHMSs (±1 kb around TSS):")
    print(f"  Up-regulated (n={n_up}): {up_hits} ({up_pct:.1f}%)")
    print(f"  Others      (n={n_oth}): {oth_hits} ({oth_pct:.1f}%)")
    with open(args.out, "w") as fo:
        fo.write("HCNE promoter overlap with DHMSs (±1 kb around TSS)\n")
        fo.write(f"Up-regulated (n={n_up}): {up_hits} ({up_pct:.1f}%)\n")
        fo.write(f"Others      (n={n_oth}): {oth_hits} ({oth_pct:.1f}%)\n")

    # Per-gene TSV (for inspection / debugging)
    os.makedirs(os.path.dirname(args.tsv) or ".", exist_ok=True)
    with open(args.tsv, "w") as ft:
        ft.write("group\tgene_symbol\toverlaps_DHMS\n")
        for g in sorted(up_genes):
            ft.write(f"UP30\t{g}\t{int(g in overlap_genes)}\n")
        for g in sorted(others):
            ft.write(f"OTHER\t{g}\t{int(g in overlap_genes)}\n")

if __name__ == "__main__":
    main()
