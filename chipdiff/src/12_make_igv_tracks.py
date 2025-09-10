#!/usr/bin/env python3
import argparse, gzip, os, sys
import re
import numpy as np

def smart_open(p): 
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")

def make_log2_bedgraph(counts_tsv, out_bg, eps=1.0, esc_over_npc=True):
    """
    Input: results/bin_counts.tsv.gz with columns: chrom start end ESC NPC
    Output: bedGraph of log2(ESC/NPC) if esc_over_npc=True, else log2(NPC/ESC)
    """
    n=0
    with smart_open(counts_tsv) as f, open(out_bg, "w") as fo:
        fo.write("track type=bedGraph name=\"log2_{}\" description=\"log2 {} per bin\"\n"
                 .format("ESCoverNPC" if esc_over_npc else "NPCoverESC",
                         "ESC/NPC" if esc_over_npc else "NPC/ESC"))
        for ln in f:
            if not ln.strip() or ln.startswith("#") or ln.lower().startswith("chrom"):
                continue
            chrom, s, e, x1, x2 = ln.strip().split("\t")[:5]
            s = int(s); e = int(e); x1 = float(x1); x2 = float(x2)
            num = (x1 + eps) if esc_over_npc else (x2 + eps)
            den = (x2 + eps) if esc_over_npc else (x1 + eps)
            val = np.log2(num/den)
            fo.write(f"{chrom}\t{s}\t{e}\t{val:.4f}\n")
            n += 1
    print(f"[ok] bedGraph written: {out_bg}  ({n} bins)")


def make_dhms_bed(dhms_bins_bed, out_bed, include_none=False):
    """
    DHMS BED writer: read bins, infer state when missing, skip malformed lines,
    then merge consecutive bins with identical state and write a colored BED.
    """
    def color_of(state):
        if state == "A1":
            return (70, 130, 180)
        if state == "A2":
            return (220, 20, 60)
        return (200, 200, 200)

    by_chr = {}
    with smart_open(dhms_bins_bed) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                s = int(parts[1]); e = int(parts[2])
            except Exception:
                continue
            name = parts[3] if len(parts) > 3 else ""
            state = parts[5].strip() if len(parts) > 5 else ""
            if not state:
                m = re.search(r"\b(A1|A2|NONE)\b", name, flags=re.I)
                state = m.group(1).upper() if m else "NONE"
            if state not in ("A1", "A2") and not include_none:
                continue
            by_chr.setdefault(chrom, []).append((s, e, state))

    merged = []
    for chrom, rows in by_chr.items():
        rows = [r for r in rows if isinstance(r, tuple) and len(r) == 3]
        if not rows:
            continue
        rows.sort()
        cur_s, cur_e, cur_state = rows[0]
        for s, e, state in rows[1:]:
            if state == cur_state and s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                merged.append((chrom, cur_s, cur_e, cur_state))
                cur_s, cur_e, cur_state = s, e, state
        merged.append((chrom, cur_s, cur_e, cur_state))

    with open(out_bed, "w") as fo:
        fo.write('track name=DHMS_states itemRgb=On description="DHMS calls (A1=ESC, A2=NPC)"\n')
        for chrom, s, e, state in sorted(merged):
            r, g, b = color_of(state)
            fo.write(f"{chrom}\t{s}\t{e}\t{state}\t0\t.\t{s}\t{e}\t{r},{g},{b}\n")

    print(f"[ok] DHMS BED written: {out_bed}")


def make_hcne_promoters_bed(xls_path, sheet, out_bed, pad=1000, fold=4.0):
    """
    Minimal reuse of XLS parser to write promoters as BED with colors:
      UP (NPC/ESC >= fold) = red, OTHER = gray
    """
    import pandas as pd
    import re
    CHR_RE   = re.compile(r"^chr[0-9XYM]+$", re.I)
    STRAND_RE= re.compile(r"^[+-]$")
    def coerce_float(s):
        if s is None: return None
        x = str(s).strip().replace(" ", "").replace(",", ".")
        try: return float(x)
        except: return None

    df = pd.read_excel(xls_path, sheet_name=sheet, header=None, dtype=str)
    header_like = df.apply(lambda r: r.astype(str).str.contains(
        r"(RefSeq|Symbol|chrom|strand|tx|expression|K27|DHMS)", case=False, regex=True, na=False
    ).any(), axis=1)
    hdr = header_like[header_like].index[0] if header_like.any() else 0
    data = df.iloc[hdr+1:].fillna("")

    out = []
    for _, r in data.iterrows():
        cells = [str(x).strip() for x in r.tolist() if str(x).strip()!=""]
        if not cells: continue
        chrom = strand = sym = None; ints=[]; floats=[]
        for c in cells:
            lc=c.lower()
            if CHR_RE.match(c): chrom=c; continue
            if STRAND_RE.match(c): strand=c; continue
            if c.isdigit(): ints.append(int(c)); continue
            v = coerce_float(c); 
            if v is not None: floats.append(v); continue
            # treat as symbol-ish fallback
            if sym is None: sym=c
        if not (chrom and strand and ints):
            continue
        txstart, txend = min(ints), max(ints)
        # expressions: take two largest floats (heuristic)
        floats.sort(reverse=True)
        esc=npc=None
        if len(floats)>=2:
            esc, npc = floats[1], floats[0]  # smaller = ESC, larger = NPC (heuristic)
        if esc is None or npc is None: 
            continue
        tss = txstart if strand=="+" else txend
        s = max(0, tss - pad); e = tss + pad
        cls = "NPC_up" if (npc/esc)>=fold else "OTHER"
        out.append((chrom,s,e,sym,cls))

    with open(out_bed,"w") as fo:
        fo.write("track name=HCNE_promoters itemRgb=On description=\"HCNE promoters (NPC_up vs OTHER)\"\n")
        for chrom,s,e,sym,cls in out:
            if cls=="NPC_up":
                rgb="220,20,60"  # red
            else:
                rgb="180,180,180"  # gray
            fo.write(f"{chrom}\t{s}\t{e}\t{sym};{cls}\t0\t.\t{s}\t{e}\t{rgb}\n")
    print(f"[ok] HCNE promoters BED written: {out_bed}")


def main():
    ap = argparse.ArgumentParser(description="Make IGV-ready tracks from ChIPDiff-faithful outputs.")
    ap.add_argument("--counts", default="chipdiff-faithful/results/bin_counts.tsv.gz")
    ap.add_argument("--dhms-bins", default="chipdiff-faithful/results/dhms_bins.bed.gz")
    ap.add_argument("--xls", default="chipdiff-faithful/refs/bioinf-2008-0517-File004.xls")
    ap.add_argument("--sheet", default="HCNE_genes")
    ap.add_argument("--outdir", default="chipdiff-faithful/igv_tracks")
    ap.add_argument("--esc-over-npc", action="store_true", help="Compute log2(ESC/NPC) instead of default log2(NPC/ESC)")
    ap.add_argument("--skip-hcne", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    # A) bedGraph log2 ratio
    out_bg = os.path.join(args.outdir, f"log2_{'ESCoverNPC' if args.esc_over_npc else 'NPCoverESC'}.bedGraph")
    make_log2_bedgraph(args.counts, out_bg, esc_over_npc=args.esc_over_npc)

    # B) DHMS colored BED
    out_dhms = os.path.join(args.outdir, "DHMS_states.bed")
    make_dhms_bed(args.dhms_bins, out_dhms)

    # C) optional HCNE promoters colored by expression call
    if not args.skip_hcne:
        out_hcne = os.path.join(args.outdir, "HCNE_promoters.bed")
        try:
            make_hcne_promoters_bed(args.xls, args.sheet, out_hcne)
        except Exception as e:
            print("[warn] HCNE promoters skipped:", e, file=sys.stderr)

    print("[done] IGV tracks in:", args.outdir)

if __name__ == "__main__":
    main()
