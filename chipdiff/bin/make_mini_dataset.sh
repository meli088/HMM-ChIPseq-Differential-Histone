#!/usr/bin/env bash
set -euo pipefail

# --- inputs you already have from the full run ---
COUNTS="chipdiff/results/bin_counts.tsv.gz"
PUTATIVE="chipdiff/results/putative_regions.bed.gz"
CHROMSIZES="chipdiff/refs/mm9.chrom.sizes"

# --- where to write the mini demo set (these WILL be committed) ---
DEMODIR="chipdiff/demo"
mkdir -p "$DEMODIR"

# Focus the demo on chr19 ~10–20 Mb (small but interesting region)
CHR="chr19"
START=10000000
END=20000000

# 1) mini chrom.sizes (only chr19 so training is instant)
awk -v c="$CHR" 'BEGIN{OFS="\t"} $1==c {print $0}' "$CHROMSIZES" > "$DEMODIR/mini.chrom.sizes"

# 2) mini bin_counts: keep header+meta+chr19 in window
zcat "$COUNTS" | awk -v c="$CHR" -v s="$START" -v e="$END" '
  BEGIN{OFS="\t"}
  /^-/ || /^#/ || NR==1 {print; next}
  $1==c && $2<e && $3>s {print}
' | gzip -c > "$DEMODIR/bin_counts.mini.tsv.gz"

# 3) mini putative regions (restrict to window; if empty, fall back to chr-wide)
zcat "$PUTATIVE" | awk -v c="$CHR" -v s="$START" -v e="$END" '
  BEGIN{OFS="\t"}
  $1==c && $2<e && $3>s {print}
' | gzip -c > "$DEMODIR/putative_regions.mini.bed.gz" || true

# if the window is too sparse, take all chr19 putatives
if [ ! -s "$DEMODIR/putative_regions.mini.bed.gz" ]; then
  echo "[info] window putatives empty; using all $CHR putatives."
  zcat "$PUTATIVE" | awk -v c="$CHR" '$1==c' | gzip -c > "$DEMODIR/putative_regions.mini.bed.gz"
fi

# 4) pick a small plotting region inside the window for the quickstart
echo -e "$CHR\t$START\t$END" > "$DEMODIR/region.tsv"

echo "[done] Wrote demo files:"
echo "  $DEMODIR/mini.chrom.sizes"
echo "  $DEMODIR/bin_counts.mini.tsv.gz"
echo "  $DEMODIR/putative_regions.mini.bed.gz"
echo "  $DEMODIR/region.tsv"
