#!/usr/bin/env bash
set -euo pipefail
W=chipdiff/results

real=$(zcat "$W/dhms_regions.bed.gz" | wc -l)
ctrl=0
if [[ -f "$W/dhms_regions.control.bed.gz" ]]; then
  ctrl=$(zcat "$W/dhms_regions.control.bed.gz" | wc -l)
fi
fpr="NA"
if [[ "${real}" -gt 0 ]]; then
  fpr=$(python3 - <<PY
real=$real; ctrl=$ctrl
print(f"{100.0*ctrl/real:.3f}%")
PY
)
fi

# HCNE detection: try to parse the % from the Up-regulated line
hcne="NA"
if [[ -f "$W/hcne_overlap_summary.txt" ]]; then
  # expects a line like: Up-regulated (n=42): 27 (64.3%)
  hcne=$(grep -Eo 'Up-regulated \(n=[0-9]+\): [0-9]+ \([0-9.]+%\)' "$W/hcne_overlap_summary.txt" | grep -Eo '\([0-9.]+%\)' | tr -d '()%')
  [[ -z "$hcne" ]] && hcne="NA"
fi

repro="NA"
if [[ -f "$W/reproducibility_summary.txt" ]]; then
  repro=$(grep -Eo 'score: [0-9.]+%' "$W/reproducibility_summary.txt" | awk '{print $2}')
fi

out="$W/performance_summary.tsv"
{
  echo -e "Method\tNo. DHMS (cell-specific)\tFPR from control\tDetection on HCNE (up-reg)\tReproducibility score"
  echo -e "ChIPDiff\t${real}\t${fpr}\t${hcne}\t${repro}"
  # If you later compute fold-change/qualitative baselines, add rows here.
} > "$out"

echo "[table] wrote $out"
cat "$out"
