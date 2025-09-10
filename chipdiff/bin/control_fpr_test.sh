#!/usr/bin/env bash
set -euo pipefail

# 1) Split each library into two halves
python3 chipdiff/src/10_split_centers.py \
  --in chipdiff/results/ESC.centers.bed.gz \
  --out1 chipdiff/results/ESC.rep1.bed.gz \
  --out2 chipdiff/results/ESC.rep2.bed.gz \
  --seed 123

python3 chipdiff/src/10_split_centers.py \
  --in chipdiff/results/NPC.centers.bed.gz \
  --out1 chipdiff/results/NPC.rep1.bed.gz \
  --out2 chipdiff/results/NPC.rep2.bed.gz \
  --seed 123

# 2) Build mixed, non–cell-specific libraries
zcat chipdiff/results/ESC.rep1.bed.gz chipdiff/results/NPC.rep1.bed.gz | gzip -c > chipdiff/results/MIXA.centers.bed.gz

zcat chipdiff/results/ESC.rep2.bed.gz chipdiff/results/NPC.rep2.bed.gz | gzip -c > chipdiff/results/MIXB.centers.bed.gz


# 3) Bin counts
python3 chipdiff/src/02_bin_counts.py \
  --esc chipdiff/results/MIXA.centers.bed.gz \
  --npc chipdiff/results/MIXB.centers.bed.gz \
  --chrom-sizes chipdiff/refs/mm9.chrom.sizes \
  --bin 1000 \
  --out chipdiff/results/bin_counts.control.tsv.gz

# 4) Rebuild LUT + train + decode
python3 chipdiff/src/04_build_emission_lut.py \
  --counts chipdiff/results/bin_counts.control.tsv.gz \
  --tau 3.0 --grid 48 \
  --out chipdiff/results/emission_lut.control.npz

python3 chipdiff/src/05_train_hmm.py \
  --obs chipdiff/results/bin_counts.control.tsv.gz \
  --putative chipdiff/results/putative_regions.bed.gz \
  --lut chipdiff/results/emission_lut.control.npz \
  --out chipdiff/results/transitions.control.npz \
  --max-iters 20 --train-regions 10000 --seed 123

python3 chipdiff/src/06_decode_posterior.py \
  --obs chipdiff/results/bin_counts.control.tsv.gz \
  --putative chipdiff/results/putative_regions.bed.gz \
  --lut chipdiff/results/emission_lut.control.npz \
  --trans chipdiff/results/transitions.control.npz \
  --rho 0.95 \
  --out-bins chipdiff/results/dhms_bins.control.bed.gz \
  --out-regions chipdiff/results/dhms_regions.control.bed.gz

# 5) Summarize
n_ctrl=$(zcat chipdiff/results/dhms_regions.control.bed.gz | wc -l)
n_real=$(zcat chipdiff/results/dhms_regions.bed.gz        | wc -l)

echo "Control DHMS (non–cell-specific): $n_ctrl"
echo "Real DHMS (ESC vs NPC):          $n_real"

if [ "$n_real" -gt 0 ]; then
  python3 - <<PY
n_ctrl = int("$n_ctrl"); n_real = int("$n_real")
fpr = 100.0 * n_ctrl / n_real
print(f"Estimated FPR = {n_ctrl}/{n_real} = {fpr:.3f}%")
PY
else
  echo "Estimated FPR = NA (no real DHMS)"
fi
