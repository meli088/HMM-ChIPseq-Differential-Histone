#!/usr/bin/env bash
set -euo pipefail

# Minimal one-shot demo runner
# Usage: bash chipdiff/bin/run_demo.sh [OUTPUT_DIR]
OUTDIR="${1:-chipdiff/demo_results}"
LOG="$OUTDIR/demo_run.log"

mkdir -p "$OUTDIR" chipdiff/results chipdiff/demo

echo "[demo] Preparing demo inputs…"
# Idempotent: creates (or verifies) mini demo inputs
bash chipdiff/bin/make_mini_dataset.sh >>"$LOG" 2>&1 || true

echo "[demo] Building LUT (this step is running)…"
python chipdiff/src/04_build_emission_lut.py \
  --counts chipdiff/demo/bin_counts.mini.tsv.gz \
  --tau 3.0 --grid 48 \
  --out chipdiff/results/emission_lut.demo.npz >>"$LOG" 2>&1

echo "[demo] Training HMM (this step is running)…"
python chipdiff/src/05_train_hmm.py \
  --obs chipdiff/demo/bin_counts.mini.tsv.gz \
  --putative chipdiff/demo/putative_regions.mini.bed.gz \
  --lut chipdiff/results/emission_lut.demo.npz \
  --out chipdiff/results/transitions.demo.npz \
  --max-iters 20 --train-regions 10000 --seed 123 >>"$LOG" 2>&1

echo "[demo] Decoding DHMS (this step is running)…"
python chipdiff/src/06_decode_posterior.py \
  --obs chipdiff/demo/bin_counts.mini.tsv.gz \
  --putative chipdiff/demo/putative_regions.mini.bed.gz \
  --lut chipdiff/results/emission_lut.demo.npz \
  --trans chipdiff/results/transitions.demo.npz \
  --rho 0.95 \
  --out-bins chipdiff/results/dhms_bins.demo.bed.gz \
  --out-regions chipdiff/results/dhms_regions.demo.bed.gz >>"$LOG" 2>&1

echo "[demo] Summarizing results…"
python chipdiff/src/07_summarize_results.py \
  --bins chipdiff/results/dhms_bins.demo.bed.gz \
  --regions chipdiff/results/dhms_regions.demo.bed.gz \
  --lut chipdiff/results/emission_lut.demo.npz \
  --trans chipdiff/results/transitions.demo.npz \
  --out-txt chipdiff/results/summary.demo.txt \
  --out-json chipdiff/results/summary.demo.json >>"$LOG" 2>&1

echo "[demo] Making a few figures…"
# Quiet plotting; warnings go to log
python chipdiff/src/11_quick_viz.py \
  --counts chipdiff/demo/bin_counts.mini.tsv.gz \
  --dhms-bins chipdiff/results/dhms_bins.demo.bed.gz \
  --ri-chr all \
  --out-dir chipdiff/results/plots >>"$LOG" 2>&1 || true

echo "[demo] Collecting outputs…"
# Bundle only what the instructor needs to see
rsync -a --include="*/" \
  --include="summary.demo.txt" \
  --include="dhms_bins.demo.bed.gz" \
  --include="dhms_regions.demo.bed.gz" \
  --include="plots/dhms_length_hist.png" \
  --include="plots/dhms_per_chrom.png" \
  --include="plots/dhms_state_bar.png" \
  --include="plots/logratio_hist_ESCoverNPC.png" \
  --include="plots/ma_all_ESCoverNPC.png" \
  --exclude="*" chipdiff/results/ "$OUTDIR"/ >>"$LOG" 2>&1

echo "[demo] Done. Results are in: $OUTDIR"
echo "       (Full log: $LOG)"
