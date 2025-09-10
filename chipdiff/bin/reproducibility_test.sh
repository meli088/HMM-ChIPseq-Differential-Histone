#!/usr/bin/env bash
set -euo pipefail

# Config (match your main run)
BIN=1000
TAU=3.0
GRID=48
RHO=0.95
SEED=123

W=chipdiff/results
S=chipdiff/src
R=chipdiff/refs

mkdir -p "$W"

# 0) Make sure split reps exist (create if missing)
if [[ ! -f "$W/ESC.rep1.bed.gz" || ! -f "$W/ESC.rep2.bed.gz" ]]; then
  python3 "$S/10_split_centers.py" --in "$W/ESC.centers.bed.gz" \
    --out1 "$W/ESC.rep1.bed.gz" --out2 "$W/ESC.rep2.bed.gz" --seed $SEED
fi

if [[ ! -f "$W/NPC.rep1.bed.gz" || ! -f "$W/NPC.rep2.bed.gz" ]]; then
  python3 "$S/10_split_centers.py" --in "$W/NPC.centers.bed.gz" \
    --out1 "$W/NPC.rep1.bed.gz" --out2 "$W/NPC.rep2.bed.gz" --seed $SEED
fi

# 1) PASS A: ESC.rep1 vs NPC.rep1
echo "[passA] bin counts"
python3 "$S/02_bin_counts.py" \
  --esc "$W/ESC.rep1.bed.gz" \
  --npc "$W/NPC.rep1.bed.gz" \
  --chrom-sizes "$R/mm9.chrom.sizes" \
  --bin $BIN \
  --out "$W/bin_counts.passA.tsv.gz"

echo "[passA] lut"
python3 "$S/04_build_emission_lut.py" \
  --counts "$W/bin_counts.passA.tsv.gz" \
  --tau $TAU --grid $GRID \
  --out "$W/emission_lut.passA.npz"

echo "[passA] train"
python3 "$S/05_train_hmm.py" \
  --obs "$W/bin_counts.passA.tsv.gz" \
  --putative "$W/putative_regions.bed.gz" \
  --lut "$W/emission_lut.passA.npz" \
  --out "$W/transitions.passA.npz" \
  --max-iters 20 --train-regions 10000 --seed $SEED

echo "[passA] decode"
python3 "$S/06_decode_posterior.py" \
  --obs "$W/bin_counts.passA.tsv.gz" \
  --putative "$W/putative_regions.bed.gz" \
  --lut "$W/emission_lut.passA.npz" \
  --trans "$W/transitions.passA.npz" \
  --rho $RHO \
  --out-bins "$W/dhms_bins.passA.bed.gz" \
  --out-regions "$W/dhms_regions.passA.bed.gz"

# 2) PASS B: ESC.rep2 vs NPC.rep2
echo "[passB] bin counts"
python3 "$S/02_bin_counts.py" \
  --esc "$W/ESC.rep2.bed.gz" \
  --npc "$W/NPC.rep2.bed.gz" \
  --chrom-sizes "$R/mm9.chrom.sizes" \
  --bin $BIN \
  --out "$W/bin_counts.passB.tsv.gz"

echo "[passB] lut"
python3 "$S/04_build_emission_lut.py" \
  --counts "$W/bin_counts.passB.tsv.gz" \
  --tau $TAU --grid $GRID \
  --out "$W/emission_lut.passB.npz"

echo "[passB] train"
python3 "$S/05_train_hmm.py" \
  --obs "$W/bin_counts.passB.tsv.gz" \
  --putative "$W/putative_regions.bed.gz" \
  --lut "$W/emission_lut.passB.npz" \
  --out "$W/transitions.passB.npz" \
  --max-iters 20 --train-regions 10000 --seed $SEED

echo "[passB] decode"
python3 "$S/06_decode_posterior.py" \
  --obs "$W/bin_counts.passB.tsv.gz" \
  --putative "$W/putative_regions.bed.gz" \
  --lut "$W/emission_lut.passB.npz" \
  --trans "$W/transitions.passB.npz" \
  --rho $RHO \
  --out-bins "$W/dhms_bins.passB.bed.gz" \
  --out-regions "$W/dhms_regions.passB.bed.gz"

# 3) Reproducibility score = overlap / average(nA, nB)
nA=$(zcat "$W/dhms_regions.passA.bed.gz" | wc -l)
nB=$(zcat "$W/dhms_regions.passB.bed.gz" | wc -l)

overlap=$(
python3 - "$W/dhms_regions.passA.bed.gz" "$W/dhms_regions.passB.bed.gz" <<'PY'
import sys, gzip
def load(p):
    op = gzip.open if p.endswith(".gz") else open
    D={}
    with op(p,'rt') as f:
        for ln in f:
            if not ln.strip() or ln[0]=='#': continue
            parts=ln.split('\t')
            if len(parts)<3: continue
            c=parts[0]; s=int(parts[1]); e=int(parts[2])
            D.setdefault(c,[]).append((s,e))
    for c in D: D[c].sort()
    return D

def count_overlaps(A,B):
    # count regions in A that overlap ANY in B (1bp overlap)
    cnt=0
    for c in set(A)|set(B):
        a=A.get(c,[]); b=B.get(c,[])
        i=j=0
        hit=0
        while i<len(a) and j<len(b):
            s1,e1=a[i]; s2,e2=b[j]
            if e1<=s2: i+=1
            elif e2<=s1: j+=1
            else:
                # overlap — count this A[i], move to next A interval
                cnt+=1
                i+=1
        # next chrom
    return cnt

A=load(sys.argv[1]); B=load(sys.argv[2])
print(count_overlaps(A,B))
PY
)

avg=$(python3 - <<PY
nA=$nA; nB=$nB
print((nA+nB)/2.0)
PY
)

score=$(python3 - <<PY
ov=$overlap
avg=$avg
print(0.0 if avg==0 else 100.0*ov/avg)
PY
)

echo "Reproducibility inputs:"
echo "  passA DHMS: $nA"
echo "  passB DHMS: $nB"
echo "  overlap(A,B): $overlap"
printf "Reproducibility score = overlap / average = %d / %.1f = %.1f%%\n" "$overlap" "$avg" "$score"

# Also write a small summary file
cat > "$W/reproducibility_summary.txt" <<EOF
Reproducibility (two independent passes):
  passA DHMS: $nA
  passB DHMS: $nB
  overlap(A,B): $overlap
  score: ${score}%
EOF

echo "[done] wrote $W/reproducibility_summary.txt"
