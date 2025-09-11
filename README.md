# HMM-ChIPseq-Differential-Histone

Tools and scripts for an HMM-based analysis pipeline (chipdiff-faithful). Includes demo data and scripts to preprocess, bin counts, find regions, train an HMM, decode posteriors, and summarize results. Based on Xu et al. (2008)

## Setup for the demo

- Install Pixi: https://pixi.sh/latest/
- From the repository root, install the environment : `pixi install`
- To run the demo: `pixi run demo`

# What this repo does

- Preprocess aligned ChIP-seq reads into 1-bp centers  
- Bin the genome into fixed windows and count fragments per bin  
- Identify putative DHMS regions  
- Build emission lookup tables  
- Train a two-state HMM (Baum–Welch / Forward–Backward)  
- Decode posterior probabilities  
- Summarize results into tables & plots  
- Annotate promoters and conserved noncoding elements  
- Evaluate reproducibility and performance


## Contents

- chipdiff/src/ — numbered pipeline steps (01–12)  
- chipdiff/demo/ — tiny demo dataset and expected shapes  
- chipdiff/results/ — generated outputs and intermediates (ignored by Git)  
- pixi.toml — environment spec managed by Pixi

## Run the pipeline on real data

The real ChIP-seq datasets used in Xu et al. (2008) can be fetched directly from GEO:
ESC (H3K27me3) → GSM307619
NPC (H3K27me3) → GSM307614
Download the aligned files into chipdiff/data/ with these names:
chipdiff/data/GSM307619_ES.H3K27me3.aligned.txt.gz
chipdiff/data/GSM307614_NP.H3K27me3.aligned.txt.gz
- Run the code step by step by checking the available tasks in pixi.toml.
    e.g : `pixi run preprocess_esc`
- Or execute the full workflow at once with: `pixi run chipdiff_all`

## Outputs for the full run (in chipdiff/results/)

*.centers.bed.gz — preprocessed tag centers
*.bins.bed.gz — binned counts
*.putative.bed.gz — candidate modification sites
lut.pkl — lookup table for emissions
hmm_model.pkl — trained HMM parameters
decoded.bed.gz — posterior-decoded states
Plots: histograms, MA plots, bar/pie charts
Tables: summary of DHMSs, promoter overlap, reproducibility, performance


## License

MIT — see `LICENSE`.










