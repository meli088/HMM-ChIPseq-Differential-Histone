# HMM-ChIPseq-Differential-Histone

Tools and scripts for an HMM-based analysis pipeline (chipdiff-faithful). Includes demo data and scripts to preprocess, bin counts, find regions, train an HMM, decode posteriors, and summarize results. Based on Xu et al. (2008)

## Contents
- `chipdiff/src/` — numbered pipeline steps (01–11).
- `chipdiff/demo/` — tiny demo dataset and expected shapes.
- `chipdiff/results/` — generated outputs and intermediates (ignored).
- `pixi.toml` — environment spec managed by Pixi.

## Setup for the demo
- Install Pixi: https://pixi.sh/latest/
- From the repository root, install the environment : `pixi install`
- To run the demo: `pixi run demo`

## Run the pipeline on real data
You can either:
- Run the code step by step by checking the available tasks in pixi.toml.
- Or execute the full workflow at once with: pixi run chipdiff_all

## Windows/WSL notes
- Work inside WSL paths (this repo is under `\\wsl.localhost\Ubuntu\...`) or linux.

## License
MIT — see `LICENSE`.








