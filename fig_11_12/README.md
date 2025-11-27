Figures 11 and 12 (April 8, 2024 eclipse)
==========================================

This directory holds the data and code needed to reproduce the eclipse layer-height time series used in Figure 11 (N5DUP) and Figure 12 (AB5YO).

Layout
- `data/recordings/`: raw WAV captures for WA5FRF->N5DUP (7.2 MHz) and WA5FRF->AB5YO (7.2 MHz, 5.3 MHz).
- `data/samples/`: template chirps used for correlation (`GerryChirp.wav`, `CorrelationV7.wav`, `CorrelationV8.wav`, `Correlation20m.wav`).
- `notebooks/April8ECLIPSE_fig11_fig12.ipynb`: cleaned notebook that generates both figures and saves outputs alongside the notebook as `fig_11_N5DUP2024.jpg` and `fig_12_2024Res.jpg`.
- `environment.yml`: minimal conda environment for running the notebook.

Quick start
1. `cd HF-TDOA-Frontiers-2025`
2. `conda env create -f fig_11_12/environment.yml && conda activate hf-tdoa`
3. Launch Jupyter in the repo root and open `fig_11_12/notebooks/April8ECLIPSE_fig11_fig12.ipynb`.
4. Run all cells; regenerated figures will overwrite the JPGs next to the notebook.
