# EDA-RSA Overlap Analysis

## Data Structure

Each `.mat` file in `data/` contains:

| Variable | Shape | Description |
|----------|-------|-------------|
| **t** | (1, N) | Temporal vector — time points (in seconds or sample indices) |
| **FS** | (1, 1) | Sampling frequency — 100 Hz |
| **eda** | (1, N) | Electrodermal Activity (EDA) waveform |
| **rsa** | (1, N) | Respiratory Sinus Arrhythmia (RSA) waveform |

- **N** = number of samples (varies per file, e.g. ~185,688 ≈ 31 min at 100 Hz)
- **t**, **eda**, **rsa** are time-aligned (same length)
- Duration in seconds = N / FS

## Implementation Steps

1. **Load data** — Read `.mat` files with `scipy.io.loadmat`
2. **Percentile normalization** — For each signal: `p(t) = rank(x(t)) / N`
3. **Binary activity** — Top 20% → 1 (active), rest → 0
4. **Overlap metrics**:
   - P(RSA active | EDA active)
   - P(EDA active | RSA active)
   - Jaccard overlap = |A ∩ B| / |A ∪ B|

## Environment

```bash
# Option 1: Conda (recommended)
conda env create -f environment.yml
conda activate autonomic

# Option 2: pip + venv
python -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Preprocessing

```bash
python preprocess.py --data-dir data --output-dir plots/preprocessed
# Optional: --eda-mode highpass (default) or lowpass
# Optional: --eda-cutoff 0.05  (highpass) or 1-3 (lowpass)
# Optional: --max-duration 300  (plot first 5 min; use 0 for full)
```

- **EDA**: Highpass 0.05 Hz (phasic component, removes drift). Spectrum: 90% power < 0.15 Hz.
- **RSA**: Bandpass 0.12–0.4 Hz (respiratory range)
- Saves comparison plots to `plots/preprocessed/`
- Use `--save-data` to save preprocessed signals (t, eda, eda_clean, rsa, rsa_clean) to `data/preprocessed/` as .mat
- See [docs/PREPROCESSING_SUMMARY.md](docs/PREPROCESSING_SUMMARY.md) for methods and rationale

## Spectrum Analysis

```bash
python analyze_spectrum.py --data-dir data
```

Analyzes frequency composition of EDA/RSA; saves plots to `plots/spectrum/`

## Run Analysis

Uses preprocessed data by default. Run `python preprocess.py --save-data` first.

```bash
python eda_rsa_overlap.py
# Uses data/preprocessed/ by default
# Optional: --data-dir data/preprocessed
# Optional: --raw  (use raw data from data/)
# Optional: --top-pct 0.20  (default: top 20% = active)
# Optional: --thresholds 0.1 0.2 0.3  (compare multiple thresholds)
# Optional: --plot  (save heatmap to plots/overlap_heatmaps.png when using --thresholds)

## Figure 1A (overlap regions)

```bash
python plot_figure1.py                    # Auto-pick best patient+threshold (20%), save draft
python plot_figure1.py --all              # All patients × thresholds [5,10,...,40]% → plots/1A/thresh{N}/
python plot_figure1.py --patient EC288_pN22 --threshold 20  # Specific
```

- EDA active: yellow (sympathetic)
- RSA active: blue (parasympathetic)
- Thresholds: 5, 10, 15, 20, 25, 30, 35, 40%
- Output: `plots/1A/Figure1A_draft.png` — see [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) for full layout

## Figure 1B (heatmaps)

```bash
python plot_figure1B.py
```

- P(EDA|RSA): participants × thresholds, bottom row = group average
- P(RSA|EDA): same structure
- Output: `plots/1B/P_EDA_given_RSA.png`, `plots/1B/P_RSA_given_EDA.png`

## Brainstorm: More Plots & Analyses

- **[Joint temporal dynamics](docs/ANALYSIS_JOINT_TEMPORAL_DYNAMICS.md)** — cross-correlation, coherence, event timing, Granger causality, HMM regimes, surrogates
- **[Plots & analyses brainstorm](docs/BRAINSTORM_PLOTS_AND_ANALYSES.md)** — summary plots, threshold sensitivity, time-course, frequency-domain, clustering, export
