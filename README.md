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
# Optional: --max-duration 300  (plot first 5 min; use 0 for full signal)
```

- **EDA**: Low-pass 3 Hz (Butterworth)
- **RSA**: Bandpass 0.12–0.4 Hz (respiratory range)
- Saves raw vs preprocessed comparison plots to `plots/preprocessed/`

## Run Analysis

```bash
python eda_rsa_overlap.py --data-dir data
# Optional: --top-pct 0.20  (default: top 20% = active)
```
