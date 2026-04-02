# Project Structure & Quick Reference

## Folder Layout

```
Autonomic/
├── data/                    # Raw .mat files (EDA, RSA waveforms)
│   └── *.mat
├── data/preprocessed/       # Preprocessed .mat (created by preprocess.py --save-data)
├── plots/                  # Python package + output directory (mixed)
│   ├── __init__.py         # Package
│   ├── overlap.py          # Heatmap plotting
│   ├── regions.py          # Figure 1A region plots
│   ├── 1A/                 # Figure 1A outputs (organized by threshold)
│   │   ├── Figure1A_draft.png
│   │   ├── thresh5/        # All patients at 5% threshold
│   │   └── ... thresh40/
│   ├── 1B/                 # Figure 1B heatmap
│   │   └── Jaccard_overlap.png
│   ├── 2/                  # Figure 2 lead-lag (amplitude peaks)
│   │   ├── {pid}.png
│   │   └── Figure2_summary.png
│   ├── 3A/                 # Slope overlap (same layout as 1A)
│   ├── 3B/                 # Jaccard_overlap_slope.png
│   ├── 4/                  # Lead–lag from max |slope|
│   │   ├── {pid}.png
│   │   └── Figure4_summary.png
│   ├── overlap_heatmaps.png
│   ├── preprocessed/       # Preprocess comparison plots
│   └── spectrum/           # Spectrum analysis plots
├── docs/
├── eda_rsa_overlap.py      # Main overlap analysis
├── eda_rsa_overlap_slope.py # Overlap on top-% |slope| (3B / shared math)
├── preprocess.py           # Preprocessing (currently passthrough)
├── plot_figure1A.py        # Figure 1A generator
├── plot_figure1B.py        # Figure 1B heatmap
├── plot_figure2.py         # Figure 2 lead-lag
├── plot_figure3A.py        # Figure 3A (slope overlap)
├── plot_figure3B.py        # Figure 3B (slope Jaccard)
├── plot_figure4.py         # Figure 4 (slope lead-lag)
└── analyze_spectrum.py     # Frequency analysis
```

---

## How to Get Figure 1A

### Step 1: Ensure preprocessed data exists
```bash
python preprocess.py --data-dir data --save-data
```
Creates `data/preprocessed/*.mat`

### Step 2: Generate Figure 1A
```bash
python plot_figure1A.py
```

**Output:** `plots/1A/Figure1A_draft.png`

---

## Where Plots Are Saved

| Script | Output Location |
|--------|-----------------|
| `plot_figure1A.py` | `plots/1A/Figure1A_draft.png` |
| `plot_figure1A.py --all` | `plots/1A/thresh{N}/{participant}.png` (by threshold) |
| `plot_figure1A.py --patient X --threshold 20` | `plots/1A/thresh20/{participant}.png` |
| `plot_figure1B.py` | `plots/1B/Jaccard_overlap.png` |
| `plot_figure2.py` | `plots/2/{pid}.png` + `Figure2_summary.png` |
| `plot_figure3A.py` | `plots/3A/Figure3A_draft.png`, `thresh*/` |
| `plot_figure3B.py` | `plots/3B/Jaccard_overlap_slope.png` |
| `plot_figure4.py` | `plots/4/{pid}.png` + `Figure4_summary.png` |
| `eda_rsa_overlap.py --plot` | `plots/overlap_heatmaps.png` |
| `preprocess.py` | `plots/preprocessed/*.png` |
| `analyze_spectrum.py` | `plots/spectrum/*.png` |

---

## Minimal Run for Figure 1A

```bash
cd /Users/wangzihan/Desktop/Projects/Autonomic
python preprocess.py --data-dir data --save-data   # if data/preprocessed/ doesn't exist
python plot_figure1A.py
# → plots/1A/Figure1A_draft.png
```
