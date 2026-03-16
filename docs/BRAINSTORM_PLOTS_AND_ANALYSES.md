# Brainstorm: Additional Plots & Analyses

Ideas for extending the EDA–RSA overlap analysis beyond Figures 1A and 1B.

---

## A. Summary & Distribution Plots

| Plot | Description |
|------|-------------|
| **Jaccard by participant** | Bar chart or violin: Jaccard overlap per participant (at fixed threshold). |
| **Metric distributions** | Histograms or KDE of P(EDA\|RSA), Jaccard, Pearson r across participants. |
| **Correlation matrix** | Heatmap of correlations between metrics (Jaccard, P(EDA\|RSA), Pearson, Spearman, duration, etc.). |
| **Participant ranking** | Ranked list or bar chart of participants by overlap strength; highlight high vs low overlap. |
| **Group summary table** | Mean ± SD of all metrics across participants; export to CSV. |

---

## B. Threshold-Sensitivity Plots

| Plot | Description |
|------|-------------|
| **Curves per participant** | Line plot: Jaccard (or P(EDA\|RSA)) vs threshold for each participant; one line per participant. |
| **Optimal threshold** | For each participant, threshold that maximizes Jaccard; histogram of optimal thresholds. |
| **Stability** | Variance or range of Jaccard across thresholds per participant; who is robust vs sensitive? |
| **Asymmetric thresholds** | Heatmap: P(EDA\|RSA) with EDA threshold on x-axis, RSA threshold on y-axis (different percentiles). |

---

## C. Time-Course & Temporal Plots

| Plot | Description |
|------|-------------|
| **Sliding-window overlap** | Jaccard or overlap % in sliding windows (e.g. 2 min) over time; line plot per participant. |
| **Time-of-recording effects** | Compare overlap in first vs middle vs last third of recording. |
| **Overlap density** | Raster: time on x-axis, participants on y-axis; color = local overlap. |
| **Event raster** | When do overlap events occur? Raster of (participant, time) for overlap onsets. |

---

## D. Signal-Level Plots

| Plot | Description |
|------|-------------|
| **EDA vs RSA scatter** | Scatter plot (EDA, RSA) with density; color by overlap (active/both/inactive). |
| **Phase portrait** | 2D trajectory: (EDA, RSA) over time; color by time or overlap state. |
| **Raw waveforms** | Side-by-side EDA and RSA for selected participants; annotate overlap regions. |
| **Amplitude distributions** | Histograms of EDA and RSA values; mark percentile cutoffs. |

---

## E. Frequency-Domain Plots

| Plot | Description |
|------|-------------|
| **Coherence** | Magnitude-squared coherence EDA–RSA vs frequency; mean ± SE across participants. |
| **Cross-spectrum** | Real/imaginary or magnitude/phase of cross-spectrum. |
| **Band-specific coupling** | Coherence or correlation in specific bands (e.g. 0.1–0.2 Hz, 0.2–0.4 Hz). |
| **Time–frequency coherence** | Spectrogram of coherence; when and at which frequencies do EDA and RSA couple? |

---

## F. Lag & Causality Plots

| Plot | Description |
|------|-------------|
| **Cross-correlation** | CCF(EDA, RSA) vs lag; one curve per participant, mean ± SE. |
| **Peak lag distribution** | Histogram of lags at which CCF peaks; does RSA lead EDA or vice versa? |
| **Lag vs overlap** | Scatter: peak lag vs Jaccard; do high-overlap participants show different lags? |
| **Granger / TE** | Bar chart: EDA→RSA vs RSA→EDA strength per participant. |

---

## G. Individual-Difference & Clustering

| Plot | Description |
|------|-------------|
| **Participant clustering** | PCA or t-SNE of metrics; color by cluster; interpret clusters. |
| **Metric correlations** | Scatter: Jaccard vs Pearson; Jaccard vs duration; etc. |
| **Subgroup comparison** | If metadata exists: overlap metrics by group (e.g. condition, age, clinical score). |
| **Outlier detection** | Flag participants with unusual metric combinations. |

---

## H. Statistical & Null-Model Plots

| Plot | Description |
|------|-------------|
| **Surrogate comparison** | Observed Jaccard vs null distribution (phase-randomized or time-shifted); p-value per participant. |
| **Bootstrap CIs** | Bar chart with error bars from bootstrap over time windows. |
| **Permutation test** | Distribution of test statistic under null; observed value marked. |

---

## I. Publication-Ready Composites

| Plot | Description |
|------|-------------|
| **Multi-panel figure** | 1A (example trace) + 1B (heatmaps) + Jaccard distribution + CCF summary. |
| **Supplementary** | One page per participant: waveform + overlap regions + metrics. |
| **Summary dashboard** | Single-page overview: key metrics, distributions, and example traces. |

---

## J. Export & Reporting

| Output | Description |
|--------|-------------|
| **CSV export** | All metrics per participant × threshold; ready for stats (R, SPSS, etc.). |
| **HTML report** | Interactive report with plots and tables; filter by participant. |
| **LaTeX table** | Formatted table of group means for manuscript. |

---

## Quick Reference: What Exists vs New

| Category | Current | New ideas |
|----------|---------|-----------|
| Overlap regions | ✓ Figure 1A | Phase portrait, scatter |
| Heatmaps | ✓ Figure 1B | Asymmetric thresholds |
| Metrics | ✓ Jaccard, P(·\|·), r | More metrics, export |
| Spectrum | ✓ Welch PSD | Coherence, time–freq |
| Temporal | — | Sliding window, CCF, lag |
| Individual | — | Clustering, subgroups |
| Null models | — | Surrogates, bootstrap |
