"""
Overlap region plots: EDA (yellow) and RSA (blue) above-threshold regions.

Figure 1A: Color-coded time course showing when each signal is active.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import rankdata


# Sympathetic (EDA) and parasympathetic (RSA) colors
COL_EDA = np.array([254, 184, 44]) / 255   # yellow
COL_RSA = np.array([74, 223, 243]) / 255   # blue
COL_OVERLAP = np.array([0.4, 0.7, 0.4])   # green for overlap


def _to_binary_active(x: np.ndarray, top_pct: float) -> np.ndarray:
    """Top top_pct → 1 (active), rest → 0."""
    x = np.asarray(x, dtype=float).flatten()
    n = len(x)
    ranks = rankdata(x, method="average")
    p = ranks / n
    return (p >= (1 - top_pct / 100)).astype(bool)  # top_pct in percent


def plot_patient_overlap_regions(
    t: np.ndarray,
    eda: np.ndarray,
    rsa: np.ndarray,
    threshold_pct: float,
    max_duration_sec: float | None = None,
    participant_id: str = "",
    figsize: tuple[float, float] = (12, 5),
) -> plt.Figure:
    """
    Single-panel plot: EDA and RSA in one graph with color-coded regions.
    - EDA only active: yellow
    - RSA only active: blue
    - Overlap (both active): green (color blend)

    Args:
        t, eda, rsa: Time and signal vectors
        threshold_pct: Top X% = active (e.g. 20 for 20%)
        max_duration_sec: Plot middle N seconds. None = full.
        participant_id: Label for title.
    """
    fs = len(t) / (t[-1] - t[0]) if len(t) > 1 else 100
    n = len(t)

    if max_duration_sec is not None:
        n_plot = min(n, int(max_duration_sec * fs))
        start = (n - n_plot) // 2  # middle window
        t_plot = t[start : start + n_plot]
        eda_plot = eda[start : start + n_plot]
        rsa_plot = rsa[start : start + n_plot]
    else:
        t_plot = t
        eda_plot = eda
        rsa_plot = rsa

    eda_active = _to_binary_active(eda_plot, threshold_pct)
    rsa_active = _to_binary_active(rsa_plot, threshold_pct)
    eda_only = eda_active & ~rsa_active
    rsa_only = rsa_active & ~eda_active
    overlap = eda_active & rsa_active

    # Z-score both for display on same axes (mean=0, std=1)
    def _zscore(x):
        x = np.asarray(x, dtype=float)
        m, s = x.mean(), x.std()
        return (x - m) / s if s > 0 else np.zeros_like(x)

    eda_n = _zscore(eda_plot)
    rsa_n = _zscore(rsa_plot)

    fig, ax = plt.subplots(figsize=figsize)

    ylo = min(eda_n.min(), rsa_n.min()) - 0.1
    yhi = max(eda_n.max(), rsa_n.max()) + 0.1

    # Full traces in gray/black (base layer)
    ax.plot(t_plot, eda_n, color="gray", linewidth=1.2, alpha=0.9, label="EDA")
    ax.plot(t_plot, rsa_n, color="black", linewidth=1.0, alpha=0.9, label="RSA")

    # Overlay colored line segments where overlap
    pad = np.concatenate(([False], overlap, [False]))
    d = np.diff(pad.astype(int))
    starts = np.where(d == 1)[0]
    ends = np.where(d == -1)[0]
    for s, e in zip(starts, ends):
        ax.plot(t_plot[s:e], eda_n[s:e], color=COL_EDA, linewidth=1.5, alpha=1, zorder=5)
        ax.plot(t_plot[s:e], rsa_n[s:e], color=COL_RSA, linewidth=1.3, alpha=1, zorder=5)

    ax.set_ylim(ylo, yhi)
    ax.set_ylabel("Z-score")
    ax.set_xlabel("Time (s)")
    ax.set_title(f"{participant_id} — Top {threshold_pct}% active (colored = overlap)")
    ax.legend(loc="upper right", ncol=2)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


def plot_patient_overlap_combined(
    t: np.ndarray,
    eda: np.ndarray,
    rsa: np.ndarray,
    threshold_pct: float,
    max_duration_sec: float | None = None,
    participant_id: str = "",
    figsize: tuple[float, float] = (12, 4),
) -> plt.Figure:
    """
    Single-panel plot: EDA (yellow) and RSA (blue) as horizontal bands,
    overlap (green) where both active. Raster-style visualization.
    """
    fs = len(t) / (t[-1] - t[0]) if len(t) > 1 else 100
    n = len(t)

    if max_duration_sec is not None:
        n_plot = min(n, int(max_duration_sec * fs))
        start = (n - n_plot) // 2
        t_plot = t[start : start + n_plot]
        eda_plot = eda[start : start + n_plot]
        rsa_plot = rsa[start : start + n_plot]
    else:
        t_plot = t
        eda_plot = eda
        rsa_plot = rsa

    eda_active = _to_binary_active(eda_plot, threshold_pct)
    rsa_active = _to_binary_active(rsa_plot, threshold_pct)
    overlap = eda_active & rsa_active

    fig, ax = plt.subplots(figsize=figsize)

    # Stacked bands: EDA row (y=2), RSA row (y=1), Overlap row (y=0)
    ax.fill_between(t_plot, 1.5, 2.5, where=eda_active, color=COL_EDA, alpha=0.8, label="EDA active")
    ax.fill_between(t_plot, 0.5, 1.5, where=rsa_active, color=COL_RSA, alpha=0.8, label="RSA active")
    ax.fill_between(t_plot, -0.5, 0.5, where=overlap, color=COL_OVERLAP, alpha=0.8, label="Overlap")

    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(["Overlap", "RSA", "EDA"])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("")
    ax.set_title(f"{participant_id} — Top {threshold_pct}% active")
    ax.legend(loc="upper right", ncol=3)
    ax.set_ylim(-0.6, 2.6)
    ax.grid(True, alpha=0.3, axis="x")
    fig.tight_layout()
    return fig


def generate_all_threshold_plots(
    data_dir: str = "data",
    thresholds: list[int] = (5, 10, 15, 20, 25, 30, 35, 40),
    max_duration_sec: float | None = None,
    output_dir: str = "plots/1A",
) -> list[dict]:
    """
    Generate overlap region plots for all patients × all thresholds.
    Output: output_dir/thresh{N}/ for each threshold, with {participant}.png per patient.
    Returns list of {participant, threshold, jaccard, path}.
    """
    from scipy.io import loadmat

    data_path = Path(data_dir)
    base_out = Path(output_dir)
    base_out.mkdir(parents=True, exist_ok=True)

    mat_files = sorted(data_path.glob("*.mat"))
    if not mat_files:
        raise FileNotFoundError(f"No .mat files in {data_dir}")

    results = []
    for f in mat_files:
        mat = loadmat(str(f))
        t = np.asarray(mat["t"]).flatten()
        eda = np.asarray(mat["eda_clean"] if "eda_clean" in mat else mat["eda"]).flatten()
        rsa = np.asarray(mat["rsa_clean"] if "rsa_clean" in mat else mat["rsa"]).flatten()

        pid = Path(f).stem.replace(" full EDA & RSA waveform_preprocessed", "").replace("_preprocessed", "")

        for thresh in thresholds:
            eda_active = _to_binary_active(eda, thresh)
            rsa_active = _to_binary_active(rsa, thresh)
            n_int = np.sum(eda_active & rsa_active)
            n_union = np.sum(eda_active | rsa_active)
            jaccard = n_int / n_union if n_union > 0 else 0

            fig = plot_patient_overlap_regions(
                t, eda, rsa, threshold_pct=thresh,
                max_duration_sec=max_duration_sec,
                participant_id=pid,
            )
            thresh_dir = base_out / f"thresh{thresh}"
            thresh_dir.mkdir(parents=True, exist_ok=True)
            save_file = thresh_dir / f"{pid}.png"
            fig.savefig(save_file, dpi=150, bbox_inches="tight")
            plt.close(fig)
            results.append({"participant": pid, "threshold": thresh, "jaccard": jaccard, "path": str(save_file)})

    return results


def pick_best_and_save_figure1a(
    data_dir: str = "data",
    thresholds: list[int] = (5, 10, 15, 20, 25, 30, 35, 40),
    max_duration_sec: float | None = None,
    output_dir: str = "plots/1A",
    prefer_threshold: int | None = 20,
) -> tuple[str, int]:
    """
    Pick best patient+threshold and save Figure 1A draft.

    Selection: Among (patient, threshold) with Jaccard in top 25%, prefer
    prefer_threshold (default 20%) for visually significant peaks.
    If prefer_threshold is None, use highest Jaccard.
    """
    from scipy.io import loadmat

    data_path = Path(data_dir)
    mat_files = sorted(data_path.glob("*.mat"))
    if not mat_files:
        raise FileNotFoundError(f"No .mat files in {data_dir}")

    # Collect all (pid, thresh, jaccard)
    candidates = []
    for f in mat_files:
        mat = loadmat(str(f))
        t = np.asarray(mat["t"]).flatten()
        eda = np.asarray(mat["eda_clean"] if "eda_clean" in mat else mat["eda"]).flatten()
        rsa = np.asarray(mat["rsa_clean"] if "rsa_clean" in mat else mat["rsa"]).flatten()
        pid = Path(f).stem.replace(" full EDA & RSA waveform_preprocessed", "").replace("_preprocessed", "")

        for thresh in thresholds:
            eda_active = _to_binary_active(eda, thresh)
            rsa_active = _to_binary_active(rsa, thresh)
            n_int = np.sum(eda_active & rsa_active)
            n_union = np.sum(eda_active | rsa_active)
            jaccard = n_int / n_union if n_union > 0 else 0
            candidates.append((pid, thresh, jaccard))

    # Select: prefer threshold near prefer_threshold among top Jaccard
    candidates.sort(key=lambda x: x[2], reverse=True)
    top_jaccard = candidates[0][2]
    top_candidates = [c for c in candidates if c[2] >= top_jaccard * 0.75]  # within 25% of best
    if prefer_threshold is not None and top_candidates:
        # Prefer closest to prefer_threshold
        best = min(top_candidates, key=lambda c: abs(c[1] - prefer_threshold))
    else:
        best = candidates[0]
    best_pid, best_thresh, best_jaccard = best

    # Generate Figure 1A
    best_file = next(f for f in mat_files if best_pid in f.stem)
    mat = loadmat(str(best_file))
    t = np.asarray(mat["t"]).flatten()
    eda = np.asarray(mat["eda_clean"] if "eda_clean" in mat else mat["eda"]).flatten()
    rsa = np.asarray(mat["rsa_clean"] if "rsa_clean" in mat else mat["rsa"]).flatten()

    fig = plot_patient_overlap_regions(
        t, eda, rsa, threshold_pct=best_thresh,
        max_duration_sec=max_duration_sec,
        participant_id=best_pid,
        figsize=(10, 5),
    )
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path / "Figure1A_draft.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure 1A draft: {best_pid}, threshold {best_thresh}% (Jaccard={best_jaccard:.4f})")
    print(f"Saved: {out_path / 'Figure1A_draft.png'}")
    return best_pid, best_thresh
