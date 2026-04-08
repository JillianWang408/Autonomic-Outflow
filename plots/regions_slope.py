"""
Figure 3A: same layout as 1A, but “active” = top X% of |d(signal)/dt| (slope magnitude).
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import rankdata

from eda_rsa_overlap_slope import abs_slope

COL_EDA = np.array([254, 184, 44]) / 255
COL_RSA = np.array([74, 223, 243]) / 255


def _to_binary_active_slope(x: np.ndarray, t: np.ndarray, top_pct: float) -> np.ndarray:
    s = abs_slope(x, t)
    n = len(s)
    ranks = rankdata(s, method="average")
    p = ranks / n
    return (p >= (1 - top_pct / 100)).astype(bool)


def plot_patient_overlap_regions_slope(
    t: np.ndarray,
    eda: np.ndarray,
    rsa: np.ndarray,
    threshold_pct: float,
    max_duration_sec: float | None = None,
    participant_id: str = "",
    figsize: tuple[float, float] = (12, 5),
) -> plt.Figure:
    """Like Figure 1A: z-scored EDA/RSA; colored segments = high-|slope| overlap where both d/dt > 0."""
    fs = len(t) / (t[-1] - t[0]) if len(t) > 1 else 100
    n = len(t)

    if max_duration_sec is not None:
        n_plot = min(n, int(max_duration_sec * fs))
        start = (n - n_plot) // 2
        t_plot = t[start : start + n_plot]
        eda_plot = eda[start : start + n_plot]
        rsa_plot = rsa[start : start + n_plot]
    else:
        t_plot, eda_plot, rsa_plot = t, eda, rsa

    eda_active = _to_binary_active_slope(eda_plot, t_plot, threshold_pct)
    rsa_active = _to_binary_active_slope(rsa_plot, t_plot, threshold_pct)
    d_eda = np.gradient(eda_plot, t_plot)
    d_rsa = np.gradient(rsa_plot, t_plot)
    overlap = eda_active & rsa_active & (d_eda > 0) & (d_rsa > 0)

    def _zscore(x):
        x = np.asarray(x, dtype=float)
        m, s = x.mean(), x.std()
        return (x - m) / s if s > 0 else np.zeros_like(x)

    eda_n = _zscore(eda_plot)
    rsa_n = _zscore(rsa_plot)

    fig, ax = plt.subplots(figsize=figsize)
    ylo = min(eda_n.min(), rsa_n.min()) - 0.1
    yhi = max(eda_n.max(), rsa_n.max()) + 0.1
    ax.plot(t_plot, eda_n, color="gray", linewidth=1.2, alpha=0.9, label="EDA")
    ax.plot(t_plot, rsa_n, color="black", linewidth=1.0, alpha=0.9, label="RSA")
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
    ax.set_title(f"{participant_id} — Top {threshold_pct}% |slope| overlap, positive slope only (3A)")
    ax.legend(loc="upper right", ncol=2)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


def generate_all_threshold_plots_slope(
    data_dir: str = "data",
    thresholds: list[int] = (5, 10, 15, 20, 25, 30, 35, 40),
    max_duration_sec: float | None = None,
    output_dir: str = "plots/3A",
) -> list[dict]:
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
        pid = Path(f).stem.split(" full")[0].replace("_preprocessed", "").strip()
        for thresh in thresholds:
            eda_active = _to_binary_active_slope(eda, t, thresh)
            rsa_active = _to_binary_active_slope(rsa, t, thresh)
            n_int = np.sum(eda_active & rsa_active)
            n_union = np.sum(eda_active | rsa_active)
            jaccard = n_int / n_union if n_union > 0 else 0
            fig = plot_patient_overlap_regions_slope(
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


def pick_best_and_save_figure3a(
    data_dir: str = "data",
    thresholds: list[int] = (5, 10, 15, 20, 25, 30, 35, 40),
    max_duration_sec: float | None = None,
    output_dir: str = "plots/3A",
    prefer_threshold: int | None = 20,
) -> tuple[str, int]:
    from scipy.io import loadmat

    data_path = Path(data_dir)
    mat_files = sorted(data_path.glob("*.mat"))
    if not mat_files:
        raise FileNotFoundError(f"No .mat files in {data_dir}")

    candidates = []
    for f in mat_files:
        mat = loadmat(str(f))
        t = np.asarray(mat["t"]).flatten()
        eda = np.asarray(mat["eda_clean"] if "eda_clean" in mat else mat["eda"]).flatten()
        rsa = np.asarray(mat["rsa_clean"] if "rsa_clean" in mat else mat["rsa"]).flatten()
        pid = Path(f).stem.split(" full")[0].replace("_preprocessed", "").strip()
        for thresh in thresholds:
            eda_active = _to_binary_active_slope(eda, t, thresh)
            rsa_active = _to_binary_active_slope(rsa, t, thresh)
            n_int = np.sum(eda_active & rsa_active)
            n_union = np.sum(eda_active | rsa_active)
            jaccard = n_int / n_union if n_union > 0 else 0
            candidates.append((pid, thresh, jaccard))

    candidates.sort(key=lambda x: x[2], reverse=True)
    top_jaccard = candidates[0][2]
    top_candidates = [c for c in candidates if c[2] >= top_jaccard * 0.75]
    if prefer_threshold is not None and top_candidates:
        best = min(top_candidates, key=lambda c: abs(c[1] - prefer_threshold))
    else:
        best = candidates[0]
    best_pid, best_thresh, best_jaccard = best

    best_file = next(f for f in mat_files if best_pid in f.stem)
    mat = loadmat(str(best_file))
    t = np.asarray(mat["t"]).flatten()
    eda = np.asarray(mat["eda_clean"] if "eda_clean" in mat else mat["eda"]).flatten()
    rsa = np.asarray(mat["rsa_clean"] if "rsa_clean" in mat else mat["rsa"]).flatten()
    fig = plot_patient_overlap_regions_slope(
        t, eda, rsa, threshold_pct=best_thresh,
        max_duration_sec=max_duration_sec,
        participant_id=best_pid,
        figsize=(10, 5),
    )
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path / "Figure3A_draft.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure 3A draft: {best_pid}, threshold {best_thresh}% (Jaccard={best_jaccard:.4f})")
    print(f"Saved: {out_path / 'Figure3A_draft.png'}")
    return best_pid, best_thresh
