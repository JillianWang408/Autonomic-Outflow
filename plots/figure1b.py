"""
Figure 1B: P(EDA|RSA) and P(RSA|EDA) heatmaps.

For each participant (y-axis) × threshold (x-axis), with group average in bottom row.
Saves two separate images with colorbar.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from eda_rsa_overlap import analyze_all_thresholds


THRESHOLDS = (0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40)  # 5–40%


def _build_matrix(
    results_by_thresh: dict[float, list[dict]],
    metric: str,
) -> tuple[np.ndarray, list[str], list[str]]:
    """Build matrix: rows = participants + group avg, cols = thresholds."""
    thresh_list = sorted(results_by_thresh.keys())
    results = results_by_thresh[thresh_list[0]]
    n_participants = len(results)

    matrix = np.zeros((n_participants + 1, len(thresh_list)))
    for j, thresh in enumerate(thresh_list):
        for i, r in enumerate(results_by_thresh[thresh]):
            val = r[metric]
            matrix[i, j] = val if not np.isnan(val) else 0
        matrix[n_participants, j] = np.nanmean([r[metric] for r in results_by_thresh[thresh]])

    row_labels = [
        Path(r["file"]).stem.replace(" full EDA & RSA waveform_preprocessed", "").replace("_preprocessed", "")
        for r in results
    ]
    row_labels.append("Group avg")
    col_labels = [f"{int(t * 100)}%" for t in thresh_list]

    return matrix, row_labels, col_labels


def plot_heatmap_single(
    matrix: np.ndarray,
    row_labels: list[str],
    col_labels: list[str],
    title: str,
    save_path: Path,
    vmin: float = 0,
    vmax: float | None = 0.4,
    cmap: str = "viridis",
    figsize: tuple[float, float] = (8, 6),
) -> plt.Figure:
    """Plot single heatmap with colorbar, save to path."""
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(matrix, aspect="auto", vmin=vmin, vmax=vmax, cmap=cmap)
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_xlabel("Threshold (top %)")
    ax.set_ylabel("Participant")
    ax.set_title(title)
    ax.axhline(len(row_labels) - 1.5, color="white", linewidth=1.5, linestyle="--")
    plt.colorbar(im, ax=ax, label="Probability", shrink=0.8)
    fig.tight_layout()

    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig


def generate_figure1b(
    data_dir: str = "data/preprocessed",
    output_dir: str = "plots/1B",
    thresholds: tuple[float, ...] = THRESHOLDS,
    use_preprocessed: bool = True,
) -> tuple[Path, Path]:
    """
    Generate Figure 1B: two heatmap images.
    - P(EDA|RSA): participants × thresholds, bottom row = group avg
    - P(RSA|EDA): same structure

    Returns (path_eda_given_rsa, path_rsa_given_eda).
    """
    results_by_thresh = analyze_all_thresholds(
        data_dir=data_dir,
        thresholds=list(thresholds),
        use_preprocessed=use_preprocessed,
    )

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # P(EDA | RSA)
    mat_eda_rsa, row_labels, col_labels = _build_matrix(results_by_thresh, "P_A_given_B")
    path_eda_rsa = out_path / "P_EDA_given_RSA.png"
    plot_heatmap_single(
        mat_eda_rsa, row_labels, col_labels,
        title="P(EDA active | RSA active)",
        save_path=path_eda_rsa,
    )
    print(f"Saved: {path_eda_rsa}")

    # P(RSA | EDA)
    mat_rsa_eda, _, _ = _build_matrix(results_by_thresh, "P_B_given_A")
    path_rsa_eda = out_path / "P_RSA_given_EDA.png"
    plot_heatmap_single(
        mat_rsa_eda, row_labels, col_labels,
        title="P(RSA active | EDA active)",
        save_path=path_rsa_eda,
    )
    print(f"Saved: {path_rsa_eda}")

    return path_eda_rsa, path_rsa_eda
