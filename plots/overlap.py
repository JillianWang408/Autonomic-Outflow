"""
Overlap analysis heatmap plots.

P(EDA|RSA) and P(RSA|EDA) by participant × threshold, with group average row.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def _build_matrix(
    results_by_thresh: dict[float, list[dict]],
    metric: str,
) -> tuple[np.ndarray, list[str], list[str]]:
    """
    Build matrix: rows = participants + group avg, cols = thresholds.
    Returns (matrix, row_labels, col_labels).
    """
    thresh_list = sorted(results_by_thresh.keys())
    results = results_by_thresh[thresh_list[0]]
    n_participants = len(results)

    matrix = np.zeros((n_participants + 1, len(thresh_list)))
    for j, thresh in enumerate(thresh_list):
        for i, r in enumerate(results_by_thresh[thresh]):
            matrix[i, j] = r[metric]
        matrix[n_participants, j] = np.mean([r[metric] for r in results_by_thresh[thresh]])

    # Short participant labels (filename stem, e.g. EC238_pN1)
    row_labels = [
        Path(r["file"]).stem.replace(" full EDA & RSA waveform_preprocessed", "").replace("_preprocessed", "")
        for r in results
    ]
    row_labels.append("Group avg")
    col_labels = [f"{int(t*100)}%" for t in thresh_list]

    return matrix, row_labels, col_labels


def plot_overlap_heatmaps(
    results_by_thresh: dict[float, list[dict]],
    save_path: str | Path | None = "plots/overlap_heatmaps.png",
    figsize: tuple[float, float] = (10, 8),
    vmin: float | None = 0,
    vmax: float | None = 0.4,
    cmap: str = "viridis",
) -> plt.Figure:
    """
    Plot P(EDA|RSA) and P(RSA|EDA) heatmaps.

    Args:
        results_by_thresh: {threshold: [results per file]} from analyze_all_thresholds()
        save_path: Where to save figure. None = don't save.
        figsize: Figure size.
        vmin, vmax: Color scale limits.
        cmap: Colormap name.

    Returns:
        matplotlib Figure
    """
    thresh_list = sorted(results_by_thresh.keys())
    if not thresh_list:
        raise ValueError("results_by_thresh is empty")

    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)

    # P(EDA|RSA)
    mat_eda_rsa, row_labels, col_labels = _build_matrix(results_by_thresh, "P_A_given_B")
    im0 = axes[0].imshow(mat_eda_rsa, aspect="auto", vmin=vmin, vmax=vmax, cmap=cmap)
    axes[0].set_xticks(range(len(col_labels)))
    axes[0].set_xticklabels(col_labels)
    axes[0].set_yticks(range(len(row_labels)))
    axes[0].set_yticklabels(row_labels, fontsize=8)
    axes[0].set_xlabel("Threshold (top %)")
    axes[0].set_ylabel("Participant")
    axes[0].set_title("P(EDA active | RSA active)")
    axes[0].axhline(len(row_labels) - 1.5, color="white", linewidth=1, linestyle="--")

    # P(RSA|EDA)
    mat_rsa_eda, _, _ = _build_matrix(results_by_thresh, "P_B_given_A")
    im1 = axes[1].imshow(mat_rsa_eda, aspect="auto", vmin=vmin, vmax=vmax, cmap=cmap)
    axes[1].set_xticks(range(len(col_labels)))
    axes[1].set_xticklabels(col_labels)
    axes[1].set_yticks([])
    axes[1].set_xlabel("Threshold (top %)")
    axes[1].set_title("P(RSA active | EDA active)")
    axes[1].axhline(len(row_labels) - 1.5, color="white", linewidth=1, linestyle="--")

    # Shared colorbar
    cbar = fig.colorbar(im1, ax=axes, shrink=0.6, label="Probability")
    fig.tight_layout()

    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig


def plot_single_heatmap(
    results_by_thresh: dict[float, list[dict]],
    metric: str = "P_B_given_A",
    title: str = "P(RSA active | EDA active)",
    save_path: str | Path | None = None,
    **kwargs,
) -> plt.Figure:
    """Plot a single heatmap for one metric."""
    matrix, row_labels, col_labels = _build_matrix(results_by_thresh, metric)
    vmin = kwargs.pop("vmin", 0)
    vmax = kwargs.pop("vmax", 0.4)
    figsize = kwargs.pop("figsize", (6, 6))

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(matrix, aspect="auto", vmin=vmin, vmax=vmax, cmap=kwargs.get("cmap", "viridis"))
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_xlabel("Threshold (top %)")
    ax.set_ylabel("Participant")
    ax.set_title(title)
    ax.axhline(len(row_labels) - 1.5, color="white", linewidth=1, linestyle="--")
    plt.colorbar(im, ax=ax, label="Probability")
    fig.tight_layout()

    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig
