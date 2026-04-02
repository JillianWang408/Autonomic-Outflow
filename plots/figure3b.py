"""
Figure 3B: Jaccard heatmap when “active” = top X% of |d(signal)/dt| (slope magnitude).
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from eda_rsa_overlap_slope import analyze_all_thresholds_slope


THRESHOLDS = (0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40)


def _build_matrix(
    results_by_thresh: dict[float, list[dict]],
    metric: str,
) -> tuple[np.ndarray, list[str], list[str]]:
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
        Path(r["file"]).stem.split(" full")[0].replace("_preprocessed", "").strip()
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
    cbar_label: str = "Jaccard",
    cmap: str = "plasma",
    figsize: tuple[float, float] = (8, 6),
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(matrix, aspect="auto", vmin=vmin, vmax=vmax, cmap=cmap, interpolation="nearest")
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_xlabel("Threshold (top % of |slope|)")
    ax.set_ylabel("Participant")
    ax.set_title(title)
    ax.axhline(len(row_labels) - 1.5, color="white", linewidth=1.5, linestyle="--")
    plt.colorbar(im, ax=ax, label=cbar_label, shrink=0.8)
    fig.tight_layout()
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig


def generate_figure3b(
    data_dir: str = "data",
    output_dir: str = "plots/3B",
    thresholds: tuple[float, ...] = THRESHOLDS,
    use_preprocessed: bool = False,
) -> Path:
    results_by_thresh = analyze_all_thresholds_slope(
        data_dir=data_dir,
        thresholds=list(thresholds),
        use_preprocessed=use_preprocessed,
    )
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    mat, row_labels, col_labels = _build_matrix(results_by_thresh, "jaccard")
    save_file = out_path / "Jaccard_overlap_slope.png"
    vmin = 0
    vmax = max(0.1, float(np.nanmax(mat)))
    plot_heatmap_single(
        mat, row_labels, col_labels,
        title="Jaccard overlap — top % of |dEDA/dt| & |dRSA/dt|",
        save_path=save_file,
        vmin=vmin,
        vmax=min(1.0, vmax * 1.05),
        cmap="plasma",
    )
    print(f"Saved: {save_file}")
    return save_file
