"""
Plotting module for EDA-RSA analysis.

Usage:
    from plots import plot_overlap_heatmaps, plot_patient_overlap_regions
    plot_overlap_heatmaps(results_by_thresh, save_path="plots/overlap_heatmap.png")
"""

from .overlap import plot_overlap_heatmaps
from .regions import (
    plot_patient_overlap_regions,
    plot_patient_overlap_combined,
    generate_all_threshold_plots,
    pick_best_and_save_figure1a,
)

__all__ = [
    "plot_overlap_heatmaps",
    "plot_patient_overlap_regions",
    "plot_patient_overlap_combined",
    "generate_all_threshold_plots",
    "pick_best_and_save_figure1a",
]
