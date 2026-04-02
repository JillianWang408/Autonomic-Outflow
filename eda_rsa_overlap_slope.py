"""
EDA–RSA overlap using |d(signal)/dt|: rank top X% of slope magnitude as “active”.

Mirrors eda_rsa_overlap metrics for Figure 3B (Jaccard heatmap).
"""

import numpy as np
from pathlib import Path

from eda_rsa_overlap import covariance_metrics, load_mat_file, overlap_metrics, to_binary_active


def abs_slope(x: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Pointwise |dx/dt| using numpy.gradient."""
    x = np.asarray(x, dtype=float).flatten()
    t = np.asarray(t, dtype=float).flatten()
    if len(x) != len(t):
        raise ValueError("x and t must match length")
    return np.abs(np.gradient(x, t))


def analyze_file_slope(
    filepath: str,
    top_pct: float = 0.20,
    use_preprocessed: bool = True,
) -> dict:
    """Overlap metrics when “active” = top top_pct of |dEDA/dt| and |dRSA/dt| (ranked separately)."""
    data = load_mat_file(filepath, use_preprocessed=use_preprocessed)
    t = data["t"]
    eda_s = abs_slope(data["eda"], t)
    rsa_s = abs_slope(data["rsa"], t)
    eda_active = to_binary_active(eda_s, top_pct=top_pct)
    rsa_active = to_binary_active(rsa_s, top_pct=top_pct)

    metrics = overlap_metrics(eda_active, rsa_active)
    cov_metrics = covariance_metrics(data["eda"], data["rsa"])
    metrics.update(cov_metrics)
    metrics["file"] = Path(filepath).name
    metrics["data_source"] = data.get("_source", "unknown")
    metrics["n_samples"] = len(t)
    metrics["FS"] = data["FS"]
    metrics["duration_sec"] = len(t) / data["FS"]
    return metrics


def analyze_all_slope(
    data_dir: str = "data",
    pattern: str = "*.mat",
    top_pct: float = 0.20,
    use_preprocessed: bool = False,
) -> list[dict]:
    data_path = Path(data_dir)
    if not data_path.exists():
        raise FileNotFoundError(f"Data directory '{data_dir}' not found.")
    results = []
    for f in sorted(data_path.glob(pattern)):
        try:
            r = analyze_file_slope(str(f), top_pct=top_pct, use_preprocessed=use_preprocessed)
            results.append(r)
        except Exception as e:
            print(f"Error processing {f}: {e}")
    return results


def analyze_all_thresholds_slope(
    data_dir: str = "data",
    pattern: str = "*.mat",
    thresholds: list[float] = (0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40),
    use_preprocessed: bool = False,
) -> dict[float, list[dict]]:
    results_by_thresh = {}
    for thresh in thresholds:
        results_by_thresh[thresh] = analyze_all_slope(
            data_dir, pattern, top_pct=thresh, use_preprocessed=use_preprocessed
        )
    return results_by_thresh


# Re-export for scripts that expect same names
__all__ = [
    "abs_slope",
    "analyze_file_slope",
    "analyze_all_slope",
    "analyze_all_thresholds_slope",
]
