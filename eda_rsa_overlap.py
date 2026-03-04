"""
EDA-RSA Overlap Analysis

Uses percentile/rank normalization to define "active" (top 20%) and computes:
- P(B active | A active)
- P(A active | B active)
- Jaccard overlap (event similarity)
"""

import numpy as np
from pathlib import Path
from scipy.io import loadmat


def load_mat_file(filepath: str) -> dict:
    """Load .mat file and return dict with t, FS, eda, rsa (flattened)."""
    mat = loadmat(filepath)
    return {
        "t": np.asarray(mat["t"]).flatten(),
        "FS": float(mat["FS"].flat[0]),
        "eda": np.asarray(mat["eda"]).flatten(),
        "rsa": np.asarray(mat["rsa"]).flatten(),
    }


def percentile_rank(x: np.ndarray) -> np.ndarray:
    """
    Empirical percentile: p(t) = rank(x(t)) / N
    Uses average rank for ties (scipy.stats.rankdata default).
    """
    from scipy.stats import rankdata
    n = len(x)
    ranks = rankdata(x, method="average")
    return ranks / n


def to_binary_active(x: np.ndarray, top_pct: float = 0.20) -> np.ndarray:
    """
    Top `top_pct` (default 20%) → 1 (active), rest → 0.
    Based on percentile rank.
    """
    p = percentile_rank(x)
    return (p >= (1 - top_pct)).astype(np.int32)


def covariance_metrics(eda: np.ndarray, rsa: np.ndarray) -> dict:
    """
    Compute covariance and correlation between EDA and RSA (raw continuous signals).
    """
    eda = np.asarray(eda, dtype=float).flatten()
    rsa = np.asarray(rsa, dtype=float).flatten()
    cov = np.cov(eda, rsa)[0, 1]
    pearson_r = np.corrcoef(eda, rsa)[0, 1] if np.std(eda) > 0 and np.std(rsa) > 0 else np.nan
    from scipy.stats import spearmanr
    spearman_r, _ = spearmanr(eda, rsa)
    return {
        "covariance": float(cov),
        "pearson_r": float(pearson_r),
        "spearman_r": float(spearman_r),
    }


def overlap_metrics(a_active: np.ndarray, b_active: np.ndarray) -> dict:
    """
    Compute overlap metrics between two binary activity vectors.
    A = first signal (e.g. EDA), B = second signal (e.g. RSA).

    Returns:
        - P_B_given_A: P(B active | A active)
        - P_A_given_B: P(A active | B active)
        - jaccard: |A ∩ B| / |A ∪ B|
    """
    a_active = np.asarray(a_active, dtype=bool)
    b_active = np.asarray(b_active, dtype=bool)

    n_a = np.sum(a_active)
    n_b = np.sum(b_active)
    n_intersection = np.sum(a_active & b_active)
    n_union = np.sum(a_active | b_active)

    p_b_given_a = n_intersection / n_a if n_a > 0 else np.nan
    p_a_given_b = n_intersection / n_b if n_b > 0 else np.nan
    jaccard = n_intersection / n_union if n_union > 0 else np.nan

    return {
        "P_B_given_A": p_b_given_a,
        "P_A_given_B": p_a_given_b,
        "jaccard": jaccard,
        "n_A_active": int(n_a),
        "n_B_active": int(n_b),
        "n_intersection": int(n_intersection),
        "n_union": int(n_union),
    }


def analyze_file(
    filepath: str,
    top_pct: float = 0.20,
) -> dict:
    """
    Run full overlap analysis on a single .mat file.
    """
    data = load_mat_file(filepath)
    eda_active = to_binary_active(data["eda"], top_pct=top_pct)
    rsa_active = to_binary_active(data["rsa"], top_pct=top_pct)

    metrics = overlap_metrics(eda_active, rsa_active)
    cov_metrics = covariance_metrics(data["eda"], data["rsa"])
    metrics.update(cov_metrics)
    metrics["file"] = Path(filepath).name
    metrics["n_samples"] = len(data["t"])
    metrics["FS"] = data["FS"]
    metrics["duration_sec"] = len(data["t"]) / data["FS"]

    return metrics


def analyze_all(
    data_dir: str = "data",
    pattern: str = "*.mat",
    top_pct: float = 0.20,
) -> list[dict]:
    """Run analysis on all .mat files in data_dir."""
    data_path = Path(data_dir)
    results = []
    for f in sorted(data_path.glob(pattern)):
        try:
            r = analyze_file(str(f), top_pct=top_pct)
            results.append(r)
        except Exception as e:
            print(f"Error processing {f}: {e}")
    return results


def print_results(results: list[dict]) -> None:
    """Pretty-print analysis results."""
    for r in results:
        print(f"\n--- {r['file']} ---")
        print(f"  N={r['n_samples']}, FS={r['FS']} Hz, duration={r['duration_sec']:.1f}s")
        print(f"  P(RSA active | EDA active) = {r['P_B_given_A']:.4f}")
        print(f"  P(EDA active | RSA active) = {r['P_A_given_B']:.4f}")
        print(f"  Jaccard overlap             = {r['jaccard']:.4f}")
        print(f"  Covariance                  = {r['covariance']:.6f}")
        print(f"  Pearson r                   = {r['pearson_r']:.4f}")
        print(f"  Spearman r                  = {r['spearman_r']:.4f}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="EDA-RSA overlap analysis")
    parser.add_argument("--data-dir", default="data", help="Directory with .mat files")
    parser.add_argument("--top-pct", type=float, default=0.20, help="Top percentile for 'active' (default 0.20)")
    args = parser.parse_args()

    results = analyze_all(args.data_dir, top_pct=args.top_pct)
    print_results(results)
