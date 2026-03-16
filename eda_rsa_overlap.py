"""
EDA-RSA Overlap Analysis

Uses preprocessed data by default. Percentile/rank normalization to define
"active" (top 20%) and computes:
- P(B active | A active)
- P(A active | B active)
- Jaccard overlap (event similarity)
"""

import numpy as np
from pathlib import Path
from scipy.io import loadmat


def load_mat_file(filepath: str, use_preprocessed: bool = True) -> dict:
    """
    Load .mat file. When use_preprocessed=True, uses eda_clean, rsa_clean if present;
    otherwise eda, rsa (raw).
    """
    mat = loadmat(filepath)
    t = np.asarray(mat["t"]).flatten()
    fs = float(mat["FS"].flat[0])
    if use_preprocessed and "eda_clean" in mat and "rsa_clean" in mat:
        eda = np.asarray(mat["eda_clean"]).flatten()
        rsa = np.asarray(mat["rsa_clean"]).flatten()
        source = "preprocessed"
    else:
        eda = np.asarray(mat["eda"]).flatten()
        rsa = np.asarray(mat["rsa"]).flatten()
        source = "raw"
    return {"t": t, "FS": fs, "eda": eda, "rsa": rsa, "_source": source}


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
    use_preprocessed: bool = True,
) -> dict:
    """
    Run full overlap analysis on a single .mat file.
    Uses preprocessed (eda_clean, rsa_clean) when available.
    """
    data = load_mat_file(filepath, use_preprocessed=use_preprocessed)
    eda_active = to_binary_active(data["eda"], top_pct=top_pct)
    rsa_active = to_binary_active(data["rsa"], top_pct=top_pct)

    metrics = overlap_metrics(eda_active, rsa_active)
    cov_metrics = covariance_metrics(data["eda"], data["rsa"])
    metrics.update(cov_metrics)
    metrics["file"] = Path(filepath).name
    metrics["data_source"] = data.get("_source", "unknown")
    metrics["n_samples"] = len(data["t"])
    metrics["FS"] = data["FS"]
    metrics["duration_sec"] = len(data["t"]) / data["FS"]

    return metrics


def analyze_all(
    data_dir: str = "data/preprocessed",
    pattern: str = "*.mat",
    top_pct: float = 0.20,
    use_preprocessed: bool = True,
) -> list[dict]:
    """Run analysis on all .mat files in data_dir. Default: preprocessed data."""
    data_path = Path(data_dir)
    if not data_path.exists():
        raise FileNotFoundError(
            f"Data directory '{data_dir}' not found. Run: python preprocess.py --data-dir data --save-data"
        )
    results = []
    for f in sorted(data_path.glob(pattern)):
        try:
            r = analyze_file(str(f), top_pct=top_pct, use_preprocessed=use_preprocessed)
            results.append(r)
        except Exception as e:
            print(f"Error processing {f}: {e}")
    return results


def print_results(results: list[dict]) -> None:
    """Pretty-print analysis results."""
    for r in results:
        src = r.get("data_source", "unknown")
        print(f"\n--- {r['file']} ({src}) ---")
        print(f"  N={r['n_samples']}, FS={r['FS']} Hz, duration={r['duration_sec']:.1f}s")
        print(f"  P(RSA active | EDA active) = {r['P_B_given_A']:.4f}")
        print(f"  P(EDA active | RSA active) = {r['P_A_given_B']:.4f}")
        print(f"  Jaccard overlap             = {r['jaccard']:.4f}")
        print(f"  Covariance                  = {r['covariance']:.6f}")
        print(f"  Pearson r                   = {r['pearson_r']:.4f}")
        print(f"  Spearman r                  = {r['spearman_r']:.4f}")


def analyze_all_thresholds(
    data_dir: str = "data/preprocessed",
    pattern: str = "*.mat",
    thresholds: list[float] = (0.10, 0.15, 0.20, 0.25, 0.30),
    use_preprocessed: bool = True,
) -> dict[float, list[dict]]:
    """Run analysis for multiple thresholds. Returns {threshold: [results per file]}."""
    results_by_thresh = {}
    for thresh in thresholds:
        results_by_thresh[thresh] = analyze_all(
            data_dir, pattern, top_pct=thresh, use_preprocessed=use_preprocessed
        )
    return results_by_thresh


def print_threshold_comparison(results_by_thresh: dict[float, list[dict]]) -> None:
    """Print comparison table across thresholds."""
    thresh_list = sorted(results_by_thresh.keys())
    if not thresh_list:
        return
    results = results_by_thresh[thresh_list[0]]
    if not results:
        return

    print("\n" + "=" * 80)
    print("P(RSA active | EDA active) by threshold (top %)")
    print("=" * 80)
    header = f"{'File':<45} " + " ".join(f"{t*100:>5.0f}%" for t in thresh_list)
    print(header)
    print("-" * len(header))
    for i, r in enumerate(results):
        fname = r["file"][:44] + "…" if len(r["file"]) > 45 else r["file"]
        row = f"{fname:<45}"
        for thresh in thresh_list:
            val = results_by_thresh[thresh][i]["P_B_given_A"]
            row += f" {val:>5.3f}"
        print(row)
    print("-" * len(header))
    mean_row = f"{'Mean':<45}"
    for thresh in thresh_list:
        vals = [r["P_B_given_A"] for r in results_by_thresh[thresh]]
        mean_row += f" {np.mean(vals):>5.3f}"
    print(mean_row)
    print()

    print("=" * 80)
    print("Jaccard overlap by threshold (top %)")
    print("=" * 80)
    print(header)
    print("-" * len(header))
    for i, r in enumerate(results):
        fname = r["file"][:44] + "…" if len(r["file"]) > 45 else r["file"]
        row = f"{fname:<45}"
        for thresh in thresh_list:
            val = results_by_thresh[thresh][i]["jaccard"]
            row += f" {val:>5.3f}"
        print(row)
    print("-" * len(header))
    mean_row = f"{'Mean':<45}"
    for thresh in thresh_list:
        vals = [r["jaccard"] for r in results_by_thresh[thresh]]
        mean_row += f" {np.mean(vals):>5.3f}"
    print(mean_row)
    print()

    # n_intersection (threshold-dependent)
    print("=" * 80)
    print("n_intersection (|EDA active ∩ RSA active|) by threshold")
    print("=" * 80)
    print(header)
    print("-" * len(header))
    for i, r in enumerate(results):
        fname = r["file"][:44] + "…" if len(r["file"]) > 45 else r["file"]
        row = f"{fname:<45}"
        for thresh in thresh_list:
            val = results_by_thresh[thresh][i]["n_intersection"]
            row += f" {val:>6d}"
        print(row)
    print("-" * len(header))
    mean_row = f"{'Mean':<45}"
    for thresh in thresh_list:
        vals = [r["n_intersection"] for r in results_by_thresh[thresh]]
        mean_row += f" {np.mean(vals):>6.0f}"
    print(mean_row)
    print()

    # Continuous metrics (threshold-independent)
    print("=" * 80)
    print("Continuous metrics (threshold-independent)")
    print("=" * 80)
    r0 = results_by_thresh[thresh_list[0]]
    print(f"{'File':<45} {'Covariance':>12} {'Pearson r':>10} {'Spearman r':>10}")
    print("-" * 80)
    for r in r0:
        fname = r["file"][:44] + "…" if len(r["file"]) > 45 else r["file"]
        print(f"{fname:<45} {r['covariance']:>12.6f} {r['pearson_r']:>10.4f} {r['spearman_r']:>10.4f}")
    print("-" * 80)
    print(f"{'Mean':<45} {np.mean([r['covariance'] for r in r0]):>12.6f} "
          f"{np.mean([r['pearson_r'] for r in r0]):>10.4f} {np.mean([r['spearman_r'] for r in r0]):>10.4f}")
    print()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="EDA-RSA overlap analysis (uses preprocessed data by default)")
    parser.add_argument("--data-dir", default="data/preprocessed",
                        help="Directory with .mat files (default: preprocessed)")
    parser.add_argument("--top-pct", type=float, default=0.20, help="Top percentile for 'active' (default 0.20)")
    parser.add_argument("--thresholds", type=float, nargs="+", default=None,
                        help="Test multiple thresholds (e.g. 0.1 0.2 0.3). Overrides --top-pct")
    parser.add_argument("--plot", action="store_true",
                        help="Save heatmap when using --thresholds")
    parser.add_argument("--plot-dir", default="plots", help="Directory for saved plots")
    parser.add_argument("--raw", action="store_true", help="Use raw data instead of preprocessed")
    args = parser.parse_args()

    if args.thresholds:
        results_by_thresh = analyze_all_thresholds(
            args.data_dir,
            thresholds=args.thresholds,
            use_preprocessed=not args.raw,
        )
        print_threshold_comparison(results_by_thresh)
        if args.plot:
            from plots import plot_overlap_heatmaps
            plot_overlap_heatmaps(
                results_by_thresh,
                save_path=Path(args.plot_dir) / "overlap_heatmaps.png",
            )
    else:
        results = analyze_all(
            args.data_dir,
            top_pct=args.top_pct,
            use_preprocessed=not args.raw,
        )
        print_results(results)
