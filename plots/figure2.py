"""
Figure 2: Lead-lag distribution (peak-to-peak).

Two histograms: EDA leads vs RSA leads.
Visual check: EDA/RSA traces with detected peaks marked.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import find_peaks
from scipy.io import loadmat


def _find_peaks(
    x: np.ndarray,
    t: np.ndarray,
    fs: float,
    min_dist_sec: float = 5,
    prominence_pct: float = 5,
    top_pct: float | None = 30,
) -> tuple[np.ndarray, np.ndarray]:
    """Find peaks; returns (peak_indices, peak_times). If top_pct, keep only peaks in top X% by amplitude."""
    x = np.asarray(x, dtype=float).flatten()
    min_dist = int(min_dist_sec * fs)
    prom = np.percentile(np.abs(np.diff(x)), prominence_pct) if prominence_pct > 0 else None
    prom = max(prom, 1e-9) if prom is not None else None
    peaks, _ = find_peaks(x, distance=min_dist, prominence=prom)
    if top_pct is not None and len(peaks) > 0:
        thresh = np.percentile(x[peaks], 100 - top_pct)  # keep peaks >= this (top top_pct)
        peaks = peaks[x[peaks] >= thresh]
    return peaks, t[peaks]


def _compute_peak_to_peak_lags(
    t: np.ndarray,
    eda: np.ndarray,
    rsa: np.ndarray,
    fs: float,
    top_peak_pct: float = 30,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    For each EDA peak, find nearest RSA peak. Lag = t_RSA - t_EDA.
    Only uses top top_peak_pct% of peaks by amplitude (default 30%).
    Returns (lags_eda_leads, lags_rsa_leads, eda_peak_times, rsa_peak_times).
    """
    eda_peaks, t_eda = _find_peaks(eda, t, fs, top_pct=top_peak_pct)
    rsa_peaks, t_rsa = _find_peaks(rsa, t, fs, top_pct=top_peak_pct)

    if len(eda_peaks) == 0 or len(rsa_peaks) == 0:
        return np.array([]), np.array([]), t_eda, t_rsa

    lags_eda_leads = []  # t_RSA - t_EDA > 0
    lags_rsa_leads = []  # t_RSA - t_EDA < 0
    for te in t_eda:
        idx = np.argmin(np.abs(t_rsa - te))
        lag = t_rsa[idx] - te
        if lag > 0:
            lags_eda_leads.append(lag)
        elif lag < 0:
            lags_rsa_leads.append(-lag)  # store as positive "RSA lead time"

    return (
        np.array(lags_eda_leads),
        np.array(lags_rsa_leads),
        t_eda,
        t_rsa,
    )


def _plot_single_patient(
    t_win: np.ndarray,
    eda_win: np.ndarray,
    rsa_win: np.ndarray,
    t_eda: np.ndarray,
    t_rsa: np.ndarray,
    lags_eda: np.ndarray,
    lags_rsa: np.ndarray,
    pid: str,
    top_peak_pct: float,
    bin_sec: float,
) -> plt.Figure:
    """Plot one patient: visual check + 2 histograms."""
    fig = plt.figure(figsize=(12, 10))

    ax0 = fig.add_subplot(3, 1, 1)
    eda_z = (eda_win - eda_win.mean()) / (eda_win.std() or 1)
    rsa_z = (rsa_win - rsa_win.mean()) / (rsa_win.std() or 1)
    ax0.plot(t_win, eda_z, color="gray", linewidth=0.8, alpha=0.8, label="EDA")
    ax0.plot(t_win, rsa_z, color="black", linewidth=0.8, alpha=0.8, label="RSA")
    ax0.scatter(t_eda, np.interp(t_eda, t_win, eda_z), color="orange", s=30, zorder=5, label="EDA peaks")
    ax0.scatter(t_rsa, np.interp(t_rsa, t_win, rsa_z), color="steelblue", s=30, zorder=5, label="RSA peaks")
    ax0.set_ylabel("Z-score")
    ax0.set_xlabel("Time (s)")
    ax0.set_title(f"{pid} — top {top_peak_pct:.0f}% peaks")
    ax0.legend(loc="upper right", ncol=2, fontsize=8)
    ax0.grid(True, alpha=0.3)

    ax1 = fig.add_subplot(3, 1, 2)
    if len(lags_eda) > 0:
        mx = max(lags_eda.max(), bin_sec)
        bins = np.arange(0, mx + bin_sec, bin_sec)
        ax1.hist(lags_eda, bins=bins, color="orange", alpha=0.7, edgecolor="darkorange")
        ax1.axvline(np.median(lags_eda), color="red", linestyle="--", linewidth=1.5, label=f"Median = {np.median(lags_eda):.1f} s")
    ax1.set_xlabel("Lag (s) — RSA peak after EDA peak")
    ax1.set_ylabel("Count")
    ax1.set_title("EDA leads")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2 = fig.add_subplot(3, 1, 3)
    if len(lags_rsa) > 0:
        mx = max(lags_rsa.max(), bin_sec)
        bins = np.arange(0, mx + bin_sec, bin_sec)
        ax2.hist(lags_rsa, bins=bins, color="steelblue", alpha=0.7, edgecolor="navy")
        ax2.axvline(np.median(lags_rsa), color="red", linestyle="--", linewidth=1.5, label=f"Median = {np.median(lags_rsa):.1f} s")
    ax2.set_xlabel("Lag (s) — EDA peak after RSA peak")
    ax2.set_ylabel("Count")
    ax2.set_title("RSA leads")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    return fig


def plot_figure2(
    data_dir: str = "data",
    output_dir: str = "plots/2",
    max_duration_sec: float = 300,
    bin_sec: float = 2,
    top_peak_pct: float = 30,
) -> list[Path]:
    """
    Generate Figure 2: per-patient plots + summary stacked bar chart.
    - Per patient: plots/2/{pid}.png (visual check + 2 histograms)
    - Summary: plots/2/Figure2_summary.png (stacked bar chart, each patient different color)
    """
    data_path = Path(data_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    mat_files = sorted(data_path.glob("*.mat"))
    if not mat_files:
        raise FileNotFoundError(f"No .mat files in {data_dir}")

    # Per-patient data
    per_patient: dict[str, tuple[np.ndarray, np.ndarray]] = {}  # pid -> (lags_eda, lags_rsa)
    saved = []

    for f in mat_files:
        mat = loadmat(str(f))
        t = np.asarray(mat["t"]).flatten()
        eda = np.asarray(mat["eda_clean"] if "eda_clean" in mat else mat["eda"]).flatten()
        rsa = np.asarray(mat["rsa_clean"] if "rsa_clean" in mat else mat["rsa"]).flatten()
        fs = float(mat["FS"].flat[0])
        pid = Path(f).stem.split(" full")[0].replace("_preprocessed", "").strip()

        n = len(t)
        n_plot = min(n, int(max_duration_sec * fs))
        start = (n - n_plot) // 2
        t_win = t[start : start + n_plot]
        eda_win = eda[start : start + n_plot]
        rsa_win = rsa[start : start + n_plot]

        lags_eda, lags_rsa, t_eda, t_rsa = _compute_peak_to_peak_lags(
            t_win, eda_win, rsa_win, fs, top_peak_pct=top_peak_pct
        )
        per_patient[pid] = (lags_eda, lags_rsa)

        # Per-patient figure
        fig = _plot_single_patient(
            t_win, eda_win, rsa_win, t_eda, t_rsa,
            lags_eda, lags_rsa, pid, top_peak_pct, bin_sec,
        )
        pfile = out_path / f"{pid}.png"
        fig.savefig(pfile, dpi=150, bbox_inches="tight")
        plt.close(fig)
        saved.append(pfile)
        print(f"Saved: {pfile}")

    # Summary: stacked bar chart (each bin stacked by patient)
    # Use gradient: orange shades for EDA, blue shades for RSA (cleaner than many distinct colors)
    pids = sorted(per_patient.keys())
    n = max(len(pids), 1)
    colors_eda = plt.colormaps["Oranges"](np.linspace(0.4, 0.95, n))
    colors_rsa = plt.colormaps["Blues"](np.linspace(0.4, 0.95, n))

    all_eda = np.concatenate([per_patient[p][0] for p in pids if len(per_patient[p][0]) > 0])
    all_rsa = np.concatenate([per_patient[p][1] for p in pids if len(per_patient[p][1]) > 0])
    max_eda = max(all_eda.max(), bin_sec) if len(all_eda) > 0 else bin_sec
    max_rsa = max(all_rsa.max(), bin_sec) if len(all_rsa) > 0 else bin_sec
    bins_eda = np.arange(0, max_eda + bin_sec, bin_sec)
    bins_rsa = np.arange(0, max_rsa + bin_sec, bin_sec)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

    # EDA leads: stacked
    if len(all_eda) > 0:
        bin_centers = (bins_eda[:-1] + bins_eda[1:]) / 2
        bottom = np.zeros(len(bin_centers))
        for i, pid in enumerate(pids):
            lags = per_patient[pid][0]
            if len(lags) > 0:
                counts, _ = np.histogram(lags, bins=bins_eda)
                ax1.bar(bin_centers, counts, width=bin_sec, bottom=bottom, color=colors_eda[i], label=pid, edgecolor="white", linewidth=0.3)
                bottom = bottom + counts
        ax1.set_xlabel("Lag (s) — RSA peak after EDA peak")
        ax1.set_ylabel("Count")
        ax1.set_title("EDA leads (stacked by patient)")
        ax1.legend(loc="upper right", ncol=2, fontsize=7)
    ax1.grid(True, alpha=0.3)

    # RSA leads: stacked
    if len(all_rsa) > 0:
        bin_centers = (bins_rsa[:-1] + bins_rsa[1:]) / 2
        bottom = np.zeros(len(bin_centers))
        for i, pid in enumerate(pids):
            lags = per_patient[pid][1]
            if len(lags) > 0:
                counts, _ = np.histogram(lags, bins=bins_rsa)
                ax2.bar(bin_centers, counts, width=bin_sec, bottom=bottom, color=colors_rsa[i], label=pid, edgecolor="white", linewidth=0.3)
                bottom = bottom + counts
        ax2.set_xlabel("Lag (s) — EDA peak after RSA peak")
        ax2.set_ylabel("Count")
        ax2.set_title("RSA leads (stacked by patient)")
        ax2.legend(loc="upper right", ncol=2, fontsize=7)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    summary_file = out_path / "Figure2_summary.png"
    fig.savefig(summary_file, dpi=150, bbox_inches="tight")
    plt.close(fig)
    saved.append(summary_file)
    print(f"Saved: {summary_file}")
    return saved
