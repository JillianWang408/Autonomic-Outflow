"""
Figure 4: Lead–lag from local maxima of |d(signal)/dt| (same pairing logic as Figure 2A).

EDA anchor vs RSA anchor; per-patient + 2×2 stacked summary.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import find_peaks
from scipy.io import loadmat

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from eda_rsa_overlap_slope import abs_slope


def _find_slope_peaks(
    t: np.ndarray,
    eda: np.ndarray,
    rsa: np.ndarray,
    fs: float,
    min_dist_sec: float = 5,
    prominence_pct: float = 5,
    top_pct: float | None = 30,
) -> tuple[np.ndarray, np.ndarray]:
    """Local maxima of |dEDA/dt| and |dRSA/dt|; returns (t_eda_events, t_rsa_events)."""
    sm_eda = abs_slope(eda, t)
    sm_rsa = abs_slope(rsa, t)
    eda_peaks, t_eda = _find_peaks_on_series(sm_eda, t, fs, min_dist_sec, prominence_pct, top_pct)
    rsa_peaks, t_rsa = _find_peaks_on_series(sm_rsa, t, fs, min_dist_sec, prominence_pct, top_pct)
    return t_eda, t_rsa


def _find_peaks_on_series(
    x: np.ndarray,
    t: np.ndarray,
    fs: float,
    min_dist_sec: float,
    prominence_pct: float,
    top_pct: float | None,
) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float).flatten()
    min_dist = int(min_dist_sec * fs)
    prom = np.percentile(np.abs(np.diff(x)), prominence_pct) if prominence_pct > 0 else None
    prom = max(prom, 1e-9) if prom is not None else None
    peaks, _ = find_peaks(x, distance=min_dist, prominence=prom)
    if top_pct is not None and len(peaks) > 0:
        thresh = np.percentile(x[peaks], 100 - top_pct)
        peaks = peaks[x[peaks] >= thresh]
    return peaks, t[peaks]


def _compute_lags_eda_anchor(
    t_eda: np.ndarray,
    t_rsa: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    if len(t_eda) == 0 or len(t_rsa) == 0:
        return np.array([]), np.array([])
    lags_eda_leads = []
    lags_rsa_leads = []
    for te in t_eda:
        idx = np.argmin(np.abs(t_rsa - te))
        lag = t_rsa[idx] - te
        if lag > 0:
            lags_eda_leads.append(lag)
        elif lag < 0:
            lags_rsa_leads.append(-lag)
    return np.array(lags_eda_leads), np.array(lags_rsa_leads)


def _compute_lags_rsa_anchor(
    t_eda: np.ndarray,
    t_rsa: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    if len(t_eda) == 0 or len(t_rsa) == 0:
        return np.array([]), np.array([])
    lags_rsa_leads = []
    lags_eda_leads = []
    for tr in t_rsa:
        idx = np.argmin(np.abs(t_eda - tr))
        lag = t_eda[idx] - tr
        if lag > 0:
            lags_rsa_leads.append(lag)
        elif lag < 0:
            lags_eda_leads.append(-lag)
    return np.array(lags_rsa_leads), np.array(lags_eda_leads)


def _compute_lags_from_slope_peaks(
    t: np.ndarray,
    eda: np.ndarray,
    rsa: np.ndarray,
    fs: float,
    top_peak_pct: float = 30,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    t_eda, t_rsa = _find_slope_peaks(t, eda, rsa, fs, top_pct=top_peak_pct)
    e_eda, e_rsa = _compute_lags_eda_anchor(t_eda, t_rsa)
    r_rsa, r_eda = _compute_lags_rsa_anchor(t_eda, t_rsa)
    return t_eda, t_rsa, e_eda, e_rsa, r_rsa, r_eda


def _hist_lag(ax, lags: np.ndarray, bin_sec: float, color: str, edge: str, xlabel: str, title: str) -> None:
    if len(lags) > 0:
        mx = max(lags.max(), bin_sec)
        bins = np.arange(0, mx + bin_sec, bin_sec)
        ax.hist(lags, bins=bins, color=color, alpha=0.7, edgecolor=edge)
        ax.axvline(np.median(lags), color="red", linestyle="--", linewidth=1.5, label=f"Median = {np.median(lags):.1f} s")
        ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)


def _plot_single_patient(
    t_win: np.ndarray,
    eda_win: np.ndarray,
    rsa_win: np.ndarray,
    t_eda_vis: np.ndarray,
    t_rsa_vis: np.ndarray,
    eda_anchor_eda_leads: np.ndarray,
    eda_anchor_rsa_leads: np.ndarray,
    rsa_anchor_rsa_leads: np.ndarray,
    rsa_anchor_eda_leads: np.ndarray,
    pid: str,
    top_peak_pct: float,
    bin_sec: float,
) -> plt.Figure:
    fig = plt.figure(figsize=(12, 14))
    ax0 = fig.add_subplot(5, 1, 1)
    eda_z = (eda_win - eda_win.mean()) / (eda_win.std() or 1)
    rsa_z = (rsa_win - rsa_win.mean()) / (rsa_win.std() or 1)
    ax0.plot(t_win, eda_z, color="gray", linewidth=0.8, alpha=0.8, label="EDA")
    ax0.plot(t_win, rsa_z, color="black", linewidth=0.8, alpha=0.8, label="RSA")
    if len(t_eda_vis) > 0:
        ax0.scatter(t_eda_vis, np.interp(t_eda_vis, t_win, eda_z), color="orange", s=30, zorder=5, label="EDA max|slope| (in window)")
    if len(t_rsa_vis) > 0:
        ax0.scatter(t_rsa_vis, np.interp(t_rsa_vis, t_win, rsa_z), color="steelblue", s=30, zorder=5, label="RSA max|slope| (in window)")
    ax0.set_ylabel("Z-score")
    ax0.set_xlabel("Time (s)")
    ax0.set_title(
        f"{pid} — top {top_peak_pct:.0f}% max |slope| (Fig 4) | trace: window only; histograms: full recording"
    )
    ax0.legend(loc="upper right", ncol=2, fontsize=8)
    ax0.grid(True, alpha=0.3)

    ax1 = fig.add_subplot(5, 1, 2)
    _hist_lag(ax1, eda_anchor_eda_leads, bin_sec, "orange", "darkorange",
                "Lag (s): RSA event after EDA event",
                "Anchor = EDA max|slope|: EDA leads / RSA follows")
    ax2 = fig.add_subplot(5, 1, 3)
    _hist_lag(ax2, eda_anchor_rsa_leads, bin_sec, "steelblue", "navy",
                "Lag (s): RSA event before EDA event",
                "Anchor = EDA max|slope|: RSA leads / EDA follows")
    ax3 = fig.add_subplot(5, 1, 4)
    _hist_lag(ax3, rsa_anchor_rsa_leads, bin_sec, "steelblue", "navy",
                "Lag (s): EDA event after RSA event",
                "Anchor = RSA max|slope|: RSA leads / EDA follows")
    ax4 = fig.add_subplot(5, 1, 5)
    _hist_lag(ax4, rsa_anchor_eda_leads, bin_sec, "orange", "darkorange",
                "Lag (s): EDA event before RSA event",
                "Anchor = RSA max|slope|: EDA leads / RSA follows")
    fig.tight_layout()
    return fig


def _stacked_bar_lags(
    ax: plt.Axes,
    per_patient: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    pids: list[str],
    lag_index: int,
    bins: np.ndarray,
    bin_sec: float,
    colors: np.ndarray,
    title: str,
    xlabel: str,
) -> None:
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bottom = np.zeros(len(bin_centers))
    for i, pid in enumerate(pids):
        lags = per_patient[pid][lag_index]
        if len(lags) > 0:
            counts, _ = np.histogram(lags, bins=bins)
            ax.bar(bin_centers, counts, width=bin_sec, bottom=bottom, color=colors[i], label=pid, edgecolor="white", linewidth=0.3)
            bottom = bottom + counts
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    ax.set_title(title)
    ax.legend(loc="upper right", ncol=2, fontsize=7)
    ax.grid(True, alpha=0.3)


def plot_figure4(
    data_dir: str = "data",
    output_dir: str = "plots/4",
    max_duration_sec: float | None = 300,
    bin_sec: float = 2,
    top_peak_pct: float = 30,
) -> list[Path]:
    data_path = Path(data_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    mat_files = sorted(data_path.glob("*.mat"))
    if not mat_files:
        raise FileNotFoundError(f"No .mat files in {data_dir}")

    per_patient: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
    saved = []

    for f in mat_files:
        mat = loadmat(str(f))
        t = np.asarray(mat["t"]).flatten()
        eda = np.asarray(mat["eda_clean"] if "eda_clean" in mat else mat["eda"]).flatten()
        rsa = np.asarray(mat["rsa_clean"] if "rsa_clean" in mat else mat["rsa"]).flatten()
        fs = float(mat["FS"].flat[0])
        pid = Path(f).stem.split(" full")[0].replace("_preprocessed", "").strip()
        n = len(t)

        t_eda, t_rsa, e_eda, e_rsa, r_rsa, r_eda = _compute_lags_from_slope_peaks(
            t, eda, rsa, fs, top_peak_pct=top_peak_pct
        )
        per_patient[pid] = (e_eda, e_rsa, r_rsa, r_eda)

        if max_duration_sec is not None and max_duration_sec > 0:
            n_plot = min(n, int(max_duration_sec * fs))
            start = (n - n_plot) // 2
            t_win = t[start : start + n_plot]
            eda_win = eda[start : start + n_plot]
            rsa_win = rsa[start : start + n_plot]
        else:
            t_win, eda_win, rsa_win = t, eda, rsa
        t_lo, t_hi = float(t_win[0]), float(t_win[-1])
        t_eda_vis = t_eda[(t_eda >= t_lo) & (t_eda <= t_hi)]
        t_rsa_vis = t_rsa[(t_rsa >= t_lo) & (t_rsa <= t_hi)]

        fig = _plot_single_patient(
            t_win, eda_win, rsa_win, t_eda_vis, t_rsa_vis,
            e_eda, e_rsa, r_rsa, r_eda,
            pid, top_peak_pct, bin_sec,
        )
        pfile = out_path / f"{pid}.png"
        fig.savefig(pfile, dpi=150, bbox_inches="tight")
        plt.close(fig)
        saved.append(pfile)
        print(f"Saved: {pfile}")

    pids = sorted(per_patient.keys())
    n_pt = max(len(pids), 1)
    colors_eda = plt.colormaps["Oranges"](np.linspace(0.4, 0.95, n_pt))
    colors_rsa = plt.colormaps["Blues"](np.linspace(0.4, 0.95, n_pt))

    def _max_lag(idx: int) -> float:
        arrs = [per_patient[p][idx] for p in pids if len(per_patient[p][idx]) > 0]
        if not arrs:
            return bin_sec
        return max(np.concatenate(arrs).max(), bin_sec)

    b0 = np.arange(0, _max_lag(0) + bin_sec, bin_sec)
    b1 = np.arange(0, _max_lag(1) + bin_sec, bin_sec)
    b2 = np.arange(0, _max_lag(2) + bin_sec, bin_sec)
    b3 = np.arange(0, _max_lag(3) + bin_sec, bin_sec)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax11, ax12 = axes[0]
    ax21, ax22 = axes[1]

    _stacked_bar_lags(ax11, per_patient, pids, 0, b0, bin_sec, colors_eda,
                      "Anchor = EDA max|slope|: EDA leads", "Lag (s): RSA after EDA")
    _stacked_bar_lags(ax12, per_patient, pids, 1, b1, bin_sec, colors_rsa,
                      "Anchor = EDA max|slope|: RSA leads", "Lag (s): RSA before EDA")
    _stacked_bar_lags(ax21, per_patient, pids, 2, b2, bin_sec, colors_rsa,
                      "Anchor = RSA max|slope|: RSA leads", "Lag (s): EDA after RSA")
    _stacked_bar_lags(ax22, per_patient, pids, 3, b3, bin_sec, colors_eda,
                      "Anchor = RSA max|slope|: EDA leads", "Lag (s): EDA before RSA")

    fig.suptitle("Figure 4: lead–lag from max |slope| events (stacked by patient)", fontsize=12, y=1.02)
    fig.tight_layout()
    summary_file = out_path / "Figure4_summary.png"
    fig.savefig(summary_file, dpi=150, bbox_inches="tight")
    plt.close(fig)
    saved.append(summary_file)
    print(f"Saved: {summary_file}")
    return saved
