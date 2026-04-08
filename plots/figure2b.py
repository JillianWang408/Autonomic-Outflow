"""
Figure 2B: Lead–lag within co‑activation episodes (Figure 1–style top‑% active, both channels).

Option 2: contiguous co‑active runs where **both** signals are active; keep a run only if it
contains ≥1 EDA peak and ≥1 RSA peak (same peak detection as 2A). Pair within each segment only.

Each raw co‑active run is extended by **±co_extend_sec** (default 3 s) at both ends (clamped to the
recording) for peak picking, pairing, and colored trace highlights.

Lag histograms use the **full** recording (all qualifying segments). Top panel: middle time window
with **extended** co‑active intervals drawn in color (EDA yellow, RSA blue) and peaks used for 2B in the window.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import loadmat

from plots.figure2 import (
    _find_peaks,
    _compute_lags_eda_anchor,
    _compute_lags_rsa_anchor,
    _hist_lag,
    _stacked_bar_lags,
)
from plots.regions import COL_EDA, COL_RSA, _to_binary_active


def _coactive_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """Contiguous True runs; returns half-open index ranges [s, e) into mask."""
    pad = np.concatenate(([False], mask, [False]))
    d = np.diff(pad.astype(int))
    starts = np.where(d == 1)[0]
    ends = np.where(d == -1)[0]
    return list(zip(starts.tolist(), ends.tolist()))


def _extended_coactive_time_mask(
    t: np.ndarray,
    co: np.ndarray,
    co_extend_sec: float,
) -> np.ndarray:
    """True where t lies inside [run_start - ext, run_end + ext] for any contiguous co-active run."""
    t = np.asarray(t, dtype=float).flatten()
    if co_extend_sec <= 0:
        return np.asarray(co, dtype=bool).copy()
    t_lo_rec, t_hi_rec = float(t[0]), float(t[-1])
    out = np.zeros(len(t), dtype=bool)
    for s, e in _coactive_runs(np.asarray(co, dtype=bool)):
        t0, t1 = float(t[s]), float(t[e - 1])
        lo = max(t_lo_rec, t0 - co_extend_sec)
        hi = min(t_hi_rec, t1 + co_extend_sec)
        out |= (t >= lo) & (t <= hi)
    return out


def _compute_lags_co_segments(
    t: np.ndarray,
    eda: np.ndarray,
    rsa: np.ndarray,
    fs: float,
    co_threshold_pct: float,
    top_peak_pct: float,
    co_extend_sec: float = 3.0,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    int,
]:
    """
    Returns (e_eda, e_rsa, r_rsa, r_eda, t_eda_used, t_rsa_used, n_qualifying_segments)
    where lag arrays concatenate within-segment pairings only; peak arrays list all peak
    times that lie in some qualifying segment (for visualization filter).
    """
    eda_a = _to_binary_active(eda, co_threshold_pct)
    rsa_a = _to_binary_active(rsa, co_threshold_pct)
    co = eda_a & rsa_a

    _, t_eda_all = _find_peaks(eda, t, fs, top_pct=top_peak_pct)
    _, t_rsa_all = _find_peaks(rsa, t, fs, top_pct=top_peak_pct)

    e_eda_list: list[np.ndarray] = []
    e_rsa_list: list[np.ndarray] = []
    r_rsa_list: list[np.ndarray] = []
    r_eda_list: list[np.ndarray] = []
    t_eda_used_list: list[float] = []
    t_rsa_used_list: list[float] = []
    n_qual = 0

    t_lo_rec, t_hi_rec = float(t[0]), float(t[-1])
    for s, e in _coactive_runs(co):
        t0, t1 = float(t[s]), float(t[e - 1])
        lo = max(t_lo_rec, t0 - co_extend_sec)
        hi = min(t_hi_rec, t1 + co_extend_sec)
        te = t_eda_all[(t_eda_all >= lo) & (t_eda_all <= hi)]
        tr = t_rsa_all[(t_rsa_all >= lo) & (t_rsa_all <= hi)]
        if len(te) < 1 or len(tr) < 1:
            continue
        n_qual += 1
        a1, a2 = _compute_lags_eda_anchor(te, tr)
        b1, b2 = _compute_lags_rsa_anchor(te, tr)
        e_eda_list.append(a1)
        e_rsa_list.append(a2)
        r_rsa_list.append(b1)
        r_eda_list.append(b2)
        t_eda_used_list.extend(te.tolist())
        t_rsa_used_list.extend(tr.tolist())

    def _cat(xs: list[np.ndarray]) -> np.ndarray:
        return np.concatenate(xs) if xs else np.array([])

    t_eda_used = np.unique(np.array(t_eda_used_list, dtype=float))
    t_rsa_used = np.unique(np.array(t_rsa_used_list, dtype=float))
    return (
        _cat(e_eda_list),
        _cat(e_rsa_list),
        _cat(r_rsa_list),
        _cat(r_eda_list),
        t_eda_used,
        t_rsa_used,
        n_qual,
    )


def _plot_single_patient_2b(
    t_win: np.ndarray,
    eda_win: np.ndarray,
    rsa_win: np.ndarray,
    co_mask_win: np.ndarray,
    t_eda_vis: np.ndarray,
    t_rsa_vis: np.ndarray,
    eda_anchor_eda_leads: np.ndarray,
    eda_anchor_rsa_leads: np.ndarray,
    rsa_anchor_rsa_leads: np.ndarray,
    rsa_anchor_eda_leads: np.ndarray,
    pid: str,
    co_threshold_pct: float,
    top_peak_pct: float,
    n_segments_full: int,
    bin_sec: float,
    co_extend_sec: float,
) -> plt.Figure:
    """Base traces + colored co‑active intervals; peaks in qualifying episodes (window)."""
    fig = plt.figure(figsize=(12, 14))

    ax0 = fig.add_subplot(5, 1, 1)
    eda_z = (eda_win - eda_win.mean()) / (eda_win.std() or 1)
    rsa_z = (rsa_win - rsa_win.mean()) / (rsa_win.std() or 1)
    ax0.plot(t_win, eda_z, color="gray", linewidth=0.8, alpha=0.8, label="EDA")
    ax0.plot(t_win, rsa_z, color="black", linewidth=0.8, alpha=0.8, label="RSA")

    pad = np.concatenate(([False], co_mask_win, [False]))
    d = np.diff(pad.astype(int))
    starts = np.where(d == 1)[0]
    ends = np.where(d == -1)[0]
    for s, e in zip(starts, ends):
        ax0.plot(t_win[s:e], eda_z[s:e], color=COL_EDA, linewidth=1.5, alpha=0.95, zorder=4)
        ax0.plot(t_win[s:e], rsa_z[s:e], color=COL_RSA, linewidth=1.3, alpha=0.95, zorder=4)

    if len(t_eda_vis) > 0:
        ax0.scatter(
            t_eda_vis, np.interp(t_eda_vis, t_win, eda_z),
            color="darkorange", s=30, zorder=5, label="EDA peaks (2B, in window)",
        )
    if len(t_rsa_vis) > 0:
        ax0.scatter(
            t_rsa_vis, np.interp(t_rsa_vis, t_win, rsa_z),
            color="navy", s=30, zorder=5, label="RSA peaks (2B, in window)",
        )

    n_eda = len(eda_anchor_eda_leads) + len(eda_anchor_rsa_leads)
    n_rsa = len(rsa_anchor_rsa_leads) + len(rsa_anchor_eda_leads)
    ax0.set_ylabel("Z-score")
    ax0.set_xlabel("Time (s)")
    ax0.set_title(
        f"{pid} — co‑active top {co_threshold_pct:.0f}% (both), ±{co_extend_sec:.0f}s pad, ≥1 peak each (extended) | "
        f"{n_segments_full} qualifying segments (full) | trace: window; "
        f"histograms: {n_eda} EDA‑anchor, {n_rsa} RSA‑anchor lags"
    )
    ax0.legend(loc="upper right", ncol=2, fontsize=8)
    ax0.grid(True, alpha=0.3)

    ax1 = fig.add_subplot(5, 1, 2)
    _hist_lag(
        ax1, eda_anchor_eda_leads, bin_sec, "orange", "darkorange",
        "Lag (s): RSA peak after EDA peak",
        "Anchor = EDA (within co‑active segment): EDA leads",
    )

    ax2 = fig.add_subplot(5, 1, 3)
    _hist_lag(
        ax2, eda_anchor_rsa_leads, bin_sec, "steelblue", "navy",
        "Lag (s): RSA peak before EDA peak",
        "Anchor = EDA: RSA leads",
    )

    ax3 = fig.add_subplot(5, 1, 4)
    _hist_lag(
        ax3, rsa_anchor_rsa_leads, bin_sec, "steelblue", "navy",
        "Lag (s): EDA peak after RSA peak",
        "Anchor = RSA (within co‑active segment): RSA leads",
    )

    ax4 = fig.add_subplot(5, 1, 5)
    _hist_lag(
        ax4, rsa_anchor_eda_leads, bin_sec, "orange", "darkorange",
        "Lag (s): EDA peak before RSA peak",
        "Anchor = RSA: EDA leads",
    )

    fig.tight_layout()
    return fig


def plot_figure2b(
    data_dir: str = "data",
    output_dir: str = "plots/2B",
    max_duration_sec: float | None = 300,
    bin_sec: float = 2,
    top_peak_pct: float = 30,
    co_threshold_pct: float = 40,
    co_extend_sec: float = 3.0,
) -> list[Path]:
    """
    Within-segment lags only (Option 2). Active = Figure 1 top co_threshold_pct% per channel.
    Each co-active run is extended by co_extend_sec before and after (clamped) for peaks and color.

    Summary file: Figure2B_summary.png
    """
    data_path = Path(data_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    mat_files = sorted(data_path.glob("*.mat"))
    if not mat_files:
        raise FileNotFoundError(f"No .mat files in {data_dir}")

    per_patient: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
    saved: list[Path] = []

    for f in mat_files:
        mat = loadmat(str(f))
        t = np.asarray(mat["t"]).flatten()
        eda = np.asarray(mat["eda_clean"] if "eda_clean" in mat else mat["eda"]).flatten()
        rsa = np.asarray(mat["rsa_clean"] if "rsa_clean" in mat else mat["rsa"]).flatten()
        fs = float(mat["FS"].flat[0])
        pid = Path(f).stem.split(" full")[0].replace("_preprocessed", "").strip()
        n = len(t)

        eda_a_full = _to_binary_active(eda, co_threshold_pct)
        rsa_a_full = _to_binary_active(rsa, co_threshold_pct)
        co_full = eda_a_full & rsa_a_full

        e_eda, e_rsa, r_rsa, r_eda, t_eda_ep, t_rsa_ep, n_seg = _compute_lags_co_segments(
            t, eda, rsa, fs,
            co_threshold_pct=co_threshold_pct,
            top_peak_pct=top_peak_pct,
            co_extend_sec=co_extend_sec,
        )
        per_patient[pid] = (e_eda, e_rsa, r_rsa, r_eda)

        ext_mask_full = _extended_coactive_time_mask(t, co_full, co_extend_sec)

        if max_duration_sec is not None and max_duration_sec > 0:
            n_plot = min(n, int(max_duration_sec * fs))
            start = (n - n_plot) // 2
            t_win = t[start : start + n_plot]
            eda_win = eda[start : start + n_plot]
            rsa_win = rsa[start : start + n_plot]
            co_win = ext_mask_full[start : start + n_plot]
        else:
            t_win, eda_win, rsa_win = t, eda, rsa
            co_win = ext_mask_full

        t_lo, t_hi = float(t_win[0]), float(t_win[-1])
        t_eda_vis = t_eda_ep[(t_eda_ep >= t_lo) & (t_eda_ep <= t_hi)]
        t_rsa_vis = t_rsa_ep[(t_rsa_ep >= t_lo) & (t_rsa_ep <= t_hi)]

        fig = _plot_single_patient_2b(
            t_win, eda_win, rsa_win, co_win, t_eda_vis, t_rsa_vis,
            e_eda, e_rsa, r_rsa, r_eda,
            pid, co_threshold_pct, top_peak_pct, n_seg, bin_sec, co_extend_sec,
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

    _stacked_bar_lags(
        ax11, per_patient, pids, 0, b0, bin_sec, colors_eda,
        "Anchor = EDA: EDA leads (within co‑active segment)",
        "Lag (s): RSA after EDA",
    )
    _stacked_bar_lags(
        ax12, per_patient, pids, 1, b1, bin_sec, colors_rsa,
        "Anchor = EDA: RSA leads",
        "Lag (s): RSA before EDA",
    )
    _stacked_bar_lags(
        ax21, per_patient, pids, 2, b2, bin_sec, colors_rsa,
        "Anchor = RSA: RSA leads",
        "Lag (s): EDA after RSA",
    )
    _stacked_bar_lags(
        ax22, per_patient, pids, 3, b3, bin_sec, colors_eda,
        "Anchor = RSA: EDA leads",
        "Lag (s): EDA before RSA",
    )

    fig.suptitle(
        "Figure 2B: co‑activation (+padding), ≥1 peak each in extended window, within‑segment lead–lag",
        fontsize=12, y=1.02,
    )
    fig.tight_layout()
    summary_file = out_path / "Figure2B_summary.png"
    fig.savefig(summary_file, dpi=150, bbox_inches="tight")
    plt.close(fig)
    saved.append(summary_file)
    print(f"Saved: {summary_file}")
    return saved
