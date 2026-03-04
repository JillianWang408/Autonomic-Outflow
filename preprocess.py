"""
Preprocessing module for EDA and RSA physiological signals.

Applies standard filters:
- EDA: Low-pass 3 Hz (NeuroKit2)
- RSA: Bandpass 0.12–0.4 Hz (respiratory range)

Saves comparison plots of raw vs preprocessed data.
"""

import numpy as np
from pathlib import Path
from scipy.io import loadmat
from scipy.signal import butter, sosfiltfilt

try:
    import neurokit2 as nk
    HAS_NEUROKIT = True
except ImportError:
    HAS_NEUROKIT = False


def _butter_lowpass(cutoff: float, fs: float, order: int = 4):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    return butter(order, normal_cutoff, btype="low", output="sos")


def _butter_bandpass(lowcut: float, highcut: float, fs: float, order: int = 4):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    return butter(order, [low, high], btype="band", output="sos")


def load_mat_file(filepath: str) -> dict:
    """Load .mat file and return dict with t, FS, eda, rsa (flattened)."""
    mat = loadmat(filepath)
    return {
        "t": np.asarray(mat["t"]).flatten(),
        "FS": float(mat["FS"].flat[0]),
        "eda": np.asarray(mat["eda"]).flatten(),
        "rsa": np.asarray(mat["rsa"]).flatten(),
    }


def preprocess_eda(eda: np.ndarray, sampling_rate: float, method: str = "neurokit") -> np.ndarray:
    """
    Clean EDA signal: low-pass filter 3 Hz (removes high-freq noise).
    Uses NeuroKit2 if available, else scipy Butterworth.
    """
    eda = np.asarray(eda, dtype=float).flatten()
    if HAS_NEUROKIT:
        return nk.eda_clean(eda, sampling_rate=sampling_rate, method=method)
    # Fallback: 4th order Butterworth lowpass at 3 Hz
    sos = _butter_lowpass(3.0, sampling_rate, order=4)
    return sosfiltfilt(sos, eda)


def preprocess_rsa(rsa: np.ndarray, sampling_rate: float, lowcut: float = 0.12, highcut: float = 0.4) -> np.ndarray:
    """
    Clean RSA waveform: bandpass filter 0.12–0.4 Hz (respiratory range).
    Removes DC drift and high-frequency noise.
    Uses NeuroKit2 if available, else scipy Butterworth.
    """
    rsa = np.asarray(rsa, dtype=float).flatten()
    if HAS_NEUROKIT:
        return nk.signal_filter(rsa, sampling_rate=sampling_rate, lowcut=lowcut, highcut=highcut, order=4)
    # Fallback: 4th order Butterworth bandpass
    sos = _butter_bandpass(lowcut, highcut, sampling_rate, order=4)
    return sosfiltfilt(sos, rsa)


def preprocess_file(
    filepath: str,
    eda_method: str = "neurokit",
    rsa_lowcut: float = 0.12,
    rsa_highcut: float = 0.4,
) -> dict:
    """
    Load and preprocess a single .mat file.
    Returns dict with raw and cleaned eda, rsa, t, FS.
    """
    data = load_mat_file(filepath)
    data["eda_clean"] = preprocess_eda(data["eda"], data["FS"], method=eda_method)
    data["rsa_clean"] = preprocess_rsa(data["rsa"], data["FS"], lowcut=rsa_lowcut, highcut=rsa_highcut)
    data["file"] = Path(filepath).name
    return data


def plot_preprocessed(
    data: dict,
    max_duration_sec: float | None = 300,
    figsize: tuple = (12, 8),
) -> "matplotlib.figure.Figure":
    """
    Plot raw vs preprocessed EDA and RSA.
    If max_duration_sec is set, only plot the first N seconds (for readability).
    """
    import matplotlib.pyplot as plt

    t = data["t"]
    fs = data["FS"]
    n = len(t)

    if max_duration_sec is not None:
        n_plot = min(n, int(max_duration_sec * fs))
        t_plot = t[:n_plot]
        eda_raw = data["eda"][:n_plot]
        eda_clean = data["eda_clean"][:n_plot]
        rsa_raw = data["rsa"][:n_plot]
        rsa_clean = data["rsa_clean"][:n_plot]
        duration_label = f" (first {max_duration_sec/60:.0f} min)"
    else:
        t_plot = t
        eda_raw = data["eda"]
        eda_clean = data["eda_clean"]
        rsa_raw = data["rsa"]
        rsa_clean = data["rsa_clean"]
        duration_label = ""

    fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=True)

    # EDA
    ax1 = axes[0]
    ax1.plot(t_plot, eda_raw, alpha=0.6, label="Raw", color="gray", linewidth=0.8)
    ax1.plot(t_plot, eda_clean, label="Preprocessed", color="tab:blue", linewidth=1)
    ax1.set_ylabel("EDA")
    ax1.set_title(f"EDA — {data['file']}{duration_label}")
    ax1.legend(loc="upper right")
    ax1.grid(True, alpha=0.3)

    # RSA
    ax2 = axes[1]
    ax2.plot(t_plot, rsa_raw, alpha=0.6, label="Raw", color="gray", linewidth=0.8)
    ax2.plot(t_plot, rsa_clean, label="Preprocessed", color="tab:orange", linewidth=1)
    ax2.set_ylabel("RSA")
    ax2.set_xlabel("Time (s)")
    ax2.set_title(f"RSA — {data['file']}{duration_label}")
    ax2.legend(loc="upper right")
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    return fig


def preprocess_and_save_plots(
    data_dir: str = "data",
    output_dir: str = "plots/preprocessed",
    pattern: str = "*.mat",
    max_duration_sec: float | None = 300,
) -> list[dict]:
    """
    Preprocess all .mat files, save comparison plots, and return preprocessed data.
    """
    data_path = Path(data_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    results = []
    for f in sorted(data_path.glob(pattern)):
        try:
            data = preprocess_file(str(f))
            fig = plot_preprocessed(data, max_duration_sec=max_duration_sec)
            out_name = Path(f).stem + "_preprocessed.png"
            out_file = output_path / out_name
            fig.savefig(out_file, dpi=150, bbox_inches="tight")
            import matplotlib.pyplot as plt
            plt.close(fig)
            results.append(data)
            print(f"  Saved: {out_file}")
        except Exception as e:
            print(f"Error processing {f}: {e}")

    return results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Preprocess EDA and RSA, save comparison plots")
    parser.add_argument("--data-dir", default="data", help="Directory with .mat files")
    parser.add_argument("--output-dir", default="plots/preprocessed", help="Directory to save plots")
    parser.add_argument("--max-duration", type=float, default=300, help="Max seconds to plot (default 300 = 5 min). 0 = full")
    args = parser.parse_args()

    max_dur = None if args.max_duration <= 0 else args.max_duration
    print(f"Preprocessing from {args.data_dir} -> {args.output_dir}")
    preprocess_and_save_plots(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        max_duration_sec=max_dur,
    )
