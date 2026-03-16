"""
Preprocessing module for EDA and RSA physiological signals.

PREPROCESSING DISABLED: Raw data is clean enough; eda_clean/rsa_clean = passthrough.
(Original: EDA highpass 0.05 Hz, RSA bandpass 0.12–0.4 Hz — see commented code.)

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


def _butter_highpass(cutoff: float, fs: float, order: int = 4):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    return butter(order, normal_cutoff, btype="high", output="sos")


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


def preprocess_eda(
    eda: np.ndarray,
    sampling_rate: float,
    mode: str = "highpass",
    cutoff_hz: float = 0.05,
) -> np.ndarray:
    """
    Preprocess EDA. Based on spectrum: 90% power < 0.15 Hz.

    mode="highpass" (default): 0.05 Hz highpass to remove slow drift,
        extract phasic (SCR-like) component — visible difference from raw.
    mode="lowpass": lowpass at cutoff_hz to remove high-freq noise.

    PREPROCESSING DISABLED: Raw data is clean enough; passthrough.
    """
    eda = np.asarray(eda, dtype=float).flatten()
    # Preprocessing disabled — raw data is clean enough
    return eda.copy()
    # if mode == "highpass":
    #     if HAS_NEUROKIT:
    #         return nk.signal_filter(eda, sampling_rate=sampling_rate, lowcut=cutoff_hz, order=4)
    #     sos = _butter_highpass(cutoff_hz, sampling_rate, order=4)
    #     return sosfiltfilt(sos, eda)
    # else:
    #     if HAS_NEUROKIT:
    #         return nk.signal_filter(eda, sampling_rate=sampling_rate, highcut=cutoff_hz, order=4)
    #     sos = _butter_lowpass(cutoff_hz, sampling_rate, order=4)
    #     return sosfiltfilt(sos, eda)


def preprocess_rsa(
    rsa: np.ndarray,
    sampling_rate: float,
    lowcut: float = 0.12,
    highcut: float = 0.4,
) -> np.ndarray:
    """
    Clean RSA waveform: bandpass filter 0.12–0.4 Hz (respiratory range).
    Removes DC drift and high-frequency noise. Pattern is preserved.

    PREPROCESSING DISABLED: Raw data is clean enough; passthrough.
    """
    rsa = np.asarray(rsa, dtype=float).flatten()
    # Preprocessing disabled — raw data is clean enough
    return rsa.copy()
    # if HAS_NEUROKIT:
    #     return nk.signal_filter(rsa, sampling_rate=sampling_rate, lowcut=lowcut, highcut=highcut, order=4)
    # sos = _butter_bandpass(lowcut, highcut, sampling_rate, order=4)
    # return sosfiltfilt(sos, rsa)


def preprocess_file(
    filepath: str,
    eda_mode: str = "highpass",
    eda_cutoff_hz: float = 0.05,
    rsa_lowcut: float = 0.12,
    rsa_highcut: float = 0.4,
) -> dict:
    """
    Load and preprocess a single .mat file.
    Returns dict with raw and cleaned eda, rsa, t, FS.
    """
    data = load_mat_file(filepath)
    data["eda_clean"] = preprocess_eda(data["eda"], data["FS"], mode=eda_mode, cutoff_hz=eda_cutoff_hz)
    data["rsa_clean"] = preprocess_rsa(data["rsa"], data["FS"], lowcut=rsa_lowcut, highcut=rsa_highcut)
    data["file"] = Path(filepath).name
    return data


def plot_preprocessed(
    data: dict,
    max_duration_sec: float | None = None,
    figsize: tuple = (14, 10),
) -> "matplotlib.figure.Figure":
    """
    Plot raw and preprocessed EDA and RSA in 4 panels (full timecourse by default).
    If max_duration_sec is set, only plot the first N seconds.
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

    fig, axes = plt.subplots(4, 1, figsize=figsize, sharex=True)

    # EDA Raw
    axes[0].plot(t_plot, eda_raw, color="tab:gray", linewidth=0.6)
    axes[0].set_ylabel("EDA")
    axes[0].set_title(f"EDA Raw — {data['file']}{duration_label}")
    axes[0].grid(True, alpha=0.3)

    # EDA Preprocessed
    axes[1].plot(t_plot, eda_clean, color="tab:blue", linewidth=0.6)
    axes[1].set_ylabel("EDA")
    axes[1].set_title(f"EDA Preprocessed — {data['file']}{duration_label}")
    axes[1].grid(True, alpha=0.3)

    # RSA Raw
    axes[2].plot(t_plot, rsa_raw, color="tab:gray", linewidth=0.6)
    axes[2].set_ylabel("RSA")
    axes[2].set_title(f"RSA Raw — {data['file']}{duration_label}")
    axes[2].grid(True, alpha=0.3)

    # RSA Preprocessed
    axes[3].plot(t_plot, rsa_clean, color="tab:orange", linewidth=0.6)
    axes[3].set_ylabel("RSA")
    axes[3].set_xlabel("Time (s)")
    axes[3].set_title(f"RSA Preprocessed — {data['file']}{duration_label}")
    axes[3].grid(True, alpha=0.3)

    fig.tight_layout()
    return fig


def save_preprocessed_data(data: dict, filepath: str, format: str = "mat") -> None:
    """Save preprocessed data to disk. format: 'mat' or 'npy'."""
    from scipy.io import savemat
    out = Path(filepath)
    out.parent.mkdir(parents=True, exist_ok=True)
    if format == "mat":
        savemat(
            str(out),
            {
                "t": data["t"],
                "FS": data["FS"],
                "eda": data["eda"],
                "eda_clean": data["eda_clean"],
                "rsa": data["rsa"],
                "rsa_clean": data["rsa_clean"],
            },
        )
    else:
        np.savez(
            str(out),
            t=data["t"],
            FS=data["FS"],
            eda=data["eda"],
            eda_clean=data["eda_clean"],
            rsa=data["rsa"],
            rsa_clean=data["rsa_clean"],
        )
    print(f"  Saved: {out}")


def preprocess_and_save_plots(
    data_dir: str = "data",
    output_dir: str = "plots/preprocessed",
    data_dir_out: str | None = None,
    pattern: str = "*.mat",
    max_duration_sec: float | None = None,
    eda_mode: str = "highpass",
    eda_cutoff_hz: float = 0.05,
    save_data: bool = False,
    data_format: str = "mat",
) -> list[dict]:
    """
    Preprocess all .mat files, save comparison plots, and return preprocessed data.
    If save_data=True, save preprocessed signals to data_dir_out (default: data/preprocessed).
    """
    data_path = Path(data_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    data_out = Path(data_dir_out) if data_dir_out else data_path / "preprocessed"
    if save_data:
        data_out.mkdir(parents=True, exist_ok=True)

    results = []
    for f in sorted(data_path.glob(pattern)):
        try:
            data = preprocess_file(str(f), eda_mode=eda_mode, eda_cutoff_hz=eda_cutoff_hz)
            fig = plot_preprocessed(data, max_duration_sec=max_duration_sec)
            out_name = Path(f).stem + "_preprocessed.png"
            fig.savefig(output_path / out_name, dpi=150, bbox_inches="tight")
            import matplotlib.pyplot as plt
            plt.close(fig)
            results.append(data)
            print(f"  Saved: {output_path / out_name}")
            if save_data:
                save_preprocessed_data(
                    data,
                    data_out / (Path(f).stem + "_preprocessed." + ("mat" if data_format == "mat" else "npz")),
                    format=data_format,
                )
        except Exception as e:
            print(f"Error processing {f}: {e}")

    return results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Preprocess EDA and RSA, save comparison plots")
    parser.add_argument("--data-dir", default="data", help="Directory with .mat files")
    parser.add_argument("--output-dir", default="plots/preprocessed", help="Directory to save plots")
    parser.add_argument("--max-duration", type=float, default=0, help="Max seconds to plot. 0 = full timecourse (default)")
    parser.add_argument("--eda-mode", choices=["highpass", "lowpass"], default="highpass",
                        help="EDA: highpass (phasic, default) or lowpass (smooth)")
    parser.add_argument("--eda-cutoff", type=float, default=0.05,
                        help="EDA cutoff in Hz. Highpass default 0.05; use 1-3 for lowpass")
    parser.add_argument("--save-data", action="store_true", help="Save preprocessed data to disk")
    parser.add_argument("--data-dir-out", default=None,
                        help="Output dir for preprocessed data (default: data/preprocessed)")
    parser.add_argument("--data-format", choices=["mat", "npy"], default="mat",
                        help="Format for saved data: mat or npy")
    args = parser.parse_args()

    max_dur = None if args.max_duration <= 0 else args.max_duration
    print(f"Preprocessing from {args.data_dir} -> {args.output_dir} (EDA: {args.eda_mode} {args.eda_cutoff} Hz)")
    preprocess_and_save_plots(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        data_dir_out=args.data_dir_out,
        max_duration_sec=max_dur,
        eda_mode=args.eda_mode,
        eda_cutoff_hz=args.eda_cutoff,
        save_data=args.save_data,
        data_format=args.data_format,
    )
