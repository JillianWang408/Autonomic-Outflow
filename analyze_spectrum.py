"""
Analyze frequency composition of EDA and RSA signals.

Uses Welch's method for power spectral density. Saves spectrum plots
to help choose appropriate filter cutoffs.
"""

import numpy as np
from pathlib import Path
from scipy.io import loadmat
from scipy.signal import welch


def load_mat_file(filepath: str) -> dict:
    """Load .mat file and return dict with t, FS, eda, rsa (flattened)."""
    mat = loadmat(filepath)
    return {
        "t": np.asarray(mat["t"]).flatten(),
        "FS": float(mat["FS"].flat[0]),
        "eda": np.asarray(mat["eda"]).flatten(),
        "rsa": np.asarray(mat["rsa"]).flatten(),
    }


def compute_spectrum(signal: np.ndarray, fs: float, nperseg: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute power spectral density using Welch's method.
    Returns frequencies (Hz) and power (V^2/Hz or equivalent).
    """
    if nperseg is None:
        nperseg = min(2**14, len(signal) // 4)  # ~16k or 1/4 of signal
    f, pxx = welch(signal, fs=fs, nperseg=nperseg)
    return f, pxx


def analyze_and_plot(
    filepath: str,
    output_dir: str = "plots/spectrum",
    max_freq_hz: float = 10.0,
) -> dict:
    """
    Analyze EDA and RSA frequency content, save spectrum plots.
    """
    import matplotlib.pyplot as plt

    data = load_mat_file(filepath)
    fs = data["FS"]
    fname = Path(filepath).stem

    # EDA spectrum
    f_eda, pxx_eda = compute_spectrum(data["eda"], fs)
    mask_eda = f_eda <= max_freq_hz
    f_eda_plot = f_eda[mask_eda]
    pxx_eda_plot = pxx_eda[mask_eda]

    # RSA spectrum
    f_rsa, pxx_rsa = compute_spectrum(data["rsa"], fs)
    mask_rsa = f_rsa <= max_freq_hz
    f_rsa_plot = f_rsa[mask_rsa]
    pxx_rsa_plot = pxx_rsa[mask_rsa]

    # Cumulative power (normalized 0-1)
    cumsum_eda = np.cumsum(pxx_eda_plot)
    cumsum_eda = cumsum_eda / cumsum_eda[-1] if cumsum_eda[-1] > 0 else cumsum_eda
    cumsum_rsa = np.cumsum(pxx_rsa_plot)
    cumsum_rsa = cumsum_rsa / cumsum_rsa[-1] if cumsum_rsa[-1] > 0 else cumsum_rsa

    # Find frequency bands containing 50%, 90%, 95%, 99% of power
    def freq_at_cumulative(cumsum, freqs, pct):
        idx = np.searchsorted(cumsum, pct / 100)
        return freqs[min(idx, len(freqs) - 1)]

    eda_50 = freq_at_cumulative(cumsum_eda, f_eda_plot, 50)
    eda_90 = freq_at_cumulative(cumsum_eda, f_eda_plot, 90)
    eda_95 = freq_at_cumulative(cumsum_eda, f_eda_plot, 95)
    eda_99 = freq_at_cumulative(cumsum_eda, f_eda_plot, 99)
    rsa_50 = freq_at_cumulative(cumsum_rsa, f_rsa_plot, 50)
    rsa_90 = freq_at_cumulative(cumsum_rsa, f_rsa_plot, 90)

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # EDA power spectrum
    ax = axes[0, 0]
    ax.semilogy(f_eda_plot, pxx_eda_plot + 1e-20, color="tab:blue")
    ax.axvline(eda_90, color="gray", linestyle="--", alpha=0.7, label=f"90% power < {eda_90:.3f} Hz")
    ax.axvline(eda_99, color="gray", linestyle=":", alpha=0.7, label=f"99% power < {eda_99:.3f} Hz")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Power")
    ax.set_title(f"EDA Power Spectrum — {fname}")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, max_freq_hz)

    # EDA cumulative power
    ax = axes[0, 1]
    ax.plot(f_eda_plot, cumsum_eda * 100, color="tab:blue")
    ax.axhline(50, color="gray", linestyle="--", alpha=0.5)
    ax.axhline(90, color="gray", linestyle="--", alpha=0.5)
    ax.axhline(99, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Cumulative power (%)")
    ax.set_title(f"EDA Cumulative Power — 50%<{eda_50:.3f}Hz, 90%<{eda_90:.3f}Hz, 99%<{eda_99:.3f}Hz")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, max_freq_hz)
    ax.set_ylim(0, 105)

    # RSA power spectrum
    ax = axes[1, 0]
    ax.semilogy(f_rsa_plot, pxx_rsa_plot + 1e-20, color="tab:orange")
    ax.axvline(rsa_90, color="gray", linestyle="--", alpha=0.7, label=f"90% power < {rsa_90:.3f} Hz")
    ax.axvspan(0.12, 0.4, alpha=0.2, color="green", label="RSA band (0.12-0.4 Hz)")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Power")
    ax.set_title(f"RSA Power Spectrum — {fname}")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, max_freq_hz)

    # RSA cumulative power
    ax = axes[1, 1]
    ax.plot(f_rsa_plot, cumsum_rsa * 100, color="tab:orange")
    ax.axvspan(0.12, 0.4, alpha=0.2, color="green")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Cumulative power (%)")
    ax.set_title(f"RSA Cumulative Power — 50%<{rsa_50:.3f}Hz, 90%<{rsa_90:.3f}Hz")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, max_freq_hz)
    ax.set_ylim(0, 105)

    fig.tight_layout()
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    out_file = out_path / f"{fname}_spectrum.png"
    fig.savefig(out_file, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return {
        "file": fname,
        "eda": {"f50": eda_50, "f90": eda_90, "f95": eda_95, "f99": eda_99},
        "rsa": {"f50": rsa_50, "f90": rsa_90},
    }


def analyze_all(
    data_dir: str = "data",
    output_dir: str = "plots/spectrum",
    pattern: str = "*.mat",
    max_freq_hz: float = 10.0,
) -> list[dict]:
    """Analyze all .mat files and print summary."""
    data_path = Path(data_dir)
    results = []
    for f in sorted(data_path.glob(pattern)):
        try:
            r = analyze_and_plot(str(f), output_dir=output_dir, max_freq_hz=max_freq_hz)
            results.append(r)
            print(f"  {r['file']}: EDA 90%<{r['eda']['f90']:.3f}Hz, 99%<{r['eda']['f99']:.3f}Hz | RSA 90%<{r['rsa']['f90']:.3f}Hz")
        except Exception as e:
            print(f"Error: {f}: {e}")
    return results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Analyze EDA/RSA frequency composition")
    parser.add_argument("--data-dir", default="data")
    parser.add_argument("--output-dir", default="plots/spectrum")
    parser.add_argument("--max-freq", type=float, default=10.0, help="Max frequency to plot (Hz)")
    args = parser.parse_args()

    print(f"Analyzing spectrum from {args.data_dir} -> {args.output_dir}\n")
    analyze_all(args.data_dir, args.output_dir, max_freq_hz=args.max_freq)
