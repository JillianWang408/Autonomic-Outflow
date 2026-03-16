#!/usr/bin/env python3
"""
Generate overlap region plots for Figure 1A.

Thresholds: [5, 10, 15, 20, 25, 30, 35, 40]%
Colors: EDA (sympathetic) yellow, RSA (parasympathetic) blue.

Usage:
    python plot_figure1A.py                    # Generate all, pick best, save Figure 1A draft
    python plot_figure1A.py --all             # Generate plots for all patients × thresholds
    python plot_figure1A.py --patient EC288 --threshold 20  # Specific patient/threshold
"""

import argparse
from pathlib import Path

from plots import (
    generate_all_threshold_plots,
    pick_best_and_save_figure1a,
    plot_patient_overlap_regions,
)
from scipy.io import loadmat


def main():
    parser = argparse.ArgumentParser(description="Plot EDA/RSA overlap regions for Figure 1A")
    parser.add_argument("--data-dir", default="data/preprocessed", help="Data directory")
    parser.add_argument("--output-dir", default="plots/1A", help="Output directory (Figure 1A folder)")
    parser.add_argument("--max-duration", type=float, default=300, help="Plot middle N seconds (0 = full time course)")
    parser.add_argument("--all", action="store_true", help="Generate all patient×threshold plots")
    parser.add_argument("--patient", type=str, help="Specific participant ID (e.g. EC288_pN22)")
    parser.add_argument("--threshold", type=int, default=20, help="Threshold percent (default 20)")
    parser.add_argument("--thresholds", type=int, nargs="+", default=(5, 10, 15, 20, 25, 30, 35, 40),
                        help="Threshold range for --all or auto-pick")
    args = parser.parse_args()
    max_dur = None if args.max_duration <= 0 else args.max_duration

    if args.patient:
        # Single patient, single threshold
        data_path = Path(args.data_dir)
        candidates = list(data_path.glob(f"*{args.patient}*.mat"))
        if not candidates:
            print(f"No file matching {args.patient} in {args.data_dir}")
            return
        mat = loadmat(str(candidates[0]))
        t = mat["t"].flatten()
        eda = (mat["eda_clean"] if "eda_clean" in mat else mat["eda"]).flatten()
        rsa = (mat["rsa_clean"] if "rsa_clean" in mat else mat["rsa"]).flatten()
        pid = Path(candidates[0].name).stem.replace(" full EDA & RSA waveform_preprocessed", "")
        fig = plot_patient_overlap_regions(
            t, eda, rsa, threshold_pct=args.threshold,
            max_duration_sec=max_dur,
            participant_id=pid,
        )
        out_dir = Path(args.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / f"thresh{args.threshold}").mkdir(parents=True, exist_ok=True)
        out = out_dir / f"thresh{args.threshold}" / f"{pid}.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        import matplotlib.pyplot as plt
        plt.close(fig)
        print(f"Saved: {out}")

    elif args.all:
        # All patients × all thresholds → plots/1A/thresh{N}/{patient}.png
        results = generate_all_threshold_plots(
            data_dir=args.data_dir,
            thresholds=args.thresholds,
            max_duration_sec=max_dur,
            output_dir=args.output_dir,
        )
        print(f"Generated {len(results)} plots in {args.output_dir}/ (by threshold)")

    else:
        # Auto-pick best and save Figure 1A draft
        pick_best_and_save_figure1a(
            data_dir=args.data_dir,
            thresholds=args.thresholds,
            max_duration_sec=max_dur,
            output_dir=args.output_dir,
        )


if __name__ == "__main__":
    main()
