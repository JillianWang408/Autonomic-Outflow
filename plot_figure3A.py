#!/usr/bin/env python3
"""
Figure 3A: same as 1A, but overlap = top X% of |d(signal)/dt| on each channel.

Usage:
    python plot_figure3A.py
    python plot_figure3A.py --all
    python plot_figure3A.py --patient EC288 --threshold 20
"""

import argparse
from pathlib import Path

from plots.regions_slope import (
    generate_all_threshold_plots_slope,
    pick_best_and_save_figure3a,
    plot_patient_overlap_regions_slope,
)
from scipy.io import loadmat


def main():
    parser = argparse.ArgumentParser(description="Figure 3A — slope-magnitude overlap")
    parser.add_argument("--data-dir", default="data", help="Data directory (default: raw)")
    parser.add_argument("--output-dir", default="plots/3A", help="Output directory")
    parser.add_argument("--max-duration", type=float, default=300, help="Middle N seconds (0 = full)")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--patient", type=str)
    parser.add_argument("--threshold", type=int, default=20)
    parser.add_argument("--thresholds", type=int, nargs="+", default=(5, 10, 15, 20, 25, 30, 35, 40))
    args = parser.parse_args()
    max_dur = None if args.max_duration <= 0 else args.max_duration

    if args.patient:
        data_path = Path(args.data_dir)
        candidates = list(data_path.glob(f"*{args.patient}*.mat"))
        if not candidates:
            print(f"No file matching {args.patient}")
            return
        mat = loadmat(str(candidates[0]))
        t = mat["t"].flatten()
        eda = (mat["eda_clean"] if "eda_clean" in mat else mat["eda"]).flatten()
        rsa = (mat["rsa_clean"] if "rsa_clean" in mat else mat["rsa"]).flatten()
        pid = Path(candidates[0].name).stem.split(" full")[0].replace("_preprocessed", "").strip()
        fig = plot_patient_overlap_regions_slope(
            t, eda, rsa, threshold_pct=args.threshold,
            max_duration_sec=max_dur, participant_id=pid,
        )
        out_dir = Path(args.output_dir)
        (out_dir / f"thresh{args.threshold}").mkdir(parents=True, exist_ok=True)
        out = out_dir / f"thresh{args.threshold}" / f"{pid}.png"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        import matplotlib.pyplot as plt
        plt.close(fig)
        print(f"Saved: {out}")
    elif args.all:
        generate_all_threshold_plots_slope(
            data_dir=args.data_dir,
            thresholds=list(args.thresholds),
            max_duration_sec=max_dur,
            output_dir=args.output_dir,
        )
        print(f"Done: {args.output_dir}/")
    else:
        pick_best_and_save_figure3a(
            data_dir=args.data_dir,
            thresholds=list(args.thresholds),
            max_duration_sec=max_dur,
            output_dir=args.output_dir,
        )


if __name__ == "__main__":
    main()
