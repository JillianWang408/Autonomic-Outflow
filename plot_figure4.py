#!/usr/bin/env python3
"""
Figure 4: Lead–lag from local maxima of |dEDA/dt| and |dRSA/dt| (same pairing as Figure 2A).

Output: plots/4/{pid}.png, plots/4/Figure4_summary.png
"""

import argparse

from plots.figure4 import plot_figure4


def main():
    parser = argparse.ArgumentParser(description="Figure 4 lead-lag (max |slope| events)")
    parser.add_argument("--data-dir", default="data")
    parser.add_argument("--output-dir", default="plots/4")
    parser.add_argument("--max-duration", type=float, default=300,
                        help="Middle N seconds for trace+markers (0 = full trace). Lags use full recording.")
    parser.add_argument("--preprocessed", action="store_true")
    parser.add_argument("--top-peak-pct", type=float, default=30,
                        help="Top X%% of max-|slope| events by height (default 30)")
    args = parser.parse_args()
    data_dir = "data/preprocessed" if args.preprocessed else args.data_dir
    max_dur = None if args.max_duration <= 0 else args.max_duration
    plot_figure4(
        data_dir=data_dir,
        output_dir=args.output_dir,
        max_duration_sec=max_dur,
        top_peak_pct=args.top_peak_pct,
    )


if __name__ == "__main__":
    main()
