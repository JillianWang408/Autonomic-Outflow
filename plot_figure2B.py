#!/usr/bin/env python3
"""
Figure 2B: Lead–lag restricted to co‑activation segments (Figure 1 top‑% active on both channels),
only segments that contain ≥1 EDA and ≥1 RSA peak; pairing within each segment.

Output: plots/2B/{pid}.png, plots/2B/Figure2B_summary.png
"""

import argparse

from plots.figure2b import plot_figure2b


def main():
    parser = argparse.ArgumentParser(description="Generate Figure 2B (co‑activation segment lags)")
    parser.add_argument("--data-dir", default="data")
    parser.add_argument("--output-dir", default="plots/2B")
    parser.add_argument("--max-duration", type=float, default=300,
                        help="Middle N seconds for trace only (0 = full trace). Lags use full recording.")
    parser.add_argument("--preprocessed", action="store_true")
    parser.add_argument("--top-peak-pct", type=float, default=30,
                        help="Peaks: keep top X%% by amplitude (same as 2A)")
    parser.add_argument("--co-threshold-pct", type=float, default=40,
                        help="Co‑activation: top X%% active per channel (Figure 1 style)")
    parser.add_argument("--co-extend-sec", type=float, default=3,
                        help="Extend each co‑active run by this many seconds before/after (clamped)")
    args = parser.parse_args()
    data_dir = "data/preprocessed" if args.preprocessed else args.data_dir
    max_dur = None if args.max_duration <= 0 else args.max_duration

    plot_figure2b(
        data_dir=data_dir,
        output_dir=args.output_dir,
        max_duration_sec=max_dur,
        top_peak_pct=args.top_peak_pct,
        co_threshold_pct=args.co_threshold_pct,
        co_extend_sec=args.co_extend_sec,
    )


if __name__ == "__main__":
    main()
