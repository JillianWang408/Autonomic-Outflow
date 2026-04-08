#!/usr/bin/env python3
"""
Generate Figure 2A: Lead-lag distribution (peak-to-peak).

Two histograms: EDA leads vs RSA leads.
Visual check: EDA/RSA traces with detected peaks marked.
Output: plots/2A/{pid}.png, plots/2A/Figure2A_summary.png

Usage:
    python plot_figure2A.py
    python plot_figure2A.py --data-dir data/preprocessed --output-dir plots/2A
"""

import argparse

from plots.figure2 import plot_figure2


def main():
    parser = argparse.ArgumentParser(description="Generate Figure 2A lead-lag histograms")
    parser.add_argument("--data-dir", default="data", help="Data directory (default: raw)")
    parser.add_argument("--output-dir", default="plots/2A", help="Output directory")
    parser.add_argument("--max-duration", type=float, default=300,
                        help="Middle N seconds for trace+peak markers only (0 = full trace). Lags always use full recording.")
    parser.add_argument("--preprocessed", action="store_true", help="Use preprocessed data (default: raw)")
    parser.add_argument("--top-peak-pct", type=float, default=30, help="Use only top X%% of peaks by amplitude (default: 30)")
    args = parser.parse_args()
    data_dir = "data/preprocessed" if args.preprocessed else args.data_dir
    max_dur = None if args.max_duration <= 0 else args.max_duration

    plot_figure2(
        data_dir=data_dir,
        output_dir=args.output_dir,
        max_duration_sec=max_dur,
        top_peak_pct=args.top_peak_pct,
    )


if __name__ == "__main__":
    main()
