#!/usr/bin/env python3
"""
Generate Figure 2: Lead-lag distribution (peak-to-peak).

Two histograms: EDA leads vs RSA leads.
Visual check: EDA/RSA traces with detected peaks marked.
Output: plots/2/Figure2_lead_lag.png

Usage:
    python plot_figure2.py
    python plot_figure2.py --data-dir data/preprocessed --output-dir plots/2
"""

import argparse

from plots.figure2 import plot_figure2


def main():
    parser = argparse.ArgumentParser(description="Generate Figure 2 lead-lag histograms")
    parser.add_argument("--data-dir", default="data", help="Data directory (default: raw)")
    parser.add_argument("--output-dir", default="plots/2", help="Output directory")
    parser.add_argument("--max-duration", type=float, default=300, help="Middle N seconds per participant")
    parser.add_argument("--preprocessed", action="store_true", help="Use preprocessed data (default: raw)")
    parser.add_argument("--top-peak-pct", type=float, default=30, help="Use only top X%% of peaks by amplitude (default: 30)")
    args = parser.parse_args()
    data_dir = "data/preprocessed" if args.preprocessed else args.data_dir

    plot_figure2(
        data_dir=data_dir,
        output_dir=args.output_dir,
        max_duration_sec=args.max_duration,
        top_peak_pct=args.top_peak_pct,
    )


if __name__ == "__main__":
    main()
