#!/usr/bin/env python3
"""
Generate Figure 1B: P(EDA|RSA) and P(RSA|EDA) heatmaps.

For each participant (y-axis) × threshold (x-axis), with group average in bottom row.
Saves to plots/1B/:
  - P_EDA_given_RSA.png
  - P_RSA_given_EDA.png

Usage:
    python plot_figure1B.py
    python plot_figure1B.py --data-dir data/preprocessed --output-dir plots/1B
"""

import argparse

from plots.figure1b import generate_figure1b, THRESHOLDS


def main():
    parser = argparse.ArgumentParser(description="Generate Figure 1B heatmaps")
    parser.add_argument("--data-dir", default="data/preprocessed", help="Data directory")
    parser.add_argument("--output-dir", default="plots/1B", help="Output directory")
    parser.add_argument("--raw", action="store_true", help="Use raw data instead of preprocessed")
    args = parser.parse_args()

    generate_figure1b(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        thresholds=THRESHOLDS,
        use_preprocessed=not args.raw,
    )


if __name__ == "__main__":
    main()
