#!/usr/bin/env python3
"""
Generate Figure 1B: Jaccard overlap heatmap.

Participants (y-axis) × thresholds (x-axis), bottom row = group average.
Output: plots/1B/Jaccard_overlap.png

Usage:
    python plot_figure1B.py
    python plot_figure1B.py --preprocessed  # use preprocessed data
"""

import argparse

from plots.figure1b import generate_figure1b, THRESHOLDS


def main():
    parser = argparse.ArgumentParser(description="Generate Figure 1B Jaccard overlap heatmap")
    parser.add_argument("--data-dir", default="data", help="Data directory (default: raw)")
    parser.add_argument("--output-dir", default="plots/1B", help="Output directory")
    parser.add_argument("--preprocessed", action="store_true", help="Use preprocessed data (default: raw)")
    args = parser.parse_args()
    data_dir = "data/preprocessed" if args.preprocessed else args.data_dir

    generate_figure1b(
        data_dir=data_dir,
        output_dir=args.output_dir,
        thresholds=THRESHOLDS,
        use_preprocessed=args.preprocessed,
    )


if __name__ == "__main__":
    main()
