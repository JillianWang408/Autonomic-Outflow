#!/usr/bin/env python3
"""Figure 3B: Jaccard heatmap for slope-based overlap (top % of |dEDA/dt| & |dRSA/dt|)."""

import argparse

from plots.figure3b import THRESHOLDS, generate_figure3b


def main():
    parser = argparse.ArgumentParser(description="Figure 3B Jaccard heatmap (slope)")
    parser.add_argument("--data-dir", default="data", help="Data directory")
    parser.add_argument("--output-dir", default="plots/3B", help="Output directory")
    parser.add_argument("--preprocessed", action="store_true")
    args = parser.parse_args()
    data_dir = "data/preprocessed" if args.preprocessed else args.data_dir
    generate_figure3b(
        data_dir=data_dir,
        output_dir=args.output_dir,
        thresholds=THRESHOLDS,
        use_preprocessed=args.preprocessed,
    )


if __name__ == "__main__":
    main()
