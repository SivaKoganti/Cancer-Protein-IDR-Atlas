#!/usr/bin/env python3
"""Summarize top-ranked residue case studies for publication-style narrative support."""

import argparse
from pathlib import Path

import pandas as pd


def summarize_case_studies(atlas_dir, output_path):
    atlas_dir = Path(atlas_dir)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for path in sorted(atlas_dir.glob("*_atlas.tsv")):
        df = pd.read_csv(path, sep="\t")
        if {"gene", "pos", "llps_proxy", "vipp_score", "virus_interaction", "virus_count"}.issubset(df.columns):
            subset = df[["gene", "pos", "llps_proxy", "vipp_score", "virus_interaction", "virus_count"]].copy()
            subset["gene"] = subset["gene"].fillna(path.stem.replace("_atlas", ""))
            subset = subset.sort_values(["vipp_score", "llps_proxy"], ascending=False)
            top = subset.head(3).copy()
            top["case_rank"] = range(1, len(top) + 1)
            rows.append(top)

    if not rows:
        raise ValueError("No atlas files with required columns found")

    summary = pd.concat(rows, ignore_index=True)
    summary = summary.sort_values(["vipp_score", "llps_proxy"], ascending=False).head(20)
    summary.to_csv(output_path, sep="\t", index=False)
    return summary


def main():
    parser = argparse.ArgumentParser(description="Summarize top-ranked residue case studies")
    parser.add_argument("--atlas-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    summarize_case_studies(args.atlas_dir, args.output)


if __name__ == "__main__":
    main()
