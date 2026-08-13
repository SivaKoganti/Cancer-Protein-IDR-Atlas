#!/usr/bin/env python3
"""Validate that the computational repository contains the core atlas, benchmark, and repurposing outputs."""

import argparse
from pathlib import Path

import pandas as pd


def validate_repository(atlas_dir, results_dir, output_path):
    atlas_dir = Path(atlas_dir)
    results_dir = Path(results_dir)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    checks = []

    atlas_files = list(atlas_dir.glob("*_atlas.tsv"))
    checks.append(("atlas_files", len(atlas_files) > 0, f"Found {len(atlas_files)} atlas files"))

    llps_library = results_dir / "llps" / "llps_compound_library.tsv"
    checks.append(("llps_library", llps_library.exists(), f"LLPS library exists: {llps_library.exists()}"))

    benchmark_metrics = results_dir / "benchmark" / "model_metrics.tsv"
    checks.append(("benchmark_metrics", benchmark_metrics.exists(), f"Benchmark metrics exist: {benchmark_metrics.exists()}"))

    ablation_matrix = results_dir / "benchmark" / "ablation_matrix.tsv"
    checks.append(("ablation_matrix", ablation_matrix.exists(), f"Ablation matrix exists: {ablation_matrix.exists()}"))

    repurposing = results_dir / "llps" / "ncats_repurposing_candidates.tsv"
    checks.append(("repurposing", repurposing.exists(), f"Repurposing candidates exist: {repurposing.exists()}"))

    if llps_library.exists():
        llps_df = pd.read_csv(llps_library, sep="\t")
        checks.append(("llps_library_not_empty", not llps_df.empty, f"LLPS library rows: {len(llps_df)}"))

    if benchmark_metrics.exists():
        bench_df = pd.read_csv(benchmark_metrics, sep="\t")
        checks.append(("benchmark_rows", not bench_df.empty, f"Benchmark rows: {len(bench_df)}"))

    if repurposing.exists():
        rep_df = pd.read_csv(repurposing, sep="\t")
        checks.append(("repurposing_rows", not rep_df.empty, f"Repurposing rows: {len(rep_df)}"))

    report = pd.DataFrame(
        [
            {"check_name": name, "status": "PASS" if passed else "FAIL", "details": detail}
            for name, passed, detail in checks
        ]
    )
    report.to_csv(output_path, sep="\t", index=False)
    return report


def main():
    parser = argparse.ArgumentParser(description="Validate the computational repository outputs")
    parser.add_argument("--atlas-dir", default="results/atlas")
    parser.add_argument("--results-dir", default="results")
    parser.add_argument("--output", default="results/validation/computational_repository_validation.tsv")
    args = parser.parse_args()

    validate_repository(args.atlas_dir, args.results_dir, args.output)


if __name__ == "__main__":
    main()
