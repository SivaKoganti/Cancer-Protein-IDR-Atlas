#!/usr/bin/env python3
"""Generate a small-molecule library for LLPS- and virus-linked residues."""

import argparse
from pathlib import Path

import pandas as pd


def _classify_llps(score):
    if score < 0.3:
        return "low"
    if score < 0.6:
        return "moderate"
    return "high"


def generate_library(atlas_dir, output_path):
    atlas_dir = Path(atlas_dir)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    frames = []
    for path in sorted(atlas_dir.glob("*_atlas.tsv")):
        df = pd.read_csv(path, sep="\t")
        if {"llps_proxy", "virus_interaction", "gene"}.issubset(df.columns):
            subset = df[["gene", "pos", "llps_proxy", "virus_interaction", "virus_count"]].copy()
            subset["llps_class"] = subset["llps_proxy"].apply(_classify_llps)
            frames.append(subset)

    if not frames:
        raise ValueError("No atlas files with LLPS and virus columns found")

    combined = pd.concat(frames, ignore_index=True)
    summary = combined.groupby("llps_class").agg(
        residues=("llps_class", "size"),
        virus_linked_residues=("virus_interaction", "sum"),
        max_llps=("llps_proxy", "max"),
    ).reset_index()

    compound_map = {
        "low": ["small-molecule scaffold stabilizer", "conformation-aware chaperone mimetic"],
        "moderate": ["1,6-hexanediol analog", "aromatic patch disruptor"],
        "high": ["RNA-binding condensate inhibitor", "phase-separation modulator"],
    }

    summary["candidate_compounds"] = summary["llps_class"].map(lambda cls: "; ".join(compound_map.get(cls, [])))
    summary["rationale"] = summary["llps_class"].map({
        "low": "Prefer scaffolds that stabilize local interactions without broadly perturbing condensates.",
        "moderate": "Target weak hydrophobic and aromatic contacts that can tip borderline LLPS assemblies.",
        "high": "Prioritize compounds that disrupt multivalent RNA/protein or protein/protein condensates at high LLPS propensity.",
    })
    summary.to_csv(output_path, sep="\t", index=False)
    return summary


def main():
    parser = argparse.ArgumentParser(description="Generate small-molecule candidates for LLPS states")
    parser.add_argument("--atlas-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    generate_library(args.atlas_dir, args.output)


if __name__ == "__main__":
    main()
