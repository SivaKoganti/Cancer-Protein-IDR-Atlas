#!/usr/bin/env python3
"""Attach ClinVar variant annotations to a per-gene atlas TSV."""

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd


def load_atlas(path):
    if not path.exists():
        raise FileNotFoundError(f"atlas file not found: {path}")
    return pd.read_csv(path, sep="\t")


def load_variants(path):
    if not path.exists():
        raise FileNotFoundError(f"variants file not found: {path}")
    return pd.read_csv(path, sep="\t", header=None, names=["gene", "pos", "ref", "alt", "clin"])


def build_variant_lookup(variants_df):
    lookup = defaultdict(list)
    for _, row in variants_df.iterrows():
        gene = str(row["gene"]).strip()
        pos = int(row["pos"])
        lookup[gene].append({
            "pos": pos,
            "ref": str(row["ref"]).strip(),
            "alt": str(row["alt"]).strip(),
            "clin": str(row["clin"]).strip(),
        })
    return lookup


def main():
    parser = argparse.ArgumentParser(description="Merge ClinVar variants into an atlas TSV")
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--variants", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    atlas_path = Path(args.atlas)
    variants_path = Path(args.variants)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    atlas = load_atlas(atlas_path)
    variants = load_variants(variants_path)
    variant_lookup = build_variant_lookup(variants)

    gene = str(atlas.iloc[0]["gene"]).strip() if not atlas.empty else None
    if gene is None:
        raise ValueError("atlas file does not contain any rows")

    merged_rows = []
    for _, row in atlas.iterrows():
        pos = int(row["pos"])
        matching_variants = [
            v for v in variant_lookup.get(str(row["gene"]).strip(), []) if v["pos"] == pos
        ]
        clinvar_sig = ";".join(v["clin"] for v in matching_variants)
        merged_rows.append({
            **row.to_dict(),
            "clinvar_count": len(matching_variants),
            "clinvar_significance": clinvar_sig,
        })

    merged_df = pd.DataFrame(merged_rows)
    output_path = outdir / f"{gene}_atlas_with_clinvar.tsv"
    merged_df.to_csv(output_path, sep="\t", index=False)
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
