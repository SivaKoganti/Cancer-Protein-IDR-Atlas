#!/usr/bin/env python3
"""Create a more informative residue-level ClinVar benchmark table from existing mappings."""

import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Create a benchmark-friendly ClinVar residue table")
    parser.add_argument("--atlas-dir", required=True)
    parser.add_argument("--variants", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    atlas_dir = Path(args.atlas_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    variants = pd.read_csv(args.variants, sep="\t", header=None, names=["gene", "pos", "ref", "alt", "clin"])
    for atlas_path in sorted(atlas_dir.glob("*_atlas.tsv")):
        atlas = pd.read_csv(atlas_path, sep="\t")
        gene = str(atlas.iloc[0]["gene"]).strip() if not atlas.empty else atlas_path.stem.replace("_atlas", "")
        gene_variants = variants[variants["gene"].astype(str).str.strip() == gene]
        if gene_variants.empty:
            continue
        rows = []
        for _, row in atlas.iterrows():
            pos = int(row["pos"])
            matches = gene_variants[gene_variants["pos"].astype(int) == pos]
            if matches.empty:
                rows.append({**row.to_dict(), "clinvar_count": 0, "clinvar_significance": ""})
            else:
                significance = ";".join(str(x) for x in matches["clin"].tolist())
                rows.append({**row.to_dict(), "clinvar_count": len(matches), "clinvar_significance": significance})
        out_path = output_dir / f"{gene}_atlas_with_clinvar.tsv"
        pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)

    print(f"Wrote ClinVar-augmented atlas files to {output_dir}")


if __name__ == "__main__":
    main()
