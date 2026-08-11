#!/usr/bin/env python3
"""Generate presentation-ready VIPP hotspot tables from per-gene atlas TSV files."""

import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Summarize top VIPP residues across atlas files.")
    parser.add_argument("--atlas-dir", required=True, help="Directory containing *_atlas.tsv files")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("--per-gene-top", type=int, default=10, help="Top residues per gene")
    parser.add_argument("--global-top", type=int, default=200, help="Top residues overall")
    args = parser.parse_args()

    atlas_dir = Path(args.atlas_dir)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    frames = []
    for tsv in sorted(atlas_dir.glob("*_atlas.tsv")):
        gene = tsv.stem.replace("_atlas", "")
        df = pd.read_csv(tsv, sep="\t")
        required = {"gene", "pos", "aa", "vipp_score"}
        if not required.issubset(df.columns):
            continue

        subset_cols = [
            "gene", "pos", "aa", "vipp_score", "virus_interaction", "virus_count",
            "iupred_score", "llps_proxy", "conservation", "variants",
        ]
        keep = [c for c in subset_cols if c in df.columns]
        top_df = df.sort_values(["vipp_score", "virus_interaction", "llps_proxy"], ascending=False).head(args.per_gene_top)
        top_df = top_df[keep].copy()
        top_df.insert(0, "rank_type", f"top_{args.per_gene_top}_per_gene")
        top_df.insert(1, "gene_rank", range(1, len(top_df) + 1))
        frames.append(top_df)

    if not frames:
        pd.DataFrame(columns=["rank_type", "gene_rank", "gene", "pos", "aa", "vipp_score"]).to_csv(
            output_path, sep="\t", index=False
        )
        print(f"No VIPP-enabled atlas files found under {atlas_dir}")
        return

    per_gene_table = pd.concat(frames, ignore_index=True)

    global_table = per_gene_table.sort_values(
        ["vipp_score", "virus_interaction", "llps_proxy"], ascending=False
    ).head(args.global_top).copy()
    global_table.insert(0, "global_rank", range(1, len(global_table) + 1))

    # Keep a single presentation table with both per-gene and global rank tags.
    global_table.to_csv(output_path, sep="\t", index=False)
    print(f"Wrote VIPP hotspots: {output_path} ({len(global_table)} rows)")


if __name__ == "__main__":
    main()
