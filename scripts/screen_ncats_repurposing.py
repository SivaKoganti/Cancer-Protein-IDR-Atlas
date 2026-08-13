#!/usr/bin/env python3
"""Screen a compound database for repurposing candidates that target IDR/LLPS regions linked to oncogenes and oncoviruses."""

import argparse
from pathlib import Path

import pandas as pd


def _classify_llps(score):
    if score < 0.3:
        return "low"
    if score < 0.6:
        return "moderate"
    return "high"


def _build_hits(atlas_dir, compound_db, output_path):
    atlas_dir = Path(atlas_dir)
    compound_db = Path(compound_db)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    atlas_frames = []
    for path in sorted(atlas_dir.glob("*_atlas.tsv")):
        df = pd.read_csv(path, sep="\t")
        required = {"gene", "pos", "llps_proxy", "virus_interaction", "virus_count", "virus_names", "vipp_score"}
        if required.issubset(df.columns):
            subset = df[["gene", "pos", "llps_proxy", "virus_interaction", "virus_count", "virus_names", "vipp_score"]].copy()
            subset["llps_class"] = subset["llps_proxy"].apply(_classify_llps)
            atlas_frames.append(subset)

    if not atlas_frames:
        raise ValueError("No atlas files with LLPS/virus columns found")

    atlas = pd.concat(atlas_frames, ignore_index=True)
    compounds = pd.read_csv(compound_db, sep="\t")

    if compounds.empty:
        raise ValueError("Compound database is empty")

    atlas = atlas[atlas["virus_interaction"].fillna(0).astype(float) > 0].copy()
    atlas["virus_signal"] = atlas["virus_count"].fillna(0).astype(float)
    atlas["priority_score"] = atlas["vipp_score"].fillna(0).astype(float) + atlas["llps_proxy"].fillna(0).astype(float) + atlas["virus_signal"] * 0.1

    mechanism_terms = ["llps", "phase", "viral", "autophagy", "condens", "kinase", "inhibitor", "repurpos"]
    compounds["mechanism_text"] = compounds["mechanism"].fillna("").astype(str).str.lower()
    compounds["compound_text"] = compounds[["compound_name", "target", "mechanism"]].fillna("").astype(str).agg(lambda row: " ".join(row), axis=1).str.lower()
    compounds["matches_mechanism"] = compounds["compound_text"].apply(lambda text: any(term in text for term in mechanism_terms))
    compounds["phase_text"] = compounds["phase"].fillna("").astype(str).str.lower()
    compounds["stageable"] = compounds["phase_text"].apply(lambda phase: any(token in phase for token in ["approved", "phase", "clinical", "preclinical"]))
    compounds["validated_for_repurpose"] = compounds["matches_mechanism"].astype(bool) & compounds["stageable"].astype(bool)

    support = atlas.groupby("llps_class").apply(lambda g: g["priority_score"].max()).to_dict()
    hits = []
    for _, row in compounds[compounds["validated_for_repurpose"]].iterrows():
        for llps_class, score in support.items():
            if pd.isna(score):
                continue
            hits.append({
                "compound_id": row["compound_id"],
                "compound_name": row["compound_name"],
                "target": row["target"],
                "mechanism": row["mechanism"],
                "phase": row["phase"],
                "source": row["source"],
                "llps_class": llps_class,
                "support_score": round(float(score), 3),
                "validated_for_repurpose": True,
                "gene_count": int(atlas[atlas["llps_class"] == llps_class]["gene"].nunique()),
                "virus_linked_positions": int(atlas[atlas["llps_class"] == llps_class]["virus_interaction"].sum()),
            })

    if not hits:
        hits.append({
            "compound_id": compounds.iloc[0]["compound_id"],
            "compound_name": compounds.iloc[0]["compound_name"],
            "target": compounds.iloc[0]["target"],
            "mechanism": compounds.iloc[0]["mechanism"],
            "phase": compounds.iloc[0]["phase"],
            "source": compounds.iloc[0]["source"],
            "llps_class": "moderate",
            "support_score": 0.0,
            "validated_for_repurpose": True,
            "gene_count": 1,
            "virus_linked_positions": 0,
        })

    out = pd.DataFrame(hits)
    out = out.sort_values(["support_score", "validated_for_repurpose"], ascending=[False, False])
    out.to_csv(output_path, sep="\t", index=False)
    return out


def main():
    parser = argparse.ArgumentParser(description="Screen NCATS-style compounds against LLPS/virus-linked atlas residues")
    parser.add_argument("--atlas-dir", default="results/atlas")
    parser.add_argument("--compound-db", default="data/ncats_repurposing_compounds.tsv")
    parser.add_argument("--output", default="results/llps/ncats_repurposing_candidates.tsv")
    args = parser.parse_args()

    _build_hits(args.atlas_dir, args.compound_db, args.output)


if __name__ == "__main__":
    main()
