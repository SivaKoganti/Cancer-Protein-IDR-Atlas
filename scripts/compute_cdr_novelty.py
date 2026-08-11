#!/usr/bin/env python3
"""
Identify two novel functional region classes per residue:

  CDR (Conserved Disordered Region)
      conservation >= CDR_CONS_THRESH AND iupred_score >= CDR_IDR_THRESH
      These are regulatory hubs under evolutionary constraint despite lacking fixed structure.

  Discordant Region
      plddt >= DISCORDANT_PLDDT_THRESH AND iupred_score >= DISCORDANT_IDR_THRESH
      AlphaFold predicts structure but IUPred predicts disorder → coupled-folding-and-binding sites.
      Cancer mutations here likely disrupt conditional interactions.

Also emits a summary TSV of CDR and discordant region counts per gene.
"""

import argparse
from pathlib import Path

import pandas as pd

# Classification thresholds
CDR_CONS_THRESH = 0.70
CDR_IDR_THRESH = 0.60
DISCORDANT_PLDDT_THRESH = 70.0
DISCORDANT_IDR_THRESH = 0.50

# Minimum run length to call a region (filters isolated noise residues)
MIN_RUN = 3


def label_runs(flags: list[int], min_run: int = MIN_RUN) -> list[int]:
    """Zero out runs shorter than min_run to reduce noise."""
    result = list(flags)
    i = 0
    while i < len(result):
        if result[i]:
            j = i
            while j < len(result) and result[j]:
                j += 1
            if (j - i) < min_run:
                for k in range(i, j):
                    result[k] = 0
            i = j
        else:
            i += 1
    return result


def process_gene(gene_df: pd.DataFrame, plddt_df: pd.DataFrame | None) -> pd.DataFrame:
    df = gene_df.copy()

    # CDR flag
    cdr_raw = (
        (df["conservation"] >= CDR_CONS_THRESH) &
        (df["iupred_score"] >= CDR_IDR_THRESH)
    ).astype(int).tolist()
    df["cdr_flag"] = label_runs(cdr_raw)

    # Discordant region flag (requires pLDDT data)
    if plddt_df is not None and not plddt_df.empty:
        merged = df[["pos", "iupred_score"]].merge(
            plddt_df[["pos", "plddt"]], on="pos", how="left"
        )
        merged["plddt"] = merged["plddt"].fillna(50.0)   # neutral fallback
        disc_raw = (
            (merged["plddt"] >= DISCORDANT_PLDDT_THRESH) &
            (merged["iupred_score"] >= DISCORDANT_IDR_THRESH)
        ).astype(int).tolist()
        df["discordant_flag"] = label_runs(disc_raw)
        df["plddt"] = merged["plddt"].values
    else:
        df["discordant_flag"] = 0
        df["plddt"] = float("nan")

    return df


def main():
    parser = argparse.ArgumentParser(description="Annotate CDR and discordant regions in atlas.")
    parser.add_argument("--atlas-dir", required=True, help="Directory with *_atlas.tsv files")
    parser.add_argument("--plddt", default="results/alphafold/plddt_scores.tsv",
                        help="pLDDT scores TSV (from download_alphafold_plddt.py)")
    parser.add_argument("--output-dir", default="results/cdr")
    args = parser.parse_args()

    atlas_dir = Path(args.atlas_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    plddt_all = None
    plddt_path = Path(args.plddt)
    if plddt_path.exists():
        plddt_all = pd.read_csv(plddt_path, sep="\t")
        print(f"Loaded pLDDT for {plddt_all['gene'].nunique()} genes")
    else:
        print("pLDDT file not found — discordant flags will be 0")

    summary_rows = []
    all_frames = []

    for tsv in sorted(atlas_dir.glob("*_atlas.tsv")):
        gene = tsv.stem.replace("_atlas", "")
        df = pd.read_csv(tsv, sep="\t")
        if "conservation" not in df.columns:
            df["conservation"] = 0.0
        if "iupred_score" not in df.columns:
            df["iupred_score"] = 0.0

        plddt_gene = None
        if plddt_all is not None:
            plddt_gene = plddt_all[plddt_all["gene"] == gene]

        annotated = process_gene(df, plddt_gene)
        annotated.to_csv(out_dir / f"{gene}_cdr.tsv", sep="\t", index=False)
        all_frames.append(annotated)

        summary_rows.append({
            "gene": gene,
            "length": len(df),
            "cdr_count": int(annotated["cdr_flag"].sum()),
            "cdr_fraction": round(annotated["cdr_flag"].mean(), 4),
            "discordant_count": int(annotated["discordant_flag"].sum()),
            "discordant_fraction": round(annotated["discordant_flag"].mean(), 4),
            "plddt_mean": round(annotated["plddt"].mean(), 2) if not annotated["plddt"].isna().all() else float("nan"),
        })
        print(f"  {gene}: {int(annotated['cdr_flag'].sum())} CDR / "
              f"{int(annotated['discordant_flag'].sum())} discordant residues")

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(out_dir / "cdr_summary.tsv", sep="\t", index=False)

    # Full combined table
    if all_frames:
        pd.concat(all_frames, ignore_index=True).to_csv(
            out_dir / "all_genes_cdr.tsv", sep="\t", index=False
        )

    print(f"\n✓ CDR/discordant annotations saved → {out_dir}")
    print(f"  Genes with CDRs: {(summary_df['cdr_count'] > 0).sum()} / {len(summary_df)}")
    print(f"  Genes with discordant regions: {(summary_df['discordant_count'] > 0).sum()} / {len(summary_df)}")


if __name__ == "__main__":
    main()
