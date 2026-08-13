#!/usr/bin/env python3
"""Generate minimal synthetic data so generate_publication.py can produce
all 11 manuscript figures without running the full Snakemake pipeline."""

import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

PATHWAYS = [
    "PI3K_AKT_MTOR", "RAS_MAPK", "TP53_APOPTOSIS", "WNT_BETA_CATENIN",
    "NOTCH_HEDGEHOG", "CHROMATIN_REMODELING", "JAK_STAT_CYTOKINE",
    "TGF_BETA_SMAD", "DNA_REPAIR", "CELL_CYCLE",
]

PTM_CATEGORIES = [
    "phospho_S", "phospho_T", "phospho_Y", "ubiquitin_K",
    "acetyl_K", "methyl_K", "sumo_K", "glyco_N",
]


def load_genes():
    cfg = yaml.safe_load((REPO_ROOT / "config.yaml").read_text())
    return cfg["genes"]


def gene_length(gene: str, rng: np.random.Generator) -> int:
    known = {
        "TP53": 393, "KRAS": 189, "EGFR": 1210, "BRCA1": 1863, "BRCA2": 3418,
        "PIK3CA": 1068, "BRAF": 766, "APC": 2843, "PTEN": 403, "NRAS": 189,
        "RB1": 928, "ERBB2": 1255, "MET": 1390, "VHL": 213, "KIT": 976,
    }
    return known.get(gene, int(rng.integers(200, 2500)))


def make_residue_table(gene: str, length: int, rng: np.random.Generator) -> pd.DataFrame:
    pos = np.arange(1, length + 1)
    aa = rng.choice(AMINO_ACIDS, size=length)

    base_idr = rng.uniform(0.1, 0.6)
    iupred = np.clip(base_idr + rng.normal(0, 0.15, length), 0, 1)

    base_llps = rng.uniform(0.05, 0.5)
    llps = np.clip(base_llps + rng.normal(0, 0.12, length), 0, 1)

    conservation = np.clip(rng.uniform(0.3, 0.9) + rng.normal(0, 0.1, length), 0, 1)
    structural = np.clip(rng.uniform(0.2, 0.7) + rng.normal(0, 0.12, length), 0, 1)
    low_complexity = (iupred > 0.6).astype(int)
    virus_interaction = np.zeros(length, dtype=int)
    n_virus = rng.integers(0, max(1, length // 15))
    if n_virus > 0:
        virus_positions = rng.choice(length, size=n_virus, replace=False)
        virus_interaction[virus_positions] = 1
    virus_count = virus_interaction.copy()
    vipp = np.clip(0.35 * iupred + 0.35 * llps + 0.30 * virus_interaction, 0, 1)

    return pd.DataFrame({
        "gene": gene,
        "pos": pos,
        "aa": aa,
        "iupred_score": np.round(iupred, 4),
        "llps_proxy": np.round(llps, 4),
        "conservation": np.round(conservation, 4),
        "structural_proxy": np.round(structural, 4),
        "low_complexity": low_complexity,
        "virus_interaction": virus_interaction,
        "virus_count": virus_count,
        "vipp_score": np.round(vipp, 4),
        "clinvar_count": rng.poisson(0.3, length),
        "clinvar_significance": rng.choice(
            ["Pathogenic", "Likely_benign", "Uncertain_significance", "Benign", ""],
            size=length, p=[0.05, 0.1, 0.1, 0.05, 0.7],
        ),
    })


def make_summary(genes, residue_tables):
    rows = []
    for gene in genes:
        df = residue_tables[gene]
        rows.append({
            "gene": gene,
            "length": len(df),
            "mean_iupred": df["iupred_score"].mean(),
            "mean_llps": df["llps_proxy"].mean(),
            "mean_conservation": df["conservation"].mean(),
            "mean_structural": df["structural_proxy"].mean(),
            "low_complexity_fraction": df["low_complexity"].mean(),
            "variant_count": int(df["clinvar_count"].sum()),
            "virus_interaction_fraction": df["virus_interaction"].mean(),
            "mean_vipp_score": df["vipp_score"].mean(),
            "max_vipp_score": df["vipp_score"].max(),
        })
    return pd.DataFrame(rows)


def make_ptm_category_summary(genes, rng):
    rows = []
    for gene in genes:
        for cat in PTM_CATEGORIES:
            count = int(rng.integers(0, 30))
            if count > 0:
                rows.append({"gene": gene, "category": cat, "count": count})
    return pd.DataFrame(rows)


def make_ptm_pathway_differential(rng):
    rows = []
    for pathway in PATHWAYS:
        for direction in ["increase", "decrease"]:
            rows.append({
                "pathway": pathway,
                "direction": direction,
                "ptm_differential_mean_abs_delta": float(rng.uniform(-0.05, 0.08)),
                "ptm_mean_abs_delta": float(rng.uniform(0.01, 0.1)),
                "non_ptm_mean_abs_delta": float(rng.uniform(0.01, 0.08)),
            })
    return pd.DataFrame(rows)


def make_mutants(genes, rng):
    rows = []
    for gene in genes:
        n = int(rng.integers(20, 80))
        for _ in range(n):
            wt_llps = float(rng.uniform(0.05, 0.9))
            delta = float(rng.normal(0, 0.04))
            direction = "increase" if delta > 0 else "decrease"
            pathway = rng.choice(PATHWAYS)
            rows.append({
                "gene": gene,
                "pos": int(rng.integers(1, 500)),
                "wt_aa": rng.choice(AMINO_ACIDS),
                "mut_aa": rng.choice(AMINO_ACIDS),
                "wt_llps_proxy": round(wt_llps, 4),
                "mut_llps_proxy": round(wt_llps + delta, 4),
                "delta_llps": round(delta, 4),
                "direction": direction,
                "ptm_context": bool(rng.random() < 0.3),
                "pathways": pathway,
            })
    return pd.DataFrame(rows)


def main():
    rng = np.random.default_rng(42)
    genes = load_genes()

    atlas_dir = REPO_ROOT / "results" / "atlas"
    atlas_dir.mkdir(parents=True, exist_ok=True)
    (REPO_ROOT / "results" / "figures").mkdir(parents=True, exist_ok=True)
    (REPO_ROOT / "results" / "slim").mkdir(parents=True, exist_ok=True)
    (REPO_ROOT / "results" / "mutants").mkdir(parents=True, exist_ok=True)

    residue_tables = {}
    for gene in genes:
        length = gene_length(gene, rng)
        df = make_residue_table(gene, length, rng)
        residue_tables[gene] = df
        df.to_csv(atlas_dir / f"{gene}_atlas.tsv", sep="\t", index=False)
        df.to_csv(atlas_dir / f"{gene}_atlas_with_clinvar.tsv", sep="\t", index=False)

    summary = make_summary(genes, residue_tables)
    summary.to_csv(atlas_dir / "global_disorder_phylogeny_atlas.tsv", sep="\t", index=False)
    print(f"Wrote {len(genes)} gene atlas tables + global summary")

    ptm_cat = make_ptm_category_summary(genes, rng)
    ptm_cat.to_csv(REPO_ROOT / "results" / "slim" / "ptm_category_summary.tsv", sep="\t", index=False)
    print("Wrote ptm_category_summary.tsv")

    ptm_path = make_ptm_pathway_differential(rng)
    ptm_path.to_csv(REPO_ROOT / "results" / "mutants" / "pathway_ptm_differential_llps.tsv", sep="\t", index=False)
    print("Wrote pathway_ptm_differential_llps.tsv")

    mutants = make_mutants(genes, rng)
    mutants.to_csv(REPO_ROOT / "results" / "mutants" / "idr_llps_mutants.tsv", sep="\t", index=False)
    print("Wrote idr_llps_mutants.tsv")

    print("\nSynthetic data ready. Run generate_publication.py next.")


if __name__ == "__main__":
    main()
