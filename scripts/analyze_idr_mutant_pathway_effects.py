#!/usr/bin/env python3
"""Generate IDR-focused LLPS mutant effects and pathway-level differential analysis."""

import argparse
from pathlib import Path
import re

import numpy as np
import pandas as pd
import yaml
from Bio import SeqIO

ALL_AA = list("ACDEFGHIKLMNPQRSTVWY")

DEFAULT_PATHWAYS = {
    "TP53_CELL_CYCLE": ["TP53", "ATM", "CHEK2", "CDKN2A", "RB1", "STK11"],
    "RTK_RAS_MAPK": ["EGFR", "ERBB2", "KRAS", "NRAS", "BRAF", "MAP2K1", "MAP2K2", "MET", "KIT", "PDGFRA", "ALK", "ROS1", "RET", "NTRK1"],
    "PI3K_AKT_MTOR": ["PIK3CA", "PTEN", "TSC1", "TSC2", "STK11"],
    "WNT_BETA_CATENIN": ["APC", "CTNNB1", "RNF43"],
    "NOTCH_HEDGEHOG": ["NOTCH1", "NOTCH2", "PTCH1"],
    "DNA_DAMAGE_REPAIR": ["BRCA1", "BRCA2", "ATM", "CHEK2", "MLH1", "MSH2", "MSH6", "PMS2", "ERCC2"],
    "CHROMATIN_REMODELING": ["ARID1A", "ARID2", "SMARCA4", "KMT2A", "KMT2D", "FBXW7"],
    "JAK_STAT_CYTOKINE": ["JAK2"],
    "TGF_BETA_SMAD": ["SMAD4"],
    "METABOLIC_HYPOXIA_REDOX": ["IDH1", "IDH2", "VHL", "NFE2L2", "GNAS"],
}


def normalize_gene_id(gene):
    if gene is None:
        return ""
    value = str(gene).strip()
    if not value:
        return ""
    return re.sub(r"_\\d+$", "", value)


def build_gene_to_pathways(pathway_cfg):
    gene_to_pathways = {}
    for pathway, genes in pathway_cfg.items():
        if not isinstance(genes, list):
            continue
        for gene in genes:
            norm_gene = normalize_gene_id(gene)
            gene_to_pathways.setdefault(norm_gene, []).append(pathway)
    return gene_to_pathways


def load_config(path):
    with open(path) as fh:
        cfg = yaml.safe_load(fh) or {}

    llps_cfg = cfg.get("llps", {})
    weights = {
        "window": int(llps_cfg.get("window", 41)),
        "qn": float(llps_cfg.get("qn_weight", 0.45)),
        "arom": float(llps_cfg.get("aromatic_weight", 0.25)),
        "hydrophobic": float(llps_cfg.get("hydrophobic_weight", 0.15)),
        "charge": float(llps_cfg.get("charge_weight", llps_cfg.get("charge_pattern_weight", 0.15))),
    }

    pathways = cfg.get("oncogenic_pathways", DEFAULT_PATHWAYS)
    if not isinstance(pathways, dict) or not pathways:
        pathways = DEFAULT_PATHWAYS

    return weights, pathways


def load_sequences(fasta_path):
    sequences = {}
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        gid = normalize_gene_id(rec.id.split("|")[0].strip())
        sequences[gid] = str(rec.seq).upper()
    return sequences


def llps_window_score(seq, center, window, weights):
    half = window // 2
    start = max(0, center - half)
    end = min(len(seq), center + half + 1)
    sub = seq[start:end]
    if not sub:
        return 0.0

    n = len(sub)
    qn = sum(1 for c in sub if c in "QN") / n
    arom = sum(1 for c in sub if c in "FYWI") / n
    charged = sum(1 for c in sub if c in "DEKRH") / n
    hydrophobic = sum(1 for c in sub if c in "AILMVFWY") / n

    return (
        weights["qn"] * qn
        + weights["arom"] * arom
        + weights["hydrophobic"] * hydrophobic
        - weights["charge"] * charged
    )


def compute_position_mutants(seq, pos_1_based, window, weights):
    idx = pos_1_based - 1
    wt_aa = seq[idx]
    wt_score = llps_window_score(seq, idx, window, weights)
    deltas = []

    for alt_aa in ALL_AA:
        if alt_aa == wt_aa:
            continue
        mut_seq = seq[:idx] + alt_aa + seq[idx + 1:]
        mut_score = llps_window_score(mut_seq, idx, window, weights)
        delta = mut_score - wt_score
        deltas.append((alt_aa, mut_score, delta))

    max_gain = max(deltas, key=lambda x: x[2])
    max_loss = min(deltas, key=lambda x: x[2])
    vulnerability = float(np.mean([abs(d[2]) for d in deltas]))

    return wt_aa, wt_score, vulnerability, max_gain, max_loss


def load_gene_idr_scores(atlas_path):
    df = pd.read_csv(atlas_path, sep="\t")
    required = {"gene", "pos", "aa", "iupred_score"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {atlas_path}: {sorted(missing)}")

    df = df.copy()
    df["gene"] = df["gene"].map(normalize_gene_id)
    for col, default in (
        ("ptm_count", 0),
        ("ptm_categories", ""),
        ("ptm_ids", ""),
        ("ptm_in_idr_count", 0),
        ("ptm_clinvar_overlap_count", 0),
    ):
        if col not in df.columns:
            df[col] = default
    return df


def build_mutant_table(atlas_dir, sequences, gene_to_pathways, window, weights, idr_threshold, impact_threshold):
    rows = []

    for atlas_file in sorted(Path(atlas_dir).glob("*_atlas.tsv")):
        if atlas_file.name == "global_disorder_phylogeny_atlas.tsv":
            continue
        atlas_df = load_gene_idr_scores(atlas_file)
        if atlas_df.empty:
            continue

        gene = normalize_gene_id(atlas_df.iloc[0]["gene"])
        seq = sequences.get(gene)
        if not seq:
            continue

        pathways = gene_to_pathways.get(gene, ["UNASSIGNED"])

        idr_df = atlas_df[pd.to_numeric(atlas_df["iupred_score"], errors="coerce") >= idr_threshold].copy()
        if idr_df.empty:
            continue

        for _, row in idr_df.iterrows():
            pos = int(row["pos"])
            if pos < 1 or pos > len(seq):
                continue

            wt_aa, wt_score, vulnerability, max_gain, max_loss = compute_position_mutants(
                seq,
                pos,
                window,
                weights,
            )

            for direction, payload in (("increase", max_gain), ("decrease", max_loss)):
                alt_aa, mut_score, delta = payload
                rows.append({
                    "gene": gene,
                    "pathways": ";".join(pathways),
                    "pathway_count": len(pathways),
                    "pos": pos,
                    "wt_aa": wt_aa,
                    "mut_aa": alt_aa,
                    "mutation": f"{wt_aa}{pos}{alt_aa}",
                    "idr_score": float(row["iupred_score"]),
                    "wt_llps_proxy": round(float(wt_score), 4),
                    "mut_llps_proxy": round(float(mut_score), 4),
                    "delta_llps": round(float(delta), 4),
                    "llps_vulnerability": round(float(vulnerability), 4),
                    "ptm_count": int(pd.to_numeric(row.get("ptm_count", 0), errors="coerce") if pd.notna(row.get("ptm_count", 0)) else 0),
                    "ptm_categories": str(row.get("ptm_categories", "")),
                    "ptm_ids": str(row.get("ptm_ids", "")),
                    "ptm_in_idr_count": int(pd.to_numeric(row.get("ptm_in_idr_count", 0), errors="coerce") if pd.notna(row.get("ptm_in_idr_count", 0)) else 0),
                    "ptm_clinvar_overlap_count": int(pd.to_numeric(row.get("ptm_clinvar_overlap_count", 0), errors="coerce") if pd.notna(row.get("ptm_clinvar_overlap_count", 0)) else 0),
                    "ptm_context": int(pd.to_numeric(row.get("ptm_count", 0), errors="coerce") if pd.notna(row.get("ptm_count", 0)) else 0) > 0,
                    "direction": direction,
                    "high_impact": abs(float(delta)) >= impact_threshold,
                })

    return pd.DataFrame(rows)


def summarize_pathways(mutants_df):
    exploded = mutants_df.assign(pathway=mutants_df["pathways"].str.split(";")).explode("pathway")

    summary = (
        exploded.groupby(["pathway", "direction"], as_index=False)
        .agg(
            n_mutants=("mutation", "size"),
            n_genes=("gene", "nunique"),
            mean_delta_llps=("delta_llps", "mean"),
            median_delta_llps=("delta_llps", "median"),
            std_delta_llps=("delta_llps", "std"),
            mean_abs_delta_llps=("delta_llps", lambda s: float(np.mean(np.abs(s)))),
            q90_abs_delta_llps=("delta_llps", lambda s: float(np.quantile(np.abs(s), 0.9))),
            high_impact_fraction=("high_impact", "mean"),
            mean_idr_score=("idr_score", "mean"),
            mean_ptm_count=("ptm_count", "mean"),
            ptm_context_fraction=("ptm_context", "mean"),
        )
        .fillna(0.0)
    )

    global_stats = (
        exploded.groupby("direction", as_index=False)
        .agg(
            global_mean_delta=("delta_llps", "mean"),
            global_mean_abs_delta=("delta_llps", lambda s: float(np.mean(np.abs(s)))),
        )
    )

    summary = summary.merge(global_stats, on="direction", how="left")
    summary["differential_mean_delta"] = summary["mean_delta_llps"] - summary["global_mean_delta"]
    summary["differential_mean_abs_delta"] = summary["mean_abs_delta_llps"] - summary["global_mean_abs_delta"]
    summary["vulnerability_index"] = summary["mean_abs_delta_llps"] / summary["global_mean_abs_delta"].replace(0, np.nan)
    summary["vulnerability_index"] = summary["vulnerability_index"].replace([np.inf, -np.inf], np.nan).fillna(0.0)

    summary = summary.sort_values(["direction", "mean_abs_delta_llps"], ascending=[True, False]).reset_index(drop=True)
    return summary


def summarize_pathway_ptm_context(mutants_df):
    exploded = mutants_df.assign(pathway=mutants_df["pathways"].str.split(";")).explode("pathway")
    exploded["ptm_context"] = exploded["ptm_context"].astype(bool)

    summary = (
        exploded.groupby(["pathway", "direction", "ptm_context"], as_index=False)
        .agg(
            n_mutants=("mutation", "size"),
            n_genes=("gene", "nunique"),
            mean_delta_llps=("delta_llps", "mean"),
            mean_abs_delta_llps=("delta_llps", lambda s: float(np.mean(np.abs(s)))),
            high_impact_fraction=("high_impact", "mean"),
            mean_ptm_count=("ptm_count", "mean"),
        )
        .fillna(0.0)
    )

    summary["ptm_context_label"] = np.where(summary["ptm_context"], "with_ptm", "without_ptm")

    paired = summary.pivot_table(
        index=["pathway", "direction"],
        columns="ptm_context_label",
        values="mean_abs_delta_llps",
        aggfunc="first",
    ).reset_index()

    if "with_ptm" not in paired.columns:
        paired["with_ptm"] = np.nan
    if "without_ptm" not in paired.columns:
        paired["without_ptm"] = np.nan

    paired["ptm_differential_mean_abs_delta"] = paired["with_ptm"].fillna(0.0) - paired["without_ptm"].fillna(0.0)
    paired = paired.rename(columns={
        "with_ptm": "mean_abs_delta_with_ptm",
        "without_ptm": "mean_abs_delta_without_ptm",
    })

    return summary.merge(paired, on=["pathway", "direction"], how="left")


def write_critical_analysis(mutants_df, pathway_summary_df, pathway_ptm_summary_df, output_path, idr_threshold):
    output_path.parent.mkdir(parents=True, exist_ok=True)

    overall = (
        pathway_summary_df.groupby("pathway", as_index=False)
        .agg(
            mean_abs_delta_llps=("mean_abs_delta_llps", "mean"),
            differential_mean_abs_delta=("differential_mean_abs_delta", "mean"),
            high_impact_fraction=("high_impact_fraction", "mean"),
        )
        .sort_values("differential_mean_abs_delta", ascending=False)
    )

    increases = pathway_summary_df[pathway_summary_df["direction"] == "increase"].sort_values("mean_delta_llps", ascending=False)
    decreases = pathway_summary_df[pathway_summary_df["direction"] == "decrease"].sort_values("mean_delta_llps", ascending=True)

    lines = []
    lines.append("# IDR Mutant LLPS Differential Analysis")
    lines.append("")
    lines.append("## Scope")
    lines.append(f"- IDR residues were defined as iupred_score >= {idr_threshold:.2f}.")
    lines.append("- For each IDR residue, two in-silico mutants were generated: maximum LLPS-increase and maximum LLPS-decrease substitution under the project LLPS proxy.")
    lines.append("- Pathway-level results aggregate these residue-level mutants; genes can contribute to multiple pathways.")
    lines.append("")

    lines.append("## Key Signals")
    top_overall = overall.head(3)
    if top_overall.empty:
        lines.append("- No pathway-level signals were detected because no IDR mutants were generated.")
    else:
        for _, row in top_overall.iterrows():
            lines.append(
                f"- {row['pathway']}: differential_mean_abs_delta={row['differential_mean_abs_delta']:.4f}, high_impact_fraction={row['high_impact_fraction']:.3f}"
            )

    lines.append("")
    lines.append("## Directional Interpretation")
    if not increases.empty:
        top_inc = increases.iloc[0]
        lines.append(
            f"- Strongest LLPS gain bias: {top_inc['pathway']} (mean_delta_llps={top_inc['mean_delta_llps']:.4f}, vulnerability_index={top_inc['vulnerability_index']:.3f})."
        )
    if not decreases.empty:
        top_dec = decreases.iloc[0]
        lines.append(
            f"- Strongest LLPS loss bias: {top_dec['pathway']} (mean_delta_llps={top_dec['mean_delta_llps']:.4f}, vulnerability_index={top_dec['vulnerability_index']:.3f})."
        )

    mean_abs = mutants_df["delta_llps"].abs().mean() if not mutants_df.empty else 0.0
    high_frac = mutants_df["high_impact"].mean() if not mutants_df.empty else 0.0
    ptm_frac = mutants_df["ptm_context"].mean() if not mutants_df.empty else 0.0
    lines.append(f"- Across all generated mutants, mean_abs_delta_llps={mean_abs:.4f} and high_impact_fraction={high_frac:.3f}.")
    lines.append(f"- PTM-context coverage among generated mutants: {ptm_frac:.3f}.")

    lines.append("")
    lines.append("## PTM Differential Signals")
    if pathway_ptm_summary_df is None or pathway_ptm_summary_df.empty:
        lines.append("- PTM-stratified pathway differentials were not available for this run.")
    else:
        dedup = pathway_ptm_summary_df.drop_duplicates(subset=["pathway", "direction"]).copy()
        top_ptm_gain = dedup.sort_values("ptm_differential_mean_abs_delta", ascending=False).head(3)
        top_ptm_loss = dedup.sort_values("ptm_differential_mean_abs_delta", ascending=True).head(3)

        lines.append("- Strongest PTM-associated increase-context signals (with PTM minus without PTM):")
        for _, row in top_ptm_gain.iterrows():
            lines.append(
                f"  - {row['pathway']} ({row['direction']}): Δmean_abs={row['ptm_differential_mean_abs_delta']:.4f} "
                f"[with_ptm={row['mean_abs_delta_with_ptm']:.4f}, without_ptm={row['mean_abs_delta_without_ptm']:.4f}]"
            )

        lines.append("- Strongest PTM-associated decrease-context signals (with PTM minus without PTM):")
        for _, row in top_ptm_loss.iterrows():
            lines.append(
                f"  - {row['pathway']} ({row['direction']}): Δmean_abs={row['ptm_differential_mean_abs_delta']:.4f} "
                f"[with_ptm={row['mean_abs_delta_with_ptm']:.4f}, without_ptm={row['mean_abs_delta_without_ptm']:.4f}]"
            )

    lines.append("")
    lines.append("## Critical Caveats")
    lines.append("- LLPS scores are proxy values derived from sequence composition windows; they are not direct biophysical phase-separation measurements.")
    lines.append("- Mutants were evaluated one residue at a time, so combinatorial epistasis and context-dependent compensation are not captured.")
    lines.append("- Pathway assignments are curated and intentionally overlapping; pathway-level contrasts should be interpreted as enrichment-like signals, not strict causal attribution.")
    lines.append("- The analysis does not model expression, post-translational regulation, tissue context, or treatment exposure, which can alter observed oncogenic behavior.")

    lines.append("")
    lines.append("## Suggested Follow-up")
    lines.append("- Prioritize high-impact mutants in pathways with elevated vulnerability_index for orthogonal structural and conservation cross-checks.")
    lines.append("- Validate top gain and loss mutants experimentally in pathway-relevant cellular models.")

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Analyze IDR mutant effects on LLPS and oncogenic pathways")
    parser.add_argument("--fasta", required=True, help="Input FASTA used for atlas generation")
    parser.add_argument("--atlas-dir", required=True, help="Directory containing *_atlas.tsv files")
    parser.add_argument("--config", default="config.yaml", help="Project config with llps and pathway settings")
    parser.add_argument("--output-dir", default="results/mutants", help="Output directory")
    parser.add_argument("--idr-threshold", type=float, default=0.5, help="IUPred threshold to define IDR residues")
    parser.add_argument("--impact-threshold", type=float, default=0.08, help="Absolute delta LLPS threshold for high-impact mutants")
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    weights, pathways = load_config(args.config)
    gene_to_pathways = build_gene_to_pathways(pathways)
    sequences = load_sequences(args.fasta)

    mutants = build_mutant_table(
        args.atlas_dir,
        sequences,
        gene_to_pathways,
        window=weights["window"],
        weights=weights,
        idr_threshold=args.idr_threshold,
        impact_threshold=args.impact_threshold,
    )

    mutants_out = outdir / "idr_llps_mutants.tsv"
    pathway_out = outdir / "pathway_differential_llps.tsv"
    pathway_ptm_out = outdir / "pathway_ptm_differential_llps.tsv"
    report_out = outdir / "pathway_critical_analysis.md"

    if mutants.empty:
        pd.DataFrame(
            columns=[
                "gene", "pathways", "pathway_count", "pos", "wt_aa", "mut_aa", "mutation", "idr_score",
                "wt_llps_proxy", "mut_llps_proxy", "delta_llps", "llps_vulnerability", "ptm_count",
                "ptm_categories", "ptm_ids", "ptm_in_idr_count", "ptm_clinvar_overlap_count", "ptm_context",
                "direction", "high_impact",
            ]
        ).to_csv(mutants_out, sep="\t", index=False)
        pd.DataFrame(
            columns=[
                "pathway", "direction", "n_mutants", "n_genes", "mean_delta_llps", "median_delta_llps", "std_delta_llps",
                "mean_abs_delta_llps", "q90_abs_delta_llps", "high_impact_fraction", "mean_idr_score", "global_mean_delta",
                "global_mean_abs_delta", "differential_mean_delta", "differential_mean_abs_delta", "vulnerability_index",
                "mean_ptm_count", "ptm_context_fraction",
            ]
        ).to_csv(pathway_out, sep="\t", index=False)
        pd.DataFrame(
            columns=[
                "pathway", "direction", "ptm_context", "n_mutants", "n_genes", "mean_delta_llps",
                "mean_abs_delta_llps", "high_impact_fraction", "mean_ptm_count", "ptm_context_label",
                "mean_abs_delta_with_ptm", "mean_abs_delta_without_ptm", "ptm_differential_mean_abs_delta",
            ]
        ).to_csv(pathway_ptm_out, sep="\t", index=False)
        write_critical_analysis(mutants, pd.DataFrame(), pd.DataFrame(), report_out, args.idr_threshold)
        print("No IDR mutants generated; wrote empty outputs.")
        return

    mutants = mutants.sort_values(["gene", "pos", "direction"]).reset_index(drop=True)
    pathway_summary = summarize_pathways(mutants)
    pathway_ptm_summary = summarize_pathway_ptm_context(mutants)

    mutants.to_csv(mutants_out, sep="\t", index=False)
    pathway_summary.to_csv(pathway_out, sep="\t", index=False)
    pathway_ptm_summary.to_csv(pathway_ptm_out, sep="\t", index=False)
    write_critical_analysis(mutants, pathway_summary, pathway_ptm_summary, report_out, args.idr_threshold)

    print(f"Wrote mutant table: {mutants_out} ({len(mutants)} rows)")
    print(f"Wrote pathway differential summary: {pathway_out} ({len(pathway_summary)} rows)")
    print(f"Wrote pathway PTM differential summary: {pathway_ptm_out} ({len(pathway_ptm_summary)} rows)")
    print(f"Wrote critical analysis: {report_out}")


if __name__ == "__main__":
    main()
