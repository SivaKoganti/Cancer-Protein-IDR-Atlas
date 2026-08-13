#!/usr/bin/env python3
"""generate_publication.py

Generate publication-quality figures and summary statistics for the
Cancer Protein IDR Atlas manuscript.
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    from scripts.render_chromosome_atlas import GENE_LOCI, CHROM_SIZES
except ModuleNotFoundError:
    from render_chromosome_atlas import GENE_LOCI, CHROM_SIZES


def load_data(atlas_dir: Path, summary_path: Path):
    """Load atlas data."""
    summary = pd.read_csv(summary_path, sep="\t")
    
    residue_data = {}
    for _, row in summary.iterrows():
        gene = str(row["gene"])
        p = atlas_dir / f"{gene}_atlas_with_clinvar.tsv"
        if p.exists():
            df = pd.read_csv(p, sep="\t")
            residue_data[gene] = df
    
    return summary, residue_data


def figure_1_overview(summary: pd.DataFrame, output_dir: Path):
    """Figure 1: Global atlas overview and statistics."""
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)
    
    # A: Genes per chromosome
    ax_a = fig.add_subplot(gs[0, :2])
    chrom_counts = {}
    chrom_loci = {
        "TP53": "17", "KRAS": "12", "EGFR": "7", "BRCA1": "17", "BRCA2": "13",
        "PIK3CA": "3", "BRAF": "7", "APC": "5", "PTEN": "10", "NRAS": "1",
        "CDKN2A": "9", "RB1": "13", "ERBB2": "17", "SMAD4": "18", "STK11": "19",
        "ATM": "11", "CHEK2": "22", "MET": "7", "VHL": "3", "GNAS": "20",
        "IDH1": "2", "IDH2": "15", "JAK2": "9", "KIT": "4", "PDGFRA": "4",
        "ALK": "2", "ROS1": "6", "RET": "10", "NFE2L2": "2", "PTCH1": "9",
        "FBXW7": "4", "CTNNB1": "3", "MAP2K1": "15", "MAP2K2": "19",
        "NOTCH1": "9", "NOTCH2": "1", "ARID1A": "1", "ARID2": "12",
        "SMARCA4": "19", "KMT2A": "11", "KMT2D": "12", "TSC1": "9",
        "TSC2": "16", "MLH1": "3", "MSH2": "2", "MSH6": "2", "PMS2": "7",
        "ERCC2": "19", "RNF43": "17", "NTRK1": "1",
    }
    for gene in summary["gene"]:
        c = chrom_loci.get(str(gene), "X")
        chrom_counts[c] = chrom_counts.get(c, 0) + 1
    
    chroms = sorted(chrom_counts.keys(), key=lambda x: (x!='X', x!='Y', int(x) if x.isdigit() else 0))
    counts = [chrom_counts.get(c, 0) for c in chroms]
    colors = plt.cm.tab20(np.linspace(0, 1, len(chroms)))
    ax_a.bar(chroms, counts, color=colors, edgecolor='black', linewidth=0.7)
    ax_a.set_xlabel("Chromosome", fontsize=11, fontweight='bold')
    ax_a.set_ylabel("Number of Genes", fontsize=11, fontweight='bold')
    ax_a.set_title("A. Gene Distribution Across Human Chromosomes", fontsize=12, fontweight='bold', loc='left')
    ax_a.grid(axis='y', alpha=0.3)
    
    # B: Protein length distribution
    ax_b = fig.add_subplot(gs[0, 2])
    ax_b.hist(summary["length"], bins=20, color='steelblue', edgecolor='black', linewidth=0.7, alpha=0.8)
    ax_b.set_xlabel("Length (aa)", fontsize=10, fontweight='bold')
    ax_b.set_ylabel("Count", fontsize=10, fontweight='bold')
    ax_b.set_title("B. Protein Length Distribution", fontsize=11, fontweight='bold', loc='left')
    ax_b.grid(axis='y', alpha=0.3)
    
    # C: IDR vs LLPS scatter
    ax_c = fig.add_subplot(gs[1, 0])
    scatter = ax_c.scatter(summary["mean_iupred"], summary["mean_llps"],
                          s=summary["variant_count"]*5 + 30, alpha=0.6, c=summary["mean_conservation"],
                          cmap='viridis', edgecolors='black', linewidth=0.5)
    ax_c.set_xlabel("Mean IDR Score", fontsize=10, fontweight='bold')
    ax_c.set_ylabel("Mean LLPS Propensity", fontsize=10, fontweight='bold')
    ax_c.set_title("C. IDR vs LLPS (colored by conservation)", fontsize=11, fontweight='bold', loc='left')
    ax_c.grid(alpha=0.3)
    cbar = plt.colorbar(scatter, ax=ax_c)
    cbar.set_label("Conservation", fontsize=9)
    
    # D: IDR score distribution
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.hist(summary["mean_iupred"], bins=15, color='#d73027', edgecolor='black', linewidth=0.7, alpha=0.8)
    ax_d.set_xlabel("Mean IDR Score", fontsize=10, fontweight='bold')
    ax_d.set_ylabel("Count", fontsize=10, fontweight='bold')
    ax_d.set_title("D. IDR Score Distribution", fontsize=11, fontweight='bold', loc='left')
    ax_d.grid(axis='y', alpha=0.3)
    
    # E: Conservation distribution
    ax_e = fig.add_subplot(gs[1, 2])
    ax_e.hist(summary["mean_conservation"], bins=15, color='#1a9850', edgecolor='black', linewidth=0.7, alpha=0.8)
    ax_e.set_xlabel("Mean Conservation Score", fontsize=10, fontweight='bold')
    ax_e.set_ylabel("Count", fontsize=10, fontweight='bold')
    ax_e.set_title("E. Conservation Distribution", fontsize=11, fontweight='bold', loc='left')
    ax_e.grid(axis='y', alpha=0.3)
    
    # F: ClinVar variant distribution
    ax_f = fig.add_subplot(gs[2, :])
    genes_sorted = summary.sort_values("variant_count", ascending=False).head(20)
    colors_cv = plt.cm.RdYlGn_r(np.linspace(0.2, 0.8, len(genes_sorted)))
    ax_f.barh(genes_sorted["gene"], genes_sorted["variant_count"], color=colors_cv, edgecolor='black', linewidth=0.7)
    ax_f.set_xlabel("ClinVar Variants", fontsize=11, fontweight='bold')
    ax_f.set_title("F. Top 20 Genes by ClinVar Variant Count", fontsize=12, fontweight='bold', loc='left')
    ax_f.grid(axis='x', alpha=0.3)
    ax_f.invert_yaxis()
    
    plt.suptitle("Cancer Protein IDR Atlas: Global Overview", fontsize=14, fontweight='bold', y=0.995)
    fig.savefig(output_dir / "figure_1_overview.png", dpi=300, bbox_inches='tight')
    print(f"✓ Saved figure_1_overview.png")


def figure_2_representative_genes(residue_data: dict, output_dir: Path):
    """Figure 2: Representative per-residue profiles (TP53, KRAS, BRCA1)."""
    fig, axes = plt.subplots(3, 1, figsize=(14, 10))
    genes = ["TP53", "KRAS", "BRCA1"]
    
    for idx, gene in enumerate(genes):
        ax = axes[idx]
        if gene not in residue_data:
            continue
        
        df = residue_data[gene]
        pos = df["pos"].values
        iupred = df["iupred_score"].values
        llps = df["llps_proxy"].values
        cons = df["conservation"].values
        
        ax.fill_between(pos, 0, iupred, alpha=0.4, color='#d73027', label='IDR (IUPred2A)')
        ax.plot(pos, iupred, color='#d73027', linewidth=1.5)
        
        ax2 = ax.twinx()
        ax2.plot(pos, llps, color='#7b2d8b', linewidth=1.5, linestyle='--', label='LLPS Propensity')
        ax2.plot(pos, cons, color='#1a9850', linewidth=1.5, linestyle=':', label='Conservation')
        
        ax.set_xlabel("Amino Acid Position", fontsize=10, fontweight='bold')
        ax.set_ylabel("IDR Score", fontsize=10, fontweight='bold', color='#d73027')
        ax2.set_ylabel("LLPS / Conservation", fontsize=10, fontweight='bold')
        ax.set_title(f"{gene} ({len(df)} residues) — Per-Residue Predictions", fontsize=11, fontweight='bold')
        ax.grid(alpha=0.3)
        ax.tick_params(axis='y', labelcolor='#d73027')
        
        # Add legend
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1+lines2, labels1+labels2, loc='upper right', fontsize=9)
    
    plt.tight_layout()
    fig.savefig(output_dir / "figure_2_representative_genes.png", dpi=300, bbox_inches='tight')
    print(f"✓ Saved figure_2_representative_genes.png")


def figure_3_heatmap(summary: pd.DataFrame, output_dir: Path):
    """Figure 3: Correlation heatmap of atlas metrics."""
    metrics = summary[["mean_iupred", "mean_llps", "mean_conservation", 
                       "mean_structural", "low_complexity_fraction"]].copy()
    metrics.columns = ["IDR", "LLPS", "Conservation", "Structural", "Low Complexity"]
    
    corr = metrics.corr()
    
    fig, ax = plt.subplots(figsize=(8, 7))
    sns.heatmap(corr, annot=True, fmt='.2f', cmap='coolwarm', center=0, square=True,
                linewidths=1, cbar_kws={"shrink": 0.8}, ax=ax, vmin=-1, vmax=1)
    ax.set_title("Correlation Matrix of Atlas Metrics Across 50 Genes", fontsize=12, fontweight='bold', pad=14)
    
    plt.tight_layout()
    fig.savefig(output_dir / "figure_3_correlation_heatmap.png", dpi=300, bbox_inches='tight')
    print(f"✓ Saved figure_3_correlation_heatmap.png")


def figure_4_metrics_comparison(summary: pd.DataFrame, output_dir: Path):
    """Figure 4: Multi-metric comparison boxplots."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    
    metrics = [
        ("mean_iupred", "IDR Score", "#d73027"),
        ("mean_llps", "LLPS Propensity", "#7b2d8b"),
        ("mean_conservation", "Conservation Score", "#1a9850"),
        ("mean_structural", "Structural Proxy", "#2166ac"),
    ]
    
    for idx, (col, label, color) in enumerate(metrics):
        ax = axes[idx // 2, idx % 2]
        data = summary[col].values
        bp = ax.boxplot(data, patch_artist=True, widths=0.5, orientation='vertical')
        for patch in bp['boxes']:
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        ax.scatter([1]*len(data), data, alpha=0.3, s=40, color='black', edgecolors='black', linewidth=0.5)
        ax.set_ylabel(label, fontsize=11, fontweight='bold')
        ax.set_title(f"{label}", fontsize=11, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        ax.set_xticklabels([])
        
        # Add statistics
        stats_text = f"μ={data.mean():.3f}\nσ={data.std():.3f}\nmin={data.min():.3f}\nmax={data.max():.3f}"
        ax.text(0.98, 0.97, stats_text, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    fig.suptitle("Distribution of Atlas Metrics Across 50 Cancer Genes", fontsize=13, fontweight='bold')
    plt.tight_layout()
    fig.savefig(output_dir / "figure_4_metrics_comparison.png", dpi=300, bbox_inches='tight')
    print(f"✓ Saved figure_4_metrics_comparison.png")


def figure_7_structural_library(atlas_dir: Path, output_dir: Path):
    """Figure 7: Structural-library overview of IDR residues."""
    rows = []
    for atlas_path in sorted(atlas_dir.glob("*_atlas.tsv")):
        if atlas_path.name == "global_disorder_phylogeny_atlas.tsv":
            continue
        df = pd.read_csv(atlas_path, sep="\t")
        if df.empty:
            continue
        df = df.copy()
        df = df[df["iupred_score"].notna() & df["llps_proxy"].notna() & df["structural_proxy"].notna()]
        if df.empty:
            continue
        df["gene"] = df["gene"].astype(str)
        df["state"] = df.apply(
            lambda r: "condensate" if r["llps_proxy"] >= 0.75 and r["iupred_score"] >= 0.65
            else "expanded" if r["llps_proxy"] >= 0.45 and r["structural_proxy"] <= 0.35
            else "compact",
            axis=1,
        )
        df["quantum_biophysical_score"] = (
            0.35 * df["iupred_score"] + 0.35 * df["llps_proxy"] + 0.15 * df["structural_proxy"] + 0.15 * df.get("quantum_proxy", 0.0)
        ).clip(0.0, 1.0)
        rows.append(df)

    if not rows:
        return

    atlas_df = pd.concat(rows, ignore_index=True)
    atlas_df = atlas_df.sample(n=min(4000, len(atlas_df)), random_state=7)
    virus_count_series = atlas_df["virus_count"].fillna(0) if "virus_count" in atlas_df.columns else pd.Series(0, index=atlas_df.index)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    ax1, ax2, ax3, ax4 = axes.flatten()

    state_counts = atlas_df["state"].value_counts().reindex(["compact", "expanded", "condensate"]).fillna(0)
    state_counts.plot(kind="bar", ax=ax1, color=["#4C78A8", "#F58518", "#54A24B"], edgecolor="black")
    ax1.set_title("A. LLPS state distribution")
    ax1.set_ylabel("Residue count")
    ax1.set_xlabel("State")

    scatter = ax2.scatter(
        atlas_df["iupred_score"],
        atlas_df["llps_proxy"],
        c=atlas_df["quantum_biophysical_score"],
        s=20 + virus_count_series * 35,
        cmap="viridis",
        alpha=0.75,
        edgecolors="black",
        linewidth=0.3,
    )
    ax2.set_title("B. Disorder vs LLPS with biophysical score")
    ax2.set_xlabel("IDR score")
    ax2.set_ylabel("LLPS propensity")
    cbar = fig.colorbar(scatter, ax=ax2, pad=0.04)
    cbar.set_label("Quantum-biophysical score")

    ax3.hist(atlas_df["quantum_biophysical_score"], bins=20, color="#6B5B95", edgecolor="black")
    ax3.set_title("C. Quantum-biophysical score distribution")
    ax3.set_xlabel("Score")
    ax3.set_ylabel("Residue count")

    top_genes = atlas_df.groupby("gene")["quantum_biophysical_score"].mean().sort_values(ascending=False).head(10)
    top_genes.plot(kind="bar", ax=ax4, color="#C44E52", edgecolor="black")
    ax4.set_title("D. Top genes by mean score")
    ax4.set_ylabel("Mean score")
    ax4.set_xlabel("Gene")
    ax4.tick_params(axis="x", rotation=45)

    fig.suptitle("Figure 7. Structural-library overview of IDR residues", fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0.02, 1, 0.98])
    fig.savefig(output_dir / "figure_7_structural_library.png", dpi=300, bbox_inches="tight")
    print("✓ Saved figure_7_structural_library.png")


def figure_6_entrez_chromosome_loci(summary: pd.DataFrame, output_dir: Path):
    """Figure 6: Entrez-derived chromosome loci with IDR and LLPS overlays."""
    chrom_order = [str(i) for i in range(1, 23)] + ["X", "Y"]
    y_map = {chrom: i for i, chrom in enumerate(chrom_order)}

    rows = []
    for _, row in summary.iterrows():
        gene = str(row["gene"])
        locus = GENE_LOCI.get(gene)
        if not locus:
            continue
        chrom = str(locus["chrom"])
        if chrom not in y_map:
            continue
        start = int(locus["start"])
        end = int(locus["end"])
        mid_mb = (start + end) / 2.0 / 1_000_000.0
        chrom_len_mb = CHROM_SIZES.get(chrom, 0) / 1_000_000.0
        rows.append({
            "gene": gene,
            "chrom": chrom,
            "mid_mb": mid_mb,
            "chrom_len_mb": chrom_len_mb,
            "mean_iupred": float(row["mean_iupred"]),
            "mean_llps": float(row["mean_llps"]),
            "variant_count": float(row.get("variant_count", 0.0)),
        })

    if not rows:
        return

    loci_df = pd.DataFrame(rows)
    max_len_mb = max(v / 1_000_000.0 for v in CHROM_SIZES.values())

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 12), sharex=True)

    for ax, color_col, cmap, panel_title in (
        (ax1, "mean_iupred", "Reds", "A. Entrez GRCh38 loci colored by mean IDR"),
        (ax2, "mean_llps", "Purples", "B. Entrez GRCh38 loci colored by mean LLPS"),
    ):
        for chrom in chrom_order:
            y = y_map[chrom]
            chrom_len_mb = CHROM_SIZES[chrom] / 1_000_000.0
            ax.plot([0, chrom_len_mb], [y, y], color="#c7ccd9", linewidth=8, solid_capstyle='round', zorder=1)

        sizes = 60 + (loci_df["mean_llps"] * 320.0)
        scatter = ax.scatter(
            loci_df["mid_mb"],
            loci_df["chrom"].map(y_map),
            c=loci_df[color_col],
            s=sizes,
            cmap=cmap,
            edgecolors="#1f2a44",
            linewidth=0.5,
            alpha=0.92,
            zorder=3,
        )

        top_labeled = loci_df.sort_values("variant_count", ascending=False).head(12)
        for _, r in top_labeled.iterrows():
            ax.annotate(
                r["gene"],
                xy=(r["mid_mb"], y_map[r["chrom"]]),
                xytext=(4, 3),
                textcoords="offset points",
                fontsize=8,
                color="#25304b",
            )

        cbar = fig.colorbar(scatter, ax=ax, shrink=0.75, pad=0.01)
        cbar.set_label("Mean score", fontsize=9)

        ax.set_yticks(list(y_map.values()))
        ax.set_yticklabels(chrom_order, fontsize=9)
        ax.set_ylabel("Chromosome", fontsize=10, fontweight='bold')
        ax.set_title(panel_title, fontsize=12, fontweight='bold', loc='left')
        ax.grid(axis='x', alpha=0.25, linestyle='--')
        ax.set_xlim(-2, max_len_mb + 6)
        ax.invert_yaxis()

    ax2.set_xlabel("Genomic position on chromosome (Mb)", fontsize=11, fontweight='bold')

    legend_handles = [
        mpatches.Patch(color="#f2f2f2", label="Locus coordinates: NCBI Entrez Gene (GRCh38)"),
    ]
    ax2.legend(handles=legend_handles, loc='lower right', fontsize=9, frameon=True)

    fig.suptitle(
        "Figure 6. Entrez Chromosome Locus Maps with IDR and LLPS Context",
        fontsize=14,
        fontweight='bold',
        y=0.995,
    )
    fig.text(
        0.01,
        0.005,
        "Gene loci are mapped from NCBI Entrez Gene coordinates (GRCh38). Marker size scales with mean LLPS propensity.",
        fontsize=9,
        color="#445",
    )
    plt.tight_layout(rect=[0, 0.02, 1, 0.98])
    fig.savefig(output_dir / "figure_6_entrez_chromosome_loci.png", dpi=300, bbox_inches='tight')
    print("✓ Saved figure_6_entrez_chromosome_loci.png")


def figure_8_ptm_category_distribution(ptm_category_path: Path, output_dir: Path):
    """Figure 8: PTM category distribution across the cancer-gene cohort."""
    if not ptm_category_path.exists():
        print(f"! Skipping Figure 8 (missing PTM category table): {ptm_category_path}")
        return

    df = pd.read_csv(ptm_category_path, sep="\t")
    if df.empty:
        print("! Skipping Figure 8 (PTM category table is empty)")
        return

    if not {"gene", "category", "count"}.issubset(df.columns):
        print("! Skipping Figure 8 (PTM category table missing required columns)")
        return

    category_order = [c for c in df["category"].dropna().unique() if str(c).strip()]
    if not category_order:
        print("! Skipping Figure 8 (no PTM categories available)")
        return

    summary = (
        df.groupby("category", as_index=False)["count"]
        .sum()
        .sort_values("count", ascending=False)
    )
    summary = summary.head(12)

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(data=summary, x="category", y="count", palette="viridis", ax=ax)
    ax.set_title("Figure 8. PTM category distribution across cancer genes", fontsize=12, fontweight="bold")
    ax.set_xlabel("PTM category")
    ax.set_ylabel("Total motif count")
    ax.grid(axis="y", alpha=0.25)
    plt.xticks(rotation=20)
    plt.tight_layout()
    fig.savefig(output_dir / "figure_9_ptm_category_distribution.png", dpi=300, bbox_inches="tight")
    print("✓ Saved figure_9_ptm_category_distribution.png")


def figure_9_ptm_pathway_differential(ptm_pathway_path: Path, output_dir: Path):
    """Figure 9: PTM-context pathway differential effects on LLPS mutant impact."""
    if not ptm_pathway_path.exists():
        print(f"! Skipping Figure 8 (missing PTM pathway table): {ptm_pathway_path}")
        return

    df = pd.read_csv(ptm_pathway_path, sep="\t")
    if df.empty:
        print("! Skipping Figure 8 (PTM pathway table is empty)")
        return

    dedup = df.drop_duplicates(subset=["pathway", "direction"]).copy()
    if dedup.empty:
        print("! Skipping Figure 8 (no unique pathway-direction rows)")
        return

    inc = dedup[dedup["direction"] == "increase"].sort_values("ptm_differential_mean_abs_delta", ascending=False)
    dec = dedup[dedup["direction"] == "decrease"].sort_values("ptm_differential_mean_abs_delta", ascending=True)

    top_n = 10
    inc_plot = inc.head(top_n).copy()
    dec_plot = dec.head(top_n).copy()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8), sharex=False)

    if not inc_plot.empty:
        ax1.barh(
            inc_plot["pathway"],
            inc_plot["ptm_differential_mean_abs_delta"],
            color="#2A9D8F",
            edgecolor="black",
            linewidth=0.7,
        )
        ax1.invert_yaxis()
    ax1.set_title("A. Top PTM-associated increase-context signals", fontsize=11, fontweight="bold")
    ax1.set_xlabel("with PTM - without PTM (mean |delta LLPS|)", fontsize=10)
    ax1.set_ylabel("Pathway", fontsize=10)
    ax1.grid(axis="x", alpha=0.25)

    if not dec_plot.empty:
        ax2.barh(
            dec_plot["pathway"],
            dec_plot["ptm_differential_mean_abs_delta"],
            color="#E76F51",
            edgecolor="black",
            linewidth=0.7,
        )
        ax2.invert_yaxis()
    ax2.set_title("B. Top PTM-associated decrease-context signals", fontsize=11, fontweight="bold")
    ax2.set_xlabel("with PTM - without PTM (mean |delta LLPS|)", fontsize=10)
    ax2.set_ylabel("Pathway", fontsize=10)
    ax2.grid(axis="x", alpha=0.25)

    fig.suptitle("Figure 9. PTM-stratified pathway differential effects", fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0.02, 1, 0.96])
    fig.savefig(output_dir / "figure_8_ptm_pathway_differential.png", dpi=300, bbox_inches="tight")
    print("✓ Saved figure_8_ptm_pathway_differential.png")


def _assign_llp_phase(series: pd.Series) -> pd.Series:
    """Assign LLPS phases from mutant baseline LLPS using robust tertile bins."""
    labels = ["Phase I (lower LLP)", "Phase II (mid LLP)", "Phase III (higher LLP)"]
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.dropna().empty:
        return pd.Series([labels[1]] * len(series), index=series.index)

    uniq = numeric.dropna().nunique()
    if uniq >= 3:
        try:
            phased = pd.qcut(numeric, q=3, labels=labels)
            return phased.astype(str).replace("nan", labels[1])
        except ValueError:
            pass

    vmin = float(numeric.min())
    vmax = float(numeric.max())
    if np.isclose(vmin, vmax):
        return pd.Series([labels[1]] * len(series), index=series.index)
    bins = np.linspace(vmin, vmax, 4)
    phased = pd.cut(numeric, bins=bins, include_lowest=True, labels=labels)
    return phased.astype(str).replace("nan", labels[1])


def figure_10_viable_mutant_llp_phase_differential(mutants_path: Path, output_dir: Path):
    """Figure 10: Viable IDR mutant effects by LLP phase and PTM/pathway differentials."""
    if not mutants_path.exists():
        print(f"! Skipping Figure 9 (missing mutant table): {mutants_path}")
        return

    mutants = pd.read_csv(mutants_path, sep="\t")
    if mutants.empty:
        print("! Skipping Figure 9 (mutant table is empty)")
        return

    required = {"wt_llps_proxy", "delta_llps", "direction", "ptm_context", "pathways"}
    if not required.issubset(mutants.columns):
        missing = sorted(required - set(mutants.columns))
        print(f"! Skipping Figure 9 (missing columns: {missing})")
        return

    df = mutants.copy()
    df["wt_llps_proxy"] = pd.to_numeric(df["wt_llps_proxy"], errors="coerce")
    df["delta_llps"] = pd.to_numeric(df["delta_llps"], errors="coerce")
    df = df.dropna(subset=["wt_llps_proxy", "delta_llps", "direction", "pathways"])
    if df.empty:
        print("! Skipping Figure 9 (no numeric mutant records)")
        return

    df["abs_delta_llps"] = df["delta_llps"].abs()
    viability_cut = float(df["abs_delta_llps"].quantile(0.90))
    df["viable_mutant"] = df["abs_delta_llps"] >= viability_cut
    df["llp_phase"] = _assign_llp_phase(df["wt_llps_proxy"])
    df["ptm_context"] = df["ptm_context"].astype(bool)

    phase_order = ["Phase I (lower LLP)", "Phase II (mid LLP)", "Phase III (higher LLP)"]
    direction_order = ["increase", "decrease"]

    fig, axes = plt.subplots(2, 2, figsize=(16, 11))
    ax1, ax2, ax3, ax4 = axes.flatten()

    # Panel A: viable fraction by LLP phase and direction
    viable_rate = (
        df.groupby(["llp_phase", "direction"], as_index=False)
        .agg(
            viable_fraction=("viable_mutant", "mean"),
            n_mutants=("viable_mutant", "size"),
        )
    )
    sns.barplot(
        data=viable_rate,
        x="llp_phase",
        y="viable_fraction",
        hue="direction",
        order=phase_order,
        hue_order=direction_order,
        palette={"increase": "#2A9D8F", "decrease": "#E76F51"},
        ax=ax1,
    )
    ax1.set_title("A. Viable mutant fraction by LLP phase and mutation direction", fontsize=11, fontweight="bold")
    ax1.set_xlabel("LLP phase (IDR residue baseline)")
    ax1.set_ylabel("Viable mutant fraction")
    ax1.grid(axis="y", alpha=0.25)
    ax1.tick_params(axis="x", rotation=12)

    # Panel B: mean absolute delta by LLP phase and PTM context
    phase_ptm = (
        df.groupby(["llp_phase", "ptm_context"], as_index=False)
        .agg(
            mean_abs_delta=("abs_delta_llps", "mean"),
            n_mutants=("abs_delta_llps", "size"),
        )
    )
    phase_ptm["ptm_context_label"] = np.where(phase_ptm["ptm_context"], "with PTM", "without PTM")
    sns.barplot(
        data=phase_ptm,
        x="llp_phase",
        y="mean_abs_delta",
        hue="ptm_context_label",
        order=phase_order,
        hue_order=["with PTM", "without PTM"],
        palette={"with PTM": "#6A4C93", "without PTM": "#577590"},
        ax=ax2,
    )
    ax2.set_title("B. Mean |delta LLPS| by LLP phase and PTM context", fontsize=11, fontweight="bold")
    ax2.set_xlabel("LLP phase (IDR residue baseline)")
    ax2.set_ylabel("Mean absolute LLPS shift")
    ax2.grid(axis="y", alpha=0.25)
    ax2.tick_params(axis="x", rotation=12)

    # Panel C: pathway PTM differential heatmap across LLP phases and directions
    exploded = df.assign(pathway=df["pathways"].astype(str).str.split(";")).explode("pathway")
    exploded = exploded[exploded["pathway"].notna() & (exploded["pathway"].astype(str).str.strip() != "")]
    diff_rows = []
    for (pathway, phase, direction), grp in exploded.groupby(["pathway", "llp_phase", "direction"]):
        with_ptm = grp[grp["ptm_context"]]["abs_delta_llps"]
        without_ptm = grp[~grp["ptm_context"]]["abs_delta_llps"]
        if len(with_ptm) == 0 or len(without_ptm) == 0:
            continue
        diff_rows.append({
            "pathway": pathway,
            "llp_phase": phase,
            "direction": direction,
            "ptm_diff": float(with_ptm.mean() - without_ptm.mean()),
        })

    diff_df = pd.DataFrame(diff_rows)
    if not diff_df.empty:
        diff_df["phase_direction"] = diff_df["llp_phase"] + " | " + diff_df["direction"].str.capitalize()
        top_pathways = (
            diff_df.groupby("pathway")["ptm_diff"]
            .apply(lambda s: float(np.max(np.abs(s))))
            .sort_values(ascending=False)
            .head(10)
            .index
        )
        heatmap_df = diff_df[diff_df["pathway"].isin(top_pathways)].pivot_table(
            index="pathway",
            columns="phase_direction",
            values="ptm_diff",
            aggfunc="mean",
        )
        desired_cols = [
            f"{phase} | {direction.capitalize()}"
            for phase in phase_order
            for direction in direction_order
        ]
        heatmap_df = heatmap_df.reindex(columns=[c for c in desired_cols if c in heatmap_df.columns])
        sns.heatmap(
            heatmap_df,
            cmap="coolwarm",
            center=0.0,
            linewidths=0.4,
            linecolor="white",
            cbar_kws={"label": "PTM differential mean |delta LLPS|\n(with PTM - without PTM)"},
            ax=ax3,
        )
        ax3.set_title("C. Pathway PTM differentials by LLP phase and direction", fontsize=11, fontweight="bold")
        ax3.set_xlabel("Phase-direction condition")
        ax3.set_ylabel("Pathway")
    else:
        ax3.axis("off")
        ax3.text(0.5, 0.5, "Insufficient PTM/non-PTM contrast\nfor pathway-phase heatmap", ha="center", va="center")

    # Panel D: viable mutant counts by pathway and LLP phase
    viable_exploded = exploded[exploded["viable_mutant"]].copy()
    if not viable_exploded.empty:
        pathway_phase_counts = (
            viable_exploded.groupby(["pathway", "llp_phase"], as_index=False)
            .size()
            .rename(columns={"size": "viable_count"})
        )
        top_viable_pathways = (
            pathway_phase_counts.groupby("pathway")["viable_count"]
            .sum()
            .sort_values(ascending=False)
            .head(8)
            .index
        )
        stacked = pathway_phase_counts[pathway_phase_counts["pathway"].isin(top_viable_pathways)].pivot_table(
            index="pathway",
            columns="llp_phase",
            values="viable_count",
            aggfunc="sum",
            fill_value=0,
        )
        stacked = stacked.reindex(columns=[c for c in phase_order if c in stacked.columns])
        stacked.plot(
            kind="bar",
            stacked=True,
            color=["#577590", "#43AA8B", "#F8961E"],
            edgecolor="black",
            linewidth=0.4,
            ax=ax4,
        )
        ax4.set_title("D. Top pathways by viable mutant burden across LLP phases", fontsize=11, fontweight="bold")
        ax4.set_xlabel("Pathway")
        ax4.set_ylabel("Viable mutant count")
        ax4.tick_params(axis="x", rotation=35)
        ax4.grid(axis="y", alpha=0.25)
        ax4.legend(title="LLP phase", fontsize=8, title_fontsize=9, loc="upper right")
    else:
        ax4.axis("off")
        ax4.text(0.5, 0.5, "No viable mutants under\ncurrent quantile threshold", ha="center", va="center")

    fig.suptitle(
        (
            "Figure 10. Viable IDR mutant effects on phase separation and LLP-phase PTM/pathway differentials\n"
            f"Viable threshold = top 10% |delta LLPS| (cutoff={viability_cut:.4f})"
        ),
        fontsize=13,
        fontweight="bold",
    )
    plt.tight_layout(rect=[0, 0.02, 1, 0.95])
    fig.savefig(output_dir / "figure_9_viable_mutant_llp_phase_differential.png", dpi=300, bbox_inches="tight")
    print("✓ Saved figure_9_viable_mutant_llp_phase_differential.png")


def figure_5_workflow_schema(output_dir: Path):
    """Figure 5: End-to-end computational workflow schema."""
    fig, ax = plt.subplots(figsize=(14, 9))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 9)
    ax.set_aspect("equal")
    ax.axis("off")

    box_style = dict(boxstyle="round,pad=0.4", linewidth=1.5)
    arrow_kw = dict(arrowstyle="-|>", lw=1.8, color="#2d5a8e")

    layers = [
        (1.5, 7.8, "INPUT LAYER", "#1e3a5f", "#fff"),
        (1.5, 5.4, "PREDICTOR LAYER", "#2d5a8e", "#fff"),
        (1.5, 3.0, "INTEGRATION LAYER", "#d73027", "#fff"),
        (1.5, 0.8, "OUTPUT LAYER", "#1a9850", "#fff"),
    ]
    for x, y, label, fc, tc in layers:
        ax.text(x, y, label, fontsize=9, fontweight="bold", color=tc,
                bbox=dict(facecolor=fc, edgecolor="none", boxstyle="round,pad=0.25", alpha=0.9),
                ha="center", va="center")

    inputs = [
        (4.0, 8.2, "Gene list\n(50 CGC genes)"),
        (6.5, 8.2, "FASTA\nsequences"),
        (9.0, 8.2, "ClinVar\nvariants"),
        (11.5, 8.2, "Oncovirus\ncurated table"),
    ]
    for x, y, label in inputs:
        ax.text(x, y, label, fontsize=8.5, ha="center", va="center",
                bbox=dict(facecolor="#e8edf3", edgecolor="#1e3a5f", **box_style))

    predictors = [
        (3.5, 5.8, "IUPred2A\n(disorder)"),
        (5.5, 5.8, "LLPS proxy\n(phase sep)"),
        (7.5, 5.8, "SEG / PLAAC\n(low complexity)"),
        (9.5, 5.8, "AlphaFold\n(pLDDT)"),
        (11.5, 5.8, "Phylogeny\n(conservation)"),
        (5.0, 4.8, "SLiM / PTM\nscanning"),
        (7.5, 4.8, "Delta-LLPS\nmutant map"),
        (10.0, 4.8, "Structural\nproxy"),
    ]
    for x, y, label in predictors:
        ax.text(x, y, label, fontsize=8, ha="center", va="center",
                bbox=dict(facecolor="#dce4f0", edgecolor="#2d5a8e", **box_style))

    integrations = [
        (4.5, 3.3, "Per-residue\natlas builder"),
        (7.5, 3.3, "VIPP scoring\n(IDR+LLPS+virus)"),
        (10.5, 3.3, "Benchmark\nablation"),
    ]
    for x, y, label in integrations:
        ax.text(x, y, label, fontsize=8.5, ha="center", va="center",
                bbox=dict(facecolor="#fde0dd", edgecolor="#d73027", **box_style))

    int2 = [
        (4.5, 2.3, "ClinVar\naugmentation"),
        (7.5, 2.3, "IDR mutant\npathway analysis"),
        (10.5, 2.3, "Structural\nlibrary"),
    ]
    for x, y, label in int2:
        ax.text(x, y, label, fontsize=8, ha="center", va="center",
                bbox=dict(facecolor="#fde0dd", edgecolor="#d73027", **box_style))

    outputs = [
        (3.5, 1.0, "Per-gene atlas\nTSV tables"),
        (6.0, 1.0, "Global summary\n+ hotspot tables"),
        (8.5, 1.0, "Interactive\nchromosome map"),
        (11.0, 1.0, "Publication\nfigures (1–10)"),
    ]
    for x, y, label in outputs:
        ax.text(x, y, label, fontsize=8.5, ha="center", va="center",
                bbox=dict(facecolor="#d4edda", edgecolor="#1a9850", **box_style))

    for inp_x, _, _ in inputs:
        ax.annotate("", xy=(inp_x, 6.3), xytext=(inp_x, 7.7),
                     arrowprops=dict(**arrow_kw, connectionstyle="arc3,rad=0"))
    for px, _, _ in integrations:
        ax.annotate("", xy=(px, 1.7), xytext=(px, 2.8),
                     arrowprops=dict(**arrow_kw, connectionstyle="arc3,rad=0"))

    fig.suptitle("Figure 5. End-to-end computational workflow schema",
                 fontsize=14, fontweight="bold", y=0.98)
    fig.text(0.5, 0.01,
             "Reproducible pipeline from sequence inputs through predictor and integration layers "
             "to atlas construction, benchmarking, and publication artifacts.",
             fontsize=9, ha="center", color="#555")
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    fig.savefig(output_dir / "figure_5_workflow_schema.png", dpi=300, bbox_inches="tight")
    print("✓ Saved figure_5_workflow_schema.png")


def figure_vipp_model(output_dir: Path):
    """Figure 8 (VIPP scoring model) as SVG."""
    svg = """<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 720 380" width="720" height="380"
     font-family="system-ui, -apple-system, sans-serif">
  <defs>
    <linearGradient id="bg" x1="0" y1="0" x2="0" y2="1">
      <stop offset="0%" stop-color="#f8f9fb"/>
      <stop offset="100%" stop-color="#eef1f5"/>
    </linearGradient>
  </defs>
  <rect width="720" height="380" rx="12" fill="url(#bg)" stroke="#c8cdd8" stroke-width="1"/>
  <text x="360" y="32" text-anchor="middle" font-size="16" font-weight="bold" fill="#1e3a5f">
    Figure 8. VIPP Scoring Model for Residue Prioritization
  </text>

  <!-- Input boxes -->
  <rect x="30" y="70" width="140" height="60" rx="8" fill="#d73027" opacity="0.85"/>
  <text x="100" y="96" text-anchor="middle" font-size="12" font-weight="bold" fill="#fff">IDR Disorder</text>
  <text x="100" y="114" text-anchor="middle" font-size="11" fill="#fdd">I(i) — IUPred2A</text>

  <rect x="200" y="70" width="140" height="60" rx="8" fill="#7b2d8b" opacity="0.85"/>
  <text x="270" y="96" text-anchor="middle" font-size="12" font-weight="bold" fill="#fff">LLPS Propensity</text>
  <text x="270" y="114" text-anchor="middle" font-size="11" fill="#ede">L(i) — proxy</text>

  <rect x="370" y="70" width="140" height="60" rx="8" fill="#2166ac" opacity="0.85"/>
  <text x="440" y="96" text-anchor="middle" font-size="12" font-weight="bold" fill="#fff">Virus Contact</text>
  <text x="440" y="114" text-anchor="middle" font-size="11" fill="#cde">V(i) — binary</text>

  <!-- Weight labels -->
  <text x="100" y="155" text-anchor="middle" font-size="13" font-weight="bold" fill="#d73027">w = 0.35</text>
  <text x="270" y="155" text-anchor="middle" font-size="13" font-weight="bold" fill="#7b2d8b">w = 0.35</text>
  <text x="440" y="155" text-anchor="middle" font-size="13" font-weight="bold" fill="#2166ac">w = 0.30</text>

  <!-- Arrows down -->
  <line x1="100" y1="160" x2="100" y2="190" stroke="#d73027" stroke-width="2.5" marker-end="url(#arrowR)"/>
  <line x1="270" y1="160" x2="270" y2="190" stroke="#7b2d8b" stroke-width="2.5"/>
  <line x1="440" y1="160" x2="440" y2="190" stroke="#2166ac" stroke-width="2.5"/>

  <!-- Summation box -->
  <rect x="60" y="190" width="450" height="50" rx="10" fill="#1e3a5f" opacity="0.9"/>
  <text x="285" y="215" text-anchor="middle" font-size="13" fill="#fff" font-weight="bold">
    Weighted sum: VIPP(i) = 0.35 · I(i) + 0.35 · L(i) + 0.30 · V(i)
  </text>
  <text x="285" y="232" text-anchor="middle" font-size="11" fill="#aac">
    Clipped to [0, 1]
  </text>

  <!-- Arrow to output -->
  <line x1="285" y1="240" x2="285" y2="270" stroke="#1e3a5f" stroke-width="2.5"/>

  <!-- Output box -->
  <rect x="155" y="270" width="260" height="50" rx="10" fill="#1a9850" opacity="0.9"/>
  <text x="285" y="296" text-anchor="middle" font-size="13" font-weight="bold" fill="#fff">
    VIPP Score (per residue)
  </text>
  <text x="285" y="312" text-anchor="middle" font-size="10" fill="#cec">
    Prioritization for experimental follow-up
  </text>

  <!-- Side note -->
  <rect x="540" y="70" width="160" height="250" rx="8" fill="#f0f4f8" stroke="#c8cdd8"/>
  <text x="620" y="92" text-anchor="middle" font-size="11" font-weight="bold" fill="#1e3a5f">
    Supporting layers
  </text>
  <text x="555" y="115" font-size="10" fill="#555">• Conservation</text>
  <text x="555" y="135" font-size="10" fill="#555">• Structural proxy</text>
  <text x="555" y="155" font-size="10" fill="#555">• ClinVar variants</text>
  <text x="555" y="175" font-size="10" fill="#555">• Low complexity</text>
  <text x="555" y="195" font-size="10" fill="#555">• AlphaFold pLDDT</text>
  <text x="555" y="215" font-size="10" fill="#555">• SLiM / PTM motifs</text>
  <text x="555" y="235" font-size="10" fill="#555">• Delta-LLPS map</text>
  <text x="555" y="270" font-size="9" fill="#888">
    These layers inform
  </text>
  <text x="555" y="284" font-size="9" fill="#888">
    atlas context but do
  </text>
  <text x="555" y="298" font-size="9" fill="#888">
    not enter the VIPP
  </text>
  <text x="555" y="312" font-size="9" fill="#888">
    weighted sum directly.
  </text>

  <text x="360" y="365" text-anchor="middle" font-size="9" fill="#888">
    Default heuristic weights (current release). Future calibrated model will replace these.
  </text>
</svg>"""
    (output_dir / "figure_vipp_model.svg").write_text(svg, encoding="utf-8")
    print("✓ Saved figure_vipp_model.svg")


def main():
    parser = argparse.ArgumentParser(description="Generate publication figures for the Cancer Protein IDR Atlas")
    parser.add_argument("--atlas-dir", default="results/atlas")
    parser.add_argument("--summary-path", default="results/atlas/global_disorder_phylogeny_atlas.tsv")
    parser.add_argument("--ptm-pathway-path", default="results/mutants/pathway_ptm_differential_llps.tsv")
    parser.add_argument("--ptm-category-path", default="results/slim/ptm_category_summary.tsv")
    parser.add_argument("--mutants-path", default="results/mutants/idr_llps_mutants.tsv")
    parser.add_argument("--output-dir", default="results/figures")
    args = parser.parse_args()

    atlas_dir = Path(args.atlas_dir)
    summary_path = Path(args.summary_path)
    output_dir = Path(args.output_dir)
    ptm_pathway_path = Path(args.ptm_pathway_path)
    ptm_category_path = Path(args.ptm_category_path)
    mutants_path = Path(args.mutants_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    summary, residue_data = load_data(atlas_dir, summary_path)

    print("Generating publication-quality figures...")
    figure_1_overview(summary, output_dir)
    figure_2_representative_genes(residue_data, output_dir)
    figure_3_heatmap(summary, output_dir)
    figure_4_metrics_comparison(summary, output_dir)
    figure_5_workflow_schema(output_dir)
    figure_6_entrez_chromosome_loci(summary, output_dir)
    figure_7_structural_library(atlas_dir, output_dir)
    figure_vipp_model(output_dir)
    figure_8_ptm_category_distribution(ptm_category_path, output_dir)
    figure_9_ptm_pathway_differential(ptm_pathway_path, output_dir)
    figure_10_viable_mutant_llp_phase_differential(mutants_path, output_dir)

    print(f"\n✓ All figures saved to {output_dir}/")


if __name__ == "__main__":
    main()
