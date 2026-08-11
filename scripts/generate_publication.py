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
    fig.savefig(output_dir / "figure_8_ptm_category_distribution.png", dpi=300, bbox_inches="tight")
    print("✓ Saved figure_8_ptm_category_distribution.png")


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
    fig.savefig(output_dir / "figure_9_ptm_pathway_differential.png", dpi=300, bbox_inches="tight")
    print("✓ Saved figure_9_ptm_pathway_differential.png")


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
    fig.savefig(output_dir / "figure_10_viable_mutant_llp_phase_differential.png", dpi=300, bbox_inches="tight")
    print("✓ Saved figure_10_viable_mutant_llp_phase_differential.png")


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
    figure_6_entrez_chromosome_loci(summary, output_dir)
    figure_7_structural_library(atlas_dir, output_dir)
    figure_8_ptm_category_distribution(ptm_category_path, output_dir)
    figure_9_ptm_pathway_differential(ptm_pathway_path, output_dir)
    figure_10_viable_mutant_llp_phase_differential(mutants_path, output_dir)

    print(f"\n✓ All figures saved to {output_dir}/")


if __name__ == "__main__":
    main()
