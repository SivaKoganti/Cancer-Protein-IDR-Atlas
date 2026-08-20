#!/usr/bin/env python3
"""
Residue-Level Docking Interaction Figures
Generates detailed visualizations of binding site residue interactions
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
from pathlib import Path

# Publication-quality settings
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.linewidth'] = 1.2

class ResidueInteractionFigures:
    """Generate residue-level docking interaction figures"""

    def __init__(self, output_dir: str = "./figures"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def figure_31_lead_coordination(self):
        """Figure 3.1 Panel D: Lead coordination energetics"""
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Left: Residue distances to lead
        residues = ['Cys-34\nThiol S', 'His-67\nImidazole N',
                   'Asp-108\nCarboxyl O', 'Asp-183\nCarboxyl O']
        distances = [2.3, 2.6, 2.8, 3.1]
        energies = [-4.2, -2.1, -1.8, -0.9]
        colors = ['#E74C3C', '#3498DB', '#2ECC71', '#F39C12']

        x = np.arange(len(residues))
        ax1 = axes[0]
        bars1 = ax1.bar(x - 0.2, distances, 0.4, label='Distance (Å)',
                        color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)
        ax1.set_ylabel('Distance (Å)', fontsize=11, fontweight='bold')
        ax1.set_xlabel('Lead-Coordinating Residue', fontsize=11, fontweight='bold')
        ax1.set_title('Lead Coordination Geometry', fontsize=12, fontweight='bold')
        ax1.set_xticks(x)
        ax1.set_xticklabels(residues, fontsize=9)
        ax1.set_ylim(0, 3.5)
        ax1.grid(axis='y', alpha=0.3)

        for bar, val in zip(bars1, distances):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                    f'{val:.1f}', ha='center', va='bottom', fontweight='bold', fontsize=9)

        # Right: Energetic contribution
        ax2 = axes[1]
        bars2 = ax2.bar(x, energies, color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)
        ax2.set_ylabel('Interaction Energy (kcal/mol)', fontsize=11, fontweight='bold')
        ax2.set_xlabel('Lead-Coordinating Residue', fontsize=11, fontweight='bold')
        ax2.set_title('Energetic Contribution to Lead Binding', fontsize=12, fontweight='bold')
        ax2.set_xticks(x)
        ax2.set_xticklabels(residues, fontsize=9)
        ax2.set_ylim(-5, 0)
        ax2.grid(axis='y', alpha=0.3)
        ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.8)

        for bar, val in zip(bars2, energies):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height - 0.2,
                    f'{val:.1f}', ha='center', va='top', fontweight='bold', fontsize=9)

        plt.tight_layout()
        plt.savefig(self.output_dir / 'Figure_3.1D_Lead_Coordination_Energetics.png',
                   dpi=300, bbox_inches='tight')
        plt.close()

    def figure_32_5fu_residue_heatmap(self):
        """Figure 3.2 Panel C: 5-FU residue interaction heatmap"""
        fig, ax = plt.subplots(figsize=(10, 6))

        residues = ['Lys-199', 'Tyr-150', 'Arg-196', 'Phe-206', 'Trp-214', 'Ile-82']
        conditions = ['Native\nOccurrence', 'Lead-Bound\nOccurrence', 'Interaction\nLoss']

        data = np.array([
            [92, 38, 54],
            [78, 22, 56],
            [68, 15, 53],
            [85, 31, 54],
            [79, 28, 51],
            [71, 25, 46]
        ])

        im = ax.imshow(data, cmap='RdYlGn_r', aspect='auto', vmin=0, vmax=100)

        ax.set_xticks(np.arange(len(conditions)))
        ax.set_yticks(np.arange(len(residues)))
        ax.set_xticklabels(conditions, fontsize=10, fontweight='bold')
        ax.set_yticklabels(residues, fontsize=10, fontweight='bold')

        ax.set_xlabel('Condition', fontsize=11, fontweight='bold')
        ax.set_ylabel('Binding Site Residue', fontsize=11, fontweight='bold')
        ax.set_title('5-FU Residue Interaction Heatmap: Native vs. Lead-Bound HSA',
                    fontsize=12, fontweight='bold')

        # Add text annotations
        for i in range(len(residues)):
            for j in range(len(conditions)):
                text = ax.text(j, i, f'{int(data[i, j])}%',
                             ha="center", va="center", color="black", fontsize=9, fontweight='bold')

        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Occurrence (%)', fontsize=10, fontweight='bold')

        plt.tight_layout()
        plt.savefig(self.output_dir / 'Figure_3.2C_5FU_Interaction_Heatmap.png',
                   dpi=300, bbox_inches='tight')
        plt.close()

    def figure_33_allosteric_pathway_dynamics(self):
        """Figure 3.3 Panel C: B-factor changes along pathway"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        pathway_residues = ['Cys-34\n(Lead Site)', 'Lys-129\n(Pathway 1)',
                           'Asp-183\n(Pathway 2)', 'Trp-214\n(Pathway 3)',
                           'Lys-199\n(Drug Site)']
        b_factor_increase = [18, 14, 16, 12, 8]
        contact_freq_native = [0, 24, 31, 38, 62]
        contact_freq_lead = [0, 68, 72, 81, 18]

        x = np.arange(len(pathway_residues))

        # Left: B-factor changes
        ax1.plot(x, b_factor_increase, 'o-', linewidth=2.5, markersize=10,
                color='#E74C3C', label='B-factor increase (Ų)')
        ax1.fill_between(x, b_factor_increase, alpha=0.3, color='#E74C3C')
        ax1.set_ylabel('B-Factor Increase (Ų)', fontsize=11, fontweight='bold')
        ax1.set_xlabel('Allosteric Pathway Position', fontsize=11, fontweight='bold')
        ax1.set_title('Conformational Flexibility Along Pathway', fontsize=12, fontweight='bold')
        ax1.set_xticks(x)
        ax1.set_xticklabels(pathway_residues, fontsize=9)
        ax1.set_ylim(0, 20)
        ax1.grid(alpha=0.3)

        for i, val in enumerate(b_factor_increase):
            ax1.text(i, val + 0.5, f'{val}', ha='center', va='bottom', fontweight='bold', fontsize=9)

        # Right: Contact frequency changes
        width = 0.35
        bars1 = ax2.bar(x - width/2, contact_freq_native, width, label='Native HSA',
                       color='#3498DB', edgecolor='black', linewidth=1.5, alpha=0.8)
        bars2 = ax2.bar(x + width/2, contact_freq_lead, width, label='Lead-Bound HSA',
                       color='#E74C3C', edgecolor='black', linewidth=1.5, alpha=0.8)

        ax2.set_ylabel('Contact Frequency (%)', fontsize=11, fontweight='bold')
        ax2.set_xlabel('Allosteric Pathway Position', fontsize=11, fontweight='bold')
        ax2.set_title('Residue Contact Frequency: Native vs. Lead-Bound', fontsize=12, fontweight='bold')
        ax2.set_xticks(x)
        ax2.set_xticklabels(pathway_residues, fontsize=9)
        ax2.set_ylim(0, 100)
        ax2.legend(fontsize=10, loc='upper center')
        ax2.grid(axis='y', alpha=0.3)

        plt.tight_layout()
        plt.savefig(self.output_dir / 'Figure_3.3C_Pathway_Dynamics.png',
                   dpi=300, bbox_inches='tight')
        plt.close()

    def figure_33_network_connectivity(self):
        """Figure 3.3 Panel D: Network centrality and connectivity"""
        fig, ax = plt.subplots(figsize=(10, 6))

        path_segments = ['Cys-34\n→ Lys-129', 'Lys-129\n→ Asp-183',
                        'Asp-183\n→ Trp-214', 'Trp-214\n→ Lys-199']
        distances = [15, 7, 6, 3]
        network_degree = [5, 4, 4, 3]
        centrality = [0.73, 0.68, 0.71, 0.62]

        x = np.arange(len(path_segments))

        # Create twin axis for centrality score
        ax1 = ax
        ax2 = ax.twinx()

        bars = ax1.bar(x, distances, width=0.6, label='Distance (Å)',
                      color='#3498DB', edgecolor='black', linewidth=1.5, alpha=0.8)
        line = ax2.plot(x, centrality, 'o-', linewidth=2.5, markersize=10,
                       color='#E74C3C', label='Centrality Score')

        ax1.set_ylabel('Path Distance (Å)', fontsize=11, fontweight='bold', color='#3498DB')
        ax2.set_ylabel('Network Centrality Score', fontsize=11, fontweight='bold', color='#E74C3C')
        ax1.set_xlabel('Allosteric Path Segment', fontsize=11, fontweight='bold')
        ax1.set_title('Network Connectivity and Path Distances', fontsize=12, fontweight='bold')
        ax1.set_xticks(x)
        ax1.set_xticklabels(path_segments, fontsize=9)
        ax1.set_ylim(0, 18)
        ax2.set_ylim(0.5, 0.85)
        ax1.grid(axis='y', alpha=0.3)
        ax1.tick_params(axis='y', labelcolor='#3498DB')
        ax2.tick_params(axis='y', labelcolor='#E74C3C')

        for i, (bar, val) in enumerate(zip(bars, distances)):
            ax1.text(bar.get_x() + bar.get_width()/2, val + 0.3, f'{int(val)}Å',
                    ha='center', va='bottom', fontweight='bold', fontsize=9)
            ax2.text(i, centrality[i] + 0.01, f'{centrality[i]:.2f}',
                    ha='center', va='bottom', fontweight='bold', fontsize=9, color='#E74C3C')

        # Add legend
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=10)

        plt.tight_layout()
        plt.savefig(self.output_dir / 'Figure_3.3D_Network_Connectivity.png',
                   dpi=300, bbox_inches='tight')
        plt.close()

    def generate_all_residue_figures(self):
        """Generate all residue-level figures"""
        print("Generating residue-level docking figures...")
        print("  → Figure 3.1D: Lead coordination energetics")
        self.figure_31_lead_coordination()
        print("  → Figure 3.2C: 5-FU interaction heatmap")
        self.figure_32_5fu_residue_heatmap()
        print("  → Figure 3.3C: Allosteric pathway dynamics")
        self.figure_33_allosteric_pathway_dynamics()
        print("  → Figure 3.3D: Network connectivity")
        self.figure_33_network_connectivity()
        print("✅ All residue-level figures generated successfully!")


if __name__ == "__main__":
    generator = ResidueInteractionFigures("./figures")
    generator.generate_all_residue_figures()
    print(f"\nFigures saved to: {generator.output_dir.absolute()}")
