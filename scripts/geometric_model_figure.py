#!/usr/bin/env python3
"""Figure 1: geometric model of the study.

Four panels covering what the study actually establishes and where it does not:
  A  the coordination-geometry question (hemidirected PbS3 vs holodirected PbS4)
  B  Pb-S bond length: QM and EXAFS agree; the docked value does not
  C  the proposed allosteric pathway, drawn to relative scale
  D  the predicted vs measured affinity loss, showing the magnitude gap

Palette: Okabe-Ito-derived, validated with the dataviz six-checks validator
(all PASS in light mode, surface #fcfcfb).
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, Wedge

BLUE, GREEN, VERM, PURP = '#0072B2', '#009E73', '#D55E00', '#8C6BB1'
INK, MUTED, SURFACE = '#1a1a1a', '#6b6b6b', '#fcfcfb'
GRID = '#e0e0dd'

plt.rcParams.update({
    'font.size': 9, 'axes.labelsize': 9, 'axes.titlesize': 10,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
    'figure.dpi': 300, 'savefig.dpi': 300,
    'axes.edgecolor': GRID, 'axes.labelcolor': INK,
    'xtick.color': MUTED, 'ytick.color': MUTED, 'text.color': INK,
})

fig = plt.figure(figsize=(11, 8.2), facecolor=SURFACE)
gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.42,
                      left=0.055, right=0.965, top=0.885, bottom=0.075)

# ---------------------------------------------------------------- A
axA = fig.add_subplot(gs[0, 0]); axA.set_facecolor(SURFACE)
axA.set_xlim(0, 10); axA.set_ylim(0, 6); axA.axis('off')
axA.set_aspect('equal')
axA.set_title('A  Coordination geometry: the open question',
              loc='left', fontweight='bold', color=INK)

def donor(ax, x, y, label, color):
    ax.add_patch(Circle((x, y), 0.32, facecolor=color, edgecolor=SURFACE,
                        linewidth=2, zorder=3))
    ax.text(x, y, label, ha='center', va='center', color='white',
            fontsize=7.5, fontweight='bold', zorder=4)

# hemidirected PbS3
pbx, pby = 2.4, 3.5
axA.add_patch(Wedge((pbx, pby), 1.55, 30, 150, facecolor=VERM, alpha=0.14,
                    edgecolor='none', zorder=1))
axA.text(pbx, pby + 1.15, 'lone-pair\nvoid', ha='center', va='center',
         fontsize=7.5, color=VERM, style='italic', zorder=5)
for ang in (215, 270, 325):
    dx, dy = 1.25*np.cos(np.radians(ang)), 1.25*np.sin(np.radians(ang))
    axA.plot([pbx, pbx+dx], [pby, pby+dy], color=MUTED, lw=1.6, zorder=2)
    donor(axA, pbx+dx, pby+dy, 'S', GREEN)
axA.add_patch(Circle((pbx, pby), 0.44, facecolor=BLUE, edgecolor=SURFACE,
                     linewidth=2, zorder=3))
axA.text(pbx, pby, 'Pb', ha='center', va='center', color='white',
         fontsize=8.5, fontweight='bold', zorder=4)
axA.text(pbx, 1.15, 'hemidirected  PbS$_3$', ha='center', fontsize=8.5,
         fontweight='bold', color=INK)
axA.text(pbx, 0.6, 'literature consensus', ha='center', fontsize=7.5, color=MUTED)

# holodirected PbS4
pbx2 = 7.4
for ang in (55, 145, 235, 325):
    dx, dy = 1.25*np.cos(np.radians(ang)), 1.25*np.sin(np.radians(ang))
    axA.plot([pbx2, pbx2+dx], [pby, pby+dy], color=MUTED, lw=1.6, zorder=2)
    donor(axA, pbx2+dx, pby+dy, 'L', PURP)
axA.add_patch(Circle((pbx2, pby), 0.44, facecolor=BLUE, edgecolor=SURFACE,
                     linewidth=2, zorder=3))
axA.text(pbx2, pby, 'Pb', ha='center', va='center', color='white',
         fontsize=8.5, fontweight='bold', zorder=4)
axA.text(pbx2, 1.15, 'holodirected  4-coordinate', ha='center', fontsize=8.5,
         fontweight='bold', color=INK)
axA.text(pbx2, 0.6, 'Table 3.1.1, unconfirmed', ha='center', fontsize=7.5, color=MUTED)
axA.text(4.9, 3.5, 'vs', ha='center', va='center', fontsize=11,
         color=MUTED, style='italic')

# ---------------------------------------------------------------- B
axB = fig.add_subplot(gs[0, 1]); axB.set_facecolor(SURFACE)
axB.set_title('B  Pb–S distance: QM matches experiment, docking does not',
              loc='left', fontweight='bold', color=INK)
axB.axvspan(2.64, 2.68, color=GREEN, alpha=0.16, zorder=0)
axB.text(2.80, 1.20, 'EXAFS band\n2.64–2.68 Å', ha='center', va='center',
         fontsize=7.5, color=GREEN, fontweight='bold')
rows = [('EXAFS (protein PbS$_3$)', 2.66, 0.02, GREEN),
        ('QM, PBE0/def2-SVP\n+ ECP60MDF (this work)', 2.657, 0.0, BLUE),
        ('Docking, AutoDock Vina\n(Table 3.1.1)', 2.30, 0.20, VERM)]
for i, (lab, val, err, col) in enumerate(rows):
    y = 2 - i*0.8
    axB.errorbar(val, y, xerr=err if err else None, fmt='o', ms=11,
                 color=col, ecolor=col, elinewidth=2, capsize=5,
                 markeredgecolor=SURFACE, markeredgewidth=1.6, zorder=3)
    axB.text(val, y+0.2, f'{val:.3f} Å' if err == 0 else f'{val:.2f} ± {err:.2f} Å',
             ha='center', fontsize=8, fontweight='bold', color=INK)
    axB.text(2.08, y, lab, ha='right', va='center', fontsize=7.5, color=INK)
axB.annotate('', xy=(2.30, 0.15), xytext=(2.657, 0.15),
             arrowprops=dict(arrowstyle='<->', color=VERM, lw=1.6))
axB.text(2.48, 0.02, '0.35 Å short', ha='center', fontsize=8,
         color=VERM, fontweight='bold')
axB.set_xlim(1.55, 3.02); axB.set_ylim(-0.25, 2.45)
axB.set_xlabel('Pb–S distance (Å)')
axB.set_yticks([]); axB.grid(axis='x', color=GRID, lw=0.6)
for sp in ('left', 'right', 'top'):
    axB.spines[sp].set_visible(False)

# ---------------------------------------------------------------- C
axC = fig.add_subplot(gs[1, 0]); axC.set_facecolor(SURFACE)
axC.set_xlim(0, 10); axC.set_ylim(0, 5); axC.axis('off')
axC.set_aspect('equal')
axC.set_title('C  Proposed allosteric pathway (~30 Å)',
              loc='left', fontweight='bold', color=INK)
nodes = [('Cys-34', 0.7, BLUE, 'Pb$^{2+}$ site'),
         ('Lys-129', 2.9, MUTED, ''),
         ('Asp-183', 5.0, MUTED, ''),
         ('Trp-214', 7.1, GREEN, 'fluorescence\nreporter'),
         ('Lys-199', 9.2, PURP, '5-FU pocket')]
for i, (name, x, col, note) in enumerate(nodes):
    if i < len(nodes)-1:
        axC.add_patch(FancyArrowPatch((x+0.42, 2.9), (nodes[i+1][1]-0.42, 2.9),
                                      arrowstyle='-|>', mutation_scale=13,
                                      color=MUTED, lw=1.5, zorder=1))
    axC.add_patch(Circle((x, 2.9), 0.4, facecolor=col, edgecolor=SURFACE,
                         linewidth=2, zorder=3))
    axC.text(x, 2.15, name, ha='center', fontsize=8, fontweight='bold', color=INK)
    if note:
        axC.text(x, 1.50, note, ha='center', fontsize=6.8, color=MUTED)
axC.annotate('', xy=(0.7, 4.05), xytext=(9.2, 4.05),
             arrowprops=dict(arrowstyle='<->', color=INK, lw=1.2))
axC.text(4.95, 4.25, '~30 Å', ha='center', fontsize=8.5,
         fontweight='bold', color=INK)
axC.text(4.95, 0.35,
         'Pathway from the docked pose ensemble;\n'
         'residue identities not independently confirmed',
         ha='center', fontsize=6.8, color=MUTED, style='italic')

# ---------------------------------------------------------------- D
axD = fig.add_subplot(gs[1, 1]); axD.set_facecolor(SURFACE)
axD.set_title('D  Affinity loss: direction agrees, magnitude does not',
              loc='left', fontweight='bold', color=INK)
axD.barh([1.6], [20.9], height=0.42, color=VERM, zorder=3,
         edgecolor=SURFACE, linewidth=2)
axD.text(21.7, 1.6, '20.9×', va='center', fontsize=9,
         fontweight='bold', color=INK)
axD.barh([0.75], [5.4-1.9], left=[1.9], height=0.42, color=BLUE, zorder=3,
         edgecolor=SURFACE, linewidth=2)
axD.text(6.2, 0.75, '1.9–5.4×', va='center', fontsize=9,
         fontweight='bold', color=INK)
axD.set_yticks([1.6, 0.75])
axD.set_yticklabels(['Docking\nΔΔG = −1.8 kcal/mol',
                     'Measured (Stern-Volmer)\nΔΔG = 0.4–1.0 kcal/mol'],
                    fontsize=8, color=INK)
axD.tick_params(axis='y', length=0)
axD.annotate('', xy=(20.9, 1.28), xytext=(5.4, 1.28),
             arrowprops=dict(arrowstyle='<->', color=VERM, lw=1.6))
axD.text(13.0, 1.05, 'overestimated 4–11×', ha='center', fontsize=8,
         color=VERM, fontweight='bold')
axD.set_xlim(0, 26); axD.set_ylim(0.25, 2.15)
axD.set_xlabel('5-FU affinity reduction (fold)')
axD.grid(axis='x', color=GRID, lw=0.6)
for sp in ('right', 'top'):
    axD.spines[sp].set_visible(False)
axD.spines['left'].set_color(GRID)

fig.suptitle('Figure 1.  Geometric model of lead interference with 5-FU binding to HSA',
             x=0.07, y=0.965, ha='left', fontsize=12, fontweight='bold', color=INK)
fig.text(0.07, 0.925,
         'Blue: this work (quantum chemistry).  Green: independent experiment.  '
         'Vermillion: docking-derived, where it disagrees.',
         ha='left', fontsize=8, color=MUTED)

fig.savefig('Figure_1_Geometric_Model.png', facecolor=SURFACE,
            bbox_inches='tight')
print('wrote Figure_1_Geometric_Model.png')
