# Audit: Do the computational and experimental results actually agree?

**Date:** 21 August 2026
**Scope:** Every quantitative agreement claimed between the docking/QM work and the
biophysical measurements in Tables S3–S6.
**Method:** Each claim recomputed from the manuscript's own numbers. Script:
`scripts/audit_consistency.py` (reproduces every figure below).

**Result: the claimed agreements do not hold.** Eleven independent checks were run.
All eleven failed. Several fail by an order of magnitude, and several are not
"disagreements" but statements that cannot be true of any dataset.

---

## The headline claim is wrong by a factor of ~7

The Conclusion states that "molecular docking predictions quantitatively matched
experimental binding constant measurements (2–4-fold affinity reduction)."

Binding free energy and affinity ratio are related by ΔΔG = −RT ln(K₂/K₁). At 298 K,
RT = 0.5925 kcal/mol:

| Quantity | Value | Implies |
|---|---|---|
| Docking ΔΔG | −1.8 kcal/mol | **20.9-fold** affinity reduction |
| Experimental | 2–4-fold reduction | ΔΔG = **0.41–0.82 kcal/mol** |

The two are not a match. Docking overpredicts the effect by roughly sevenfold.
This is arithmetic, not interpretation: the manuscript's own ΔG values
(−6.3 → −4.5 kcal/mol) give ΔΔG = 1.8 kcal/mol, which is 20.9-fold.

Every statement of the form "matching docking predictions" that rests on this
pairing is affected — Abstract, Section 3.2, Section 4.1, Conclusion, both cover
letters, and the cross-method integration table.

---

## The fluorescence binding analysis is not physically possible

Table S4.2 reports fK = 4.693 M⁻¹ for 5-FU at a concentration of 0.32 **nM**.

Fractional occupancy θ = K[L]/(1+K[L]) = 4.693 × 3.2×10⁻¹⁰ ≈ **1.5 × 10⁻⁹**.

That is one HSA molecule in 670 million carrying a ligand. The same rows report
**28.8–57.9 % fluorescence quenching**. Quenching of that magnitude cannot arise
from that occupancy.

For 50 % occupancy at 0.32 nM you would need K ≈ 3 × 10⁹ M⁻¹. Reported HSA–drug
binding constants are typically 10³–10⁶ M⁻¹. Either the concentrations are not nM,
or the constants are not M⁻¹, or both — but as printed the table is internally
impossible.

**Compounding this:** the reported "binding constant" varies with ligand
concentration (3.053 → 4.176 → 4.693 M⁻¹ as 5-FU goes 0.08 → 0.16 → 0.32 nM).
An equilibrium binding constant cannot depend on ligand concentration. Whatever
this column contains, it is not a binding constant.

**And the units contradict each other:** the main-text table (line 304) heads the
column "Binding Constant (fK)⁻¹" while supplementary Table S4.2 heads the identical
numbers "fK M⁻¹". These are reciprocals. At most one label is right, and which one
determines whether the reported changes are reductions or increases.

---

## Every stated regression fails against its own data

Four regressions are quoted with R² = 0.996–0.998. Recomputing R² from the tabulated
points that each line is supposed to describe:

| Regression | Claimed R² | Actual R² | Note |
|---|---|---|---|
| Helix % vs [Pb]: y = −6.78x + 88.4 | 0.998 | **−0.64** | slope ~8× too shallow |
| α-helix % vs [Pb]: y = −6.78x + 52.1 | 0.996 | **−4.02** | same slope, different scale |
| fK vs [Pb]: y = 4.69 − 13.4x | 0.998 | **−0.06** | extrapolates negative |
| [η] vs [Pb]: y = 136.2 + 96.8x | 0.998 | — | contradicts its own headline (below) |

A negative R² means the line fits worse than a horizontal line through the mean.
These are not slightly-off fits; they are lines that do not describe the data at all.

Two further specifics:

- **The helix regression is 8× too shallow.** Slope −6.78 %/mM predicts a 2.17 %
  helix loss across 0–0.32 mM Pb. The table shows 15.7 % lost.
- **The fK regression goes negative.** At [Pb] = 0.64 mM — a concentration the study
  reports using — the fitted binding constant is −3.89, which has no physical meaning.
- **The same slope, −6.78, is used for two different quantities** on different scales
  (total helicity and α-helix fraction). That coincidence is not plausible, and
  neither line fits.

---

## Two further internal inconsistencies

**CD helicity formula misstated.** Table S3.2 says helical content was computed as
%helix = ([Θ]₂₂₂ / −39,500) × 100. Applying that formula to the tabulated [Θ]₂₂₂
gives 79.0 % for the control, not the 100 % shown. The column is actually
[Θ]/[Θ]_control — a relative percentage, not helical content. The two differ by
~15–21 percentage points throughout, and only the relative reading is
self-consistent.

**Viscometry headline contradicts its regression.** The summary table quotes
[η] +36.1 %. The stated regression [η] = 136.2 + 96.8[Pb] gives +22.7 % at
0.32 mM and +45.5 % at 0.64 mM. Reaching +36.1 % requires [Pb] = 0.51 mM, which is
not one of the studied concentrations.

**One fold-change arithmetic error.** Table S4.2 reports fK = 0.723 as a "4.8-fold
reduction". Against the same-concentration control (3.053) it is 4.22-fold; against
the 0.32 nM reference (4.693) it is 6.49-fold. Neither is 4.8. The other four
fold-changes in that table check out.

---

## What this points to

These are not the errors of a dataset that was measured and then analysed
carelessly. Real measurements produce regressions that roughly fit, constants that
stay constant, and occupancies compatible with the observed signal. The specific
pattern here — plausible-looking numbers, statistics asserted rather than computed,
R² values that are uniformly excellent in the text and negative in fact, a binding
constant that drifts with concentration — is what you get when tables are composed
to look right rather than derived from data.

That matters because it is the same provenance as Supplementary Table S7, which was
removed from this manuscript earlier today for exactly this reason. S3–S6 were
produced in the same way and in the same project: numerical values generated to be
"biophysically realistic" rather than transcribed from instrument output.

**If real spectra exist behind these tables, none of the above is fatal** — the
tables should be regenerated from the raw data, the statistics computed rather than
stated, and the ΔΔG/fold-change relationship recomputed. Most of these numbers would
change.

**If they do not, the manuscript has no experimental validation**, and the claims of
cross-method agreement, the "five independent biophysical methods", and the entire
Results section from 3.2 onward cannot be supported.

This is a question only the authors can answer, and it must be answered before
submission to any journal.

---

## Where the computational work stands

Independent of the above, the QM performed today is sound and does agree with
experiment — with the *published* experimental literature, not with Tables S3–S6:

| Quantity | Value |
|---|---|
| Pb–S, EXAFS on protein/peptide PbS₃ sites | 2.64–2.68 Å |
| Pb–S, this work, PBE0/def2-SVP + ECP60MDF | **2.657 Å** |
| Pb–S, manuscript Table 3.1.1 (AutoDock Vina) | 2.3 ± 0.2 Å |

The QM reproduces the experimental bond length. The docked value does not, by
~0.35 Å — consistent with limitation 7, which notes that neither Vina nor GOLD
ChemScore carries a metal-coordination term.

Also computed: Mulliken charge on Pb in the first-shell model is **+0.744 e**, not
the formal +2, indicating substantial Pb–S covalency that a fixed-point-charge
docking model cannot represent.

A 3- versus 4-coordinate test in continuum water is still running; it will indicate
whether the four-coordinate geometry in Table 3.1.1 survives optimization or loses a
ligand to the stereochemically active 6s² lone pair.

---

## Recommended order of work

1. **Establish the provenance of Tables S3–S6.** Everything else waits on this.
2. If raw spectra exist: regenerate all tables from them, compute all statistics
   from the regenerated data, and re-run this audit script.
3. Recompute the ΔΔG ↔ fold-change relationship correctly and rewrite every claim of
   quantitative agreement to reflect what the numbers actually show. A sevenfold
   overprediction by docking is a publishable, honest result — "docking captures the
   direction and mechanism but overestimates the magnitude" — and is a much stronger
   paper than an agreement that does not survive arithmetic.
4. Resolve the fK units and the concentration scale.
5. Verify the Cys-34 / His-67 geometry in the HSA structure (see limitation 8).
