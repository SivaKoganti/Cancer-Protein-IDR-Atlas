# Anticipated Reviewer Comments and Prepared Responses

**Manuscript:** Molecular Docking, Quantum-Chemical and Biophysical Studies of Lead
Interference with Human Serum Albumin and 5-Fluorouracil Binding
**Prepared:** 21 August 2026
**Purpose:** Internal preparation. Not for submission. Reviewer comments are
paraphrased in the form they are most likely to arrive.

---

## How to use this document

Comments are grouped by how dangerous they are to the manuscript, not by section.
The first group contains the objections that could result in rejection; they are
the ones to rehearse. Group 2 will likely appear as required revisions. Group 3 are
routine.

Two of these — R1 and R2 — **cannot be answered from the current draft**. They
depend on data only the authors hold. Prepare those answers before submitting, not
after the review arrives.

---

# GROUP 1 — Potentially fatal. Answer these before submission.

## R1. "The fluorescence binding constants are not physically possible as reported."

> *Likely wording:* "The authors report binding constants of ~4.7 M⁻¹ with 5-FU at
> 0.32 nM. At these values fewer than one albumin molecule in 10⁸ would carry a
> ligand, yet 12–58% quenching is reported. Either the concentrations or the
> constants are wrong. I cannot evaluate the central claim of the paper until this
> is resolved."

**This is the most likely cause of rejection.** A referee who checks one number will
check this one.

**Prepared response:**

> We thank the reviewer for identifying this. The reviewer is correct that the
> combination as tabulated is not physically possible, and we addressed it in
> Section 4.5, limitation 10, of the submitted manuscript.
>
> [IF THE RAW DATA CONFIRM A SCALE ERROR:] On re-examining the original
> Stern-Volmer plots we confirm that the 5-FU concentrations are [micromolar /
> millimolar] and that the tabulated constants carry a factor of [10^n]. Corrected
> values are given in the revised Table S4.2, and the raw fluorescence spectra and
> Stern-Volmer plots are now deposited at [DOI]. The binding constant for 5-FU in
> native HSA becomes [value] M⁻¹, consistent with the weak albumin binding reported
> for fluoropyrimidines. No conclusion changes, because every claim we draw from
> these data is a ratio between conditions measured in the same session, and a
> uniform scale error cancels in that ratio.

**What must be done first:** recover the original Stern-Volmer plots and establish
the true concentration axis and units. If they cannot be recovered, the fluorescence
section must be repeated or removed. **Do not submit without settling this.**

---

## R2. "What is the provenance of the data in Tables S3–S6?"

> *Likely wording:* "The supplementary tables report coefficients of determination
> of 0.996–0.998 for dose-response fits. When I refit the tabulated values I obtain
> substantially lower values. I would like to see the raw spectra."

**Prepared response:**

> The reviewer is right, and we corrected this before submission: all four
> regressions were refitted by least squares from the tabulated data for the current
> version, giving R² = 0.54–0.87, and the manuscript now states that a linear model
> is inappropriate for a saturating dose-response (Section 4.5, limitation 11). The
> earlier coefficients were not obtained from these data and have been withdrawn.
> The raw spectra are deposited at [DOI], and `scripts/audit_consistency.py`
> recomputes every relationship we claim.

**What must be done first:** deposit the raw data. A referee who asks this question
and receives no data will recommend rejection. If the underlying spectra do not
exist in recoverable form, the affected sections cannot be submitted.

---

## R3. "The docking overestimates the effect by an order of magnitude. Why should I believe the rest of the computational work?"

> *Likely wording:* "The authors' own analysis shows the docking ΔΔG corresponds to
> a 21-fold affinity change against a measured 2–5-fold. This is a substantial
> failure. On what basis do they retain the allosteric pathway, which comes from the
> same calculation?"

**This is the strongest scientific objection**, and the answer is genuinely
defensible — but only if given precisely.

**Prepared response:**

> The reviewer raises the right question, and we agree that a scoring function which
> fails by an order of magnitude on energetics cannot be trusted for energetics. Our
> claim is narrower. The Vina scoring function is a five-term empirical expression
> with no electrostatic, desolvation or metal-coordination term [Trott & Olson 2010],
> and benchmarking on metalloproteins has found that docking programs pose such
> complexes acceptably while failing to rank affinities [Chen et al. 2019]. That is
> precisely the pattern we observe: the pose and the site are informative, the energy
> is not. We therefore use the docking only to generate a hypothesis about *where*
> a coupling might run, and we test that hypothesis experimentally — Trp-214 lies on
> the proposed path and responds to lead, and four further methods report structural
> change over the same concentration range. The magnitude of the effect is taken
> entirely from the measurements. We have rewritten the Discussion (Section 4.1) to
> make this division explicit rather than presenting the two as mutually confirming.

---

## R4. "Cys-34 as the primary lead site contradicts the only study that measured it."

> *Likely wording:* "Belatik et al. (PLoS ONE 2012) localized Pb(II) on albumin to
> nitrogen and oxygen donors by XPS, not to the Cys-34 thiol. The authors assert a
> thiolate site without engaging this result."

**Prepared response:**

> We address this in Section 4.5, limitation 9. We agree the assignment is not
> settled, and we do not claim otherwise in the revised manuscript: our site
> assignment derives from docking with the scoring-function limitations described
> above, and Belatik et al. remains the only study to have probed the question
> spectroscopically. We note also the precedent that the strong Cd(II) sites on
> albumin do not involve Cys-34 [Sadler & Viles 1996]. We have restated the Cys-34
> assignment as a hypothesis requiring independent structural confirmation — EXAFS at
> the Pb L₃ edge, ²⁰⁷Pb NMR, or crystallography — rather than as an established
> result. Importantly, the allosteric argument does not depend on which remote site
> lead occupies, only that it is remote from the drug pocket.

---

## R5. "Can Cys-34 and His-67 coordinate the same metal ion?"

> *Likely wording:* "His-67 is a ligand of albumin's interdomain site A. Cys-34 is
> in subdomain IA. Table 3.1.1 places both around one Pb(II). What is the
> Sγ···Nε distance in the structure?"

**A referee familiar with albumin will ask this.** It is a factual question with a
factual answer, and the answer is not currently in the manuscript.

**Prepared response:**

> [MEASURE THIS BEFORE SUBMITTING.] The Sγ(Cys-34)···Nε(His-67) distance in [PDB
> ID] is [X] Å. [IF COMPATIBLE:] This is consistent with joint coordination of a
> single Pb(II), whose ionic radius and typical Pb–S/Pb–N distances are [ ]. [IF
> NOT:] The reviewer is correct that these residues cannot coordinate one metal ion.
> We have restricted Table 3.1.1 to the Cys-34 contact and reassigned the remaining
> residues to a separate site, which does not affect the allosteric argument.

**Action:** open the structure in PyMOL and measure it. One command:
`dist d1, /1ao6//A/CYS`34/SG, /1ao6//A/HIS`67/NE2`

---

# GROUP 2 — Expect these as required revisions.

## R6. "Why cluster models rather than QM/MM?"

**Prepared response:**

> Cluster models were chosen deliberately for the question at hand — the geometry
> and electronic structure of the first coordination shell — which is well posed at
> this level and for which the relevant experimental comparison (EXAFS Pb–S
> distances) is also a first-shell quantity. QM/MM would be required to treat the
> allosteric pathway itself, and we agree this is the natural next step; we note that
> converged QM/MM for a metal site may require QM regions of several hundred atoms
> [Kulik et al. 2016] together with conformational sampling [Mehmood & Kulik 2020],
> which places it beyond the scope of the present work. We have stated this in
> Section 4.6.

## R7. "The coordination-number question is left unresolved."

**Prepared response:**

> Correct, and we report it as unresolved rather than overstating it. Our
> Pb(SCH₃)₄²⁻ calculation in continuum water did not converge to a defined
> geometry: all four Pb–S distances lengthened together while the energy flattened,
> the behaviour of a marginally bound dianion, and continuum solvation alone does not
> stabilize a −2 thiolate complex. We report this as inconclusive in Section 3.1.2
> rather than presenting the endpoint as a result. Resolving it properly requires
> explicit first-shell waters or counterions. The published literature favours
> three-coordinate hemidirected PbS₃ [Shimoni-Livny 1998; Magyar 2005; Gourlaouen &
> Parisel 2007], and we defer to it while noting that our own calculations do not
> independently establish the point.

## R8. "Lead concentrations are far above physiological."

**Prepared response:**

> We agree, and state this in Section 4.5, limitation 2. Concentrations of
> 0.032–0.64 mM exceed blood lead levels in even heavily exposed workers (1–2 µM) by
> two to three orders of magnitude. They were chosen so that structural changes would
> be resolvable by the methods used. The mechanistic inference — that lead binding at
> a remote site perturbs drug binding — does not require the concentrations to be
> physiological, but the clinical extrapolation does, and we have restricted our
> clinical statements accordingly. Establishing whether the effect operates at
> occupational exposure levels requires measurements at those concentrations, which
> we identify as the priority next experiment.

## R9. "The allosteric pathway residues are not experimentally validated."

**Prepared response:**

> Correct. The pathway is a computational proposal supported by one direct
> observation — Trp-214 lies on the predicted path and responds to lead — and by
> aggregate structural changes consistent with, but not diagnostic of, that
> particular route. We say so in Section 3.3 and in the Conclusion. Site-directed
> mutagenesis of Lys-129 and Asp-183, with the fluorescence and CD measurements
> repeated on the mutants, is the decisive experiment and is identified as such in
> Section 4.6.

## R10. "The DSC result contradicts the CD and FTIR results."

**Prepared response:**

> The two are reconcilable and we now explain the reconciliation in Section 4.1.
> Unfolding enthalpy rises 2.8-fold while secondary structure content falls. A metal
> ion bridging carboxylate side chains introduces electrostatic cross-links that must
> be broken on unfolding, raising the enthalpic cost, while locally disordering the
> helices it perturbs. Thermodynamic stability and secondary-structure content are
> distinct quantities, and lead appears to raise one while lowering the other.

---

# GROUP 3 — Routine.

## R11. Generative-AI use
**Response:** Disclosed in full in Section 6.1, following ICMJE and ACS policy,
with the specific tasks itemized. The assistant does not meet authorship criteria
and the authors retain responsibility for all content.

## R12. Only one docking program is described in detail
**Response:** Two were used, AutoDock Vina and GOLD ChemScore, with comparable
results (Section 3.1). Given the scoring-function limitations we now document, we
would not argue that agreement between two empirical functions constitutes
validation — both share the same structural deficiency at metal centres.

## R13. Statistical treatment
**Response:** All fits were recomputed by least squares for this version;
`scripts/audit_consistency.py` reproduces each one. Where a linear model fits
poorly we say so and give the log-linear alternative rather than quoting the linear
fit alone.

## R14. Figure quality
**Response:** All figures are 300 dpi. Figure 1 uses a colourblind-safe palette
validated for CVD separation and contrast against the page surface.

---

# Pre-submission checklist

Derived from the above. The first two are blocking.

- [ ] **R1** — Resolve the fK units and concentration scale against the original
      Stern-Volmer plots. **Blocking.**
- [ ] **R2** — Deposit the raw spectroscopic data and insert the DOI in Section 6.3.
      **Blocking.**
- [ ] **R5** — Measure the Cys-34/His-67 distance; adjust Table 3.1.1 if required.
- [ ] Complete author names, affiliations, and the corresponding author's details.
- [ ] Complete funding sources and grant numbers (Section 6, Section 6.1).
- [ ] Complete the CRediT author contributions (Section 6.2).
- [ ] Declare conflicts of interest.
- [ ] Confirm the PDB accession used for docking is stated in Section 2.1.
- [ ] Run `python3 scripts/audit_consistency.py` and confirm 0 failures.
- [ ] Read Section 4.5 in full — it is long, and it is what protects the manuscript
      under review.

---

## A note on strategy

The manuscript's strongest feature is not the mechanism; it is that the
computational predictions were checked against experiment and against the
coordination-chemistry literature, and that the failures are reported rather than
smoothed over. Two independent lines — the quantum-chemical Pb–S distance matching
EXAFS, and the docking ΔΔG failing by an order of magnitude — together give a clear
account of what empirical docking can and cannot establish at a metal site.

That is worth foregrounding in the cover letter. Referees are used to manuscripts in
which every method agrees with every other; a paper that says where its methods
disagree, and why, tends to be read as more careful rather than less convincing.

The corresponding risk is that the honesty in Section 4.5 will be read as weakness
if the two blocking items above are unresolved. A manuscript that documents a
problem it has fixed is strong. One that documents a problem it has not fixed
invites the reviewer to fix it by rejection.
