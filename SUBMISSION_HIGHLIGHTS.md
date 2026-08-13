# Submission Highlights

## Title

A residue-resolved atlas of disorder, phase separation, and virus-interface signals in cancer proteins

## One-Sentence Summary

We present a computational atlas that integrates intrinsic disorder, LLPS propensity, structural-confidence, conservation, and oncovirus-host interaction evidence to prioritize cancer-protein residues potentially sensitive to viral rewiring and regulatory disruption, with benchmark interpretation explicitly limited to internal evaluation.

## Significance Statement

Most cancer-protein resources remain centered on sequence variation or folded-domain interpretation and do not explicitly model how disordered regulatory sequence, condensate propensity, and viral interface burden intersect at residue resolution. This study addresses that gap by introducing a virus-aware residue-level atlas and a unified prioritization framework that can rank candidate hotspots across 50 canonical cancer genes. The resulting platform is designed for fully computational hypothesis generation, comparative analysis, and future benchmark-driven refinement.

## Key Advances

1. Integrates intrinsic disorder, LLPS propensity, conservation, structural-confidence, and oncovirus interaction evidence in one residue-level framework.
2. Introduces VIPP, a composite residue prioritization score for identifying virus-sensitive regulatory hotspots.
3. Adds PTM-stratified mutant differential analysis to quantify pathway-specific LLPS sensitivity in PTM-context versus non-PTM-context residues.
4. Adds viable IDR-mutant LLP-phase differential analysis to quantify phase-conditioned mutant impact and phase-specific PTM/pathway coupling.
5. Activates novelty layers beyond standard atlas summaries, including discordant pLDDT-disorder regions and conserved disordered regions.
6. Produces publication-ready figures, ranked hotspot tables, and interactive atlas outputs from a single computational workflow.

## Current Quantitative Snapshot

1. 50 cancer genes analyzed.
2. 64,903 residue records in the virus-aware mapping layer.
3. 873 virus-interacting residues across 7 genes.
4. Internal benchmark baselines were near chance (IDR-only AUROC 0.502534; IDR+LLPS AUROC 0.488739).
5. Virus-aware and full formulations showed strong separation in the internal benchmark setting (IDR+LLPS+virus AUROC 0.999992; full VIPP AUROC 0.999999).
6. Top-10 enrichment in the full VIPP configuration reached 10.000462 in the internal benchmark.
7. PTM-context pathway differential outputs are included in the manuscript package (Figure 8).
8. Viable-mutant LLP-phase differential outputs are included in the manuscript package (Figure 9).

## Key Caveat

The current build is computationally mature but not yet biologically final: in the absence of parsed predictor-backed IUPred output, disorder scores fall back to a deterministic heuristic, and benchmark positives are internally curated. Accordingly, benchmark values should be interpreted as controlled internal ranking behavior rather than external predictive generalization.

## Suggested Cover-Letter Framing

This manuscript introduces a computationally novel atlas framework for cancer proteins in which host-virus interface evidence is modeled together with intrinsic disorder and phase-separation-related sequence behavior at residue resolution. The study is particularly distinctive in its explicit residue-level integration of viral interaction burden with regulatory biophysics, yielding a platform that can generate mechanistically focused hypotheses not accessible from static structure-centric annotations alone.

## Suggested Highlights

1. Virus-aware residue-level atlas across 50 Cancer Gene Census genes.
2. Unified disorder, LLPS, conservation, pLDDT, and host-virus interaction framework.
3. VIPP prioritization highlights candidate virus-sensitive regulatory hotspots.
4. Active computational novelty layers include discordant and conserved disordered regions.
5. End-to-end reproducible workflow visualization supports method transparency and re-execution (Figure 5).
6. Entrez chromosome-locus visualization links genomic position with locus-level IDR and LLPS context (Figure 6).
7. PTM-context pathway differentials provide direction-aware LLPS vulnerability contrasts across oncogenic pathways (Figure 8).
8. Viable IDR-mutant LLP-phase analysis shows phase-specific LLPS effect patterns and downstream PTM/pathway coupling (Figure 9).