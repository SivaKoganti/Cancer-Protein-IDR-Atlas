# Cancer Protein IDR-LLPS Atlas

**A residue-resolved atlas of intrinsic disorder, phase separation, and virus-interface signals in cancer proteins**

> **CONFIDENTIAL** -- This repository is under embargo until journal publication. Unauthorized access, copying, or distribution is prohibited. Contact skoganti@ufl.edu for access.

---

## Overview

The Cancer Protein IDR-LLPS Atlas is a reproducible bioinformatics pipeline that constructs a residue-level atlas of intrinsic disorder regions (IDRs), liquid-liquid phase separation (LLPS) propensity, phylogenetic conservation, structural properties, clinical variant annotations, and oncovirus interaction mapping for **50 Cancer Gene Census (CGC) genes**.

The pipeline's central contribution is the **VIPP score** (Virus-Interaction + IDR + LLPS Prioritization) -- a weighted composite score that ranks each residue by its convergent disorder, phase-separation, and viral-interface evidence:

$$
\mathrm{VIPP}(i) = 0.35 \cdot I(i) + 0.35 \cdot L(i) + 0.30 \cdot V(i)
$$

where *I(i)* is the normalized disorder score, *L(i)* is the normalized LLPS proxy score, and *V(i)* is a binary virus-contact indicator. Weights are configurable in [`config.yaml`](config.yaml).

---

## Key Features

- **Residue-level multi-evidence atlas** -- Integrates disorder (IUPred2A), low-complexity (SEG), prion-like composition (PLAAC), LLPS propensity, AlphaFold pLDDT, phylogenetic conservation, ClinVar variants, SLiM/PTM motifs, and oncovirus interaction data into a single per-residue table per gene
- **VIPP prioritization score** -- Transparent, configurable composite score with ablation benchmarking against ClinVar pathogenic/benign labels
- **Delta-LLPS vulnerability** -- Simulates all 19 amino acid substitutions at each position to quantify phase-separation sensitivity to mutation
- **Conserved Disordered Regions (CDRs)** -- Identifies residues that are both conserved (>= 0.70) and disordered (>= 0.60), plus discordant regions where AlphaFold confidence conflicts with disorder predictions
- **SLiM and PTM motif scanning** -- 27 curated Short Linear Motif patterns from the ELM database covering phosphorylation, degradation, docking, nuclear signals, and SUMO modification
- **Oncovirus interaction mapping** -- Curated interaction regions for 7 oncoviruses (HPV16, HPV18, EBV, HBV, HCV, KSHV, MCPyV) mapped to per-residue flags
- **In-silico mutant pathway analysis** -- Generates IDR mutants, computes LLPS gain/loss, and maps effects across 10 oncogenic signaling pathways with PTM stratification
- **Drug repurposing screen** -- Screens NCATS-style compound databases against LLPS/virus-linked atlas residues
- **Interactive visualizations** -- D3.js chromosome karyotype atlas and landing-page dashboard with per-gene multi-track charts
- **Fully reproducible** -- Snakemake workflow, Docker container, and CI pipeline

---

## Gene Panel

50 high-confidence Cancer Gene Census genes:

| | | | | |
|---|---|---|---|---|
| TP53 | KRAS | EGFR | BRCA1 | BRCA2 |
| PIK3CA | BRAF | APC | PTEN | NRAS |
| CDKN2A | RB1 | ERBB2 | SMAD4 | STK11 |
| ATM | CHEK2 | MET | VHL | GNAS |
| IDH1 | IDH2 | JAK2 | KIT | PDGFRA |
| ALK | ROS1 | RET | NFE2L2 | PTCH1 |
| FBXW7 | CTNNB1 | MAP2K1 | MAP2K2 | NOTCH1 |
| NOTCH2 | ARID1A | ARID2 | SMARCA4 | KMT2A |
| KMT2D | TSC1 | TSC2 | MLH1 | MSH2 |
| MSH6 | PMS2 | ERCC2 | RNF43 | NTRK1 |

---

## Repository Structure

```
Cancer-Protein-IDR-Atlas/
|-- Snakefile                  # Workflow (17 rules)
|-- config.yaml                # Pipeline configuration (genes, weights, species)
|-- Dockerfile                 # Reproducible container (Ubuntu 22.04)
|-- requirements.txt           # Python dependencies
|-- scripts/                   # Pipeline scripts (27 modules)
|   |-- download_human_gene_fasta.py
|   |-- download_uniprot_orthologs.py
|   |-- download_alphafold_plddt.py
|   |-- run_iupred_batch.py
|   |-- run_seg_batch.py
|   |-- run_plaac_batch.py
|   |-- run_llps_proxy.py
|   |-- run_structure_proxy.py
|   |-- map_clinvar_to_proteins.py
|   |-- compute_phylo_conservation.py
|   |-- map_oncovirus_interactions.py
|   |-- build_disorder_atlas.py
|   |-- build_clinvar_view.py
|   |-- compute_delta_llps.py
|   |-- compute_cdr_novelty.py
|   |-- scan_slim_motifs.py
|   |-- analyze_idr_mutant_pathway_effects.py
|   |-- benchmark_vipp.py
|   |-- screen_ncats_repurposing.py
|   |-- build_idr_structural_library.py
|   |-- generate_llps_compound_library.py
|   |-- render_chromosome_atlas.py
|   |-- render_dashboard.py
|   |-- generate_publication.py
|   |-- summarize_case_studies.py
|   |-- summarize_vipp_hotspots.py
|   |-- validate_computational_repository.py
|-- tests/                     # Unit and integration tests (10 modules)
|-- .github/workflows/ci.yml   # CI smoke test
|-- data/                      # Downloaded data (gitignored, generated at runtime)
|-- results/                   # Pipeline outputs (gitignored, generated at runtime)
```

---

## Pipeline Architecture

The Snakemake workflow is organized into five stages:

### 1. Data Acquisition

| Script | Source | Output |
|--------|--------|--------|
| `download_human_gene_fasta.py` | UniProt REST API | `data/fasta/selected_proteins.fasta` |
| `download_uniprot_orthologs.py` | UniProt (6 species) | `data/orthologs/{gene}.fasta` |
| `download_alphafold_plddt.py` | EBI AlphaFold DB v4 | `results/alphafold/plddt_scores.tsv` |
| ClinVar download (wget) | NCBI FTP | `data/clinvar/variant_summary.txt.gz` |

### 2. Predictors

| Script | Method | Output |
|--------|--------|--------|
| `run_iupred_batch.py` | IUPred2A or amino-acid propensity heuristic (21-residue window) | `results/iupred/iupred_scores.tsv` |
| `run_seg_batch.py` | Shannon entropy low-complexity detection (12-residue window) | `results/seg/seg_regions.tsv` |
| `run_plaac_batch.py` | Prion-like Q/N composition (41-residue window) | `results/plaac/plaac_scores.tsv` |
| `run_llps_proxy.py` | Weighted LLPS propensity (Q/N + aromatic + charge + hydrophobic) | `results/llps/llps_scores.tsv` |
| `run_structure_proxy.py` | Structural and quantum-inspired biophysical scores | `results/structure/structure_scores.tsv` |

### 3. Mapping and Conservation

| Script | Description | Output |
|--------|-------------|--------|
| `map_clinvar_to_proteins.py` | Maps ClinVar variants to proteins by gene symbol and protein change | `results/clinvar/mapped_variants.tsv` |
| `compute_phylo_conservation.py` | Per-residue conservation from pairwise ortholog alignments (6 species) | `results/phylogeny/conservation.tsv` |
| `map_oncovirus_interactions.py` | Maps curated oncovirus-host interaction regions to residues | `results/oncovirus/oncovirus_per_residue.tsv` |

### 4. Novelty Layers

| Script | Description | Output |
|--------|-------------|--------|
| `compute_delta_llps.py` | Mutational vulnerability: mean absolute delta-LLPS across all 19 substitutions | `results/delta_llps/delta_llps_scores.tsv` |
| `compute_cdr_novelty.py` | CDR (conserved + disordered) and discordant region (pLDDT vs disorder) annotation | `results/cdr/cdr_summary.tsv` |
| `scan_slim_motifs.py` | 27 ELM-derived SLiM/PTM motif patterns with ClinVar overlap flags | `results/slim/slim_hits.tsv` + 4 more |
| `analyze_idr_mutant_pathway_effects.py` | In-silico IDR mutants mapped to 10 oncogenic pathways with PTM stratification | `results/mutants/` (4 files) |

### 5. Integration and Visualization

| Script | Description | Output |
|--------|-------------|--------|
| `build_disorder_atlas.py` | Merges all layers into per-gene atlas tables and global summary; computes VIPP | `results/atlas/{gene}_atlas.tsv` |
| `build_clinvar_view.py` | Adds ClinVar annotations to atlas tables | `results/atlas/{gene}_atlas_with_clinvar.tsv` |
| `benchmark_vipp.py` | Ablation benchmarking (AUROC, average precision) against ClinVar labels | `results/benchmark/model_metrics.tsv` |
| `render_chromosome_atlas.py` | Interactive D3.js chromosome karyotype with per-gene multi-track charts | `results/map/chromosome_atlas.html` |
| `render_dashboard.py` | Landing-page dashboard with summary statistics and gene table | `results/index.html` |
| `generate_publication.py` | Publication-quality matplotlib/seaborn figures (300 DPI) | `results/figures/` |

---

## Quickstart

### Prerequisites

- Python 3.10+
- Dependencies listed in [`requirements.txt`](requirements.txt): biopython, pandas, numpy, matplotlib, seaborn, pyfaidx, tqdm, pysam, requests, pyyaml

### Installation

```bash
git clone https://github.com/SivaKoganti/Cancer-Protein-IDR-Atlas.git
cd Cancer-Protein-IDR-Atlas
python -m venv .venv
source .venv/bin/activate   # On Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

### Run the Full Pipeline

```bash
snakemake --cores 4
```

Or with Docker:

```bash
docker build -t idr-atlas .
docker run --rm idr-atlas
```

### Run Tests

```bash
pytest -q tests/
```

### Validate Outputs

```bash
python scripts/validate_computational_repository.py
```

---

## Configuration

All pipeline parameters are centralized in [`config.yaml`](config.yaml):

| Section | Parameters |
|---------|------------|
| `genes` | List of 50 CGC gene symbols |
| `llps` | Window size (41), Q/N weight (0.45), aromatic weight (0.25), charge pattern weight (0.15), hydrophobic weight (0.15) |
| `structural.quantum_proxy` | Dielectric constant (4.0), polarizability factor (1.0), stacking coefficient (1.2), charge penalty (0.6) |
| `phylogeny` | 6 species (human, mouse, rat, zebrafish, chicken, fruitfly), source: UniProt |
| `oncovirus` | 7 viruses: HPV16, HPV18, EBV, HBV, HCV, KSHV, MCPyV |
| `vipp` | IDR weight (0.35), LLPS weight (0.35), virus weight (0.30) |

---

## Oncogenic Pathways Analyzed

The mutant pathway analysis maps genes to 10 canonical cancer signaling pathways:

1. TP53 / Cell Cycle
2. RTK / RAS / MAPK
3. PI3K / AKT / mTOR
4. WNT / Beta-Catenin
5. Notch / Hedgehog
6. DNA Damage Repair
7. Chromatin Remodeling
8. JAK / STAT / Cytokine
9. TGF-Beta / SMAD
10. Metabolic / Hypoxia / Redox

---

## Core Outputs

| Output | Path | Description |
|--------|------|-------------|
| Per-gene atlas | `results/atlas/{gene}_atlas.tsv` | Residue-level disorder, LLPS, conservation, structure, VIPP, virus, SLiM, PTM, CDR, delta-LLPS |
| Global atlas | `results/atlas/global_disorder_phylogeny_atlas.tsv` | Summary statistics across all 50 genes |
| ClinVar-enriched atlas | `results/atlas/{gene}_atlas_with_clinvar.tsv` | Atlas with variant count and clinical significance |
| VIPP benchmark | `results/benchmark/model_metrics.tsv` | AUROC and average precision for ablated models |
| Ablation matrix | `results/benchmark/ablation_matrix.tsv` | Feature contribution analysis |
| Compound library | `results/llps/llps_compound_library.tsv` | LLPS-stratified small-molecule candidates |
| Repurposing candidates | `results/llps/ncats_repurposing_candidates.tsv` | NCATS-screened compound hits |
| Chromosome atlas | `results/map/chromosome_atlas.html` | Interactive D3.js karyotype visualization |
| Dashboard | `results/index.html` | Summary landing page |
| Publication figures | `results/figures/` | Matplotlib/seaborn figures at 300 DPI |
| Pathway analysis | `results/mutants/pathway_differential_llps.tsv` | Pathway-level LLPS differential statistics |
| CDR summary | `results/cdr/cdr_summary.tsv` | Conserved disordered and discordant regions |
| SLiM hits | `results/slim/slim_hits.tsv` | Short linear motif matches with ClinVar overlaps |
| VIPP hotspots | `results/vipp/top_vipp_hotspots.tsv` | Top-ranked VIPP residues across all genes |

---

## Tests

10 test modules covering core pipeline components:

| Test | Covers |
|------|--------|
| `test_build_disorder_atlas.py` | Atlas construction and VIPP computation |
| `test_clinvar_integration.py` | ClinVar variant mapping |
| `test_analyze_idr_mutant_pathway_effects.py` | Mutant generation and pathway analysis |
| `test_benchmark_vipp.py` | VIPP ablation benchmarking |
| `test_generate_llps_compound_library.py` | Compound library generation |
| `test_generate_publication.py` | Figure generation |
| `test_idr_structural_library.py` | Structural library construction |
| `test_screen_ncats_repurposing.py` | NCATS repurposing screen |
| `test_summarize_case_studies.py` | Case study extraction |
| `test_validate_computational_repository.py` | Output file validation |

---

## CI/CD

GitHub Actions runs a smoke test on every push and pull request to `main`:

- **Runner:** `ubuntu-latest`, Python 3.10
- **Test:** `scripts/ci_smoke.py`

---

## Docker

The Dockerfile provides a fully reproducible execution environment:

- **Base:** Ubuntu 22.04
- **Includes:** Python 3, build-essential, wget, git, JRE, Perl
- **Default command:** `snakemake --cores 4`

```bash
docker build -t idr-atlas .
docker run --rm idr-atlas
```

---

## Documentation

| Document | Description |
|----------|-------------|
| [MANUSCRIPT.md](MANUSCRIPT.md) | Full peer-review-ready manuscript |
| [SUPPLEMENTARY_INFORMATION.md](SUPPLEMENTARY_INFORMATION.md) | Extended methods, validation metrics, per-gene analysis |
| [INDEX.md](INDEX.md) | Publication package table of contents |
| [PUBLICATION_README.md](PUBLICATION_README.md) | Publication package quick links |
| [CONTRIBUTING.md](CONTRIBUTING.md) | Contribution guidelines (post-publication) |
| [CONTRIBUTORS.md](CONTRIBUTORS.md) | Author and acknowledgment details |
| [CONFIDENTIALITY.md](CONFIDENTIALITY.md) | Confidentiality and data use agreement |
| [CONFIDENTIALITY_TRANSITION.md](CONFIDENTIALITY_TRANSITION.md) | Plan for public release after publication |

---

## Licensing and Intellectual Property

**Copyright 2026, Siva Koganti and Co-Authors. All rights reserved.**

This repository is currently **proprietary and confidential**. No license to use, copy, modify, or distribute is granted. The project may be released under an open-source license (MIT, CC-BY-4.0, or similar) after publication.

See [LICENSE](LICENSE) and [CONFIDENTIALITY.md](CONFIDENTIALITY.md) for full terms.

---

## Acknowledgments

This project builds on the following tools and databases:

- [BioPython](https://biopython.org/) -- sequence analysis and alignment
- [IUPred2A](https://iupred2a.elte.hu/) -- intrinsic disorder prediction
- [PLAAC](http://plaac.wi.mit.edu/) -- prion-like amino acid composition
- [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/) -- predicted structure confidence (pLDDT)
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) -- clinical variant annotations
- [ELM](http://elm.eu.org/) -- eukaryotic linear motif patterns
- [UniProt](https://www.uniprot.org/) -- protein sequences and ortholog data
- [Snakemake](https://snakemake.readthedocs.io/) -- workflow management
- [Pandas](https://pandas.pydata.org/), [NumPy](https://numpy.org/) -- data processing
- [Matplotlib](https://matplotlib.org/), [Seaborn](https://seaborn.pydata.org/) -- publication figures
- [D3.js](https://d3js.org/) -- interactive chromosome visualization

---

**Contact:** skoganti@ufl.edu
