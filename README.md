# ⚠️ CONFIDENTIAL - Cancer Protein IDR & LLMS Atlas

**⚠️ CONFIDENTIAL MATERIAL - NOT FOR PUBLIC DISTRIBUTION**

This repository is currently under confidential development and embargo. Unauthorized access, copying, or distribution is strictly prohibited.

**Authorized Personnel Only** - Contact s.koganti@[institution].edu for access

---

## Project Description (Internal)

This repository contains a computationally validated pipeline for intrinsic disorder region (IDR), low-complexity region (LCR), and liquid–liquid phase separation (LLPS) analyses in cancer-associated human proteins. It integrates ClinVar evidence, phylogenetic conservation, oncovirus interaction mapping, VIPP prioritization, a sequence-only spin-glass-inspired heuristic layer, and repurposing screening into a reproducible computational framework.

The current build is validated through generated artifacts and a repository-level validation report in [results/validation/computational_repository_validation.tsv](results/validation/computational_repository_validation.tsv).

### Core computational outputs

- Atlas tables per gene in [results/atlas](results/atlas)
- Spin-glass-inspired residue heuristics in [results/spin_glass/spin_glass_scores.tsv](results/spin_glass/spin_glass_scores.tsv)
- LLPS compound library in [results/llps/llps_compound_library.tsv](results/llps/llps_compound_library.tsv)
- VIPP benchmark metrics in [results/benchmark/model_metrics.tsv](results/benchmark/model_metrics.tsv) with virus-aware models showing near-perfect ranking performance on the current benchmark set
- Repurposing candidates in [results/llps/ncats_repurposing_candidates.tsv](results/llps/ncats_repurposing_candidates.tsv)
- Repository validation report in [results/validation/computational_repository_validation.tsv](results/validation/computational_repository_validation.tsv)

### VIPP scoring model

The repository now computes a formal, configurable VIPP score for each residue as a weighted combination of normalized evidence channels:

### Novelty relative to existing pipelines

This project is not merely a wrapper around standard disorder or LLPS predictors. It is a cancer-gene-specific atlas that combines residue-level disorder, LLPS, conservation, structural proxies, ClinVar context, curated oncovirus interaction evidence, and a transparent spin-glass-inspired comparative score in one reproducible framework, then ranks residues with a transparent VIPP score and benchmarks the contribution of each evidence layer through ablation analyses.

### Spin-glass-inspired analysis layer

The spin-glass layer is a deterministic, sequence-only comparative analysis. For each residue, the pipeline encodes local amino-acid biophysical categories, computes local couplings over a configurable window, and reports interpretable heuristic features including:

- `local_energy`
- `local_frustration`
- `coupling_variance`
- `susceptibility_like`
- `spin_glass_score`
- shuffled-window null expectations for the main features

The resulting `spin_glass_score` is intended as a residue-level comparative surrogate for heterogeneous local sequence organization. It does **not** replace IUPred, SEG, PLAAC, LLPS proxies, AlphaFold/pLDDT, or other existing layers.

### Scientific limitations

- This is a **spin-glass-inspired statistical surrogate**, not molecular dynamics, not a validated phase diagram, and not proof that a protein region is a physical spin glass.
- The model uses **sequence-only categorical encodings** and local shuffled controls; it should be interpreted as a reproducible comparative heuristic, not direct biophysical validation.
- Any apparent signal must be evaluated alongside the existing atlas layers rather than treated as stand-alone evidence for disease mechanism, LLPS, or experimental phenotype.

### Comparison with representative model classes

| Model class / example | Primary input | Typical output | Basis | What this repository adds |
| --- | --- | --- | --- | --- |
| IUPred2A-style disorder predictors | Sequence | Residue disorder propensity | Statistical energetics for disorder tendency | Preserves disorder scores and adds cancer-focused multi-layer atlas context |
| SEG-style low-complexity detectors | Sequence | Low-complexity spans | Compositional complexity filtering | Keeps low-complexity calls as one layer rather than final interpretation |
| PLAAC-style prion-like region models | Sequence | Prion-like/domain scores | Composition and enrichment heuristics | Combines prion-like signal with atlas-wide residue evidence |
| LLPS/FuzDrop-like sequence predictors | Sequence | LLPS propensity or driver-like scores | Sequence-derived phase-separation heuristics | Adds conservation, variants, oncovirus context, and the spin-glass-inspired comparative layer |
| AlphaFold / pLDDT proxies | Sequence/model inference | Structure confidence proxy | Structure-prediction confidence | Uses pLDDT as a complementary override/context layer, not a disorder replacement |
| Sequence protein language models | Sequence | Embeddings or task-specific scores | Large-scale learned statistical representations | Keeps the current pipeline explicit and interpretable, with tunable hand-specified couplings and null controls |

### Reproducibility and interpretation

The new layer is configured in [config.yaml](config.yaml):

```yaml
spin_glass:
  enabled: true
  window_size: 15
  temperature: 1.5
  coupling_scale: 1.0
  field_scale: 0.35
  num_shuffles: 16
  encoding: biophysical_v1
```

It can be reproduced directly with:

```bash
python3 scripts/compute_spin_glass_features.py \
  --fasta data/fasta/selected_proteins.fasta \
  --config config.yaml \
  --output results/spin_glass/spin_glass_scores.tsv
```

Interpret `spin_glass_score` as a relative, null-referenced ranking within this repository's sequence-analysis framework:

- values near 0.5 indicate little deviation from the shuffled local baseline
- higher values indicate stronger combined deviation in energy/frustration/variance/response features
- lower values indicate weaker or more null-like local organization
- final biological interpretation should use the score together with disorder, LLPS, conservation, structural, variant, and virus-aware layers

$$
\mathrm{VIPP}(i)=w_{\mathrm{IDR}}\,I(i)+w_{\mathrm{LLPS}}\,L(i)+w_{\mathrm{virus}}\,V(i)
$$

where $I(i)$ is the normalized disorder score, $L(i)$ is the normalized LLPS proxy score, and $V(i)$ is a binary virus-contact indicator for residue $i$. The default weights are declared in [config.yaml](config.yaml) as $w_{\mathrm{IDR}}=0.35$, $w_{\mathrm{LLPS}}=0.35$, and $w_{\mathrm{virus}}=0.30$, and the resulting score is clipped to $[0,1]$ for interpretability.

### Validation workflow

Run the full validation stack with:

```bash
source .venv-2/bin/activate
pytest -q tests/test_validate_computational_repository.py
python scripts/validate_computational_repository.py
```

This project now includes:

- Snakemake pipeline to download data and run predictors (IUPred2A, SEG, PLAAC) and map ClinVar variants.
- Dockerfile and requirements for reproducible execution.
- Scripts to run batch predictors and aggregate results.
- Integration of phylogeny-aware sequence conservation using orthologs.
- Integration of a deterministic spin-glass-inspired sequence heuristic with shuffled null controls.
- LLPS and structural-proxy analysis with physics-inspired quantum biophysics rules.
- Curated oncovirus-interaction and repurposing screening layers for computational hypothesis generation.
- A curated starter list of 50 high-confidence Cancer Gene Census genes for the MVP.

---

## Licensing & Intellectual Property

**⚠️ CONFIDENTIAL & PROPRIETARY**

- **Copyright © 2026** - All rights reserved
- **Current Status:** Confidential draft (NOT open-source)
- **Embargo Period:** Until publication and copyright/license registration
- **Public Release Date:** TBD (following journal publication)

**For authorization or licensing questions:** s.koganti@[institution].edu

See [LICENSE](./LICENSE) and [CONFIDENTIALITY.md](./CONFIDENTIALITY.md) for full terms.

---

## Access Restrictions

**Access is RESTRICTED to authorized personnel only**

Prohibited without explicit written permission:
- ✗ Public distribution
- ✗ Public repository posting (GitHub, Zenodo, etc.)
- ✗ Sharing with unauthorized persons
- ✗ Commercial use
- ✗ Any use prior to formal publication

---

**License:** Proprietary and Confidential (See LICENSE)
