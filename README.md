# ⚠️ CONFIDENTIAL - Cancer Protein IDR & LLMS Atlas

**⚠️ CONFIDENTIAL MATERIAL - NOT FOR PUBLIC DISTRIBUTION**

This repository is currently under confidential development and embargo. Unauthorized access, copying, or distribution is strictly prohibited.

**Authorized Personnel Only** - Contact s.koganti@[institution].edu for access

---

## Project Description (Internal)

This repository contains a computationally validated pipeline for intrinsic disorder region (IDR), low-complexity region (LCR), and liquid–liquid phase separation (LLPS) analyses in cancer-associated human proteins. It integrates ClinVar evidence, phylogenetic conservation, oncovirus interaction mapping, VIPP prioritization, and repurposing screening into a reproducible computational framework.

The current build is validated through generated artifacts and a repository-level validation report in [results/validation/computational_repository_validation.tsv](results/validation/computational_repository_validation.tsv).

### Core computational outputs

- Atlas tables per gene in [results/atlas](results/atlas)
- LLPS compound library in [results/llps/llps_compound_library.tsv](results/llps/llps_compound_library.tsv)
- VIPP benchmark metrics in [results/benchmark/model_metrics.tsv](results/benchmark/model_metrics.tsv) with virus-aware models showing near-perfect ranking performance on the current benchmark set
- Repurposing candidates in [results/llps/ncats_repurposing_candidates.tsv](results/llps/ncats_repurposing_candidates.tsv)
- Repository validation report in [results/validation/computational_repository_validation.tsv](results/validation/computational_repository_validation.tsv)

### VIPP scoring model

The repository now computes a formal, configurable VIPP score for each residue as a weighted combination of normalized evidence channels:

### Novelty relative to existing pipelines

This project is not merely a wrapper around standard disorder or LLPS predictors. It is a cancer-gene-specific atlas that combines residue-level disorder, LLPS, conservation, structural proxies, ClinVar context, and curated oncovirus interaction evidence in one reproducible framework, then ranks residues with a transparent VIPP score and bench-marks the contribution of each evidence layer through ablation analyses.

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
