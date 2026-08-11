# Computational Novelty Sprint Plan

This document is a machine-runnable milestone plan for making the study computationally unique, robust, and publication-ready.

## Scope

- Fully computational only (no wet-lab steps)
- Deterministic, reproducible outputs
- Quantitative novelty evidence through benchmarking and ablations

## Global Success Gates

1. Method novelty gate:
   - Virus + IDR + LLPS + conservation + discordance all active in final scoring table.
2. Biological confidence gate:
   - Non-placeholder IDR values, non-empty pLDDT, and non-zero CDR/discordant regions.
3. Reproducibility gate:
   - Full rerun reproduces summary metrics and top-ranked residues within tolerance.

## Environment Setup

Run from repository root:

```bash
set -euo pipefail

# Activate project Python environment
source .venv-2/bin/activate

# Quick health checks
python --version
snakemake --version
```

Expected outcome:
- Python and Snakemake available without errors.

---

## Milestone M1: Activate Non-Placeholder Core Signals

Objective:
- Ensure IDR and pLDDT are real-valued and usable for novelty layers.

Commands:

```bash
set -euo pipefail
source .venv-2/bin/activate

# Rebuild core upstream layers
snakemake --cores 1 results/iupred/iupred_scores.tsv results/alphafold/plddt_scores.tsv

# Rebuild novelty-dependent layers
snakemake --cores 1 results/cdr/cdr_summary.tsv results/cdr/all_genes_cdr.tsv

# Rebuild atlas summaries
snakemake --cores 1 results/atlas/global_disorder_phylogeny_atlas.tsv
```

Verification:

```bash
python - <<'PY'
import pandas as pd
from pathlib import Path

iup = pd.read_csv('results/iupred/iupred_scores.tsv', sep='\t')
print('iupred_unique_scores', iup['score'].nunique())

pl = Path('results/alphafold/plddt_scores.tsv')
print('plddt_exists', pl.exists(), 'size', pl.stat().st_size if pl.exists() else -1)

cdr = pd.read_csv('results/cdr/cdr_summary.tsv', sep='\t')
print('genes_with_cdr', int((cdr['cdr_count'] > 0).sum()))
print('genes_with_discordant', int((cdr['discordant_count'] > 0).sum()))
PY
```

Pass criteria:
- iupred_unique_scores > 1
- plddt file size > 1 byte
- genes_with_cdr > 0
- genes_with_discordant > 0

Primary outputs:
- results/iupred/iupred_scores.tsv
- results/alphafold/plddt_scores.tsv
- results/cdr/cdr_summary.tsv
- results/cdr/all_genes_cdr.tsv

---

## Milestone M2: Virus Layer Depth and Evidence Weighting

Objective:
- Strengthen virus-host signal quality and computational traceability.

Commands:

```bash
set -euo pipefail
source .venv-2/bin/activate

# Recompute virus mapping layer
snakemake --cores 1 results/oncovirus/oncovirus_per_residue.tsv

# Rebuild atlas with updated virus layer
snakemake --cores 1 results/atlas/global_disorder_phylogeny_atlas.tsv
```

Verification:

```bash
python - <<'PY'
import pandas as pd

onc = pd.read_csv('results/oncovirus/oncovirus_per_residue.tsv', sep='\t')
print('rows', len(onc))
print('genes', onc['gene'].nunique())
print('interaction_residues', int(onc['virus_interaction'].sum()))
print('genes_with_interactions', int(onc.groupby('gene')['virus_interaction'].max().sum()))

summary = pd.read_csv('results/atlas/global_disorder_phylogeny_atlas.tsv', sep='\t')
for c in ['virus_interaction_fraction','mean_virus_count','mean_vipp_score','max_vipp_score']:
    print(c, c in summary.columns)
PY
```

Pass criteria:
- interaction_residues > 0
- genes_with_interactions > 0
- virus and VIPP columns present in global summary

Primary outputs:
- results/oncovirus/oncovirus_per_residue.tsv
- results/atlas/global_disorder_phylogeny_atlas.tsv

---

## Milestone M3: VIPP Ranking and Hotspot Stability

Objective:
- Produce robust, presentation-ready computational prioritization.

Commands:

```bash
set -euo pipefail
source .venv-2/bin/activate

# Build all VIPP-relevant artifacts
snakemake --cores 1 results/vipp/top_vipp_residues.tsv results/map/chromosome_atlas.html results/index.html
```

Verification:

```bash
python - <<'PY'
import pandas as pd

v = pd.read_csv('results/vipp/top_vipp_residues.tsv', sep='\t')
print('rows', len(v))
print('unique_genes', v['gene'].nunique())
print('max_vipp', float(v['vipp_score'].max()))
print(v.head(5).to_string(index=False))
PY
```

Pass criteria:
- rows == configured global top (default 200)
- unique_genes > 5
- max_vipp > 0

Primary outputs:
- results/vipp/top_vipp_residues.tsv
- results/map/chromosome_atlas.html
- results/index.html

---

## Milestone M4: Benchmark and Ablation Matrix

Objective:
- Quantify novelty contribution of each computational layer.

Required protocol:
1. Baseline model: IDR only
2. Add LLPS
3. Add virus layer
4. Add conservation
5. Add CDR/discordant
6. Full model

Suggested command template:

```bash
# Example placeholder for each model variant
# python scripts/benchmark_vipp_model.py --variant idr_only --out results/benchmark/idr_only.tsv
# python scripts/benchmark_vipp_model.py --variant idr_llps --out results/benchmark/idr_llps.tsv
# ...
```

Required outputs:
- results/benchmark/ablation_matrix.tsv
- results/benchmark/model_metrics.tsv
- results/benchmark/ranking_stability.tsv

Required metrics columns:
- model_name
- auroc
- auprc
- brier_score
- calibration_error
- topk_enrichment
- delta_vs_baseline

Pass criteria:
- Full model outperforms baseline on at least one primary endpoint.
- Effect sizes and confidence intervals are reported.

---

## Milestone M5: Deterministic Release Bundle

Objective:
- Freeze computational release with complete provenance.

Commands:

```bash
set -euo pipefail
source .venv-2/bin/activate

# Full target regeneration
snakemake --cores 1 \
  results/atlas/global_disorder_phylogeny_atlas.tsv \
  results/vipp/top_vipp_residues.tsv \
  results/map/chromosome_atlas.html \
  results/index.html

# Checksums for key release artifacts
mkdir -p results/release
sha256sum \
  results/atlas/global_disorder_phylogeny_atlas.tsv \
  results/vipp/top_vipp_residues.tsv \
  results/map/chromosome_atlas.html \
  results/index.html \
  > results/release/release_checksums.sha256
```

Verification:

```bash
snakemake --cores 1 \
  results/atlas/global_disorder_phylogeny_atlas.tsv \
  results/vipp/top_vipp_residues.tsv \
  results/map/chromosome_atlas.html \
  results/index.html
```

Pass criteria:
- Snakemake reports all requested files up to date.
- Checksums file exists and is non-empty.

Primary outputs:
- results/release/release_checksums.sha256

---

## Final Acceptance Checklist

- [ ] M1 passed
- [ ] M2 passed
- [ ] M3 passed
- [ ] M4 passed
- [ ] M5 passed
- [ ] Novelty scorecard updated in manuscript and readme from current run outputs

## Current Known Risk Flags (from latest run)

1. Uniform IDR score currently detected (placeholder behavior).
2. pLDDT file currently empty in this environment.
3. CDR and discordant layers currently inactive due missing pLDDT and placeholder IDR state.

These must be resolved before claiming high biological-confidence novelty.
