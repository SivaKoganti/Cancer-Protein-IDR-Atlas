from pathlib import Path
import subprocess
import sys

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_analyze_idr_mutant_pathway_effects_generates_outputs(tmp_path):
    fasta_path = tmp_path / "selected_proteins.fasta"
    fasta_path.write_text(
        ">TP53\nQQQQQQQQ\n>KRAS\nDEDEDEDE\n",
        encoding="utf-8",
    )

    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        """
llps:
  window: 5
  qn_weight: 0.45
  aromatic_weight: 0.25
  hydrophobic_weight: 0.15
  charge_pattern_weight: 0.15
oncogenic_pathways:
  TP53_AXIS:
    - TP53
  RTK_RAS_MAPK:
    - KRAS
""".strip()
        + "\n",
        encoding="utf-8",
    )

    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()

    tp53_atlas = pd.DataFrame([
      {"gene": "TP53", "pos": 1, "aa": "Q", "iupred_score": 0.8, "ptm_count": 1, "ptm_categories": "phospho", "ptm_ids": "MOD_CDK_1", "ptm_in_idr_count": 1, "ptm_clinvar_overlap_count": 0},
      {"gene": "TP53", "pos": 2, "aa": "Q", "iupred_score": 0.7, "ptm_count": 0, "ptm_categories": "", "ptm_ids": "", "ptm_in_idr_count": 0, "ptm_clinvar_overlap_count": 0},
    ])
    tp53_atlas.to_csv(atlas_dir / "TP53_atlas.tsv", sep="\t", index=False)

    kras_atlas = pd.DataFrame([
      {"gene": "KRAS", "pos": 1, "aa": "D", "iupred_score": 0.85, "ptm_count": 0, "ptm_categories": "", "ptm_ids": "", "ptm_in_idr_count": 0, "ptm_clinvar_overlap_count": 0},
      {"gene": "KRAS", "pos": 2, "aa": "E", "iupred_score": 0.9, "ptm_count": 2, "ptm_categories": "sumo", "ptm_ids": "MOD_SUMO_1", "ptm_in_idr_count": 2, "ptm_clinvar_overlap_count": 1},
    ])
    kras_atlas.to_csv(atlas_dir / "KRAS_atlas.tsv", sep="\t", index=False)

    output_dir = tmp_path / "mutants"

    completed = subprocess.run([
        sys.executable,
        str(REPO_ROOT / "scripts" / "analyze_idr_mutant_pathway_effects.py"),
        "--fasta",
        str(fasta_path),
        "--atlas-dir",
        str(atlas_dir),
        "--config",
        str(config_path),
        "--output-dir",
        str(output_dir),
        "--idr-threshold",
        "0.5",
    ], capture_output=True, text=True, check=False)

    assert completed.returncode == 0, completed.stderr

    mutants_path = output_dir / "idr_llps_mutants.tsv"
    pathway_path = output_dir / "pathway_differential_llps.tsv"
    pathway_ptm_path = output_dir / "pathway_ptm_differential_llps.tsv"
    report_path = output_dir / "pathway_critical_analysis.md"

    assert mutants_path.exists()
    assert pathway_path.exists()
    assert pathway_ptm_path.exists()
    assert report_path.exists()

    mutants = pd.read_csv(mutants_path, sep="\t")
    assert not mutants.empty
    assert {"gene", "pathways", "mutation", "direction", "delta_llps", "high_impact", "ptm_count", "ptm_context"}.issubset(mutants.columns)
    assert set(mutants["direction"].unique()) == {"increase", "decrease"}

    pathway_summary = pd.read_csv(pathway_path, sep="\t")
    assert not pathway_summary.empty
    assert {"pathway", "direction", "differential_mean_delta", "vulnerability_index"}.issubset(pathway_summary.columns)
    assert {"TP53_AXIS", "RTK_RAS_MAPK"}.issubset(set(pathway_summary["pathway"]))

    pathway_ptm_summary = pd.read_csv(pathway_ptm_path, sep="\t")
    assert not pathway_ptm_summary.empty
    assert {"pathway", "direction", "ptm_context", "ptm_differential_mean_abs_delta"}.issubset(pathway_ptm_summary.columns)

    report_text = report_path.read_text(encoding="utf-8")
    assert "# IDR Mutant LLPS Differential Analysis" in report_text
    assert "## Critical Caveats" in report_text
