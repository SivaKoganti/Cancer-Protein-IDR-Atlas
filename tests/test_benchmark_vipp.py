from pathlib import Path
import importlib.util

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("benchmark_vipp", REPO_ROOT / "scripts" / "benchmark_vipp.py")
benchmark_vipp = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(benchmark_vipp)
benchmark_models = benchmark_vipp.benchmark_models


def test_benchmark_models_writes_summary_tables(tmp_path):
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()

    tp53 = pd.DataFrame([
        {"gene": "TP53", "pos": 1, "aa": "M", "iupred_score": 0.8, "llps_proxy": 0.2, "virus_interaction": 0, "conservation": 0.1, "vipp_score": 0.25, "clinvar_significance": "Pathogenic"},
        {"gene": "TP53", "pos": 2, "aa": "A", "iupred_score": 0.4, "llps_proxy": 0.4, "virus_interaction": 1, "conservation": 0.2, "vipp_score": 0.55, "clinvar_significance": "Benign"},
        {"gene": "TP53", "pos": 3, "aa": "C", "iupred_score": 0.6, "llps_proxy": 0.3, "virus_interaction": 0, "conservation": 0.3, "vipp_score": 0.45, "clinvar_significance": "Pathogenic/Likely pathogenic"},
    ])
    kras = pd.DataFrame([
        {"gene": "KRAS", "pos": 12, "aa": "G", "iupred_score": 0.5, "llps_proxy": 0.5, "virus_interaction": 0, "conservation": 0.4, "vipp_score": 0.35, "clinvar_significance": "Benign"},
        {"gene": "KRAS", "pos": 13, "aa": "D", "iupred_score": 0.7, "llps_proxy": 0.3, "virus_interaction": 1, "conservation": 0.6, "vipp_score": 0.6, "clinvar_significance": "Likely pathogenic"},
    ])

    tp53.to_csv(atlas_dir / "TP53_atlas_with_clinvar.tsv", sep="\t", index=False)
    kras.to_csv(atlas_dir / "KRAS_atlas_with_clinvar.tsv", sep="\t", index=False)

    output_dir = tmp_path / "benchmark"
    metrics, ablation = benchmark_models(atlas_dir, output_dir)

    assert not metrics.empty
    assert {"model_name", "auroc", "auprc", "topk_enrichment"}.issubset(metrics.columns)
    assert not ablation.empty
    assert "model_name" in ablation.columns
    assert (metrics["model_name"] == "full_vipp").any()
