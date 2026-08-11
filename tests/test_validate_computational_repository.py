from pathlib import Path
import subprocess
import sys

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_validate_computational_repository_writes_report(tmp_path):
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    pd.DataFrame([
        {"gene": "TP53", "pos": 1, "llps_proxy": 0.25, "virus_interaction": 1, "virus_count": 2, "vipp_score": 0.7},
    ]).to_csv(atlas_dir / "TP53_atlas.tsv", sep="\t", index=False)

    results_dir = tmp_path / "results"
    llps_dir = results_dir / "llps"
    benchmark_dir = results_dir / "benchmark"
    llps_dir.mkdir(parents=True)
    benchmark_dir.mkdir(parents=True)

    pd.DataFrame([
        {"llps_class": "low", "residues": 10, "virus_linked_residues": 1, "max_llps": 0.3, "candidate_compounds": "test", "rationale": "test"},
    ]).to_csv(llps_dir / "llps_compound_library.tsv", sep="\t", index=False)

    pd.DataFrame([
        {"model_name": "full_vipp", "auroc": 0.8, "auprc": 0.7, "topk_enrichment": 2.0},
    ]).to_csv(benchmark_dir / "model_metrics.tsv", sep="\t", index=False)

    pd.DataFrame([
        {"model_name": "full_vipp", "ablation": "full", "score": 0.8},
    ]).to_csv(benchmark_dir / "ablation_matrix.tsv", sep="\t", index=False)

    repurposing = results_dir / "llps" / "ncats_repurposing_candidates.tsv"
    repurposing.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([
        {"compound_id": "CMP001", "compound_name": "Test", "validated_for_repurpose": True},
    ]).to_csv(repurposing, sep="\t", index=False)

    output = tmp_path / "validation.tsv"
    completed = subprocess.run([
        sys.executable,
        str(REPO_ROOT / "scripts" / "validate_computational_repository.py"),
        "--atlas-dir",
        str(atlas_dir),
        "--results-dir",
        str(results_dir),
        "--output",
        str(output),
    ], capture_output=True, text=True, check=False)

    assert completed.returncode == 0, completed.stderr
    report = pd.read_csv(output, sep="\t")
    assert not report.empty
    assert report["status"].isin(["PASS", "FAIL"]).all()
    assert (report["status"] == "PASS").any()
