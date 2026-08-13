from pathlib import Path
import subprocess
import sys

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_summarize_case_studies_writes_table(tmp_path):
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()

    tp53 = pd.DataFrame([
        {"gene": "TP53", "pos": 1, "llps_proxy": 0.2, "vipp_score": 0.15, "virus_interaction": 0, "virus_count": 0},
        {"gene": "TP53", "pos": 20, "llps_proxy": 0.7, "vipp_score": 0.72, "virus_interaction": 1, "virus_count": 2},
    ])
    kras = pd.DataFrame([
        {"gene": "KRAS", "pos": 12, "llps_proxy": 0.65, "vipp_score": 0.63, "virus_interaction": 1, "virus_count": 1},
    ])

    tp53.to_csv(atlas_dir / "TP53_atlas.tsv", sep="\t", index=False)
    kras.to_csv(atlas_dir / "KRAS_atlas.tsv", sep="\t", index=False)

    output = tmp_path / "case_studies.tsv"
    completed = subprocess.run([
        sys.executable,
        str(REPO_ROOT / "scripts" / "summarize_case_studies.py"),
        "--atlas-dir",
        str(atlas_dir),
        "--output",
        str(output),
    ], capture_output=True, text=True, check=False)

    assert completed.returncode == 0, completed.stderr
    table = pd.read_csv(output, sep="\t")
    assert not table.empty
    assert {"gene", "pos", "llps_proxy", "vipp_score", "virus_interaction"}.issubset(table.columns)
    assert table["gene"].isin(["TP53", "KRAS"]).all()
