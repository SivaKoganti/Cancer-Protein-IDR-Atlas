from pathlib import Path
import subprocess
import sys

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_generate_llps_compound_library_writes_summary(tmp_path):
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()

    tp53 = pd.DataFrame([
        {"gene": "TP53", "pos": 1, "aa": "M", "llps_proxy": 0.2, "virus_interaction": 0, "virus_count": 0},
        {"gene": "TP53", "pos": 2, "aa": "A", "llps_proxy": 0.7, "virus_interaction": 1, "virus_count": 2},
    ])
    tp53.to_csv(atlas_dir / "TP53_atlas.tsv", sep="\t", index=False)

    output = tmp_path / "library.tsv"
    completed = subprocess.run([
        sys.executable,
        str(REPO_ROOT / "scripts" / "generate_llps_compound_library.py"),
        "--atlas-dir",
        str(atlas_dir),
        "--output",
        str(output),
    ], capture_output=True, text=True, check=False)

    assert completed.returncode == 0, completed.stderr
    library = pd.read_csv(output, sep="\t")
    assert not library.empty
    assert {"llps_class", "candidate_compounds"}.issubset(library.columns)
    assert library["llps_class"].isin(["low", "moderate", "high"]).all()
