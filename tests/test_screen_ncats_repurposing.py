from pathlib import Path
import subprocess
import sys

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_screen_ncats_repurposing_generates_candidates(tmp_path):
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()

    tp53 = pd.DataFrame([
        {"gene": "TP53", "pos": 1, "aa": "M", "llps_proxy": 0.25, "virus_interaction": 0, "virus_count": 0, "virus_names": "", "vipp_score": 0.2},
        {"gene": "TP53", "pos": 2, "aa": "A", "llps_proxy": 0.62, "virus_interaction": 1, "virus_count": 2, "virus_names": "HPV16;EBV", "vipp_score": 0.7},
        {"gene": "TP53", "pos": 3, "aa": "P", "llps_proxy": 0.81, "virus_interaction": 1, "virus_count": 1, "virus_names": "HPV16", "vipp_score": 0.8},
    ])
    tp53.to_csv(atlas_dir / "TP53_atlas.tsv", sep="\t", index=False)

    compounds = pd.DataFrame([
        {"compound_id": "CMP001", "compound_name": "Hexanediol", "target": "phase separation", "mechanism": "LLPS disruptor", "phase": "preclinical", "source": "example"},
        {"compound_id": "CMP002", "compound_name": "Chloroquine", "target": "viral replication", "mechanism": "autophagy modulator", "phase": "approved", "source": "example"},
        {"compound_id": "CMP003", "compound_name": "Dasatinib", "target": "kinase", "mechanism": "tyrosine kinase inhibitor", "phase": "approved", "source": "example"},
    ])
    compounds.to_csv(tmp_path / "ncats_compounds.tsv", sep="\t", index=False)

    output = tmp_path / "repurposing.tsv"
    completed = subprocess.run([
        sys.executable,
        str(REPO_ROOT / "scripts" / "screen_ncats_repurposing.py"),
        "--atlas-dir",
        str(atlas_dir),
        "--compound-db",
        str(tmp_path / "ncats_compounds.tsv"),
        "--output",
        str(output),
    ], capture_output=True, text=True, check=False)

    assert completed.returncode == 0, completed.stderr
    library = pd.read_csv(output, sep="\t")
    assert not library.empty
    assert {"compound_id", "compound_name", "validated_for_repurpose"}.issubset(library.columns)
    assert library["validated_for_repurpose"].astype(bool).any()
