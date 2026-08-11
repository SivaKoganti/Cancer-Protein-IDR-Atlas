from pathlib import Path
import pandas as pd
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_build_idr_structural_library_writes_residue_state_library(tmp_path):
    atlas_path = tmp_path / "sample_atlas.tsv"
    atlas_path.write_text(
        "gene\tpos\taa\tiupred_score\tllps_proxy\tstructural_proxy\tquantum_proxy\n"
        "TP53\t1\tM\t0.82\t0.71\t0.21\t0.53\n"
        "TP53\t2\tA\t0.64\t0.36\t0.48\t0.41\n",
        encoding="utf-8",
    )

    output_path = tmp_path / "library.tsv"
    completed = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "build_idr_structural_library.py"),
            str(atlas_path),
            str(output_path),
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    library_df = pd.read_csv(output_path, sep="\t")
    assert not library_df.empty
    assert {"gene", "pos", "structural_library_json", "coarse_model_xyz", "simulated_crystallogram", "quantum_biophysical_score", "llps_state"}.issubset(library_df.columns)
    assert library_df.loc[0, "llps_state"] in {"compact", "expanded", "condensate"}
    assert library_df.loc[0, "quantum_biophysical_score"] >= 0.0
