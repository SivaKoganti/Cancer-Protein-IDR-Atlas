from pathlib import Path
import subprocess
import sys

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_build_clinvar_view_integrates_variants(tmp_path):
    atlas_path = tmp_path / "TP53_atlas.tsv"
    atlas = pd.DataFrame(
        [
            {"gene": "TP53", "pos": 12, "aa": "A", "iupred_score": 0.2},
            {"gene": "TP53", "pos": 13, "aa": "B", "iupred_score": 0.4},
        ]
    )
    atlas.to_csv(atlas_path, sep="\t", index=False)

    variants_path = tmp_path / "variants.tsv"
    pd.DataFrame(
        [
            ["TP53", 12, "R", "H", "Pathogenic"],
            ["TP53", 12, "R", "C", "Likely pathogenic"],
        ]
    ).to_csv(variants_path, sep="\t", header=False, index=False)

    outdir = tmp_path / "out"
    cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "build_clinvar_view.py"),
        "--atlas",
        str(atlas_path),
        "--variants",
        str(variants_path),
        "--outdir",
        str(outdir),
    ]

    completed = subprocess.run(cmd, capture_output=True, text=True, check=False)

    assert completed.returncode == 0, completed.stderr
    merged_path = outdir / "TP53_atlas_with_clinvar.tsv"
    assert merged_path.exists(), "expected merged ClinVar atlas output"

    merged = pd.read_csv(merged_path, sep="\t")
    row = merged.loc[merged["pos"] == 12].iloc[0]
    assert row["clinvar_count"] == 2
    assert "Pathogenic" in row["clinvar_significance"]
    assert "Likely pathogenic" in row["clinvar_significance"]
