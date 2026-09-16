from pathlib import Path
import subprocess
import sys

import pandas as pd

from scripts.compute_spin_glass_features import build_output_dataframe, load_spin_glass_config, residue_vector


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_residue_vector_unknown_defaults_to_zero_state():
    assert residue_vector("X") == (0.0, 0.0, 0.0, 0.0)


def test_load_spin_glass_config_normalizes_even_windows(tmp_path):
    config_path = tmp_path / "config.yaml"
    config_path.write_text("spin_glass:\n  window_size: 8\n  num_shuffles: 4\n", encoding="utf-8")

    cfg = load_spin_glass_config(config_path)

    assert cfg["window_size"] == 9
    assert cfg["num_shuffles"] == 4


def test_compute_spin_glass_script_is_deterministic_and_handles_short_unknown_sequences(tmp_path):
    fasta_path = tmp_path / "sample.fasta"
    fasta_path.write_text(">TP53_1\nAXD\n>KRAS\nG\n", encoding="utf-8")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "spin_glass:\n"
        "  enabled: true\n"
        "  window_size: 5\n"
        "  temperature: 1.2\n"
        "  coupling_scale: 0.8\n"
        "  field_scale: 0.25\n"
        "  num_shuffles: 6\n",
        encoding="utf-8",
    )
    output_a = tmp_path / "spin_a.tsv"
    output_b = tmp_path / "spin_b.tsv"

    for output_path in (output_a, output_b):
        completed = subprocess.run([
            sys.executable,
            str(REPO_ROOT / "scripts" / "compute_spin_glass_features.py"),
            "--fasta",
            str(fasta_path),
            "--config",
            str(config_path),
            "--output",
            str(output_path),
        ], capture_output=True, text=True, check=False)
        assert completed.returncode == 0, completed.stderr

    df_a = pd.read_csv(output_a, sep="\t")
    df_b = pd.read_csv(output_b, sep="\t")

    pd.testing.assert_frame_equal(df_a, df_b)
    assert len(df_a) == 4
    assert df_a["gene"].tolist() == ["TP53", "TP53", "TP53", "KRAS"]
    assert df_a["window_length"].tolist() == [3, 3, 3, 1]
    assert df_a["spin_glass_score"].between(0, 1).all()
    numeric_columns = [
        "spin_state",
        "local_energy",
        "local_frustration",
        "coupling_mean",
        "coupling_variance",
        "susceptibility_like",
        "null_energy_mean",
        "null_energy_std",
        "null_frustration_mean",
        "null_frustration_std",
        "null_coupling_variance_mean",
        "null_coupling_variance_std",
        "null_susceptibility_mean",
        "null_susceptibility_std",
        "spin_glass_score",
    ]
    assert df_a[numeric_columns].notna().all().all()


def test_compute_spin_glass_disabled_writes_empty_schema(tmp_path):
    fasta_path = tmp_path / "sample.fasta"
    fasta_path.write_text(">TP53\nACDE\n", encoding="utf-8")
    config_path = tmp_path / "config.yaml"
    config_path.write_text("spin_glass:\n  enabled: false\n", encoding="utf-8")

    df = build_output_dataframe(fasta_path, config_path)

    assert df.empty
    assert {"gene", "spin_glass_score", "spin_glass_method", "spin_glass_version"}.issubset(df.columns)
