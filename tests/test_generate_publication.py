from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts import generate_publication


def test_figure_7_structural_library_is_generated(tmp_path):
    output_dir = tmp_path / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    generate_publication.figure_7_structural_library(
        REPO_ROOT / "results" / "atlas",
        output_dir,
    )

    assert (output_dir / "figure_7_structural_library.png").exists()


def test_figure_9_viable_mutant_llp_phase_differential_is_generated(tmp_path):
    output_dir = tmp_path / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    generate_publication.figure_9_viable_mutant_llp_phase_differential(
        REPO_ROOT / "results" / "mutants" / "idr_llps_mutants.tsv",
        output_dir,
    )

    assert (output_dir / "figure_9_viable_mutant_llp_phase_differential.png").exists()
