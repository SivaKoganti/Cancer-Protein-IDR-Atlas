from pathlib import Path
import subprocess
import sys

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_build_disorder_atlas_generates_tables_from_available_inputs(tmp_path):
    outdir = tmp_path / "atlas"
    cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "build_disorder_atlas.py"),
        "--iupred",
        str(REPO_ROOT / "results" / "iupred" / "iupred_scores.tsv"),
        "--seg",
        str(REPO_ROOT / "results" / "seg" / "seg_regions.tsv"),
        "--plaac",
        str(REPO_ROOT / "results" / "plaac" / "plaac_scores.tsv"),
        "--llps",
        str(REPO_ROOT / "results" / "llps" / "llps_scores.tsv"),
        "--structure",
        str(REPO_ROOT / "results" / "structure" / "structure_scores.tsv"),
        "--variants",
        str(REPO_ROOT / "results" / "clinvar" / "mapped_variants.tsv"),
        "--conservation",
        str(REPO_ROOT / "results" / "phylogeny" / "conservation.tsv"),
        "--fasta",
        str(REPO_ROOT / "data" / "fasta" / "selected_proteins.fasta"),
        "--genes-file",
        str(REPO_ROOT / "config.yaml"),
        "--outdir",
        str(outdir),
    ]

    completed = subprocess.run(cmd, capture_output=True, text=True, check=False)

    assert completed.returncode == 0, completed.stderr
    atlas_files = sorted(path.name for path in outdir.glob("*_atlas.tsv"))
    assert atlas_files, "expected atlas TSVs to be generated"

    summary = pd.read_csv(outdir / "global_disorder_phylogeny_atlas.tsv", sep="\t")
    assert not summary.empty
    assert "gene" in summary.columns


def test_build_disorder_atlas_normalizes_isoform_ids_and_ignores_invalid_optional_rows(tmp_path):
    outdir = tmp_path / "atlas"
    config_path = tmp_path / "config.yaml"
    config_path.write_text("genes:\n  - TP53\n", encoding="utf-8")

    fasta_path = tmp_path / "sample.fasta"
    fasta_path.write_text(">TP53_1\nACDEFGHIKLMN\n", encoding="utf-8")

    (tmp_path / "iupred.tsv").write_text("gene\tpos\tscore\nTP53_1\t1\t0.9\n", encoding="utf-8")
    (tmp_path / "seg.tsv").write_text("gene\tstart\tend\nTP53_1\t1\t2\n", encoding="utf-8")
    (tmp_path / "plaac.tsv").write_text("gene\tstart\tend\tq_n_fraction\nTP53_1\t1\t2\t0.5\n", encoding="utf-8")
    (tmp_path / "llps.tsv").write_text("gene\tstart\tend\tllps_proxy\nTP53_1\t1\t2\t0.6\n", encoding="utf-8")
    (tmp_path / "structure.tsv").write_text("gene\tpos\tstructural_proxy\tquantum_proxy\nTP53_1\t1\t0.1\t0.2\n", encoding="utf-8")
    (tmp_path / "conservation.tsv").write_text("gene\tpos\tconservation_score\nTP53_1\t1\t0.7\n", encoding="utf-8")
    (tmp_path / "plddt.tsv").write_text("", encoding="utf-8")
    (tmp_path / "delta_llps.tsv").write_text("gene\tpos\tllps_vulnerability\tdelta_llps_max\tdelta_llps_min\nTP53_1\tNaN\t0.1\t0.2\t0.3\n", encoding="utf-8")
    (tmp_path / "cdr.tsv").write_text("", encoding="utf-8")
    (tmp_path / "slim.tsv").write_text("", encoding="utf-8")
    (tmp_path / "ptm.tsv").write_text("gene\tpos\tptm_count\tptm_categories\tptm_ids\tptm_in_idr_count\tptm_clinvar_overlap_count\nTP53_1\t1\t1\tphospho\tMOD_CDK_1\t1\t0\n", encoding="utf-8")

    completed = subprocess.run([
        sys.executable,
        str(REPO_ROOT / "scripts" / "build_disorder_atlas.py"),
        "--iupred",
        str(tmp_path / "iupred.tsv"),
        "--seg",
        str(tmp_path / "seg.tsv"),
        "--plaac",
        str(tmp_path / "plaac.tsv"),
        "--llps",
        str(tmp_path / "llps.tsv"),
        "--structure",
        str(tmp_path / "structure.tsv"),
        "--conservation",
        str(tmp_path / "conservation.tsv"),
        "--fasta",
        str(fasta_path),
        "--genes-file",
        str(config_path),
        "--plddt",
        str(tmp_path / "plddt.tsv"),
        "--delta-llps",
        str(tmp_path / "delta_llps.tsv"),
        "--cdr",
        str(tmp_path / "cdr.tsv"),
        "--slim",
        str(tmp_path / "slim.tsv"),
        "--ptm",
        str(tmp_path / "ptm.tsv"),
        "--outdir",
        str(outdir),
    ], capture_output=True, text=True, check=False)

    assert completed.returncode == 0, completed.stderr
    atlas_path = outdir / "TP53_atlas.tsv"
    assert atlas_path.exists(), "expected atlas table for normalized gene"

    atlas_df = pd.read_csv(atlas_path, sep="\t")
    assert atlas_df["gene"].eq("TP53").all()
    assert atlas_df["pos"].tolist()[:2] == [1, 2]
    assert {"ptm_count", "ptm_categories", "ptm_ids"}.issubset(atlas_df.columns)
    assert {"vipp_score"}.issubset(atlas_df.columns)
    assert atlas_df["vipp_score"].between(0, 1).all()

    summary_df = pd.read_csv(outdir / "global_disorder_phylogeny_atlas.tsv", sep="\t")
    assert {"mean_vipp_score", "max_vipp_score", "virus_interaction_fraction", "mean_virus_count"}.issubset(summary_df.columns)
    assert summary_df["mean_vipp_score"].between(0, 1).all()


def test_scan_slim_motifs_writes_empty_outputs_when_no_hits_are_found(tmp_path):
    fasta_path = tmp_path / "sample.fasta"
    fasta_path.write_text(">TP53_1\nACDEFGHIKLMN\n", encoding="utf-8")
    iupred_path = tmp_path / "iupred.tsv"
    iupred_path.write_text("gene\tpos\tscore\nTP53_1\t1\t0.9\n", encoding="utf-8")

    outdir = tmp_path / "slim"
    completed = subprocess.run([
        sys.executable,
        str(REPO_ROOT / "scripts" / "scan_slim_motifs.py"),
        "--fasta",
        str(fasta_path),
        "--iupred",
        str(iupred_path),
        "--output-dir",
        str(outdir),
    ], capture_output=True, text=True, check=False)

    assert completed.returncode == 0, completed.stderr
    assert (outdir / "slim_hits.tsv").exists()
    assert (outdir / "slim_per_residue.tsv").exists()
    assert (outdir / "ptm_hits.tsv").exists()
    assert (outdir / "ptm_per_residue.tsv").exists()
    assert (outdir / "ptm_category_summary.tsv").exists()


def test_scan_slim_motifs_normalizes_isoform_gene_ids(tmp_path):
    fasta_path = tmp_path / "sample.fasta"
    fasta_path.write_text(">TP53_1\nSPVA\n", encoding="utf-8")
    iupred_path = tmp_path / "iupred.tsv"
    iupred_path.write_text("gene\tpos\tscore\nTP53_1\t1\t0.9\nTP53_1\t2\t0.8\n", encoding="utf-8")
    variants_path = tmp_path / "variants.tsv"
    variants_path.write_text("TP53_1\t1\tA\tG\tPathogenic\n", encoding="utf-8")

    outdir = tmp_path / "slim"
    completed = subprocess.run([
        sys.executable,
        str(REPO_ROOT / "scripts" / "scan_slim_motifs.py"),
        "--fasta",
        str(fasta_path),
        "--iupred",
        str(iupred_path),
        "--variants",
        str(variants_path),
        "--output-dir",
        str(outdir),
    ], capture_output=True, text=True, check=False)

    assert completed.returncode == 0, completed.stderr
    hits_df = pd.read_csv(outdir / "slim_hits.tsv", sep="\t")
    assert not hits_df.empty
    assert hits_df["mean_idr"].max() > 0.0
    assert hits_df["clinvar_overlap"].sum() == 1

    ptm_df = pd.read_csv(outdir / "ptm_hits.tsv", sep="\t")
    assert not ptm_df.empty
    assert ptm_df["motif_id"].str.startswith("MOD_").all()
