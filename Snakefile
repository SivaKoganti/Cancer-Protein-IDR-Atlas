# Snakemake workflow for the Cancer Protein IDR atlas with phylogeny-aware conservation

configfile: "config.yaml"

rule all:
  input:
    expand("results/atlas/{gene}_atlas.tsv", gene=config["genes"]),
    expand("results/atlas/{gene}_atlas_with_clinvar.tsv", gene=config["genes"]),
    "results/atlas/global_disorder_phylogeny_atlas.tsv",
    "results/alphafold/plddt_scores.tsv",
    "results/delta_llps/delta_llps_scores.tsv",
    "results/mutants/idr_llps_mutants.tsv",
    "results/mutants/pathway_differential_llps.tsv",
    "results/mutants/pathway_critical_analysis.md",
    "results/cdr/cdr_summary.tsv",
    "results/slim/slim_hits.tsv",
    "results/slim/slim_per_residue.tsv",
    "results/slim/ptm_hits.tsv",
    "results/slim/ptm_per_residue.tsv",
    "results/slim/ptm_category_summary.tsv",
    "results/map/chromosome_atlas.html",
    "results/index.html"

rule download_clinvar:
  output:
    "data/clinvar/variant_summary.txt.gz"
  shell:
    """
    mkdir -p data/clinvar
    wget -q -O {output} https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
    """

rule download_uniprot:
  output:
    "data/fasta/selected_proteins.fasta"
  shell:
    """
    mkdir -p data/fasta
    python3 scripts/download_human_gene_fasta.py --config config.yaml --output {output}
    """

rule download_orthologs:
  input:
    fasta="data/fasta/selected_proteins.fasta"
  output:
    expand("data/orthologs/{gene}.fasta", gene=config["genes"])
  run:
    import os, subprocess
    os.makedirs("data/orthologs", exist_ok=True)
    subprocess.run([
      "python3",
      "scripts/download_uniprot_orthologs.py",
      "--config", "config.yaml",
      "--outdir", "data/orthologs"
    ], check=True)
    missing = [str(path) for path in output if not os.path.exists(str(path))]
    if missing:
      raise RuntimeError("Missing ortholog outputs: " + ", ".join(missing))

rule run_iupred:
  input:
    "data/fasta/selected_proteins.fasta"
  output:
    "results/iupred/iupred_scores.tsv"
  shell:
    """
    mkdir -p results/iupred
    python3 scripts/run_iupred_batch.py {input} {output}
    """

rule run_seg:
  input:
    "data/fasta/selected_proteins.fasta"
  output:
    "results/seg/seg_regions.tsv"
  shell:
    """
    mkdir -p results/seg
    python3 scripts/run_seg_batch.py {input} {output}
    """

rule run_plaac:
  input:
    "data/fasta/selected_proteins.fasta"
  output:
    "results/plaac/plaac_scores.tsv"
  shell:
    """
    mkdir -p results/plaac
    python3 scripts/run_plaac_batch.py {input} {output}
    """

rule run_llps:
  input:
    "data/fasta/selected_proteins.fasta"
  output:
    "results/llps/llps_scores.tsv"
  shell:
    """
    mkdir -p results/llps
    python3 scripts/run_llps_proxy.py {input} {output}
    """

rule run_structural_proxy:
  input:
    "data/fasta/selected_proteins.fasta"
  output:
    "results/structure/structure_scores.tsv"
  shell:
    """
    mkdir -p results/structure
    python3 scripts/run_structure_proxy.py {input} {output}
    """

rule map_clinvar:
  input:
    clinvar="data/clinvar/variant_summary.txt.gz",
    fasta="data/fasta/selected_proteins.fasta"
  output:
    "results/clinvar/mapped_variants.tsv"
  shell:
    """
    mkdir -p results/clinvar
    python3 scripts/map_clinvar_to_proteins.py --clinvar {input.clinvar} --fasta {input.fasta} > {output}
    """

rule compute_conservation:
  input:
    expand("data/orthologs/{gene}.fasta", gene=config["genes"])
  output:
    "results/phylogeny/conservation.tsv"
  shell:
    """
    mkdir -p results/phylogeny
    python3 scripts/compute_phylo_conservation.py --ortholog-dir data/orthologs --output {output}
    """

rule build_atlas:
  input:
    iupred="results/iupred/iupred_scores.tsv",
    seg="results/seg/seg_regions.tsv",
    plaac="results/plaac/plaac_scores.tsv",
    llps="results/llps/llps_scores.tsv",
    structure="results/structure/structure_scores.tsv",
    variants="results/clinvar/mapped_variants.tsv",
    conservation="results/phylogeny/conservation.tsv",
    plddt="results/alphafold/plddt_scores.tsv",
    delta_llps="results/delta_llps/delta_llps_scores.tsv",
    cdr="results/cdr/all_genes_cdr.tsv",
    slim="results/slim/slim_per_residue.tsv",
    ptm="results/slim/ptm_per_residue.tsv",
    oncovirus="results/oncovirus/oncovirus_per_residue.tsv"
  output:
    expand("results/atlas/{gene}_atlas.tsv", gene=config["genes"]),
    "results/atlas/global_disorder_phylogeny_atlas.tsv"
  shell:
    """
    mkdir -p results/atlas
    python3 scripts/build_disorder_atlas.py \
      --iupred {input.iupred} \
      --seg {input.seg} \
      --plaac {input.plaac} \
      --llps {input.llps} \
      --structure {input.structure} \
      --variants {input.variants} \
      --conservation {input.conservation} \
      --plddt {input.plddt} \
      --delta-llps {input.delta_llps} \
      --cdr {input.cdr} \
      --slim {input.slim} \
      --ptm {input.ptm} \
      --oncovirus {input.oncovirus} \
      --fasta data/fasta/selected_proteins.fasta \
      --genes-file config.yaml \
      --outdir results/atlas
    """

rule build_atlas_with_clinvar:
  input:
    atlas="results/atlas/{gene}_atlas.tsv",
    variants="results/clinvar/mapped_variants.tsv"
  output:
    "results/atlas/{gene}_atlas_with_clinvar.tsv"
  shell:
    """
    mkdir -p results/atlas
    python3 scripts/build_clinvar_view.py \
      --atlas {input.atlas} \
      --variants {input.variants} \
      --outdir results/atlas
    """

rule render_map:
  input:
    summary="results/atlas/global_disorder_phylogeny_atlas.tsv",
    clinvar=expand("results/atlas/{gene}_atlas_with_clinvar.tsv", gene=config["genes"])
  output:
    "results/map/chromosome_atlas.html"
  shell:
    """
    mkdir -p results/map
    python3 scripts/render_chromosome_atlas.py \
      --atlas-dir results/atlas \
      --global-summary {input.summary} \
      --output {output}
    """

rule render_dashboard:
  input:
    summary="results/atlas/global_disorder_phylogeny_atlas.tsv",
    map_html="results/map/chromosome_atlas.html"
  output:
    "results/index.html"
  shell:
    """
    mkdir -p results
    python3 scripts/render_dashboard.py \
      --atlas-dir results/atlas \
      --global-summary {input.summary} \
      --output {output}
    """

# ── Novelty Layer 1: AlphaFold2 pLDDT ──────────────────────────────────────
rule download_alphafold:
  output:
    "results/alphafold/plddt_scores.tsv"
  shell:
    """
    mkdir -p results/alphafold
    python3 scripts/download_alphafold_plddt.py --output {output}
    """

# ── Novelty Layer 2: Delta-LLPS vulnerability map ──────────────────────────
rule compute_delta_llps:
  input:
    fasta="data/fasta/selected_proteins.fasta"
  output:
    "results/delta_llps/delta_llps_scores.tsv"
  shell:
    """
    mkdir -p results/delta_llps
    python3 scripts/compute_delta_llps.py --fasta {input.fasta} --output {output}
    """

rule analyze_idr_mutants_pathways:
  input:
    fasta="data/fasta/selected_proteins.fasta",
    config="config.yaml",
    atlases=expand("results/atlas/{gene}_atlas.tsv", gene=config["genes"])
  output:
    "results/mutants/idr_llps_mutants.tsv",
    "results/mutants/pathway_differential_llps.tsv",
    "results/mutants/pathway_ptm_differential_llps.tsv",
    "results/mutants/pathway_critical_analysis.md"
  shell:
    """
    mkdir -p results/mutants
    python3 scripts/analyze_idr_mutant_pathway_effects.py \
      --fasta {input.fasta} \
      --atlas-dir results/atlas \
      --config {input.config} \
      --output-dir results/mutants
    """

# ── Novelty Layer 3: CDR + AlphaFold2 discordant regions ───────────────────
rule compute_cdr:
  input:
    atlases=expand("results/atlas/{gene}_atlas.tsv", gene=config["genes"]),
    plddt="results/alphafold/plddt_scores.tsv"
  output:
    "results/cdr/cdr_summary.tsv"
  shell:
    """
    mkdir -p results/cdr
    python3 scripts/compute_cdr_novelty.py \
      --atlas-dir results/atlas \
      --plddt {input.plddt} \
      --output-dir results/cdr
    """

# ── Novelty Layer 4: SLiM scanning ─────────────────────────────────────────
rule scan_slims:
  input:
    fasta="data/fasta/selected_proteins.fasta",
    iupred="results/iupred/iupred_scores.tsv",
    variants="results/clinvar/mapped_variants.tsv"
  output:
    "results/slim/slim_hits.tsv",
    "results/slim/slim_per_residue.tsv",
    "results/slim/ptm_hits.tsv",
    "results/slim/ptm_per_residue.tsv",
    "results/slim/ptm_category_summary.tsv"
  shell:
    """
    mkdir -p results/slim
    python3 scripts/scan_slim_motifs.py \
      --fasta {input.fasta} \
      --iupred {input.iupred} \
      --variants {input.variants} \
      --output-dir results/slim
    """
