# Snakemake workflow skeleton for MVP

configfile: "config.yaml"

rule all:
  input:
    expand("results/{{gene}}/final_table.tsv", gene=[g for g in config["genes"]])

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
    "data/uniprot/human_proteome.fasta.gz"
  shell:
    """
    mkdir -p data/uniprot
    wget -q -O {output} https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/proteomes/UP000005640_9606.fasta.gz
    """

rule extract_proteins:
  input:
    fasta="data/uniprot/human_proteome.fasta.gz",
    genes="config.yaml"
  output:
    "data/fasta/selected_proteins.fasta"
  run:
    import gzip
    genes = set(config["genes"])
    outpath = output[0]
    with gzip.open(input.fasta, "rt") as fh, open(outpath, "w") as out:
      write = False
      header = None
      seq = []
      for line in fh:
        if line.startswith(">"):
          if header and seq and any(g in header for g in genes):
            out.write(header)
            out.write(''.join(seq))
          header = line
          seq = []
        else:
          seq.append(line)
      if header and seq and any(g in header for g in genes):
        out.write(header)
        out.write(''.join(seq))

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

rule aggregate:
  input:
    iupred="results/iupred/iupred_scores.tsv",
    seg="results/seg/seg_regions.tsv",
    plaac="results/plaac/plaac_scores.tsv",
    variants="results/clinvar/mapped_variants.tsv"
  output:
    "results/FINAL_AGGREGATED_PER_PROTEIN.tar.gz"
  shell:
    """
    mkdir -p results/aggregate
    python3 scripts/aggregate_annotations.py --iupred {input.iupred} --seg {input.seg} --plaac {input.plaac} --variants {input.variants} --outdir results/aggregate
    tar -czf {output} results/aggregate
    """
