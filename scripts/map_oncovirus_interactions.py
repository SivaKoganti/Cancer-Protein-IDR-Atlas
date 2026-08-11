#!/usr/bin/env python3
"""
Map curated oncovirus-host interaction regions to per-residue flags.

Input TSV columns:
  gene, virus, start, end, evidence, source

Output TSV columns:
  gene, pos, virus_interaction, virus_count, virus_names, virus_evidence
"""

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO
import yaml


def normalise_gene_id(gene_id: str) -> str:
    token = gene_id.strip()
    if "_" in token:
        left, right = token.rsplit("_", 1)
        if right.isdigit():
            return left
    return token


def load_gene_list(config_path: Path | None) -> set[str] | None:
    if config_path is None or not config_path.exists():
        return None
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)
    return set(cfg.get("genes", []))


def load_lengths(fasta_path: Path, allowed_genes: set[str] | None = None) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        gid = normalise_gene_id(rec.id.split("|")[0].strip())
        if allowed_genes is not None and gid not in allowed_genes:
            continue
        seq_len = len(str(rec.seq))
        if gid not in lengths or seq_len > lengths[gid]:
            lengths[gid] = seq_len
    return lengths


def main():
    parser = argparse.ArgumentParser(description="Map oncovirus interactions to per-residue table.")
    parser.add_argument("--fasta", required=True, help="FASTA with selected proteins")
    parser.add_argument("--interactions", required=True, help="Curated interaction region TSV")
    parser.add_argument("--output", required=True, help="Output per-residue TSV")
    parser.add_argument("--genes-file", default=None, help="Optional config.yaml with genes list")
    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    interactions_path = Path(args.interactions)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    allowed_genes = load_gene_list(Path(args.genes_file)) if args.genes_file else None
    lengths = load_lengths(fasta_path, allowed_genes=allowed_genes)
    interactions = pd.read_csv(interactions_path, sep="\t")
    interactions["gene"] = interactions["gene"].astype(str).map(normalise_gene_id)

    rows = []
    for gene, length in sorted(lengths.items()):
        gene_hits = interactions[interactions["gene"].astype(str).str.strip() == gene]

        by_pos = {
            pos: {"viruses": set(), "evidence": set()}
            for pos in range(1, length + 1)
        }

        for _, row in gene_hits.iterrows():
            try:
                start = int(row["start"])
                end = int(row["end"])
            except (TypeError, ValueError):
                continue

            if end < start:
                start, end = end, start

            start = max(1, start)
            end = min(length, end)
            if start > end:
                continue

            virus = str(row.get("virus", "")).strip()
            evidence = str(row.get("evidence", "")).strip()
            source = str(row.get("source", "")).strip()
            evidence_label = evidence if not source else f"{evidence} [{source}]".strip()

            for pos in range(start, end + 1):
                if virus:
                    by_pos[pos]["viruses"].add(virus)
                if evidence_label:
                    by_pos[pos]["evidence"].add(evidence_label)

        for pos in range(1, length + 1):
            viruses = sorted(by_pos[pos]["viruses"])
            evidence_labels = sorted(by_pos[pos]["evidence"])
            rows.append(
                {
                    "gene": gene,
                    "pos": pos,
                    "virus_interaction": 1 if viruses else 0,
                    "virus_count": len(viruses),
                    "virus_names": ";".join(viruses),
                    "virus_evidence": " | ".join(evidence_labels),
                }
            )

    out_df = pd.DataFrame(rows)
    out_df.to_csv(output_path, sep="\t", index=False)

    interacting_genes = out_df.groupby("gene")["virus_interaction"].max()
    print(f"Wrote {len(out_df)} residue rows to {output_path}")
    print(f"Genes with at least one mapped virus interaction: {int(interacting_genes.sum())} / {len(interacting_genes)}")


if __name__ == "__main__":
    main()
