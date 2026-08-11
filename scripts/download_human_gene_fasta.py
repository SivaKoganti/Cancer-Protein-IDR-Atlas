#!/usr/bin/env python3
import argparse
import sys
import time
from io import StringIO
from pathlib import Path

import requests
import yaml
from Bio import SeqIO

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/stream"


def fetch_gene_fasta(gene, taxon_id=9606):
    query = f"gene_exact:{gene} AND organism_id:{taxon_id} AND reviewed:true"
    headers = {"User-Agent": "IDR-Atlas/1.0"}
    params = {
        "query": query,
        "format": "fasta",
    }
    resp = requests.get(UNIPROT_URL, params=params, headers=headers, timeout=60)
    if resp.status_code != 200:
        raise RuntimeError(f"Failed to fetch {gene}: status {resp.status_code} - {resp.text.strip()}")
    return resp.text.strip()


def parse_fasta_records(text):
    return list(SeqIO.parse(StringIO(text), "fasta"))


def main():
    parser = argparse.ArgumentParser(description="Download human gene FASTA entries from UniProt.")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--output", default="data/fasta/selected_proteins.fasta")
    parser.add_argument("--taxon-id", default=9606, type=int)
    args = parser.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    genes = cfg.get("genes", [])
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    with open(args.output, "w") as out:
        for gene in genes:
            try:
                fasta_text = fetch_gene_fasta(gene, taxon_id=args.taxon_id)
                if fasta_text:
                    records = parse_fasta_records(fasta_text)
                    if not records:
                        print(f"Warning: no FASTA records returned for {gene}", file=sys.stderr)
                        continue
                    for idx, rec in enumerate(records, start=1):
                        rec.id = gene if len(records) == 1 else f"{gene}_{idx}"
                        rec.name = ""
                        rec.description = rec.description
                        SeqIO.write(rec, out, "fasta")
                else:
                    print(f"Warning: no FASTA returned for {gene}", file=sys.stderr)
            except Exception as exc:
                print(f"Warning: failed to fetch {gene}: {exc}", file=sys.stderr)
            time.sleep(0.1)

    if Path(args.output).stat().st_size == 0:
        raise RuntimeError("No sequences were downloaded for the configured gene list.")


if __name__ == '__main__':
    main()
