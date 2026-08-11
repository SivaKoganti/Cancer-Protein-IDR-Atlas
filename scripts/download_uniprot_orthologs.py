#!/usr/bin/env python3
import argparse
import os
import re
import sys
import time
from io import StringIO
from pathlib import Path

import requests
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/stream"


def get_gene_symbol_from_record(rec):
    m = re.search(r"GN=([A-Za-z0-9_-]+)", rec.description)
    if m:
        return m.group(1)
    if "|" in rec.id:
        parts = rec.id.split("|")
        if len(parts) >= 3:
            return parts[2].split("_")[0]
    return rec.id


def parse_fasta_records(text):
    if not text:
        return []
    return list(SeqIO.parse(StringIO(text), "fasta"))


def fetch_fasta(gene, taxon_id, exact=True):
    if exact:
        query = f"gene_exact:{gene} AND organism_id:{taxon_id} AND reviewed:true"
    else:
        query = f"gene:{gene} AND organism_id:{taxon_id} AND reviewed:true"
    headers = {"User-Agent": "IDR-Atlas/1.0"}
    params = {
        "query": query,
        "format": "fasta",
    }
    resp = requests.get(UNIPROT_URL, params=params, headers=headers, timeout=30)
    if resp.status_code != 200:
        raise RuntimeError(f"Uniprot fetch failed {resp.status_code}: {resp.text}")
    return parse_fasta_records(resp.text)


def load_fasta_entry(gene, fasta_path):
    try:
        for rec in SeqIO.parse(str(fasta_path), "fasta"):
            if rec.id == gene or rec.id.startswith(f"{gene}_") or gene in rec.description:
                return rec
    except FileNotFoundError:
        return None
    return None


def main():
    parser = argparse.ArgumentParser(description="Download ortholog FASTA sets for cancer genes.")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--outdir", default="data/orthologs")
    args = parser.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    species = cfg.get("phylogeny", {}).get("species", [])
    genes = cfg.get("genes", [])
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    seed_fasta = Path("data/fasta/selected_proteins.fasta")

    for gene in genes:
        out_path = outdir / f"{gene}.fasta"
        records = []
        for spec in species:
            taxon_id = spec.get("taxon_id")
            name = spec.get("name")
            try:
                fetched_records = fetch_fasta(gene, taxon_id)
                if not fetched_records:
                    fetched_records = fetch_fasta(gene, taxon_id, exact=False)
                    if fetched_records:
                        print(
                            f"Warning: gene_exact search returned no results; using gene search for {gene} {name} ({taxon_id})",
                            file=sys.stderr,
                        )
                if fetched_records:
                    for idx, rec in enumerate(fetched_records, start=1):
                        gene_name = get_gene_symbol_from_record(rec)
                        rec.id = f"{gene_name}_{name}_{idx}"
                        rec.name = ""
                        records.append(rec)
                    time.sleep(0.2)
                else:
                    print(
                        f"Warning: no sequence returned for {gene} {name} ({taxon_id}); skipping.",
                        file=sys.stderr,
                    )
            except Exception as exc:
                print(f"Warning: failed to fetch {gene} {name} ({taxon_id}): {exc}", file=sys.stderr)

        if not records and seed_fasta.exists():
            print(f"Warning: no orthologs downloaded for {gene}; falling back to human sequence.", file=sys.stderr)
            fallback = load_fasta_entry(gene, seed_fasta)
            if fallback:
                fallback.id = f"{gene}_human_fallback"
                fallback.name = ""
                records.append(fallback)

        if not records:
            print(f"Warning: no ortholog information available for {gene}; writing placeholder entry.", file=sys.stderr)
            placeholder = SeqRecord(Seq("X"), id=f"{gene}_placeholder", description="placeholder")
            records.append(placeholder)

        SeqIO.write(records, out_path, "fasta")


if __name__ == '__main__':
    main()
