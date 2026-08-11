#!/usr/bin/env python3
"""map_clinvar_to_proteins.py

Simple mapper that parses ClinVar variant_summary and matches by GeneSymbol and ProteinChange (if present).
Outputs: gene\tprotein_pos\tref\talt\tclinical_significance\t

This is a best-effort lightweight mapper for the MVP; for production use, use VEP/Ensembl mapping for canonical transcripts.
"""
import gzip
import argparse
import re


def normalize_gene_id(gene):
    value = str(gene).strip()
    if not value:
        return ""
    return re.sub(r"_\d+$", "", value)


def extract_protein_change(text):
    """Extract (ref, pos, alt) from free-text protein-change fields.

    Accepts forms like:
      p.Arg248Trp
      p.R248W
      p.Trp24*
      p.Gly12=
      p.Arg337fs
    """
    if text is None:
        return None
    token = str(text).strip()
    if not token:
        return None

    # Prefer explicit parenthesized protein change in ClinVar Name, e.g. (...(p.Arg248Trp))
    paren_match = re.search(r"\(p\.([A-Za-z*]{1,3})(\d+)([A-Za-z*=]{1,10})\)", token)
    if paren_match:
        return paren_match.groups()

    direct_match = re.search(r"p\.([A-Za-z*]{1,3})(\d+)([A-Za-z*=]{1,10})", token)
    if direct_match:
        return direct_match.groups()

    return None


def row_gene_symbols(raw_gene_field):
    if raw_gene_field is None:
        return []
    raw = str(raw_gene_field).strip()
    if not raw:
        return []
    return [normalize_gene_id(part) for part in re.split(r"[;,]", raw) if part.strip()]

p = argparse.ArgumentParser()
p.add_argument('--clinvar', required=True)
p.add_argument('--fasta', required=True)
args = p.parse_args()

# load gene list from fasta headers
genes = set()
with open(args.fasta) as fh:
    for line in fh:
        if line.startswith('>'):
            header = line[1:].strip()
            gid = header.split()[0]
            genes.add(normalize_gene_id(gid))

# parse clinvar summary
infile = args.clinvar
openf = gzip.open if infile.endswith('.gz') else open
with openf(infile, 'rt') as fh:
    hdr = next(fh).strip().split('\t')
    idx = {k:i for i,k in enumerate(hdr)}
    out = []
    protein_change_key = None
    if 'ProteinChange' in idx:
        protein_change_key = 'ProteinChange'
    elif 'Name' in idx:
        protein_change_key = 'Name'
    for line in fh:
        cols = line.strip().split('\t')
        raw_gene = cols[idx.get('GeneSymbol', '-')] if 'GeneSymbol' in idx else ''
        row_genes = row_gene_symbols(raw_gene)
        target_genes = [g for g in row_genes if g in genes]
        if not target_genes:
            continue

        prot_raw = cols[idx.get(protein_change_key, '-')] if protein_change_key else '-'
        protein_change = extract_protein_change(prot_raw)
        if protein_change is None and 'Name' in idx:
            protein_change = extract_protein_change(cols[idx['Name']])
        if protein_change is None:
            continue

        ref, pos, alt = protein_change
        cs = cols[idx.get('ClinicalSignificance','-')] if 'ClinicalSignificance' in idx else ''
        for gene in target_genes:
            out.append((gene, pos, ref, alt, cs))

# write to stdout
for row in out:
    print('\t'.join(map(str,row)))
