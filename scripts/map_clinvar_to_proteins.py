#!/usr/bin/env python3
"""map_clinvar_to_proteins.py

Simple mapper that parses ClinVar variant_summary and matches by GeneSymbol and ProteinChange (if present).
Outputs: gene\tprotein_pos\tref\talt\tclinical_significance\t

This is a best-effort lightweight mapper for the MVP; for production use, use VEP/Ensembl mapping for canonical transcripts.
"""
import gzip
import argparse
import re
from collections import defaultdict

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
            genes.add(gid)

# parse clinvar summary
infile = args.clinvar
openf = gzip.open if infile.endswith('.gz') else open
with openf(infile, 'rt') as fh:
    hdr = next(fh).strip().split('\t')
    idx = {k:i for i,k in enumerate(hdr)}
    out = []
    for line in fh:
        cols = line.strip().split('\t')
        gene = cols[idx.get('GeneSymbol','-')]
        if gene in genes:
            prot = cols[idx.get('ProteinChange','-')]
            if prot and prot.startswith('p.'):
                m = re.match(r'p\.([A-Za-z]+)(\d+)([A-Za-z*]+)', prot)
                if m:
                    ref, pos, alt = m.groups()
                    cs = cols[idx.get('ClinicalSignificance','-')]
                    out.append((gene, pos, ref, alt, cs))

# write to stdout
for row in out:
    print('\t'.join(map(str,row)))
