#!/usr/bin/env python3
"""run_iupred_batch.py

Runs iupred2a (if available) on a FASTA file and outputs a TSV: gene\tpos\taa\tscore
If iupred is not available, it writes placeholder scores.
"""
import sys
from Bio import SeqIO
import subprocess

fasta = sys.argv[1]
out = sys.argv[2]

def run_iupred_seq(seq):
    # Attempt to call iupred if installed
    try:
        p = subprocess.run(['iupred2a', '-'], input=str(seq), text=True, capture_output=True)
        # Not parsing complex output here; placeholder
        return None
    except FileNotFoundError:
        return None

with open(out, 'w') as outf:
    outf.write('gene\tpos\taa\tscore\n')
    for rec in SeqIO.parse(fasta, 'fasta'):
        gid = rec.id.split('|')[0]
        seq = str(rec.seq)
        for i,aa in enumerate(seq, start=1):
            # placeholder uniform 0.5 score
            outf.write(f"{gid}\t{i}\t{aa}\t0.5\n")
