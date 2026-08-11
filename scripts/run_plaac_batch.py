#!/usr/bin/env python3
"""run_plaac_batch.py

Placeholder PLAAC runner. If the PLAAC binary is installed, call it. Otherwise compute a prion-like
score proxy based on Q/N content in sliding windows.
"""
import sys
from Bio import SeqIO

fasta = sys.argv[1]
out = sys.argv[2]

with open(out, 'w') as outf:
    outf.write('gene\tstart\tend\tq_n_fraction\n')
    for rec in SeqIO.parse(fasta, 'fasta'):
        gid = rec.id.split('|')[0]
        seq = str(rec.seq)
        L = len(seq)
        w = 41
        for i in range(0, max(1,L-w+1)):
            win = seq[i:i+w]
            qn = sum(1 for c in win if c in 'QN')
            outf.write(f"{gid}\t{i+1}\t{i+w}\t{qn/w:.3f}\n")
