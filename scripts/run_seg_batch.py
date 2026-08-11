#!/usr/bin/env python3
"""run_seg_batch.py

Runs SEG for low complexity detection. Placeholder implementation that flags low-complexity
windows with simple Shannon entropy if SEG isn't available.
"""
import sys
from Bio import SeqIO
import math

fasta = sys.argv[1]
out = sys.argv[2]

def shannon_entropy(s):
    from collections import Counter
    L = len(s)
    freqs = Counter(s)
    H = 0.0
    for v in freqs.values():
        p = v/L
        H -= p * math.log2(p)
    return H

with open(out, 'w') as outf:
    outf.write('gene\tstart\tend\tscore\n')
    for rec in SeqIO.parse(fasta, 'fasta'):
        gid = rec.id.split('|')[0]
        seq = str(rec.seq)
        L = len(seq)
        w = 12
        for i in range(0, max(1,L-w+1)):
            win = seq[i:i+w]
            H = shannon_entropy(win)
            if H < 3.0:
                outf.write(f"{gid}\t{i+1}\t{i+w}\t{H:.3f}\n")
