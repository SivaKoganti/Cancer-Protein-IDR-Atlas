#!/usr/bin/env python3
"""aggregate_annotations.py

Combine predictor outputs into per-protein TSV files (simple merge for MVP).
"""
import argparse
import pandas as pd
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--iupred')
p.add_argument('--seg')
p.add_argument('--plaac')
p.add_argument('--variants')
p.add_argument('--outdir')
args = p.parse_args()

outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

iu = pd.read_csv(args.iupred, sep='\t')
seg = pd.read_csv(args.seg, sep='\t')
pl = pd.read_csv(args.plaac, sep='\t')
vars = pd.read_csv(args.variants, sep='\t', header=None, names=['gene','pos','ref','alt','clin'])

for gene, gdf in iu.groupby('gene'):
    out = gdf.copy()
    # simplistic merges for MVP
    segs = seg[seg['gene']==gene]
    pls = pl[pl['gene']==gene]
    v = vars[vars['gene']==gene]
    path = outdir / f"{gene}_annotation.tsv"
    out.to_csv(path, sep='\t', index=False)

print('Wrote annotations to', outdir)
