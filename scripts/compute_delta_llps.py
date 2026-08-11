#!/usr/bin/env python3
"""
Compute delta-LLPS scores: how each amino acid substitution changes phase-separation propensity.

Produces three outputs per residue:
  llps_vulnerability  — mean |delta_LLPS| across all 19 possible substitutions (positional sensitivity)
  delta_llps_max      — largest gain in LLPS propensity achievable by any single substitution
  delta_llps_min      — largest loss (most negative delta) achievable
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

# Amino acid group memberships used in LLPS proxy (mirrors run_llps_proxy.py)
QN = set("QN")
AROMATIC = set("FYWI")
HYDROPHOBIC = set("LVAM")
CHARGED = set("RKDE")

WINDOW = 41      # must match config llps window
QN_W = 0.45
AR_W = 0.25
HY_W = 0.15
CH_W = 0.15

ALL_AA = list("ACDEFGHIKLMNPQRSTVWY")


def llps_window_score(seq: str, center: int, window: int = WINDOW) -> float:
    half = window // 2
    start = max(0, center - half)
    end = min(len(seq), center + half + 1)
    sub = seq[start:end]
    n = len(sub)
    if n == 0:
        return 0.0
    qn = sum(1 for a in sub if a in QN) / n
    ar = sum(1 for a in sub if a in AROMATIC) / n
    hy = sum(1 for a in sub if a in HYDROPHOBIC) / n
    ch = sum(1 for a in sub if a in CHARGED) / n
    return QN_W * qn + AR_W * ar + HY_W * hy - CH_W * ch


def compute_delta_llps_profile(seq: str) -> pd.DataFrame:
    """Return per-residue delta-LLPS statistics over all possible single substitutions."""
    seq = seq.upper()
    rows = []
    wt_scores = [llps_window_score(seq, i) for i in range(len(seq))]

    for i, wt_aa in enumerate(seq):
        wt_s = wt_scores[i]
        deltas = []
        for alt_aa in ALL_AA:
            if alt_aa == wt_aa:
                continue
            mut_seq = seq[:i] + alt_aa + seq[i+1:]
            mut_s = llps_window_score(mut_seq, i)
            deltas.append(mut_s - wt_s)
        rows.append({
            "pos": i + 1,
            "aa": wt_aa,
            "wt_llps": round(wt_s, 4),
            "llps_vulnerability": round(float(np.mean(np.abs(deltas))), 4),
            "delta_llps_max": round(float(max(deltas)), 4),   # max gain
            "delta_llps_min": round(float(min(deltas)), 4),   # max loss
        })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description="Compute delta-LLPS vulnerability per residue.")
    parser.add_argument("--fasta", required=True, help="Selected proteins FASTA")
    parser.add_argument("--output", default="results/delta_llps/delta_llps_scores.tsv")
    args = parser.parse_args()

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    frames = []
    for rec in SeqIO.parse(args.fasta, "fasta"):
        gid = rec.id.split("|")[0].strip()
        # Normalise isoform names (e.g. CDKN2A_1 → CDKN2A)
        base = gid.rsplit("_", 1)[0] if "_" in gid and gid.rsplit("_", 1)[1].isdigit() else gid
        seq = str(rec.seq).upper()
        print(f"  {base}: {len(seq)} aa ...", end=" ", flush=True)
        df = compute_delta_llps_profile(seq)
        df.insert(0, "gene", base)
        frames.append(df)
        print("done")

    result = pd.concat(frames, ignore_index=True)
    result.to_csv(out_path, sep="\t", index=False)
    print(f"\n✓ Saved delta-LLPS scores → {out_path} ({len(result)} residues)")


if __name__ == "__main__":
    main()
