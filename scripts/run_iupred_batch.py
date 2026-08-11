#!/usr/bin/env python3
"""Run IUPred if available, otherwise emit a deterministic disorder heuristic."""

import shutil
import subprocess
import sys

from Bio import SeqIO

fasta = sys.argv[1]
out = sys.argv[2]

DISORDER_PROPENSITY = {
    "A": 0.38, "C": 0.35, "D": 0.72, "E": 0.74, "F": 0.27,
    "G": 0.62, "H": 0.56, "I": 0.23, "K": 0.71, "L": 0.24,
    "M": 0.31, "N": 0.64, "P": 0.78, "Q": 0.68, "R": 0.66,
    "S": 0.65, "T": 0.53, "V": 0.25, "W": 0.22, "Y": 0.33,
}
WINDOW = 21


def parse_iupred_output(stdout, expected_len):
    scores = []
    for line in stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            int(parts[0])
            score = float(parts[2])
        except ValueError:
            continue
        scores.append(score)
    return scores if len(scores) == expected_len else None


def heuristic_scores(seq):
    half = WINDOW // 2
    scores = []
    for idx, aa in enumerate(seq):
        start = max(0, idx - half)
        end = min(len(seq), idx + half + 1)
        window = seq[start:end]
        base = DISORDER_PROPENSITY.get(aa, 0.5)
        mean_prop = sum(DISORDER_PROPENSITY.get(res, 0.5) for res in window) / len(window)
        charged = sum(1 for res in window if res in "DEKR") / len(window)
        aromatic = sum(1 for res in window if res in "FWYIAVLMC") / len(window)
        gly_pro = sum(1 for res in window if res in "GPQSENKRD") / len(window)
        score = 0.45 * base + 0.35 * mean_prop + 0.15 * gly_pro + 0.10 * charged - 0.20 * aromatic
        scores.append(min(0.95, max(0.05, round(score, 4))))
    return scores

def run_iupred_seq(seq):
    for executable in ("iupred2a", "iupred"):
        if not shutil.which(executable):
            continue
        try:
            proc = subprocess.run(
                [executable, "-"],
                input=str(seq),
                text=True,
                capture_output=True,
                check=False,
            )
        except OSError:
            continue
        if proc.returncode == 0:
            scores = parse_iupred_output(proc.stdout, len(seq))
            if scores is not None:
                return scores
    return heuristic_scores(seq)

with open(out, 'w') as outf:
    outf.write('gene\tpos\taa\tscore\n')
    for rec in SeqIO.parse(fasta, 'fasta'):
        gid = rec.id.split('|')[0]
        seq = str(rec.seq)
        scores = run_iupred_seq(seq)
        for i, aa in enumerate(seq, start=1):
            outf.write(f"{gid}\t{i}\t{aa}\t{scores[i-1]:.4f}\n")
