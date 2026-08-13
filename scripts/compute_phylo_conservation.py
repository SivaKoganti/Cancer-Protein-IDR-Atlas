#!/usr/bin/env python3
import argparse
from pathlib import Path
from collections import Counter

from Bio import SeqIO
from Bio import pairwise2


def align_reference(reference, seq):
    aln = pairwise2.align.globalxx(reference, seq, one_alignment_only=True)[0]
    return aln.seqA, aln.seqB


def conservation_from_alignments(reference, aligned_seqs):
    ref_len = len(reference)
    scores = [0.0] * ref_len
    count = [0] * ref_len
    for aligned in aligned_seqs:
        for i, (a, b) in enumerate(zip(reference, aligned)):
            if a == '-':
                continue
            if b == '-':
                continue
            count[i] += 1
            if a == b:
                scores[i] += 1
    result = []
    for i in range(ref_len):
        if count[i] > 0:
            result.append(scores[i] / count[i])
        else:
            result.append(0.0)
    return result


def parse_ortholog_file(path):
    seqs = [rec for rec in SeqIO.parse(str(path), "fasta")]
    return seqs


def main():
    parser = argparse.ArgumentParser(description="Compute simple conservation scores from ortholog alignments.")
    parser.add_argument("--ortholog-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    ortholog_dir = Path(args.ortholog_dir)
    out_path = Path(args.output)
    rows = []

    for fasta_path in sorted(ortholog_dir.glob("*.fasta")):
        seqs = parse_ortholog_file(fasta_path)
        if not seqs:
            continue
        reference = str(seqs[0].seq)
        aligned_sets = []
        for seq in seqs[1:]:
            aligned_ref, aligned_seq = align_reference(reference, str(seq.seq))
            aligned_sets.append(aligned_seq)
        if not aligned_sets:
            scores = [1.0] * len(reference)
        else:
            scores = conservation_from_alignments(reference, aligned_sets)

        gene = fasta_path.stem
        for pos, (aa, score) in enumerate(zip(reference, scores), start=1):
            rows.append((gene, pos, aa, round(score, 3), len(seqs)))

    with open(out_path, "w") as out:
        out.write("gene\tpos\tref_aa\tconservation_score\tortholog_count\n")
        for row in rows:
            out.write("\t".join(map(str, row)) + "\n")


if __name__ == '__main__':
    main()
