#!/usr/bin/env python3
"""run_structure_proxy.py

Compute structural and quantum-inspired residue scores using simplified physics rules.
"""
import argparse
from pathlib import Path
from Bio import SeqIO

AMINO_ACID_PROPERTIES = {
    "A": {"polarity": 1.8, "size": 1.0, "pi": 0.0},
    "R": {"polarity": 10.5, "size": 2.0, "pi": 1.0},
    "N": {"polarity": 11.6, "size": 1.5, "pi": 0.5},
    "D": {"polarity": 13.0, "size": 1.5, "pi": 0.5},
    "C": {"polarity": 5.5, "size": 1.7, "pi": 1.0},
    "Q": {"polarity": 10.0, "size": 1.8, "pi": 0.5},
    "E": {"polarity": 12.3, "size": 1.8, "pi": 0.5},
    "G": {"polarity": 9.0, "size": 0.8, "pi": 0.0},
    "H": {"polarity": 10.4, "size": 2.0, "pi": 0.8},
    "I": {"polarity": 5.2, "size": 2.0, "pi": 0.0},
    "L": {"polarity": 4.9, "size": 2.0, "pi": 0.0},
    "K": {"polarity": 11.3, "size": 2.0, "pi": 1.0},
    "M": {"polarity": 5.7, "size": 2.0, "pi": 0.5},
    "F": {"polarity": 5.2, "size": 2.2, "pi": 1.5},
    "P": {"polarity": 8.0, "size": 1.9, "pi": 0.0},
    "S": {"polarity": 9.2, "size": 1.4, "pi": 0.0},
    "T": {"polarity": 8.6, "size": 1.5, "pi": 0.0},
    "W": {"polarity": 5.4, "size": 2.5, "pi": 2.0},
    "Y": {"polarity": 6.2, "size": 2.3, "pi": 1.5},
    "V": {"polarity": 5.9, "size": 1.9, "pi": 0.0},
}


def quantum_proxy(aa, config):
    properties = AMINO_ACID_PROPERTIES.get(aa, {"polarity": 9.0, "size": 1.5, "pi": 0.0})
    dielectric = config["dielectric_constant"]
    polarizability = properties["polarity"] / dielectric
    stacking = properties["pi"] * config["stacking_coefficient"]
    charge_penalty = config["charge_penalty"] if aa in "DEKRH" else 0.0
    return max(0.0, polarizability + stacking - charge_penalty)


def main():
    parser = argparse.ArgumentParser(description="Compute structural provenance and quantum-inspired scores per residue.")
    parser.add_argument("fasta")
    parser.add_argument("output")
    args = parser.parse_args()

    config = {
        "dielectric_constant": 4.0,
        "polarizability_factor": 1.0,
        "stacking_coefficient": 1.2,
        "charge_penalty": 0.6,
    }

    with open(args.output, "w") as outf:
        outf.write("gene\tpos\taa\tstructural_proxy\tquantum_proxy\n")
        for rec in SeqIO.parse(args.fasta, "fasta"):
            gid = rec.id.split("|")[0]
            seq = str(rec.seq)
            for pos, aa in enumerate(seq, start=1):
                qp = quantum_proxy(aa, config)
                properties = AMINO_ACID_PROPERTIES.get(aa, {"polarity": 9.0, "size": 1.5, "pi": 0.0})
                spu = 1.0 / (1.0 + properties["size"])
                outf.write(f"{gid}\t{pos}\t{aa}\t{spu:.4f}\t{qp:.4f}\n")


if __name__ == '__main__':
    main()
