#!/usr/bin/env python3
"""run_llps_proxy.py

Compute an LLPS propensity proxy using prion-like and physicochemical features.
"""
import argparse
import math
from pathlib import Path
from Bio import SeqIO
import yaml


def score_window(seq, weights):
    qn = sum(1 for c in seq if c in "QN") / len(seq)
    arom = sum(1 for c in seq if c in "FYWI") / len(seq)
    charged = sum(1 for c in seq if c in "DEKRH") / len(seq)
    hydrophobic = sum(1 for c in seq if c in "AILMVFWY") / len(seq)
    return (
        weights["qn"] * qn
        + weights["arom"] * arom
        + weights["hydrophobic"] * hydrophobic
        - weights["charge"] * charged
    )


def main():
    parser = argparse.ArgumentParser(description="Compute an LLPS proxy score per residue.")
    parser.add_argument("fasta")
    parser.add_argument("output")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    config = {
        "window": 41,
        "qn_weight": 0.45,
        "aromatic_weight": 0.25,
        "hydrophobic_weight": 0.15,
        "charge_weight": 0.15,
    }

    try:
        with open(args.config) as fh:
            user_cfg = yaml.safe_load(fh) or {}
        llps_cfg = user_cfg.get("llps", {})
        config["window"] = llps_cfg.get("window", config["window"])
        config["qn_weight"] = llps_cfg.get("qn_weight", config["qn_weight"])
        config["aromatic_weight"] = llps_cfg.get("aromatic_weight", config["aromatic_weight"])
        config["hydrophobic_weight"] = llps_cfg.get("hydrophobic_weight", config["hydrophobic_weight"])
        config["charge_weight"] = llps_cfg.get("charge_weight", llps_cfg.get("charge_pattern_weight", config["charge_weight"]))
    except FileNotFoundError:
        pass

    with open(args.output, "w") as outf:
        outf.write("gene\tstart\tend\tllps_proxy\n")
        for rec in SeqIO.parse(args.fasta, "fasta"):
            gid = rec.id.split("|")[0]
            seq = str(rec.seq)
            L = len(seq)
            if L == 0:
                continue
            w = config["window"]
            if L < w:
                windows = [(0, seq)]
            else:
                windows = [(i, seq[i:i+w]) for i in range(0, L - w + 1)]

            for i, win in windows:
                if len(win) == 0:
                    continue
                score = score_window(win, {
                    "qn": config["qn_weight"],
                    "arom": config["aromatic_weight"],
                    "hydrophobic": config["hydrophobic_weight"],
                    "charge": config["charge_weight"],
                })
                outf.write(f"{gid}\t{i+1}\t{i+len(win)}\t{score:.4f}\n")


if __name__ == '__main__':
    main()
