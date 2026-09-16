#!/usr/bin/env python3
"""Compute a deterministic spin-glass-inspired residue heuristic from sequence windows.

This module intentionally provides a transparent sequence-only surrogate. It does
not claim that sequence-derived features establish physical spin-glass behavior.
"""

import argparse
import hashlib
import math
import random
import re
from pathlib import Path

import pandas as pd
import yaml
from Bio import SeqIO


METHOD = "spin_glass_inspired_sequence_heuristic"
VERSION = "1.0.0"
EPSILON = 1e-9
FEATURE_ORDER = ("charge", "polarity", "aromaticity", "hydrophobicity")
SPIN_WEIGHTS = (1.25, 1.0, 0.75, 0.75)
DEFAULT_CONFIG = {
    "enabled": True,
    "window_size": 15,
    "temperature": 1.5,
    "coupling_scale": 1.0,
    "field_scale": 0.35,
    "num_shuffles": 16,
    "encoding": "biophysical_v1",
}

RESIDUE_FEATURES = {
    "A": (0.0, 0.0, 0.0, 1.0),
    "C": (0.0, 0.0, 0.0, 1.0),
    "D": (-1.0, 1.0, 0.0, -1.0),
    "E": (-1.0, 1.0, 0.0, -1.0),
    "F": (0.0, 0.0, 1.0, 1.0),
    "G": (0.0, 0.0, 0.0, -1.0),
    "H": (0.5, 1.0, 1.0, 0.0),
    "I": (0.0, 0.0, 0.0, 1.0),
    "K": (1.0, 1.0, 0.0, -1.0),
    "L": (0.0, 0.0, 0.0, 1.0),
    "M": (0.0, 0.0, 0.0, 1.0),
    "N": (0.0, 1.0, 0.0, -1.0),
    "P": (0.0, 0.0, 0.0, -1.0),
    "Q": (0.0, 1.0, 0.0, -1.0),
    "R": (1.0, 1.0, 0.0, -1.0),
    "S": (0.0, 1.0, 0.0, -1.0),
    "T": (0.0, 1.0, 0.0, 0.0),
    "V": (0.0, 0.0, 0.0, 1.0),
    "W": (0.0, 0.0, 1.0, 1.0),
    "Y": (0.0, 1.0, 1.0, 0.0),
}


def normalize_gene_id(gene):
    if gene is None:
        return ""
    value = str(gene).strip()
    if not value:
        return ""
    return re.sub(r"_\d+$", "", value)


def load_config(config_path):
    if not config_path or not Path(config_path).exists():
        return {}
    with open(config_path, encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def load_spin_glass_config(config_path):
    cfg = dict(DEFAULT_CONFIG)
    user_cfg = load_config(config_path).get("spin_glass", {})
    if isinstance(user_cfg, dict):
        cfg.update(user_cfg)

    cfg["enabled"] = bool(cfg.get("enabled", True))
    cfg["window_size"] = max(1, int(cfg.get("window_size", DEFAULT_CONFIG["window_size"])))
    if cfg["window_size"] % 2 == 0:
        cfg["window_size"] += 1
    cfg["temperature"] = max(float(cfg.get("temperature", DEFAULT_CONFIG["temperature"])), 1e-6)
    cfg["coupling_scale"] = float(cfg.get("coupling_scale", DEFAULT_CONFIG["coupling_scale"]))
    cfg["field_scale"] = float(cfg.get("field_scale", DEFAULT_CONFIG["field_scale"]))
    cfg["num_shuffles"] = max(0, int(cfg.get("num_shuffles", DEFAULT_CONFIG["num_shuffles"])))
    cfg["encoding"] = str(cfg.get("encoding", DEFAULT_CONFIG["encoding"]))
    return cfg


def load_fasta_sequences(path):
    sequences = []
    for record in SeqIO.parse(path, "fasta"):
        gene = normalize_gene_id(record.id.split("|")[0].strip())
        sequences.append((gene, str(record.seq).upper()))
    return sequences


def residue_vector(amino_acid, encoding="biophysical_v1"):
    if encoding != "biophysical_v1":
        raise ValueError(f"Unsupported encoding: {encoding}")
    return RESIDUE_FEATURES.get(str(amino_acid).upper(), (0.0, 0.0, 0.0, 0.0))


def effective_spin(vector):
    numerator = sum(weight * value for weight, value in zip(SPIN_WEIGHTS, vector))
    denominator = sum(abs(weight) for weight in SPIN_WEIGHTS)
    return numerator / denominator if denominator else 0.0


def mean(values):
    if not values:
        return 0.0
    return sum(values) / len(values)


def variance(values):
    if len(values) < 2:
        return 0.0
    avg = mean(values)
    return sum((value - avg) ** 2 for value in values) / len(values)


def stable_seed(*parts):
    joined = "||".join(str(part) for part in parts)
    return int(hashlib.sha256(joined.encode("utf-8")).hexdigest()[:16], 16)


def safe_tanh(value):
    return math.tanh(max(-60.0, min(60.0, value)))


def safe_sigmoid(value):
    clipped = max(-60.0, min(60.0, value))
    return 1.0 / (1.0 + math.exp(-clipped))


def safe_zscore(value, reference_mean, reference_std):
    if reference_std <= EPSILON:
        return 0.0
    return (value - reference_mean) / reference_std


def summarize_null(values):
    if not values:
        return 0.0, 0.0
    return mean(values), math.sqrt(max(variance(values), 0.0))


def compute_window_metrics(residues, center_index, coupling_scale, field_scale, temperature, encoding):
    vectors = [residue_vector(residue, encoding=encoding) for residue in residues]
    spins = [effective_spin(vector) for vector in vectors]
    center_vector = vectors[center_index]
    center_spin = spins[center_index]
    local_field = field_scale * mean(spins)

    couplings = []
    signed_products = []
    effective_field = local_field

    for idx, spin in enumerate(spins):
        if idx == center_index:
            continue
        distance = abs(idx - center_index)
        coupling = coupling_scale * (
            sum(a * b for a, b in zip(center_vector, vectors[idx])) / (len(FEATURE_ORDER) * max(distance, 1))
        )
        couplings.append(coupling)
        signed_product = coupling * center_spin * spin
        signed_products.append(signed_product)
        effective_field += coupling * spin

    local_energy = -(center_spin * effective_field)
    total_coupling = sum(abs(coupling) for coupling in couplings)
    local_frustration = (
        sum(abs(coupling) for coupling, signed_product in zip(couplings, signed_products) if signed_product < 0.0)
        / (total_coupling + EPSILON)
    )
    coupling_mean = mean(couplings)
    coupling_variance = variance(couplings)
    perturbation = max(abs(field_scale) * 0.1, 1e-3)
    susceptibility_like = abs(
        safe_tanh((effective_field + perturbation) / temperature)
        - safe_tanh((effective_field - perturbation) / temperature)
    ) / (2.0 * perturbation)

    return {
        "spin_state": center_spin,
        "local_energy": local_energy,
        "local_frustration": local_frustration,
        "coupling_mean": coupling_mean,
        "coupling_variance": coupling_variance,
        "susceptibility_like": susceptibility_like,
    }


def compute_null_metrics(window_residues, center_index, config, gene, position):
    null_energies = []
    null_frustrations = []
    null_coupling_variances = []
    null_susceptibilities = []
    if config["num_shuffles"] <= 0:
        return {
            "null_energy_mean": 0.0,
            "null_energy_std": 0.0,
            "null_frustration_mean": 0.0,
            "null_frustration_std": 0.0,
            "null_coupling_variance_mean": 0.0,
            "null_coupling_variance_std": 0.0,
            "null_susceptibility_mean": 0.0,
            "null_susceptibility_std": 0.0,
        }

    rng = random.Random(stable_seed(gene, position, "".join(window_residues), config["num_shuffles"]))
    for shuffle_index in range(config["num_shuffles"]):
        shuffled = list(window_residues)
        random.Random(rng.randint(0, 10**9) + shuffle_index).shuffle(shuffled)
        metrics = compute_window_metrics(
            shuffled,
            center_index,
            coupling_scale=config["coupling_scale"],
            field_scale=config["field_scale"],
            temperature=config["temperature"],
            encoding=config["encoding"],
        )
        null_energies.append(metrics["local_energy"])
        null_frustrations.append(metrics["local_frustration"])
        null_coupling_variances.append(metrics["coupling_variance"])
        null_susceptibilities.append(metrics["susceptibility_like"])

    null_energy_mean, null_energy_std = summarize_null(null_energies)
    null_frustration_mean, null_frustration_std = summarize_null(null_frustrations)
    null_coupling_variance_mean, null_coupling_variance_std = summarize_null(null_coupling_variances)
    null_susceptibility_mean, null_susceptibility_std = summarize_null(null_susceptibilities)

    return {
        "null_energy_mean": null_energy_mean,
        "null_energy_std": null_energy_std,
        "null_frustration_mean": null_frustration_mean,
        "null_frustration_std": null_frustration_std,
        "null_coupling_variance_mean": null_coupling_variance_mean,
        "null_coupling_variance_std": null_coupling_variance_std,
        "null_susceptibility_mean": null_susceptibility_mean,
        "null_susceptibility_std": null_susceptibility_std,
    }


def compute_spin_glass_score(actual_metrics, null_metrics):
    component_values = [
        (null_metrics["null_energy_mean"] - actual_metrics["local_energy"], 0.0, null_metrics["null_energy_std"]),
        (actual_metrics["local_frustration"], null_metrics["null_frustration_mean"], null_metrics["null_frustration_std"]),
        (
            actual_metrics["coupling_variance"],
            null_metrics["null_coupling_variance_mean"],
            null_metrics["null_coupling_variance_std"],
        ),
        (
            actual_metrics["susceptibility_like"],
            null_metrics["null_susceptibility_mean"],
            null_metrics["null_susceptibility_std"],
        ),
    ]
    zscores = [safe_zscore(value, reference_mean, reference_std) for value, reference_mean, reference_std in component_values]
    return safe_sigmoid(mean(zscores))


def build_output_dataframe(fasta_path, config_path):
    config = load_spin_glass_config(config_path)
    columns = [
        "gene",
        "pos",
        "aa",
        "window_start",
        "window_end",
        "window_length",
        "spin_state",
        "local_energy",
        "local_frustration",
        "coupling_mean",
        "coupling_variance",
        "susceptibility_like",
        "null_energy_mean",
        "null_energy_std",
        "null_frustration_mean",
        "null_frustration_std",
        "null_coupling_variance_mean",
        "null_coupling_variance_std",
        "null_susceptibility_mean",
        "null_susceptibility_std",
        "spin_glass_score",
        "spin_glass_method",
        "spin_glass_version",
    ]
    if not config["enabled"]:
        return pd.DataFrame(columns=columns)

    rows = []
    radius = config["window_size"] // 2
    for gene, sequence in load_fasta_sequences(fasta_path):
        if not gene or not sequence:
            continue
        for index, amino_acid in enumerate(sequence):
            window_start = max(0, index - radius)
            window_end = min(len(sequence), index + radius + 1)
            window_residues = list(sequence[window_start:window_end])
            center_index = index - window_start
            actual_metrics = compute_window_metrics(
                window_residues,
                center_index,
                coupling_scale=config["coupling_scale"],
                field_scale=config["field_scale"],
                temperature=config["temperature"],
                encoding=config["encoding"],
            )
            null_metrics = compute_null_metrics(window_residues, center_index, config, gene, index + 1)
            rows.append(
                {
                    "gene": gene,
                    "pos": index + 1,
                    "aa": amino_acid,
                    "window_start": window_start + 1,
                    "window_end": window_end,
                    "window_length": len(window_residues),
                    **actual_metrics,
                    **null_metrics,
                    "spin_glass_score": compute_spin_glass_score(actual_metrics, null_metrics),
                    "spin_glass_method": METHOD,
                    "spin_glass_version": VERSION,
                }
            )

    return pd.DataFrame(rows, columns=columns)


def main():
    parser = argparse.ArgumentParser(
        description="Compute deterministic spin-glass-inspired residue features from FASTA windows."
    )
    parser.add_argument("--fasta", required=True, help="Input FASTA file")
    parser.add_argument("--config", required=True, help="YAML config containing the spin_glass section")
    parser.add_argument("--output", required=True, help="Output TSV path")
    args = parser.parse_args()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    dataframe = build_output_dataframe(args.fasta, args.config)
    dataframe.to_csv(output_path, sep="\t", index=False)
    print(f"Wrote {len(dataframe)} spin-glass-inspired rows to {output_path}")


if __name__ == "__main__":
    main()
