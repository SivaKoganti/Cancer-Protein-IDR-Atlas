#!/usr/bin/env python3
"""Build a residue-level structural library for IDR/LLPS states.

The output is a compact, interpretable representation that can be used to
summarize each residue in terms of an LLPS structural state, a coarse-model
coordinate sketch, a simulated crystallogram-style profile, and a quantum
biophysical proxy.
"""
import argparse
import json
from pathlib import Path

import pandas as pd


def llps_state_from_scores(iupred_score, llps_proxy, structural_proxy, quantum_proxy):
    if llps_proxy >= 0.75 and iupred_score >= 0.65:
        return "condensate"
    if llps_proxy >= 0.45 and structural_proxy <= 0.35:
        return "expanded"
    return "compact"


def coarse_model_xyz(iupred_score, llps_proxy, structural_proxy, quantum_proxy):
    x = round(0.5 + 0.35 * iupred_score, 3)
    y = round(0.5 + 0.35 * llps_proxy, 3)
    z = round(0.5 + 0.25 * structural_proxy + 0.15 * quantum_proxy, 3)
    return f"[{x:.3f}, {y:.3f}, {z:.3f}]"


def simulated_crystallogram(iupred_score, llps_proxy, structural_proxy, quantum_proxy):
    peak = round(max(iupred_score, llps_proxy, structural_proxy, quantum_proxy), 3)
    width = round(0.3 + 0.2 * llps_proxy + 0.1 * iupred_score, 3)
    return {"peak": peak, "width": width, "shape": "broad" if peak >= 0.6 else "narrow"}


def quantum_biophysical_score(iupred_score, llps_proxy, structural_proxy, quantum_proxy):
    return round(max(0.0, min(1.0, 0.35 * iupred_score + 0.35 * llps_proxy + 0.15 * structural_proxy + 0.15 * quantum_proxy)), 3)


def main():
    parser = argparse.ArgumentParser(description="Build a structural library for IDR residues")
    parser.add_argument("input_atlas", help="Residue atlas TSV with IDR and LLPS-related columns")
    parser.add_argument("output", help="Output TSV containing residue structural libraries")
    args = parser.parse_args()

    atlas = pd.read_csv(args.input_atlas, sep="\t")
    required = {"gene", "pos", "iupred_score", "llps_proxy", "structural_proxy", "quantum_proxy"}
    missing = required - set(atlas.columns)
    if missing:
        raise SystemExit(f"Input atlas is missing required columns: {sorted(missing)}")

    rows = []
    for _, row in atlas.iterrows():
        iupred_score = float(row.get("iupred_score", 0.0))
        llps_proxy = float(row.get("llps_proxy", 0.0))
        structural_proxy = float(row.get("structural_proxy", 0.0))
        quantum_proxy = float(row.get("quantum_proxy", 0.0))

        state = llps_state_from_scores(iupred_score, llps_proxy, structural_proxy, quantum_proxy)
        crystal = simulated_crystallogram(iupred_score, llps_proxy, structural_proxy, quantum_proxy)
        rows.append({
            "gene": row.get("gene", ""),
            "pos": int(row.get("pos", 0)),
            "aa": row.get("aa", ""),
            "iupred_score": iupred_score,
            "llps_proxy": llps_proxy,
            "structural_proxy": structural_proxy,
            "quantum_proxy": quantum_proxy,
            "llps_state": state,
            "coarse_model_xyz": coarse_model_xyz(iupred_score, llps_proxy, structural_proxy, quantum_proxy),
            "simulated_crystallogram": json.dumps(crystal),
            "quantum_biophysical_score": quantum_biophysical_score(iupred_score, llps_proxy, structural_proxy, quantum_proxy),
            "structural_library_json": json.dumps({
                "state": state,
                "coarse_model_xyz": coarse_model_xyz(iupred_score, llps_proxy, structural_proxy, quantum_proxy),
                "crystallogram": crystal,
                "quantum_biophysical_score": quantum_biophysical_score(iupred_score, llps_proxy, structural_proxy, quantum_proxy),
            }),
        })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
