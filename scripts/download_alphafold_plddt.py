#!/usr/bin/env python3
"""Download per-residue pLDDT confidence scores from EBI AlphaFold database (v4)."""

import argparse
import time
from pathlib import Path

import pandas as pd
import requests

# Canonical UniProt accessions for all 50 CGC genes
UNIPROT_IDS = {
    "ALK": "Q9UM73", "APC": "P25054", "ARID1A": "O14497",
    "ARID2": "Q68CP9", "ATM": "Q13315", "BRAF": "P15056",
    "BRCA1": "P38398", "BRCA2": "P51587", "CDKN2A": "P42771",
    "CHEK2": "O96017", "CTNNB1": "P35222", "EGFR": "P00533",
    "ERBB2": "P04626", "ERCC2": "P18074", "FBXW7": "Q969H0",
    "GNAS": "P63092", "IDH1": "O75874", "IDH2": "P48735",
    "JAK2": "O60674", "KIT": "P10721", "KMT2A": "Q03164",
    "KMT2D": "O14686", "KRAS": "P01116", "MAP2K1": "Q02750",
    "MAP2K2": "P36507", "MET": "P08581", "MLH1": "P40692",
    "MSH2": "P43246", "MSH6": "P52701", "NFE2L2": "Q16236",
    "NOTCH1": "P46531", "NOTCH2": "Q04721", "NRAS": "P01111",
    "NTRK1": "P04629", "PDGFRA": "P16234", "PIK3CA": "P42336",
    "PMS2": "P54278", "PTCH1": "Q13635", "PTEN": "P60484",
    "RB1": "P06400", "RET": "P07949", "RNF43": "Q68DV7",
    "ROS1": "P08922", "SMAD4": "Q13485", "SMARCA4": "P51532",
    "STK11": "Q15831", "TP53": "P04637", "TSC1": "Q92574",
    "TSC2": "P49815", "VHL": "P40337",
}

AF_PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.pdb"


def fetch_plddt_from_pdb(uniprot_id: str, retries: int = 3) -> dict[int, float]:
    """Fetch per-residue pLDDT from AlphaFold PDB file (B-factor column)."""
    url = AF_PDB_URL.format(uid=uniprot_id)
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=60)
            if r.status_code == 200:
                return _parse_pdb_bfactors(r.text)
            if r.status_code == 404:
                return {}          # protein not in AlphaFold DB
        except requests.RequestException:
            pass
        time.sleep(2 ** attempt)
    return {}


def _parse_pdb_bfactors(pdb_text: str) -> dict[int, float]:
    """Extract one pLDDT value per residue from ATOM CA records."""
    plddt = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        try:
            res_seq = int(line[22:26])
            b_factor = float(line[60:66])
            plddt[res_seq] = b_factor
        except ValueError:
            continue
    return plddt


def main():
    parser = argparse.ArgumentParser(description="Download AlphaFold2 pLDDT scores for all CGC genes.")
    parser.add_argument("--output", default="results/alphafold/plddt_scores.tsv")
    args = parser.parse_args()

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for gene, uid in UNIPROT_IDS.items():
        print(f"  {gene} ({uid}) ...", end=" ", flush=True)
        plddt = fetch_plddt_from_pdb(uid)
        if plddt:
            for pos, score in sorted(plddt.items()):
                rows.append({"gene": gene, "pos": pos, "plddt": round(score, 2)})
            print(f"{len(plddt)} residues")
        else:
            print("not found — skipped")
        time.sleep(0.5)   # polite rate limiting

    df = pd.DataFrame(rows)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\n✓ Saved {len(df)} residue pLDDT scores → {out_path}")


if __name__ == "__main__":
    main()
