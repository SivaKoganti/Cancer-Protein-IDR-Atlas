#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from collections import defaultdict

import pandas as pd
import yaml
from Bio import SeqIO


def normalize_gene_id(gene):
    if gene is None:
        return ""
    value = str(gene).strip()
    if not value:
        return ""
    return re.sub(r"_\d+$", "", value)


def gene_matches(gene_value, target_gene):
    return normalize_gene_id(gene_value) == normalize_gene_id(target_gene)


def load_config(config_path):
    if not config_path or not Path(config_path).exists():
        return {}
    with open(config_path) as fh:
        return yaml.safe_load(fh) or {}


def load_gene_list(config_path):
    return load_config(config_path).get("genes", [])


def load_vipp_weights(config_path):
    cfg = load_config(config_path)
    vipp_cfg = cfg.get("vipp", {})
    return {
        "idr_weight": float(vipp_cfg.get("idr_weight", 0.35)),
        "llps_weight": float(vipp_cfg.get("llps_weight", 0.35)),
        "virus_weight": float(vipp_cfg.get("virus_weight", 0.30)),
    }


def load_fasta_sequences(path):
    if not path or not Path(path).exists():
        return {}
    sequences = {}
    for rec in SeqIO.parse(path, "fasta"):
        gid = normalize_gene_id(rec.id.split("|")[0].strip())
        sequences[gid] = str(rec.seq)
    return sequences


def load_optional_tsv(path, *args, **kwargs):
    if not path or not Path(path).exists():
        return None
    try:
        return pd.read_csv(path, sep='\t', *args, **kwargs)
    except pd.errors.EmptyDataError:
        return None


def bits_from_regions(regions_df, length):
    flags = [0] * length
    if regions_df is None:
        return flags
    for _, row in regions_df.iterrows():
        start = int(row["start"]) - 1
        end = int(row["end"])
        for pos in range(start, min(end, length)):
            flags[pos] = 1
    return flags


def window_scores(plaac_df, length):
    scores = [0.0] * length
    if plaac_df is None:
        return scores
    for _, row in plaac_df.iterrows():
        start = int(row["start"]) - 1
        end = int(row["end"])
        qscore = float(row.get("q_n_fraction", 0.0))
        for pos in range(start, min(end, length)):
            scores[pos] = max(scores[pos], qscore)
    return scores


def group_by_gene(df, gene):
    if df is None or df.empty:
        return df
    return df[df["gene"].astype(str).str.strip().map(normalize_gene_id) == normalize_gene_id(gene)]


def clip01(value):
    try:
        score = float(value)
    except (TypeError, ValueError):
        return 0.0
    if pd.isna(score):
        return 0.0
    return max(0.0, min(1.0, score))


def pos_lookup(df, gene):
    if df is None or df.empty:
        return {}
    sub = df[df["gene"].astype(str).str.strip().map(normalize_gene_id) == normalize_gene_id(gene)]
    index = {}
    for _, row in sub.iterrows():
        try:
            pos = int(row["pos"])
        except (KeyError, TypeError, ValueError):
            continue
        if pd.notna(pos):
            index[pos] = row
    return index


def map_variants(variants_df):
    variant_map = defaultdict(list)
    if variants_df is None:
        return variant_map
    for _, row in variants_df.iterrows():
        gene = normalize_gene_id(row["gene"])
        variant_map[gene].append(f"{row['pos']}{row['ref']}>{row['alt']}({row['clin']})")
    return variant_map


def main():
    parser = argparse.ArgumentParser(description="Build a phylogeny-aware disorder atlas per gene.")
    parser.add_argument("--iupred", required=True)
    parser.add_argument("--seg", required=True)
    parser.add_argument("--plaac", required=True)
    parser.add_argument("--llps", required=True)
    parser.add_argument("--structure", required=True)
    parser.add_argument("--variants")
    parser.add_argument("--conservation", required=True)
    parser.add_argument("--fasta")
    parser.add_argument("--genes-file", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--plddt", default=None)
    parser.add_argument("--delta-llps", dest="delta_llps", default=None)
    parser.add_argument("--cdr", default=None)
    parser.add_argument("--slim", default=None)
    parser.add_argument("--ptm", default=None)
    parser.add_argument("--oncovirus", default=None)
    args = parser.parse_args()

    gene_list = load_gene_list(args.genes_file)
    vipp_weights = load_vipp_weights(args.genes_file)
    sequences = load_fasta_sequences(args.fasta)
    iupred = pd.read_csv(args.iupred, sep='\t')
    seg = pd.read_csv(args.seg, sep='\t')
    plaac = pd.read_csv(args.plaac, sep='\t')
    llps = pd.read_csv(args.llps, sep='\t')
    structure = pd.read_csv(args.structure, sep='\t')
    variants = None
    if args.variants and Path(args.variants).exists():
        variants = pd.read_csv(args.variants, sep='\t', header=None, names=['gene', 'pos', 'ref', 'alt', 'clin'])
    conservation = pd.read_csv(args.conservation, sep='\t')

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    plddt      = load_optional_tsv(args.plddt)
    delta_llps = load_optional_tsv(args.delta_llps)
    cdr_data   = load_optional_tsv(args.cdr)
    slim_data  = load_optional_tsv(args.slim)
    ptm_data   = load_optional_tsv(args.ptm)
    oncovirus  = load_optional_tsv(args.oncovirus)

    variant_map = map_variants(variants)
    summary_rows = []

    for gene in gene_list:
        gene_iupred = iupred[iupred["gene"].astype(str).str.strip().map(normalize_gene_id) == normalize_gene_id(gene)].copy()
        seq = sequences.get(normalize_gene_id(gene))
        if seq is None:
            if gene_iupred.empty or "aa" not in gene_iupred.columns:
                continue
            seq = ''.join(str(aa) for aa in gene_iupred["aa"].tolist() if str(aa).strip())
        if not seq:
            continue
        length = len(seq)

        low_complexity = bits_from_regions(group_by_gene(seg, gene), length)
        plaac_scores = window_scores(group_by_gene(plaac, gene), length)
        conservation_scores = [0.0] * length
        gene_cons = conservation[conservation["gene"].astype(str).str.strip().map(normalize_gene_id) == normalize_gene_id(gene)]
        for _, row in gene_cons.iterrows():
            try:
                idx = int(row["pos"]) - 1
            except (KeyError, TypeError, ValueError):
                continue
            if 0 <= idx < length:
                conservation_scores[idx] = float(row.get("conservation_score", 0.0))

        variant_list = [v for v in variant_map.get(normalize_gene_id(gene), [])]
        gene_variants_by_pos = defaultdict(list)
        for variant in variant_list:
            pos = int(''.join(ch for ch in variant if ch.isdigit()))
            gene_variants_by_pos[pos].append(variant)

        table = []
        llps_scores = [0.0] * length
        gene_llps = llps[llps["gene"].astype(str).str.strip().map(normalize_gene_id) == normalize_gene_id(gene)]
        for _, row in gene_llps.iterrows():
            try:
                start = int(row["start"]) - 1
                end = int(row["end"])
                score = float(row.get("llps_proxy", 0.0))
            except (KeyError, TypeError, ValueError):
                continue
            if start >= end or end <= 0:
                continue
            for idx in range(max(0, start), min(end, length)):
                llps_scores[idx] = max(llps_scores[idx], score)

        structure_scores = [0.0] * length
        quantum_scores = [0.0] * length
        gene_structure = structure[structure["gene"].astype(str).str.strip().map(normalize_gene_id) == normalize_gene_id(gene)]
        for _, row in gene_structure.iterrows():
            try:
                idx = int(row["pos"]) - 1
            except (KeyError, TypeError, ValueError):
                continue
            if 0 <= idx < length:
                structure_scores[idx] = float(row.get("structural_proxy", 0.0))
                quantum_scores[idx] = float(row.get("quantum_proxy", 0.0))

        # ── index novelty layers by position for O(1) lookup ──────────────
        plddt_idx      = pos_lookup(plddt, gene)
        delta_llps_idx = pos_lookup(delta_llps, gene)
        cdr_idx        = pos_lookup(cdr_data, gene)
        slim_idx       = pos_lookup(slim_data, gene)
        ptm_idx        = pos_lookup(ptm_data, gene)
        oncovirus_idx  = pos_lookup(oncovirus, gene)

        for pos, aa in enumerate(seq, start=1):
            iupred_score = 0.0
            pos_matches = gene_iupred[gene_iupred["pos"] == pos]
            if not pos_matches.empty:
                iupred_score = float(pos_matches.iloc[0].get("score", 0.0))

            pr = plddt_idx.get(pos)
            dr = delta_llps_idx.get(pos)
            cr = cdr_idx.get(pos)
            sr = slim_idx.get(pos)
            pr_ptm = ptm_idx.get(pos)
            ov = oncovirus_idx.get(pos)

            idr_component = clip01(iupred_score)
            llps_component = clip01(llps_scores[pos-1] if pos-1 < len(llps_scores) else 0.0)
            virus_component = 1.0 if (ov is not None and int(ov.get("virus_interaction", 0)) > 0) else 0.0
            vipp_score = (
                vipp_weights["idr_weight"] * idr_component
                + vipp_weights["llps_weight"] * llps_component
                + vipp_weights["virus_weight"] * virus_component
            )
            vipp_score = clip01(vipp_score)

            table.append({
                "gene": gene,
                "pos": pos,
                "aa": aa,
                "iupred_score": iupred_score,
                "low_complexity": low_complexity[pos-1] if pos-1 < len(low_complexity) else 0,
                "plaac_qn": plaac_scores[pos-1] if pos-1 < len(plaac_scores) else 0.0,
                "llps_proxy": llps_scores[pos-1] if pos-1 < len(llps_scores) else 0.0,
                "structural_proxy": structure_scores[pos-1] if pos-1 < len(structure_scores) else 0.0,
                "quantum_proxy": quantum_scores[pos-1] if pos-1 < len(quantum_scores) else 0.0,
                "conservation": conservation_scores[pos-1] if pos-1 < len(conservation_scores) else 0.0,
                "variants": ";".join(gene_variants_by_pos.get(pos, [])),
                "virus_interaction": int(ov.get("virus_interaction", 0)) if ov is not None else 0,
                "virus_count": int(ov.get("virus_count", 0)) if ov is not None else 0,
                "virus_names": str(ov.get("virus_names", "")) if ov is not None else "",
                "virus_evidence": str(ov.get("virus_evidence", "")) if ov is not None else "",
                "vipp_score": vipp_score,
                # ── novelty layers ──────────────────────────────────────────
                "plddt":            float(pr["plddt"]) if pr is not None and "plddt" in pr else float("nan"),
                "llps_vulnerability": float(dr["llps_vulnerability"]) if dr is not None else float("nan"),
                "delta_llps_max":   float(dr["delta_llps_max"]) if dr is not None else float("nan"),
                "delta_llps_min":   float(dr["delta_llps_min"]) if dr is not None else float("nan"),
                "cdr_flag":         int(cr["cdr_flag"]) if cr is not None and "cdr_flag" in cr else 0,
                "discordant_flag":  int(cr["discordant_flag"]) if cr is not None and "discordant_flag" in cr else 0,
                "slim_count":       int(sr["slim_count"]) if sr is not None and "slim_count" in sr else 0,
                "slim_ids":         str(sr["slim_ids"]) if sr is not None and "slim_ids" in sr else "",
                "ptm_count":        int(pr_ptm["ptm_count"]) if pr_ptm is not None and "ptm_count" in pr_ptm else 0,
                "ptm_categories":   str(pr_ptm["ptm_categories"]) if pr_ptm is not None and "ptm_categories" in pr_ptm else "",
                "ptm_ids":          str(pr_ptm["ptm_ids"]) if pr_ptm is not None and "ptm_ids" in pr_ptm else "",
                "ptm_in_idr_count": int(pr_ptm["ptm_in_idr_count"]) if pr_ptm is not None and "ptm_in_idr_count" in pr_ptm else 0,
                "ptm_clinvar_overlap_count": int(pr_ptm["ptm_clinvar_overlap_count"]) if pr_ptm is not None and "ptm_clinvar_overlap_count" in pr_ptm else 0,
            })

        gene_df = pd.DataFrame(table)
        if not gene_df.empty:
            gene_df.to_csv(outdir / f"{gene}_atlas.tsv", sep='\t', index=False)
            summary_rows.append({
                "gene": gene,
                "length": length,
                "mean_iupred": gene_df["iupred_score"].mean(),
                "mean_llps": gene_df["llps_proxy"].mean(),
                "mean_structural": gene_df["structural_proxy"].mean(),
                "mean_quantum": gene_df["quantum_proxy"].mean(),
                "mean_conservation": gene_df["conservation"].mean(),
                "low_complexity_fraction": gene_df["low_complexity"].mean(),
                "variant_count": len(variant_list),
                "mean_plddt": gene_df["plddt"].mean() if "plddt" in gene_df else float("nan"),
                "mean_llps_vulnerability": gene_df["llps_vulnerability"].mean() if "llps_vulnerability" in gene_df else float("nan"),
                "cdr_fraction": gene_df["cdr_flag"].mean() if "cdr_flag" in gene_df else 0.0,
                "discordant_fraction": gene_df["discordant_flag"].mean() if "discordant_flag" in gene_df else 0.0,
                "slim_density": gene_df["slim_count"].mean() if "slim_count" in gene_df else 0.0,
                "ptm_density": gene_df["ptm_count"].mean() if "ptm_count" in gene_df else 0.0,
                "ptm_idr_density": gene_df["ptm_in_idr_count"].mean() if "ptm_in_idr_count" in gene_df else 0.0,
                "mean_vipp_score": gene_df["vipp_score"].mean() if "vipp_score" in gene_df else float("nan"),
                "max_vipp_score": gene_df["vipp_score"].max() if "vipp_score" in gene_df else float("nan"),
                "virus_interaction_fraction": gene_df["virus_interaction"].mean() if "virus_interaction" in gene_df else 0.0,
                "mean_virus_count": gene_df["virus_count"].mean() if "virus_count" in gene_df else 0.0,
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(outdir / "global_disorder_phylogeny_atlas.tsv", sep='\t', index=False)
    print(f"Wrote atlas summaries to {outdir}")


if __name__ == '__main__':
    main()
