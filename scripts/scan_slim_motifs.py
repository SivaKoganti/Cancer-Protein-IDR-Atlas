#!/usr/bin/env python3
"""
Scan cancer protein sequences for cancer-relevant Short Linear Motifs (SLiMs).

Uses a curated set of ELM-derived regex patterns covering:
  - Cell-cycle kinase substrates (CDK, PLK1, ATM, CK2, GSK3)
  - Protein degradation signals (D-box, KEN-box, SCF-degron, PEST)
  - Protein interaction docking motifs (SH2, SH3, 14-3-3, BRCT, WW, PDZ)
  - Nuclear localisation / export signals (NLS, NES)
  - SUMO modification consensus

For each motif match:
  - Records the genomic start/end position (1-based)
  - Flags whether a ClinVar variant overlaps the motif (potential gain/loss)
  - Assigns an IDR context score (mean iupred_score across match span)
"""

import argparse
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def normalize_gene_id(gene):
    if gene is None:
        return ""
    value = str(gene).strip()
    if not value:
        return ""
    return re.sub(r"_\d+$", "", value)

# ── curated SLiM catalogue ──────────────────────────────────────────────────
# Each entry: (motif_id, category, description, regex_pattern)
SLIM_CATALOGUE = [
    # Cell-cycle kinase substrates
    ("MOD_CDK_1",    "phospho", "CDK substrate [ST]-P-x-[RK]",      r"[ST]P.[RK]"),
    ("MOD_CDK_SPxK", "phospho", "CDK minimal [ST]-P",               r"[ST]P"),
    ("MOD_MAPK_ERK", "phospho", "ERK docking P-x-[ST]-P",           r"P.[ST]P"),
    ("MOD_PLK1_1",   "phospho", "PLK1 substrate [DE]-x-[ST]-[FY]",  r"[DE].[ST][FY]"),
    ("MOD_ATM_1",    "phospho", "ATM/ATR [ST]-Q",                   r"[ST]Q"),
    ("MOD_CK2_1",    "phospho", "CK2 [ST]-x-x-[DE]",               r"[ST]..[DE]"),
    ("MOD_GSK3_1",   "phospho", "GSK3 primed [ST]-x-x-x-[ST]",     r"[ST]...[ST]"),
    ("MOD_AURORA",   "phospho", "Aurora kinase [RK]-x-[ST]",        r"[RK].[ST]"),
    # Degradation signals
    ("DEG_APCC_DBOX","degron",  "APC/C D-box R-x-x-L",             r"R..[LI]"),
    ("DEG_APCC_KEN", "degron",  "APC/C KEN-box K-E-N",             r"KEN"),
    ("DEG_SCF_TRCP", "degron",  "SCF-βTrCP D[ST]Gxx[ST]",         r"D[ST]G..[ST]"),
    ("DEG_SCF_SKP2", "degron",  "SCF-Skp2 [LI]-P",                r"[LI]P"),
    ("DEG_PEST",     "degron",  "PEST region (rapid turnover)",     r"P[EDS].{2,8}[ST]"),
    # Interaction docking motifs
    ("LIG_SH3_1",    "docking", "SH3 ligand P-x-x-P",              r"P..P"),
    ("LIG_SH3_RxxP", "docking", "SH3 type II R-x-x-P",            r"[RK]..P"),
    ("LIG_14-3-3_1", "docking", "14-3-3 [RK]xx[ST]xP",            r"[RK]..[ST].P"),
    ("LIG_14-3-3_2", "docking", "14-3-3 [ST]xxx[LI]",             r"[ST]...[LI]"),
    ("LIG_BRCT_1",   "docking", "BRCT [ST]-P-x-[FY]",             r"[ST]P.[FY]"),
    ("LIG_WW_1",     "docking", "WW domain P-P-x-Y",              r"PP.Y"),
    ("LIG_PDZ_1",    "docking", "PDZ [ST]-x-[VLI] C-terminal",     r"[ST].[VLI]$"),
    ("LIG_PCNA_PIP", "docking", "PCNA PIP-box Qxx[^DE][^DE]DA",    r"Q..[^DE][^DE]D[^DE]"),
    # Nuclear signals
    ("TRG_NLS_mono", "localise","Monopartite NLS K-K/R-x-K/R",     r"K[KR].[KR]"),
    ("TRG_NLS_bi",   "localise","Bipartite NLS [KR]{2}.{9,12}[KR]{2}", r"[KR]{2}.{9,12}[KR]{2}"),
    ("TRG_NES_1",    "localise","CRM1 NES L-x{2,3}-[LIVFM]-x{2,3}-L", r"L.{2,3}[LIVFM].{2,3}L"),
    # SUMO
    ("MOD_SUMO_1",   "sumo",   "SUMO consensus ψ-K-x-[ED]",        r"[LIVMF]K.[ED]"),
    ("MOD_SUMO_rev", "sumo",   "Inverted SUMO [ED]-x-K-ψ",         r"[ED].K[LIVMF]"),
]

PTM_CATEGORY_PREFIX = "MOD_"


def scan_sequence(seq: str, gene: str, iupred_scores: list[float]) -> list[dict]:
    """Return one row per SLiM hit with IDR context annotation."""
    hits = []
    seq_upper = seq.upper()
    n = len(seq_upper)

    for motif_id, category, description, pattern in SLIM_CATALOGUE:
        for m in re.finditer(pattern, seq_upper):
            start_1 = m.start() + 1   # 1-based
            end_1 = m.end()            # inclusive
            span_idr = [iupred_scores[i] for i in range(m.start(), m.end())
                        if i < len(iupred_scores)]
            mean_idr = round(sum(span_idr) / len(span_idr), 4) if span_idr else 0.0
            hits.append({
                "gene":        gene,
                "motif_id":    motif_id,
                "category":    category,
                "description": description,
                "start":       start_1,
                "end":         end_1,
                "match_seq":   m.group(),
                "mean_idr":    mean_idr,
                "in_idr":      int(mean_idr >= 0.5),
            })
    return hits


def annotate_variant_overlap(hits_df: pd.DataFrame, variants_tsv: str | None) -> pd.DataFrame:
    """Flag SLiM hits that overlap a ClinVar variant position."""
    hits_df["clinvar_overlap"] = 0
    if not variants_tsv or not Path(variants_tsv).exists():
        return hits_df

    var_df = pd.read_csv(variants_tsv, sep="\t", header=None,
                         names=["gene", "pos", "ref", "alt", "clin"])
    var_positions = set(zip(
        var_df["gene"].astype(str).map(normalize_gene_id),
        var_df["pos"].astype(int),
    ))

    def overlaps(row):
        gene_key = normalize_gene_id(row["gene"])
        for p in range(int(row["start"]), int(row["end"]) + 1):
            if (gene_key, p) in var_positions:
                return 1
        return 0

    hits_df["clinvar_overlap"] = hits_df.apply(overlaps, axis=1)
    return hits_df


def build_per_residue_slim_counts(hits_df: pd.DataFrame,
                                  gene: str, length: int) -> pd.DataFrame:
    """Aggregate SLiM hits into per-residue count and IDs columns."""
    count = [0] * length
    ids: list[list[str]] = [[] for _ in range(length)]
    gene_hits = hits_df[hits_df["gene"].astype(str).str.strip().map(normalize_gene_id) == normalize_gene_id(gene)]
    for _, row in gene_hits.iterrows():
        for p in range(int(row["start"]) - 1, min(int(row["end"]), length)):
            count[p] += 1
            ids[p].append(row["motif_id"])

    return pd.DataFrame({
        "gene": gene,
        "pos": range(1, length + 1),
        "slim_count": count,
        "slim_ids": [";".join(sorted(set(x))) for x in ids],
    })


def extract_ptm_hits(hits_df: pd.DataFrame) -> pd.DataFrame:
    """Keep SLiM hits that correspond to PTM motifs."""
    if hits_df.empty:
        return pd.DataFrame(columns=[
            "gene", "motif_id", "category", "description", "start", "end",
            "match_seq", "mean_idr", "in_idr", "clinvar_overlap",
        ])
    return hits_df[hits_df["motif_id"].astype(str).str.startswith(PTM_CATEGORY_PREFIX)].copy()


def build_per_residue_ptm_annotations(ptm_hits_df: pd.DataFrame,
                                      gene: str, length: int) -> pd.DataFrame:
    """Aggregate PTM motif hits into per-residue PTM-centric annotations."""
    count = [0] * length
    ids: list[list[str]] = [[] for _ in range(length)]
    categories: list[list[str]] = [[] for _ in range(length)]
    idr_hits = [0] * length
    clinvar_hits = [0] * length

    gene_hits = ptm_hits_df[
        ptm_hits_df["gene"].astype(str).str.strip().map(normalize_gene_id) == normalize_gene_id(gene)
    ]

    for _, row in gene_hits.iterrows():
        is_idr = int(row.get("in_idr", 0))
        has_clinvar_overlap = int(row.get("clinvar_overlap", 0))
        for p in range(int(row["start"]) - 1, min(int(row["end"]), length)):
            count[p] += 1
            ids[p].append(str(row["motif_id"]))
            categories[p].append(str(row.get("category", "")))
            idr_hits[p] += is_idr
            clinvar_hits[p] += has_clinvar_overlap

    return pd.DataFrame({
        "gene": gene,
        "pos": range(1, length + 1),
        "ptm_count": count,
        "ptm_categories": [";".join(sorted(set([c for c in vals if c]))) for vals in categories],
        "ptm_ids": [";".join(sorted(set(vals))) for vals in ids],
        "ptm_in_idr_count": idr_hits,
        "ptm_clinvar_overlap_count": clinvar_hits,
    })


def main():
    parser = argparse.ArgumentParser(description="Scan cancer proteins for SLiM motifs.")
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--iupred", required=True,
                        help="IUPred scores TSV (gene, pos, score)")
    parser.add_argument("--variants", default="results/clinvar/mapped_variants.tsv")
    parser.add_argument("--output-dir", default="results/slim")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    iupred_df = pd.read_csv(args.iupred, sep="\t")

    all_hits = []
    residue_frames = []

    for rec in SeqIO.parse(args.fasta, "fasta"):
        gid = rec.id.split("|")[0].strip()
        base = normalize_gene_id(gid)
        seq = str(rec.seq).upper()
        length = len(seq)

        gene_iupred = iupred_df[
            iupred_df["gene"].astype(str).str.strip().map(normalize_gene_id) == base
        ]
        iupred_scores = [0.0] * length
        for _, row in gene_iupred.iterrows():
            idx = int(row["pos"]) - 1
            if 0 <= idx < length:
                iupred_scores[idx] = float(row.get("score", 0.0))

        hits = scan_sequence(seq, base, iupred_scores)
        all_hits.extend(hits)
        print(f"  {base}: {len(hits)} SLiM hits")

    hits_df = pd.DataFrame(all_hits)
    if not hits_df.empty:
        hits_df = annotate_variant_overlap(hits_df, args.variants)

    hits_df.to_csv(out_dir / "slim_hits.tsv", sep="\t", index=False)
    ptm_hits_df = extract_ptm_hits(hits_df)
    ptm_hits_df.to_csv(out_dir / "ptm_hits.tsv", sep="\t", index=False)

    if not hits_df.empty:
        cat_summary = (
            hits_df.groupby(["gene", "category"])
                   .size()
                   .reset_index(name="count")
        )
        cat_summary.to_csv(out_dir / "slim_category_summary.tsv", sep="\t", index=False)

        for rec in SeqIO.parse(args.fasta, "fasta"):
            gid = rec.id.split("|")[0].strip()
            base = normalize_gene_id(gid)
            rdf = build_per_residue_slim_counts(hits_df, base, len(rec.seq))
            residue_frames.append(rdf)

        residue_df = pd.concat(residue_frames, ignore_index=True)
    else:
        residue_df = pd.DataFrame(columns=["gene", "pos", "slim_count", "slim_ids"])

    ptm_residue_frames = []
    for rec in SeqIO.parse(args.fasta, "fasta"):
        gid = rec.id.split("|")[0].strip()
        base = normalize_gene_id(gid)
        ptm_rdf = build_per_residue_ptm_annotations(ptm_hits_df, base, len(rec.seq))
        ptm_residue_frames.append(ptm_rdf)

    ptm_residue_df = pd.concat(ptm_residue_frames, ignore_index=True) if ptm_residue_frames else pd.DataFrame(
        columns=[
            "gene", "pos", "ptm_count", "ptm_categories", "ptm_ids",
            "ptm_in_idr_count", "ptm_clinvar_overlap_count",
        ]
    )

    residue_df.to_csv(out_dir / "slim_per_residue.tsv", sep="\t", index=False)
    ptm_residue_df.to_csv(out_dir / "ptm_per_residue.tsv", sep="\t", index=False)

    if not ptm_hits_df.empty:
        ptm_summary = (
            ptm_hits_df.groupby(["gene", "category"], as_index=False)
            .size()
            .rename(columns={"size": "count"})
        )
    else:
        ptm_summary = pd.DataFrame(columns=["gene", "category", "count"])
    ptm_summary.to_csv(out_dir / "ptm_category_summary.tsv", sep="\t", index=False)

    print(f"\n✓ SLiM scan complete → {out_dir}")
    print(f"  Total hits: {len(hits_df)}")
    if not hits_df.empty:
        print(f"  Hits overlapping ClinVar variants: {hits_df['clinvar_overlap'].sum()}")
        print(f"  Hits in IDR regions (mean_idr ≥ 0.5): {hits_df['in_idr'].sum()}")
    print(f"  PTM motif hits: {len(ptm_hits_df)}")


if __name__ == "__main__":
    main()
