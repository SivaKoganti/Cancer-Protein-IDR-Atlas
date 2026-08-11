#!/usr/bin/env python3
"""Benchmark VIPP against simpler baselines using ClinVar-derived labels."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


POSITIVE_KEYWORDS = (
    "pathogenic",
    "likely pathogenic",
)

NEGATIVE_KEYWORDS = (
    "benign",
    "likely benign",
)


def _target_from_clinvar(clinvar_value):
    """Return 1 for pathogenic, 0 for benign, and None when label is ambiguous/unusable."""
    if pd.isna(clinvar_value):
        return None

    normalized = str(clinvar_value).strip().lower()
    if not normalized:
        return None

    has_positive = any(keyword in normalized for keyword in POSITIVE_KEYWORDS)
    has_negative = any(keyword in normalized for keyword in NEGATIVE_KEYWORDS)

    # Mixed or conflicting classes are excluded to avoid noisy supervision.
    if has_positive and not has_negative:
        return 1
    if has_negative and not has_positive:
        return 0
    return None


def _prepare_benchmark_frame(atlas_dir):
    frames = []
    for path in sorted(atlas_dir.glob("*_atlas_with_clinvar.tsv")):
        df = pd.read_csv(path, sep="\t")
        required = {"clinvar_significance", "vipp_score", "virus_interaction", "iupred_score", "llps_proxy", "conservation"}
        if required.issubset(df.columns):
            df = df.copy()
            df["target"] = df["clinvar_significance"].map(_target_from_clinvar)
            df = df[df["target"].isin([0, 1])].copy()
            if not df.empty:
                frames.append(df)

    if not frames:
        return pd.DataFrame()

    combined = pd.concat(frames, ignore_index=True)
    combined = combined[combined["target"].isin([0, 1])].copy()
    return combined


def _score_row(row, model_name):
    if model_name == "idr_only":
        return float(row.get("iupred_score", 0.0))
    if model_name == "idr_llps":
        return 0.5 * float(row.get("iupred_score", 0.0)) + 0.5 * float(row.get("llps_proxy", 0.0))
    if model_name == "idr_llps_virus":
        return 0.4 * float(row.get("iupred_score", 0.0)) + 0.3 * float(row.get("llps_proxy", 0.0)) + 0.3 * float(row.get("virus_interaction", 0.0))
    if model_name == "idr_llps_virus_conservation":
        return 0.35 * float(row.get("iupred_score", 0.0)) + 0.25 * float(row.get("llps_proxy", 0.0)) + 0.2 * float(row.get("virus_interaction", 0.0)) + 0.2 * float(row.get("conservation", 0.0))
    return float(row.get("vipp_score", 0.0))


def _rankdata(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[order[end]] == values[order[start]]:
            end += 1
        avg_rank = 0.5 * (start + 1 + end)
        ranks[order[start:end]] = avg_rank
        start = end
    return ranks


def _auroc(y_true, y_score):
    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)
    pos = y_true == 1
    neg = ~pos
    pos_count = int(pos.sum())
    neg_count = int(neg.sum())
    if pos_count == 0 or neg_count == 0:
        return float("nan")

    ranks = _rankdata(y_score)
    rank_sum_pos = float(ranks[pos].sum())
    u_stat = rank_sum_pos - pos_count * (pos_count + 1) / 2.0
    return u_stat / (pos_count * neg_count)


def _average_precision(y_true, y_score):
    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)
    order = np.argsort(-y_score, kind="mergesort")
    ranked_true = y_true[order]
    positives = int(ranked_true.sum())
    if positives == 0:
        return float("nan")

    true_positives = 0
    precision_sum = 0.0
    for idx, label in enumerate(ranked_true, start=1):
        if label == 1:
            true_positives += 1
            precision_sum += true_positives / idx
    return precision_sum / positives


def benchmark_models(atlas_dir, output_dir):
    atlas_dir = Path(atlas_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    combined = _prepare_benchmark_frame(atlas_dir)

    models = ["idr_only", "idr_llps", "idr_llps_virus", "idr_llps_virus_conservation", "full_vipp"]
    rows = []

    if combined.empty:
        for model_name in models:
            rows.append({
                "model_name": model_name,
                "auroc": float("nan"),
                "auprc": float("nan"),
                "topk_enrichment": float("nan"),
            })
        metrics_df = pd.DataFrame(rows, columns=["model_name", "auroc", "auprc", "topk_enrichment"])
        metrics_df.to_csv(output_dir / "model_metrics.tsv", sep="\t", index=False)
        metrics_df.to_csv(output_dir / "ablation_matrix.tsv", sep="\t", index=False)
        return metrics_df, metrics_df.copy()

    for model_name in models:
        scores = combined.apply(lambda row: _score_row(row, model_name), axis=1)
        y_true = combined["target"].to_numpy()
        y_score = np.asarray(scores, dtype=float)

        if len(np.unique(y_true)) < 2:
            rows.append({
                "model_name": model_name,
                "auroc": float("nan"),
                "auprc": float("nan"),
                "topk_enrichment": float("nan"),
            })
            continue

        topk = int(max(1, round(0.1 * len(y_true))))
        if y_true.sum() == 0:
            topk_enrichment = float("nan")
            auroc = float("nan")
            auprc = float("nan")
            rows.append({
                "model_name": model_name,
                "auroc": auroc,
                "auprc": auprc,
                "topk_enrichment": topk_enrichment,
            })
            continue
        order = np.argsort(-y_score)
        top_hits = y_true[order[:topk]].sum() / max(1, topk)
        baseline_hits = y_true.mean()
        topk_enrichment = top_hits / baseline_hits if baseline_hits > 0 else float("nan")

        pos = y_true == 1
        neg = ~pos
        if pos.sum() and neg.sum():
            auroc = _auroc(y_true, y_score)
            auprc = _average_precision(y_true, y_score)
        else:
            auroc = float("nan")
            auprc = float("nan")

        rows.append({
            "model_name": model_name,
            "auroc": auroc,
            "auprc": auprc,
            "topk_enrichment": topk_enrichment,
        })

    metrics_df = pd.DataFrame(rows, columns=["model_name", "auroc", "auprc", "topk_enrichment"])
    for column in ["auroc", "auprc", "topk_enrichment"]:
        metrics_df[column] = pd.to_numeric(metrics_df[column], errors="coerce").round(6)
    metrics_df.to_csv(output_dir / "model_metrics.tsv", sep="\t", index=False)

    ablation_df = metrics_df[["model_name", "auroc", "auprc", "topk_enrichment"]].copy()
    ablation_df.to_csv(output_dir / "ablation_matrix.tsv", sep="\t", index=False)

    return metrics_df, ablation_df


def main():
    parser = argparse.ArgumentParser(description="Benchmark VIPP and ablated models")
    parser.add_argument("--atlas-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    benchmark_models(args.atlas_dir, args.output_dir)


if __name__ == "__main__":
    main()
