#!/usr/bin/env python3
"""Evaluate prespecified downstream model targets in Cosgrove soft-vs-stiff A549 RNA-seq.

The primary inferential unit is the biological replicate, not the gene. For each soft/stiff
sample, the script z-scores log2(TPM + 1) for each direct target across the six samples,
orients the z-score by the model-predicted direction, and averages the four oriented values.
It reports the stiff-minus-soft effect, a Welch 95% confidence interval, and the exact
one-sided label-permutation p-value (20 possible 3-vs-3 allocations).

Per-gene DESeq2 estimates and concordance remain descriptive. E2F1 and CCNE1 are reported as
a secondary downstream-support score and are not pooled into the primary test.
"""

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t


# gene -> (model prediction on stiff, which input drives it, the rate law)
TARGETS = {
    "CCND1": ("UP", "AKT,ERK", "CycD synthesis driven by AKT & ERK"),
    "MYC": ("UP", "AKT,ERK", "Myc synthesis driven by AKT & ERK"),
    "SKP2": ("UP", "FAK", "Skp2 synthesis driven by FAK"),
    "CDKN1A": ("DOWN", "FAK", "p21 synthesis repressed by FAK (FAK in denominator)"),
}

# Secondary genes downstream of CycD -> Rb -> E2F; not direct input targets.
SUPPORT = {
    "E2F1": ("UP", "downstream", "CycD->pRb phosphorylation releases E2F"),
    "CCNE1": ("UP", "downstream", "CycE is an E2F target downstream of CycD/Rb/E2F"),
}


def fmt(de, gene):
    row = de[de["symbol"] == gene]
    if row.empty:
        return None
    r = row.iloc[0]
    return {
        "symbol": gene,
        "baseMean": r.get("baseMean", np.nan),
        "log2FC": r["log2FoldChange"],
        "lfcSE": r.get("lfcSE", np.nan),
        "pvalue": r.get("pvalue", np.nan),
        "padj": r.get("padj", np.nan),
    }


def verdict(pred, lfc, padj):
    correct_dir = (lfc > 0 and pred == "UP") or (lfc < 0 and pred == "DOWN")
    sig = pd.notna(padj) and padj < 0.05
    if correct_dir and sig:
        return "correct direction + gene-level FDR < 0.05"
    if correct_dir:
        return "correct direction; gene-level FDR >= 0.05"
    return "opposite direction"


def gene_lfc_table(de):
    rows = []
    for set_name, genes in (("primary_direct", TARGETS), ("secondary_support", SUPPORT)):
        for gene, (pred, driver, law) in genes.items():
            values = fmt(de, gene)
            if values is None:
                rows.append(
                    {
                        "set": set_name,
                        "gene": gene,
                        "prediction": pred,
                        "driver": driver,
                        "rate_law_or_rationale": law,
                        "verdict": "not found",
                    }
                )
                continue
            rows.append(
                {
                    "set": set_name,
                    "gene": gene,
                    "prediction": pred,
                    "driver": driver,
                    "rate_law_or_rationale": law,
                    **values,
                    "verdict": verdict(pred, values["log2FC"], values["padj"]),
                }
            )
    return pd.DataFrame(rows)


def target_log_expression(tpm, metadata, de):
    sample_ids = metadata.index[metadata["condition"].isin(["soft", "stiff"])].tolist()
    condition_counts = metadata.loc[sample_ids, "condition"].value_counts().to_dict()
    if condition_counts != {"soft": 3, "stiff": 3}:
        raise ValueError(
            "expected exactly three soft and three stiff biological replicates; "
            f"found {condition_counts}"
        )
    missing_samples = sorted(set(sample_ids) - set(tpm.columns))
    if missing_samples:
        raise ValueError(f"TPM matrix is missing samples: {missing_samples}")

    expression = {}
    for gene in [*TARGETS, *SUPPORT]:
        ensembl_ids = de.index[de["symbol"] == gene].intersection(tpm.index)
        if ensembl_ids.empty:
            raise ValueError(f"no TPM row mapped to prespecified target {gene}")
        # Sum rare duplicate Ensembl-to-symbol mappings before transformation.
        expression[gene] = tpm.loc[ensembl_ids, sample_ids].sum(axis=0)
    return np.log2(pd.DataFrame(expression).loc[sample_ids] + 1), metadata.loc[sample_ids].copy()


def oriented_score(log_expression, target_dict):
    genes = list(target_dict)
    values = log_expression[genes]
    scales = values.std(axis=0, ddof=1)
    if (scales == 0).any():
        raise ValueError(f"constant expression for target(s): {scales.index[scales == 0].tolist()}")
    z = values.sub(values.mean(axis=0), axis=1).div(scales, axis=1)
    signs = pd.Series({gene: 1 if target_dict[gene][0] == "UP" else -1 for gene in genes})
    oriented = z.mul(signs, axis=1)
    return oriented.mean(axis=1), oriented


def exact_permutation_p(scores, conditions):
    scores = np.asarray(scores, dtype=float)
    conditions = np.asarray(conditions)
    n_stiff = int((conditions == "stiff").sum())
    observed = scores[conditions == "stiff"].mean() - scores[conditions == "soft"].mean()
    differences = []
    for stiff_indices in combinations(range(len(scores)), n_stiff):
        stiff_mask = np.zeros(len(scores), dtype=bool)
        stiff_mask[list(stiff_indices)] = True
        differences.append(scores[stiff_mask].mean() - scores[~stiff_mask].mean())
    pvalue = np.mean(np.asarray(differences) >= observed - 1e-12)
    return observed, float(pvalue), len(differences)


def effect_summary(scores, conditions, role):
    scores = np.asarray(scores, dtype=float)
    conditions = np.asarray(conditions)
    soft = scores[conditions == "soft"]
    stiff = scores[conditions == "stiff"]
    effect, permutation_p, n_permutations = exact_permutation_p(scores, conditions)

    soft_var = soft.var(ddof=1)
    stiff_var = stiff.var(ddof=1)
    se = np.sqrt(soft_var / len(soft) + stiff_var / len(stiff))
    numerator = (soft_var / len(soft) + stiff_var / len(stiff)) ** 2
    denominator = (soft_var / len(soft)) ** 2 / (len(soft) - 1)
    denominator += (stiff_var / len(stiff)) ** 2 / (len(stiff) - 1)
    welch_df = numerator / denominator
    critical = t.ppf(0.975, welch_df)

    pooled_sd = np.sqrt(
        ((len(soft) - 1) * soft_var + (len(stiff) - 1) * stiff_var)
        / (len(soft) + len(stiff) - 2)
    )
    cohen_d = effect / pooled_sd
    correction = 1 - 3 / (4 * (len(soft) + len(stiff) - 2) - 1)

    return {
        "role": role,
        "n_soft": len(soft),
        "n_stiff": len(stiff),
        "soft_mean": soft.mean(),
        "soft_sd": soft.std(ddof=1),
        "stiff_mean": stiff.mean(),
        "stiff_sd": stiff.std(ddof=1),
        "stiff_minus_soft": effect,
        "welch_95ci_low": effect - critical * se,
        "welch_95ci_high": effect + critical * se,
        "hedges_g": cohen_d * correction,
        "exact_permutation_p_one_sided": permutation_p if role == "primary_prespecified" else np.nan,
        "n_label_permutations": n_permutations if role == "primary_prespecified" else np.nan,
    }


def plot_primary_scores(sample_scores, stats, output_path):
    """Write a compact dependency-free SVG with all six biological replicates."""
    width, height = 720, 560
    left, right, top, bottom = 105, 680, 75, 390
    x_positions = {"soft": 255, "stiff": 525}
    colors = {"soft": "#4C78A8", "stiff": "#E45756"}
    values = sample_scores["primary_oriented_score"].to_numpy(dtype=float)
    value_range = values.max() - values.min()
    padding = max(0.18 * value_range, 0.15)
    y_min, y_max = values.min() - padding, values.max() + padding

    def y_pixel(value):
        return bottom - (value - y_min) / (y_max - y_min) * (bottom - top)

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#222}.tick{font-size:12px}'
        '.label{font-size:14px}.small{font-size:12px}.title{font-size:18px;font-weight:600}</style>',
        '<text x="360" y="30" text-anchor="middle" class="title">'
        'Prespecified direct targets across biological replicates</text>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="#222"/>',
        f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="#222"/>',
    ]

    for tick in np.linspace(y_min, y_max, 5):
        y = y_pixel(tick)
        svg.extend(
            [
                f'<line x1="{left - 6}" y1="{y:.1f}" x2="{left}" y2="{y:.1f}" stroke="#222"/>',
                f'<text x="{left - 10}" y="{y + 4:.1f}" text-anchor="end" class="tick">{tick:+.2f}</text>',
            ]
        )
    if y_min <= 0 <= y_max:
        zero_y = y_pixel(0)
        svg.append(
            f'<line x1="{left}" y1="{zero_y:.1f}" x2="{right}" y2="{zero_y:.1f}" '
            'stroke="#999" stroke-dasharray="4 4"/>'
        )

    for condition in ("soft", "stiff"):
        group = sample_scores[sample_scores["condition"] == condition].sort_values("rep")
        x_center = x_positions[condition]
        offsets = np.linspace(-28, 28, len(group))
        for offset, (_, row) in zip(offsets, group.iterrows()):
            x, y = x_center + offset, y_pixel(row["primary_oriented_score"])
            svg.extend(
                [
                    f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="{colors[condition]}" '
                    'stroke="white" stroke-width="1.2"/>',
                    f'<text x="{x:.1f}" y="{y - 12:.1f}" text-anchor="middle" class="small">'
                    f'rep{row["rep"]}</text>',
                ]
            )
        mean = group["primary_oriented_score"].mean()
        mean_y = y_pixel(mean)
        svg.extend(
            [
                f'<line x1="{x_center - 55}" y1="{mean_y:.1f}" x2="{x_center + 55}" '
                f'y2="{mean_y:.1f}" stroke="#111" stroke-width="3"/>',
                f'<text x="{x_center}" y="{bottom + 28}" text-anchor="middle" class="label">'
                f'{"Soft (1 kPa)" if condition == "soft" else "Stiff (50 kPa)"}</text>',
            ]
        )

    svg.extend(
        [
            '<text x="22" y="235" text-anchor="middle" class="label" '
            'transform="rotate(-90 22 235)">Primary oriented model score</text>',
            '<text x="22" y="235" text-anchor="middle" class="label" '
            'transform="rotate(-90 22 235) translate(0 18)">(mean oriented gene z-score)</text>',
            f'<text x="105" y="455" class="label">stiff - soft = {stats["stiff_minus_soft"]:+.2f}</text>',
            f'<text x="105" y="480" class="label">Welch 95% CI '
            f'[{stats["welch_95ci_low"]:+.2f}, {stats["welch_95ci_high"]:+.2f}]</text>',
            f'<text x="105" y="505" class="label">exact one-sided permutation '
            f'p = {stats["exact_permutation_p_one_sided"]:.3f}</text>',
            '<text x="105" y="535" class="small">Points are biological replicates; black bars are group means.</text>',
            '</svg>',
        ]
    )
    Path(output_path).write_text("\n".join(svg), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--deseq", default="results_cosgrove/deseq2_stiff_vs_soft.csv")
    parser.add_argument("--tpm", default="results_cosgrove/a549_tpm.csv")
    parser.add_argument("--samples", default="results_cosgrove/a549_samples.csv")
    parser.add_argument("--outdir", default="results_cosgrove")
    args = parser.parse_args()

    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    de = pd.read_csv(args.deseq, index_col=0)
    if "symbol" not in de.columns:
        raise SystemExit("expected a 'symbol' column in the DESeq2 CSV")
    if "log2FoldChange" not in de.columns:
        for alternative in ("log2FC", "LFC", "log2_fold_change"):
            if alternative in de.columns:
                de = de.rename(columns={alternative: "log2FoldChange"})
                break

    gene_table = gene_lfc_table(de)
    gene_table.to_csv(output_dir / "model_target_gene_lfc.csv", index=False)
    print("Per-gene DESeq2 estimates (descriptive; positive log2FC = up on stiff):")
    print(gene_table[["set", "gene", "prediction", "log2FC", "padj", "verdict"]].to_string(index=False))

    if not Path(args.tpm).exists():
        raise SystemExit(f"missing {args.tpm}; rerun ingest_cosgrove.py to write the TPM matrix")
    tpm = pd.read_csv(args.tpm, index_col=0)
    metadata = pd.read_csv(args.samples, index_col=0)
    log_expression, metadata = target_log_expression(tpm, metadata, de)

    primary_score, _ = oriented_score(log_expression, TARGETS)
    support_score, _ = oriented_score(log_expression, SUPPORT)
    combined_score, _ = oriented_score(log_expression, {**TARGETS, **SUPPORT})
    sample_scores = metadata[["condition", "rep", "gsm"]].copy()
    sample_scores.index.name = "sample_id"
    sample_scores["primary_oriented_score"] = primary_score
    sample_scores["secondary_support_score"] = support_score
    sample_scores["secondary_combined_score"] = combined_score
    sample_scores.to_csv(output_dir / "model_target_scores_by_sample.csv")

    summaries = []
    for score_name, role in (
        ("primary_oriented_score", "primary_prespecified"),
        ("secondary_support_score", "secondary_descriptive"),
        ("secondary_combined_score", "secondary_descriptive"),
    ):
        summaries.append(
            {
                "score": score_name,
                **effect_summary(sample_scores[score_name], sample_scores["condition"], role),
            }
        )
    summary = pd.DataFrame(summaries)
    summary.to_csv(output_dir / "model_target_score_summary.csv", index=False)

    primary_stats = summary.iloc[0].to_dict()
    plot_primary_scores(
        sample_scores.reset_index(),
        primary_stats,
        output_dir / "model_target_score_by_replicate.svg",
    )

    print("\nPrimary replicate-level result:")
    print(sample_scores[["condition", "rep", "primary_oriented_score"]].to_string())
    print(
        f"\nstiff - soft = {primary_stats['stiff_minus_soft']:+.3f}; "
        f"Welch 95% CI [{primary_stats['welch_95ci_low']:+.3f}, "
        f"{primary_stats['welch_95ci_high']:+.3f}]; "
        f"Hedges g = {primary_stats['hedges_g']:+.3f}; "
        f"exact one-sided permutation p = {primary_stats['exact_permutation_p_one_sided']:.4f} "
        f"({int(primary_stats['n_label_permutations'])} allocations)"
    )
    print("Secondary support and combined scores are descriptive and receive no inferential p-value.")


if __name__ == "__main__":
    main()
