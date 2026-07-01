#!/usr/bin/env python3
# Per-gene check of the model's DOWNSTREAM cell-cycle targets in the A549 stiff-vs-soft contrast.
#
# Usage: python check_model_targets.py --deseq results_cosgrove/deseq2_stiff_vs_soft.csv

import argparse
import numpy as np
import pandas as pd

# gene -> (model prediction on stiff, which input drives it, the rate law)
TARGETS = {
    "CCND1":  ("UP",   "AKT,ERK", "CycD synthesis driven by AKT & ERK"),
    "MYC":    ("UP",   "AKT,ERK", "Myc synthesis driven by AKT & ERK"),
    "SKP2":   ("UP",   "FAK",     "Skp2 synthesis driven by FAK"),
    "CDKN1A": ("DOWN", "FAK",     "p21 synthesis repressed by FAK (FAK in denominator)"),
}
# secondary / supporting cell-cycle genes the model also moves (E2F target program downstream of CycD->Rb->E2F)
SUPPORT = {
    "E2F1":   ("UP",  "downstream", "CycD->pRb phosphorylation releases E2F"),
    "CCNE1":  ("UP",  "downstream", "CycE, an E2F target, downstream of CycD/Rb/E2F"),
    "RB1":    ("ns",  "downstream", "Rb level itself not directly driven; phospho-state is (not in RNA)"),
    "CDKN2A": ("ns",  "n/a",       "p16 not an explicit input target in the model"),
}


def fmt(de, gene):
    row = de[de["symbol"] == gene]
    if row.empty:
        return None
    r = row.iloc[0]
    return dict(symbol=gene, baseMean=r.get("baseMean", np.nan),
                log2FC=r["log2FoldChange"], lfcSE=r.get("lfcSE", np.nan),
                pvalue=r.get("pvalue", np.nan), padj=r.get("padj", np.nan))


def verdict(pred, lfc, padj):
    if pred == "ns":
        return "—"
    correct_dir = (lfc > 0 and pred == "UP") or (lfc < 0 and pred == "DOWN")
    sig = (padj == padj) and (padj < 0.05)
    if correct_dir and sig:   return "✓ correct + sig"
    if correct_dir:           return "✓ correct dir (ns)"
    return "✗ WRONG direction"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deseq", default="results_cosgrove/deseq2_stiff_vs_soft.csv")
    args = ap.parse_args()

    de = pd.read_csv(args.deseq)
    # the ingest wrote a 'symbol' column; index is Ensembl. Be robust to header names.
    if "symbol" not in de.columns:
        raise SystemExit("expected a 'symbol' column in the DESeq2 csv (written by ingest_cosgrove.py)")
    if "log2FoldChange" not in de.columns:
        # some versions name it differently
        for alt in ["log2FC", "LFC", "log2_fold_change"]:
            if alt in de.columns:
                de = de.rename(columns={alt: "log2FoldChange"}); break

    print("Positive log2FC = UP on stiff (50 kPa) vs soft (1 kPa).  Stiff has higher AKT/ERK/FAK.\n")

    print("=" * 92)
    print("PRIMARY model targets (directly driven by an input in the ODE rate laws)")
    print("=" * 92)
    print(f"{'gene':8} {'pred':5} {'driver':8} {'log2FC':>8} {'padj':>10}   verdict")
    print("-" * 92)
    hits = 0; tested = 0
    for gene, (pred, driver, law) in TARGETS.items():
        v = fmt(de, gene)
        if v is None:
            print(f"{gene:8} {pred:5} {driver:8} {'NOT FOUND':>8}"); continue
        tested += 1
        vd = verdict(pred, v["log2FC"], v["padj"])
        if vd.startswith("✓"): hits += 1
        padj_s = f"{v['padj']:.2e}" if v["padj"] == v["padj"] else "NA"
        print(f"{gene:8} {pred:5} {driver:8} {v['log2FC']:>8.3f} {padj_s:>10}   {vd}")
        print(f"         └ {law}")
    print("-" * 92)
    print(f"PRIMARY: {hits}/{tested} genes in the model-predicted direction\n")

    print("=" * 92)
    print("SUPPORTING cell-cycle genes (downstream E2F program / not direct input targets)")
    print("=" * 92)
    print(f"{'gene':8} {'pred':5} {'log2FC':>8} {'padj':>10}   verdict")
    print("-" * 92)
    for gene, (pred, driver, law) in SUPPORT.items():
        v = fmt(de, gene)
        if v is None:
            print(f"{gene:8} {pred:5} {'NOT FOUND':>8}"); continue
        vd = verdict(pred, v["log2FC"], v["padj"])
        padj_s = f"{v['padj']:.2e}" if v["padj"] == v["padj"] else "NA"
        print(f"{gene:8} {pred:5} {v['log2FC']:>8.3f} {padj_s:>10}   {vd}")
        print(f"         └ {law}")

    # genome-wide context: where do these LFCs sit relative to all genes?
    alllfc = de["log2FoldChange"].dropna()
    print("\n" + "=" * 92)
    print("genome-wide context")
    print("=" * 92)
    print(f"all genes: median log2FC = {alllfc.median():.3f}, "
          f"IQR [{alllfc.quantile(.25):.3f}, {alllfc.quantile(.75):.3f}]")
    for gene in TARGETS:
        v = fmt(de, gene)
        if v is not None:
            pct = (alllfc < v["log2FC"]).mean() * 100
            print(f"  {gene:8} log2FC {v['log2FC']:>7.3f}  -> {pct:5.1f} percentile of all genes")

    # ---- COLLECTIVE TESTS: do the model's targets move coordinately as predicted? ----
    # No single gene is significant at 3v3, but the right question is whether the SET is shifted
    # in the model-predicted direction. Two complementary tests:
    #   (1) sign test  : how many of N predicted genes match their predicted direction (vs 50/50)
    #   (2) Wilcoxon   : orient each LFC by its prediction (LFC * sign), then ask if the oriented
    #                    target LFCs are shifted positive vs the genome-wide background (~0).
    from scipy.stats import binomtest, mannwhitneyu

    def oriented_set(target_dict):
        out = []
        for gene, (pred, _, _) in target_dict.items():
            if pred == "ns":
                continue
            v = fmt(de, gene)
            if v is None or v["log2FC"] != v["log2FC"]:   # skip missing / NaN LFC
                continue
            sign = +1 if pred == "UP" else -1
            out.append((gene, sign, v["log2FC"], sign * v["log2FC"]))
        return out

    print("\n" + "=" * 92)
    print("COLLECTIVE TEST — do the model's targets shift coordinately in the predicted direction?")
    print("=" * 92)

    for label, tdict in [("PRIMARY only (CCND1/MYC/SKP2/CDKN1A)", TARGETS),
                         ("PRIMARY + SUPPORT (adds E2F1/CCNE1)", {**TARGETS, **SUPPORT})]:
        ori = oriented_set(tdict)
        if not ori:
            continue
        n = len(ori)
        k = sum(1 for *_, o in ori if o > 0)           # concordant = oriented LFC positive
        sign_p = binomtest(k, n, 0.5, alternative="greater").pvalue
        oriented_vals = [o for *_, o in ori]
        # Wilcoxon: oriented target LFCs vs genome-wide background LFCs (one-sided, targets > background)
        u, wil_p = mannwhitneyu(oriented_vals, alllfc.values, alternative="greater")
        mean_ori = float(np.mean(oriented_vals))
        print(f"\n  {label}")
        print(f"    genes tested      : {n}  [{', '.join(g for g, *_ in ori)}]")
        print(f"    concordant        : {k}/{n} match predicted direction")
        print(f"    sign test         : p = {sign_p:.4f} (one-sided vs 50/50)")
        print(f"    mean oriented LFC : {mean_ori:+.3f}  (background median {alllfc.median():+.3f})")
        print(f"    Wilcoxon vs bkgd  : p = {wil_p:.4f} (one-sided, targets shifted positive)")

    print("\n  Interpretation: individual genes are underpowered at 3v3, but a coordinated shift of the")
    print("  model's downstream targets in the predicted direction is the claim the data can support.")


if __name__ == "__main__":
    main()