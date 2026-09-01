#!/usr/bin/env python3
"""
Fibroblast-restricted pseudobulk DE (tumor vs normal CAFs) — Kim GSE131907.

Separates broad tissue-composition effects from expression differences within the fibroblast
compartment, using the same design as v1. A score
difference between tumor and normal fibroblasts has two possible causes and cannot
distinguish them on its own: tumors may simply hold more pathologic fibroblasts, or the
fibroblasts present may each express the program harder. Running the contrast twice does
distinguish them:
  - whole-tissue DE : tumor vs normal across all cells per donor, which deliberately
                      includes cell-type composition shifts ("tumors have more fibroblasts")
  - fibroblast DE   : tumor vs normal within fibroblasts only per donor, which holds
                      composition fixed and leaves what the compartment is doing

Reading the two together:
  * a gene up in whole-tissue that drops out under restriction was compositional, driven
    by cell-type proportions rather than by CAF activation.
  * a gene up in both, and sharper (larger LFC or more significant) in the fibroblast
    analysis, supports increased expression within the fibroblast compartment. This restriction
    does not control for differences among fibroblast states.

Reads raw-count full-gene per-donor files from results_kim/qc/*.h5ad, never the 2k-HVG
annotated.h5ad — CNV and DE both need genome-wide coverage the HVG subset does not have.
Pseudobulk by hand in pandas; round to int for PyDESeq2.

Run (in an env with pydeseq2 — e.g. pip install pydeseq2 into the scanpy env):
    python workflow/scripts/caf_de.py
"""
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

MIN_CELLS = 20          # per-donor min cells to enter a pseudobulk (stable summation)
MIN_GENE_COUNTS = 10    # drop genes below this total across the pseudobulk
GROUP_A, GROUP_B = "tumor", "normal"
# the v1 stiffening / pathologic-CAF program + a couple of activation controls
FOCUS = ["CTHRC1", "POSTN", "COMP", "COL1A1", "COL3A1", "COL1A2",
         "FN1", "ELN", "FBN1", "LOX", "LOXL1", "LOXL2", "ACTA2", "TAGLN"]
FIB_LABEL = "Fibroblasts"   # author label in the QC files (NOT the coarse "Fibroblast" from annotated.h5ad)

OUT_DIR = "results_kim/de"
FIG_DIR = "results_kim/figures"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)


def build_pseudobulk(files, restrict_celltype=None):
    """Sum raw counts per donor (optionally within one cell type).
    Returns (counts donors x genes int, meta donors x [condition])."""
    per_donor, cond, gene_index = {}, {}, None
    for f in files:
        a = sc.read_h5ad(f)
        if gene_index is None:
            gene_index = a.var_names
        if restrict_celltype is not None:
            a = a[a.obs["cell_type"] == restrict_celltype]
        if a.n_obs < MIN_CELLS:
            continue
        donor = str(a.obs["donor_id"].iloc[0])
        per_donor[donor] = np.asarray(a.X.sum(axis=0)).ravel()
        cond[donor] = str(a.obs["condition"].iloc[0])
    counts = pd.DataFrame(per_donor, index=gene_index).T          # donors x genes
    counts = counts.round().astype(int)
    counts = counts.loc[:, counts.sum(axis=0) >= MIN_GENE_COUNTS]  # gene filter
    meta = pd.DataFrame({"condition": pd.Series(cond)}).loc[counts.index]
    return counts, meta


def run_deseq(counts, meta):
    dds = DeseqDataSet(counts=counts, metadata=meta, design="~condition", quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", GROUP_A, GROUP_B], quiet=True)
    st.summary()
    return st.results_df


def main():
    files = sorted(glob.glob("results_kim/qc/LUNG_*.h5ad"))
    print(f"{len(files)} per-donor QC files\n")

    # sanity: confirm qc X is raw integer counts (pseudobulk needs raw, not normalized)
    chk = sc.read_h5ad(files[0])
    xmax = chk.X.max()
    print(f"X sanity (first file): max={xmax:.1f}  -> "
          f"{'looks like raw counts' if xmax > 30 else 'WARNING: may be normalized!'}\n")

    print("=== whole-tissue pseudobulk (all cells) ===")
    ct_all, meta_all = build_pseudobulk(files, restrict_celltype=None)
    print(f"{ct_all.shape[0]} donors x {ct_all.shape[1]} genes | "
          f"{dict(meta_all.condition.value_counts())}")
    res_all = run_deseq(ct_all, meta_all)

    print("\n=== fibroblast-restricted pseudobulk ===")
    ct_fib, meta_fib = build_pseudobulk(files, restrict_celltype=FIB_LABEL)
    print(f"{ct_fib.shape[0]} donors x {ct_fib.shape[1]} genes | "
          f"{dict(meta_fib.condition.value_counts())}")
    res_fib = run_deseq(ct_fib, meta_fib)

    # ---- whole-tissue vs fibroblast-restricted: focus genes side by side ----
    comp = pd.DataFrame({
        "LFC_whole":  res_all["log2FoldChange"].reindex(FOCUS),
        "padj_whole": res_all["padj"].reindex(FOCUS),
        "LFC_fib":    res_fib["log2FoldChange"].reindex(FOCUS),
        "padj_fib":   res_fib["padj"].reindex(FOCUS),
    })
    comp["sharpens_in_fib"] = (comp["LFC_fib"] > comp["LFC_whole"]) & (comp["padj_fib"] < 0.05)
    print("\n=== whole-tissue vs fibroblast-restricted DE (focus genes) ===")
    print("  LFC>0 = up in tumor; 'sharpens_in_fib' = stronger within the fibroblast compartment")
    print(comp.round(3).to_string())

    res_all.to_csv(f"{OUT_DIR}/caf_wholetissue_de.csv")
    res_fib.to_csv(f"{OUT_DIR}/caf_fibroblast_de.csv")
    comp.to_csv(f"{OUT_DIR}/caf_focus_comparison.csv")

    # ---- volcano of the fibroblast-restricted DE, focus genes highlighted ----
    rf = res_fib.dropna(subset=["padj"]).copy()
    rf["nlp"] = -np.log10(rf["padj"].clip(lower=1e-300))
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.4))

    ax[0].scatter(rf["log2FoldChange"], rf["nlp"], s=7, c="#cccccc", alpha=.5)
    for g in FOCUS:
        if g in rf.index:
            r = rf.loc[g]
            ax[0].scatter(r["log2FoldChange"], r["nlp"], s=70, c="#2a6f74", zorder=3, edgecolor="white")
            ax[0].annotate(g, (r["log2FoldChange"], r["nlp"]), fontsize=8,
                           xytext=(4, 2), textcoords="offset points")
    ax[0].axhline(-np.log10(0.05), ls="--", lw=.7, c="#888")
    ax[0].axvline(0, lw=.5, c="#888")
    ax[0].set_xlabel("log2 FC (tumor / normal CAF)"); ax[0].set_ylabel("-log10 padj")
    ax[0].set_title("Fibroblast-restricted DE: tumor vs normal CAFs")

    # whole-tissue vs fibroblast LFC for focus genes — the composition→intrinsic shift
    c = comp.dropna(subset=["LFC_whole", "LFC_fib"])
    ax[1].scatter(c["LFC_whole"], c["LFC_fib"], s=70, c="#2a6f74", zorder=3)
    for g, r in c.iterrows():
        ax[1].annotate(g, (r["LFC_whole"], r["LFC_fib"]), fontsize=8,
                       xytext=(4, 2), textcoords="offset points")
    lim = np.nanmax(np.abs([c["LFC_whole"], c["LFC_fib"]])) * 1.15
    ax[1].plot([-lim, lim], [-lim, lim], ls="--", lw=.7, c="#888")  # y=x: above => sharpens in fib
    ax[1].axhline(0, lw=.4, c="#bbb"); ax[1].axvline(0, lw=.4, c="#bbb")
    ax[1].set_xlabel("LFC whole-tissue"); ax[1].set_ylabel("LFC fibroblast-restricted")
    ax[1].set_title("Above y=x = sharper in fibroblast-restricted DE\n(within-compartment evidence)")
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/caf_de_volcano.png", dpi=150, bbox_inches="tight")
    print(f"\nwrote {FIG_DIR}/caf_de_volcano.png and DE CSVs to {OUT_DIR}/")


if __name__ == "__main__":
    main()
