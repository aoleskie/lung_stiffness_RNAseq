#!/usr/bin/env python3
# Cosgrove et al. (Science 2024; GSE243763) A549 lung adenocarcinoma on soft (1 kPa) vs stiff (50 kPa)
# polyacrylamide hydrogels. THE clean soft-vs-stiff cancer-cell contrast the patient datasets (Kim,
# Maynard) structurally could not provide. Tests the manuscript's central INPUT claim:
#   stiffness -> FAK/AKT up -> proliferation up  (and CycD-CDK4/6 as the therapeutic lever).
# Corrected prediction (per the model): STIFF should be HIGHER for AKT / proliferation / CDK4-6.
#
# Caveat carried in the writeup: A549 is KRAS-mutant, so this confirms the MECHANISM in lung
# adenocarcinoma, not the EGFR context specifically (PC9/HCC827 cover EGFR as literature).
#
# Input : per-sample RSEM gene results from GSE243763_RAW.tar (9 A549 files), e.g.
#         GSM9224457_A549.rnaseq.1kPa.rep1.star2.rsem.genes.results.txt.gz
#         columns: gene_id (versioned ENSG), expected_count, TPM, FPKM, ...
# Output: results_cosgrove/  -> counts/TPM matrices, DESeq2 table (50kPa vs 1kPa),
#                                signature scores, figure.
#
# Usage:
#   python ingest_cosgrove.py --raw-dir geo_tmp/raw --hff-table geo_tmp/GSE243763_SupplementaryTable2.csv.gz \
#       --outdir results_cosgrove
# (the HFF table is used only as an ENSG->symbol fallback map; optional)

import os, re, glob, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- the model's nodes, as gene sets (symbols) ---------------------------------------------------
# The manuscript routes mechanotransduction through FAK/AKT/ERK -> cell cycle (NOT via YAP/TAZ),
# so FAK/focal-adhesion and AKT/proliferation/CDK4-6 are the load-bearing signatures here.
SIGS = {
    "fak_focal_adhesion": ["PTK2", "SRC", "PXN", "BCAR1", "TLN1", "VCL", "ZYX", "TNS1"],   # FAK axis
    "pi3k_akt_axis":      ["AKT1", "PIK3CA", "PDPK1", "FOXO3", "GSK3B", "RPS6KB1"],
    "proliferation":      ["MKI67", "TOP2A", "PCNA", "CCNB1", "CDK1", "BIRC5"],
    "cdk46_cellcycle":    ["CCND1", "CDK4", "CDK6", "RB1", "CDKN1A", "CDKN2A", "E2F1"],     # the lever
    "yap_taz_targets":    ["CTGF", "CCN2", "CYR61", "CCN1", "ANKRD1", "AMOTL2", "THBS1", "CAV1", "TEAD1"],  # paper's readout
    "stiffening_ecm":     ["COL1A1", "COL1A2", "COL3A1", "FN1", "ELN", "FBN1"],             # control
}

COND_FROM_NAME = [("1kPa", "soft"), ("50kPa", "stiff"), ("TCP", "tcp")]


def parse_sample(path):
    """From the RSEM filename pull (sample_id, condition, rep). e.g. ...A549.rnaseq.1kPa.rep1..."""
    base = os.path.basename(path)
    gsm = base.split("_")[0]
    cond = next((c for tok, c in COND_FROM_NAME if tok.lower() in base.lower()), "NA")
    m = re.search(r"rep(\d+)", base, re.I)
    rep = m.group(1) if m else "1"
    return f"{cond}_rep{rep}", cond, rep, gsm


def load_counts(raw_dir):
    """Assemble per-sample RSEM expected_count + TPM into matrices (genes x samples)."""
    files = sorted(glob.glob(os.path.join(raw_dir, "*A549*genes.results*")))
    if not files:
        files = sorted(glob.glob(os.path.join(raw_dir, "*9224*genes.results*")))
    assert files, f"no A549 RSEM files found in {raw_dir}"
    counts, tpm, meta = {}, {}, []
    for f in files:
        sid, cond, rep, gsm = parse_sample(f)
        df = pd.read_csv(f, sep="\t", index_col=0)
        # strip Ensembl version suffix (ENSG00000000003.13 -> ENSG00000000003)
        df.index = df.index.str.replace(r"\.\d+$", "", regex=True)
        counts[sid] = df["expected_count"].round().astype(int)   # RSEM estimates -> integers for DESeq2
        tpm[sid]    = df["TPM"]
        meta.append({"sample_id": sid, "condition": cond, "rep": rep, "gsm": gsm, "file": os.path.basename(f)})
    C = pd.DataFrame(counts); T = pd.DataFrame(tpm)
    M = pd.DataFrame(meta).set_index("sample_id").loc[C.columns]
    print(f"  loaded {C.shape[1]} samples x {C.shape[0]:,} genes")
    print("  conditions:", dict(M["condition"].value_counts()))
    return C, T, M


def ensembl_symbol_map(genes_ensembl, hff_table):
    """ENSG->symbol. Prefer pyensembl (already installed); fall back to the HFF table's Gene Name col."""
    mapping = {}
    try:
        from pyensembl import EnsemblRelease
        data = EnsemblRelease(100)
        for g in genes_ensembl:
            try:
                mapping[g] = data.gene_by_id(g).gene_name
            except Exception:
                pass
        if mapping:
            print(f"  pyensembl mapped {len(mapping):,}/{len(genes_ensembl):,} ENSG->symbol")
            return mapping
    except Exception as e:
        print(f"  pyensembl unavailable ({e}); using HFF-table map")
    if hff_table and os.path.exists(hff_table):
        h = pd.read_csv(hff_table, index_col=0)
        h.index = h.index.str.replace(r"\.\d+$", "", regex=True)
        if "Gene Name" in h.columns:
            mapping = h["Gene Name"].dropna().to_dict()
            print(f"  HFF-table mapped {len(mapping):,} ENSG->symbol")
    return mapping


def score_signature(expr_log, sym_index, genes):
    """Mean z-scored log-expression across the signature's present genes (per sample).
    expr_log: genes(ensembl) x samples log2 matrix; sym_index: ensembl->symbol Series."""
    present_ens = sym_index.index[sym_index.isin(genes)]
    present_ens = [e for e in present_ens if e in expr_log.index]
    if not present_ens:
        return None, []
    sub = expr_log.loc[present_ens]
    # z-score each gene across samples, then average -> signature score per sample
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    score = z.mean(axis=0)
    syms = sorted(set(sym_index.loc[present_ens]))
    return score, syms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", default="geo_tmp/raw")
    ap.add_argument("--hff-table", default="geo_tmp/GSE243763_SupplementaryTable2.csv.gz")
    ap.add_argument("--outdir", default="results_cosgrove")
    ap.add_argument("--run-deseq", action="store_true", default=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print("1) assembling A549 RSEM counts ...")
    C, T, M = load_counts(args.raw_dir)
    C.to_csv(os.path.join(args.outdir, "a549_counts.csv"))
    T.to_csv(os.path.join(args.outdir, "a549_tpm.csv"))
    M.to_csv(os.path.join(args.outdir, "a549_samples.csv"))

    print("2) ENSG -> symbol map ...")
    smap = ensembl_symbol_map(C.index.tolist(), args.hff_table)
    sym = pd.Series(smap).reindex(C.index)   # ensembl-indexed symbol Series (NaN where unmapped)

    # ---- soft vs stiff contrast: 1kPa vs 50kPa (drop TCP from the primary test) ----
    keep = M.index[M["condition"].isin(["soft", "stiff"])]
    Csub = C[keep]; Msub = M.loc[keep]
    print(f"\n3) DESeq2 on soft(1kPa) vs stiff(50kPa): {dict(Msub['condition'].value_counts())}")

    de = None
    if args.run_deseq:
        try:
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats
            # genes x samples -> samples x genes; filter near-zero genes
            counts_t = Csub.T
            counts_t = counts_t.loc[:, counts_t.sum(axis=0) >= 10]
            meta_d = Msub[["condition"]].copy()
            dds = DeseqDataSet(counts=counts_t, metadata=meta_d, design="~condition")
            dds.deseq2()
            # contrast: stiff vs soft  -> positive LFC = UP on stiff (the model's predicted direction)
            st = DeseqStats(dds, contrast=["condition", "stiff", "soft"])
            st.summary()
            de = st.results_df.copy()
            de["symbol"] = sym.reindex(de.index).values
            de.to_csv(os.path.join(args.outdir, "deseq2_stiff_vs_soft.csv"))
            print(f"  DESeq2 done: {de.shape[0]:,} genes; "
                  f"{(de['padj']<0.05).sum():,} sig at padj<0.05")
        except Exception as e:
            print(f"  [DESeq2 skipped: {e}] — falling back to TPM means")

    # ---- signature scoring on log2 TPM (all conditions, for the soft->stiff->TCP trend) ----
    print("\n4) scoring model-node signatures (log2 TPM, z-scored) ...")
    logT = np.log2(T + 1)
    rows = []
    score_table = {}
    for name, genes in SIGS.items():
        score, syms = score_signature(logT, sym, genes)
        if score is None:
            print(f"  {name}: 0 genes present — skipped"); continue
        score_table[name] = score
        soft = score[M.index[M.condition == "soft"]].mean()
        stiff = score[M.index[M.condition == "stiff"]].mean()
        # DESeq2 mean LFC across the signature genes (stiff vs soft), if available
        if de is not None:
            sig_lfc = de.loc[de["symbol"].isin(genes), "log2FoldChange"].mean()
            # Descriptive summary of constituent gene-level tests; this is not a
            # pathway- or signature-level adjusted p-value.
            median_gene_level_padj = de.loc[de["symbol"].isin(genes), "padj"].median()
        else:
            sig_lfc, median_gene_level_padj = np.nan, np.nan
        rows.append({"signature": name, "n_genes": len(syms),
                     "soft_score": round(soft, 3), "stiff_score": round(stiff, 3),
                     "stiff_minus_soft": round(stiff - soft, 3),
                     "mean_LFC_stiff_vs_soft": round(sig_lfc, 3) if sig_lfc == sig_lfc else np.nan,
                     "median_gene_level_padj": (f"{median_gene_level_padj:.2e}"
                                                  if median_gene_level_padj == median_gene_level_padj
                                                  else "NA"),
                     "predicted": "UP on stiff" if name in
                         ("fak_focal_adhesion","pi3k_akt_axis","proliferation","cdk46_cellcycle") else "(readout/ctrl)"})
    summary = pd.DataFrame(rows).set_index("signature")
    summary.to_csv(os.path.join(args.outdir, "signature_summary.csv"))
    print("\n=== signature scores: soft (1kPa) vs stiff (50kPa) ===")
    print(summary.to_string())

    # ---- figure: per-signature soft vs stiff (z-scored), with the model's predicted direction ----
    order = ["soft", "stiff", "tcp"]
    present_conds = [c for c in order if c in set(M.condition)]
    names = list(score_table)
    fig, axes = plt.subplots(1, len(names), figsize=(2.7*len(names), 4.0), sharey=True)
    if len(names) == 1: axes = [axes]
    for ax, name in zip(axes, names):
        sc_ = score_table[name]
        data = [sc_[M.index[M.condition == c]].values for c in present_conds]
        xs = np.arange(len(present_conds))
        for x, vals in zip(xs, data):
            ax.scatter([x]*len(vals), vals, s=28, alpha=.8)
            ax.hlines(np.mean(vals), x-0.2, x+0.2, color="k", lw=2)
        ax.set_xticks(xs); ax.set_xticklabels([c.replace("soft","1kPa").replace("stiff","50kPa").upper()
                                               for c in present_conds], fontsize=8)
        pred = name in ("fak_focal_adhesion","pi3k_akt_axis","proliferation","cdk46_cellcycle")
        ax.set_title(f"{name}\n{'model: UP on stiff' if pred else 'readout/ctrl'}", fontsize=8)
        ax.axhline(0, color="grey", lw=.6, ls=":")
    axes[0].set_ylabel("signature score (z-scored log2 TPM)")
    fig.suptitle("A549 soft (1 kPa) vs stiff (50 kPa) — model-node signatures", y=1.03)
    plt.tight_layout(); plt.savefig(os.path.join(args.outdir, "figures_signatures.png"),
                                    dpi=150, bbox_inches="tight")
    print(f"\nwrote {args.outdir}/ (counts/TPM, DESeq2 table, signature_summary.csv, figures_signatures.png)")


if __name__ == "__main__":
    main()
