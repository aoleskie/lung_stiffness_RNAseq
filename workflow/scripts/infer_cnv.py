#!/usr/bin/env python3
"""
CNV-based malignant-cell calling (inferCNVpy) for the Kim NSCLC atlas.

WHY full-gene counts: CNV detection smooths expression across genes ordered along
each chromosome, so it needs dense genome-wide coverage. annotated.h5ad is the 2,000
HVG subset (wrong sampling), so we concatenate the per-sample resources/kim/*.h5ad
(all ~29,634 genes) and attach cell_type / leiden from annotated.h5ad by barcode.

Reference = immune cells (diploid karyotype). Malignant = aneuploid epithelial cells.
Normal-lung epithelium is the built-in sanity check: it should score near-diploid.

Outputs:
  - h5ad with obs['cnv_score'] and obs['cnv_status'] (malignant/normal/reference/other)
  - per-cell cnv_score table (csv)
  - diagnostic figures: CNV score distribution + score by cell type / condition
The malignant threshold is PROVISIONAL (data-driven from the epithelial bimodality);
inspect the distribution plot and finalize cnv_score_threshold in config.
"""
import re
import glob
import numpy as np
import pandas as pd
import matplotlib
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["LOKY_MAX_CPU_COUNT"] = "1"
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import infercnvpy as cnv

sm = snakemake  # noqa: F821
cfg = sm.config["cnv"]
ccfg = sm.config["cluster"]

# ---- 1. full-gene matrix: concatenate per-sample files ----
files = sorted(glob.glob(f"{sm.config['resources_dir']}/*.h5ad"))
print(f"concatenating {len(files)} per-sample files (full gene set) ...")
parts = [sc.read_h5ad(f) for f in files]
adata = ad.concat(parts, join="inner", index_unique="-", label="batch")
adata.var_names_make_unique()
print(f"  full-gene matrix: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

# ---- 2. attach cell_type + leiden from annotated.h5ad by barcode ----
ann = sc.read_h5ad(sm.input.annotated, backed="r").obs
strip = lambda s: re.sub(r"-\d+$", "", str(s))
# match annotated barcodes (which also carry a -batch suffix) to ours via stripped id
cols = ["cell_type", "leiden", "donor_id", "condition"]
have = [c for c in cols if c in ann.columns]
ann_by_bc = {strip(bc): row for bc, row in zip(ann.index, ann[have].astype(str).to_dict("records"))}
keys = [strip(b) for b in adata.obs_names]
matched = np.mean([k in ann_by_bc for k in keys])
# match rate ≈ QC survival rate: concat has all pre-QC cells, annotated.h5ad only survivors.
# Maynard QC dropped ~18% failed wells, so ~80% is correct here (Kim was ~99.6%).
assert matched > 0.6, f"only {matched:.1%} of cells matched annotated labels — too low, check the join"
for c in have:
    adata.obs[c] = [ann_by_bc.get(k, {}).get(c, "NA") for k in keys]
# keep only cells that survived QC/clustering (present in annotated)
adata = adata[[k in ann_by_bc for k in keys]].copy()
print(f"  matched labels for {adata.n_obs:,} cells ({matched:.1%})")

# ---- 3. normalize for CNV (infercnvpy wants log-normalized X) ----
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ---- 4. attach gene positions ----
pos = pd.read_csv(cfg["gene_order_file"], sep="\t").set_index("gene")
common = adata.var_names.intersection(pos.index)
print(f"  gene positions: {len(common):,}/{adata.n_vars:,} genes placed")
adata = adata[:, list(common)].copy()
adata.var["chromosome"] = "chr" + pos.loc[adata.var_names, "chromosome"].astype(str)  # infercnvpy wants 'chr' prefix
adata.var["start"] = pos.loc[adata.var_names, "start"].astype(int).values
adata.var["end"]   = pos.loc[adata.var_names, "end"].astype(int).values

# ---- 5. run inferCNV with immune cells as the diploid reference ----
# config names author labels; map to the coarse cell_type values actually present
ref_map = {"T lymphocytes": "T/NK", "NK cells": "T/NK",
           "Myeloid cells": "Myeloid", "B lymphocytes": "B"}
ref_cats = sorted({ref_map.get(c, c) for c in cfg["reference_cats"]}
                  & set(adata.obs["cell_type"].unique()))
print(f"  diploid reference categories: {ref_cats}")

# infercnvpy ignores n_jobs and uses tqdm.process_map (spawns subprocesses ->
# on Windows each re-imports this script and dies at `sm = snakemake`).
# Force it sequential by replacing process_map with a plain serial map.
import infercnvpy.tl._infercnv as _icnv
def _serial_map(fn, *iterables, **kw):
    return [fn(*args) for args in zip(*iterables)]
_icnv.process_map = _serial_map

cnv.tl.infercnv(
    adata,
    reference_key="cell_type",
    reference_cat=ref_cats,
    window_size=cfg.get("window_size", 100),
    step=cfg.get("step", 10),
)
cnv.tl.pca(adata)
cnv.pp.neighbors(adata)                    # builds the cnv_neighbors graph leiden needs
cnv.tl.leiden(adata)                       # CNV-space clusters -> obs['cnv_leiden']
cnv.tl.cnv_score(adata)                    # mean CNV magnitude per cnv_leiden cluster -> obs['cnv_score']# ---- 6. provisional malignant call ----
# threshold from config if set; else data-driven from epithelial bimodality
thr = cfg.get("cnv_score_threshold")
epi = adata.obs.loc[adata.obs["cell_type"] == cfg["malignant_candidate"].replace(" cells", ""), "cnv_score"]
if thr is None:
    # midpoint between reference median and epithelial 75th pct as a first cut
    ref_med = adata.obs.loc[adata.obs["cell_type"].isin(ref_cats), "cnv_score"].median()
    thr = float((ref_med + epi.quantile(0.75)) / 2) if len(epi) else float(adata.obs["cnv_score"].median())
print(f"  provisional CNV-score threshold: {thr:.4f}")

is_epi = adata.obs["cell_type"] == cfg["malignant_candidate"].replace(" cells", "")
status = np.where(adata.obs["cell_type"].isin(ref_cats), "reference",
          np.where(is_epi & (adata.obs["cnv_score"] > thr), "malignant",
          np.where(is_epi, "normal_epithelium", "other")))
adata.obs["cnv_status"] = status
print("\n  cnv_status counts:\n" + adata.obs["cnv_status"].value_counts().to_string())

# ---- 7. outputs ----
adata.obs[["cell_type", "leiden", "condition", "donor_id", "cnv_score", "cnv_status"]] \
     .to_csv(sm.output["scores"])
adata.write(sm.output["h5ad"])

# diagnostics
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for ct, sub in adata.obs.groupby("cell_type", observed=True):
    ax[0].hist(sub["cnv_score"], bins=60, alpha=.5, label=ct, density=True)
ax[0].axvline(thr, color="k", ls="--", lw=1, label=f"threshold {thr:.3f}")
ax[0].set_xlabel("CNV score"); ax[0].set_ylabel("density")
ax[0].set_title("CNV score by cell type"); ax[0].legend(fontsize=7)
epi_obs = adata.obs[is_epi]
for cond, sub in epi_obs.groupby("condition", observed=True):
    ax[1].hist(sub["cnv_score"], bins=50, alpha=.55, label=f"epithelial — {cond}", density=True)
ax[1].axvline(thr, color="k", ls="--", lw=1)
ax[1].set_xlabel("CNV score"); ax[1].set_title("Epithelial: tumor vs normal")
ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig(sm.output["fig"], dpi=150, bbox_inches="tight")
print(f"\nwrote {sm.output['h5ad']}, {sm.output['scores']}, {sm.output['fig']}")
