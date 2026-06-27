"""Concatenate QC'd samples, normalize, pick HVGs, PCA -> merged.h5ad."""
import scanpy as sc
import anndata as ad
import numpy as np

sm = snakemake  # noqa: F821
cfg = sm.config["norm"]

adata = ad.concat([sc.read_h5ad(p) for p in sm.input],
                  join="outer", label="batch", index_unique="-")
adata.layers["counts"] = adata.X.copy()  # X is raw counts coming out of qc.py

# log1p CP10k (default). For analytic Pearson residuals set norm.method: pearson.
if cfg["method"] == "pearson":
    sc.experimental.pp.normalize_pearson_residuals(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=cfg["n_top_genes"],
                                flavor="seurat_v3", layer="counts")
else:
    sc.pp.normalize_total(adata, target_sum=cfg["target_sum"])
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=cfg["n_top_genes"], flavor="seurat")

adata = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=cfg["n_pcs"])
adata.write(sm.output[0])
print(f"merged: {adata.n_obs} cells, {adata.n_vars} HVGs, {cfg['n_pcs']} PCs")
