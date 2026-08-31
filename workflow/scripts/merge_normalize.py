"""Concatenate QC'd samples, normalize, pick HVGs, PCA -> merged.h5ad."""
import scanpy as sc
import anndata as ad
import numpy as np

sm = snakemake  # noqa: F821
cfg = sm.config["normalization"]
hvg_flavor = cfg["hvg_flavor"]

adata = ad.concat([sc.read_h5ad(p) for p in sm.input],
                  join="outer", label="batch", index_unique="-")
adata.layers["counts"] = adata.X.copy()  # X is raw counts coming out of qc.py

# log1p CP10k (default). For analytic Pearson residuals set normalization.method: pearson.
if cfg["method"] == "pearson":
    sc.experimental.pp.normalize_pearson_residuals(adata)
else:
    sc.pp.normalize_total(adata, target_sum=cfg["target_sum"])
    sc.pp.log1p(adata)
    adata.raw = adata

# seurat_v3 flavors model raw counts; dispersion-based flavors model log-normalized X.
hvg_layer = "counts" if hvg_flavor.startswith("seurat_v3") else None
print(f"normalization: method={cfg['method']}, hvg_flavor={hvg_flavor}, hvg_layer={hvg_layer or 'X'}")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=cfg["n_top_genes"],
    flavor=hvg_flavor,
    layer=hvg_layer,
)

adata = adata[:, adata.var["highly_variable"]].copy()
print(f"highly variable genes retained: {adata.n_vars}")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=cfg["n_pcs"])
adata.write(sm.output[0])
print(f"merged: {adata.n_obs} cells, {adata.n_vars} HVGs, {cfg['n_pcs']} PCs")
