"""Per-sample QC: load -> metrics -> filter -> doublets -> filtered .h5ad.
Run by Snakemake (the `snakemake` object is injected)."""
import os
import numpy as np
import scanpy as sc

sm = snakemake  # noqa: F821  (injected by Snakemake)
cfg = sm.config["qc"]


def load_counts(path):
    if path.endswith(".h5ad"):
        return sc.read_h5ad(path)
    if path.endswith(".h5"):
        return sc.read_10x_h5(path)
    if os.path.isdir(path):
        return sc.read_10x_mtx(path, var_names="gene_symbols", cache=False)
    raise ValueError(f"Unrecognized matrix input: {path}")


adata = load_counts(sm.input[0])
adata.var_names_make_unique()

# carry sample-level metadata through the whole pipeline
adata.obs["sample_id"] = sm.params.sample
adata.obs["condition"] = sm.params.condition
adata.obs["donor_id"] = sm.params.donor

# QC metrics
adata.var["mito"] = adata.var_names.str.startswith(cfg["mito_prefix"])
sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], percent_top=None,
                           log1p=False, inplace=True)

# cell filtering
if cfg["method"] == "mad":
    def mad_lo(x, n):
        med = np.median(x)
        mad = np.median(np.abs(x - med)) * 1.4826
        return med - n * mad
    n = cfg["n_mads"]
    g = np.log1p(adata.obs["n_genes_by_counts"].values)
    c = np.log1p(adata.obs["total_counts"].values)
    keep = (g > mad_lo(g, n)) & (c > mad_lo(c, n)) & \
           (adata.obs["pct_counts_mito"].values < cfg["max_pct_mito"])
else:
    keep = (adata.obs["n_genes_by_counts"] >= cfg["min_genes"]) & \
           (adata.obs["total_counts"] >= cfg["min_counts"]) & \
           (adata.obs["pct_counts_mito"] < cfg["max_pct_mito"])
adata = adata[np.asarray(keep)].copy()

# doublets
if cfg.get("doublet_method") == "scrublet":
    try:
        sc.pp.scrublet(adata)
        adata = adata[~adata.obs["predicted_doublet"]].copy()
    except Exception as e:  # version/skip-safe
        print("scrublet skipped:", e)

# stash raw counts for later (pseudobulk + seurat_v3 HVG)
adata.layers["counts"] = adata.X.copy()
adata.write(sm.output[0])
print(f"{sm.params.sample}: {adata.n_obs} cells x {adata.n_vars} genes after QC")
