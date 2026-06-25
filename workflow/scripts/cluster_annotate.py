"""Neighbors -> UMAP -> Leiden -> markers -> optional annotation."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

sm = snakemake  # noqa: F821
cfg = sm.config["cluster"]
adata = sc.read_h5ad(sm.input[0])
rep = adata.uns.get("use_rep", "X_pca")

sc.pp.neighbors(adata, n_neighbors=cfg["n_neighbors"], use_rep=rep)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=cfg["resolution"], key_added="leiden")
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True)

# top markers per cluster -> CSV
res = adata.uns["rank_genes_groups"]
groups = res["names"].dtype.names
mk = pd.concat(
    {g: pd.DataFrame({"gene": res["names"][g],
                      "lfc": res["logfoldchanges"][g],
                      "pval_adj": res["pvals_adj"][g]}).head(25) for g in groups},
    names=["cluster", "rank"]).reset_index(level=0)
mk.to_csv(sm.output["markers"], index=False)

# optional manual annotation map from config; otherwise label by cluster id
ann = cfg.get("annotation") or {}
adata.obs["cell_type"] = (
    adata.obs["leiden"].map(lambda c: ann.get(str(c), f"cluster_{c}"))
    if ann else "cluster_" + adata.obs["leiden"].astype(str)
).astype("category")

sc.pl.umap(adata, color=["leiden", "condition", "cell_type"], ncols=3, show=False)
plt.savefig(sm.output["umap"], dpi=150, bbox_inches="tight")
adata.write(sm.output["h5ad"])
print(f"clusters: {adata.obs['leiden'].nunique()}  cells: {adata.n_obs}")
