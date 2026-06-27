"""Batch integration (Harmony by default). Records which rep to cluster on."""
import scanpy as sc
import numpy as np

sm = snakemake  # noqa: F821
cfg = sm.config["integration"]
adata = sc.read_h5ad(sm.input[0])

if cfg["method"] == "harmony":
    import harmonypy
    ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, [cfg["batch_key"]])
    Z = np.asarray(ho.Z_corr).T                      # -> (cells, PCs)
    if Z.shape[0] != adata.n_obs:                    # guard against orientation flips
        Z = Z.T
    assert Z.shape == adata.obsm["X_pca"].shape, f"harmony shape {Z.shape} != PCA {adata.obsm['X_pca'].shape}"
    adata.obsm["X_pca_harmony"] = Z
    adata.uns["use_rep"] = "X_pca_harmony"
else:
    adata.uns["use_rep"] = "X_pca"

adata.write(sm.output[0])
print(f"integration={cfg['method']}  use_rep={adata.uns['use_rep']}  embedding shape={adata.obsm.get('X_pca_harmony', adata.obsm['X_pca']).shape}")
