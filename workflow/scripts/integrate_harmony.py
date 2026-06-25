"""Batch integration (Harmony by default). Records which rep to cluster on."""
import scanpy as sc

sm = snakemake  # noqa: F821
cfg = sm.config["integration"]
adata = sc.read_h5ad(sm.input[0])

if cfg["method"] == "harmony":
    sc.external.pp.harmony_integrate(adata, key=cfg["batch_key"])
    adata.uns["use_rep"] = "X_pca_harmony"
else:
    adata.uns["use_rep"] = "X_pca"

adata.write(sm.output[0])
print(f"integration={cfg['method']}  use_rep={adata.uns['use_rep']}")
