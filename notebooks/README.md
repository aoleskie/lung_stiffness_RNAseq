# Notebooks

Interactive exploration that sits alongside the Snakemake pipeline. The pipeline
produces the reproducible objects (`results/.../*.h5ad`); notebooks are where you
poke at them cell-by-cell.

Suggested per-stage notebooks (keep a clean + an executed copy of each, per the
project convention):

- `01_v1_explore.ipynb`  - load a QC'd object, sanity-check metrics, eyeball UMAP.
- `02_v1_fibroblasts.ipynb` - subcluster mesenchyme, find the CTHRC1+ subset, stiffening scores.
- `03_v2_tumor.ipynb` - annotate lineages, CNV-call malignant cells, model-node signatures.

Open any `.h5ad` with:  `import scanpy as sc; adata = sc.read_h5ad("results/cluster/annotated.h5ad")`
