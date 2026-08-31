# Config and statistical-change rerun - 2026-08-31

- Branch: `main`
- Starting commit: `d3e2ae2`
- Workflow: Snakemake 9.23.1, Python 3.12.13
- Kim config MD5 after annotation update: `9a8c14ad314bfad8cfecbaf21d260351`

## Kim normalization and downstream outputs

The Kim normalization-through-signature chain was forced to rerun with the unified
`normalization:` schema. The run logged:

```text
normalization: method=lognorm, hvg_flavor=seurat_v3, hvg_layer=counts
highly variable genes retained: 2000
merged: 87832 cells, 2000 HVGs, 30 PCs
```

Regenerated outputs:

- `results_kim/integrate/merged.h5ad`
- `results_kim/integrate/integrated.h5ad`
- `results_kim/cluster/annotated.h5ad`
- `results_kim/cluster/markers.csv`
- `results_kim/figures/umap_clusters.png`
- `results_kim/biology/signature_scores.csv`
- `results_kim/figures/signature_umaps.png`

The normalization change altered the unsupervised partition, so the numeric cluster annotation
map in `config/config_kim.yaml` was updated before the final annotation/signature rerun. New
clusters were checked against the author cell-type labels retained in the integrated object and
against marker genes. Cluster 13 is a cycling immune cluster; cluster 24 is pDC-like/myeloid based
on `LILRA4`, `GZMB`, and `IRF7` markers.

| Composition | Before (`seurat`) | After (`seurat_v3`) |
|---|---:|---:|
| Cells | 87,832 | 87,832 |
| Leiden clusters | 24 | 26 |
| B | 5,665 | 5,752 |
| Endothelial | 2,004 | 2,021 |
| Epithelial | 10,871 | 10,948 |
| Fibroblast | 3,392 | 3,392 |
| Mast | 2,875 | 2,862 |
| Myeloid | 26,169 | 24,726 |
| Cycling immune | not separated | 1,514 |
| T/NK | 36,856 | 36,617 |

The broad lineage composition and the writeup's biological conclusions do not materially change.
The v2 writeup's CNV, CAF differential-expression, and donor-level analyses use the full-gene QC
objects and preserved author identities rather than the HVG partition, so they were not recomputed.

All five curated signatures completed using every configured gene. The optional PROGENy branch was
skipped because the installed Decoupler/Numba combination cannot compile its GSVA helper; this does
not affect the curated signature outputs used here.

Final Kim log: `.snakemake/log/2026-08-31T170338.163247.snakemake.log`

## Cosgrove statistical outputs

The Cosgrove ingest and replicate-level target analysis were rerun from the source RSEM files. The
regenerated primary result is:

- stiff minus soft oriented score: `+0.969`
- Welch 95% CI: `+0.026` to `+1.912`
- Hedges g: `+2.419`
- exact one-sided permutation p-value: `0.050` from 20 label allocations

`results_cosgrove/signature_summary.csv` now labels the descriptive column
`median_gene_level_padj`; the YAP/TAZ value remains `0.00166` and is not represented as a
pathway-level adjusted p-value.

All Cosgrove CSV, PNG, and SVG outputs were regenerated. The notebook
`notebooks/cosgrove_a549_writeup.ipynb` was then executed with the activated analysis environment;
it completed without cell errors and embeds the revised tables and figures.

## Outputs not rerun

Maynard and Tsukui were not rerun. Their configured HVG flavor is `seurat`, which is the same flavor
the pre-fix log-normalized branch hard-coded, so the config-key cleanup did not change their
effective normalization behavior. Dry-runs for all three scRNA-seq configs successfully built the
expected DAG; the requested Kim signature target was fully up to date after regeneration.
