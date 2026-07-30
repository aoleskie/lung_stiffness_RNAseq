# Notebooks

The Snakemake pipeline produces the reproducible objects: `results_<name>/.../*.h5ad` and the CSV
tables that sit beside them. These notebooks are where those objects get interrogated cell by cell
and turned into an argument.

Four writeups, one per dataset, each self-contained and each reading only from its own
`results_<name>/`:

- `v1_writeup.ipynb` — Tsukui fibrotic lung. Recovers the CTHRC1⁺ pathologic fibroblast, then
  separates its expansion in number from its activation per cell.
- `v2_kim_writeup.ipynb` — Kim NSCLC. CNV-based malignant calling, the sequencing-depth confound
  that caps what the calling can support, and the CAF bridge showing tumor stroma rebuilding the v1
  niche.
- `v3_maynard_writeup.ipynb` — Maynard treatment trajectory. Persister-state reprogramming, plus a
  long look at why a treatment axis cannot test a stiffness prediction no matter how it is analyzed.
- `cosgrove_a549_writeup.ipynb` — Cosgrove A549 hydrogels. The one controlled soft-versus-stiff
  contrast in the project, and the only place the model's central prediction gets a fair test.

`01_v1_explore.ipynb` is a scratch notebook rather than a writeup: load an object, check the metrics
look sane, eyeball the UMAP. It earns its keep right after a rule finishes, when the only question
is whether the output is worth building on.

Each writeup walks up from the working directory to the repository root on startup, so they run from
anywhere inside the repo. To open an object by hand:

```python
import scanpy as sc
adata = sc.read_h5ad("results_kim/cluster/annotated.h5ad")
```

One thing worth knowing before poking around: `annotated.h5ad` carries only the 2,000 HVGs. If a
gene you expected is missing, that is why, and the full-gene counts live in the per-sample
`resources/*.h5ad` files.
