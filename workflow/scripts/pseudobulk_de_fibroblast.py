#!/usr/bin/env python3
"""
Fibroblast-restricted pseudobulk DE (fibrotic vs normal).

Removes the immune/epithelial composition confound from the whole-tissue DE by
summing ONLY fibroblast-lineage cells per donor. Cell-type labels come from the
annotated (HVG-subset) object; counts come from the full-gene qc files. The two
are joined by cell barcode, undoing the `-{batch}` suffix that anndata.concat
added during merge_normalize.

NOTE: this restricts to the fibroblast *lineage* (all clusters whose label
contains "fibroblast"). It removes the immune/epithelial confound but still
contains the within-fibroblast pathologic<->homeostatic shift, so it is "what the
fibroblast compartment does", not strictly single-cell-intrinsic. For the latter,
compare only the pathologic cluster fibrotic-vs-normal (the signature analysis
already showed ~8x there).
"""
import re
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import issparse

sm = snakemake  # noqa: F821
cfg = sm.config["de"]
sk, ck = cfg["sample_key"], cfg["condition_key"]
LINEAGE = "fibroblast"   # substring matched against cell_type labels

# ---- cell_type lookup keyed by ORIGINAL barcode (strip the -batch suffix) ----
ann = sc.read_h5ad(sm.input.annotated).obs
if "cell_type" not in ann.columns:
    raise ValueError("annotated object has no 'cell_type' column; set the annotation map in config and re-run cluster_annotate")
strip = lambda s: re.sub(r"-\d+$", "", str(s))      # remove single trailing -<int> (the batch tag)
ct_map = dict(zip((strip(x) for x in ann.index), ann["cell_type"].astype(str)))

# ---- per-donor pseudobulk over fibroblast cells, from FULL-gene qc counts ----
profiles, conds = {}, {}
var_names = None
matched = unmatched = 0
for f in sm.input.qc:                                # one qc file = one donor, all genes
    a = sc.read_h5ad(f)
    donor = str(a.obs[sk].iloc[0])
    labels = np.array([ct_map.get(bc) for bc in a.obs_names], dtype=object)
    matched += int(pd.notna(labels).sum())
    unmatched += int(pd.isna(labels).sum())
    is_fib = pd.Series(labels).str.contains(LINEAGE, case=False, na=False).values
    if is_fib.sum() == 0:
        print(f"WARNING: {donor} matched 0 fibroblasts — skipping")
        continue
    c = a.layers["counts"] if "counts" in a.layers else a.X
    c = c.toarray() if issparse(c) else np.asarray(c)
    profiles[donor] = c[is_fib].sum(axis=0).ravel()
    conds[donor] = str(a.obs[ck].iloc[0])
    var_names = a.var_names
    print(f"  {donor}: {is_fib.sum()} fibroblasts / {a.n_obs} cells  [{conds[donor]}]")

rate = matched / max(matched + unmatched, 1)
assert rate > 0.95, (f"only {rate:.1%} of qc cells matched annotated labels — "
                     f"barcode join is wrong (check the -batch suffix assumption)")
print(f"barcode match rate: {rate:.1%} ({matched} cells)")

# ---- build integer pseudobulk matrix, filter, run DESeq2 ----
pb = pd.DataFrame(profiles, index=var_names).T          # donors x ALL genes
pb = np.round(pb).astype(int)                            # 0.5 multimapping counts -> int
pb = pb.loc[:, pb.sum(axis=0) >= cfg["min_counts_per_gene"]]
obs = pd.DataFrame({ck: [conds[d] for d in pb.index]}, index=pb.index)
pdata = ad.AnnData(X=pb.values.astype(float), obs=obs, var=pd.DataFrame(index=pb.columns))

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
pdata.obs[ck] = pdata.obs[ck].astype("category")
dds = DeseqDataSet(adata=pdata, design_factors=ck)
dds.deseq2()
stat = DeseqStats(dds, contrast=[ck, cfg["group_a"], cfg["group_b"]])
stat.summary()
stat.results_df.sort_values("padj").to_csv(sm.output[0])
n_donor = pdata.n_obs
print(f"fibroblast-restricted DE: {n_donor} donors, {pdata.n_vars} genes, "
      f"{cfg['group_a']} vs {cfg['group_b']}")
