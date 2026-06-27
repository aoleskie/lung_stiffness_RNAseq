import scanpy as sc, numpy as np, pandas as pd, anndata as ad
from scipy.sparse import issparse

sm = snakemake  # noqa: F821
cfg = sm.config["de"]
sk, ck = cfg["sample_key"], cfg["condition_key"]

profiles, conds, var_names = {}, {}, None
for f in sm.input:                                  # one QC file = one donor, all genes
    a = sc.read_h5ad(f)
    donor = str(a.obs[sk].iloc[0])
    c = a.layers["counts"] if "counts" in a.layers else a.X
    c = c.toarray() if issparse(c) else np.asarray(c)
    profiles[donor] = np.asarray(c).sum(axis=0).ravel()
    conds[donor] = str(a.obs[ck].iloc[0])
    var_names = a.var_names if var_names is None else var_names

pb = pd.DataFrame(profiles, index=var_names).T       # donors x ALL genes
pb = np.round(pb).astype(int)                         # 0.5 multimapping counts -> int
pb = pb.loc[:, pb.sum(axis=0) >= cfg["min_counts_per_gene"]]
obs = pd.DataFrame({ck: [conds[d] for d in pb.index]}, index=pb.index)
pdata = ad.AnnData(X=pb.values.astype(float), obs=obs,
                   var=pd.DataFrame(index=pb.columns))

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
pdata.obs[ck] = pdata.obs[ck].astype("category")
dds = DeseqDataSet(adata=pdata, design_factors=ck)
dds.deseq2()
stat = DeseqStats(dds, contrast=[ck, cfg["group_a"], cfg["group_b"]])
stat.summary()
stat.results_df.sort_values("padj").to_csv(sm.output[0])
print(f"pseudobulk DE: {pdata.n_obs} donors, {pdata.n_vars} genes tested, "
      f"{cfg['group_a']} vs {cfg['group_b']}")