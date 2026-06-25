"""Per-donor pseudobulk + PyDESeq2 for condition contrast (avoids pseudoreplication).
NOTE: decoupler/PyDESeq2 APIs vary by version; adjust arg names if needed."""
import scanpy as sc
import decoupler as dc

sm = snakemake  # noqa: F821
cfg = sm.config["de"]
adata = sc.read_h5ad(sm.input[0])

# sum raw counts per donor -> pseudobulk profiles (the real replication unit)
pdata = dc.get_pseudobulk(
    adata, sample_col=cfg["sample_key"], groups_col=None, layer="counts",
    mode="sum", min_cells=cfg["min_cells_per_sample"],
    min_counts=cfg["min_counts_per_gene"],
)

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

pdata.obs[cfg["condition_key"]] = pdata.obs[cfg["condition_key"]].astype("category")
dds = DeseqDataSet(adata=pdata, design_factors=cfg["condition_key"])
dds.deseq2()
stat = DeseqStats(dds, contrast=[cfg["condition_key"], cfg["group_a"], cfg["group_b"]])
stat.summary()
stat.results_df.sort_values("padj").to_csv(sm.output[0])
print(f"pseudobulk DE: {pdata.n_obs} samples, contrast {cfg['group_a']} vs {cfg['group_b']}")
