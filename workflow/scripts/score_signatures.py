"""Score curated signatures + PROGENy pathway activity; export per-cell scores."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

sm = snakemake  # noqa: F821
adata = sc.read_h5ad(sm.input[0])
var_names = adata.raw.var_names if adata.raw is not None else adata.var_names

score_cols = []
for name, genes in sm.config["signatures"].items():
    present = [g for g in genes if g in var_names]
    if len(present) >= 2:
        sc.tl.score_genes(adata, present, score_name=name, use_raw=adata.raw is not None)
        score_cols.append(name)
        print(f"{name}: scored on {len(present)}/{len(genes)} genes")
    else:
        print(f"{name}: SKIPPED (<2 genes present)")

# PROGENy pathway activity (MAPK/EGFR ~ ERK, PI3K ~ AKT, TGFb, ...). API can
# differ across decoupler versions; wrapped so a miss doesn't fail the rule.
if sm.config["pathways"].get("use_progeny"):
    try:
        import decoupler as dc
        prog = dc.get_progeny(organism=sm.config["pathways"]["progeny_organism"],
                              top=sm.config["pathways"]["progeny_top"])
        dc.run_mlm(adata, prog, source="source", target="target",
                   weight="weight", use_raw=True)
        acts = adata.obsm["mlm_estimate"]
        for c in acts.columns:
            col = f"progeny_{c}"
            adata.obs[col] = acts[c].values
            score_cols.append(col)
    except Exception as e:
        print("PROGENy skipped:", e)

keys = ["sample_id", "donor_id", "condition", "cell_type"]
adata.obs[keys + score_cols].to_csv(sm.output["scores"])

curated = [c for c in score_cols if not c.startswith("progeny_")][:6]
if curated:
    sc.pl.umap(adata, color=curated, ncols=3, show=False)
    plt.savefig(sm.output["fig"], dpi=150, bbox_inches="tight")
else:
    plt.figure(); plt.savefig(sm.output["fig"])
print("scored signatures:", score_cols)
