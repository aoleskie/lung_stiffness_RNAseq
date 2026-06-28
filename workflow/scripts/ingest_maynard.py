#!/usr/bin/env python3
"""
Ingest Maynard 2020 (BioProject PRJNA591860) Smart-seq2 lung adenocarcinoma into the pipeline.

Smart-seq2, NOT 10x: plate-based full-length, read counts (not UMIs), ERCC spike-ins present,
one cell per well (no droplet doublets). Treatment axis lives in metadata column `analysis`:
  naive (TN) / grouped_pr (persister, residual disease) / grouped_pd (progressive disease).

Inputs (from the study's Google Drive Data_input/csv_files/):
  S01_datafinal.csv  : counts, genes x cells (dense; ~27k x 27k -> chunked load)
  S01_metacells.csv  : per-cell metadata (41 cols; analysis / sample_name / patient_id / driver_gene ...)

Writes per-sample h5ad files + a samples sheet matching common.smk's schema
(sample_id / input_type=matrix / path / donor_id / condition / plus stage & driver for reference).

Usage:
  python ingest_maynard.py --indir Data_input/csv_files --outdir resources/maynard --chunksize 2000
"""
import os
import argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

# analysis -> clean stage label. Keep BOTH a fine stage and a coarse tumor/normal-style `condition`.
STAGE_MAP = {"naive": "TN", "grouped_pr": "RD", "grouped_pd": "PD"}   # RD = residual disease (persister)


def load_counts_chunked(path, keep_cells, chunksize):
    """Stream a dense genes x cells CSV into a sparse cells x genes matrix, keeping only `keep_cells`
    columns and dropping ERCC spike-in rows. Returns (csr cells x genes, gene_names, cell_names)."""
    # header first: which columns (cells) to keep, in file order
    header = pd.read_csv(path, index_col=0, nrows=0)
    all_cells = header.columns.tolist()
    keep_mask = np.array([c in keep_cells for c in all_cells])
    kept_cells = [c for c in all_cells if c in keep_cells]
    print(f"  matrix has {len(all_cells):,} cells; keeping {len(kept_cells):,} present in metadata")

    gene_names, blocks, n_ercc = [], [], 0
    reader = pd.read_csv(path, index_col=0, chunksize=chunksize)
    for i, chunk in enumerate(reader):
        # drop ERCC spike-in rows
        is_ercc = chunk.index.astype(str).str.startswith("ERCC")
        n_ercc += int(is_ercc.sum())
        chunk = chunk.loc[~is_ercc]
        gene_names.extend(chunk.index.tolist())
        # subset to kept cells, store as sparse (genes_chunk x cells_kept)
        blocks.append(sp.csr_matrix(chunk.loc[:, keep_mask].to_numpy(dtype=np.float32)))
        if (i + 1) % 5 == 0:
            print(f"    processed {sum(b.shape[0] for b in blocks):,} genes ...")
    X_genes_by_cells = sp.vstack(blocks, format="csr")
    print(f"  dropped {n_ercc} ERCC rows; final {X_genes_by_cells.shape[0]:,} genes")
    # transpose to cells x genes
    return X_genes_by_cells.T.tocsr(), gene_names, kept_cells


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="Data_input/csv_files")
    ap.add_argument("--outdir", default="resources/maynard")
    ap.add_argument("--chunksize", type=int, default=2000)
    ap.add_argument("--min-cells-per-sample", type=int, default=20)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    meta = pd.read_csv(os.path.join(args.indir, "S01_metacells.csv"), index_col=0)
    print(f"metadata: {meta.shape[0]:,} cells x {meta.shape[1]} cols")
    # cell identifier: the matrix columns match the metadata index (well/cell ids)
    meta = meta.set_index("cell_id")
    meta.index = meta.index.astype(str)
    meta["stage"] = meta["analysis"].map(STAGE_MAP)
    if meta["stage"].isna().any():
        print(f"  WARNING: {meta['stage'].isna().sum()} cells with unmapped `analysis` value(s): "
              f"{sorted(meta.loc[meta['stage'].isna(),'analysis'].unique())}")

    counts_path = os.path.join(args.indir, "S01_datafinal.csv")
    X, genes, cells = load_counts_chunked(counts_path, set(meta.index), args.chunksize)

    # align metadata to the kept cell order
    obs = meta.loc[cells].copy()
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    adata.var_names_make_unique()
    for col in adata.obs.columns:
        if adata.obs[col].dtype == object or str(adata.obs[col].dtype) == "category":
            adata.obs[col] = adata.obs[col].astype(str).replace({"nan": "", "None": ""})

    print(f"\nassembled: {adata.n_obs:,} cells x {adata.n_vars:,} genes "
          f"(sparsity {1 - adata.X.nnz/(adata.n_obs*adata.n_vars):.1%})")
    print("stage distribution:", dict(adata.obs['stage'].value_counts()))

    # write one file per sample (biopsy) + a samples sheet matching common.smk
    rows = []
    for s, idx in adata.obs.groupby("sample_name", observed=True).groups.items():
        sub = adata[idx].copy()
        if sub.n_obs < args.min_cells_per_sample:
            print(f"  skip {s}: only {sub.n_obs} cells (< {args.min_cells_per_sample})")
            continue
        sid = str(s).replace(" ", "_").replace("/", "_")
        stage = sub.obs["stage"].iloc[0]          # stage is a biopsy-level attribute
        driver = str(sub.obs["driver_gene"].iloc[0]) if "driver_gene" in sub.obs else "NA"
        path = os.path.join(args.outdir, f"{sid}.h5ad").replace("\\", "/")
        sub.write(path)
        rows.append({
            "sample_id": sid,
            "input_type": "matrix",
            "path": path,
            "donor_id": str(sub.obs["patient_id"].iloc[0]),  # patient = donor
            "condition": stage,                              # TN / RD / PD (the treatment axis)
            "stage": stage,
            "driver_gene": driver,
            "n_cells": sub.n_obs,
        })
        print(f"  wrote {sid}.h5ad  ({sub.n_obs} cells, stage={stage}, driver={driver})")

    samples = pd.DataFrame(rows).sort_values(["condition", "sample_id"])
    os.makedirs("config", exist_ok=True)
    samples.to_csv("config/samples_maynard.tsv", sep="\t", index=False)
    samples.to_csv(os.path.join(args.outdir, "samples_maynard.tsv"), sep="\t", index=False)
    print(f"\nDone. {len(samples)} samples. Sheet -> config/samples_maynard.tsv\n")
    print(samples[["sample_id", "donor_id", "condition", "driver_gene", "n_cells"]].to_string(index=False))
    print("\nstage x sample counts:")
    print(samples.groupby("condition")["n_cells"].agg(["count", "sum"]).to_string())


if __name__ == "__main__":
    main()
