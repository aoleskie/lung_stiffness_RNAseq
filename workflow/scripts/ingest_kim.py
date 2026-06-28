#!/usr/bin/env python3
"""
Ingest Kim et al. 2020 (GSE131907) NSCLC atlas into per-sample .h5ad files.

The raw UMI matrix is a DENSE genes x cells text table (~180k cells x 29,634 genes),
so a full dense load is ~21 GB. This streams it in gene-row chunks, keeps only the
cells in the selected tissues (column subsetting during read), and builds a sparse
matrix — then attaches the author annotation (sample / tissue / cell type) and writes
one .h5ad per sample, plus a samples_kim.tsv for the Snakemake pipeline.

Usage:
    python ingest_kim.py \
        --matrix /path/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz \
        --anno   /path/GSE131907_Lung_Cancer_cell_annotation.txt.gz \
        --outdir resources/kim \
        --tissues tLung,nLung          # default; add mLN,mBrain,nLN,PE to include mets
        --chunksize 1000               # gene-rows per chunk; lower if memory is tight
"""
import argparse, gzip, os, sys
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, vstack
import anndata as ad

# tissue -> coarse condition (edit to taste; drives the pipeline's contrasts)
TISSUE_CONDITION = {
    "tLung": "tumor", "nLung": "normal",
    "mLN": "metastasis", "mBrain": "metastasis",
    "nLN": "normal_LN", "PE": "effusion",
}

def detect_columns(anno: pd.DataFrame, matrix_cells: list) -> dict:
    """Find barcode / sample / tissue / cell-type columns without hardcoding names."""
    cells = set(matrix_cells)
    # barcode col = the column whose values best overlap the matrix header
    overlaps = {c: anno[c].astype(str).isin(cells).mean() for c in anno.columns}
    barcode_col = max(overlaps, key=overlaps.get)
    if overlaps[barcode_col] < 0.5:
        sys.exit(f"ERROR: no annotation column matches matrix barcodes "
                 f"(best={barcode_col} at {overlaps[barcode_col]:.0%}). "
                 f"Columns: {list(anno.columns)}")
    low = {c: c.lower() for c in anno.columns}
    def find(*pats, exclude=()):
        for c in anno.columns:
            n = low[c]
            if any(p in n for p in pats) and not any(e in n for e in exclude):
                return c
        return None
    sample_col   = find("sample", exclude=("origin", "type")) or find("patient", "donor")
    tissue_col   = find("origin", "tissue", "site")
    celltype_col = find("cell_type", "celltype", "cell.type") or find("type", exclude=("sub",))
    return {"barcode": barcode_col, "sample": sample_col,
            "tissue": tissue_col, "celltype": celltype_col}

def load_matrix_subset(matrix_path, keep_barcodes, chunksize):
    """Stream the dense genes x cells matrix, keep only `keep_barcodes`, return (X cells x genes, genes, cells)."""
    op = gzip.open if matrix_path.endswith(".gz") else open
    with op(matrix_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    all_cells = header[1:]
    keep_pos = [i for i, c in enumerate(all_cells) if c in keep_barcodes]
    if not keep_pos:
        sys.exit("ERROR: none of the kept barcodes were found in the matrix header.")
    usecols = [0] + [i + 1 for i in keep_pos]
    kept_cells = [all_cells[i] for i in keep_pos]
    print(f"  matrix: {len(all_cells):,} cells total -> reading {len(kept_cells):,} after tissue filter")

    blocks, genes, seen = [], [], 0
    reader = pd.read_csv(matrix_path, sep="\t", header=0, usecols=usecols,
                         dtype={0: str}, chunksize=chunksize, compression="infer")
    for chunk in reader:
        genes.extend(chunk.iloc[:, 0].tolist())
        blocks.append(csr_matrix(chunk.iloc[:, 1:].to_numpy(dtype=np.float32)))
        seen += len(chunk)
        print(f"\r  genes read: {seen:,}", end="", flush=True)
    print()
    X = vstack(blocks).T.tocsr()        # genes x cells -> cells x genes
    # reorder columns of usecols-read frame back to kept_cells order (pandas preserves file order = kept_cells order)
    return X, genes, kept_cells

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--anno", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tissues", default="tLung,nLung")
    ap.add_argument("--chunksize", type=int, default=1000)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    tissues = [t.strip() for t in args.tissues.split(",") if t.strip()]

    print("Reading annotation ...")
    anno = pd.read_csv(args.anno, sep="\t", compression="infer", low_memory=False)
    print(f"  {len(anno):,} annotated cells, columns: {list(anno.columns)}")

    # peek matrix header to resolve the barcode column
    op = gzip.open if args.matrix.endswith(".gz") else open
    with op(args.matrix, "rt") as fh:
        matrix_cells = fh.readline().rstrip("\n").split("\t")[1:]
    cols = detect_columns(anno, matrix_cells)
    print(f"  detected columns -> {cols}")
    if cols["tissue"] is None:
        sys.exit("ERROR: could not find a tissue/origin column; pass through manually.")

    anno = anno.rename(columns={cols["barcode"]: "barcode", cols["sample"]: "sample",
                                cols["tissue"]: "tissue", cols["celltype"]: "cell_type"})
    print("\n  cells per tissue:\n" + anno["tissue"].value_counts().to_string())
    anno = anno[anno["tissue"].isin(tissues)].copy()
    keep = set(anno["barcode"].astype(str))
    print(f"\n  keeping {len(keep):,} cells from tissues {tissues}")

    print("\nLoading matrix (chunked, sparse) ...")
    X, genes, kept_cells = load_matrix_subset(args.matrix, keep, args.chunksize)
    anno = anno.set_index("barcode").loc[kept_cells]          # align to matrix order
    adata = ad.AnnData(X=X, obs=anno.reset_index(drop=False).set_index("barcode"),
                       var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    adata.var_names_make_unique()
    print(f"  assembled AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} genes "
          f"(sparsity {1 - adata.X.nnz/(adata.n_obs*adata.n_vars):.1%})")

    # write one file per sample + a samples sheet matching the pipeline schema
    # (common.smk indexes on sample_id and filters input_type=="matrix"; the qc rule
    #  reads donor_id + condition from the sheet, and sample_path uses the path column)
    rows = []
    for s, idx in adata.obs.groupby("sample", observed=True).groups.items():
        sub = adata[idx].copy()
        tissue = sub.obs["tissue"].iloc[0]
        cond = TISSUE_CONDITION.get(tissue, tissue)
        sub.obs["condition"] = cond
        path = os.path.join(args.outdir, f"{s}.h5ad").replace("\\", "/")
        sub.write(path)
        rows.append({
            "sample_id": s,            # index column common.smk expects
            "input_type": "matrix",    # so the sample lands in MATRIX_SAMPLES
            "path": path,              # resolved by sample_path()
            "donor_id": s,             # one sample = one patient-tissue unit
            "condition": cond,         # tumor / normal (contrast key)
            "tissue": tissue,          # kept for reference / optional covariate
            "n_cells": sub.n_obs,
        })
        print(f"  wrote {s}.h5ad  ({sub.n_obs:,} cells, {tissue} -> {cond})")
    samples = pd.DataFrame(rows).sort_values(["condition", "sample_id"])
    # write next to the pipeline's other sheets (config/) and alongside the data
    out_sheet = os.path.join("config", "samples_kim.tsv")
    os.makedirs("config", exist_ok=True)
    samples.to_csv(out_sheet, sep="\t", index=False)
    samples.to_csv(os.path.join(args.outdir, "samples_kim.tsv"), sep="\t", index=False)
    print(f"\nDone. {len(samples)} samples. Sheet -> {out_sheet}\n")
    print(samples.to_string(index=False))

if __name__ == "__main__":
    main()