#!/usr/bin/env python3
"""
Ingest an EBI Single Cell Expression Atlas (SCEA) MatrixMarket bundle into
per-donor .h5ad files that the Snakemake pipeline can consume.

SCEA gives you ONE aggregated matrix (genes x cells) plus two sidecar files and
a separate experiment-design TSV. This script:
  1. loads the .mtx + .mtx_rows (genes) + .mtx_cols (cells),
  2. transposes to cells x genes (AnnData orientation),
  3. resolves gene symbols (from the rows file, or via mygene as a fallback),
  4. joins the experiment-design TSV for condition (normal/fibrotic) + donor,
  5. writes one .h5ad per donor into --outdir and prints samples.tsv rows.

Run it once, in the scanpy conda env, before the pipeline:

  python workflow/scripts/ingest_scea.py \
      --mtx    resources/tsukui/E-CURD-126.aggregated_filtered_counts.mtx \
      --rows   resources/tsukui/E-CURD-126.aggregated_filtered_counts.mtx_rows \
      --cols   resources/tsukui/E-CURD-126.aggregated_filtered_counts.mtx_cols \
      --design resources/tsukui/ExpDesign-E-CURD-126.tsv \
      --outdir resources/tsukui

If the column auto-detection picks the wrong field, override with
--condition-col / --donor-col / --id-col and re-run. The script prints exactly
what it detected so you can check.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad


def log(msg):
    print(f"[ingest] {msg}", flush=True)


def die(msg):
    print(f"[ingest] ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def detect_col(df, keywords):
    """First column whose name contains any keyword (case-insensitive)."""
    for col in df.columns:
        low = col.lower()
        if any(k in low for k in keywords):
            return col
    return None


def sanitize(name):
    return "".join(c if c.isalnum() else "_" for c in str(name)).strip("_")


def normalize_condition(value):
    """Map a free-text disease value to the pipeline's 'normal' / 'fibrotic'."""
    v = str(value).lower()
    normal_keys = ["normal", "control", "healthy", "uninvolved", "non-fibrotic",
                   "non fibrotic", "unaffected"]
    return "normal" if any(k in v for k in normal_keys) else "fibrotic"


def load_genes(rows_path):
    """Return (ensembl_ids, symbols_or_None) from a .mtx_rows file."""
    rows = pd.read_csv(rows_path, sep="\t", header=None, dtype=str).fillna("")
    ensembl = rows.iloc[:, 0].astype(str).values
    symbols = None
    if rows.shape[1] >= 2:
        cand = rows.iloc[:, 1].astype(str).values
        # SCEA sometimes duplicates the ID into col 2; only use real symbols
        if not np.array_equal(cand, ensembl) and any(c for c in cand):
            symbols = np.array([
                c if (c and c.lower() != "nan" and c != e) else e
                for c, e in zip(cand, ensembl)
            ])
            log(f"using gene symbols from column 2 of {Path(rows_path).name}")
    return ensembl, symbols


def map_symbols_mygene(ensembl):
    """Best-effort Ensembl->symbol via mygene (needs network). None on failure."""
    try:
        import mygene
        mg = mygene.MyGeneInfo()
        res = mg.querymany(list(ensembl), scopes="ensembl.gene", fields="symbol",
                           species="human", as_dataframe=True, verbose=False)
        res = res[~res.index.duplicated(keep="first")]
        sym = res.reindex(ensembl)["symbol"]
        out = np.array([s if isinstance(s, str) and s else e
                        for s, e in zip(sym.values, ensembl)])
        log(f"mygene resolved {int((out != ensembl).sum())}/{len(ensembl)} symbols")
        return out
    except Exception as e:  # noqa: BLE001
        log(f"mygene mapping unavailable ({e}); keeping Ensembl IDs as var_names")
        return None


def main():
    ap = argparse.ArgumentParser(description="Ingest an SCEA MatrixMarket bundle.")
    ap.add_argument("--mtx", required=True)
    ap.add_argument("--rows", required=True)
    ap.add_argument("--cols", required=True)
    ap.add_argument("--design", required=True, help="experiment-design TSV")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--id-col", default=None, help="design column matching cell IDs (default: first column)")
    ap.add_argument("--condition-col", default=None, help="design column for disease/condition")
    ap.add_argument("--donor-col", default=None, help="design column for donor/individual")
    ap.add_argument("--samples-out", default="config/samples_tsukui.tsv")
    ap.add_argument("--use-mygene", action="store_true",
                    help="try mygene if the rows file has no symbols")
    args = ap.parse_args()

    for p in [args.mtx, args.rows, args.cols, args.design]:
        if not Path(p).exists():
            die(f"file not found: {p}")

    # ---- matrix + sidecars -------------------------------------------------
    log(f"reading matrix {Path(args.mtx).name} ...")
    M = scipy.io.mmread(args.mtx)          # SCEA orientation: genes x cells
    M = sp.csr_matrix(M)
    ensembl, symbols = load_genes(args.rows)
    cells = pd.read_csv(args.cols, sep="\t", header=None, dtype=str).iloc[:, 0].astype(str).values
    log(f"matrix {M.shape}  | genes(rows file)={len(ensembl)}  cells(cols file)={len(cells)}")

    # orient to cells x genes
    if M.shape == (len(ensembl), len(cells)):
        X = M.T.tocsr()
    elif M.shape == (len(cells), len(ensembl)):
        log("matrix already cells x genes; not transposing")
        X = M
    else:
        die(f"matrix shape {M.shape} matches neither genes x cells "
            f"({len(ensembl)} x {len(cells)}) nor its transpose. "
            f"Check that rows/cols files belong to this matrix.")

    # gene symbols
    if symbols is None and args.use_mygene:
        symbols = map_symbols_mygene(ensembl)
    var_names = symbols if symbols is not None else ensembl
    if symbols is None:
        log("WARNING: no gene symbols resolved; var_names are Ensembl IDs. "
            "Mitochondrial QC (MT- prefix) and symbol-based signatures will not "
            "match until symbols are added. Re-run with --use-mygene (needs net).")

    var = pd.DataFrame({"gene_ids": ensembl}, index=pd.Index(var_names, name=None))
    obs = pd.DataFrame(index=pd.Index(cells, name=None))
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    log(f"built AnnData: {adata.n_obs} cells x {adata.n_vars} genes")

    # ---- experiment design metadata ---------------------------------------
    design = pd.read_csv(args.design, sep="\t", header=0, dtype=str)
    id_col = args.id_col or design.columns[0]
    if id_col not in design.columns:
        die(f"id column '{id_col}' not in design; columns: {list(design.columns)}")
    cond_col = args.condition_col or detect_col(design, ["disease", "condition", "phenotype"])
    donor_col = args.donor_col or detect_col(design, ["individual", "donor", "patient", "subject"])
    log(f"design columns ({len(design.columns)}): {list(design.columns)}")
    log(f"using id_col='{id_col}'  condition_col='{cond_col}'  donor_col='{donor_col}'")
    if cond_col is None or donor_col is None:
        die("could not auto-detect condition/donor columns; pass --condition-col / --donor-col")

    design = design.drop_duplicates(subset=[id_col]).set_index(id_col)
    common = adata.obs_names.intersection(design.index)
    rate = len(common) / max(adata.n_obs, 1)
    log(f"cell-ID match to design: {len(common)}/{adata.n_obs} ({rate:.1%})")
    if len(common) == 0:
        log("first 3 matrix cell IDs:  " + ", ".join(map(str, adata.obs_names[:3])))
        log("first 3 design cell IDs:  " + ", ".join(map(str, design.index[:3])))
        die("no cell IDs matched the design file. Wrong --id-col, or mismatched bundle.")
    if rate < 1.0:
        log(f"dropping {adata.n_obs - len(common)} cells with no metadata")
    adata = adata[common].copy()

    adata.obs["disease_raw"] = design.loc[adata.obs_names, cond_col].astype(str).values
    adata.obs["donor_id"] = design.loc[adata.obs_names, donor_col].astype(str).values
    adata.obs["condition"] = [normalize_condition(v) for v in adata.obs["disease_raw"]]

    # show the disease -> condition mapping so it can be verified
    mapping = (adata.obs[["disease_raw", "condition"]]
               .drop_duplicates().sort_values("condition"))
    log("disease_raw -> condition mapping:")
    for _, r in mapping.iterrows():
        log(f"    {r['disease_raw']!r:40s} -> {r['condition']}")

    # ---- per-donor split + samples.tsv ------------------------------------
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows_out = []
    for donor, sub in adata.obs.groupby("donor_id", observed=True):
        conds = sub["condition"].unique()
        if len(conds) > 1:
            log(f"WARNING: donor {donor} has mixed conditions {list(conds)}; using majority")
            cond = sub["condition"].value_counts().idxmax()
        else:
            cond = conds[0]
        sid = sanitize(donor)
        sub_adata = adata[adata.obs_names[adata.obs["donor_id"] == donor]].copy()
        path = outdir / f"{sid}.h5ad"
        sub_adata.write(path)
        rows_out.append((sid, sid, cond, "matrix", str(path).replace("\\", "/")))
        log(f"  wrote {path.name}: {sub_adata.n_obs} cells  [{cond}]")

    samples_df = pd.DataFrame(rows_out,
                              columns=["sample_id", "donor_id", "condition", "input_type", "path"])
    Path(args.samples_out).parent.mkdir(parents=True, exist_ok=True)
    samples_df.to_csv(args.samples_out, sep="\t", index=False)
    log(f"wrote {args.samples_out} ({len(samples_df)} donors). Copy its rows into config/samples.tsv:")
    print()
    print(samples_df.to_string(index=False))
    print()
    n_norm = (samples_df["condition"] == "normal").sum()
    n_fib = (samples_df["condition"] == "fibrotic").sum()
    log(f"done: {len(samples_df)} donors  ({n_norm} normal, {n_fib} fibrotic)")


if __name__ == "__main__":
    main()
