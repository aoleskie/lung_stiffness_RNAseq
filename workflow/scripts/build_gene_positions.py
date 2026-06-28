#!/usr/bin/env python3
"""
Build a gene-position table for inferCNVpy, aligned to a dataset's var_names.

inferCNVpy orders genes along the genome to detect large-scale CNV, so every gene
needs chromosome + start + end. This pulls coordinates from pyensembl (release 100,
GRCh38 — Kim's era) for the genes actually present in the AnnData, writes a TSV, and
reports how many genes resolved (unresolved genes are simply dropped from CNV, which
is fine — CNV is a smoothed large-window signal).

Usage:
    python build_gene_positions.py \
        --h5ad results_kim/cluster/annotated.h5ad \
        --out  resources/kim/gene_positions.tsv \
        --release 100
"""
import argparse
import pandas as pd
import scanpy as sc
from pyensembl import EnsemblRelease

# autosomes + sex chromosomes only (drop scaffolds/MT — not used for CNV)
CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
CHROM_ORDER = {c: i for i, c in enumerate(CHROMS)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--release", type=int, default=100)
    args = ap.parse_args()

    data = EnsemblRelease(args.release)
    var_names = sc.read_h5ad(args.h5ad, backed="r").var_names.tolist()
    print(f"{len(var_names):,} genes in dataset; querying Ensembl release {args.release} ...")

    rows, miss, ambig = [], 0, 0
    for g in var_names:
        try:
            hits = data.genes_by_name(g)          # match by gene symbol
        except Exception:
            miss += 1
            continue
        # keep only hits on a standard chromosome; if several, take the first standard one
        hits = [h for h in hits if h.contig in CHROM_ORDER]
        if not hits:
            miss += 1
            continue
        if len(hits) > 1:
            ambig += 1
        h = hits[0]
        rows.append({"gene": g, "chromosome": h.contig,
                     "start": int(h.start), "end": int(h.end)})

    pos = pd.DataFrame(rows)
    # order along the genome: chromosome, then start
    pos["chr_order"] = pos["chromosome"].map(CHROM_ORDER)
    pos = pos.sort_values(["chr_order", "start"]).drop(columns="chr_order")
    pos.to_csv(args.out, sep="\t", index=False)

    print(f"resolved {len(pos):,}/{len(var_names):,} genes "
          f"({len(pos)/len(var_names):.1%}); {miss:,} unmapped, {ambig:,} multi-locus")
    print("per-chromosome gene counts:")
    print(pos["chromosome"].value_counts().reindex(CHROMS).to_string())
    print(f"\nwrote {args.out}")

if __name__ == "__main__":
    main()
