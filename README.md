# scRNA-seq: the stiffening lung and NSCLC drug tolerance

An end-to-end single-cell pipeline (raw reads -> biology) built as a Snakemake workflow.
It teaches the canonical scRNA-seq workflow **and** tests predictions from a dynamical-systems
NSCLC model in real patient data, in two installments that share one pipeline:

- **v1 (fibrosis):** Tsukui 2020 collagen-producing cell atlas (E-CURD-126, GEO GSE132771).
  Recover the CTHRC1+ pathologic fibroblast population and quantify the matrix-stiffening
  program (collagens + LOX/LOXL2 crosslinkers) that expands in fibrosis.
- **v2 (NSCLC):** Kim 2020 LUAD TME (GSE131907; Maynard 2020 for TKI timepoints, if accessible).
  Call malignant cells by CNV, then ask whether drug-tolerant-looking tumor cells carry the
  mechanosensing (YAP/TAZ) and PI3K-AKT signatures the model predicts.

> scRNA-seq sees mRNA, not phospho-FAK or limit cycles, so the model's nodes are tested via their
> **transcriptional shadows** (PROGENy MAPK/PI3K activity, YAP/TAZ target genes, slow-cycling state).

## Layout

```
config/        config.yaml + samples.tsv  (edit these)
workflow/
  Snakefile    targets + rule includes
  rules/       refs, quant, qc, integrate, cluster, biology
  envs/        one conda env per tool group
  scripts/     scanpy / decoupler / PyDESeq2 steps
resources/     genomes, whitelists, downloaded matrices (gitignored)
results/       all outputs (gitignored)
notebooks/     interactive exploration (clean + executed)
```

## Quick start

```bash
# 1. driver env (runs Snakemake itself)
mamba env create -f environment.yml && conda activate scrna-driver

# 2. tell the pipeline about your data
#    edit config/samples.tsv  -> one row per sample (matrix path or fastq dir)
#    edit config/config.yaml  -> QC thresholds, contrast, signatures

# 3. see the plan without running anything (works before data is downloaded)
snakemake -n
snakemake --dag | dot -Tpng > dag.png      # picture of the DAG

# 4. run, letting Snakemake build each rule's conda env
snakemake --use-conda --cores 8

# build just one target
snakemake --use-conda --cores 8 results/cluster/annotated.h5ad
```

## Getting the data

- **v1 matrices:** Tsukui via the SCEA download tab (E-CURD-126) or GEO GSE132771. Put per-sample
  matrices under `resources/tsukui/` and point `samples.tsv` at them.
- **The FASTQ lesson (optional, once):** patient FASTQs are usually controlled-access, so use a
  single open 10x sample. Add it to `samples.tsv` as `input_type=fastq`, list its id under
  `fastq_lesson_samples`, and fill in the reference URLs in `config.yaml`.
- **v2 matrices:** Kim GSE131907 processed matrix from GEO (open). Maynard 2020 - verify access first.

### Ingesting the SCEA MatrixMarket bundle (v1)

SCEA hands you one aggregated matrix (genes x cells) + two sidecar files + a separate
experiment-design TSV - not per-sample matrices. `scripts/ingest_scea.py` converts that bundle
into the per-donor `.h5ad` files this pipeline expects. Download the **raw counts archive** and the
**experiment-design TSV** from E-CURD-126, unzip into `resources/tsukui/`, then run (in the scanpy env):

```bash
python workflow/scripts/ingest_scea.py \
    --mtx    resources/tsukui/E-CURD-126.aggregated_filtered_counts.mtx \
    --rows   resources/tsukui/E-CURD-126.aggregated_filtered_counts.mtx_rows \
    --cols   resources/tsukui/E-CURD-126.aggregated_filtered_counts.mtx_cols \
    --design resources/tsukui/ExpDesign-E-CURD-126.tsv \
    --outdir resources/tsukui
```

It transposes to cells x genes, resolves gene symbols, joins the design TSV for condition + donor,
writes one `.h5ad` per donor, and emits `config/samples_tsukui.tsv`. It prints exactly which columns
it detected and the disease->condition mapping - check those, and override with
`--condition-col` / `--donor-col` / `--id-col` if needed. Then copy the generated rows into
`config/samples.tsv` and run `snakemake -n`.

**Gene symbols.** SCEA's `.mtx_rows` is often Ensembl-only, which breaks `MT-` mito QC and
symbol-based signatures. Resolve symbols one of two ways:
- `--gtf path/to/annotation.gtf` - **offline and reproducible** (recommended). Use the GENCODE/Ensembl
  GTF for the release the data was built on; the script reads `gene_id`->`gene_name` and strips version suffixes.
- `--use-mygene` - queries mygene.info over the network. Quick, but maps against current annotations.

Genes with no symbol keep their Ensembl ID; the original IDs are always preserved in `adata.var['gene_ids']`.

## How to work with it

Each `scripts/*.py` is a Snakemake `script:` step (the `snakemake` object is injected). They are
**starting points** meant to be run rule-by-rule and refined as you inspect outputs - matching an
iterative, output-driven workflow. The notebooks are for interactive poking at the `.h5ad` objects.

## Statistics conventions (baked in)

- DE and composition tests run on **per-donor pseudobulk** (`scripts/pseudobulk_de.py`), never per
  cell - the donor is the replication unit.
- Signature scores export with donor/condition keys so you can aggregate to donor level and bootstrap.
- Treat any subpopulation difference in a small donor subset as provisional until bootstrapped.

## Switching v1 -> v2

1. Repoint `samples.tsv` at the NSCLC matrices; set `de.group_a/group_b` to the relevant contrast.
2. Annotate all lineages (consider CellTypist / HLCA label transfer) instead of fibroblast-only.
3. Enable the CNV rule in `rules/biology.smk` (scaffolded) to separate malignant from normal epithelium.
4. The v2 signature sets (yap_taz_targets, pi3k_akt_axis, cdk46_cellcycle) are already in `config.yaml`.
