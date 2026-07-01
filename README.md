# The stiffening lung and NSCLC drug tolerance

An end-to-end single-cell pipeline built as a Snakemake workflow.
I wanted to learn scRNA-seq workflow and test predictions from a dynamical-systems
NSCLC model with real data. The central thesis is that a stiff/fibrotic microenvironment increases 
drug tolerance by increasing FAK/AKT/ERK activity in tumor cells, sustaining proliferation 
(cyclin D–CDK4/6) and promoting EGFR-TKI tolerance. The project tests this across four datasets 
that share one pipeline.

## Background Model

The mitosis-apoptosis model (soon to be submitted!) takes **AKT, ERK, and FAK as three 
independent measured inputs** where kinase activities are set per (stiffness, drug) condition from 
phospho-Westerns, feeding into a cell-cycle/apoptosis ODE system. Two considerations for this analysis 
are:

1. **The kinase inputs are phospho-regulated, not transcriptional.** scRNA-seq / RNA-seq cannot test
   AKT/ERK/FAK activity. Those are measured by Westerns in the manuscript. What
   transcriptomics can test is the model's downstream species: the cell-cycle components the ODE
   rate laws drive when the inputs rise (CCND1, MYC, SKP2 up; CDKN1A down; the E2F program following).
2. **YAP/TAZ is not a node in this model.** It routes mechanotransduction through FAK/AKT/ERK → cell
   cycle. YAP/TAZ target genes are scored only as an external mechanobiology readout (and because the
   Cosgrove source paper centers on them), never as evidence for this model.

## The four data sets and analyses

- **v1 - fibrosis (Tsukui 2020, E-CURD-126 / GSE132771).** Recover the CTHRC1+ pathologic fibroblast
  population and quantify the matrix-stiffening program (collagens + LOX/LOXL2 crosslinkers) that
  expands in fibrosis. *Result: the stiffening niche exists and is cell-intrinsic (pathologic
  fraction separates fibrotic vs normal, p=0.018; pseudobulk DE separates compositional from
  cell-intrinsic effects).*
- **v2 - NSCLC TME (Kim 2020, GSE131907).** Call malignant cells by CNV, characterize the tumor
  microenvironment. *Result: the CAF bridge - tumor fibroblasts cell-intrinsically
  reactivate the CTHRC1+ stiffening program (vs normal lung fibroblasts, p≈0; fibroblast-restricted
  PyDESeq2 confirms cell-intrinsic activation). Tumors build the stiff niche. Malignant-cell
  mechanotransduction signatures are flat-to-lower, and malignant calling is depth-confounded
  (r=−0.79), a droplet artifact that caps the tumor-intrinsic claim.*
- **v3 - treatment trajectory (Maynard 2020, PRJNA591860, Smart-seq2).** EGFR/ALK lung adeno across
  treatment-naive → residual-disease (persister) → progression. Cleanest CNV calling in the project
  (uniform Smart-seq2 depth, no confound). *Result: genuine persister biology, the alveolar-
  regenerative (AT2) program activates at residual disease, cyclin D–CDK4/6 collapses (persister
  quiescence). The treatment axis confounds drug/selection/time and never varies stiffness, so
  it cannot test the stiffness prediction; and the cell-level mechano-signatures are pseudoreplication. 
  They do not survive aggregation to the patient (YAP RD-vs-TN p≈1.0 across N=5). Reframed as a
  treatment-trajectory result, not a stiffness test.*
- **v4 - controlled stiffness experiments (Cosgrove 2024, *Science*; GSE243763, bulk RNA-seq).** A549
  lung adenocarcinoma on soft (1 kPa) vs stiff (50 kPa) hydrogels, 3 replicates/condition.
  *Result: the model's downstream cell-cycle targets shift coordinately in the predicted direction
  on stiff substrates (CCND1/MYC/SKP2 up, CDKN1A down, E2F program following; collective sign test
  p=0.016, oriented Wilcoxon vs genome-wide background p=0.013). No single gene significant at 3v3,
  the result is a coordinated set-level shift. A separate, strong YAP/TAZ response (padj 0.0017)
  reproduces the source paper but lies outside the model. Caveat: A549 is KRAS-mutant, the mechanism
  in lung adenocarcinoma, not the EGFR context.*

**Project arc:** tumors build the stiff niche (v1/v2, patient data) → cancer cells respond to it with
the cell-cycle program the model predicts (v4, controlled stiffness). The treatment-trajectory data
(v3) is presented as a different axis → the kinase-input layer is tested by the manuscript's
own phospho-Westerns, not by transcriptomics.

## Layout

```
config/        config.yaml (+ per-dataset config_<name>.yaml) and samples_<name>.tsv
workflow/
  Snakefile    targets + rule includes
  rules/       common, qc, integrate, cluster, biology
  scripts/     ingest_*.py, build_gene_positions.py, infer_cnv.py, caf_de.py,
               ingest_cosgrove.py, check_model_targets.py
resources/     downloaded matrices, per-dataset subdirs (gitignored)
results_<name>/  per-dataset outputs (results_kim, results_maynard, results_cosgrove, ...) (gitignored)
notebooks/     writeups: v1_writeup, v2_kim_writeup, v3_maynard_writeup, cosgrove_a549_writeup
               + stiffness_model_synthesis.md (cross-dataset synthesis)
```

Each dataset is isolated by its own `config_<name>.yaml` (own `resources_dir` / `results_dir` /
`samples`), so installments won't conflict. Run with `--configfile config/config_<name>.yaml`.

## Environment notes (Windows 11 + Anaconda. Read before running!)

The pipeline was developed on Windows; a few platform specific issues exist and are worth knowing
before starting, you should probably just run this on Linux :)

- **`--use-conda` is unreliable here.** The mamba activation banner pollutes Snakemake's Python
  version-probe JSON. Workaround: `pip install snakemake` into the hashed scanpy conda env, activate
  that env, and run Snakemake without `--use-conda` (all rules through CNV/DE share that one env,
  which also has scanpy / infercnvpy / pydeseq2 / pyensembl).
- **inferCNVpy on Windows** spawns subprocesses (tqdm `process_map`) that re-import the Snakemake
  script and crash. `infer_cnv.py` patches `process_map` with a serial map; `n_jobs=1` alone
  does not fix it. CNV scoring needs the full chain `cnv.tl.pca → cnv.pp.neighbors →
  cnv.tl.leiden → cnv.tl.cnv_score`.
- **HVG-subset** `merge_normalize` writes `annotated.h5ad` with only the 2,000 HVGs. CNV /
  gene-positions / pseudobulk-DE must read full-gene counts from the per-sample `resources/*.h5ad`,
  never from `annotated.h5ad`.
- **annotation fix.** An empty `annotation: {}` cannot override a populated map through
  Snakemake's deep-merge. Provide a complete cluster→label map under `cluster:` (every cluster id
  present). YAML: spaces not tabs, no semicolons.
- **Barcode joins across `ad.concat(index_unique="-")`:** strip the `-{batch}` suffix; the match rate
  ≈ the QC survival rate (≈99% for Kim, ≈80% for Maynard since its QC drops more failed wells. The
  CNV assert threshold is set accordingly).

## Quick start

```bash
# driver env (runs Snakemake itself)
mamba env create -f environment.yml && conda activate scrna-driver

# OR the working Windows path: activate the hashed scanpy env that has snakemake pip-installed,
# and run without --use-conda (see Environment notes)

# see the plan without running anything
snakemake -n --configfile config/config_kim.yaml

# run a target for one dataset
snakemake --configfile config/config_maynard.yaml --cores 4 results_maynard/cluster/annotated.h5ad
```

## Per-dataset entry points

- **v1 Tsukui:** `ingest_scea.py` converts the SCEA MatrixMarket bundle (one aggregated matrix +
  sidecars + experiment-design TSV) into per-donor `.h5ad`. SCEA counts are fractional → HVG flavor
  forced to `seurat`. (Details in the SCEA section below.)
- **v2 Kim:** `ingest_kim.py` (chunked loader for the big UMI matrix; integer UMIs → `seurat_v3`).
  Then `build_gene_positions.py` (pyensembl release 100) and `infer_cnv.py`. Fibroblast-restricted DE
  via `caf_de.py` (PyDESeq2; note author label is `Fibroblasts`, plural).
- **v3 Maynard:** `ingest_maynard.py` (chunked dense-CSV loader; strips ERCC spike-ins; maps the
  `analysis` column naive/grouped_pr/grouped_pd → TN/RD/PD; cell id is the `cell_id` column, not the
  row index). Smart-seq2 → `doublet_method: none`, `hvg_flavor: seurat`; mito genes are stripped from
  the matrix so the mito QC gate is disabled. Then the standard CNV path.
- **v4 Cosgrove:** `ingest_cosgrove.py` assembles the 9 A549 RSEM `*.genes.results` files
  (from `GSE243763_RAW.tar`), strips Ensembl version suffixes, maps ENSG→symbol (pyensembl, HFF-table
  fallback), runs PyDESeq2 stiff-vs-soft, and scores the signatures. Then `check_model_targets.py`
  runs the per-gene + collective tests on the model's downstream cell-cycle targets. Bulk + 3v3, so
  this is the simplest installment, no CNV, no integration.

## Getting the data

- **v1 Tsukui:** SCEA download tab (E-CURD-126) or GEO GSE132771 → `resources/tsukui/`.
- **v2 Kim:** GSE131907 processed matrix (`..._raw_UMI_matrix.txt.gz` + `cell_annotation.txt.gz`) →
  `resources/kim/`.
- **v3 Maynard:** processed CSVs from the study's Google Drive (`Data_input/csv_files/S01_datafinal.csv`
  + `S01_metacells.csv`), linked from `github.com/czbiohub-sf/scell_lung_adenocarcinoma`. Pull with
  `gdown --folder`. (Raw SRA under PRJNA591860 is too heavy; use the processed CSVs.)
- **v4 Cosgrove:** GEO GSE243763. `GSE243763_RAW.tar` holds the per-sample RSEM files; the 9 A549
  samples are GSM9224457–9224465. `GSE243763_SupplementaryTable2.csv.gz` is the HFF DE table
  (useful only as an ENSG→symbol fallback). Pull via `GEOparse` or direct HTTPS from the GEO FTP path.

### Ingesting the SCEA MatrixMarket bundle (v1)

SCEA gives you one aggregated matrix (genes × cells) + two sidecar files + a separate
experiment-design TSV. `scripts/ingest_scea.py` converts that bundle into per-donor `.h5ad`:

```bash
python workflow/scripts/ingest_scea.py \
    --mtx    resources/tsukui/E-CURD-126.aggregated_filtered_counts.mtx \
    --rows   resources/tsukui/E-CURD-126.aggregated_filtered_counts.mtx_rows \
    --cols   resources/tsukui/E-CURD-126.aggregated_filtered_counts.mtx_cols \
    --design resources/tsukui/ExpDesign-E-CURD-126.tsv \
    --outdir resources/tsukui
```

It transposes to cells × genes, resolves gene symbols, joins the design TSV for condition + donor,
writes one `.h5ad` per donor, and emits `config/samples_tsukui.tsv`. It prints the detected columns
and disease→condition mapping — check those, override with `--condition-col` / `--donor-col` /
`--id-col` if needed.

**Gene symbols.** SCEA's `.mtx_rows` is often Ensembl-only, which breaks `MT-` mito QC and
symbol-based signatures. Resolve via `--gtf path/to/annotation.gtf` (offline, reproducible,
recommended; use the release the data was built on) or `--use-mygene` (network, current annotations).
Unmapped genes keep their Ensembl ID; originals are preserved in `adata.var['gene_ids']`.

## How to work with it

Each `scripts/*.py` is a Snakemake `script:` step (the `snakemake` object is injected) or a
standalone CLI (the v4 scripts). They are starting points, run rule-by-rule and refined as you
inspect outputs, matching an iterative, output-driven workflow. The notebooks regenerate each
installment's figures from the `results_<name>/` outputs.

## Statistics conventions

- **DE and composition tests run on per-donor pseudobulk.** The donor is the replication unit, never
  the cell. Cell-level p-values are pseudoreplication.
- **Always pair cell-level signature claims with a per-donor check.** This is not optional: in v3, the
  dramatic cell-level signature shifts (KW p < 1e-100) almost entirely vanished at the patient level
  (N=5). The per-donor view is what the result actually is.
- **For multi-gene model predictions, test the set, not each gene.** With small replicate counts no
  single gene reaches significance; the right test is whether the predicted gene set is coordinately
  shifted (sign test + oriented-LFC Wilcoxon vs background, as in `check_model_targets.py`).
- **Match the axis of variation to the hypothesis.** A treatment axis (v3) cannot test a stiffness
  prediction however well-powered; only a stiffness contrast (v4) can.
- Signature scores export with donor/condition keys for aggregation and bootstrapping. Treat any
  subpopulation difference in a small donor subset as provisional until checked at the donor level.

## Open next steps

- **EGFR-context repeat of v4:** a transcriptome-wide soft-vs-stiff dataset in PC9 or HCC827
  (EGFR exon19del) would put the model's own population on the footing A549 (KRAS) currently provides.
  Martino et al. 2024 show the right direction (YAP/CTGF/CYR61/Ki67 up with stiffness) but only as
  targeted qPCR/protein.
- **The bimodality prediction (single-cell, still untested):** the model predicts a bimodal G1
  distribution (cycling + quiescent attractors) on soft + EGFRi substrates. The Cosgrove scRNA-seq
  CRISPRi series (GSE243756) could test whether the soft/stiff contrast produces that two-population
  cell-state structure, using non-targeting guides for the clean contrast, targeted guides as a
  perturbation test. More parsing work (guide demultiplexing), but it's the most distinctive remaining
  claim.
- **Spatial NSCLC:** define stiff/soft regions by local CAF/ECM density and compare the same tumor's
  cancer cells across them, the in-situ version of the v4 contrast.

See `notebooks/stiffness_model_synthesis.md` for the full cross-dataset synthesis with all numbers
and caveats.
