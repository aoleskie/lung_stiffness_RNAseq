# Transcriptomic evidence for stromal remodeling and stiffness-associated responses in NSCLC

This repository evaluates which components of a proposed stiffness-to-drug-tolerance mechanism are
supported by four public lung fibrosis and lung adenocarcinoma datasets. The proposed mechanism is
that a fibrotic, mechanically altered microenvironment increases FAK/AKT/ERK activity in tumor
cells, sustains cyclin D-CDK4/6-associated proliferation, and promotes EGFR-TKI tolerance.

The datasets address distinct components of this mechanism. Fibrosis and NSCLC single-cell data
characterize stromal ECM-remodeling programs; a controlled hydrogel experiment tests
stiffness-associated transcriptional responses in lung adenocarcinoma cells; and longitudinal
patient data characterize treatment-associated persister states. No single dataset tests the full
causal chain.

## The model being tested

The mitosis-apoptosis model takes **AKT, ERK, and FAK as three independent
measured inputs**: kinase activities set per (stiffness, drug) condition from phospho-Westerns,
feeding a cell-cycle/apoptosis ODE system. Two aspects of that structure define what the
transcriptomic analyses can test.

1. **The kinase inputs are phospho-regulated, not transcriptional.** No amount of scRNA-seq will
   measure AKT/ERK/FAK activity; the manuscript measures those by Western. What transcriptomics can
   test is the layer underneath, the cell-cycle species the ODE rate laws drive when the inputs
   rise (CCND1, MYC, SKP2 up; CDKN1A down; the E2F program following a step later).
2. **YAP/TAZ is not a node in this model.** Mechanotransduction routes through FAK/AKT/ERK into the
   cell cycle. YAP/TAZ targets are scored throughout as an external mechanobiology readout, partly
   because the Cosgrove paper centers on them, and never as evidence for this model.

## The four datasets

- **v1, fibrosis (Tsukui 2020, E-CURD-126 / GSE132771).** Recovers the CTHRC1+ pathologic fibroblast
  population and quantifies an ECM-remodeling program that includes collagens and LOX-family
  crosslinkers. *Result: CTHRC1+ fibroblasts expand in fibrotic lung, with complete donor-level
  separation from normal lung (5 fibrotic versus 3 normal donors; exact p = 0.018). The
  fibroblast-restricted pseudobulk analysis also identifies increased ECM-remodeling expression
  within the fibroblast compartment. Because fibroblast states remain heterogeneous, this result
  is not interpreted as strictly cell-intrinsic activation.*
- **v2, the NSCLC microenvironment (Kim 2020, GSE131907).** Call malignant cells by CNV and
  characterize the tumor microenvironment. *Result: tumor-associated fibroblasts recapitulate a
  CTHRC1-associated fibrotic ECM-remodeling program linked to matrix stiffening. The signal persists
  in fibroblast-restricted pseudobulk, distinguishing it from immune or epithelial composition;
  however, differences among fibroblast states remain a possible contributor. Malignant-cell
  mechanotransduction signatures are flat or lower, and CNV-based malignant calling is confounded
  by sequencing depth (r = -0.79), limiting tumor-cell inference.*
- **v3, the treatment trajectory (Maynard 2020, PRJNA591860, Smart-seq2).** EGFR/ALK lung
  adenocarcinoma sampled from treatment-naive through residual disease (the persister state) to
  progression. Smart-seq2 depth reduces the CNV-scoring confound observed in Kim. *Result: residual
  disease is associated with an increased alveolar-regenerative/AT2 program and decreased cyclin
  D-CDK4/6-associated expression, consistent with a treatment-associated persister state. The
  design combines drug exposure, elapsed time, and clonal selection and contains no stiffness axis;
  it therefore provides orthogonal treatment-state context rather than a direct stiffness test.
  Cell-level mechanotransduction differences do not persist after aggregation to patients (YAP/TAZ,
  RD versus TN, p approximately 1.0; N = 5 per group).*
- **v4, controlled stiffness (Cosgrove 2024, *Science*; GSE243763, bulk RNA-seq).** A549 lung
  adenocarcinoma grown on soft (1 kPa) and stiff (50 kPa) hydrogels, 3 replicates per condition.
  *Result: a prespecified score of the four direct downstream targets (CCND1/MYC/SKP2 up, CDKN1A
  down) is higher in all 3 stiff replicates than in all 3 soft replicates (stiff - soft = +0.97
  oriented z-score units, Welch 95% CI +0.03 to +1.91, exact one-sided permutation p = 0.050).
  E2F1/CCNE1 are retained as descriptive downstream support rather than pooled into the primary
  test. No single gene clears significance at 3 vs 3.
  A separate, pronounced YAP/TAZ response (mean LFC +0.58; median constituent gene-level adjusted
  p-value 0.0017) reproduces the source paper but sits outside the model. The 0.0017 value is a
  descriptive median, not a pathway-level adjusted p-value. Because A549 is KRAS-mutant, this
  result supports a stiffness-associated response in lung adenocarcinoma but does not establish the
  response in the EGFR-mutant context addressed by the model.*

**Evidence synthesis.** Fibrotic lung and NSCLC stroma share a CTHRC1-associated ECM-remodeling
program consistent with generation of a mechanically altered microenvironment (v1 and v2). In a
controlled stiffness experiment, lung adenocarcinoma cells show coordinated transcriptional changes
in prespecified cell-cycle targets downstream of the model's mechanosensitive inputs (v4). The
treatment trajectory provides separate evidence about persister-state biology but cannot connect
that state to stiffness (v3). Transcriptomics does not directly assay the phospho-regulated kinase
inputs or validate the complete stiffness-to-drug-tolerance causal chain.

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
```

Each dataset is walled off by its own `config_<name>.yaml`, with its own `resources_dir`,
`results_dir`, and `samples`, so the installments cannot collide. Run with
`--configfile config/config_<name>.yaml`.

## Platform notes (Windows 11 and Anaconda)

The following Windows-specific issues affect reproducibility and environment activation.

- **`--use-conda` is unreliable here.** The mamba activation banner pollutes the JSON that
  Snakemake's Python version-probe expects to parse. Workaround: `pip install snakemake` into the
  hashed scanpy conda env, activate that env, and run Snakemake without `--use-conda`. Every rule
  from QC through CNV and DE shares that one env, which already carries scanpy, infercnvpy,
  pydeseq2, and pyensembl.
- **inferCNVpy on Windows** spawns subprocesses (tqdm's `process_map`) that re-import the Snakemake
  script and crash. `infer_cnv.py` patches `process_map` with a serial map; setting `n_jobs=1` on
  its own does not fix it. CNV scoring needs the full chain `cnv.tl.pca → cnv.pp.neighbors →
  cnv.tl.leiden → cnv.tl.cnv_score`.
- **Full-gene analyses must not use the HVG subset.** `merge_normalize` writes `annotated.h5ad` with only the 2,000 HVGs.
  CNV, gene positions, and pseudobulk DE all need full-gene counts and must read the per-sample
  `resources/*.h5ad`, never `annotated.h5ad`.
- **Annotation overrides.** An empty `annotation: {}` cannot override a populated map through
  Snakemake's deep-merge, so provide a complete cluster→label map under `cluster:` with every
  cluster id present. YAML wants spaces, not tabs, and no semicolons.
- **Barcode joins across `ad.concat(index_unique="-")`.** Strip the `-{batch}` suffix; the match
  rate lands at roughly the QC survival rate (about 99% for Kim, about 80% for Maynard, whose QC
  drops more failed wells). The CNV assert threshold is set accordingly.

## Quick start

```bash
# driver env (runs Snakemake itself)
mamba env create -f environment.yml && conda activate scrna-driver
```

Or take the working Windows path: activate the hashed scanpy env that has snakemake pip-installed
and run without `--use-conda`, as described above.

```bash
# see the plan without running anything
snakemake -n --configfile config/config_kim.yaml
```

```bash
# run a target for one dataset
snakemake --configfile config/config_maynard.yaml --cores 4 results_maynard/cluster/annotated.h5ad
```

## Per-dataset entry points

- **v1 Tsukui:** `ingest_scea.py` converts the SCEA MatrixMarket bundle (one aggregated matrix plus
  sidecars plus an experiment-design TSV) into per-donor `.h5ad`. SCEA counts are fractional, so
  HVG flavor is forced to `seurat`. Details in the SCEA section below.
- **v2 Kim:** `ingest_kim.py` chunk-loads the large UMI matrix; the counts are integer UMIs, so
  `seurat_v3`. Then `build_gene_positions.py` (pyensembl release 100) and `infer_cnv.py`.
  Fibroblast-restricted DE runs through `caf_de.py`; the author-provided label is the plural
  `Fibroblasts`.
- **v3 Maynard:** `ingest_maynard.py` chunk-loads the dense CSV, strips ERCC spike-ins, and maps the
  `analysis` column naive/grouped_pr/grouped_pd onto TN/RD/PD. Cell identity lives in the `cell_id`
  column, not the row index. Smart-seq2 means `doublet_method: none` and `hvg_flavor: seurat`, and
  since mito genes are stripped from the matrix the mito QC gate is disabled. Then the standard CNV
  path.
- **v4 Cosgrove:** `ingest_cosgrove.py` assembles the 9 A549 RSEM `*.genes.results` files from
  `GSE243763_RAW.tar`, strips Ensembl version suffixes, maps ENSG→symbol (pyensembl, with an
  HFF-table fallback), runs PyDESeq2 stiff-vs-soft, and exports counts, TPM, and signature scores.
  Then `check_model_targets.py` preserves the per-gene LFC table descriptively and runs the
  prespecified oriented target score with biological replicate as the inferential unit. Bulk data
  at 3 versus 3 and requires no CNV or integration step.

## Getting the data

- **v1 Tsukui:** the SCEA download tab (E-CURD-126) or GEO GSE132771, into `resources/tsukui/`.
- **v2 Kim:** the GSE131907 processed matrix (`..._raw_UMI_matrix.txt.gz` plus
  `cell_annotation.txt.gz`), into `resources/kim/`.
- **v3 Maynard:** processed CSVs from the study's Google Drive (`Data_input/csv_files/S01_datafinal.csv`
  and `S01_metacells.csv`), linked from `github.com/czbiohub-sf/scell_lung_adenocarcinoma`. Pull them
  with `gdown --folder`. The processed CSVs are the recommended inputs; the raw SRA under
  PRJNA591860 is substantially larger and is not required for this workflow.
- **v4 Cosgrove:** GEO GSE243763. `GSE243763_RAW.tar` holds the per-sample RSEM files, and the 9
  A549 samples are GSM9224457–9224465. `GSE243763_SupplementaryTable2.csv.gz` is the HFF DE table,
  useful here only as an ENSG→symbol fallback. Pull via `GEOparse` or direct HTTPS from the GEO FTP
  path.

### Ingesting the SCEA MatrixMarket bundle (v1)

SCEA hands you one aggregated matrix (genes × cells), two sidecar files, and a separate
experiment-design TSV. `scripts/ingest_scea.py` turns that bundle into per-donor `.h5ad`:

```bash
python workflow/scripts/ingest_scea.py \
    --mtx    resources/tsukui/E-CURD-126.aggregated_filtered_counts.mtx \
    --rows   resources/tsukui/E-CURD-126.aggregated_filtered_counts.mtx_rows \
    --cols   resources/tsukui/E-CURD-126.aggregated_filtered_counts.mtx_cols \
    --design resources/tsukui/ExpDesign-E-CURD-126.tsv \
    --outdir resources/tsukui
```

It transposes to cells × genes, resolves gene symbols, joins the design TSV for condition and
donor, writes one `.h5ad` per donor, and emits `config/samples_tsukui.tsv`. It also prints the
columns it detected and the inferred disease-to-condition mapping. Review these diagnostics and
override with `--condition-col`, `--donor-col`, or `--id-col` when necessary.

**Gene symbols.** SCEA's `.mtx_rows` is often Ensembl-only, which quietly breaks both `MT-` mito QC
and every symbol-based signature. Resolve them with `--gtf path/to/annotation.gtf` (offline,
reproducible; use the release on which the data were built) or `--use-mygene`
(network, current annotations). Unmapped genes keep their Ensembl ID, and the originals are
preserved in `adata.var['gene_ids']` either way.

## How to work with it

Each `scripts/*.py` is either a Snakemake `script:` step, with the `snakemake` object injected, or a
standalone CLI (the v4 scripts). Treat them as starting points. I ran them rule by rule and refined
them while looking at the outputs, which is the workflow they are shaped for. The notebooks
regenerate each installment's figures from the `results_<name>/` outputs.

## Statistics conventions

The analyses use the following statistical conventions.

- **DE and composition tests run on per-donor pseudobulk.** The donor is the replication unit, never
  the cell. Cell-level tests across a small number of donors constitute pseudoreplication.
- **Pair every cell-level signature claim with a per-donor check.** This is not optional. In v3 the
  large cell-level signature differences (KW p < 1e-100) were not retained at the patient level
  across N = 5 donors. The per-donor result is the relevant inferential result; the cell-level test
  counted non-independent observations.
- **For a multi-gene prediction, score the set within each biological replicate.** Correlated genes
  are not independent trials. `check_model_targets.py` orients and averages the prespecified direct
  targets into one score per sample, then reports the stiff-minus-soft effect, uncertainty, and an
  exact label-permutation p-value across the 3-vs-3 design.
- **Match the axis of variation to the hypothesis.** A treatment axis (v3) cannot test a stiffness
  prediction; a controlled stiffness contrast (v4) is required. Increased sample size does not
  compensate for the absence of the relevant experimental variable.
- **Signature scores export with donor and condition keys** so they can be aggregated and
  bootstrapped. Treat any difference found in a small subpopulation as provisional until it has
  survived a look at the donor level.

## Open next steps

- **Repeat v4 in an EGFR context.** A transcriptome-wide soft-vs-stiff dataset in PC9 or HCC827
  (EGFR exon19del) would put the model's own population on the footing A549 currently provides for
  KRAS. Martino et al. 2024 point the right way (YAP, CTGF, CYR61, and Ki67 all rising with
  stiffness) but only through targeted qPCR and protein, not transcriptome-wide.
- **The bimodality prediction, still untested.** The model predicts a bimodal G1 distribution, a
  cycling attractor and a quiescent one, on soft substrates under EGFR inhibition. The Cosgrove
  scRNA-seq CRISPRi series (GSE243756) could test whether the soft/stiff contrast produces that
  two-population structure, using the non-targeting guides for a clean contrast and the targeted
  guides as a perturbation. It needs real parsing work (guide demultiplexing), and it is the most
  distinctive claim the model has left.
- **Spatial NSCLC.** Define stiff and soft regions by local CAF and ECM density, then compare the
  same tumor's cancer cells across them. That is the v4 contrast run in situ, without the hydrogel.

## Summary of evidence

Fibrotic lung and NSCLC stroma share a CTHRC1-associated ECM-remodeling program consistent with a
mechanically altered microenvironment. The fibroblast-restricted analyses show that this signal is
not explained solely by immune or epithelial composition, although differences among fibroblast
states remain unresolved. In a controlled hydrogel experiment, A549 cells show concordant changes
in prespecified cell-cycle targets downstream of the mechanosensitive inputs represented in the
model.

The Maynard treatment trajectory independently characterizes a drug-tolerant persister state but
does not contain a stiffness contrast. The Kim malignant-cell analysis is limited by depth-dependent
CNV scoring. Together, the datasets support distinct stromal and tumor-cell components of the
proposed mechanism without directly validating the complete causal chain from stiffness through
FAK/AKT/ERK activity to EGFR-TKI tolerance. EGFR-specific transcriptomic confirmation under
controlled stiffness remains a central unresolved test.
