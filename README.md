# The stiffening lung and NSCLC drug tolerance

A lung tumor can build the ground it stands on. Fibroblasts in and around the tumor lay down collagen
and crosslink it until the tissue is measurably stiffer than the lung it replaced, and the working
claim is that cancer cells living in that stiffened tissue are harder to kill with an EGFR
inhibitor. This repository is where I went to find out how much of that four public datasets will
actually support.

There were two reasons to build it. I wanted a single-cell workflow I had assembled myself rather
than read about, and I wanted to put a dynamical-systems model of NSCLC drug tolerance in front of
data it had never seen. The thesis under test: a stiff, fibrotic microenvironment increases drug
tolerance by raising FAK/AKT/ERK activity in tumor cells, which sustains proliferation through
cyclin D–CDK4/6 and promotes EGFR-TKI tolerance. Four datasets, one Snakemake pipeline.

## The model being tested

The mitosis-apoptosis model (soon to be submitted!) takes **AKT, ERK, and FAK as three independent
measured inputs**: kinase activities set per (stiffness, drug) condition from phospho-Westerns,
feeding a cell-cycle/apoptosis ODE system. Two things follow from that structure, and both of them
limit what transcriptomics is allowed to say here.

1. **The kinase inputs are phospho-regulated, not transcriptional.** No amount of scRNA-seq will
   measure AKT/ERK/FAK activity; the manuscript measures those by Western. What transcriptomics can
   test is the layer underneath, the cell-cycle species the ODE rate laws drive when the inputs
   rise (CCND1, MYC, SKP2 up; CDKN1A down; the E2F program following a step later).
2. **YAP/TAZ is not a node in this model.** Mechanotransduction routes through FAK/AKT/ERK into the
   cell cycle. YAP/TAZ targets are scored throughout as an external mechanobiology readout, partly
   because the Cosgrove paper centers on them, and never as evidence for this model.

## The four datasets

- **v1, fibrosis (Tsukui 2020, E-CURD-126 / GSE132771).** Recover the CTHRC1+ pathologic fibroblast
  population and quantify the matrix-stiffening program (collagens plus the LOX/LOXL2 crosslinkers)
  that expands in fibrotic lung. *Result: the stiffening niche is there, and it is cell-intrinsic.
  The pathologic fraction separates fibrotic from normal donors at p = 0.018, which is the floor a
  5-vs-3 comparison can report; every fibrotic donor sits above every normal one. Pseudobulk DE
  then separates the compositional effect from the cell-intrinsic one.*
- **v2, the NSCLC microenvironment (Kim 2020, GSE131907).** Call malignant cells by CNV and
  characterize the tumor microenvironment. *Result: the CAF bridge. Tumor fibroblasts
  cell-intrinsically reactivate the CTHRC1+ stiffening program relative to normal lung fibroblasts
  (p ≈ 0, below what this sample size can resolve), and fibroblast-restricted PyDESeq2 confirms the
  activation is per-cell rather than compositional. Tumors build the stiff niche. The malignant
  cells themselves are a different story: their mechanotransduction signatures run flat to lower,
  and malignant calling is confounded by sequencing depth (r = −0.79), a droplet artifact that caps
  how far the tumor-intrinsic claim can go.*
- **v3, the treatment trajectory (Maynard 2020, PRJNA591860, Smart-seq2).** EGFR/ALK lung
  adenocarcinoma sampled from treatment-naive through residual disease (the persister state) to
  progression. The cleanest CNV calling in the project, since Smart-seq2 sequences every cell to
  comparable depth and the confound that limited Kim never appears. *Result: real persister
  biology. The alveolar-regenerative (AT2) program switches on at residual disease and cyclin
  D–CDK4/6 collapses, which is persister quiescence. The treatment axis, though, bundles drug
  exposure with clonal selection with elapsed time, and it never varies stiffness at all, so it
  cannot test the stiffness prediction. The cell-level mechano-signatures turn out to be
  pseudoreplication and do not survive aggregation to the patient (YAP, RD vs TN, p ≈ 1.0 across
  N = 5). Reframed as a treatment-trajectory result rather than a stiffness test.*
- **v4, controlled stiffness (Cosgrove 2024, *Science*; GSE243763, bulk RNA-seq).** A549 lung
  adenocarcinoma grown on soft (1 kPa) and stiff (50 kPa) hydrogels, 3 replicates per condition.
  *Result: the model's downstream cell-cycle targets shift together in the predicted direction on
  stiff substrates (CCND1/MYC/SKP2 up, CDKN1A down, the E2F program following; collective sign test
  p = 0.016, oriented Wilcoxon against a genome-wide background p = 0.013). No single gene clears
  significance at 3 vs 3, so what this is measuring is a coordinated set-level shift, not a gene.
  A separate and much stronger YAP/TAZ response (padj 0.0017) reproduces the source paper but sits
  outside the model. Caveat worth its weight: A549 is KRAS-mutant, so this places the mechanism in
  lung adenocarcinoma generally, not in the EGFR context the model is ultimately about.*

**The arc.** Tumors build the stiff niche, and patient data shows it (v1 and v2). Cancer cells
respond to stiffness with the cell-cycle program the model predicts, and a controlled contrast
shows that (v4). The treatment trajectory (v3) belongs to a different axis and is presented as one.
The kinase-input layer is tested by the manuscript's own phospho-Westerns; transcriptomics never
gets a vote on it.

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

## Environment notes (Windows 11 + Anaconda; read this part first)

I developed the pipeline on Windows, which I can now report was a decision. A handful of
platform-specific problems are worth knowing about before you start, and the shortest path around
every one of them is to run this on Linux :)

- **`--use-conda` is unreliable here.** The mamba activation banner pollutes the JSON that
  Snakemake's Python version-probe expects to parse. Workaround: `pip install snakemake` into the
  hashed scanpy conda env, activate that env, and run Snakemake without `--use-conda`. Every rule
  from QC through CNV and DE shares that one env, which already carries scanpy, infercnvpy,
  pydeseq2, and pyensembl.
- **inferCNVpy on Windows** spawns subprocesses (tqdm's `process_map`) that re-import the Snakemake
  script and crash. `infer_cnv.py` patches `process_map` with a serial map; setting `n_jobs=1` on
  its own does not fix it. CNV scoring needs the full chain `cnv.tl.pca → cnv.pp.neighbors →
  cnv.tl.leiden → cnv.tl.cnv_score`.
- **The HVG subset is a trap.** `merge_normalize` writes `annotated.h5ad` with only the 2,000 HVGs.
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
  Fibroblast-restricted DE runs through `caf_de.py` (PyDESeq2; the author's label is `Fibroblasts`,
  plural, which will bite you once).
- **v3 Maynard:** `ingest_maynard.py` chunk-loads the dense CSV, strips ERCC spike-ins, and maps the
  `analysis` column naive/grouped_pr/grouped_pd onto TN/RD/PD. Cell identity lives in the `cell_id`
  column, not the row index. Smart-seq2 means `doublet_method: none` and `hvg_flavor: seurat`, and
  since mito genes are stripped from the matrix the mito QC gate is disabled. Then the standard CNV
  path.
- **v4 Cosgrove:** `ingest_cosgrove.py` assembles the 9 A549 RSEM `*.genes.results` files from
  `GSE243763_RAW.tar`, strips Ensembl version suffixes, maps ENSG→symbol (pyensembl, with an
  HFF-table fallback), runs PyDESeq2 stiff-vs-soft, and scores the signatures. Then
  `check_model_targets.py` runs the per-gene and collective tests on the model's downstream
  cell-cycle targets. Bulk data at 3 vs 3, so this is the simplest installment by a wide margin:
  no CNV, no integration.

## Getting the data

- **v1 Tsukui:** the SCEA download tab (E-CURD-126) or GEO GSE132771, into `resources/tsukui/`.
- **v2 Kim:** the GSE131907 processed matrix (`..._raw_UMI_matrix.txt.gz` plus
  `cell_annotation.txt.gz`), into `resources/kim/`.
- **v3 Maynard:** processed CSVs from the study's Google Drive (`Data_input/csv_files/S01_datafinal.csv`
  and `S01_metacells.csv`), linked from `github.com/czbiohub-sf/scell_lung_adenocarcinoma`. Pull them
  with `gdown --folder`. The raw SRA under PRJNA591860 is far too heavy for what you get; use the
  processed CSVs.
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
columns it detected and the disease→condition mapping it inferred, which are worth reading before
you trust them; override with `--condition-col`, `--donor-col`, or `--id-col` if it guessed wrong.

**Gene symbols.** SCEA's `.mtx_rows` is often Ensembl-only, which quietly breaks both `MT-` mito QC
and every symbol-based signature. Resolve them with `--gtf path/to/annotation.gtf` (offline,
reproducible, and my recommendation; use the release the data was built on) or `--use-mygene`
(network, current annotations). Unmapped genes keep their Ensembl ID, and the originals are
preserved in `adata.var['gene_ids']` either way.

## How to work with it

Each `scripts/*.py` is either a Snakemake `script:` step, with the `snakemake` object injected, or a
standalone CLI (the v4 scripts). Treat them as starting points. I ran them rule by rule and refined
them while looking at the outputs, which is the workflow they are shaped for. The notebooks
regenerate each installment's figures from the `results_<name>/` outputs.

## Statistics conventions

These are the rules the project runs on, and most of them exist because ignoring one of them cost me
a result I briefly believed.

- **DE and composition tests run on per-donor pseudobulk.** The donor is the replication unit, never
  the cell. A cell-level p-value across a handful of donors is pseudoreplication wearing a lab coat.
- **Pair every cell-level signature claim with a per-donor check.** This is not optional. In v3 the
  dramatic cell-level signature shifts (KW p < 1e-100, which looks unanswerable) very nearly
  vanished at the patient level across N = 5 donors. The per-donor view is what the result actually
  was; the cell-level view was an artifact of counting thousands of non-independent cells.
- **For a multi-gene prediction, test the set, not each gene.** With small replicate counts no
  single gene will reach significance, and demanding that it does asks the wrong question. Test
  whether the predicted gene set moved coordinately: a sign test plus an oriented-LFC Wilcoxon
  against background, as in `check_model_targets.py`.
- **Match the axis of variation to the hypothesis.** A treatment axis (v3) cannot test a stiffness
  prediction however beautifully powered it is; only a stiffness contrast (v4) can. Sample size
  does not rescue the wrong instrument.
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

## Where this landed

Across four datasets, the stromal half of the model holds up and the tumor-cell half holds up in a
controlled contrast but not in patient tissue. Fibrotic lung and NSCLC stroma both carry the same
CTHRC1+ collagen program, and in both it survives the test that separates "more of these cells" from
"these cells doing more." Cancer cells on stiff hydrogels shift the model's downstream cell-cycle
targets coordinately in the predicted direction. The two patient datasets could never have shown the
second thing, because neither one contains a soft-versus-stiff comparison, and saying so is most of
what v2 and v3 contribute.

That is a narrower result than the project set out to get, and I think it is the honest shape of it.
The value is in knowing which datasets can answer which question, which is a cheaper lesson to learn
here than in a manuscript.
