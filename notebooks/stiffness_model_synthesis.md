# Microenvironment stiffness and EGFR-mutant NSCLC drug tolerance: what the data support

*Synthesis across three patient scRNA-seq analyses (Tsukui fibrosis, Kim NSCLC TME, Maynard treatment-trajectory) and one controlled in-vitro stiffness contrast (Cosgrove A549 soft-vs-stiff RNA-seq).*

---

## 1. The model, and its two separable halves

The dynamical-systems model proposes that a stiff/fibrotic microenvironment raises FAK/AKT activity in EGFR-mutant tumor cells, which sustains proliferation (cyclin D–CDK4/6), shrinks the drug-sensitive slow-divider population, and thereby promotes tolerance to EGFR-TKIs. Cyclin D–CDK4/6 is framed as the lever: suppressing it should restore erlotinib sensitivity.

This single narrative contains **two logically distinct claims**, and they require different data to test:

- **(A) The source of stiffness.** Do tumors *build* a stiffer microenvironment? This is a statement about the stromal compartment — fibroblast activation, ECM deposition, collagen crosslinking.
- **(B) The cancer-cell response to stiffness.** Do tumor cells, when in a *stiffer* environment, raise FAK/AKT/YAP and proliferation relative to tumor cells in a *softer* one? This is a statement about a stiffness *contrast* in the malignant compartment, holding cell identity fixed.

Claim (A) is testable with dissociated patient scRNA-seq. **Claim (B) is not** — it requires comparing cancer cells in a *known soft* versus *known stiff* condition, with stiffness as the manipulated or measured variable, which no patient-biopsy dataset provides. Recognizing that boundary is the central methodological result of the patient-data work — and it is why a controlled in-vitro stiffness contrast (Section 2.4) was needed to close the gap. In patient tissue, stiffness is neither measured per cell nor varied independently of every other tumor property; in a hydrogel experiment it is the one thing that varies.

---

## 2. What each dataset established

### 2.1 Tsukui (fibrosis reference) — the pathologic stiffening niche exists and is cell-intrinsic

Eight donors (fibrotic vs normal lung), ~101K cells. Established four results, all per-donor / pseudobulk:

- A **CTHRC1+ pathologic fibroblast** population (cluster 4) carrying the collagen / ECM-stiffening program.
- **Fibrosis-associated expansion** of this population with reciprocal contraction of homeostatic fibroblasts.
- The **stiffening-ECM signature concentrates** in the pathologic fibroblast cluster.
- **Pseudobulk DE separates compositional from cell-intrinsic** effects: an immune confound (IGHG3) collapses to non-significance under fibroblast-restricted testing, while POSTN sharpens — i.e., the pathologic program is genuinely cell-intrinsic activation, not just a shift in cell proportions.

This is the reference definition of the stiffening niche. It is about fibroblasts, not cancer cells — relevant to claim (A), silent on claim (B).

### 2.2 Kim (NSCLC TME) — tumors rebuild the stiffening niche; cancer-cell mechanotransduction is *not* elevated

11 tumor + 11 normal samples, ~88K cells, treatment-naive. Cancer cells called by inferCNV.

**The robust, headline result — the CAF bridge (claim A, in cancer):** Tumor-associated fibroblasts cell-intrinsically reactivate the pathologic program. Compared to normal-lung fibroblasts, tumor CAFs score far higher on `pathologic_caf` (~0.76 vs −0.03, every tumor donor above nearly every normal donor) and on `stiffening_ecm` (0.52 vs 0.17, p≈0.003). Fibroblast-restricted pseudobulk DE (PyDESeq2) confirms cell-intrinsic CAF activation: CTHRC1 (LFC~3, padj~6e-5), COL3A1, LOXL2, COL1A1 all up, with the genes lying on/above the y=x line of whole-tissue-vs-fibroblast LFC (i.e., not a compositional artifact). **The Tsukui CTHRC1+ stiffening niche reappears in NSCLC tumor stroma.** Tumors build the stiff niche.

**The cancer-cell finding (bearing on claim B, and negative):** On malignant cells specifically, the mechanotransduction signatures are *not* elevated. `yap_taz_targets` is **lower** in malignant cells (p=1e-55), `pi3k_akt_axis` is flat, and apparent `stiffening_ecm` signal in malignant cells tracks stromal/EMT contamination rather than a tumor-cell-intrinsic program. So in treatment-naive bulk tumor, cancer cells do not show the elevated FAK/AKT/YAP the model predicts.

**Critical limitation:** Malignant calling in Kim is confounded by sequencing depth (per-donor malignant-rate vs depth r=−0.79); shallow normal samples were falsely called malignant because inferCNV reads low-depth noise as CNV. This is a droplet-scRNA-seq depth artifact, not fixable post hoc, and it caps the per-donor tumor-intrinsic claim — though it does not affect the CAF (stromal) result, which is the solid one.

### 2.3 Maynard (treatment trajectory) — real persister biology, but the *wrong axis* for stiffness

49 biopsies / 30 patients, EGFR/ALK/etc., Smart-seq2, sampled before treatment (TN), at residual disease (RD, the drug-tolerant persister state), and at progression (PD). Cancer cells called by inferCNV — and here the CNV diagnostics were the **cleanest in the project** (uniform Smart-seq2 depth, no depth confound; immune reference tight and low, malignant epithelium widely separated and high). EGFR-restricted malignant cells: TN 400, RD 478, PD 1053.

Scoring the model nodes across TN→RD→PD on EGFR malignant cells gave large, monotone cell-level effects:

| signature | TN | RD | PD | cell-level KW |
|---|---|---|---|---|
| at2_regenerative *(positive control)* | −0.15 | **+1.19** | −0.45 | p≈1e-294 |
| cdk46_cellcycle | 0.338 | −0.038 | −0.315 | p≈1e-232 |
| yap_taz_targets | 0.121 | 0.165 | −0.081 | p≈1e-110 |
| pi3k_akt_axis | −0.078 | −0.031 | 0.019 | p≈1e-21 |

The positive control fires decisively (residual-disease persisters adopt the alveolar-regenerative AT2-like state the original authors described), validating that the pipeline captures real persister biology.

**But two things must be stated plainly:**

**(i) The CDK4/6 collapse is the *opposite* of the stiffness prediction, and is persister biology, not a stiffness readout.** The model predicts stiffness *raises* cyclin D–CDK4/6 (more dividing cells). What Maynard shows is CDK4/6 *collapsing* in persisters — because cells under TKI exit the cell cycle to survive. That quiescence is the defining feature of drug-tolerant persistence and is orthogonal to matrix stiffness. Reading the CDK4/6 drop as relevant to the stiffness axis was an error of interpretation; the direction and the cause both belong to treatment response, not to a soft→stiff transition.

**(ii) The cell-level significances are pseudoreplication; almost nothing survives per-donor.** Aggregating cells to patient-stage means (RD n=5, TN n=5, PD n=3 EGFR patients) and testing across patients:

| signature | TN donors (median) | RD donors (median) | RD vs TN per-donor |
|---|---|---|---|
| cdk46_cellcycle | 0.255 | −0.041 | p=0.151 (largest effect, direction holds) |
| at2_regenerative | 0.595 | 1.214 | p=0.841 (direction holds, underpowered) |
| pi3k_akt_axis | −0.073 | −0.036 | p=0.421 |
| yap_taz_targets | 0.155 | 0.134 | **p=1.000 (peak vanishes)** |

The YAP peak at RD — the most model-specific signal in the cell-level view — **does not survive aggregation to the patient level.** Even the authors' own AT2 finding is directional-but-not-significant at N=5. Only the CDK4/6 drop holds directionally per-donor, and that, as noted, is persister quiescence rather than a stiffness effect.

**Why Maynard cannot test claim (B) even in principle:** its axis is *treatment stage*, which confounds drug exposure, clonal selection, persister survival, and time. No two stages differ *only* in stiffness; stiffness is never measured. So "signatures change across TN/RD/PD" — whatever they do — cannot be attributed to stiffness. Maynard is a treatment-trajectory dataset, and it answers treatment-trajectory questions well; it is simply not a stiffness experiment.

### 2.4 Cosgrove A549 (controlled in-vitro stiffness contrast) — the cancer-cell response, tested directly

Cosgrove, Bounds, Taylor et al. (*Science* 2024; GEO GSE243763), bulk RNA-seq. A549 lung adenocarcinoma on **soft (1 kPa) vs stiff (50 kPa)** polyacrylamide hydrogels (plus TCP), 3 biological replicates per condition. This is the contrast claim (B) actually requires: same cells, stiffness manipulated, everything else fixed — the controlled test the patient data structurally cannot provide.

**The assay boundary that frames the test.** The model (`cell_cycle_death_model.py`) takes **AKT, ERK, FAK as measured inputs** — kinase activities set per (stiffness, drug) condition from phospho-Westerns, *not* transcript levels. So bulk RNA-seq cannot test the input layer. What it *can* test is the model's **downstream species** — the cell-cycle players the ODE rate laws produce when the input knobs rise. The rate laws give a specific, pre-registered prediction set for stiff (higher AKT/ERK/FAK): **CCND1 up** (CycD synthesis driven by AKT & ERK), **MYC up** (Myc synthesis, AKT & ERK), **SKP2 up** (Skp2 synthesis, FAK), **CDKN1A down** (p21 repressed by FAK), with the downstream **E2F program (E2F1, CCNE1)** following.

**Result — the model's downstream targets shift coordinately as predicted (PyDESeq2, stiff vs soft):** all six predicted genes move the predicted way. CCND1 (LFC +0.39, 88th percentile genome-wide), SKP2 (+0.27, 82nd), MYC (+0.13, 69th) up; CDKN1A down (weakest, just below center); E2F1 and CCNE1 up. No single gene is significant at 3v3 (all padj > 0.05) — expected, and not the right test for a multi-gene prediction. The correct **set-level** tests are significant: **sign test 6/6 concordant p = 0.016**, and an **oriented-LFC Wilcoxon vs the genome-wide background** (each LFC multiplied by its predicted sign, then tested against the ~20K-gene background centered at ~0) **p = 0.013**. The result is robust to weighting direction alone (sign test) or direction-and-magnitude (Wilcoxon). The pattern also tracks the model's wiring: the most-directly-driven genes (CCND1, SKP2) show the largest shifts, the two-steps-downstream E2F genes the smallest — attenuation down the cascade, as predicted. CCND1 being the standout matters: it is the upstream node of the model's therapeutic lever (cyclin D–CDK4/6).

**A separate, strong finding that is *not* model support:** the **YAP/TAZ target signature** is significantly up on stiff (mean LFC +0.58, padj 0.0017), reproducing the source paper's headline. But YAP/TAZ is **not a node in this model**, which routes mechanotransduction through FAK/AKT/ERK → cell cycle. So this is a real, parallel mechanotransduction response, reported as its own result, not folded into the model test. (FAK→YAP crosstalk is documented in the mechanobiology literature and would be where a future model *extension* might connect them — but that is not a claim this dataset establishes.)

**Caveat:** A549 is **KRAS-mutant**, so this confirms the *mechanism* in lung adenocarcinoma, not the EGFR context specifically. PC9/HCC827 hydrogel data (Martino et al. 2024; raised YAP/CTGF/CYR61/Ki67 with stiffness) provide EGFR-specific support as literature.

---

## 3. The honest status of the model

- **Claim (A) — tumors build the stiff niche — is supported (patient data).** The CTHRC1+/collagen/crosslinking program is cell-intrinsically reactivated in NSCLC CAFs (Kim), reproducing the fibrosis reference (Tsukui). The project's solid, defensible patient-data result.

- **Claim (B) — cancer cells raise the proliferative/cell-cycle program *in response to stiffness* — is now supported by a controlled contrast (Cosgrove A549), at the level the assay can test.** The two patient datasets could *not* test it: in treatment-naive tumor (Kim) malignant-cell mechanotransduction signatures are flat-to-lower, and in the treatment trajectory (Maynard) the apparent signal does not survive per-donor — but neither contains a soft-vs-stiff contrast, so those are non-tests, not evidence against. The controlled in-vitro contrast (Cosgrove) does contain it, and there the model's downstream cell-cycle targets shift coordinately in the predicted direction (collective p ≈ 0.013). The kinase-*input* layer (AKT/ERK/FAK) is phospho-regulated and remains outside what RNA-seq can adjudicate — it is tested by the manuscript's own Western data, not by transcriptomics.

So the boundary has moved: **the model's stromal half is corroborated in patient tissue, and its cancer-cell-response half — the downstream cell-cycle consequences the ODEs predict — is now corroborated in a controlled stiffness contrast.** Two caveats bound the second claim: it is a *set-level* coordinated shift (no single gene significant at 3v3), and it is in **KRAS-mutant** A549 (the mechanism in lung adenocarcinoma; EGFR-specific support is literature, PC9/HCC827). A separate strong YAP/TAZ stiffness response in the same data is real mechanobiology but lies outside the model's wiring.

---

## 4. Remaining gaps and what would further strengthen the model

The controlled contrast closes the central gap (does the cancer-cell program respond to stiffness as predicted), but three things would sharpen it further, in decreasing priority:

### 4.1 The EGFR context (the standing caveat)
The Cosgrove test is in KRAS-mutant A549. The model is about EGFR-mutant disease. PC9/HCC827 hydrogel data (Martino et al. 2024; YAP/CTGF/CYR61/Ki67 up with stiffness across 0.5–32 kPa) show the same direction in the canonical EGFR exon19del lines, but appear to be targeted qPCR/protein, not transcriptome-wide — strong literature support, not a dataset scored here. A transcriptome-wide soft-vs-stiff dataset in an EGFR-mutant line (± EGFR-TKI) would put the model's own population on the same footing as the A549 result, and would let the *downstream-target* test be repeated in-context.

### 4.2 The heterogeneity / bimodality prediction (single-cell, untested)
The manuscript's most distinctive claim is not just a mean shift but a **bimodal G1 distribution** on soft + EGFRi substrates — a fast-cycling population and a slow-divider/quiescent peak, attributed to noise-driven switching between attractors. This is intrinsically a between-cell-variability claim that bulk RNA-seq averages away. The Cosgrove **scRNA-seq CRISPRi series (GSE243756)** could test it: cluster the cancer cells, score proliferation/CDK4-6 per cell, and ask whether the soft/stiff contrast produces the predicted bimodal cell-state structure (using non-targeting control guides for the clean contrast; targeted guides as a perturbation test of whether disrupting a mechanoenhancer collapses the proliferative sub-population). This is meaningfully more parsing work (guide demultiplexing) but is the natural next single-cell test and the one aligned with the model's most novel prediction.

### 4.3 Spatial transcriptomics of NSCLC (the in-vivo version)
Visium/Xenium/MERFISH on NSCLC tissue, defining stiff vs soft regions by *local ECM/CAF density* (dense CTHRC1+ collagen stroma = stiff proxy; loose stroma = soft), then asking whether cancer cells in stiff neighborhoods show higher proliferation/cell-cycle program than the *same tumor's* cells in soft neighborhoods. Keeps the malignant clone fixed and varies the local niche — the in-situ version of the contrast, bridging the controlled in-vitro result back to patient tissue.

### 4.4 Pathway perturbation (mechanistic complement)
FAK inhibition/knockdown vs control in EGFR-mutant cells ± TKI, asking whether disrupting the input shifts the drug-tolerant / cell-cycle state — directly relevant to the model's therapeutic claim (CDK4/6 as the lever). The Cosgrove CRISPRi screen is a ready-made version of this logic for the mechanoenhancer layer.

---

## 5. Where the project stands

The four datasets now assemble into a complete, honest arc, each claim stated at the confidence its data support:

- **Tumors build the stiff niche (patient data, supported).** NSCLC CAFs cell-intrinsically reactivate the CTHRC1+/collagen/crosslinking program (Kim), reproducing the fibrosis reference (Tsukui).
- **The treatment-trajectory data is orthogonal to stiffness (patient data, correctly bracketed).** Maynard shows genuine persister biology (AT2 reprogramming, CDK4/6 collapse) but on a treatment axis that confounds drug/selection/time and never measures stiffness; its mechano-signatures do not survive per-donor. Not a stiffness test, and not read as one.
- **Cancer cells respond to stiffness as the model's downstream predictions require (controlled contrast, supported at the set level).** In A549 on soft vs stiff hydrogels, the model's downstream cell-cycle targets shift coordinately in the predicted direction (CCND1/MYC/SKP2 up, CDKN1A down, E2F program following; collective sign test p = 0.016, oriented Wilcoxon vs background p = 0.013). CCND1 — the upstream node of the therapeutic lever — is the standout.
- **Separately, a strong YAP/TAZ stiffness response** (padj 0.0017) reproduces the source paper but lies outside the model's wiring; reported as parallel mechanobiology, not model support.

For a manuscript, this reads as: **tumors build the stiff niche (patient data); cancer cells respond to it with the cell-cycle program the model predicts (controlled in-vitro contrast); the kinase-input layer is tested by the manuscript's own phospho-Westerns, not by transcriptomics; and the treatment-trajectory data is honestly bracketed as a different axis.** The standing caveat — A549 is KRAS-mutant, so the controlled result is the mechanism in lung adenocarcinoma rather than the EGFR context — is the clearest target for the next experiment (Section 4.1), with the single-cell heterogeneity/bimodality prediction (Section 4.2) the most distinctive remaining untested claim.

> **One-line summary.** The model's stromal half is corroborated in patient tissue and its cancer-cell-response half (downstream cell-cycle targets) in a controlled soft-vs-stiff contrast; what remains open is the EGFR-specific repeat and the single-cell bimodality prediction — both with identified, analyzable datasets.
