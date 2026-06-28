# The biology: signature scoring + pathway activity, and pseudobulk DE.
rule signatures:
    input:  f"{RESULTS}/cluster/annotated.h5ad"
    output:
        scores=f"{RESULTS}/biology/signature_scores.csv",
        fig=f"{RESULTS}/figures/signature_umaps.png",
    conda:  "../envs/scanpy.yaml"
    script: "../scripts/score_signatures.py"

rule pseudobulk_de:
    input:  expand(f"{RESULTS}/qc/{{s}}.h5ad", s=MATRIX_SAMPLES)
    output: f"{RESULTS}/biology/pseudobulk_de.csv"
    conda:  "../envs/pydeseq2.yaml"
    script: "../scripts/pseudobulk_de.py"

rule pseudobulk_de_fibroblast:
    input: qc=expand(f"{RESULTS}/qc/{{s}}.h5ad", s=MATRIX_SAMPLES), annotated=f"{RESULTS}/cluster/annotated.h5ad",
    output: f"{RESULTS}/biology/pseudobulk_de_fibroblast.csv"
    conda:  "../envs/pydeseq2.yaml"
    script: "../scripts/pseudobulk_de_fibroblast.py"

rule infer_cnv:
    input: annotated = f"{RESULTS}/cluster/annotated.h5ad",
    output:
        h5ad   = f"{RESULTS}/cnv/cnv.h5ad",
        scores = f"{RESULTS}/cnv/cnv_scores.csv",
        fig    = f"{RESULTS}/figures/cnv_diagnostics.png",
    script: "../scripts/infer_cnv.py"

# v2 stretch (scaffold): malignant-cell calling by CNV. Wire into `all` for v2.
# rule infercnv:
#     input:  f"{RESULTS}/cluster/annotated.h5ad"
#     output: f"{RESULTS}/biology/cnv.h5ad"
#     conda:  "../envs/infercnv.yaml"
#     script: "../scripts/infercnv.py"
