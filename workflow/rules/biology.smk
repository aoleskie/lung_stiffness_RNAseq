# The biology: signature scoring + pathway activity, and pseudobulk DE.
rule signatures:
    input:  f"{RESULTS}/cluster/annotated.h5ad"
    output:
        scores=f"{RESULTS}/biology/signature_scores.csv",
        fig=f"{RESULTS}/figures/signature_umaps.png",
    conda:  "../envs/scanpy.yaml"
    script: "../scripts/score_signatures.py"

rule pseudobulk_de:
    input:  f"{RESULTS}/cluster/annotated.h5ad"
    output: f"{RESULTS}/biology/pseudobulk_de.csv"
    conda:  "../envs/pydeseq2.yaml"
    script: "../scripts/pseudobulk_de.py"

# v2 stretch (scaffold): malignant-cell calling by CNV. Wire into `all` for v2.
# rule infercnv:
#     input:  f"{RESULTS}/cluster/annotated.h5ad"
#     output: f"{RESULTS}/biology/cnv.h5ad"
#     conda:  "../envs/infercnv.yaml"
#     script: "../scripts/infercnv.py"
