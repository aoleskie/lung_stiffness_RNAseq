# Merge all QC'd samples -> normalize/HVG/PCA -> batch integration.
rule merge_normalize:
    input:  expand(f"{RESULTS}/qc/{{s}}.h5ad", s=MATRIX_SAMPLES)
    output: f"{RESULTS}/integrate/merged.h5ad"
    conda:  "../envs/scanpy.yaml"
    script: "../scripts/merge_normalize.py"

rule integrate:
    input:  f"{RESULTS}/integrate/merged.h5ad"
    output: f"{RESULTS}/integrate/integrated.h5ad"
    conda:  "../envs/scanpy.yaml"     # for scVI: switch to ../envs/scvi.yaml + a scvi script
    script: "../scripts/integrate_harmony.py"
