# Per-sample QC: load matrix, metrics, filter, doublets -> filtered .h5ad
rule qc:
    input:  lambda wc: samples.loc[wc.sample, "path"]
    output: f"{RESULTS}/qc/{{sample}}.h5ad"
    params:
        sample=lambda wc: wc.sample,
        condition=lambda wc: samples.loc[wc.sample, "condition"],
        donor=lambda wc: samples.loc[wc.sample, "donor_id"],
    log:    f"{RESULTS}/qc/logs/{{sample}}.log"
    conda:  "../envs/scanpy.yaml"
    script: "../scripts/qc.py"
