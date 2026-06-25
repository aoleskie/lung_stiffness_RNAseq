# FASTQ -> counts with STARsolo (the "from raw data" lesson, one sample).
rule starsolo:
    input:
        idx=config["reference"]["star_index"],
        wl=config["reference"]["cb_whitelist"],
    output:
        f"{RESULTS}/quant/{{sample}}/Solo.out/Gene/filtered/matrix.mtx.gz",
    params:
        fastq_dir=lambda wc: samples.loc[wc.sample, "path"],
        out_prefix=lambda wc: f"{RESULTS}/quant/{wc.sample}/",
        cb_len=config["reference"]["cb_len"],
        umi_len=config["reference"]["umi_len"],
    threads: 8
    conda: "../envs/star.yaml"
    shell:
        "bash workflow/scripts/run_starsolo.sh "
        "{params.fastq_dir} {input.idx} {input.wl} {params.out_prefix} "
        "{params.cb_len} {params.umi_len} {threads}"
