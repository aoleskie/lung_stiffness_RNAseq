# Reference build for the STARsolo lesson only (skipped if no fastq lesson).
rule download_genome:
    output: f"{RESOURCES}/genome/genome.fa"
    params: url=config["reference"]["genome_fasta_url"]
    shell:  "mkdir -p $(dirname {output}) && curl -L {params.url} | gunzip -c > {output}"

rule download_gtf:
    output: f"{RESOURCES}/genome/annotation.gtf"
    params: url=config["reference"]["gtf_url"]
    shell:  "mkdir -p $(dirname {output}) && curl -L {params.url} | gunzip -c > {output}"

rule star_index:
    input:
        fa=f"{RESOURCES}/genome/genome.fa",
        gtf=f"{RESOURCES}/genome/annotation.gtf",
    output: directory(config["reference"]["star_index"])
    threads: 8
    conda:  "../envs/star.yaml"
    shell:
        "mkdir -p {output} && "
        "STAR --runMode genomeGenerate --genomeDir {output} "
        "--genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} "
        "--runThreadN {threads} --sjdbOverhang 90"
