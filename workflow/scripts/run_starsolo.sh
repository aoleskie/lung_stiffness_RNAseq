#!/usr/bin/env bash
# Usage: run_starsolo.sh FASTQ_DIR INDEX WHITELIST OUT_PREFIX CB_LEN UMI_LEN THREADS
# 10x convention: R1 = cell barcode + UMI, R2 = cDNA. Adjust the globs if your
# files are named differently (e.g. *_1.fastq.gz / *_2.fastq.gz).
set -euo pipefail
FASTQ_DIR=$1; IDX=$2; WL=$3; OUT=$4; CBLEN=$5; UMILEN=$6; THREADS=$7

R1=$(ls ${FASTQ_DIR}/*_R1_*.fastq.gz | sort | paste -sd, -)
R2=$(ls ${FASTQ_DIR}/*_R2_*.fastq.gz | sort | paste -sd, -)
mkdir -p ${OUT}

STAR --runMode alignReads --genomeDir ${IDX} --runThreadN ${THREADS} \
     --readFilesCommand zcat --outFileNamePrefix ${OUT} \
     --readFilesIn ${R2} ${R1} \
     --soloType CB_UMI_Simple --soloCBwhitelist ${WL} \
     --soloCBstart 1 --soloCBlen ${CBLEN} \
     --soloUMIstart $((CBLEN+1)) --soloUMIlen ${UMILEN} \
     --soloFeatures Gene --outSAMtype BAM Unsorted

# Compress matrix outputs so scanpy.read_10x_mtx finds the expected .gz files.
gzip -f ${OUT}Solo.out/Gene/filtered/*.tsv ${OUT}Solo.out/Gene/filtered/*.mtx 2>/dev/null || true
echo "STARsolo done -> ${OUT}Solo.out/Gene/filtered/"
