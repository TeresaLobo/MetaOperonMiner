#!/bin/bash

set -euo pipefail

THREADS=4
NREADS=5000000
GROUP="river"

CLEAN_DIR="data/clean/${GROUP}"
SUB_DIR="data/subsample/${GROUP}"
ASSEMBLY_DIR="results/assembly_subsample/${GROUP}"
LOG_DIR="logs"

mkdir -p "$SUB_DIR" "$ASSEMBLY_DIR" "$LOG_DIR"

for R1 in ${CLEAN_DIR}/*_R1.clean.fastq.gz
do
    SAMPLE=$(basename "$R1" _R1.clean.fastq.gz)
    R2="${CLEAN_DIR}/${SAMPLE}_R2.clean.fastq.gz"

    echo "===================================="
    echo "Subamostrando e montando: $SAMPLE"
    echo "===================================="

    if [ ! -f "$R2" ]; then
        echo "R2 ausente para $SAMPLE. Pulando."
        continue
    fi

    if [ ! -f "${SUB_DIR}/${SAMPLE}_R1.sub.fastq.gz" ]; then
        seqtk sample -s100 "$R1" "$NREADS" > "${SUB_DIR}/${SAMPLE}_R1.sub.fastq"
        seqtk sample -s100 "$R2" "$NREADS" > "${SUB_DIR}/${SAMPLE}_R2.sub.fastq"

        pigz -p "$THREADS" "${SUB_DIR}/${SAMPLE}_R1.sub.fastq"
        pigz -p "$THREADS" "${SUB_DIR}/${SAMPLE}_R2.sub.fastq"
    fi

    if [ -d "${ASSEMBLY_DIR}/${SAMPLE}_megahit" ]; then
        echo "Assembly já existe para $SAMPLE. Pulando."
        continue
    fi

    megahit \
      -1 "${SUB_DIR}/${SAMPLE}_R1.sub.fastq.gz" \
      -2 "${SUB_DIR}/${SAMPLE}_R2.sub.fastq.gz" \
      -o "${ASSEMBLY_DIR}/${SAMPLE}_megahit" \
      --min-contig-len 1000 \
      --memory 0.5 \
      -t "$THREADS"

    echo "Finalizado: $SAMPLE"
done
