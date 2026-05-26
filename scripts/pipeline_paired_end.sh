#!/bin/bash

set -euo pipefail

THREADS=8
GROUP="river"

RAW_DIR="data/raw/${GROUP}"
CLEAN_DIR="data/clean/${GROUP}"
QC_DIR="results/qc/${GROUP}"
ASSEMBLY_DIR="results/assembly/${GROUP}"
LOG_DIR="logs"

mkdir -p "$CLEAN_DIR" "$QC_DIR" "$ASSEMBLY_DIR" "$LOG_DIR"

for R1 in ${RAW_DIR}/*_1.fastq.gz
do
    SAMPLE=$(basename "$R1" _1.fastq.gz)
    R2="${RAW_DIR}/${SAMPLE}_2.fastq.gz"

    echo "======================================"
    echo "Processando amostra: $SAMPLE"
    echo "======================================"

    if [ ! -f "$R2" ]; then
        echo "ERRO: arquivo R2 não encontrado para $SAMPLE"
        echo "Pulando amostra..."
        continue
    fi

    echo "Rodando fastp..."

    fastp \
      -i "$R1" \
      -I "$R2" \
      -o "${CLEAN_DIR}/${SAMPLE}_R1.clean.fastq.gz" \
      -O "${CLEAN_DIR}/${SAMPLE}_R2.clean.fastq.gz" \
      --thread "$THREADS" \
      --detect_adapter_for_pe \
      --qualified_quality_phred 20 \
      --unqualified_percent_limit 30 \
      --length_required 50 \
      --html "${QC_DIR}/${SAMPLE}_fastp.html" \
      --json "${QC_DIR}/${SAMPLE}_fastp.json"

    echo "Rodando MEGAHIT..."

    megahit \
      -1 "${CLEAN_DIR}/${SAMPLE}_R1.clean.fastq.gz" \
      -2 "${CLEAN_DIR}/${SAMPLE}_R2.clean.fastq.gz" \
      -o "${ASSEMBLY_DIR}/${SAMPLE}_megahit" \
      --min-contig-len 1000 \
      -t "$THREADS"

    echo "Finalizado: $SAMPLE"
done

echo "Pipeline concluído para o grupo $GROUP."
