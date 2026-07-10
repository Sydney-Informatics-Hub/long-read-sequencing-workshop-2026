#!/bin/bash

# Usage: ./qc_before.sh samples.txt

SAMPLE_LIST="$1"

if [[ ! -f "$SAMPLE_LIST" ]]; then
    echo "ERROR: Sample list file '$SAMPLE_LIST' not found."
    exit 1
fi

mkdir -p qc_before/fastqc
mkdir -p qc_before/nanoplot
mkdir -p qc_before/multiqc

while read -r FASTQ; do
    # Skip empty lines
    [[ -z "$FASTQ" ]] && continue

    if [[ ! -f "$FASTQ" ]]; then
        echo "WARNING: File '$FASTQ' not found. Skipping."
        continue
    fi

    ID=$(basename "$FASTQ" .fastq.gz)
    ID=$(basename "$ID" .fq.gz)   # also handle .fq.gz

    echo "Processing sample: $ID"

    # FASTQC
    mkdir -p qc_before/fastqc/"$ID"
    fastqc \
        -f fastq \
        -o qc_before/fastqc/"$ID" \
        "$FASTQ"

    # NANOPLOT
    mkdir -p qc_before/nanoplot/"$ID"
    NanoPlot \
        --fastq "$FASTQ" \
        -p "${ID}_" \
        --loglength \
        -o qc_before/nanoplot/"$ID"/

done < "$SAMPLE_LIST"

# MULTIQC
multiqc \
    -o qc_before/multiqc \
    -f \
    --fullnames \
    qc_before

echo "QC pipeline completed."
