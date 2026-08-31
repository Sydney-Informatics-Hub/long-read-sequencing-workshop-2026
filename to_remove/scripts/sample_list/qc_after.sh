#!/bin/bash

# Usage: ./qc_after.sh samples.txt

SAMPLE_LIST="$1"

if [[ ! -f "$SAMPLE_LIST" ]]; then
    echo "ERROR: Sample list file '$SAMPLE_LIST' not found."
    exit 1
fi

mkdir -p qc_after/fastqc
mkdir -p qc_after/nanoplot
mkdir -p qc_after/multiqc

while read -r FASTQ; do
    # Skip empty lines
    [[ -z "$FASTQ" ]] && continue

    if [[ ! -f "$FASTQ" ]]; then
        echo "WARNING: File '$FASTQ' not found. Skipping."
        continue
    fi

    ID=$(basename "$FASTQ")
    ID=${ID%.fastq.gz}
    ID=${ID%.fastq}

    echo "Processing sample: $ID"

    # FASTQC
    mkdir -p qc_after/fastqc/"$ID"
    fastqc \
        -f fastq \
        -o qc_after/fastqc/"$ID" \
        "$FASTQ"

    # NANOPLOT
    mkdir -p qc_after/nanoplot/"$ID"
    NanoPlot \
        --fastq "$FASTQ" \
        -p "${ID}_" \
        --loglength \
        --N50 \
        -o qc_after/nanoplot/"$ID"/

done < "$SAMPLE_LIST"

# MULTIQC
multiqc \
    -o qc_after/multiqc \
    -f \
    --fullnames \
    qc_after

echo "QC pipeline completed."
