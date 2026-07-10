#!/bin/bash

# Usage: ./qc_after.sh samples.txt

SAMPLE_LIST="$1"

if [[ ! -f "$SAMPLE_LIST" ]]; then
    echo "ERROR: Sample list file '$SAMPLE_LIST' not found."
    exit 1
fi

mkdir -p qc_after

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

    for PREFIX in filtered_ subsampled_; do
        NEW_ID="${PREFIX}${ID}"

        # FASTQC
        mkdir -p qc_after/fastqc/"$NEW_ID"
        fastqc \
            -f fastq \
            -o qc_after/fastqc/"$NEW_ID" \
            "filter/$NEW_ID.fastq.gz"

        # NANOPLOT
        mkdir -p qc_after/nanoplot/"$NEW_ID"
        NanoPlot \
            --fastq "filter/$NEW_ID.fastq.gz" \
            -p "${NEW_ID}_" \
            --loglength \
            -o qc_after/nanoplot/"$NEW_ID"/
    done


done < "$SAMPLE_LIST"

# MULTIQC
multiqc \
    -o qc_after/multiqc \
    -f \
    --fullnames \
    qc_after

echo "QC pipeline completed."
