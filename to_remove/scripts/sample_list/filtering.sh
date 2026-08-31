#!/bin/bash

# Usage: ./qc_before.sh samples.txt

SAMPLE_LIST="$1"

if [[ ! -f "$SAMPLE_LIST" ]]; then
    echo "ERROR: Sample list file '$SAMPLE_LIST' not found."
    exit 1
fi

mkdir -p filter

while read -r FASTQ; do
    # Skip empty lines
    [[ -z "$FASTQ" ]] && continue

    if [[ ! -f "$FASTQ" ]]; then
        echo "WARNING: File '$FASTQ' not found. Skipping."
        continue
    fi

    ID=$(basename "$FASTQ" .fastq.gz)
    FOLD=$(dirname "$FASTQ")

    gunzip "$FASTQ"
  
    echo "Processing sample: $ID"

    mkdir -p filter
    # NanoFilt -q 7 --length 1000 < "$FOLD/$ID.fastq" > "filter/filtered_$ID.fastq"
    filtlong \
        --min_length 1kb \
        --keep_percent 50 \
        "$FOLD/$ID.fastq" > "filter/$ID.filtered.fastq"

    mkdir -p trim_filter
    mkdir -p fastplong_filter

    fastplong \
        --in ${FASTQ} \
        --qualified_quality_phred 13 \
        --mean_qual 13 \
        --low_complexity_filter \
        --complexity_threshold 10 \
        --out trim_filter/${ID}.trimmed.filtered.fastq.gz \
        --html fastplong_filter/${ID}.fastplong_report.html \
        --json fastplong_filter/${ID}.fastplong_report.json
    
    gzip "filter/$ID.filtered.fastq"
    gzip "$FOLD/$ID.fastq"

done < "$SAMPLE_LIST"
