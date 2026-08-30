#!/bin/bash

set -euo pipefail

SAMPLE_LIST="$1"

if [[ ! -f "$SAMPLE_LIST" ]]; then
    echo "ERROR: Sample list file '$SAMPLE_LIST' not found."
    exit 1
fi

mkdir -p species_id

while read -r FASTQ; do
    # Skip empty lines
    [[ -z "$FASTQ" ]] && continue

    if [[ ! -f "$FASTQ" ]]; then
        echo "WARNING: File '$FASTQ' not found. Skipping."
        continue
    fi

    ID=$(basename ${FASTQ} .fastq.gz)
    K2DB=/home/tdev2/data/ref/kalamari

    mkdir -p species_id

    # Bind /cvmfs so singularity can find the kraken2 database
    export SINGULARITY_COMMAND_OPTS="-B /cvmfs"

    echo "Processing sample: $ID with Kraken2"
    kraken2 \
        ${FASTQ} \
        --db ${K2DB} \
        --report species_id/${ID}.k2report \
        --report-minimizer-data \
        --minimum-hit-groups 3 \
        --threads 2 \
        --output species_id/${ID}_k2_out.txt

done < "$SAMPLE_LIST"
