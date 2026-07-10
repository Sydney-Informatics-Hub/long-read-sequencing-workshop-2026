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

    echo "Processing sample: $ID with Flye"

    mkdir -p assembly/${ID}/flye

    flye \
        --nano-hq ${FASTQ} \
        --threads 2 \
        --out-dir assembly/${ID}/flye

    DB=/home/tdev2/data/ref/plassembler/plasmid_db_plassembler

    echo "Processing sample: $ID with Plassembler"

    mkdir -p assembly/${ID}/plassembler

    plassembler-exec plassembler long \
        -d ${DB} \
        -l ${FASTQ} \
        --flye_assembly assembly/${ID}/flye/assembly.fasta \
        --flye_info assembly/${ID}/flye/assembly_info.txt \
        -t 2 \
        -f \
        -o assembly/${ID}/plassembler

done < "$SAMPLE_LIST"
