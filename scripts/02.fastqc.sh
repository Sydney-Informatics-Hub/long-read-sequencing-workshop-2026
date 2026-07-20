#!/bin/bash
shopt -s expand_aliases
module load fastqc

set -euo pipefail

FASTQ="${1}"
ID=$(basename ${FASTQ} .fastq)

mkdir -p qc/fastqc/${ID}

fastqc \
    -f fastq \
    -o qc/fastqc/${ID} \
    --memory 2048 \
    ${FASTQ}

