#!/bin/bash
shopt -s expand_aliases
module load filtlong

set -euo pipefail

FASTQ="${1}"
ID=$(basename ${FASTQ} .fastq.gz)

mkdir -p filter

filtlong \
    --min_length 1kb \
    --keep_percent 50 \
    ${FASTQ} \
    | gzip -c > filter/${ID}.filtered.fastq.gz
