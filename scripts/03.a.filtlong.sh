#!/bin/bash
shopt -s expand_aliases
module load filtlong

set -euo pipefail

FASTQ="${1}"
ID=$(basename ${FASTQ} .fastq)

mkdir -p filter

filtlong \
    --min_length 1kb \
    --target_bases 50mb \
    ${FASTQ} \
    > filter/${ID}.filtered.fastq
