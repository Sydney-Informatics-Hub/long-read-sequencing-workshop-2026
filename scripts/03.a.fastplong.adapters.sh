#!/bin/bash
shopt -s expand_aliases
module load fastplong

set -euo pipefail

FASTQ="${1}"
ID=$(basename ${FASTQ} .fastq.gz)

mkdir -p trim
mkdir -p qc/fastplong

fastplong \
    -i ${FASTQ} \
    --disable_quality_filtering \
    --disable_length_filtering \
    --trim_poly_x \
    -o trim/${ID}.trimmed.fastq.gz \
    -h qc/fastplong/${ID}.fastplong_report.html \
    -j qc/fastplong/${ID}.fastplong_report.json
