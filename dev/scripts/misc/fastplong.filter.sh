#!/bin/bash
shopt -s expand_aliases
module load fastplong

set -euo pipefail

FASTQ="${1}"
ID=$(basename ${FASTQ} .fastq)

mkdir -p trim_filter
mkdir -p qc/fastplong_filter

fastplong \
    -i ${FASTQ} \
    -q 13 \
    --mean_qual 13 \
    --trim_poly_x \
    -o trim_filter/${ID}.trimmed.filtered.fastq \
    -h qc/fastplong_filter/${ID}.fastplong_report.html \
    -j qc/fastplong_filter/${ID}.fastplong_report.json
