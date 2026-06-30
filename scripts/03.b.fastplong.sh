#!/bin/bash
shopt -s expand_aliases

set -euo pipefail

FASTQ="${1}"
ID=$(basename ${FASTQ} .fastq)

mkdir -p trim
mkdir -p qc/fastplong

fastplong \
    -i ${FASTQ} \
    --disable_quality_filtering \
    --disable_length_filtering \
    -o trim/${ID}.trimmed.fastq \
    -h qc/fastplong/${ID}.fastplong_report.html \
    -j qc/fastplong/${ID}.fastplong_report.json
