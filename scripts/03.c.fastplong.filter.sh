#!/bin/bash
shopt -s expand_aliases
module load fastplong

set -euo pipefail

FASTQ="${1}"
ID=$(basename ${FASTQ} .fastq)

mkdir -p trim_filter
mkdir -p qc/fastplong

fastplong \
    -i ${FASTQ} \
    -q 13 \
    -o trim_filter/${ID}.trimmed.filtered.fastq \
    -h qc/fastplong/${ID}.fastplong_filter_report.html \
    -j qc/fastplong/${ID}.fastplong_filter_report.json
