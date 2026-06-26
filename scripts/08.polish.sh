#!/bin/bash

set -euo pipefail

FASTQ="${1}"
ASSEMBLY="${2}"
ID=$(basename ${FASTQ} .fastq)

mkdir -p assembly/${ID}/medaka

medaka_consensus \
    -i ${FASTQ} \
    -d ${ASSEMBLY} \
    -o assembly/${ID}/medaka \
    -t 3
