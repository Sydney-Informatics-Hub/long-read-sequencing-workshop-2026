#!/bin/bash
shopt -s expand_aliases
module load flye

set -euo pipefail

FASTQ="${1}"
ID=$(basename ${FASTQ} .fastq)

mkdir -p assembly/${ID}/flye

flye \
    --nano-raw ${FASTQ} \
    --threads 2 \
    --out-dir assembly/${ID}/flye