#!/bin/bash

set -euo pipefail

FASTQ="${1}"
ID=$(basename ${FASTQ} .fastq)

mkdir -p assembly/${ID}/flye

flye \
    --nano-hq ${FASTQ} \
    --threads 2 \
    --out-dir assembly/${ID}/flye