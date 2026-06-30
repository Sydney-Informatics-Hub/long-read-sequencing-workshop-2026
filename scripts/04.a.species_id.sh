#!/bin/bash
shopt -s expand_aliases

set -euo pipefail

FASTQ="${1}"
ID=$(basename ${FASTQ} .fastq)
K2DB=data/ref/kalamari

mkdir -p qc/species_id/report

# Bind /cvmfs so singularity can find the kraken2 database
export SINGULARITY_COMMAND_OPTS="-B /cvmfs"

kraken2 \
    ${FASTQ} \
    --db ${K2DB} \
    --report qc/species_id/report/${ID}.k2report \
    --report-minimizer-data \
    --minimum-hit-groups 3 \
    --threads 2 \
    --output qc/species_id/${ID}_k2_out.txt
