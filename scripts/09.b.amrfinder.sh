#!/bin/bash
shopt -s expand_aliases
module load ncbi-amrfinderplus

set -euo pipefail

ASSEMBLY="${1}"
ID="${2}"
mkdir -p annotate/${ID}/amrfinder

amrfinder \
    -n ${ASSEMBLY} \
    -d data/ref/amrfinderplus_db/2026-05-15.1 \
    --threads 3 > annotate/${ID}/amrfinder/${ID}.tsv
