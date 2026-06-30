#!/bin/bash
shopt -s expand_aliases

set -euo pipefail

ASSEMBLY="${1}"
ID="${2}"

mkdir -p qc/assembly/busco/
mkdir qc/assembly/busco/${ID}

busco \
    --in ${ASSEMBLY} \
    --out qc/assembly/busco/${ID} \
    --mode genome \
    --lineage_dataset data/ref/busco/bacteria_odb12.2 \
    --force \
    --offline \
    --cpu 3
