#!/bin/bash
shopt -s expand_aliases
module load quast

set -euo pipefail

ASSEMBLY="${1}"
ID="${2}"

mkdir -p qc/assembly/quast/
mkdir -p qc/assembly/quast/${ID}

quast \
    ${ASSEMBLY} \
    --labels ${ID} \
    --output-dir qc/assembly/quast/${ID} \
    --threads 3
