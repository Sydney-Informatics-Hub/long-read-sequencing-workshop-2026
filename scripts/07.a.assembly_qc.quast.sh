#!/bin/bash

set -euo pipefail

ASSEMBLY="${1}"
ID="${2}"

mkdir -p qc/assembly/quast/
mkdir qc/assembly/quast/${ID}

quast \
    ${ASSEMBLY} \
    --labels ${ID} \
    --output-dir qc/assembly/quast/${ID} \
    --threads 3
