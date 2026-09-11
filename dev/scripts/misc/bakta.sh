#!/bin/bash
shopt -s expand_aliases
module load bakta

set -euo pipefail

ASSEMBLY="${1}"
ID="${2}"
mkdir -p annotate/${ID}

bakta ${ASSEMBLY} \
    --db data/ref/bakta \
    --output annotate/${ID}/ \
    --prefix ${ID} \
    --force \
    --threads 3
