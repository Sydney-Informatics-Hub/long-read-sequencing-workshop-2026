#!/bin/bash
shopt -s expand_aliases
module load bakta

set -euo pipefail

ASSEMBLY="${1}"
ID="${2}"
mkdir -p annotate/${ID}

bakta ${ASSEMBLY} \
    --db /cvmfs/data.galaxyproject.org/byhand/bakta_database/10522951 \
    --output annotate/${ID}/ \
    --prefix ${ID} \
    --force \
    --threads 3\
