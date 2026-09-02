#!/bin/bash
shopt -s expand_aliases
module load bandage

set -euo pipefail

ASSEMBLY="${1}"
ID="${2}"

mkdir -p qc/assembly/bandage/
mkdir -p qc/assembly/bandage/${ID}

bandage-exec \
Bandage image \
    ${ASSEMBLY} \
    qc/assembly/bandage/${ID}/${ID}.assembly.svg
