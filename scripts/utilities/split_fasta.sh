#!/bin/bash
shopt -s expand_aliases
module load seqkit

ASSEMBLY="${1}"
ID="${2}"

mkdir -p assembly/${ID}/split/

seqkit split -i "${ASSEMBLY}" --by-id-prefix assembly/${ID}/split/
