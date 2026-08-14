#!/bin/bash
shopt -s expand_aliases
module load seqkit

ASSEMBLY="${1}"
ID="${2}"

seqkit split -i "${ASSEMBLY}"
