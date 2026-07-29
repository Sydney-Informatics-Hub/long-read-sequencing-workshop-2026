#!/bin/bash

set -euo pipefail

ASSEMBLY="${1}"
ID="${2}"
CHROMID="${3}"
mkdir -p assembly/${ID}/plasmids

awk -v exclude="^>${CHROMID}" '$0 ~ /^>/ { p = 1 } $0 ~ exclude { p = 0 } p == 1 { print $0 }' ${ASSEMBLY} > assembly/${ID}/plasmids/${ID}.plasmids.fasta
