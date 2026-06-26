#!/bin/bash

set -euo pipefail

ASSEMBLY="${1}"
ID="${2}"

mkdir -p qc/assembly/bandage/
mkdir qc/assembly/bandage/${ID}

Bandage image \
    ${ASSEMBLY} \
    qc/assembly/bandage/${ID}/${ID}.assembly.svg
