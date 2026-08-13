#!/bin/bash

module load seqkit

ASSEMBLY="${1}"
ID="${2}"

mkdir -p assembly/split/${ID}

seqkit split -i "${ASSEMBLY}" --by-id-prefix assembly/split/${ID}
