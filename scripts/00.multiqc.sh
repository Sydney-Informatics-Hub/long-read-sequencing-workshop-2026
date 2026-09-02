#!/bin/bash
shopt -s expand_aliases
module load multiqc

mkdir -p multiqc

multiqc \
    -o multiqc \
    -f \
    --fullnames \
    qc
