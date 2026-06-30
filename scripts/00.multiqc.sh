#!/bin/bash
shopt -s expand_aliases

mkdir -p multiqc

multiqc \
    -o multiqc \
    -f \
    --fullnames \
    qc
