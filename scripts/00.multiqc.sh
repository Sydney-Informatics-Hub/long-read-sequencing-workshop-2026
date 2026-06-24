#!/bin/bash

mkdir -p multiqc

multiqc \
    -o multiqc \
    -f \
    --fullnames \
    qc
