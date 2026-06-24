#!/bin/bash

SEQ_SUMMARY=data/sequencing_summary.txt

mkdir -p qc/nanoplot

NanoPlot \
    --summary ${SEQ_SUMMARY} \
    --loglength \
    -o qc/nanoplot/summary
