#!/bin/bash
shopt -s expand_aliases

FASTQ_IN=trim_filter/SRR11909877.trim_filter.fastq
ID=$(basename ${FASTQ_IN} .trim_filter.fastq)

# TODO
