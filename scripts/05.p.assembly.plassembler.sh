#!/bin/bash

FASTQ_IN=data/fastqs/SRR11909877.sampled.2k.fastq
ID=$(basename ${FASTQ_IN} .sampled.2k.fastq)
DB=ref/plassembler

mkdir -p assembly/${ID}/plassembler


plassembler long \
    -d ${DB} \
    -l ${FASTQ_IN} \
    --flye_assembly assembly/${ID}/flye/assembly.fasta \
    --flye_info assembly/${ID}/flye/assembly.info \
    -t 2 \
    -f \
    -o assembly/${ID}/plassembler