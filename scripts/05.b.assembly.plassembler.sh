#!/bin/bash
shopt -s expand_aliases
module load plassembler

FLYE="${2}"

set -euo pipefail

FASTQ="${1}"

ID=$(basename ${FASTQ} .fastq)

DB=data/ref/plassembler

mkdir -p assembly/${ID}/plassembler

FLYE_ARGS=""
if [ -d "${FLYE}" ] && [ -f "${FLYE}/assembly.fasta" ] && [ -f "${FLYE}/assembly_info.txt" ]; then
    FLYE_ARGS="--flye_assembly ${FLYE}/assembly.fasta --flye_info ${FLYE}/assembly_info.txt"
fi

plassembler long \
    -d ${DB} \
    -l ${FASTQ} \
    ${FLYE_ARGS} \
    -t 2 \
    -f \
    -o assembly/${ID}/plassembler