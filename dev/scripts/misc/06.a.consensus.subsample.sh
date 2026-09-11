#!/bin/bash
shopt -s expand_aliases
module load autocycler

# Subsampling your reads prior to assembly (with multiple assemblers)
# is recommended when using Autocycler to construct a consensus assembly

FASTQ="${1}"
ID=$(basename ${FASTQ} .fastq)
OUTDIR=subsample/${ID}
mkdir -p ${OUTDIR}

set -euo pipefail

autocycler subsample \
    --reads ${FASTQ} \
    --out_dir ${OUTDIR} \
    --count 4 \
    --genome_size 5000000

# Rename FASTQs from sample_xx.fastq
# to ${ID}.sample_xx.fastq
for f in ${OUTDIR}/*.fastq; do
    mv ${f} ${OUTDIR}/${ID}.$(basename ${f})
done
