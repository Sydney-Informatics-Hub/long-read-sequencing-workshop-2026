#!/bin/bash
shopt -s expand_aliases
module load autocycler

set -euo pipefail

SAMPLE="${1}"

# Gather assemblies into a single directory
CONSENSUSDIR=assembly/${SAMPLE}/consensus
ASSEBMLYDIR=${CONSENSUSDIR}/assemblies
mkdir -p ${ASSEBMLYDIR}
for ASSEMBLER in flye plassembler; do  # Add others if they were run
    # Assume there is only one FASTA in the folder
    ASSEMBLY=$(find assembly/${SAMPLE}/${ASSEMBLER}/ -maxdepth 1 -type f -name "*.fasta" | head -n 1)
    # Skip FASTA if empty
    if [ ! -s ${ASSEMBLY} ]; then continue; fi
    ln -s ${PWD}/${ASSEMBLY} ${ASSEBMLYDIR}/${SAMPLE}.${ASSEMBLER}.fasta
done

# Step 1: compress
AUTOCYCLERDIR=${CONSENSUSDIR}/autocycler
mkdir -p ${AUTOCYCLERDIR}
autocycler compress \
    --assemblies_dir ${ASSEBMLYDIR} \
    --autocycler_dir ${AUTOCYCLERDIR} \
    --threads 3

# Step 2: cluster
autocycler cluster \
    --autocycler_dir ${AUTOCYCLERDIR} \
    --min_assemblies 1

for d in ${AUTOCYCLERDIR}/clustering/qc_pass/cluster_*; do
    # Step 3: trim
    autocycler trim \
        --cluster_dir ${d} \
        --threads 3

    # Step 4: resolve
    autocycler resolve \
        --cluster_dir ${d}
done

# Step 5: combine
autocycler combine \
    --autocycler_dir ${AUTOCYCLERDIR} \
    --in_gfas ${AUTOCYCLERDIR}/clustering/qc_pass/cluster_*/5_final.gfa

# Step 6: summary table
autocycler table > ${CONSENSUSDIR}/autocycler_metrics.tsv  # Header row
autocycler table \
    --autocycler_dir ${AUTOCYCLERDIR} \
    --name consensus >> ${CONSENSUSDIR}/autocycler_metrics.tsv
