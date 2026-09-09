#!/bin/bash

shopt -s expand_aliases

set -euo pipefail

# --- Defaults -----------------------------------------------------------
FLYE=""
PLASSEMBLER=""
KEEP=""
EXCLUDE=""
CHROM=""
OUTPUT=""

# --- Flag parsing ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --flye)
            FLYE="$2"
            shift
            shift
            ;;
        --plassembler)
            PLASSEMBLER="$2"
            shift
            shift
            ;;
        --keep)
            KEEP="$2"
            shift
            shift
            ;;
        --exclude)
            EXCLUDE="$2"
            shift
            shift
            ;;
        --chromosome)
            CHROM="$2"
            shift
            shift
            ;;
        --output)
            OUTPUT="$2"
            shift
            shift
            ;;
        -h|--help)
            echo "Usage: $0 --flye /path/to/flye/assembly.fasta --plassembler /path/to/plassembler/plassembler_plasmids.fasta --chromosome contig_1 [ --[keep|exclude] contig_1[,contig_2] ] --output /path/to/output.fasta"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

if [ ! -f "${FLYE}" ]; then
    echo "ERROR: Must provide a flye assembly."
    exit 1
fi

if [ ! -f "${PLASSEMBLER}" ]; then
    echo "ERROR: Must provide a plassembler assembly."
    exit 1
fi

if [ -z "${CHROM}" ]; then
    echo "ERROR: Must provide a name of a contig to mark as the chromosome."
    exit 1
fi

if [ -n "${KEEP}" ] && [ -n "${EXCLUDE}" ]; then
    echo "ERROR: --keep and --exclude cannot be used together."
    exit 1
fi

if [ -z "${OUTPUT}" ]; then
    echo "ERROR: Must provide an output FASTA file."
    exit 1
fi

if [ -f "${OUTPUT}" ]; then
    echo "ERROR: Output file '${OUTPUT}' already exists."
    exit 1
fi

TEMP_DIR=".merge_contigs.tmp"
mkdir "${TEMP_DIR}"

# Step 1 - concatenate FASTA files
echo "Step 1: Concatenate FASTA files"
cat "${FLYE}" "${PLASSEMBLER}" > "${TEMP_DIR}/concat.fasta"

# Step 2 - strip unwanted contigs
if [ -n "${KEEP}" ]; then
    echo "Step 2: Strip unwanted contigs"
    # Split by commas
    echo "${KEEP}" | tr , '\n' > "${TEMP_DIR}/contigs_to_keep.txt"
    # Extract contigs
    seqkit grep -f "${TEMP_DIR}/contigs_to_keep.txt" "${TEMP_DIR}/concat.fasta" > "${TEMP_DIR}/filtered.fasta"
elif [ -n "${EXCLUDE}" ]; then
    echo "Step 2: Strip unwanted contigs"
    # Split by commas
    echo "${EXCLUDE}" | tr , '\n' > "${TEMP_DIR}/contigs_to_exclude.txt"
    # Extract contigs
    seqkit grep -v -f "${TEMP_DIR}/contigs_to_exclude.txt" "${TEMP_DIR}/concat.fasta" > "${TEMP_DIR}/filtered.fasta"
else
    echo "SKIPPING: Step 2: Strip unwanted contigs"
    mv "${TEMP_DIR}/concat.fasta" "${TEMP_DIR}/filtered.fasta"
fi

# Step 3 - extract chromosome
echo "Step 3: Extract chromosome"
seqkit grep -p "${CHROM}" "${TEMP_DIR}/filtered.fasta" > "${TEMP_DIR}/chromosome.fasta"
seqkit grep -v -p "${CHROM}" "${TEMP_DIR}/filtered.fasta" > "${TEMP_DIR}/plasmids.fasta"

# Step 4 - rename contigs
echo "Step 4: Rename contigs"
seqkit replace -p "(.+)" -r "chromosome" "${TEMP_DIR}/chromosome.fasta" > "${TEMP_DIR}/chromosome.renamed.fasta"
seqkit replace -p "(.+)" -r "plasmid_{nr}" "${TEMP_DIR}/plasmids.fasta" > "${TEMP_DIR}/plasmids.renamed.fasta"

# Step 5 - assemble final FASTA
echo "Step 5: Assemble final FASTA"
cat "${TEMP_DIR}/chromosome.renamed.fasta" "${TEMP_DIR}/plasmids.renamed.fasta" > "${OUTPUT}"

# Step 6 - generate map of old to new IDs
echo "Step 6: Generate map of old to new IDs"
MAPFILE="${OUTPUT}.id_map.txt"
set +e
grep '^>' "${TEMP_DIR}/chromosome.fasta" | sed -E -e 's|^>||g' -e 's|[ \t].*||g' > "${TEMP_DIR}/chromosome.old"
grep '^>' "${TEMP_DIR}/chromosome.renamed.fasta" | sed -E -e 's|^>||g' -e 's|[ \t].*||g' > "${TEMP_DIR}/chromosome.new"
grep '^>' "${TEMP_DIR}/plasmids.fasta" | sed -E -e 's|^>||g' -e 's|[ \t].*||g' > "${TEMP_DIR}/plasmids.old"
grep '^>' "${TEMP_DIR}/plasmids.renamed.fasta" | sed -E -e 's|^>||g' -e 's|[ \t].*||g' > "${TEMP_DIR}/plasmids.new"
set -e
paste "${TEMP_DIR}/chromosome.old" "${TEMP_DIR}/chromosome.new" > "${MAPFILE}"
paste "${TEMP_DIR}/plasmids.old" "${TEMP_DIR}/plasmids.new" >> "${MAPFILE}"

rm -r "${TEMP_DIR}"