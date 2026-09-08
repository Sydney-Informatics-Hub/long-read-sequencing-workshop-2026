#!/bin/bash

shopt -s expand_aliases

set -euo pipefail

# --- Defaults -----------------------------------------------------------
DB=""
FASTA=""
OUTPUT=""

# --- Flag parsing ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --db)
            DB="$2"
            shift
            shift
            ;;
        --fasta)
            FASTA="$2"
            shift
            shift
            ;;
        --output)
            OUTPUT="$2"
            shift
            shift
            ;;
        -h|--help)
            echo "Usage: $0 --fasta /path/to/assembly.fasta --db /path/to/plassembler/db/dir [--output /path/to/output.tsv]"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

if [ -z "${FASTA}" ]; then
    echo "ERROR: Must provide a FASTA file."
    exit 1
fi

if [ -z "${DB}" ] || [ ! -d "${DB}" ]; then
    echo "ERROR: Must provide a plassembler database."
    exit 1
fi

DB_MSH_FILE=$(find "${DB}" -type f -name "plsdb_*.msh" | head -n 1)
DB_TSV_FILE="${DB}/"$(basename "${DB_MSH_FILE}" .msh).tsv
if [ ! -f "${DB_MSH_FILE}" ] || [ ! -f "${DB_TSV_FILE}" ]; then
    echo "ERROR: Must provide a valid plassembler database."
    exit 1
fi

if [ -z "${OUTPUT}" ]; then
    OUTPUT="${FASTA}.plsdb_results.tsv"
fi

if [ -f "${OUTPUT}" ]; then
    echo "ERROR: Output file '${OUTPUT}' already exists."
    exit 1
fi

# Split FASTA by contig ID
seqkit split "${FASTA}" --by-id
SPLIT_CONTIG_DIR="${FASTA}.split"

# Run mash screen
rm -f "${OUTPUT}"
for CONTIG in ${SPLIT_CONTIG_DIR}/*.fasta; do
    rm -f "${OUTPUT}.tmp"{1,2}
    singularity exec /cvmfs/singularity.galaxyproject.org/all/plassembler:1.8.2--pyhdfd78af_0 \
    mash screen \
        -i 0.99 \
        -v 0.1 \
        ${DB_MSH_FILE} \
        ${CONTIG} >> "${OUTPUT}.tmp1"
    cut -f 5 "${OUTPUT}.tmp1" | while read PLASMID; do
        grep -P "\t${PLASMID}\t" ${DB_TSV_FILE} | cut -f 2-9 >> "${OUTPUT}.tmp2"
    done
    echo "===== "$(basename "${CONTIG}")" =====" >> "${OUTPUT}"
    echo -e "Identity\tShared Hashes\tMedian Multiplicity\tP-value\tNUCCORE_ACC\tDescription" >> "${OUTPUT}"
    cat "${OUTPUT}.tmp1" >> "${OUTPUT}"
    echo "====================" >> "${OUTPUT}"
    head -n 1 ${DB_TSV_FILE} | cut -f 2-9 >> "${OUTPUT}"
    cat "${OUTPUT}.tmp2" >> "${OUTPUT}"
done
