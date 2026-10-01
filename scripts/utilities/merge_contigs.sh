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

# Step 7 - merge graphs
echo "Step 7: Merge graphs"
## contig_N / bare Plassembler N -> final name, e.g. chromosome, plasmid_1
FLYE_GFA="$(dirname "${FLYE}")/assembly_graph.gfa"
PLASSEMBLER_GFA="$(dirname "${PLASSEMBLER}")/plassembler_plasmids.gfa"
# edge_N -> contig_N, read straight off the Flye graph's own P lines
awk 'BEGIN{FS="\t"; OFS="\t"} $1=="P"{
    n = split($3, edges, ",");
    for (i=1; i<=n; i++) {
        edge = edges[i];
        sub(/[+-]$/, "", edge);
        print edge, $2;
    }
}' "${FLYE_GFA}" > "${TEMP_DIR}/edge_to_contig.tsv"

rename_awk='
    BEGIN {
        while ((getline line < idmap) > 0) { split(line, a, "\t"); final[a[1]] = a[2] }
        while ((getline line < edgemap) > 0) { split(line, a, "\t"); contig[a[1]] = a[2] }
    }
    function known(id,    c) {
        if (id in final) return 1;
        c = (id in contig) ? contig[id] : id;
        return (c in final);
    }
    function rename(id,    c) {
        if (id in final) return final[id];
        c = (id in contig) ? contig[id] : id;
        return final[c];
    }
    $1=="S" {
        if (!known($2)) next;
        $2=rename($2); print; next
    }
    $1=="L" {
        if (!known($2) || !known($4)) next;
        $2=rename($2); $4=rename($4); print; next
    }
    $1=="P" {
        n=split($3, arr, ",");
        out="";
        drop=0;
        for (i=1; i<=n; i++) {
            orient=substr(arr[i], length(arr[i]), 1);
            id=substr(arr[i], 1, length(arr[i])-1);
            if (!known(id)) { drop=1; break }
            out=out (i>1?",":"") rename(id) orient;
        }
        if (drop) next;
        $2=rename($2);
        $3=out;
        print;
        next
    }
    $1=="H" { next }
    { print }
'

{
    echo -e "H\tVN:Z:1.0"
    awk -F'\t' -v OFS='\t' -v idmap="${OUTPUT}.id_map.txt" -v edgemap="${TEMP_DIR}/edge_to_contig.tsv" "${rename_awk}" "${FLYE_GFA}"
    awk -F'\t' -v OFS='\t' -v idmap="${OUTPUT}.id_map.txt" -v edgemap=/dev/null "${rename_awk}" "${PLASSEMBLER_GFA}"
} > "${OUTPUT}.gfa"


rm -r "${TEMP_DIR}"