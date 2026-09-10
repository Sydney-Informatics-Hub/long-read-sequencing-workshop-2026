#!/bin/bash

set -euo pipefail

# Runs the full single-sample workflow used across this training, start to finish,
# for one FASTQ file: QC -> filter -> species ID -> assembly -> plasmid recovery ->
# assembly QC -> polish -> assembly QC again -> AMR gene detection.
#
# Usage: run_sample_pipeline.sh <sample.fastq.gz>

# ─── Reference data paths ─────────────────────────────────────────────────────
K2DB=/home/tdev3/data/ref/kalamari                          # Kraken2 Kalamari database
PLASSEMBLER_DB=/home/tdev3/data/ref/plasmid_db_plassembler   # Plassembler plasmid database
BUSCO_DB=/home/tdev3/data/ref/busco/bacteria_odb12.2         # BUSCO lineage dataset (offline)
AMRFINDER_DB=/home/tdev3/data/ref/amrfinderplus_db/2026-08-07.1
MEDAKA_IMAGE_PATH=/home/tdev3/sing_images/medaka_1.3.3--py38h130def0_0
MEDAKA_MODEL=r941_min_high_g360
AMR_CASSETTE_DIAGRAM_SCRIPT="$(dirname "$(realpath "$0")")/utilities/amr_cassette_diagram.py"
COMBINE_BANDAGE_GRAPHS_SCRIPT="$(dirname "$(realpath "$0")")/utilities/combine_bandage_graphs.sh"

# ─── Thread count ────────────────────────────────────────────────────────────────
THREADS=4

# ─── Helpers ──────────────────────────────────────────────────────────────────
log() {
    local msg
    msg=$(printf '\n[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*")
    printf '%s' "${msg}"
    printf '%s' "${msg}" >&2
}

die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

# ─── Argument parsing ─────────────────────────────────────────────────────────
usage="Usage: $(basename "$0") <sample.fastq.gz>"

[[ $# -eq 1 ]] || die "${usage}"

input_fastq=$(realpath "$1")
[[ -f "${input_fastq}" ]] || die "FASTQ not found: ${input_fastq}"

sample_id=$(basename "${input_fastq}" .fastq.gz)
sample_id=$(basename "${sample_id}" .fastq)

log "Sample ID   : ${sample_id}"
log "Input FASTQ : ${input_fastq}"

# ─── Step 1 · QC on raw reads ──────────────────────────────────────────────────
log "Step 1: QC on raw reads (FastQC + NanoPlot + MultiQC)"

mkdir -p "read_qc/${sample_id}/fastqc/raw" "read_qc/${sample_id}/nanoplot/raw"

fastqc \
    -f fastq \
    -o "read_qc/${sample_id}/fastqc/raw/" \
    "${input_fastq}"

NanoPlot \
    --fastq "${input_fastq}" \
    -p "${sample_id}_raw_" \
    --loglength \
    --N50 \
    -o "read_qc/${sample_id}/nanoplot/raw/"

# ─── Step 2 · Filter reads ─────────────────────────────────────────────────────
log "Step 2: Filter reads with Filtlong"

mkdir -p filtered
filtered_fastq="filtered/${sample_id}.filtered.fastq.gz"

filtlong \
    --min_length 1kb \
    --target_bases 150mb \
    "${input_fastq}" \
    | gzip > "${filtered_fastq}"

log "Step 2b: QC on filtered reads (FastQC + NanoPlot + MultiQC)"

mkdir -p "read_qc/${sample_id}/fastqc/filtered" "read_qc/${sample_id}/nanoplot/filtered"

fastqc \
    -f fastq \
    -o "read_qc/${sample_id}/fastqc/filtered" \
    "${filtered_fastq}"

NanoPlot \
    --fastq "${filtered_fastq}" \
    -p "${sample_id}_filtered_" \
    --loglength \
    --N50 \
    -o "read_qc/${sample_id}/nanoplot/filtered/"

multiqc \
    -o "read_qc" \
    -f \
    --fullnames \
    read_qc

# ─── Step 3 · Species identification (Kraken2) ─────────────────────────────────
log "Step 3: Species identification with Kraken2"

mkdir -p kraken2

# Bind /cvmfs so Singularity can find the database
export SINGULARITY_COMMAND_OPTS="-B /cvmfs"

kraken2 \
    --db "${K2DB}" \
    --report "kraken2/${sample_id}.k2report" \
    --output "kraken2/${sample_id}.k2_out.txt" \
    --threads "${THREADS}" \
    "${filtered_fastq}"

multiqc \
    -o kraken2 \
    -f \
    kraken2

# ─── Step 4 · De novo assembly (Flye) ──────────────────────────────────────────
log "Step 4: De novo assembly with Flye"

mkdir -p flye

flye \
    --nano-raw "${filtered_fastq}" \
    --out-dir flye \
    --threads "${THREADS}"

[[ -s flye/assembly.fasta ]] || die "Flye assembly missing: flye/assembly.fasta"

# ─── Step 5 · Plasmid recovery (Plassembler) ───────────────────────────────────
log "Step 5: Plasmid recovery with Plassembler"

plassembler-exec plassembler long \
    -l "${filtered_fastq}" \
    -d "${PLASSEMBLER_DB}" \
    --flye_assembly flye/assembly.fasta \
    --flye_info flye/assembly_info.txt \
    -t "${THREADS}" \
    -o plassembler

# Fold any additional plasmid contigs Plassembler found into a single draft assembly
draft_assembly="draft_assembly.fasta"
plassembler_plasmids=(plassembler/*_plasmids.fasta)

if [[ -s "${plassembler_plasmids[0]}" ]]; then
    cat flye/assembly.fasta "${plassembler_plasmids[0]}" > "${draft_assembly}"
else
    cp flye/assembly.fasta "${draft_assembly}"
fi

# ─── Step 6 · Assembly QC on the draft assembly ────────────────────────────────
log "Step 6: Assembly QC (QUAST + BUSCO + Bandage) on the draft assembly"

mkdir -p assembly_qc/quast/draft assembly_qc/busco assembly_qc/bandage

quast \
    "${draft_assembly}" \
    --labels "${sample_id}.draft" \
    --output-dir assembly_qc/quast/draft \
    --threads "${THREADS}"

busco \
    --in "${draft_assembly}" \
    --lineage_dataset "${BUSCO_DB}" \
    --out assembly_qc/busco/draft \
    --mode genome \
    --offline \
    --cpu "${THREADS}"

bandage-exec Bandage image \
    flye/assembly_graph.gfa \
    "assembly_qc/bandage/${sample_id}.flye_assembly_graph.svg"

if [[ -s "${plassembler_plasmids[0]}" ]]; then
    plassembler_gfa=(plassembler/*_plasmids.gfa)

    bandage-exec Bandage image \
        "${plassembler_gfa[0]}" \
        "assembly_qc/bandage/${sample_id}.plassembler_plasmids_graph.svg"

    # Combined Flye + Plassembler graph, contigs labelled by source program,
    # name, and length (see utilities/combine_bandage_graphs.sh for how). Best-effort:
    # a failure here shouldn't take down the rest of the pipeline over what
    # is just an extra visualisation on top of the two Bandage images above.
    "${COMBINE_BANDAGE_GRAPHS_SCRIPT}" \
        flye/assembly_graph.gfa \
        "${plassembler_gfa[0]}" \
        "assembly_qc/bandage/${sample_id}.combined_assembly_graph.svg" \
        || log "Warning: combined Bandage graph failed, continuing without it"
fi

# ─── Step 7 · Polish the assembly (Medaka) ─────────────────────────────────────
log "Step 7: Polish assembly with Medaka"

singularity exec "${MEDAKA_IMAGE_PATH}" medaka_consensus \
    -i "${filtered_fastq}" \
    -d "${draft_assembly}" \
    -m "${MEDAKA_MODEL}" \
    -o medaka \
    -t "${THREADS}" \
    -b 50

polished_assembly="medaka/consensus.fasta"
[[ -s "${polished_assembly}" ]] || die "Polished assembly missing: ${polished_assembly}"

log "Step 7b: Assembly QC (QUAST + BUSCO) on the polished assembly"

quast \
    "${polished_assembly}" \
    --labels "${sample_id}.polished" \
    --output-dir assembly_qc/quast/polished \
    --threads "${THREADS}"

busco \
    --in "${polished_assembly}" \
    --lineage_dataset "${BUSCO_DB}" \
    --out assembly_qc/busco/polished \
    --mode genome \
    --offline \
    --cpu "${THREADS}"

multiqc \
    -o assembly_qc \
    -f \
    --fullnames \
    assembly_qc

# ─── Step 8 · AMR gene detection (AMRFinderPlus) ───────────────────────────────
log "Step 8: AMR gene detection with AMRFinderPlus"

mkdir -p amrfinder

amrfinder \
    -n "${polished_assembly}" \
    -d "${AMRFINDER_DB}" \
    --plus \
    --threads "${THREADS}" \
    > "amrfinder/${sample_id}.amrfinder_plus.tsv"

# ─── Step 9 · AMR gene-cassette diagram (pyGenomeViz) ──────────────────────────
log "Step 9: AMR gene-cassette diagram"

pygenomeviz-exec python3 "${AMR_CASSETTE_DIAGRAM_SCRIPT}" \
    --amrfinder-tsv "amrfinder/${sample_id}.amrfinder_plus.tsv"

# ─── Done ───────────────────────────────────────────────────────────────────
log "Pipeline completed for ${sample_id}"
log "Draft assembly    : ${draft_assembly}"
log "Polished assembly : ${polished_assembly}"
log "AMRFinderPlus      : amrfinder/${sample_id}.amrfinder_plus.tsv"
log "AMR cassette plots : amrfinder/${sample_id}.amrfinder_plus.<contig>.png"
