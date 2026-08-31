#!/bin/bash

set -euo pipefail

# Runs the full single-sample workflow used across this training, start to finish,
# for one FASTQ file: QC -> filter -> species ID -> assembly -> plasmid recovery ->
# assembly QC -> polish -> assembly QC again -> AMR gene detection.
#
# Usage: run_sample_pipeline.sh <sample.fastq.gz>

# ─── Reference data paths ─────────────────────────────────────────────────────
K2DB=/home/tdev2/data/ref/kalamari                          # Kraken2 Kalamari database
PLASSEMBLER_DB=/home/tdev2/data/ref/plasmid_db_plassembler   # Plassembler plasmid database
BUSCO_DB=/home/tdev2/data/ref/busco/bacteria_odb12.2         # BUSCO lineage dataset (offline)
AMRFINDER_DB=/home/tdev2/data/ref/amrfinderplus_db/2026-05-15.1_4.2.7

# ─── Thread count ─────────────────────────────────────────────────────────────
THREADS=4

# ─── Helpers ──────────────────────────────────────────────────────────────────
log() {
    printf '\n[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

# ─── Argument parsing ─────────────────────────────────────────────────────────
[[ $# -eq 1 ]] || die "Usage: $(basename "$0") <sample.fastq.gz>"

input_fastq=$(realpath "$1")
[[ -f "${input_fastq}" ]] || die "FASTQ not found: ${input_fastq}"

sample_id=$(basename "${input_fastq}" .fastq.gz)
sample_id=$(basename "${sample_id}" .fastq)

log "Sample ID   : ${sample_id}"
log "Input FASTQ : ${input_fastq}"

# ─── Step 1 · QC on raw reads ──────────────────────────────────────────────────
log "Step 1: QC on raw reads (FastQC + NanoPlot + MultiQC)"

mkdir -p "qc_before/fastqc/${sample_id}" "qc_before/nanoplot/${sample_id}"

fastqc \
    -f fastq \
    -o "qc_before/fastqc/${sample_id}" \
    "${input_fastq}"

NanoPlot \
    --fastq "${input_fastq}" \
    -p "${sample_id}_" \
    --loglength \
    --N50 \
    -o "qc_before/nanoplot/${sample_id}/"

mkdir -p qc_before/multiqc
multiqc \
    -o qc_before/multiqc \
    -f \
    --fullnames \
    qc_before

# ─── Step 2 · Filter reads ─────────────────────────────────────────────────────
log "Step 2: Filter reads with Filtlong"

mkdir -p filtered
filtered_fastq="filtered/${sample_id}.filtered.fastq.gz"

filtlong \
    --min_length 1kb \
    --keep_percent 50 \
    "${input_fastq}" \
    | gzip > "${filtered_fastq}"

log "Step 2b: QC on filtered reads (FastQC + NanoPlot + MultiQC)"

mkdir -p "qc_after/fastqc/${sample_id}" "qc_after/nanoplot/${sample_id}"

fastqc \
    -f fastq \
    -o "qc_after/fastqc/${sample_id}" \
    "${filtered_fastq}"

NanoPlot \
    --fastq "${filtered_fastq}" \
    -p "${sample_id}_" \
    --loglength \
    --N50 \
    -o "qc_after/nanoplot/${sample_id}/"

mkdir -p qc_after/multiqc
multiqc \
    -o qc_after/multiqc \
    -f \
    --fullnames \
    qc_after

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
    --nano-hq "${filtered_fastq}" \
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

mkdir -p quast/draft busco bandage

quast \
    "${draft_assembly}" \
    --labels "${sample_id}.draft" \
    --output-dir quast/draft \
    --threads "${THREADS}"

busco \
    --in "${draft_assembly}" \
    --lineage_dataset "${BUSCO_DB}" \
    --out busco/draft \
    --mode genome \
    --offline \
    --cpu "${THREADS}"

bandage-exec Bandage image \
    flye/assembly_graph.gfa \
    "bandage/${sample_id}.flye_assembly_graph.svg"

if [[ -s "${plassembler_plasmids[0]}" ]]; then
    bandage-exec Bandage image \
        plassembler/*_plasmids.gfa \
        "bandage/${sample_id}.plassembler_plasmids_graph.svg"
fi

# ─── Step 7 · Polish the assembly (Medaka) ─────────────────────────────────────
log "Step 7: Polish assembly with Medaka"

medaka_consensus \
    -i "${filtered_fastq}" \
    -d "${draft_assembly}" \
    -o medaka \
    -t "${THREADS}" \
    --bacteria

polished_assembly="medaka/consensus.fasta"
[[ -s "${polished_assembly}" ]] || die "Polished assembly missing: ${polished_assembly}"

log "Step 7b: Assembly QC (QUAST + BUSCO) on the polished assembly"

mkdir -p busco/polished

quast \
    "${polished_assembly}" \
    --labels "${sample_id}.polished" \
    --output-dir quast/polished \
    --threads "${THREADS}"

busco \
    --in "${polished_assembly}" \
    --lineage_dataset "${BUSCO_DB}" \
    --out busco/polished \
    --mode genome \
    --offline \
    --cpu "${THREADS}"

# ─── Step 8 · AMR gene detection (AMRFinderPlus) ───────────────────────────────
log "Step 8: AMR gene detection with AMRFinderPlus"

mkdir -p amrfinder

amrfinder \
    -n "${polished_assembly}" \
    -d "${AMRFINDER_DB}" \
    --threads "${THREADS}" \
    > "amrfinder/${sample_id}.amrfinder.tsv"

# contig_diagram.py is a separate, manual visualisation step (see 12.annotation.md)
# that combines the Flye contig info with the AMRFinderPlus summary above; its
# command-line interface isn't documented in this training's materials, so it
# isn't automated here.

# ─── Done ───────────────────────────────────────────────────────────────────
log "Pipeline completed for ${sample_id}"
log "Draft assembly    : ${draft_assembly}"
log "Polished assembly : ${polished_assembly}"
log "AMRFinderPlus      : amrfinder/${sample_id}.amrfinder.tsv"
