#!/bin/bash

set -euo pipefail

# ─── Reference data paths ─────────────────────────────────────────────────────
K2DB=/home/tdev2/data/ref/kalamari                          # Kraken2 Kalamari database
PLASSEMBLER_DB=/home/tdev2/data/ref/plassembler/plasmid_db_plassembler             # Plassembler plasmid database
BUSCO_DB=/home/tdev2/data/ref/busco/bacteria_odb12.2        # BUSCO lineage dataset (offline)
AMRFINDER_DB=/home/tdev2/data/ref/amrfinderplus_db/amrfinderplus_V3.12_2024-05-02.2
BAKTA_DB=/home/tdev2/data/ref/bakta_database/7669534

# ─── Thread count ─────────────────────────────────────────────────────────────
THREADS=1

# ─── Skip flags ───────────────────────────────────────────────────────────────
skip_fastqc=0
skip_trim=0
skip_species_id=0
skip_flye=0
skip_plassembler=0
skip_consensus=0
skip_polish=0
skip_assembly_qc=0
skip_amrfinder=0
skip_bakta=0
skip_multiqc=0
chromosome_contig=""
OUTPUT_DIR="$(pwd)"

# ─── Helpers ──────────────────────────────────────────────────────────────────
usage() {
    echo "Usage: run_sample_pipeline.sh [options] <fastq|fastq.gz>"
    echo ""
    echo "Run the full analysis pipeline for one sample: QC -> trim -> species ID ->"
    echo "assembly -> consensus -> polish -> assembly QC -> annotation -> report."
    echo ""
    echo "Options:"
    echo "  --output-dir DIR      Output working directory (default: current directory)"
    echo "  --skip-fastqc         Skip FastQC"
    echo "  --skip-trim           Skip read trimming with fastplong"
    echo "  --skip-species-id     Skip Kraken2 species identification"
    echo "  --skip-flye           Skip Flye assembly"
    echo "  --skip-plassembler    Skip Plassembler assembly"
    echo "  --skip-consensus      Skip Autocycler consensus assembly"
    echo "  --skip-polish         Skip Medaka polishing"
    echo "  --skip-assembly-qc    Skip QUAST, BUSCO, and Bandage QC"
    echo "  --skip-amrfinder      Skip AMRFinderPlus"
    echo "  --skip-bakta          Skip Bakta annotation"
    echo "  --skip-multiqc        Skip final MultiQC summary"
    echo "  --chromosome-contig NAME"
    echo "                        Extract plasmid contigs by excluding this contig header"
    echo "  -h, --help            Show this help"
}

log() {
    printf '\n[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

find_first_file() {
    local search_dir="$1"
    shift
    local pattern
    for pattern in "$@"; do
        while IFS= read -r candidate; do
            if [[ -n "${candidate}" && -s "${candidate}" ]]; then
                printf '%s\n' "${candidate}"
                return 0
            fi
        done < <(find "${search_dir}" -maxdepth 3 -type f -name "${pattern}" | sort)
    done
    return 1
}

# ─── Argument parsing ─────────────────────────────────────────────────────────
args=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir)
            [[ $# -ge 2 ]] || die "--output-dir requires a value"
            OUTPUT_DIR="$2"; shift ;;
        --skip-fastqc)       skip_fastqc=1 ;;
        --skip-trim)         skip_trim=1 ;;
        --skip-species-id)   skip_species_id=1 ;;
        --skip-flye)         skip_flye=1 ;;
        --skip-plassembler)  skip_plassembler=1 ;;
        --skip-consensus)    skip_consensus=1 ;;
        --skip-polish)       skip_polish=1 ;;
        --skip-assembly-qc)  skip_assembly_qc=1 ;;
        --skip-amrfinder)    skip_amrfinder=1 ;;
        --skip-bakta)        skip_bakta=1 ;;
        --skip-multiqc)      skip_multiqc=1 ;;
        --chromosome-contig)
            [[ $# -ge 2 ]] || die "--chromosome-contig requires a value"
            chromosome_contig="$2"; shift ;;
        -h|--help) usage; exit 0 ;;
        --*) die "Unknown option: $1" ;;
        *) args+=("$1") ;;
    esac
    shift
done

[[ ${#args[@]} -eq 1 ]] || { usage >&2; exit 1; }

input_fastq=$(realpath "${args[0]}")
[[ -f "${input_fastq}" ]] || die "FASTQ not found: ${input_fastq}"

mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")
cd "${OUTPUT_DIR}"

# Strip .fastq.gz or .fastq to get a bare sample ID
sample_id=$(basename "${input_fastq}" .fastq.gz)
sample_id=$(basename "${sample_id}" .fastq)
trimmed_fastq="trim/${sample_id}.trimmed.fastq"

current_fastq="${input_fastq}"
current_id="${sample_id}"

log "Input FASTQ  : ${current_fastq}"
log "Sample ID    : ${current_id}"
log "Output dir   : ${OUTPUT_DIR}"

# ─── Step 1 · FastQC ──────────────────────────────────────────────────────────
if [[ ${skip_fastqc} -eq 0 ]]; then
    log "Step 1: FastQC"

    mkdir -p qc/fastqc/${current_id}

    fastqc \
        -f fastq \
        -o qc/fastqc/${current_id} \
        ${current_fastq}
fi

# ─── Step 2 · Read trimming (fastplong) ───────────────────────────────────────
if [[ ${skip_trim} -eq 0 ]]; then
    log "Step 2: fastplong trimming"
    

    mkdir -p trim
    mkdir -p qc/fastplong

    fastplong \
        -i ${current_fastq} \
        --disable_quality_filtering \
        --disable_length_filtering \
        -o ${trimmed_fastq} \
        -h qc/fastplong/${current_id}.fastplong_report.html \
        -j qc/fastplong/${current_id}.fastplong_report.json

    [[ -s "${trimmed_fastq}" ]] || die "Trimmed FASTQ missing: ${trimmed_fastq}"
fi

# Always prefer trimmed reads when available, even if trimming is skipped.
if [[ -s "${trimmed_fastq}" ]]; then
    current_fastq="${trimmed_fastq}"
    current_id=$(basename "${current_fastq}" .fastq)
    log "Using trimmed FASTQ: ${current_fastq}"
else
    current_fastq="${input_fastq}"
    current_id="${sample_id}"
fi

# ─── Step 3 · Species identification (Kraken2) ────────────────────────────────
if [[ ${skip_species_id} -eq 0 ]]; then
    log "Step 3: Kraken2 species identification"

    mkdir -p qc/species_id/report

    # Bind /cvmfs so Singularity can find the database
    export SINGULARITY_COMMAND_OPTS="-B /cvmfs"

    kraken2 \
        ${current_fastq} \
        --db ${K2DB} \
        --report qc/species_id/report/${current_id}.k2report \
        --report-minimizer-data \
        --minimum-hit-groups 3 \
        --threads ${THREADS} \
        --output qc/species_id/${current_id}_k2_out.txt
fi

# ─── Step 4a · Flye assembly ──────────────────────────────────────────────────
flye_dir="assembly/${current_id}/flye"

if [[ ${skip_flye} -eq 0 ]]; then
    log "Step 4a: Flye assembly"

    mkdir -p ${flye_dir}

    flye \
        --nano-hq ${current_fastq} \
        --threads ${THREADS} \
        --out-dir ${flye_dir}

    [[ -s "${flye_dir}/assembly.fasta" ]] || \
        die "Flye assembly missing: ${flye_dir}/assembly.fasta"
fi

# ─── Step 4b · Plassembler ────────────────────────────────────────────────────
if [[ ${skip_plassembler} -eq 0 ]]; then
    log "Step 4b: Plassembler"

    mkdir -p assembly/${current_id}/plassembler

    FLYE_ARGS=""
    if [[ -d "${flye_dir}" && -f "${flye_dir}/assembly.fasta" && -f "${flye_dir}/assembly_info.txt" ]]; then
        FLYE_ARGS="--flye_assembly ${flye_dir}/assembly.fasta --flye_info ${flye_dir}/assembly_info.txt"
    fi

    plassembler-exec plassembler long \
        -d ${PLASSEMBLER_DB} \
        -l ${current_fastq} \
        ${FLYE_ARGS} \
        -t ${THREADS} \
        -f \
        -o assembly/${current_id}/plassembler
fi

# ─── Step 5 · Autocycler subsampling + consensus ──────────────────────────────
assembly_for_polish=""

if [[ ${skip_consensus} -eq 0 ]]; then
    log "Step 5a: Autocycler subsampling"
    

    SUBSAMPLE_DIR=subsample/${current_id}
    mkdir -p ${SUBSAMPLE_DIR}

    autocycler subsample \
        --reads ${current_fastq} \
        --out_dir ${SUBSAMPLE_DIR} \
        --count 4 \
        --genome_size 5000000

    # Prefix subsampled FASTQs with the sample ID
    for f in ${SUBSAMPLE_DIR}/*.fastq; do
        mv "${f}" "${SUBSAMPLE_DIR}/${current_id}.$(basename "${f}")"
    done

    log "Step 5b: Autocycler consensus"

    CONSENSUSDIR=assembly/${current_id}/consensus
    ASSEMBLYDIR=${CONSENSUSDIR}/assemblies
    mkdir -p ${ASSEMBLYDIR}

    for ASSEMBLER in flye plassembler; do
        ASSEMBLY=$(find assembly/${current_id}/${ASSEMBLER}/ -maxdepth 1 -type f -name "*.fasta" | head -n 1)
        [[ -s "${ASSEMBLY}" ]] || continue
        ln -sf "${PWD}/${ASSEMBLY}" "${ASSEMBLYDIR}/${current_id}.${ASSEMBLER}.fasta"
    done

    AUTOCYCLERDIR=${CONSENSUSDIR}/autocycler
    mkdir -p ${AUTOCYCLERDIR}

    autocycler compress \
        --assemblies_dir ${ASSEMBLYDIR} \
        --autocycler_dir ${AUTOCYCLERDIR} \
        --threads ${THREADS}

    autocycler cluster \
        --autocycler_dir ${AUTOCYCLERDIR} \
        --min_assemblies 1

    for d in ${AUTOCYCLERDIR}/clustering/qc_pass/cluster_*; do
        autocycler trim \
            --cluster_dir ${d} \
            --threads ${THREADS}

        autocycler resolve \
            --cluster_dir ${d}
    done

    autocycler combine \
        --autocycler_dir ${AUTOCYCLERDIR} \
        --in_gfas ${AUTOCYCLERDIR}/clustering/qc_pass/cluster_*/5_final.gfa

    autocycler table > ${CONSENSUSDIR}/autocycler_metrics.tsv
    autocycler table \
        --autocycler_dir ${AUTOCYCLERDIR} \
        --name consensus >> ${CONSENSUSDIR}/autocycler_metrics.tsv

    if assembly_for_polish=$(find_first_file "${CONSENSUSDIR}" "consensus*.fasta" "*.fasta"); then
        log "Using consensus assembly: ${assembly_for_polish}"
    else
        die "Consensus assembly not found under ${CONSENSUSDIR}"
    fi
fi

if [[ -z "${assembly_for_polish}" ]]; then
    assembly_for_polish="${flye_dir}/assembly.fasta"
    [[ -s "${assembly_for_polish}" ]] || die "No assembly available for polishing/QC"
fi

# ─── Step 6 · Medaka polishing ────────────────────────────────────────────────
final_assembly="${assembly_for_polish}"

if [[ ${skip_polish} -eq 0 ]]; then
    log "Step 6: Medaka polishing"
    

    medaka_dir=assembly/${current_id}/medaka
    mkdir -p ${medaka_dir}

    medaka_consensus \
        -i ${current_fastq} \
        -d ${assembly_for_polish} \
        -o ${medaka_dir} \
        -t ${THREADS}

    if final_assembly=$(find_first_file "${medaka_dir}" "consensus*.fasta" "*.fasta"); then
        log "Using polished assembly: ${final_assembly}"
    else
        die "Medaka assembly not found under ${medaka_dir}"
    fi
fi

# ─── Step 7 · Assembly QC ─────────────────────────────────────────────────────
if [[ ${skip_assembly_qc} -eq 0 ]]; then
    qc_id="${current_id}.final"

    log "Step 7a: QUAST"
    

    mkdir -p qc/assembly/quast/${qc_id}

    quast \
        ${final_assembly} \
        --labels ${qc_id} \
        --output-dir qc/assembly/quast/${qc_id} \
        --threads ${THREADS}

    log "Step 7b: BUSCO"
    

    mkdir -p qc/assembly/busco/${qc_id}

    busco \
        --in ${final_assembly} \
        --out qc/assembly/busco/${qc_id} \
        --mode genome \
        --lineage_dataset ${BUSCO_DB} \
        --force \
        --offline \
        --cpu ${THREADS}

    bandage_graph=""
    if [[ -s "${flye_dir}/assembly_graph.gfa" ]]; then
        bandage_graph="${flye_dir}/assembly_graph.gfa"
    elif [[ ${skip_consensus} -eq 0 ]]; then
        bandage_graph=$(find_first_file "assembly/${current_id}/consensus" "*.gfa" || true)
    fi

    if [[ -n "${bandage_graph}" ]]; then
        log "Step 7c: Bandage"
        

        mkdir -p qc/assembly/bandage/${qc_id}

        bandage_exec \
        Bandage image \
            ${bandage_graph} \
            qc/assembly/bandage/${qc_id}/${qc_id}.assembly.svg
    else
        log "Step 7c: Bandage skipped (no assembly graph found)"
    fi
fi

# ─── Step 8 · Extract plasmid contigs ─────────────────────────────────────────
if [[ -n "${chromosome_contig}" ]]; then
    log "Step 8: Extracting plasmid contigs (excluding ${chromosome_contig})"

    mkdir -p assembly/${current_id}/plasmids

    awk -v exclude="^>${chromosome_contig}" \
        '$0 ~ /^>/ { p = 1 } $0 ~ exclude { p = 0 } p == 1 { print $0 }' \
        ${final_assembly} \
        > assembly/${current_id}/plasmids/${current_id}.plasmids.fasta
fi

# ─── Step 9 · AMRFinderPlus ───────────────────────────────────────────────────
if [[ ${skip_amrfinder} -eq 0 ]]; then
    log "Step 9: AMRFinderPlus"
    

    mkdir -p annotate/${current_id}/amrfinder

    amrfinder \
        -n ${final_assembly} \
        -d ${AMRFINDER_DB} \
        --threads ${THREADS} \
        > annotate/${current_id}/amrfinder/${current_id}.tsv
fi

# ─── Step 10 · Bakta annotation ───────────────────────────────────────────────
if [[ ${skip_bakta} -eq 0 ]]; then
    log "Step 10: Bakta annotation"
    

    mkdir -p annotate/${current_id}

    bakta ${final_assembly} \
        --db ${BAKTA_DB} \
        --output annotate/${current_id}/ \
        --prefix ${current_id} \
        --force \
        --threads ${THREADS}
fi

# ─── Step 11 · MultiQC summary ────────────────────────────────────────────────
if [[ ${skip_multiqc} -eq 0 ]]; then
    log "Step 11: MultiQC"

    mkdir -p multiqc

    multiqc \
        -o multiqc \
        -f \
        --fullnames \
        qc
fi

# ─── Done ─────────────────────────────────────────────────────────────────────
log "Pipeline completed for ${current_id}"
log "Final assembly: ${final_assembly}"
