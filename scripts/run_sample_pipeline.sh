#!/bin/bash

set -euo pipefail

# ─── Reference data paths ─────────────────────────────────────────────────────
K2DB=/home/tdev2/data/ref/kalamari                          # Kraken2 Kalamari database
PLASSEMBLER_DB=/home/tdev2/data/ref/plassembler/plasmid_db_plassembler             # Plassembler plasmid database
BUSCO_DB=/home/tdev2/data/ref/busco/bacteria_odb12.2        # BUSCO lineage dataset (offline)
AMRFINDER_DB=/home/tdev2/data/ref/amrfinderplus_db/amrfinderplus_V3.11_2022-12-19.1
BAKTA_DB=/home/tdev2/data/ref/bakta_database/7669534

# ─── Thread count ─────────────────────────────────────────────────────────────
THREADS=2

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
# --skip-fastqc --skip-trim --skip-species-id --skip-flye --skip-plassembler --skip-consensus --skip-polish --skip-assembly-qc --skip-amrfinder --skip-bakta
# bash run_sample_pipeline.sh --output-dir /home/tdev2/analysis --skip-fastqc --skip-trim --skip-species-id --skip-flye --skip-plassembler --skip-consensus --skip-bakta /home/tdev2/data/SAMEA5226451_A.baumannii/ERR8282753.fastq.gz > /home/tdev2/analysis/logs/ERR8282753/log04 2> /home/tdev2/analysis/logs/ERR8282753/err04
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

run_assembly_qc_stage() {
    local qc_id="$1"
    local assembly_path="$2"
    local bandage_graph="${3:-}"

    [[ ${skip_assembly_qc} -eq 0 ]] || return 0

    if [[ ! -s "${assembly_path}" ]]; then
        log "QC (${qc_id}) skipped (missing assembly: ${assembly_path})"
        return 0
    fi

    log "QC (${qc_id}): QUAST"
    mkdir -p "qc/assembly/quast/${qc_id}"
    quast \
        "${assembly_path}" \
        --labels "${qc_id}" \
        --output-dir "qc/assembly/quast/${qc_id}" \
        --threads "${THREADS}"

    log "QC (${qc_id}): BUSCO"
    mkdir -p "qc/assembly/busco/${qc_id}"

    # Bind /cvmfs so Singularity can find the database
    export SINGULARITY_COMMAND_OPTS="-B /cvmfs"

    busco \
        --in "${assembly_path}" \
        --out "qc/assembly/busco/${qc_id}" \
        --mode genome \
        --lineage_dataset "${BUSCO_DB}" \
        --force \
        --offline \
        --cpu "${THREADS}"

    if [[ -n "${bandage_graph}" && -s "${bandage_graph}" ]]; then
        log "QC (${qc_id}): Bandage"
        mkdir -p "qc/assembly/bandage/${qc_id}"
        bandage-exec Bandage image \
            "${bandage_graph}" \
            "qc/assembly/bandage/${qc_id}/${qc_id}.assembly.svg"
    else
        log "QC (${qc_id}): Bandage skipped (no assembly graph found)"
    fi
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
    log "Step 1: FastQC (pre-trim)"

    mkdir -p qc/fastqc/${sample_id}/pre_trim

    fastqc \
        -f fastq \
        -o qc/fastqc/${sample_id}/pre_trim \
        ${input_fastq}
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

if [[ ${skip_fastqc} -eq 0 ]]; then
    log "Step 2b: FastQC (post-trim)"

    mkdir -p qc/fastqc/${sample_id}/post_trim

    fastqc \
        -f fastq \
        -o qc/fastqc/${sample_id}/post_trim \
        ${current_fastq}
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

    run_assembly_qc_stage \
        "${current_id}.flye" \
        "${flye_dir}/assembly.fasta" \
        "${flye_dir}/assembly_graph.gfa"
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

    plassembler_qc_assembly=""
    for candidate in \
        "assembly/${current_id}/plassembler/plassembler_chromosome.fasta" \
        "assembly/${current_id}/plassembler/assembly.fasta" \
        "assembly/${current_id}/plassembler/plassembler_plasmids.fasta"; do
        if [[ -s "${candidate}" ]]; then
            plassembler_qc_assembly="${candidate}"
            break
        fi
    done

    run_assembly_qc_stage "${current_id}.plassembler" "${plassembler_qc_assembly}"
fi

# ─── Step 5 · Autocycler subsampling + consensus ──────────────────────────────
assembly_for_polish=""
CONSENSUSDIR=assembly/${current_id}/consensus
AUTOCYCLERDIR=${CONSENSUSDIR}/autocycler
existing_consensus_assembly="${AUTOCYCLERDIR}/consensus_assembly.fasta"

if [[ ${skip_consensus} -eq 0 ]]; then
    log "Step 5a: Autocycler subsampling"
    

    SUBSAMPLE_DIR=subsample/${current_id}
    mkdir -p ${SUBSAMPLE_DIR}

    autocycler subsample \
        --reads ${current_fastq} \
        --out_dir ${SUBSAMPLE_DIR} \
        --count 4 \
        --genome_size 5000000
        #  \
        # --min_read_depth 10

    # Prefix subsampled FASTQs with the sample ID
    for f in ${SUBSAMPLE_DIR}/*.fastq; do
        mv "${f}" "${SUBSAMPLE_DIR}/${current_id}.$(basename "${f}")"
    done

    log "Step 5b: Autocycler consensus"

    ASSEMBLYDIR=${CONSENSUSDIR}/assemblies
    mkdir -p ${ASSEMBLYDIR}

    for ASSEMBLER in flye plassembler; do        
        if [[ "${ASSEMBLER}" == "flye" ]]; then
            ASSEMBLY="assembly/${current_id}/flye/assembly.fasta"
        elif [[ "${ASSEMBLER}" == "plassembler" ]]; then
            ASSEMBLY="assembly/${current_id}/plassembler/plassembler_plasmids.fasta"
        fi
        [[ -s "${ASSEMBLY}" ]] || continue
        ln -sf "${PWD}/${ASSEMBLY}" "${ASSEMBLYDIR}/${current_id}.${ASSEMBLER}.fasta"
    done

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

    assembly_for_polish="${AUTOCYCLERDIR}/consensus_assembly.fasta"

    if [[ -s "${assembly_for_polish}" ]]; then
        log "Using consensus assembly: ${assembly_for_polish}"
        run_assembly_qc_stage \
            "${current_id}.consensus" \
            "${assembly_for_polish}" \
            "${AUTOCYCLERDIR}/consensus_assembly.gfa"
    else
            die "Consensus assembly ${assembly_for_polish} not found"
    fi
fi

if [[ -z "${assembly_for_polish}" && -s "${existing_consensus_assembly}" ]]; then
    assembly_for_polish="${existing_consensus_assembly}"
    log "Using existing consensus assembly: ${assembly_for_polish}"
fi

if [[ -z "${assembly_for_polish}" ]]; then
    flye_assembly="${flye_dir}/assembly.fasta"
    plassembler_assembly=""

    for candidate in \
        "assembly/${current_id}/plassembler/plassembler_chromosome.fasta" \
        "assembly/${current_id}/plassembler/assembly.fasta"; do
        if [[ -s "${candidate}" ]]; then
            plassembler_assembly="${candidate}"
            break
        fi
    done

    if [[ -z "${plassembler_assembly}" && -d "assembly/${current_id}/plassembler" ]]; then
        restore_nullglob=0
        shopt -q nullglob && restore_nullglob=1
        shopt -s nullglob

        for candidate in \
            assembly/${current_id}/plassembler/*.fasta \
            assembly/${current_id}/plassembler/*/*.fasta; do
            if [[ -s "${candidate}" ]]; then
                plassembler_assembly="${candidate}"
                break
            fi
        done

        if [[ ${restore_nullglob} -eq 0 ]]; then
            shopt -u nullglob
        fi
    fi

    if [[ -s "${flye_assembly}" ]]; then
        assembly_for_polish="${flye_assembly}"
        log "Using Flye assembly: ${assembly_for_polish}"
    elif [[ -n "${plassembler_assembly}" && -s "${plassembler_assembly}" ]]; then
        assembly_for_polish="${plassembler_assembly}"
        log "Using Plassembler assembly: ${assembly_for_polish}"
    else
        die "No assembly available for polishing/QC"
    fi
fi

# ─── Step 6 · Medaka polishing ────────────────────────────────────────────────
final_assembly="${assembly_for_polish}"
medaka_dir=assembly/${current_id}/medaka
existing_polished_assembly="${medaka_dir}/consensus.fasta"

if [[ ${skip_polish} -eq 0 ]]; then
    log "Step 6: Medaka polishing"

    medaka_consensus \
        -i ${current_fastq} \
        -d ${assembly_for_polish} \
        -o ${medaka_dir} \
        -t ${THREADS}

    final_assembly="${medaka_dir}/consensus.fasta"
    if [[ -s "${final_assembly}" ]]; then
        log "Using polished assembly: ${final_assembly}"
        run_assembly_qc_stage "${current_id}.medaka" "${final_assembly}"
    else
        die "Medaka assembly not found under ${medaka_dir}"
    fi
elif [[ -s "${existing_polished_assembly}" ]]; then
    final_assembly="${existing_polished_assembly}"
    log "Using existing polished assembly: ${final_assembly}"
fi

# ─── Step 7 · Assembly QC ─────────────────────────────────────────────────────
if [[ ${skip_assembly_qc} -eq 0 ]]; then
    qc_id="${current_id}.final"

    bandage_graph=""
    if [[ -s "${flye_dir}/assembly_graph.gfa" ]]; then
        bandage_graph="${flye_dir}/assembly_graph.gfa"
    elif [[ -d "${AUTOCYCLERDIR}" ]]; then
        restore_nullglob=0
        shopt -q nullglob && restore_nullglob=1
        shopt -s nullglob

        for candidate in \
            ${AUTOCYCLERDIR}/consensus_assembly.gfa \
            ${AUTOCYCLERDIR}/clustering/qc_pass/cluster_*/5_final.gfa \
            ${AUTOCYCLERDIR}/*.gfa; do
            if [[ -s "${candidate}" ]]; then
                bandage_graph="${candidate}"
                break
            fi
        done

        if [[ ${restore_nullglob} -eq 0 ]]; then
            shopt -u nullglob
        fi
    fi

    run_assembly_qc_stage "${qc_id}" "${final_assembly}" "${bandage_graph}"
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

    # Bind /cvmfs so Singularity can find the database
    export SINGULARITY_COMMAND_OPTS="-B /cvmfs"

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

    # Bind /cvmfs so Singularity can find the database
    export SINGULARITY_COMMAND_OPTS="-B /cvmfs"

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

# ─── Step 12 · Compress uncompressed FASTQ files ─────────────────────────────
log "Step 12: Compressing uncompressed FASTQ files"

mapfile -d '' fastq_files < <(find . -type f -name "*.fastq" -print0)

if [[ ${#fastq_files[@]} -eq 0 ]]; then
    log "No uncompressed FASTQ files found"
else
    for fastq in "${fastq_files[@]}"; do
        gzip -f "${fastq}"
    done
    log "Compressed ${#fastq_files[@]} FASTQ file(s)"
fi

# ─── Done ─────────────────────────────────────────────────────────────────────
log "Pipeline completed for ${current_id}"
log "Final assembly: ${final_assembly}"
