#!/bin/bash

set -euo pipefail

# ─── Reference data paths ─────────────────────────────────────────────────────
K2DB=/home/tdev2/data/ref/kalamari                          # Kraken2 Kalamari database
PLASSEMBLER_DB=/home/tdev2/data/ref/plasmid_db_plassembler             # Plassembler plasmid database
BUSCO_DB=/home/tdev2/data/ref/busco/bacteria_odb12.2        # BUSCO lineage dataset (offline)
AMRFINDER_DB=/home/tdev2/data/ref/amrfinderplus_db/2026-05-15.1_4.2.7

# ─── Thread count ─────────────────────────────────────────────────────────────
THREADS=2

# ─── Pipeline steps (execution order; also the vocabulary for --start-from) ───
STEP_ORDER=(fastqc trim species-id flye plassembler consensus polish assembly-qc amrfinder multiqc)

declare -A explicit_skip=()
start_from=""
chromosome_contig=""
OUTPUT_DIR="$(pwd)"

# Example: resume polishing from existing Flye/Plassembler outputs, treating consensus as excluded:
# bash run_sample_pipeline.sh --output-dir /home/tdev2/analysis --start-from polish --skip-consensus /home/tdev2/data/SAMEA5226451_A.baumannii/ERR8282753.fastq.gz > log 2> err

# ─── Helpers ──────────────────────────────────────────────────────────────────
usage() {
    echo "Usage: run_sample_pipeline.sh [options] <fastq|fastq.gz>"
    echo ""
    echo "Run the full analysis pipeline for one sample: QC -> trim -> species ID ->"
    echo "assembly -> consensus -> polish -> assembly QC -> AMRFinderPlus -> report."
    echo ""
    echo "Options:"
    echo "  --output-dir DIR      Output working directory (default: current directory)"
    echo "  --start-from STEP     Resume from STEP, treating all earlier steps as already"
    echo "                        run in a previous invocation: their outputs must already"
    echo "                        exist on disk and will be reused (the pipeline dies with"
    echo "                        a clear error if they don't). STEP is one of:"
    echo "                        ${STEP_ORDER[*]}"
    echo "  --skip-fastqc         Exclude FastQC entirely"
    echo "  --skip-trim           Exclude read trimming with fastplong"
    echo "  --skip-species-id     Exclude Kraken2 species identification"
    echo "  --skip-flye           Exclude Flye assembly"
    echo "  --skip-plassembler    Exclude Plassembler assembly"
    echo "  --skip-consensus      Exclude Autocycler consensus; Flye and Plassembler"
    echo "                        assemblies are polished with Medaka separately instead,"
    echo "                        each followed by its own AMRFinderPlus run"
    echo "  --skip-polish         Exclude Medaka polishing"
    echo "  --skip-assembly-qc    Exclude QUAST, BUSCO, and Bandage QC"
    echo "  --skip-amrfinder      Exclude AMRFinderPlus"
    echo "  --skip-multiqc        Exclude final MultiQC summary"
    echo "  --chromosome-contig NAME"
    echo "                        Extract plasmid contigs by excluding this contig header"
    echo "  -h, --help            Show this help"
    echo ""
    echo "--skip-X means the step is left out of the pipeline entirely (its absence is"
    echo "expected downstream). --start-from means everything before STEP was already"
    echo "run and its output is required; combine the two to say e.g. \"resume from"
    echo "polishing, and I never want consensus at all\"."
}

log() {
    printf '\n[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

step_index() {
    local target="$1" i
    for i in "${!STEP_ORDER[@]}"; do
        if [[ "${STEP_ORDER[$i]}" == "${target}" ]]; then
            printf '%s' "$i"
            return 0
        fi
    done
    return 1
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

find_plassembler_assembly() {
    local id="$1"
    local found=""

    for candidate in \
        "assembly/${id}/plassembler/plassembler_chromosome.fasta" \
        "assembly/${id}/plassembler/assembly.fasta" \
        "assembly/${id}/plassembler/plassembler_plasmids.fasta"; do
        if [[ -s "${candidate}" ]]; then
            found="${candidate}"
            break
        fi
    done

    if [[ -z "${found}" && -d "assembly/${id}/plassembler" ]]; then
        local restore_nullglob=0
        shopt -q nullglob && restore_nullglob=1
        shopt -s nullglob

        for candidate in \
            assembly/${id}/plassembler/*.fasta \
            assembly/${id}/plassembler/*/*.fasta; do
            if [[ -s "${candidate}" ]]; then
                found="${candidate}"
                break
            fi
        done

        if [[ ${restore_nullglob} -eq 0 ]]; then
            shopt -u nullglob
        fi
    fi

    printf '%s' "${found}"
}

find_fastq() {
    local base="$1"

    if [[ -s "${base}.fastq" ]]; then
        printf '%s' "${base}.fastq"
    elif [[ -s "${base}.fastq.gz" ]]; then
        printf '%s' "${base}.fastq.gz"
    fi
}

# ─── Argument parsing ─────────────────────────────────────────────────────────
args=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir)
            [[ $# -ge 2 ]] || die "--output-dir requires a value"
            OUTPUT_DIR="$2"; shift ;;
        --start-from)
            [[ $# -ge 2 ]] || die "--start-from requires a value"
            start_from="$2"; shift ;;
        --skip-fastqc)       explicit_skip[fastqc]=1 ;;
        --skip-trim)         explicit_skip[trim]=1 ;;
        --skip-species-id)   explicit_skip[species-id]=1 ;;
        --skip-flye)         explicit_skip[flye]=1 ;;
        --skip-plassembler)  explicit_skip[plassembler]=1 ;;
        --skip-consensus)    explicit_skip[consensus]=1 ;;
        --skip-polish)       explicit_skip[polish]=1 ;;
        --skip-assembly-qc)  explicit_skip[assembly-qc]=1 ;;
        --skip-amrfinder)    explicit_skip[amrfinder]=1 ;;
        --skip-multiqc)      explicit_skip[multiqc]=1 ;;
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

# ─── Resolve --start-from / --skip-* into per-step skip & require flags ───────
start_from_idx=-1
if [[ -n "${start_from}" ]]; then
    start_from_idx=$(step_index "${start_from}") || \
        die "Unknown --start-from step: ${start_from} (expected one of: ${STEP_ORDER[*]})"
fi

declare -A skip=()
declare -A required=()
for i in "${!STEP_ORDER[@]}"; do
    name="${STEP_ORDER[$i]}"
    is_excluded=${explicit_skip[${name}]:-0}
    is_prior=0
    if [[ ${start_from_idx} -ge 0 && $i -lt ${start_from_idx} ]]; then
        is_prior=1
    fi

    if [[ ${is_excluded} -eq 1 || ${is_prior} -eq 1 ]]; then
        skip[${name}]=1
    else
        skip[${name}]=0
    fi

    if [[ ${is_prior} -eq 1 && ${is_excluded} -eq 0 ]]; then
        required[${name}]=1
    else
        required[${name}]=0
    fi
done

skip_fastqc=${skip[fastqc]}
skip_trim=${skip[trim]}
skip_species_id=${skip[species-id]}
skip_flye=${skip[flye]}
skip_plassembler=${skip[plassembler]}
skip_consensus=${skip[consensus]}
skip_polish=${skip[polish]}
skip_assembly_qc=${skip[assembly-qc]}
skip_amrfinder=${skip[amrfinder]}
skip_multiqc=${skip[multiqc]}

require_trim=${required[trim]}
require_flye=${required[flye]}
require_plassembler=${required[plassembler]}
require_consensus=${required[consensus]}
require_polish=${required[polish]}

# --skip-consensus (an explicit exclusion, not a --start-from side effect) always
# means: polish Flye and Plassembler assemblies separately, one AMRFinderPlus run each.
separate_mode=${explicit_skip[consensus]:-0}

input_fastq=$(realpath "${args[0]}")
[[ -f "${input_fastq}" ]] || die "FASTQ not found: ${input_fastq}"

mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")
cd "${OUTPUT_DIR}"

# Strip .fastq.gz or .fastq to get a bare sample ID
sample_id=$(basename "${input_fastq}" .fastq.gz)
sample_id=$(basename "${sample_id}" .fastq)
trimmed_fastq_base="trim/${sample_id}.trimmed"
trimmed_fastq="${trimmed_fastq_base}.fastq"

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
elif [[ -z "$(find_fastq "${trimmed_fastq_base}")" && ${require_trim} -eq 1 ]]; then
    die "Trimmed FASTQ not found: ${trimmed_fastq_base}.fastq(.gz) (required because --start-from assumes trimming already ran)"
fi

# Always prefer trimmed reads when available, even if trimming is skipped.
# (Trimmed reads may have been gzipped by a prior run's compression step.)
existing_trimmed_fastq=$(find_fastq "${trimmed_fastq_base}")
if [[ -n "${existing_trimmed_fastq}" ]]; then
    current_fastq="${existing_trimmed_fastq}"
    current_id=$(basename "${current_fastq}" .fastq.gz)
    current_id=$(basename "${current_id}" .fastq)
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
elif [[ -s "${flye_dir}/assembly.fasta" ]]; then
    log "Step 4a: Flye assembly excluded; using existing assembly: ${flye_dir}/assembly.fasta"
elif [[ ${require_flye} -eq 1 ]]; then
    die "Flye assembly not found: ${flye_dir}/assembly.fasta (required because --start-from assumes Flye already ran)"
else
    log "Step 4a: Flye assembly excluded; no assembly available from this stage"
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

    plassembler_qc_assembly=$(find_plassembler_assembly "${current_id}")
else
    plassembler_qc_assembly=$(find_plassembler_assembly "${current_id}")
    if [[ -n "${plassembler_qc_assembly}" ]]; then
        log "Step 4b: Plassembler excluded; using existing assembly: ${plassembler_qc_assembly}"
    elif [[ ${require_plassembler} -eq 1 && ! -d "assembly/${current_id}/plassembler" ]]; then
        die "Plassembler output not found under assembly/${current_id}/plassembler (required because --start-from assumes Plassembler already ran)"
    else
        log "Step 4b: Plassembler excluded; no existing assembly found (sample may have no plasmids)"
    fi
fi

# ─── Step 5 · Consensus assembly (Autocycler) ─────────────────────────────────
assembly_for_polish=""
declare -a denovo_assemblies=()
declare -a amrfinder_targets=()
CONSENSUSDIR=assembly/${current_id}/consensus
AUTOCYCLERDIR=${CONSENSUSDIR}/autocycler
existing_consensus_assembly="${AUTOCYCLERDIR}/consensus_assembly.fasta"

if [[ ${separate_mode} -eq 1 ]]; then
    log "Step 5: Consensus excluded (--skip-consensus); Flye/Plassembler assemblies will be polished separately"

    flye_assembly="${flye_dir}/assembly.fasta"
    plassembler_assembly="${plassembler_qc_assembly}"

    [[ -s "${flye_assembly}" ]] && denovo_assemblies+=("flye::${flye_assembly}")
    [[ -n "${plassembler_assembly}" && -s "${plassembler_assembly}" ]] && denovo_assemblies+=("plassembler::${plassembler_assembly}")

    [[ ${#denovo_assemblies[@]} -gt 0 ]] || die "No assembly available for polishing/QC"

elif [[ ${skip_consensus} -eq 0 ]]; then
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

    ASSEMBLYDIR=${CONSENSUSDIR}/assemblies
    mkdir -p ${ASSEMBLYDIR}

    for ASSEMBLER in flye plassembler; do
        if [[ "${ASSEMBLER}" == "flye" ]]; then
            ASSEMBLY="assembly/${current_id}/flye/assembly.fasta"
        elif [[ "${ASSEMBLER}" == "plassembler" ]]; then
            ASSEMBLY="${plassembler_qc_assembly}"
        fi
        [[ -n "${ASSEMBLY}" && -s "${ASSEMBLY}" ]] || continue
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

    [[ -s "${assembly_for_polish}" ]] || die "Consensus assembly ${assembly_for_polish} not found"

    log "Using consensus assembly: ${assembly_for_polish}"
    run_assembly_qc_stage \
        "${current_id}.consensus" \
        "${assembly_for_polish}" \
        "${AUTOCYCLERDIR}/consensus_assembly.gfa"

elif [[ -s "${existing_consensus_assembly}" ]]; then
    assembly_for_polish="${existing_consensus_assembly}"
    log "Step 5: Consensus excluded; using existing consensus assembly: ${assembly_for_polish}"

else
    die "Consensus assembly not found: ${existing_consensus_assembly} (required because --start-from assumes consensus already ran)"
fi

# ─── Step 6 · Medaka polishing ────────────────────────────────────────────────
final_assembly="${assembly_for_polish}"
medaka_dir=assembly/${current_id}/medaka
existing_polished_assembly="${medaka_dir}/consensus.fasta"

if [[ ${skip_polish} -eq 0 ]]; then
    if [[ ${separate_mode} -eq 1 ]]; then
        log "Step 6: Medaka polishing (separate Flye/Plassembler outputs)"

        final_assembly=""
        for entry in "${denovo_assemblies[@]}"; do
            assembler_name="${entry%%::*}"
            assembler_assembly="${entry#*::}"
            assembler_medaka_dir="assembly/${current_id}/medaka_${assembler_name}"

            medaka_consensus \
                -i ${current_fastq} \
                -d ${assembler_assembly} \
                -o ${assembler_medaka_dir} \
                -t ${THREADS} -b 25 -f

            polished_assembly="${assembler_medaka_dir}/consensus.fasta"
            [[ -s "${polished_assembly}" ]] || die "Medaka assembly not found under ${assembler_medaka_dir}"

            [[ -z "${final_assembly}" ]] && final_assembly="${polished_assembly}"
            amrfinder_targets+=("${assembler_name}::${polished_assembly}")
            log "Using polished ${assembler_name} assembly: ${polished_assembly}"
            run_assembly_qc_stage "${current_id}.medaka.${assembler_name}" "${polished_assembly}"
        done
    else
        log "Step 6: Medaka polishing"

        medaka_consensus \
            -i ${current_fastq} \
            -d ${assembly_for_polish} \
            -o ${medaka_dir} \
            -t ${THREADS} -b 25 -f

        final_assembly="${medaka_dir}/consensus.fasta"
        [[ -s "${final_assembly}" ]] || die "Medaka assembly not found under ${medaka_dir}"

        log "Using polished assembly: ${final_assembly}"
        run_assembly_qc_stage "${current_id}.medaka" "${final_assembly}"
        amrfinder_targets+=("final::${final_assembly}")
    fi
elif [[ ${separate_mode} -eq 1 ]]; then
    log "Step 6: Medaka polishing excluded; looking for existing per-assembler polished assemblies"

    final_assembly=""
    for entry in "${denovo_assemblies[@]}"; do
        assembler_name="${entry%%::*}"
        assembler_assembly="${entry#*::}"
        assembler_medaka_dir="assembly/${current_id}/medaka_${assembler_name}"
        polished_assembly="${assembler_medaka_dir}/consensus.fasta"

        if [[ -s "${polished_assembly}" ]]; then
            [[ -z "${final_assembly}" ]] && final_assembly="${polished_assembly}"
            amrfinder_targets+=("${assembler_name}::${polished_assembly}")
            log "Using existing polished ${assembler_name} assembly: ${polished_assembly}"
        elif [[ ${require_polish} -eq 1 ]]; then
            die "Medaka assembly not found under ${assembler_medaka_dir} (required because --start-from assumes polishing already ran)"
        else
            [[ -z "${final_assembly}" ]] && final_assembly="${assembler_assembly}"
            amrfinder_targets+=("${assembler_name}::${assembler_assembly}")
            log "No polished ${assembler_name} assembly found; using unpolished assembly: ${assembler_assembly}"
        fi
    done
elif [[ -s "${existing_polished_assembly}" ]]; then
    final_assembly="${existing_polished_assembly}"
    log "Using existing polished assembly: ${final_assembly}"
    amrfinder_targets+=("final::${final_assembly}")
elif [[ ${require_polish} -eq 1 ]]; then
    die "Medaka assembly not found under ${medaka_dir} (required because --start-from assumes polishing already ran)"
else
    amrfinder_targets+=("final::${assembly_for_polish}")
    log "No polished assembly found; using unpolished assembly: ${assembly_for_polish}"
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

    [[ ${#amrfinder_targets[@]} -gt 0 ]] || die "No assembly available for AMRFinderPlus"

    for entry in "${amrfinder_targets[@]}"; do
        assembly_name="${entry%%::*}"
        assembly_path="${entry#*::}"
        report_path="annotate/${current_id}/amrfinder/${current_id}.${assembly_name}.tsv"

        amrfinder \
            -n ${assembly_path} \
            -d ${AMRFINDER_DB} \
            --threads ${THREADS} \
            > ${report_path}
    done
fi

# ─── Step 10 · MultiQC summary ─────────────────────────────────────────────────
if [[ ${skip_multiqc} -eq 0 ]]; then
    log "Step 10: MultiQC"

    mkdir -p multiqc

    multiqc \
        -o multiqc \
        -f \
        --fullnames \
        qc
fi

# ─── Step 11 · Compress uncompressed FASTQ files ──────────────────────────────
log "Step 11: Compressing uncompressed FASTQ files"

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
