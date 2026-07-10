#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)

usage() {
    cat <<'EOF'
Usage: scripts/run_sample_pipeline.sh [options] <fastq>

Run the workshop analysis scripts for a single FASTQ from preprocessing through
assembly, annotation, and reporting.

Options:
  --skip-fastqc         Skip FastQC on the current FASTQ
  --skip-trim           Skip read trimming with fastplong
  --skip-species-id     Skip Kraken2 species identification
  --skip-flye           Skip Flye assembly
  --skip-plassembler    Skip Plassembler assembly
  --skip-consensus      Skip Autocycler consensus assembly
  --skip-polish         Skip Medaka polishing
  --skip-assembly-qc    Skip QUAST, BUSCO, and Bandage QC
  --skip-amrfinder      Skip AMRFinderPlus
  --skip-bakta          Skip Bakta annotation
  --skip-multiqc        Skip final MultiQC summary
  --chromosome-contig   Extract plasmid contigs by excluding this contig header
  -h, --help            Show this help

Notes:
  - Run this from anywhere; it will cd into the repository root.
  - The wrapper preserves each script's basename-based sample ID logic, so the
    sample ID usually changes after trimming (for example sample -> sample.trimmed).
EOF
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

run_stage() {
    local script_name="$1"
    shift

    log "Running ${script_name} $*"
    bash "${SCRIPT_DIR}/${script_name}" "$@"
}

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

args=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-fastqc)
            skip_fastqc=1
            ;;
        --skip-trim)
            skip_trim=1
            ;;
        --skip-species-id)
            skip_species_id=1
            ;;
        --skip-flye)
            skip_flye=1
            ;;
        --skip-plassembler)
            skip_plassembler=1
            ;;
        --skip-consensus)
            skip_consensus=1
            ;;
        --skip-polish)
            skip_polish=1
            ;;
        --skip-assembly-qc)
            skip_assembly_qc=1
            ;;
        --skip-amrfinder)
            skip_amrfinder=1
            ;;
        --skip-bakta)
            skip_bakta=1
            ;;
        --skip-multiqc)
            skip_multiqc=1
            ;;
        --chromosome-contig)
            [[ $# -ge 2 ]] || die "--chromosome-contig requires a value"
            chromosome_contig="$2"
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --*)
            die "Unknown option: $1"
            ;;
        *)
            args+=("$1")
            ;;
    esac
    shift
done

[[ ${#args[@]} -eq 1 ]] || {
    usage >&2
    exit 1
}

input_fastq="${args[0]}"
cd "${REPO_ROOT}"

[[ -f "${input_fastq}" ]] || die "FASTQ not found: ${input_fastq}"

current_fastq="${input_fastq}"
current_id=$(basename "${current_fastq}" .fastq)

if [[ ${skip_fastqc} -eq 0 ]]; then
    run_stage 02.fastqc.sh "${current_fastq}"
fi

if [[ ${skip_trim} -eq 0 ]]; then
    run_stage 03.b.fastplong.sh "${current_fastq}"
    current_fastq="trim/${current_id}.trimmed.fastq"
    [[ -s "${current_fastq}" ]] || die "Trimmed FASTQ missing: ${current_fastq}"
    current_id=$(basename "${current_fastq}" .fastq)
fi

if [[ ${skip_species_id} -eq 0 ]]; then
    run_stage 04.a.species_id.sh "${current_fastq}"
fi

flye_dir="assembly/${current_id}/flye"
if [[ ${skip_flye} -eq 0 ]]; then
    run_stage 05.a.assembly.flye.sh "${current_fastq}"
    [[ -s "${flye_dir}/assembly.fasta" ]] || die "Flye assembly missing: ${flye_dir}/assembly.fasta"
fi

if [[ ${skip_plassembler} -eq 0 ]]; then
    run_stage 05.b.assembly.plassembler.sh "${current_fastq}" "${flye_dir}"
fi

assembly_for_polish=""
if [[ ${skip_consensus} -eq 0 ]]; then
    run_stage 06.a.consensus.subsample.sh "${current_fastq}"
    run_stage 06.b.consensus.autocycler.sh "${current_id}"
    consensus_dir="assembly/${current_id}/consensus"
    if assembly_for_polish=$(find_first_file "${consensus_dir}" "consensus*.fasta" "*.fasta"); then
        log "Using consensus assembly: ${assembly_for_polish}"
    else
        die "Consensus assembly not found under ${consensus_dir}"
    fi
fi

if [[ -z "${assembly_for_polish}" ]]; then
    assembly_for_polish="${flye_dir}/assembly.fasta"
    [[ -s "${assembly_for_polish}" ]] || die "No assembly available for polishing/QC"
fi

final_assembly="${assembly_for_polish}"
if [[ ${skip_polish} -eq 0 ]]; then
    run_stage 08.polish.sh "${current_fastq}" "${assembly_for_polish}"
    medaka_dir="assembly/${current_id}/medaka"
    if final_assembly=$(find_first_file "${medaka_dir}" "consensus*.fasta" "*.fasta"); then
        log "Using polished assembly: ${final_assembly}"
    else
        die "Medaka assembly not found under ${medaka_dir}"
    fi
fi

if [[ ${skip_assembly_qc} -eq 0 ]]; then
    qc_id="${current_id}.final"
    run_stage 07.a.assembly_qc.quast.sh "${final_assembly}" "${qc_id}"
    run_stage 07.b.assembly_qc.busco.sh "${final_assembly}" "${qc_id}"

    bandage_graph=""
    if [[ -s "${flye_dir}/assembly_graph.gfa" ]]; then
        bandage_graph="${flye_dir}/assembly_graph.gfa"
    elif [[ ${skip_consensus} -eq 0 ]]; then
        consensus_graph_dir="assembly/${current_id}/consensus"
        bandage_graph=$(find_first_file "${consensus_graph_dir}" "*.gfa" || true)
    fi

    if [[ -n "${bandage_graph}" ]]; then
        run_stage 07.c.assembly_qc.bandage.sh "${bandage_graph}" "${qc_id}"
    else
        log "Skipping Bandage because no assembly graph was found"
    fi
fi

if [[ -n "${chromosome_contig}" ]]; then
    run_stage 09.a.get_plasmids.sh "${final_assembly}" "${current_id}" "${chromosome_contig}"
fi

if [[ ${skip_amrfinder} -eq 0 ]]; then
    run_stage 09.b.amrfinder.sh "${final_assembly}" "${current_id}"
fi

if [[ ${skip_bakta} -eq 0 ]]; then
    run_stage 09.c.bakta.sh "${final_assembly}" "${current_id}"
fi

if [[ ${skip_multiqc} -eq 0 ]]; then
    run_stage 00.multiqc.sh
fi

log "Pipeline completed for ${current_id}"
log "Final assembly: ${final_assembly}"