#!/bin/bash

set -euo pipefail

# Merge a Flye assembly-graph GFA and a Plassembler plasmid-graph GFA into
# one Bandage image, with segments labelled using the same chromosome/
# plasmid_N names merge_contigs.sh assigned when it built the draft assembly
# (see its <output>.id_map.txt), instead of raw per-program IDs -- so this
# image lines up with QUAST/BUSCO/AMRFinderPlus reports on the same contigs.
#
# The id_map only knows contig_N (Flye) and the bare Plassembler numbers, but
# Flye's own graph segments are edge_N, not contig_N. To bridge that, we read
# the Flye GFA's P lines (each one is a contig_N and the edge path that makes
# it up) and use them to translate edge_N to contig_N before the id_map
# lookup. Plassembler's segment IDs already match its side of the id_map
# directly. Anything with no id_map entry -- on either side -- is dropped
# from the combined image entirely rather than kept under a fallback name,
# so the image only ever shows contigs that are actually in
# flye/assembly_info.txt (or the draft assembly, if merge_contigs.sh's
# --exclude was used to drop one of those from it).
#
# On the Flye side, that filters out raw assembly-graph edges that never
# made it into any final contig's path (repeat copies, alternative/bubble
# paths, low-coverage tips Flye's repeat resolution didn't select) -- these
# were never in assembly.fasta to begin with, they're graph noise upstream
# of the final contigs, not discarded FASTA sequences.
#
# On the Plassembler side, plassembler_plasmids.gfa is not "the final
# plasmids as a graph" -- it's Plassembler's internal Unicycler/SPAdes
# hybrid-assembly graph, i.e. every contig produced before Plassembler's
# plasmid-calling filters (copy number, circularity, depth vs. chromosome,
# PLSDB hits -- see plassembler_summary.tsv) picked which ones were real
# plasmids and wrote only those to plassembler_plasmids.fasta. So the .gfa
# routinely has far more segments than the .fasta has sequences; the extra
# ones are leftover low-depth/small fragments that never made the cut and
# have no id_map entry either.
#
# --names/--lengths label rendering needs a Qt platform plugin unavailable in
# this headless container, hence QT_QPA_PLATFORM=offscreen below (plain,
# unlabelled Bandage image calls don't hit that code path and don't need it).
#
# Usage: combine_bandage_graphs.sh <flye_assembly_graph.gfa> <plassembler_plasmids.gfa> <id_map.txt> <output.svg>

if [[ $# -ne 4 ]]; then
    echo "Usage: $(basename "$0") <flye_assembly_graph.gfa> <plassembler_plasmids.gfa> <id_map.txt> <output.svg>" >&2
    exit 1
fi

flye_gfa="${1}"
plassembler_gfa="${2}"
id_map="${3}"
out_svg="${4}"
out_gfa="${out_svg%.svg}.gfa"

work_dir=$(mktemp -d)
trap 'rm -rf "${work_dir}"' EXIT

# contig_N / bare Plassembler N -> final name, e.g. chromosome, plasmid_1
cp "${id_map}" "${work_dir}/id_map.tsv"

# edge_N -> contig_N, read straight off the Flye graph's own P lines
awk 'BEGIN{FS="\t"; OFS="\t"} $1=="P"{
    n = split($3, edges, ",");
    for (i=1; i<=n; i++) {
        edge = edges[i];
        sub(/[+-]$/, "", edge);
        print edge, $2;
    }
}' "${flye_gfa}" > "${work_dir}/edge_to_contig.tsv"

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
    awk -F'\t' -v OFS='\t' -v idmap="${work_dir}/id_map.tsv" -v edgemap="${work_dir}/edge_to_contig.tsv" "${rename_awk}" "${flye_gfa}"
    awk -F'\t' -v OFS='\t' -v idmap="${work_dir}/id_map.tsv" -v edgemap=/dev/null "${rename_awk}" "${plassembler_gfa}"
} > "${out_gfa}"

QT_QPA_PLATFORM=offscreen bandage-exec Bandage image \
    "${out_gfa}" \
    "${out_svg}" \
    --names --lengths
