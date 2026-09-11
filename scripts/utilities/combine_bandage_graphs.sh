#!/bin/bash

set -euo pipefail

# Merge a Flye assembly-graph GFA and a Plassembler plasmid-graph GFA into
# one Bandage image, each contig labelled with its name and length. Segment/
# link/path IDs are prefixed by source program (flye_/plassembler_) -- the
# two graphs' own ID namespaces don't collide today (Flye: edge_1, edge_2,
# ...; Plassembler: 1, 2, ...) but this keeps it safe if that ever changes,
# and means the --names label below reads as "flye_edge_1" / "plassembler_3"
# rather than a bare, ambiguous number. Bandage draws disconnected graphs
# fine, each component laid out separately, so the two just appear side by
# side in one image.
#
# Bandage's --colour custom/--color <csv> flags, for colouring nodes by
# source instead, are documented but broken in the biocontainers Bandage
# 0.9.0 build -- hence encoding origin in the name prefix instead of colour.
# --names/--lengths label rendering needs a Qt platform plugin unavailable
# in this headless container, hence QT_QPA_PLATFORM=offscreen below (plain,
# unlabelled Bandage image calls don't hit that code path and don't need it).
#
# Usage: combine_bandage_graphs.sh <flye_assembly_graph.gfa> <plassembler_plasmids.gfa> <output.svg>

flye_gfa="${1}"
plassembler_gfa="${2}"
out_svg="${3}"
out_gfa="${out_svg%.svg}.gfa"

{
    echo -e "H\tVN:Z:1.0"
    awk 'BEGIN{OFS="\t"}
        $1=="S"{$2="flye_"$2}
        $1=="L"{$2="flye_"$2; $4="flye_"$4}
        $1=="P"{
            n=split($3, arr, ",");
            out="";
            for (i=1; i<=n; i++) {
                orient=substr(arr[i], length(arr[i]), 1);
                name=substr(arr[i], 1, length(arr[i])-1);
                out=out (i>1?",":"") "flye_" name orient;
            }
            $3=out;
        }
        $1!="H"{print}' "${flye_gfa}"
    awk 'BEGIN{OFS="\t"}
        $1=="S"{$2="plassembler_"$2}
        $1=="L"{$2="plassembler_"$2; $4="plassembler_"$4}
        $1!="H"{print}' "${plassembler_gfa}"
} > "${out_gfa}"

QT_QPA_PLATFORM=offscreen bandage-exec Bandage image \
    "${out_gfa}" \
    "${out_svg}" \
    --names --lengths
