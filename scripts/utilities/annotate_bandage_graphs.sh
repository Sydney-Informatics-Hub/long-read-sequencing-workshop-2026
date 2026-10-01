#!/bin/bash

set -euo pipefail

# Usage: annotate_bandage_graphs.sh <assembly_graph.gfa> <output.svg>

if [[ $# -ne 2 ]]; then
    echo "Usage: $(basename "$0") <assembly_graph.gfa> <output.svg>" >&2
    exit 1
fi

assembly_gfa="${1}"
out_svg="${2}"

QT_QPA_PLATFORM=offscreen bandage-exec Bandage image \
    "${assembly_gfa}" \
    "${out_svg}" \
    --names --lengths
