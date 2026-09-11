#!/bin/bash

set -euo pipefail

awk -v FS="\t" -v OFS="," '
    NR == 1 {
        print "contig", "name", "start", "end", "strand", "coverage", "pident"
    }
    NR > 1 {
        strand = "1";
        if ( $5 == "-" ) { strand = "-1" };
        print $2, $6, $3, $4, strand, $16, $17
    }' "${1}" > "${2}"
