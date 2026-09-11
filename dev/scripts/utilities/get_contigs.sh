#!/bin/bash

set -euo pipefail

awk -v FS="\t" -v OFS="," '
    NR == 1 {
        print "name", "length", "topology"
    }
    NR > 1 {
        topology = "linear";
        if ( $4 == "Y" ) { topology = "circular" };
        print $1, $2, topology
    }' "${1}" > "${2}"
