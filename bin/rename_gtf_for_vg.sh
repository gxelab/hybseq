#!/usr/bin/env bash
# Rewrite contig names in a GTF so they match the haplotype-path names inside the
# Minigraph-Cactus GBZ (<sample>#0#<contig>, where #0 is the reference-haplotype tag
# produced by cactus --haplo). Only lines starting with NC_ or NW_ are rewritten.
#
# Usage: rename_gtf_for_vg.sh <input.gtf> <sample> <output.gtf>
set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input.gtf> <sample> <output.gtf>" >&2
    exit 1
fi

input=$1
sample=$2
output=$3

sed -E "s/^(NC_|NW_)/${sample}#0#\1/" "$input" > "$output"
