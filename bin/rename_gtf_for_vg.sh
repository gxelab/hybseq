#!/usr/bin/env bash
# Rewrite contig names in a GTF so they match the haplotype-path names inside the
# Minigraph-Cactus GBZ (<sample>#0#<contig>, where #0 is the reference-haplotype tag
# produced by cactus --haplo).
#
# The contig field (first tab-delimited field) of every non-comment line is
# prefixed, so the rule works for any contig naming scheme (NC_/NW_ accessions,
# chr*, scaffold*, Ensembl-style names, ...). Lines starting with '#' (comments),
# blank lines, and lines with no tab (not GTF) are passed through unchanged.
# Do not feed an already-renamed GTF: it would be prefixed a second time.
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

# Escape the sample name for use as a sed replacement (backslash first).
escaped_sample=${sample//\\/\\\\}
escaped_sample=${escaped_sample//&/\\&}
escaped_sample=${escaped_sample//|/\\|}

# Anchor on the whole first field plus its tab, so the prefix lands at the start of
# the contig field and '#'-comments, blank lines and non-tab lines are left alone.
sed -E "s/^([^#][^\t]*\t)/${escaped_sample}#0#\1/" "$input" > "$output"
