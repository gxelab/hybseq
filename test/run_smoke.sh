#!/usr/bin/env bash
# Smoke test for the hybseq pipeline.
#
# Validates pipeline syntax, channel wiring, module I/O contracts and the bin script
# with tiny synthetic inputs (neutral names). Everything runs with -stub-run (process
# stubs), so no containers or real tools are needed and NO biological results are
# produced -- this is a connectivity/syntax test, not a data analysis. Real-data
# validation is documented in the README.
#
# Usage: bash test/run_smoke.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

NF=./nextflow
PROFILE=test
RESULTS=test/results

# Rewrite the cactus seqfile with absolute paths so a (hypothetical) real cactus run
# would also resolve them; the content is not read under -stub-run.
printf 'ref_a\t%s/test/data/ref_a.fa\nref_b\t%s/test/data/ref_b.fa\n' "$REPO_ROOT" "$REPO_ROOT" > test/assemblies.txt

clean() {
    rm -rf work "$RESULTS"
}
clean
trap clean EXIT

pass() { echo "  [ok] $1"; }
fail() { echo "  [FAIL] $1" >&2; exit 1; }

echo "== 1. bin script: syntax + unit check =="
bash -n bin/rename_gtf_for_vg.sh
tmp_gtf=$(mktemp)
bash bin/rename_gtf_for_vg.sh test/data/ref_a.gtf ref_a "$tmp_gtf"
diff "$tmp_gtf" test/data/expected_ref_a.gtf
rm -f "$tmp_gtf"
pass "rename_gtf_for_vg.sh rewrites NC_/NW_ prefixes to ref_a#0#"

echo "== 2. config parses =="
$NF config . -profile $PROFILE > /dev/null
pass "nextflow config"

echo "== 3. DNA-seq workflow (stub-run) =="
$NF run . -profile $PROFILE -stub-run --run dnaseq
for f in \
    "$RESULTS/gfa/test.gbz" "$RESULTS/gfa/test.dist" "$RESULTS/gfa/test.dist2" \
    "$RESULTS/gfa/test.shortread.withzip.min" "$RESULTS/gfa/test.shortread.zipcodes" \
    "$RESULTS/gfa/test.snarls" "$RESULTS/gfa/test.log" \
    "$RESULTS/gam/combined.vcf.gz" "$RESULTS/gam/combined.vcf.tsv.gz"
do
    [ -f "$f" ] || fail "missing expected output: $f"
done
for s in dna_a dna_b dna_hybrid dna_c; do
    for ext in gam gam.log pack vcf.gz vcf.gz.csi; do
        [ -f "$RESULTS/gam/$s.$ext" ] || fail "missing expected output: $RESULTS/gam/$s.$ext"
    done
done
[ ! -d "$RESULTS/rna" ] || fail "dnaseq workflow must not create results/rna"
pass "dnaseq workflow outputs complete"
clean

echo "== 4. RNA-seq workflow (stub-run) =="
$NF run . -profile $PROFILE -stub-run --run rnaseq
for f in \
    "$RESULTS/gfa/ref_a.gtf" \
    "$RESULTS/gfa/test_spliced.pg" "$RESULTS/gfa/test_spliced.xg" \
    "$RESULTS/gfa/test_spliced.gcsa" "$RESULTS/gfa/test_spliced.gcsa.lcp" \
    "$RESULTS/gfa/test_spliced.snarls" "$RESULTS/gfa/test_spliced.dist"
do
    [ -f "$f" ] || fail "missing expected output: $f"
done
for s in rna_a rna_b; do
    for ext in gam pack vcf.gz vcf.tsv.gz; do
        [ -f "$RESULTS/rna/$s.$ext" ] || fail "missing expected output: $RESULTS/rna/$s.$ext"
    done
done
[ -z "$(find "$RESULTS" -name '*pruned.pg*' -print -quit)" ] || fail "pruned.pg must not be published (ephemeral intermediate)"
[ ! -d "$RESULTS/gam" ] || fail "rnaseq workflow must not create results/gam"
pass "rnaseq workflow outputs complete"
clean

echo "== 5. default entry: both assays, shared reference (stub-run) =="
# .dot needs no graphviz (unlike .svg); written into the gitignored test/results/
$NF run . -profile $PROFILE -stub-run -with-dag "$RESULTS/dag.dot"
[ -f "$RESULTS/gam/combined.vcf.tsv.gz" ] || fail "missing DNA output in default entry"
[ -f "$RESULTS/rna/rna_a.vcf.tsv.gz" ] || fail "missing RNA output in default entry"
[ -s "$RESULTS/dag.dot" ] || fail "missing DAG file"
grep -q 'CACTUS_PANGENOME' "$RESULTS/dag.dot" || fail "CACTUS_PANGENOME missing from DAG"
pass "default entry outputs complete (reference workflow invoked once for both assays)"
clean

echo "== 6. negative tests =="
tmp_ss=$(mktemp)
{
    echo "sample,assay,fastq_1,fastq_2"
    echo "dna_x,bogus,test/data/fastq/dna_a_R1.fq.gz,test/data/fastq/dna_a_R2.fq.gz"
} > "$tmp_ss"
if out=$($NF run . -profile $PROFILE -stub-run --run dnaseq --samplesheet "$tmp_ss" 2>&1); then
    fail "bogus assay was accepted"
fi
echo "$out" | grep -q "Unknown assay" || fail "bogus assay: error message missing"
clean

{
    echo "sample,assay,fastq_1,fastq_2"
    echo "dna_x,dnaseq,,test/data/fastq/dna_a_R2.fq.gz"
} > "$tmp_ss"
if out=$($NF run . -profile $PROFILE -stub-run --run dnaseq --samplesheet "$tmp_ss" 2>&1); then
    fail "row with missing fastq_1 was accepted"
fi
echo "$out" | grep -q "fastq_1 missing" || fail "missing fastq_1: error message missing"
clean
rm -f "$tmp_ss"
pass "negative tests"

echo "ALL SMOKE TESTS PASSED"
