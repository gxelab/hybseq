#!/usr/bin/env bash
# Smoke test for the hybseq pipeline.
#
# Validates pipeline syntax, channel wiring, module I/O contracts and the bin script
# with small inputs simulated from the committed references (neutral names). The
# reads are generated first (test/simulate_reads.py) because they are not committed.
# Everything then runs with -stub-run (process stubs), so no containers or real tools
# are needed and NO biological results are produced -- this is a connectivity/syntax
# test, not a data analysis. Real-data validation is documented in development.md.
#
# Usage: bash test/run_smoke.sh
#        READS=500 bash test/run_smoke.sh    # smaller inputs, quicker run
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

echo "== 0. test data: generate the reads from the committed references =="
# test/data/fastq/ is generated, not committed: the reads are simulated from
# test/data/ref_a.fa, ref_b.fa and ref_a.gtf, so a reference change can never
# leave stale reads behind. Re-run this step alone with:
#     python3 test/simulate_reads.py --reads N
READS=${READS:-5000}
command -v python3 >/dev/null 2>&1 || fail "python3 is required to generate the test reads"
python3 test/simulate_reads.py --reads "$READS"
python3 test/simulate_reads.py --check --reads "$READS"
pass "reads generated from ref_a/ref_b/ref_a.gtf and verified against them"

echo "== 1. bin script: syntax + unit check =="
bash -n bin/rename_gtf_for_vg.sh
tmp_gtf=$(mktemp)
bash bin/rename_gtf_for_vg.sh test/data/ref_a.gtf ref_a "$tmp_gtf"
diff "$tmp_gtf" test/data/expected_ref_a.gtf
rm -f "$tmp_gtf"
pass "rename_gtf_for_vg.sh reproduces the NC_/NW_ fixture (ref_a#0#)"

# Contig-agnostic rule: any contig name is prefixed; comment, blank and non-tab
# lines pass through unchanged. Tab-delimited fields matter, so the fixtures use
# literal tabs (printf).
tmp_in=$(mktemp)
printf '#!genome-build test\n\
scaffold_b1\ttest\texon\t1\t2\t.\t+\t.\ttranscript_id "tx1";\n\
chr1\ttest\texon\t3\t4\t.\t+\t.\ttranscript_id "tx1";\n\
CM012345.1\ttest\texon\t5\t6\t.\t+\t.\ttranscript_id "tx1";\n\
1\ttest\texon\t7\t8\t.\t+\t.\ttranscript_id "tx1";\n\
\n\
no tabs here\n' > "$tmp_in"
tmp_out=$(mktemp)
bash bin/rename_gtf_for_vg.sh "$tmp_in" ref_b "$tmp_out"
expected=$(mktemp)
printf '#!genome-build test\n\
ref_b#0#scaffold_b1\ttest\texon\t1\t2\t.\t+\t.\ttranscript_id "tx1";\n\
ref_b#0#chr1\ttest\texon\t3\t4\t.\t+\t.\ttranscript_id "tx1";\n\
ref_b#0#CM012345.1\ttest\texon\t5\t6\t.\t+\t.\ttranscript_id "tx1";\n\
ref_b#0#1\ttest\texon\t7\t8\t.\t+\t.\ttranscript_id "tx1";\n\
\n\
no tabs here\n' > "$expected"
diff "$tmp_out" "$expected"
rm -f "$tmp_in" "$tmp_out" "$expected"
pass "rename_gtf_for_vg.sh prefixes any contig name; comments/blank/non-tab lines pass through"

# A sample name with sed-special characters must be emitted literally.
tmp_in=$(mktemp)
printf 'chr1\ttest\texon\t1\t2\t.\t+\t.\ttranscript_id "tx1";\n' > "$tmp_in"
tmp_out=$(mktemp)
bash bin/rename_gtf_for_vg.sh "$tmp_in" 'a&b' "$tmp_out"
expected=$(mktemp)
printf 'a&b#0#chr1\ttest\texon\t1\t2\t.\t+\t.\ttranscript_id "tx1";\n' > "$expected"
diff "$tmp_out" "$expected"
rm -f "$tmp_in" "$tmp_out" "$expected"
pass "rename_gtf_for_vg.sh handles sed-special characters in the sample name"

echo "== 2. config parses; profiles resolve and directives stay paired =="
$NF config . -profile $PROFILE > /dev/null

# Each execution profile must resolve to `enabled = true` in its own config scope.
# This is a config parse only: no container engine is invoked and no Conda
# environment is created.
assert_scope_enabled() {
    local prof="$1" scope="$2" key="$3" out
    if ! out=$($NF config . -profile "test,$prof" 2>&1); then
        echo "$out" | tail -20
        fail "nextflow config -profile test,$prof failed"
    fi
    printf '%s\n' "$out" | awk -v s="$scope" -v k="$key" '
        $0 == s " {" { inb=1; next }
        inb && $0 == "}" { inb=0 }
        inb && $0 ~ ("^[[:space:]]*" k " = true$") { found=1 }
        END { exit found ? 0 : 1 }' \
        || fail "profile '$prof' does not resolve to $scope.$key = true"
}
assert_scope_enabled docker         docker         enabled
assert_scope_enabled podman         podman         enabled
assert_scope_enabled apptainer      apptainer      enabled

# Retired execution paths must not reappear: cactus needs an AVX2-capable x86-64
# CPU (no macOS emulation provides it) and bioconda's cactus is not usable for
# this pipeline, so the applecontainer and conda profiles were removed. An
# unknown profile name makes `nextflow config` exit non-zero, so this asserts
# both that they are gone and that they stay gone.
for retired in applecontainer conda; do
    if $NF config . -profile "test,$retired" >/dev/null 2>&1; then
        fail "retired profile '$retired' still resolves"
    fi
done

# Every tool-running process declares a container image; there is no conda
# fallback, and RENAME_GTF (host sed) declares neither directive.
mods_with_container=$(grep -rlE '^[[:space:]]*container[[:space:]]' modules/local || true)
[ -n "$mods_with_container" ] || fail "no module declares a container directive"
if grep -qE '^[[:space:]]*(container|conda)[[:space:]]' modules/local/reference/rename_gtf.nf; then
    fail "RENAME_GTF must stay container-less (host sed)"
fi
if grep -rqE '^[[:space:]]*conda[[:space:]]' modules/local; then
    fail "conda directives are retired; containers are the only execution path"
fi
pass "nextflow config; execution profiles resolve and the directive contract holds"

echo "== 3. DNA-seq workflow (stub-run) =="
$NF run . -profile $PROFILE -stub-run --run dnaseq
for f in \
    "$RESULTS/ref/test.gbz" "$RESULTS/ref/test.dist" "$RESULTS/ref/test.dist.bak" \
    "$RESULTS/ref/test.shortread.withzip.min" "$RESULTS/ref/test.shortread.zipcodes" \
    "$RESULTS/ref/test.snarls" "$RESULTS/ref/test.log" \
    "$RESULTS/dna/combined.vcf.gz" "$RESULTS/dna/combined.vcf.tsv.gz"
do
    [ -f "$f" ] || fail "missing expected output: $f"
done
for s in dna_a dna_b dna_hybrid dna_c; do
    for ext in gam gam.log pack vcf.gz vcf.gz.csi; do
        [ -f "$RESULTS/dna/$s.$ext" ] || fail "missing expected output: $RESULTS/dna/$s.$ext"
    done
done
[ -z "$(find "$RESULTS" -name '*.dist2' -print -quit)" ] || fail "dist2 must not be published (intermediate)"
# Retired directory names must never reappear (guards a partial rename in any module).
for d in gfa gam; do
    [ ! -d "$RESULTS/$d" ] || fail "retired output dir results/$d was created"
done
[ ! -d "$RESULTS/rna" ] || fail "dnaseq workflow must not create results/rna"
pass "dnaseq workflow outputs complete"
clean

echo "== 4. RNA-seq workflow (stub-run) =="
$NF run . -profile $PROFILE -stub-run --run rnaseq
for f in \
    "$RESULTS/ref/ref_a.gtf" \
    "$RESULTS/ref/test_spliced.pg" "$RESULTS/ref/test_spliced.xg" \
    "$RESULTS/ref/test_spliced.gcsa" "$RESULTS/ref/test_spliced.gcsa.lcp" \
    "$RESULTS/ref/test_spliced.snarls" "$RESULTS/ref/test_spliced.dist" \
    "$RESULTS/rna/combined.vcf.gz" "$RESULTS/rna/combined.vcf.tsv.gz"
do
    [ -f "$f" ] || fail "missing expected output: $f"
done
for s in rna_a rna_b; do
    for ext in gam pack vcf.gz vcf.gz.csi; do
        [ -f "$RESULTS/rna/$s.$ext" ] || fail "missing expected output: $RESULTS/rna/$s.$ext"
    done
    # RNA VCFs are indexed and merged like the DNA ones; per-sample TSVs are retired.
    [ ! -f "$RESULTS/rna/$s.vcf.tsv.gz" ] || fail "per-sample RNA TSV must not be published: $RESULTS/rna/$s.vcf.tsv.gz"
done
[ -z "$(find "$RESULTS" -name '*pruned.pg*' -print -quit)" ] || fail "pruned.pg must not be published (ephemeral intermediate)"
# Retired directory names must never reappear (guards a partial rename in any module).
for d in gfa gam; do
    [ ! -d "$RESULTS/$d" ] || fail "retired output dir results/$d was created"
done
[ ! -d "$RESULTS/dna" ] || fail "rnaseq workflow must not create results/dna"
pass "rnaseq workflow outputs complete"
clean

echo "== 5. default entry: both assays, shared reference (stub-run) =="
# .dot needs no graphviz (unlike .svg); written into the gitignored test/results/
$NF run . -profile $PROFILE -stub-run -with-dag "$RESULTS/dag.dot"
[ -f "$RESULTS/dna/combined.vcf.tsv.gz" ] || fail "missing DNA output in default entry"
[ -f "$RESULTS/rna/combined.vcf.tsv.gz" ] || fail "missing RNA output in default entry"
[ -s "$RESULTS/dag.dot" ] || fail "missing DAG file"
grep -q 'CACTUS_PANGENOME' "$RESULTS/dag.dot" || fail "CACTUS_PANGENOME missing from DAG"
pass "default entry outputs complete (reference workflow invoked once for both assays)"
clean

echo "== 6. single-sample assays: no index, no merge (stub-run) =="
# One row per assay: each run sees exactly one sample for its assay, so there is
# nothing to merge and the sample VCF must be queried directly.
tmp_ss=$(mktemp)
{
    echo "sample,assay,fastq_1,fastq_2"
    echo "dna_solo,dnaseq,test/data/fastq/dna_a_R1.fq.gz,test/data/fastq/dna_a_R2.fq.gz"
    echo "rna_solo,rnaseq,test/data/fastq/rna_a.fastq.gz"
} > "$tmp_ss"

if out=$($NF run . -profile $PROFILE -stub-run --run dnaseq --samplesheet "$tmp_ss" 2>&1); then :; else
    echo "$out" | tail -20
    fail "single-sample dnaseq run failed"
fi
for f in "$RESULTS/dna/dna_solo.vcf.gz" "$RESULTS/dna/combined.vcf.tsv.gz"; do
    [ -f "$f" ] || fail "missing expected output: $f"
done
for f in "$RESULTS/dna/dna_solo.vcf.gz.csi" "$RESULTS/dna/combined.vcf.gz"; do
    [ ! -f "$f" ] || fail "single sample must not publish: $f"
done
if echo "$out" | grep -q -E 'BCFTOOLS_INDEX|BCFTOOLS_MERGE'; then
    fail "index/merge ran for a single dnaseq sample"
fi
clean

if out=$($NF run . -profile $PROFILE -stub-run --run rnaseq --samplesheet "$tmp_ss" 2>&1); then :; else
    echo "$out" | tail -20
    fail "single-sample rnaseq run failed"
fi
for f in "$RESULTS/rna/rna_solo.vcf.gz" "$RESULTS/rna/combined.vcf.tsv.gz"; do
    [ -f "$f" ] || fail "missing expected output: $f"
done
for f in "$RESULTS/rna/rna_solo.vcf.gz.csi" "$RESULTS/rna/combined.vcf.gz"; do
    [ ! -f "$f" ] || fail "single sample must not publish: $f"
done
if echo "$out" | grep -q -E 'BCFTOOLS_INDEX|BCFTOOLS_MERGE'; then
    fail "index/merge ran for a single rnaseq sample"
fi
clean
rm -f "$tmp_ss"
pass "single-sample assays publish the queried table without index or merge"

echo "== 7. negative tests =="
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
