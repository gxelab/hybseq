# hybseq — development notes

Workflow internals, provenance, design decisions, known limitations, and the smoke test for the hybseq Nextflow pipeline. User-facing documentation lives in [`README.md`](README.md); session rules for AI agents live in [`AGENTS.md`](AGENTS.md).

## Repository layout

| Path | Purpose |
|---|---|
| `main.nf` | Entry workflow: samplesheet parsing/validation, `REFERENCE` once, then `DNASEQ` and/or `RNASEQ` by `params.run` |
| `workflows/reference.nf` | Graph construction and indexing, shared by both assays |
| `workflows/dnaseq.nf` | DNA-seq / pool-seq workflow |
| `workflows/rnaseq.nf` | RNA-seq workflow |
| `modules/local/reference/` | `CACTUS_PANGENOME`, `VG_INDEX_DIST2` |
| `modules/local/dnaseq/` | `VG_GIRAFFE`, `VG_PACK_DNA`, `VG_CALL_DNA`, `BCFTOOLS_INDEX`, `BCFTOOLS_MERGE_DNA`, `BCFTOOLS_QUERY_DNA` |
| `modules/local/rnaseq/` | `RENAME_GTF`, `VG_RNA`, `VG_INDEX_XG`, `VG_PRUNE`, `VG_INDEX_GCSA`, `VG_SNARLS`, `VG_INDEX_DIST_RNA`, `VG_MPMAP`, `VG_PACK_RNA`, `VG_CALL_RNA`, `BCFTOOLS_QUERY_RNA` |
| `bin/rename_gtf_for_vg.sh` | Rewrites `NC_`/`NW_` GTF contig prefixes to `<sample>#0#` |
| `assets/` | Placeholder samplesheet and cactus seqfile |
| `test/` | Synthetic stub-run inputs and `run_smoke.sh` |
| `nextflow.config` | All params, per-process resources, `docker`/`apptainer`/`test` profiles |

## Workflow architecture

```
main.nf
  ├── makeSamplesChannel()            // parse + validate samplesheet
  ├── REFERENCE(assemblies)           // once per run
  │     └── CACTUS_PANGENOME → VG_INDEX_DIST2
  ├── DNASEQ(...)   if params.run in ['both','dnaseq']
  └── RNASEQ(...)   if params.run in ['both','rnaseq']   // requires params.gtf
```

**I/O contracts**

| Workflow | take | emit |
|---|---|---|
| `REFERENCE` | `assemblies` (seqfile path) | `gbz`, `dist`, `min`, `zipcodes`, `snarls`, `dist2` |
| `DNASEQ` | `gbz`, `dist`, `min`, `zipcodes`, `snarls`, `samples` (`[sample,assay,fq1,fq2,idx]`) | `combined_vcf`, `tsv` |
| `RNASEQ` | `gbz`, `gtf`, `samples` | `tsv` |

Each workflow filters `samples` by assay and errors out with `No dnaseq/rnaseq samples found in the samplesheet` when the filtered channel is empty. `params.run` is validated in `main.nf` (`both|dnaseq|rnaseq`), and `--run rnaseq|both` without `--gtf` errors before any task starts.

## Reference workflow

`CACTUS_PANGENOME` runs `cactus-pangenome` as a managed foreground task:

```bash
cactus-pangenome js <assemblies> \
    --outDir gfa --outName <outname> --reference <ref_name> \
    --giraffe clip filter --gbz clip filter full --gfa clip filter full --vcf \
    --permissiveContigFilter --haplo --chrom-vg clip filter --chrom-og full --viz \
    --consCores 16 --indexCores 16 --mgCores 16 --mapCores 4 \
    --batchSystem single_machine --logFile gfa/<outname>.log
gunzip -c gfa/<outname>.gfa.gz > gfa/<outname>.gfa
gunzip -c gfa/<outname>.vcf.gz > gfa/<outname>.vcf
```

The Toil jobstore (`js/`) stays in the Nextflow work directory and is never published.

`VG_INDEX_DIST2` then runs `vg index -j gfa/<outname>.dist2 gfa/<outname>.gbz`. This rebuilds a distance index (as the notebook did after cactus), but **mapping uses the cactus-produced `<outname>.dist`**; `dist2` is published and consumed by nothing downstream.

## DNA-seq (pool-seq) workflow

Per sample (single-end uses one `-f`, paired-end two, both inferred from the samplesheet):

```bash
vg giraffe -Z <outname>.gbz -m <outname>.shortread.withzip.min \
    -z <outname>.shortread.zipcodes -d <outname>.dist -f <fq1> [-f <fq2>] \
    -t <cpus> > gam/<sample>.gam 2> gam/<sample>.gam.log

vg pack -x <outname>.gbz -g gam/<sample>.gam -o gam/<sample>.pack -t <cpus> -Q 5

vg call <outname>.gbz -r <outname>.snarls -k gam/<sample>.pack \
    -s <sample> -z -a -t <cpus> | bgzip -c > gam/<sample>.vcf.gz
```

Then, across all DNA samples:

```bash
bcftools index -f gam/<sample>.vcf.gz
bcftools merge <all vcfs> -O z -o gam/combined.vcf.gz
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT\t%DP\t%AD{0}\t%AD{1}\t%GQ]\n' \
    gam/combined.vcf.gz | gzip -c > gam/combined.vcf.tsv.gz
```

**Merge order.** Each samplesheet row carries a zero-based `idx` through the whole workflow. Before `BCFTOOLS_MERGE_DNA` the `(sample, vcf, csi, idx)` tuples are collected with `toSortedList { a, b -> a[3] <=> b[3] }` and reshaped into one tuple of lists, so the merged VCF sample columns follow samplesheet order regardless of task completion order. The notebook relied on shell-glob order (`for i in *.vcf.gz`), which is not reproducible across filesystems.

`BCFTOOLS_INDEX` re-emits the VCF purely to carry it into the merge; its `publishDir` pattern is `*.csi`, so only the index is published there (the VCF itself was already published by `VG_CALL_DNA`).

## RNA-seq workflow

```bash
# 1. rename GTF contigs to match haplotype-path names in the GBZ
sed -E 's/^(NC_|NW_)/<ref_name>#0#\1/' <gtf> > gfa/<ref_name>.gtf

# 2. spliced graph + indexes
vg rna -p --threads <cpus> --transcripts gfa/<ref_name>.gtf --use-hap-ref --gbz-format <gbz> \
    > gfa/<outname>_spliced.pg
vg index -x gfa/<outname>_spliced.xg gfa/<outname>_spliced.pg
vg prune gfa/<outname>_spliced.pg > <outname>_spliced.pruned.pg
vg index -t <cpus> -g gfa/<outname>_spliced.gcsa -b tmp <outname>_spliced.pruned.pg
vg snarls gfa/<outname>_spliced.pg > gfa/<outname>_spliced.snarls
vg index -j gfa/<outname>_spliced.dist gfa/<outname>_spliced.xg

# 3. per sample
vg mpmap -x gfa/<outname>_spliced.xg -g gfa/<outname>_spliced.gcsa -d gfa/<outname>_spliced.dist \
    -n RNA -l short -F GAM -t <cpus> -f <fq1> [-f <fq2>] > rna/<sample>.gam
vg pack -x gfa/<outname>_spliced.xg -g rna/<sample>.gam -o rna/<sample>.pack -t <cpus> -Q 5
vg call gfa/<outname>_spliced.xg -r gfa/<outname>_spliced.snarls -k rna/<sample>.pack \
    -s <sample|-rna_call_sample> -z -a -t <cpus> | bgzip -c > rna/<sample>.vcf.gz
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT\t%DP\t%AD{0}\t%AD{1}\t%GQ]\n' \
    rna/<sample>.vcf.gz | gzip -c > rna/<sample>.vcf.tsv.gz
```

Notes:

- `vg rna` output is a PackedGraph; `--gbz-format` is passed because the input is GBZ.
- `vg prune` output is an ephemeral intermediate for the GCSA build: `VG_PRUNE` declares no `publishDir`, so `*_spliced.pruned.pg` never reaches `--outdir`. `VG_INDEX_GCSA` emits both `.gcsa` and `.gcsa.lcp` because `vg mpmap` needs them adjacent.
- Snarls are not needed by `vg mpmap` when a dist index is supplied, but they are required by `vg call`; `trivial.snarls` is not needed to build the distance index.
- No merge is performed for RNA-seq: each sample yields its own `.vcf.gz` + `.vcf.tsv.gz`, and cross-sample comparison is a downstream concern.

## Processes and published outputs

| Process | Tool | Published (relative to `--outdir`) |
|---|---|---|
| `CACTUS_PANGENOME` | `cactus-pangenome` | `gfa/<outname>.gbz/.dist/.shortread.withzip.min/.shortread.zipcodes/.snarls/.log/.gfa/.vcf` |
| `VG_INDEX_DIST2` | `vg index -j` | `gfa/<outname>.dist2` |
| `VG_GIRAFFE` | `vg giraffe` | `gam/<sample>.gam`, `gam/<sample>.gam.log` |
| `VG_PACK_DNA` | `vg pack` | `gam/<sample>.pack` |
| `VG_CALL_DNA` | `vg call \| bgzip` | `gam/<sample>.vcf.gz` |
| `BCFTOOLS_INDEX` | `bcftools index` | `gam/<sample>.vcf.gz.csi` |
| `BCFTOOLS_MERGE_DNA` | `bcftools merge` | `gam/combined.vcf.gz` |
| `BCFTOOLS_QUERY_DNA` | `bcftools query` | `gam/combined.vcf.tsv.gz` |
| `RENAME_GTF` | host `sed` (`bin/rename_gtf_for_vg.sh`) | `gfa/<ref_name>.gtf` |
| `VG_RNA` | `vg rna` | `gfa/<outname>_spliced.pg` |
| `VG_INDEX_XG` | `vg index -x` | `gfa/<outname>_spliced.xg` |
| `VG_PRUNE` | `vg prune` | **nothing** (work dir only) |
| `VG_INDEX_GCSA` | `vg index -g -b` | `gfa/<outname>_spliced.gcsa`, `.gcsa.lcp` |
| `VG_SNARLS` | `vg snarls` | `gfa/<outname>_spliced.snarls` |
| `VG_INDEX_DIST_RNA` | `vg index -j` | `gfa/<outname>_spliced.dist` |
| `VG_MPMAP` | `vg mpmap` | `rna/<sample>.gam` |
| `VG_PACK_RNA` | `vg pack` | `rna/<sample>.pack` |
| `VG_CALL_RNA` | `vg call \| bgzip` | `rna/<sample>.vcf.gz` |
| `BCFTOOLS_QUERY_RNA` | `bcftools query` | `rna/<sample>.vcf.tsv.gz` |

## Provenance: pipeline stage → source notebook

`notebooks/` is gitignored but is the source of truth for every stage. Line numbers refer to the current copies.

| Pipeline stage | Source |
|---|---|
| Pangenome graph build (`cactus-pangenome`) | [`notebooks/S05_pangenome_poolseq.qmd`](notebooks/S05_pangenome_poolseq.qmd) L23 |
| `<outname>.dist2` rebuild | S05 L26 |
| DNA mapping (`vg giraffe`) | S05 L30–65; clean one-line version L208 |
| Read support (`vg pack -Q 5`) | S05 L70–94; L211 |
| Variant calling (`vg call -z -a`) | S05 L99–101; L214 |
| `bcftools index` → `merge` → `query` (DNA) | S05 L104–108; L217–221 |
| GTF contig rename | [`notebooks/S06_pangenome_rnaseq.qmd`](notebooks/S06_pangenome_rnaseq.qmd) L22 |
| Spliced graph (`vg rna`) | S06 L25 |
| xg / prune / GCSA indexes | S06 L28–34 |
| Snarls (`vg snarls`) | S06 L37 |
| Spliced dist index | S06 L40 |
| RNA mapping (`vg mpmap`) | S06 L44–51 |
| RNA pack / call / query | S06 L54–58 |
| Long-read RNA-seq | [`notebooks/S07_pangenome_lrs.qmd`](notebooks/S07_pangenome_lrs.qmd) — empty (stub notebook); not implemented |
| Downstream R analysis (site filtering, grenedalf comparison, plots) | S05 L111–202 — exploratory, deliberately **not** pipelined |

## Design decisions and deviations

- **Both layouts in both assays.** The notebooks show paired-end DNA (S05) and single-end RNA (S06); the pipeline infers single-end vs paired-end per row for both assays from an empty `fastq_2`.
- **Deterministic merge order.** DNA merge follows samplesheet order via the carried row index instead of shell-glob order (see DNA-seq section).
- **Explicit mapper format.** `vg mpmap` passes `-F GAM` explicitly; the notebook showed it inconsistently.
- **Graph artifacts as plain files.** The notebook ran cactus manually under apptainer and used the `.gz` artifacts; the pipeline `gunzip`s `<outname>.gfa.gz` → `.gfa` and `<outname>.vcf.gz` → `.vcf` so the published layout is uniform and predictable.
- **`ref_name` couples cactus and the GTF.** `params.ref_name` is both `cactus --reference` and the GTF rename prefix, so it must equal the reference sample name in the seqfile, or the renamed GTF contig names will not match GBZ haplotype paths.
- **No replicate and no per-sample reference.** Derived from what the notebooks actually require: every sample maps to the one shared graph, and pooled samples are a single sample (pooling happens upstream).
- **Threshold exposure.** The hard-coded `-Q 5` pack filter became `params.min_mapq` (default 5, i.e. the notebook value); thread counts, cactus cores, and `gcsa_tmpdir` are likewise params.
- **Pass-through flags.** `--permissiveContigFilter`, `--haplo`, `--chrom-vg clip filter`, `--chrom-og full`, and `--viz` are passed to cactus exactly as in the notebook.
- **`RENAME_GTF` is container-less.** The rename is a one-line `sed`; it uses the helper in `bin/` (available on `PATH` in the Nextflow script environment) rather than pulling a container.

## Known limitations

- `vg call | bgzip -c` requires `bgzip` inside the vg image. vgteam images ship htslib tools; verify at the first real run or override `--vg_container`.
- `<outname>.dist2` is built as a separate step but consumed by nothing: mapping uses the cactus-produced `<outname>.dist`, and `dist2` is only published.
- `--viz` is passed to `cactus-pangenome` as is; the cactus image must provide its dependencies (graphviz/odgi).
- `bcftools_container` uses a placeholder tag (`1.19--h3ea31c5_0`). Pin it to a tag available on your cluster; no particular version is required.
- Cactus resources are estimates (64 GB RAM / 48 h per process) and must be tuned to the assemblies; all other processes default to 1 CPU, 8 GB, 8 h, with per-process overrides in `nextflow.config`.
- `publish_dir_mode` defaults to `copy`, which will copy large GAM/pack files. `symlink` is recommended for real runs.
- `RENAME_GTF` runs on the host with `sed` and only rewrites `NC_` / `NW_` prefixes. Other contig naming schemes need a different rename rule, and `bin/` must be reachable on `PATH`.
- RNA-seq has no merge step by design, so cross-sample RNA comparison must be done downstream of the per-sample TSVs.
- No biological validation is automated in this repository: the smoke test is stub-only, and real-data correctness is established by the notebook comparisons (e.g. S05's grenedalf cross-check), not by CI.
- The `test` profile's `withName` overrides exist because profile-level params do not propagate into the base `process` block; keep them in sync when adding processes.

## Smoke test

```bash
bash test/run_smoke.sh
```

Requirements: the vendored Nextflow launcher at `./nextflow` and Nextflow ≥ 24.10. Everything runs with `-stub-run` — process stubs only, no containers, **no biological results are produced**; this is a connectivity/syntax test.

Inputs are tiny synthetic files with neutral names (`ref_a`/`ref_b`, `dna_a`/`dna_b`/`dna_hybrid`/`dna_c`, `rna_a`/`rna_b`) covering all four assay/layout combinations. The script rewrites `test/assemblies.txt` with absolute paths (the content is not read under `-stub-run`) and removes `work/` and `test/results/` before and after the run (`trap clean EXIT`).

Steps:

1. **bin script unit check** — `bash -n bin/rename_gtf_for_vg.sh`, then run it on `test/data/ref_a.gtf` and `diff` the result against `test/data/expected_ref_a.gtf` (verifies `NC_`/`NW_` → `ref_a#0#` rewriting).
2. **Config parse** — `./nextflow config . -profile test`.
3. **`--run dnaseq`** — asserts `gfa/test.gbz`, `.dist`, `.dist2`, `.shortread.withzip.min`, `.shortread.zipcodes`, `.snarls`, `.log`, `gam/combined.vcf.gz`, `gam/combined.vcf.tsv.gz`, and `gam/<sample>.{gam,gam.log,pack,vcf.gz,vcf.gz.csi}` for `dna_a`, `dna_b`, `dna_hybrid`, `dna_c`; also asserts that `results/rna` is **not** created.
4. **`--run rnaseq`** — asserts `gfa/ref_a.gtf`, `gfa/test_spliced.{pg,xg,gcsa,gcsa.lcp,snarls,dist}` and `rna/<sample>.{gam,pack,vcf.gz,vcf.tsv.gz}` for `rna_a`, `rna_b`; asserts no `*pruned.pg*` is published (ephemeral intermediate) and that `results/gam` is **not** created.
5. **Default entry (`both`)** — one run with `-with-dag`, asserting both DNA and RNA outputs exist, the DAG file is non-empty, and `CACTUS_PANGENOME` appears exactly once in the DAG (single shared reference).
6. **Negative tests** — samplesheets with a bogus assay and with an empty `fastq_1` must fail, and the error messages must contain `Unknown assay` / `fastq_1 missing`.

Success ends with `ALL SMOKE TESTS PASSED`.

## Real-data runs

See the Usage section of [`README.md`](README.md) for command lines. Additional practical notes:

- **Apptainer.** Pass local SIFs via `--vg_sif`, `--cactus_sif`, `--bcftools_sif`; they override the registry images. `apptainer.autoMounts = true` is set, so bind-mount data directories the same way you would with `apptainer exec`, e.g. `-c 'apptainer.runOptions = "--bind /nfs_data"'`.
- **Resources.** Raise `cactus_cons_cores` / `cactus_index_cores` / `cactus_mg_cores` / `cactus_map_cores` and the `CACTUS_PANGENOME` memory/time to match the assemblies.
- **Publishing.** Use `--publish_dir_mode symlink` when GAM and pack files are large.
- **Naming.** `--ref_name` must match the reference name in the seqfile; `--outname` prefixes every graph artifact.

## Open items

- Long-read RNA-seq (S07) is an empty notebook stub and has no implementation.
- `TODO.md` (gitignored) holds the original build brief; its requirements are satisfied except for long-read support. The provenance table above fills the documentation-provenance requirement.
