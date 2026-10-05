# hybseq — development notes

Workflow internals, provenance, design decisions, known limitations, and the smoke test for the hybseq Nextflow pipeline. User-facing documentation lives in [`README.md`](README.md); session rules for AI agents live in [`AGENTS.md`](AGENTS.md).

## Repository layout

| Path | Purpose |
|---|---|
| `main.nf` | Entry workflow: samplesheet parsing/validation, `REFERENCE` once, then `DNASEQ` and/or `RNASEQ` by `params.run` |
| `workflows/reference.nf` | Graph construction and indexing, shared by both assays |
| `workflows/dnaseq.nf` | DNA-seq / pool-seq workflow |
| `workflows/rnaseq.nf` | RNA-seq workflow |
| `modules/local/reference/` | Reference-graph steps, whichever workflow invokes them: `CACTUS_PANGENOME`, `VG_INDEX_DIST_UPDATE`, `RENAME_GTF`, `VG_RNA`, `VG_INDEX_XG`, `VG_PRUNE`, `VG_INDEX_GCSA`, `VG_SNARLS`, `VG_INDEX_DIST_RNA` |
| `modules/local/variation/` | `VG_PACK`, `VG_CALL`, `BCFTOOLS_INDEX`, `BCFTOOLS_MERGE`, `BCFTOOLS_QUERY` — shared by both assays (everything from `vg pack` onwards) |
| `modules/local/dnaseq/` | `VG_GIRAFFE` (per-sample DNA mapping) |
| `modules/local/rnaseq/` | `VG_MPMAP` (per-sample RNA mapping) |
| `bin/rename_gtf_for_vg.sh` | Prefixes every non-comment GTF contig with `<sample>#0#`, whatever the contig naming scheme |
| `assets/` | Placeholder samplesheet and cactus seqfile |
| `test/` | Synthetic stub-run inputs and `run_smoke.sh` |
| `nextflow.config` | All params, per-process resources, `docker`/`apptainer`/`test` profiles |

## Workflow architecture

```
main.nf
  ├── makeSamplesChannel()            // parse + validate samplesheet
  ├── REFERENCE(assemblies)           // once per run
  │     └── CACTUS_PANGENOME → VG_INDEX_DIST_UPDATE
  ├── DNASEQ(...)   if params.run in ['both','dnaseq']
  └── RNASEQ(...)   if params.run in ['both','rnaseq']   // requires params.gtf
```

**I/O contracts**

| Workflow | take | emit |
|---|---|---|
| `REFERENCE` | `assemblies` (seqfile path) | `gbz`, `dist`, `min`, `zipcodes`, `snarls` |
| `DNASEQ` | `gbz`, `dist`, `min`, `zipcodes`, `snarls`, `samples` (`[sample,assay,fq1,fq2,idx]`) | `combined_vcf`, `tsv` |
| `RNASEQ` | `gbz`, `gtf`, `samples` | `combined_vcf`, `tsv` |

Each workflow filters `samples` by assay and errors out with `No dnaseq/rnaseq samples found in the samplesheet` when the filtered channel is empty. `params.run` is validated in `main.nf` (`both|dnaseq|rnaseq`), and `--run rnaseq|both` without `--gtf` errors before any task starts. `combined_vcf` is an empty channel when the assay has a single sample (nothing is merged); `tsv` is always emitted.

## Reference workflow

`CACTUS_PANGENOME` runs `cactus-pangenome` as a managed foreground task:

```bash
mkdir -p <outdir>/toil_work
cactus-pangenome <outdir>/toil_work/js <assemblies> \
    --outDir ref --outName <outname> --reference <ref_name> \
    --giraffe clip filter --gbz clip filter full --gfa clip filter full --vcf \
    --permissiveContigFilter --haplo --chrom-vg clip filter --chrom-og full --viz \
    --consCores 16 --indexCores 16 --mgCores 16 --mapCores 4 \
    --batchSystem single_machine --logFile ref/<outname>.log
gunzip -c ref/<outname>.gfa.gz > ref/<outname>.gfa
gunzip -c ref/<outname>.vcf.gz > ref/<outname>.vcf
mv ref/<outname>.dist ref/<outname>.dist.bak
rm -rf <outdir>/toil_work
```

The Toil jobstore is the first positional argument and lives at `<outdir>/toil_work/js`, so all Toil scratch stays under `--outdir`; the task removes `toil_work/` when it finishes. Nothing in `toil_work/` is published. The distance index cactus produced is archived as `<outname>.dist.bak` inside the task, so `CACTUS_PANGENOME` publishes no file named `<outname>.dist` (which avoids racing the promoted index, see below).

`VG_INDEX_DIST_UPDATE` then rebuilds the distance index and promotes it to the primary filename: `vg index -j ref/<outname>.dist2 ref/<outname>.gbz`, then `mv ref/<outname>.dist2 ref/<outname>.dist`. The cactus index was already renamed to `<outname>.dist.bak`, so the published `<outname>.dist` is the rebuilt index and every `dist` consumer (`vg giraffe -d`) uses it — exactly the notebook sequence (S05 L26–28). `<outname>.dist2` is an intermediate and is never published.

## Shared variant-calling steps (pack → call → index → merge → query)

Everything from `vg pack` onwards is common to both assays and exists once, under `modules/local/variation/`. Each module takes `val assay_dir` (`dna`|`rna`), which selects the published subdirectory; `VG_CALL` additionally takes `graph_is_gbz` and `call_sample_name` (see the RNA-seq section).

Per sample (`<assay_dir>` is `dna` for DNA-seq and `rna` for RNA-seq; the graph is the GBZ for DNA-seq and the spliced `xg` for RNA-seq):

```bash
vg pack -x <graph> -g <assay_dir>/<sample>.gam -o <assay_dir>/<sample>.pack -t <cpus> -Q 5

# DNA-seq: -z restricts calling to the GBZ haplotypes
vg call <outname>.gbz -r <outname>.snarls -k <assay_dir>/<sample>.pack \
    -s <sample> -z -a -t <cpus> | bgzip -c > <assay_dir>/<sample>.vcf.gz

# RNA-seq: no -z (it applies only to GBZ input), -s may be overridden
vg call <outname>_spliced.xg -r <outname>_spliced.snarls -k <assay_dir>/<sample>.pack \
    -s <sample> -a -t <cpus> | bgzip -c > <assay_dir>/<sample>.vcf.gz

bcftools index -f <assay_dir>/<sample>.vcf.gz
```

Then, when the assay has more than one sample, across all of them:

```bash
bcftools merge <all vcfs> -O z -o <assay_dir>/combined.vcf.gz
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT\t%DP\t%AD{0}\t%AD{1}\t%GQ]\n' \
    <assay_dir>/combined.vcf.gz | gzip -c > <assay_dir>/combined.vcf.tsv.gz
```

**Single-sample assays skip index and merge.** With exactly one sample there is nothing to merge, and the index exists only to enable merging, so neither `bcftools index` nor `bcftools merge` runs: `BCFTOOLS_QUERY` reads the sample's own VCF directly (`bcftools query` needs no index unless a region is requested) and writes the same `<assay_dir>/combined.vcf.tsv.gz` table. Such an assay publishes no `<sample>.vcf.gz.csi` and no `combined.vcf.gz`. The branch lives in each workflow, not in a module: both sample VCFs are collected once with `toSortedList { a, b -> a[3] <=> b[3] }` (samplesheet order), split into the single-sample and multi-sample cases, the multi-sample list is turned back into per-sample items for `BCFTOOLS_INDEX` with `flatMap { items -> items }` (not `flatten`, which recurses into the tuples), the indexed tuples are re-collected for `BCFTOOLS_MERGE`, and the merged VCF is `mix`ed with the single sample's VCF as the input of `BCFTOOLS_QUERY`. A collecting operator on an empty channel emits an **empty list** (verified on the installed Nextflow), so the multi-sample branch drops it with `filter { !it.isEmpty() }` before reshaping; without that guard the merge would receive an empty tuple of lists.

**Merge order.** With more than one sample, each samplesheet row carries a zero-based `idx` through the whole workflow. Before `BCFTOOLS_MERGE` the `(sample, vcf, csi, idx)` tuples are collected with `toSortedList { a, b -> a[3] <=> b[3] }` and reshaped into one tuple of lists, so the merged VCF sample columns follow samplesheet order regardless of task completion order. The notebook relied on shell-glob order (`for i in *.vcf.gz`), which is not reproducible across filesystems. Keep the explicit `items.collect { it[n] }` reshape: the channel `transpose()` operator is a Nextflow operator (it transposes, it does not split a list into items), not Groovy `List.transpose()`.

`BCFTOOLS_INDEX` re-emits the VCF purely to carry it into the merge; its `publishDir` pattern is `*.csi`, so only the index is published there (the VCF itself was already published by `VG_CALL` under `<assay_dir>/`). Its publish path is the closure `{ "${params.outdir}/${assay_dir}" }`: a closure is the form Nextflow re-evaluates per task for a directive argument, while a `${...}` string inside a `publishDir` attribute is resolved once, when the process is defined, and fails on an input variable.

## DNA-seq (pool-seq) workflow

Per sample (single-end uses one `-f`, paired-end two, both inferred from the samplesheet):

```bash
vg giraffe -Z <outname>.gbz -m <outname>.shortread.withzip.min \
    -z <outname>.shortread.zipcodes -d <outname>.dist -f <fq1> [-f <fq2>] \
    -t <cpus> > dna/<sample>.gam 2> dna/<sample>.gam.log
```

The giraffe log is published as its own process output (`dna/<sample>.gam.log`), not carried inside the sample tuple, so both assays feed the shared chain the same `(sample, file, idx)` shape.

The shared steps then run with `assay_dir = 'dna'`, the GBZ as pack/call graph (`VG_CALL` needs `graph_is_gbz = true`, hence `-z`) and no `vg call -s` override.

## RNA-seq workflow

```bash
# 1. rename GTF contigs to match haplotype-path names in the GBZ (any contig naming scheme)
sed -E 's/^([^#][^\t]*\t)/<ref_name>#0#\1/' <gtf> > ref/<ref_name>.gtf

# 2. spliced graph + indexes
vg rna -p --threads <cpus> --transcripts ref/<ref_name>.gtf --use-hap-ref --gbz-format <gbz> \
    > ref/<outname>_spliced.pg
vg index -x ref/<outname>_spliced.xg ref/<outname>_spliced.pg
vg prune ref/<outname>_spliced.pg > <outname>_spliced.pruned.pg
vg index -t <cpus> -g ref/<outname>_spliced.gcsa -b tmp <outname>_spliced.pruned.pg
vg snarls ref/<outname>_spliced.pg > ref/<outname>_spliced.snarls
vg index -j ref/<outname>_spliced.dist ref/<outname>_spliced.xg

# 3. per sample
vg mpmap -x ref/<outname>_spliced.xg -g ref/<outname>_spliced.gcsa -d ref/<outname>_spliced.dist \
    -n RNA -l short -F GAM -t <cpus> -f <fq1> [-f <fq2>] > rna/<sample>.gam
```

The shared steps then run with `assay_dir = 'rna'`, the spliced `xg` as pack/call graph (`graph_is_gbz = false`, so no `-z`), and `vg call -s` = the sample id unless `--rna_call_sample` is set.

Notes:

- The contig rename is contig-agnostic: it prefixes the first tab-delimited field of every non-comment line, so `NC_`/`NW_` accessions, `chr*`, `scaffold*` and other schemes all work. Lines starting with `#`, blank lines and lines without a tab are passed through unchanged.
- `vg rna` output is a PackedGraph; `--gbz-format` is passed because the input is GBZ.
- `vg prune` output is an ephemeral intermediate for the GCSA build: `VG_PRUNE` declares no `publishDir`, so `*_spliced.pruned.pg` never reaches `--outdir`. `VG_INDEX_GCSA` emits both `.gcsa` and `.gcsa.lcp` because `vg mpmap` needs them adjacent.
- Snarls are not needed by `vg mpmap` when a dist index is supplied, but they are required by `vg call`; `trivial.snarls` is not needed to build the distance index.
- **`vg call -z` is DNA-only.** S06 L56 passes `-z` together with the spliced `xg`, but `-z` restricts calling to the GBZ haplotypes and `vg call` rejects it when the input graph is not a GBZ (see the tool's own check: *"-z can only be used when input graph is in GBZ format"*). The pipeline therefore passes it only for the GBZ (`VG_CALL`'s `graph_is_gbz`), which is also why S06's `-s dsim` is exposed as an opt-in `--rna_call_sample`.
- **`-s` and merging.** `vg call -s` defaults to the sample id (unique per row, which is what `bcftools merge` needs for the sample columns). A fixed `--rna_call_sample` gives every RNA VCF the same sample column, so `bcftools merge` fails on duplicate sample names; use it only for single-sample RNA runs.
- RNA-seq is indexed, merged and queried exactly like DNA-seq (one `combined.vcf.gz`/`combined.vcf.tsv.gz` per assay); cross-assay comparison remains a downstream concern.

## Processes and published outputs

| Process | Tool | Published (relative to `--outdir`) |
|---|---|---|
| `CACTUS_PANGENOME` | `cactus-pangenome` | `ref/<outname>.gbz/.dist.bak/.shortread.withzip.min/.shortread.zipcodes/.snarls/.log/.gfa/.vcf` |
| `VG_INDEX_DIST_UPDATE` | `vg index -j` | `ref/<outname>.dist` |
| `VG_GIRAFFE` | `vg giraffe` | `dna/<sample>.gam`, `dna/<sample>.gam.log` |
| `VG_PACK` | `vg pack` | `dna/<sample>.pack` or `rna/<sample>.pack` |
| `VG_CALL` | `vg call \| bgzip` | `dna/<sample>.vcf.gz` or `rna/<sample>.vcf.gz` |
| `BCFTOOLS_INDEX` | `bcftools index` | `dna/<sample>.vcf.gz.csi` or `rna/<sample>.vcf.gz.csi` — only when that assay has ≥2 samples |
| `BCFTOOLS_MERGE` | `bcftools merge` | `dna/combined.vcf.gz` or `rna/combined.vcf.gz` — only when that assay has ≥2 samples |
| `BCFTOOLS_QUERY` | `bcftools query` | `dna/combined.vcf.tsv.gz` or `rna/combined.vcf.tsv.gz` |
| `RENAME_GTF` | host `sed` (`bin/rename_gtf_for_vg.sh`) | `ref/<ref_name>.gtf` |
| `VG_RNA` | `vg rna` | `ref/<outname>_spliced.pg` |
| `VG_INDEX_XG` | `vg index -x` | `ref/<outname>_spliced.xg` |
| `VG_PRUNE` | `vg prune` | **nothing** (work dir only) |
| `VG_INDEX_GCSA` | `vg index -g -b` | `ref/<outname>_spliced.gcsa`, `.gcsa.lcp` |
| `VG_SNARLS` | `vg snarls` | `ref/<outname>_spliced.snarls` |
| `VG_INDEX_DIST_RNA` | `vg index -j` | `ref/<outname>_spliced.dist` |
| `VG_MPMAP` | `vg mpmap` | `rna/<sample>.gam` |

`VG_PACK`, `VG_CALL`, `BCFTOOLS_INDEX`, `BCFTOOLS_MERGE` and `BCFTOOLS_QUERY` are defined once in `modules/local/variation/` and invoked by both workflows; `assay_dir` picks the published subdirectory.

## Provenance: pipeline stage → source notebook

`notebooks/` is gitignored but is the source of truth for every stage. Line numbers refer to the current copies.

| Pipeline stage | Source |
|---|---|
| Pangenome graph build (`cactus-pangenome`) | [`notebooks/S05_pangenome_poolseq.qmd`](notebooks/S05_pangenome_poolseq.qmd) L23 |
| Distance index rebuild + promotion (`.dist` → `.dist.bak`, `.dist2` → `.dist`) | S05 L26–28 |
| DNA mapping (`vg giraffe`) | S05 L30–65; clean one-line version L208 |
| Read support (`vg pack -Q 5`) | S05 L70–94; L211 (shared: also S06 L54) |
| Variant calling (`vg call -z -a`) | S05 L99–101; L214 (shared: also S06 L56, without `-z`) |
| `bcftools index` → `merge` → `query` (DNA) | S05 L104–108; L217–221 |
| GTF contig rename | [`notebooks/S06_pangenome_rnaseq.qmd`](notebooks/S06_pangenome_rnaseq.qmd) L22 |
| Spliced graph (`vg rna`) | S06 L25 |
| xg / prune / GCSA indexes | S06 L28–34 |
| Snarls (`vg snarls`) | S06 L37 |
| Spliced dist index | S06 L40 |
| RNA mapping (`vg mpmap`) | S06 L44–51 |
| `bcftools index` → `merge` → `query` (RNA) | **not in S06** (it queries per sample); S05 L104–108 applied to the RNA VCFs by request |
| Long-read RNA-seq | [`notebooks/S07_pangenome_lrs.qmd`](notebooks/S07_pangenome_lrs.qmd) — empty (stub notebook); not implemented |
| Downstream R analysis (site filtering, grenedalf comparison, plots) | S05 L111–202 — exploratory, deliberately **not** pipelined |

## Design decisions and deviations

- **Both layouts in both assays.** The notebooks show paired-end DNA (S05) and single-end RNA (S06); the pipeline infers single-end vs paired-end per row for both assays from an empty `fastq_2`.
- **Deterministic merge order.** Both assays merge in samplesheet order via the carried row index instead of shell-glob order (see the shared variant-calling section).
- **Distance index promotion.** S05 L26–28 rebuilds the distance index after cactus and renames files on the fly. The pipeline encodes the same three commands, split by ownership: `CACTUS_PANGENOME` archives its index as `<outname>.dist.bak`, and `VG_INDEX_DIST_UPDATE` builds `<outname>.dist2` and promotes it to `<outname>.dist`. Only one process ever publishes a file named `<outname>.dist` (`publish_dir_mode` is `copy`, so two publishers would race). The published layout has `<outname>.dist` (rebuilt, consumed by `vg giraffe -d`) and `<outname>.dist.bak` (archived cactus index); `<outname>.dist2` never reaches `--outdir`.
- **Explicit mapper format.** `vg mpmap` passes `-F GAM` explicitly; the notebook showed it inconsistently.
- **Graph artifacts as plain files.** The notebook ran cactus manually under apptainer and used the `.gz` artifacts; the pipeline `gunzip`s `<outname>.gfa.gz` → `.gfa` and `<outname>.vcf.gz` → `.vcf` so the published layout is uniform and predictable.
- **Directory names follow the data, not the file format.** Published output is grouped as `ref/` (reference graph, its indexes, the RNA spliced graph and the renamed GTF), `dna/` (DNA per-sample and combined VCF/TSV) and `rna/` (RNA per-sample and combined VCF/TSV). The notebooks wrote graph artifacts to `gfa/` and everything vg-mapped to `gam/`, but `gfa`/`gam` name file formats rather than content categories. The directory names are the literal `assay_dir` values (`'dna'`/`'rna'`) each workflow passes to the shared modules — layout choices, not params — and the smoke test asserts that the retired `gfa/` and `gam/` directories are never created.
- **Shared post-pack modules, one copy each.** From `vg pack` onwards the DNA and RNA steps differ only in the graph they read, the output subdirectory and two `vg call` details, so each step exists once under `modules/local/variation/`: `VG_PACK`, `VG_CALL`, `BCFTOOLS_INDEX`, `BCFTOOLS_MERGE`, `BCFTOOLS_QUERY`. `assay_dir` (`'dna'`/`'rna'`) selects the published subdirectory; `VG_CALL` also takes `graph_is_gbz` (the `-z` flag, GBZ input only) and `call_sample_name` (`''` → the sample id, RNA-seq may pass `--rna_call_sample`). `VG_GIRAFFE` publishes its `.gam.log` as a separate output so both assays hand the shared chain the same `(sample, file, idx)` tuple shape.
- **One directory, one meaning.** `modules/local/reference/` holds every reference-side step: GBZ construction and its distance index, the GTF contig rename, and the spliced graph with its xg/gcsa/snarls/dist indexes. `modules/local/dnaseq/` and `modules/local/rnaseq/` hold only the per-sample mappers (`VG_GIRAFFE`, `VG_MPMAP`), and `modules/local/variation/` holds the shared pack → query chain. The spliced-graph modules live in `reference/` but are invoked by `RNASEQ` rather than `REFERENCE`, because they need `--gtf`, which a DNA-only run must not require.
- **Single-sample assays skip index and merge.** `bcftools index` exists only to enable `bcftools merge`, and a single sample has nothing to merge, so both steps are skipped and the sample's VCF goes straight to `BCFTOOLS_QUERY` (which needs no index without a region). The consequences are deliberate: that assay publishes `<sample>.vcf.gz` and the assay-level `combined.vcf.tsv.gz` (the same columns as a merged run), but no `.csi` and no `combined.vcf.gz`; `combined_vcf` is an empty channel for such a run. The gating is channel-level in each workflow, so no task is even submitted and no param or module changes.
- **`ref_name` couples cactus and the GTF.** `params.ref_name` is both `cactus --reference` and the GTF rename prefix, so it must equal the reference sample name in the seqfile, or the renamed GTF contig names will not match GBZ haplotype paths.
- **Contig-agnostic GTF rename.** The notebook's `sed 's/^(NC_|NW_)/<ref>#0#\1/'` (S06 L22) only works for RefSeq-style accessions; reference assemblies routinely use `chr*`, `scaffold*` or other names, so the pipeline prefixes the first tab-delimited field of every non-comment line instead. Comment/blank lines and lines without a tab (not GTF) are passed through unchanged rather than guessed at, and an already-renamed GTF must not be fed in again.
- **No replicate and no per-sample reference.** Derived from what the notebooks actually require: every sample maps to the one shared graph, and pooled samples are a single sample (pooling happens upstream).
- **Threshold exposure.** The hard-coded `-Q 5` pack filter became `params.min_mapq` (default 5, i.e. the notebook value); thread counts, cactus cores, and `gcsa_tmpdir` are likewise params.
- **Pass-through flags.** `--permissiveContigFilter`, `--haplo`, `--chrom-vg clip filter`, `--chrom-og full`, and `--viz` are passed to cactus exactly as in the notebook.
- **`RENAME_GTF` is container-less.** The rename is a single contig-agnostic `sed` substitution; it uses the helper in `bin/` (available on `PATH` in the Nextflow script environment) rather than pulling a container.
- **`vg call -z` is DNA-only.** S06 L56 copies `-z` from the DNA command, but `-z` restricts calling to the GBZ haplotypes and `vg call` rejects it for a non-GBZ graph, so `VG_CALL` passes it only when the graph is the GBZ (`graph_is_gbz`); the RNA (`-a` only) and DNA (`-z -a`) flag sets are otherwise unchanged from the notebooks.
- **RNA cross-sample merge (requested extension).** S06 stops at a per-sample `bcftools query`; at the user's request the RNA VCFs are now indexed, merged and queried per assay exactly like the DNA ones, using the shared modules, so the pipeline produces `rna/combined.vcf.gz` and `rna/combined.vcf.tsv.gz` instead of per-sample TSVs. There is no notebook source for this step (see provenance).

## Known limitations

- `vg call | bgzip -c` requires `bgzip` inside the vg image. vgteam images ship htslib tools; verify at the first real run or override `--vg_container`.
- `--viz` is passed to `cactus-pangenome` as is; the cactus image must provide its dependencies (graphviz/odgi).
- `bcftools_container` uses a placeholder tag (`1.19--h3ea31c5_0`). Pin it to a tag available on your cluster; no particular version is required.
- Cactus resources are estimates (64 GB RAM / 48 h per process) and must be tuned to the assemblies; all other processes default to 1 CPU, 8 GB, 8 h, with per-process overrides in `nextflow.config`.
- `publish_dir_mode` defaults to `copy`, which will copy large GAM/pack files. `symlink` is recommended for real runs.
- `CACTUS_PANGENOME` deletes `<outdir>/toil_work` only when it succeeds. A failed or interrupted cactus run leaves the Toil jobstore behind, and rerunning into the same `--outdir` will hand Toil an existing jobstore; remove `<outdir>/toil_work` manually before such a rerun.
- `RENAME_GTF` runs on the host with `sed`, so `bin/` must be reachable on `PATH`. It assumes a tab-delimited GTF: lines without a tab are passed through unprefixed rather than guessed at.
- RNA and DNA are merged separately, one `combined.vcf.gz`/`combined.vcf.tsv.gz` per assay: there is no cross-assay merge, and cross-assay comparison stays downstream.
- `--rna_call_sample` sets a single fixed `vg call -s` for every RNA sample, so with more than one RNA sample the per-sample VCFs share a sample column and `bcftools merge` fails on duplicate sample names. Leave it unset (the default: sample id) unless the run has a single RNA sample.
- No biological validation is automated in this repository: the smoke test is stub-only, and real-data correctness is established by the notebook comparisons (e.g. S05's grenedalf cross-check), not by CI.
- The `test` profile's `withName` overrides exist because profile-level params do not propagate into the base `process` block; keep them in sync when adding processes.

## Smoke test

```bash
bash test/run_smoke.sh
```

Requirements: the vendored Nextflow launcher at `./nextflow` and Nextflow ≥ 24.10. Everything runs with `-stub-run` — process stubs only, no containers, **no biological results are produced**; this is a connectivity/syntax test.

Inputs are tiny synthetic files with neutral names (`ref_a`/`ref_b`, `dna_a`/`dna_b`/`dna_hybrid`/`dna_c`, `rna_a`/`rna_b`) covering all four assay/layout combinations. The script rewrites `test/assemblies.txt` with absolute paths (the content is not read under `-stub-run`) and removes `work/` and `test/results/` before and after the run (`trap clean EXIT`).

Steps:

1. **bin script unit check** — `bash -n bin/rename_gtf_for_vg.sh`, then run it on `test/data/ref_a.gtf` and `diff` the result against `test/data/expected_ref_a.gtf` (regression guard for the `NC_`/`NW_` case → `ref_a#0#`), plus inline cases asserting that any contig name (`scaffold_b1`, `chr1`, `CM012345.1`, `1`) is prefixed, that comment/blank/non-tab lines pass through unchanged, and that sed-special characters in the sample name are emitted literally.
2. **Config parse** — `./nextflow config . -profile test`.
3. **`--run dnaseq`** — asserts `ref/test.gbz`, `.dist`, `.dist.bak`, `.shortread.withzip.min`, `.shortread.zipcodes`, `.snarls`, `.log`, `dna/combined.vcf.gz`, `dna/combined.vcf.tsv.gz`, and `dna/<sample>.{gam,gam.log,pack,vcf.gz,vcf.gz.csi}` for `dna_a`, `dna_b`, `dna_hybrid`, `dna_c`; also asserts that no `*.dist2` is published (intermediate), that `results/rna` is **not** created, and that the retired `results/gfa` and `results/gam` directories do not exist.
4. **`--run rnaseq`** — asserts `ref/ref_a.gtf`, `ref/test_spliced.{pg,xg,gcsa,gcsa.lcp,snarls,dist}`, `rna/combined.vcf.gz`, `rna/combined.vcf.tsv.gz` and `rna/<sample>.{gam,pack,vcf.gz,vcf.gz.csi}` for `rna_a`, `rna_b`; asserts no per-sample `rna/<sample>.vcf.tsv.gz` is published (retired), no `*pruned.pg*` is published (ephemeral intermediate), that `results/dna` is **not** created, and that the retired `results/gfa` and `results/gam` directories do not exist.
5. **Default entry (`both`)** — one run with `-with-dag`, asserting both DNA and RNA outputs exist, the DAG file is non-empty, and `CACTUS_PANGENOME` appears exactly once in the DAG (single shared reference).
6. **Single-sample assays** — one `--run dnaseq` and one `--run rnaseq` with a samplesheet holding exactly one row per assay: asserts the run succeeds, publishes `<assay_dir>/<sample>.vcf.gz` and `<assay_dir>/combined.vcf.tsv.gz`, publishes **no** `<sample>.vcf.gz.csi` and **no** `combined.vcf.gz`, and does not submit `BCFTOOLS_INDEX`/`BCFTOOLS_MERGE` (checked in the captured run log).
7. **Negative tests** — samplesheets with a bogus assay and with an empty `fastq_1` must fail, and the error messages must contain `Unknown assay` / `fastq_1 missing`.

Success ends with `ALL SMOKE TESTS PASSED`.

## Real-data runs

See the Usage section of [`README.md`](README.md) for command lines. Additional practical notes:

- **Apptainer.** Pass local SIFs via `--vg_sif`, `--cactus_sif`, `--bcftools_sif`; they override the registry images. `apptainer.autoMounts = true` is set, so bind-mount data directories the same way you would with `apptainer exec`, e.g. `-c 'apptainer.runOptions = "--bind /nfs_data"'`.
- **Resources.** Raise `cactus_cons_cores` / `cactus_index_cores` / `cactus_mg_cores` / `cactus_map_cores` and the `CACTUS_PANGENOME` memory/time to match the assemblies.
- **Publishing.** Use `--publish_dir_mode symlink` when GAM and pack files are large.
- **Naming.** `--ref_name` must match the reference name in the seqfile; `--outname` prefixes every graph artifact.

## Open items

- Long-read RNA-seq (S07) is an empty notebook stub and has no implementation.
