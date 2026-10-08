# hybseq  <img src="logo_light.svg" align="right" height="138" alt="HybSeq Logo" />

A Nextflow DSL2 pipeline for **pangenome-graph based processing of interspecies hybrid genomic data** (DNA-seq and RNA-seq) built with [Minigraph-Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) and [vg](https://github.com/vgteam/vg).

The pipeline builds a pangenome graph from a reference assembly and one or more additional assemblies. It maps short reads to this graph, using a spliced graph derived from it for RNA-seq, then computes read support, calls variants for each sample, and tabulates the results by assay.

DNA-seq and RNA-seq are implemented as separate workflows, but share the same reference graph and the same downstream steps for read-support calculation, variant calling, and result tabulation.

## Workflows

The pangenome graph is built once and reused by both assays; after mapping, DNA-seq and RNA-seq run the same read-support, variant-calling, and tabulation steps. `--run` selects which workflows execute (`both` by default).

### Pangenome graph construction

**Built once per run, for both assays** (reference workflow):

1. Build the pangenome graph with **Minigraph-Cactus 3.1.4** (`cactus-pangenome`, GBZ/GFA/VCF outputs, giraffe mapping indexes)
2. Rebuild the distance index after cactus: archive the cactus-produced `<outname>.dist` as `<outname>.dist.bak`, build the index with `vg index -j`, and promote the rebuilt index to `<outname>.dist`

**Needed for RNA-seq only** (skipped entirely with `--run dnaseq`; requires `--gtf`):

3. Prefix every GTF contig with `<ref>#0#` so the names match GBZ haplotype-path names — any contig naming scheme (`bin/rename_gtf_for_vg.sh`)
4. Build the spliced pangenome graph with **vg rna** (PackedGraph, `--use-hap-ref --gbz-format`)
5. Index the spliced graph: `vg index -x` (xg) → `vg prune` → `vg index -g` (GCSA) → `vg snarls` → `vg index -j` (dist)

A DNA-seq analysis uses steps 1–2 only; an RNA-seq analysis uses all five. Step 2 runs whichever assay is selected (it belongs to the shared reference workflow), but only DNA-seq consumes it: `vg giraffe -d` takes the rebuilt `<outname>.dist`, whereas `vg mpmap -d` takes the spliced `<outname>_spliced.dist` from step 5.

### Read mapping

One mapper run per sample, single-end or paired-end (layout inferred from the samplesheet, one `-f` per read file — see [Samplesheet format](#samplesheet-format)). Only the graph mapped against and the mapper differ between the assays:

| | DNA-seq | RNA-seq |
|---|---|---|
| Graph mapped against | GBZ (`<outname>.gbz`) | spliced PackedGraph built by `vg rna` (`<outname>_spliced.pg`) |
| Mapper | `vg giraffe -Z <outname>.gbz -m <outname>.shortread.withzip.min -z <outname>.shortread.zipcodes -d <outname>.dist` | `vg mpmap -x <outname>_spliced.xg -g <outname>_spliced.gcsa -d <outname>_spliced.dist -n RNA -l short` |
| Indexes needed | steps 1–2 (giraffe min/zipcodes + rebuilt distance index) | step 5 (spliced xg/gcsa/dist) |
| Reads per sample | single-end or paired-end | single-end or paired-end |

### Variants calling

Shared by both assays, per sample:

1. Compute read support with **vg pack** (`-Q 5`, i.e. `min_mapq`)
2. Call variants with **vg call** (`-a`, plus `-z` for DNA-seq only)
3. Index, merge and tabulate the assay's samples with **bcftools** (`index -f` → `merge` → `query`; skipped for a single sample, which is queried directly into `combined.vcf.tsv.gz`)

Graph-dependent differences:

| | DNA-seq | RNA-seq |
|---|---|---|
| Graph for `vg pack` / `vg call` | `<outname>.gbz` | `<outname>_spliced.xg` |
| Snarls | cactus `<outname>.snarls` | `<outname>_spliced.snarls` (step 5) |
| `vg call` sample and flags | `-z -a`, `-s` = sample id (`-z` uses the GBZ haplotypes) | `-a` only (`-z` is only valid for a GBZ); `-s` = sample id or `--rna_call_sample` |

Each assay aggregates to its own `combined.vcf.gz` / `combined.vcf.tsv.gz` under `dna/` or `rna/`; see [Outputs](#outputs) for the published files.

## Inputs

- **Assemblies seqfile** (`--assemblies`): Minigraph-Cactus seqfile, one `<name><TAB><fasta path>` per line; the reference must match `--ref_name`. See [`assets/assemblies.txt`](assets/assemblies.txt).
- **Samplesheet** (`--samplesheet`): see below. See [`assets/samplesheet.csv`](assets/samplesheet.csv) for an example.
- **GTF** (`--gtf`, required for `--run rnaseq|both`): transcript annotation of the reference assembly.

## Samplesheet format

CSV with header and four columns:

```csv
sample,assay,fastq_1,fastq_2
dna_a,dnaseq,/path/to/dna_a_R1.fq.gz,/path/to/dna_a_R2.fq.gz
dna_b,dnaseq,/path/to/dna_b.fq.gz
dna_hybrid,dnaseq,/path/to/dna_hybrid_R1.fq.gz,/path/to/dna_hybrid_R2.fq.gz
rna_a,rnaseq,/path/to/rna_a.fastq.gz
rna_b,rnaseq,/path/to/rna_b_R1.fq.gz,/path/to/rna_b_R2.fq.gz
```

- `sample`: sample id (used as `vg call -s` and as the output file prefix)
- `assay`: `dnaseq` or `rnaseq`
- `fastq_1` / `fastq_2`: read paths. **Library layout is inferred** for both assays: `fastq_2` empty → single-end, non-empty → paired-end. For single-end rows simply omit the trailing column (missing columns are padded with empty strings).

Validation: unknown assays and rows without `fastq_1` are rejected with explicit errors.

There is no replicate column and no per-sample reference column: every sample maps to the *same* graph. Pooled samples are simply one sample — pooling happens upstream of the pipeline.

## Configuration

All settings live in [`nextflow.config`](nextflow.config):

| param | default | meaning |
|---|---|---|
| `ref_name` | `ref` | cactus `--reference`; GTF rename prefix (applied to every contig name, whatever the naming scheme); must match the seqfile reference name |
| `outname` | `pangenome` | cactus `--outName`; prefixes all graph artifacts |
| `giraffe_threads` / `pack_threads` | 16 | `vg giraffe` / `vg pack` `-t` |
| `call_threads` | 4 | `vg call` `-t` |
| `mpmap_threads` | 8 | `vg mpmap` `-t` |
| `rna_threads` | 8 | `vg rna --threads` |
| `gcsa_threads` | 4 | `vg index -t` (GCSA build) |
| `cactus_cons/index/mg_cores`, `cactus_map_cores` | 16 / 16 / 16 / 4 | `--consCores/--indexCores/--mgCores/--mapCores` |
| `min_mapq` | 5 | `vg pack -Q` (ignore reads below this MAPQ) |
| `gcsa_tmpdir` | `tmp` | `vg index -b` temporary directory |
| `rna_call_sample` | `null` | `vg call -s` for RNA; null → sample id, or a fixed name (a fixed name is only valid for a single RNA sample: it makes the sample column identical in every RNA VCF, which `bcftools merge` rejects) |
| `run` | `both` | `both` \| `dnaseq` \| `rnaseq` |
| `outdir` / `publish_dir_mode` | `results` / `copy` | outputs; `symlink` recommended for large GAMs |

Container images (versions the pipeline is developed and tested against):

- `vg_container` = `quay.io/vgteam/vg:v1.73.0` (vg 1.73)
- `cactus_container` = `quay.io/biocontainers/cactus:3.1.4` (cactus 3.1.4)
- `bcftools_container` = `quay.io/biocontainers/bcftools:1.19--h3ea31c5_0` (placeholder tag — pin to whatever is available on your cluster; no particular version is required)

**Profiles**

- `docker` — run with Docker images.
- `apptainer` — run with Apptainer. Local SIF paths override the registry images:

  ```bash
  -profile apptainer \
    --vg_sif /path/to/vg.sif \
    --cactus_sif /path/to/cactus.sif \
    --bcftools_sif /path/to/bcftools.sif
  ```

  To bind-mount host directories into Apptainer containers, provide a Nextflow configuration file containing:

  ```groovy
  apptainer.runOptions = '--bind /my/data'
  ```

  and pass it with `-c`:

  ```bash
  nextflow run ... -profile apptainer -c apptainer.config
  ```
- `test` — tiny synthetic inputs for the smoke test (see [development.md](development.md#smoke-test))

Resource note: cactus RAM/time are config defaults (64 GB / 48 h) and should be adjusted to your assemblies. All cactus Toil scratch stays under `--outdir` (`toil_work/`, with the jobstore at `toil_work/js`) and is removed when the cactus task finishes; nothing in `toil_work/` is published.

## Usage

```bash
# Full pipeline (reference graph once, then both assays), Docker
./nextflow run . -profile docker --assemblies assemblies.txt --samplesheet samplesheet.csv --gtf annotation.gtf

# One assay only
./nextflow run . -profile docker --run dnaseq ...
./nextflow run . -profile docker --run rnaseq --gtf annotation.gtf ...

# Apptainer with local SIFs
./nextflow run . -profile apptainer \
    --vg_sif /path/to/vg.sif --cactus_sif /path/to/cactus.sif \
    --assemblies assemblies.txt --samplesheet samplesheet.csv --gtf annotation.gtf
```

## Outputs

Published under `--outdir` (`results/` by default):

| dir | contents |
|---|---|
| `ref/` | reference graph + index artifacts: `<outname>.gbz/.dist/.dist.bak/.shortread.withzip.min/.shortread.zipcodes/.snarls/.log/.gfa/.vcf`, the renamed `<ref_name>.gtf`, and the spliced-graph artifacts `<outname>_spliced.pg/.xg/.gcsa/.gcsa.lcp/.snarls/.dist` |
| `dna/` | per DNA sample: `<sample>.gam`, `<sample>.gam.log`, `<sample>.pack`, `<sample>.vcf.gz`; with ≥2 DNA samples also `<sample>.vcf.gz.csi` per sample and combined: `combined.vcf.gz`, `combined.vcf.tsv.gz`; with a single DNA sample: `combined.vcf.tsv.gz` only |
| `rna/` | per RNA sample: `<sample>.gam`, `<sample>.pack`, `<sample>.vcf.gz`; with ≥2 RNA samples also `<sample>.vcf.gz.csi` per sample and combined: `combined.vcf.gz`, `combined.vcf.tsv.gz`; with a single RNA sample: `combined.vcf.tsv.gz` only |

`<outname>_spliced.pruned.pg` and `<outname>.dist2` are intermediates and are never published; `<outname>.dist.bak` is the archived cactus distance index and `<outname>.dist` is the rebuilt index used for mapping. The TSVs carry one row per VCF record (CHROM, POS, ID, REF, ALT, QUAL, FILTER, per-sample GT, DP, AD{0}, AD{1}, GQ) and are suited as inputs to downstream tabular analyses.

`<sample>.vcf.gz.csi` and `combined.vcf.gz` are produced only when an assay has more than one sample: an index only enables merging, and a single sample has nothing to merge, so the sample VCF is queried directly and only the assay-level table (`combined.vcf.tsv.gz`) is published. The table's columns are the same either way.

## AI assistance

This pipeline and its documentation were developed with AI assistance. AI-generated content was reviewed and revised by the authors, who are responsible for the final content.

| | |
|---|---|
| Models | DeepSeek V4 Flash / Pro |
| Harness | deepseek-harness + Claude Code |


## Citation

Citation for hybseq and for every dependency it runs (Minigraph-Cactus, vg, giraffe, mpmap, bcftools, Nextflow, …) are collected in [`citations.md`](citations.md).

## License

This project is released under the MIT License. See [`LICENSE`](LICENSE) for the full text.

Copyright (c) 2026 RNA GxE lab @ LZU.
