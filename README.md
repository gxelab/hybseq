# hybseq

A Nextflow DSL2 pipeline for **pangenome-graph based processing of interspecies hybrid genomic data** (DNA-seq and RNA-seq) built with [Minigraph-Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) and [vg](https://github.com/vgteam/vg).

The pipeline builds one pangenome graph from a reference and one or more additional assemblies, maps short reads to it, computes read support, calls variants per sample, and tabulates the results. DNA-seq and RNA-seq are separate workflows that share the same reference graph.

## Workflows

### DNA-seq / pool-seq

1. Build the pangenome graph with **Minigraph-Cactus 3.1.4** (`cactus-pangenome`, GBZ/GFA/VCF outputs, giraffe mapping indexes)
2. Rebuild a distance index as `<outname>.dist2` with `vg index -j`
3. Map reads to the graph with **vg giraffe** (vg 1.73; single-end or paired-end)
4. Compute read support with **vg pack** (`-Q 5`)
5. Call variants per sample with **vg call** (`-z -a`)
6. Index, merge and tabulate all samples with **bcftools** (`index -f` → `merge` → `query`)

### RNA-seq

1. Rewrite `NC_`/`NW_` contig prefixes in the GTF to `<ref>#0#` so they match GBZ haplotype-path names (`bin/rename_gtf_for_vg.sh`)
2. Build the spliced pangenome graph with **vg rna** (PackedGraph, `--use-hap-ref --gbz-format`)
3. Index: `vg index -x` (xg) → `vg prune` → `vg index -g` (GCSA) → `vg snarls` → `vg index -j` (dist)
4. Map reads with **vg mpmap** (`-n RNA -l short`; single-end or paired-end)
5. Read support (**vg pack**, `-Q 5`) and variant calling (**vg call**, `-z -a`) on the spliced graph
6. Per-sample tabulation with **bcftools query** (no merge is performed for RNA-seq)

## DNA-seq vs RNA-seq

| | DNA-seq | RNA-seq |
|---|---|---|
| Graph mapped against | GBZ (`<outname>.gbz`) | Spliced PackedGraph built by `vg rna` from the same GBZ |
| Mapper | `vg giraffe -Z gbz -m min -z zipcodes -d dist`, one `-f` per read file | `vg mpmap -x xg -g gcsa -d dist -n RNA -l short`, one `-f` per read file |
| Read support / variant call inputs | `-x <outname>.gbz`, cactus snarls | `-x <outname>_spliced.xg`, `<outname>_spliced.snarls` |
| VCF aggregation | `bcftools index` + **`bcftools merge`** of all samples, then one `query` | **per-sample** `bcftools query`, no merge |
| Reads per sample | single-end or paired-end (inferred from the samplesheet) | single-end or paired-end (inferred from the samplesheet) |

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
| `ref_name` | `ref` | cactus `--reference`; GTF rename prefix; must match the seqfile reference name |
| `outname` | `pangenome` | cactus `--outName`; prefixes all graph artifacts |
| `giraffe_threads` / `pack_threads` | 16 | `vg giraffe` / `vg pack` `-t` |
| `call_threads` | 4 | `vg call` `-t` |
| `mpmap_threads` | 8 | `vg mpmap` `-t` |
| `rna_threads` | 8 | `vg rna --threads` |
| `gcsa_threads` | 4 | `vg index -t` (GCSA build) |
| `cactus_cons/index/mg_cores`, `cactus_map_cores` | 16 / 16 / 16 / 4 | `--consCores/--indexCores/--mgCores/--mapCores` |
| `min_mapq` | 5 | `vg pack -Q` (ignore reads below this MAPQ) |
| `gcsa_tmpdir` | `tmp` | `vg index -b` temporary directory |
| `rna_call_sample` | `null` | `vg call -s` for RNA; null → sample id, or a fixed name |
| `run` | `both` | `both` \| `dnaseq` \| `rnaseq` |
| `outdir` / `publish_dir_mode` | `results` / `copy` | outputs; `symlink` recommended for large GAMs |

Container images (versions the pipeline is developed and tested against):

- `vg_container` = `quay.io/vgteam/vg:v1.73.0` (vg 1.73)
- `cactus_container` = `quay.io/biocontainers/cactus:3.1.4` (cactus 3.1.4)
- `bcftools_container` = `quay.io/biocontainers/bcftools:1.19--h3ea31c5_0` (placeholder tag — pin to whatever is available on your cluster; no particular version is required)

**Profiles**

- `docker` — run with Docker images
- `apptainer` — run with Apptainer; local SIF paths override the registry images:
  `-profile apptainer --vg_sif /path/to/vg.sif --cactus_sif /path/to/cactus.sif --bcftools_sif /path/to/bcftools.sif`
  Bind-mount your data directories the same way you would with `apptainer exec`, e.g. `-c 'apptainer.runOptions = "--bind /my/data"'`.
- `test` — tiny synthetic inputs for the smoke test (see [development.md](development.md#smoke-test))

Resource note: cactus RAM/time are config defaults (64 GB / 48 h) and should be adjusted to your assemblies. The cactus Toil jobstore (`js/`) stays in the Nextflow work directory and is never published.

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
| `gfa/` | graph + index artifacts: `<outname>.gbz/.dist/.dist2/.shortread.withzip.min/.shortread.zipcodes/.snarls/.log/.gfa/.vcf`, `<ref_name>.gtf`, `<outname>_spliced.pg/.xg/.gcsa/.gcsa.lcp/.snarls/.dist` |
| `gam/` | per DNA sample: `<sample>.gam`, `<sample>.gam.log`, `<sample>.pack`, `<sample>.vcf.gz`, `<sample>.vcf.gz.csi`; combined: `combined.vcf.gz`, `combined.vcf.tsv.gz` |
| `rna/` | per RNA sample: `<sample>.gam`, `<sample>.pack`, `<sample>.vcf.gz`, `<sample>.vcf.tsv.gz` |

`<outname>_spliced.pruned.pg` is an ephemeral intermediate and is never published. The TSVs carry one row per VCF record (CHROM, POS, ID, REF, ALT, QUAL, FILTER, per-sample GT, DP, AD{0}, AD{1}, GQ) and are suited as inputs to downstream tabular analyses.

## AI assistance

This pipeline and its documentation were developed with AI assistance. AI-generated content was reviewed and revised by the authors, who are responsible for the final content.

| | |
|---|---|
| Models | DeepSeek V4 Flash / Pro |
| Harness | deepseek-harness + Claude Code |


## License

This project is released under the MIT License. See [`LICENSE`](LICENSE) for the full text.

Copyright (c) 2026 RNA GxE lab @ LZU.
