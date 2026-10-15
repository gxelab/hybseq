# AGENTS.md

Conventions and rules for AI coding agents working in this repository. See [`README.md`](README.md) for user documentation and [`development.md`](development.md) for design notes. Keep this file short — only rules that matter in every session.

## Project and layout

- `hybseq` is a Nextflow DSL2 pipeline (Minigraph-Cactus + vg) for pangenome-graph analysis of interspecies hybrid pooled DNA-seq and RNA-seq. Entry point `main.nf`, workflows in `workflows/`, one process per file under `modules/local/<reference|variation|dnaseq|rnaseq>/`. Keep the split by meaning: `reference/` = reference-graph construction and indexing (GBZ, distance index, GTF rename, spliced graph + its indexes), `dnaseq/`/`rnaseq/` = per-sample mapping only, `variation/` = the steps shared by both assays from `vg pack` onwards (never copy a shared step into `dnaseq/` or `rnaseq/`).
- Run everything through the vendored `./nextflow` launcher (pins Nextflow 26.04.6); never assume a system-installed `nextflow`.
- `workflows/reference.nf` builds the graph once, then `workflows/dnaseq.nf` and `workflows/rnaseq.nf` consume it. Keep them separate; never fold the two assays into one workflow.
- `bin/rename_gtf_for_vg.sh` is a pipeline helper called by the `RENAME_GTF` process; keep it in `bin/` and keep `test/data/expected_ref_a.gtf` in sync if you change it.

## Documentation ownership

- `README.md` is user-facing only: intro, usage, inputs, samplesheet, configuration, outputs, AI-assistance declaration, license. Keep architecture, rationale, limitations, and test notes out of it.
- `development.md` holds workflow internals, provenance, design decisions, limitations, and smoke-test instructions; `AGENTS.md` holds session rules. Update all affected files when behavior changes.

## Configuration

- Never hardcode paths or cluster specifics; every knob is a documented `params` entry in `nextflow.config` with a safe default.
- Thread counts are used both in the tool flags (`-t ${task.cpus}`) and in the `process { withName: ... }` block. Update both, including the `test` profile overrides.
- Do not add a macOS/Apple-silicon execution profile: cactus needs an AVX2-capable x86-64 CPU, so macOS hosts are for `-stub-run` smoke tests only. There is no conda path either — do not re-add `conda` directives or `*_conda_env` params.

## Contracts

- The samplesheet schema is fixed at `sample,assay,fastq_1,fastq_2`; library layout (single-end vs paired-end) is inferred from an empty `fastq_2`. Do not add columns (replicate, per-sample reference) without an explicit request.
- Every process keeps a `stub:` block and a `container` directive. There is no conda fallback: the `conda` profile and the `*_conda_env` params were removed. `RENAME_GTF` is the only process without a `container` directive (host `sed` via `bin/rename_gtf_for_vg.sh`).
- Publish outputs as work-dir-relative paths via `publishDir "${params.outdir}"`. Give ephemeral intermediates (e.g. `VG_PRUNE`'s `*_spliced.pruned.pg`) no `publishDir` at all.
- Dynamic directive values must be closures: a `${...}` string inside a `publishDir` attribute is resolved when the process is defined, so `BCFTOOLS_INDEX` (which needs its `assay_dir` input in the publish path) uses `path: { ... }`.

## Testing and safety

- Verify changes with `bash test/run_smoke.sh` (synthetic inputs, `-profile test -stub-run`, no containers). It checks connectivity, never biology.
- Never launch a real cactus/vg analysis, pull containers, or run without `-stub-run`; the default `--outdir` is `${launchDir}/results` and can be huge.
- Do not modify files outside this repository, and do not delete fixtures under `test/data/` or helper scripts under `bin/`.
- Do not commit or push; leave the working tree for the user to review.

## Style

- Docs and comments in English, matching the existing tone. Use tables for params and outputs, and keep diffs minimal and surgical.
