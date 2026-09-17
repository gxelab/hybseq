# AGENTS.md

Conventions and rules for AI coding agents working in this repository. See [`README.md`](README.md) for user documentation and [`development.md`](development.md) for design notes. Keep this file short — only rules that matter in every session.

## Project and layout

- `hybseq` is a Nextflow DSL2 pipeline (Minigraph-Cactus + vg) for pangenome-graph analysis of interspecies hybrid pooled DNA-seq and RNA-seq. Entry point `main.nf`, workflows in `workflows/`, one process per file under `modules/local/<reference|dnaseq|rnaseq>/`.
- Run everything through the vendored `./nextflow` launcher (Nextflow ≥ 24.10); never assume a system-installed `nextflow`.
- `workflows/reference.nf` builds the graph once, then `workflows/dnaseq.nf` and `workflows/rnaseq.nf` consume it. Keep them separate; never fold the two assays into one workflow.
- `bin/rename_gtf_for_vg.sh` is a pipeline helper called by the `RENAME_GTF` process; keep it in `bin/` and keep `test/data/expected_ref_a.gtf` in sync if you change it.

## Source of truth

- `notebooks/S05_pangenome_poolseq.qmd` and `notebooks/S06_pangenome_rnaseq.qmd` record the real in-house workflows and are the source of truth. Translate them faithfully; never edit the notebooks.
- Do not modernize tools, flags, thresholds, or versions (vg 1.73, cactus 3.1.4) unless explicitly asked; `S07_pangenome_lrs.qmd` is empty and long-read support is not implemented.

## Documentation ownership

- `README.md` is user-facing only: intro, usage, inputs, samplesheet, configuration, outputs, AI-assistance declaration, license. Keep architecture, rationale, limitations, and test notes out of it.
- `development.md` holds workflow internals, provenance, design decisions, limitations, and smoke-test instructions; `AGENTS.md` holds session rules. Update all affected files when behavior changes.

## Configuration

- Never hardcode paths or cluster specifics; every knob is a documented `params` entry in `nextflow.config` with a safe default.
- Thread counts are used both in the tool flags (`-t ${task.cpus}`) and in the `process { withName: ... }` block. Update both, including the `test` profile overrides.

## Contracts

- The samplesheet schema is fixed at `sample,assay,fastq_1,fastq_2`; library layout (single-end vs paired-end) is inferred from an empty `fastq_2`. Do not add columns (replicate, per-sample reference) without an explicit request.
- Every process keeps a `stub:` block and a `container` directive so `-stub-run` works. `RENAME_GTF` is the only container-less process (host `sed` via `bin/rename_gtf_for_vg.sh`).
- Publish outputs as work-dir-relative paths via `publishDir "${params.outdir}"`. Give ephemeral intermediates (e.g. `VG_PRUNE`'s `*_spliced.pruned.pg`) no `publishDir` at all.

## Testing and safety

- Verify changes with `bash test/run_smoke.sh` (synthetic inputs, `-profile test -stub-run`, no containers). It checks connectivity, never biology.
- Never launch a real cactus/vg analysis, pull containers, or run without `-stub-run`; the default `--outdir` is `${launchDir}/results` and can be huge.
- Do not modify files outside this repository, and do not delete fixtures under `test/data/` or helper scripts under `bin/`.
- Do not commit or push; leave the working tree for the user to review.

## Style

- Docs and comments in English, matching the existing tone. Use tables for params and outputs, and keep diffs minimal and surgical.
