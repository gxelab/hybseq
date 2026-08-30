#!/usr/bin/env nextflow
// hybseq: pangenome-graph processing of interspecies hybrid genomic data
// (pooled DNA-seq and RNA-seq).
//
// DSL2 pipeline: Minigraph-Cactus pangenome graph construction and indexing,
// vg mapping / read support / graph-based genotyping, and bcftools tabulation.
// Separate DNA-seq and RNA-seq workflows share one reference graph.
nextflow.enable.dsl = 2

include { REFERENCE } from './workflows/reference'
include { DNASEQ }    from './workflows/dnaseq'
include { RNASEQ }    from './workflows/rnaseq'

// Parse and validate the samplesheet.
// Schema: sample, assay, fastq_1, fastq_2.
// Library layout (single-end vs paired-end) is inferred per row for BOTH assays:
// fastq_2 empty -> single-end, fastq_2 non-empty -> paired-end.
// No replicate column and no per-sample reference column: every sample maps to
// one shared graph, and pooled samples are simply one sample (pooling happens
// upstream of the pipeline).
def makeSamplesChannel() {
    def idx = 0
    Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def sample = row.sample
            def assay  = row.assay
            def fq1    = row.fastq_1
            // Missing trailing columns are padded with empty strings; strip any
            // surrounding quotes so `""` is also read as empty.
            def fq2    = row.fastq_2 ? row.fastq_2.toString().replaceAll(/^"+|"+$/, '') : ''
            if (!fq1) { error "fastq_1 missing for sample ${sample} in ${params.samplesheet}" }
            if (!(assay in ['dnaseq', 'rnaseq'])) {
                error "Unknown assay '${assay}' for sample ${sample} in ${params.samplesheet} (expected dnaseq|rnaseq)"
            }
            def row_idx = idx;
            idx = idx + 1;
            [sample, assay, file(fq1, checkIfExists: true), fq2 ? file(fq2, checkIfExists: true) : '', row_idx]
        }
}

workflow {
    // Entry workflow: build the reference graph once, then run the DNA-seq and/or
    // RNA-seq workflows against it, selected by params.run (both|dnaseq|rnaseq).
    main:
    if (!(params.run in ['both', 'dnaseq', 'rnaseq'])) {
        error "Invalid --run '${params.run}' (expected both|dnaseq|rnaseq)"
    }
    ch_samples = makeSamplesChannel()
    ref = REFERENCE(file(params.assemblies, checkIfExists: true))

    if (params.run in ['both', 'dnaseq']) {
        DNASEQ(ref.gbz, ref.dist, ref.min, ref.zipcodes, ref.snarls, ch_samples)
    }
    if (params.run in ['both', 'rnaseq']) {
        if (!params.gtf) { error 'params.gtf is required to run the RNA-seq workflow' }
        RNASEQ(ref.gbz, Channel.value(file(params.gtf, checkIfExists: true)), ch_samples)
    }
}
