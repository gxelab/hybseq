// DNA-seq (pool-seq) workflow: per sample giraffe -> pack -> call, then
// bcftools index -> merge -> query. Everything from pack onwards is shared with
// the RNA-seq workflow (modules/local/variation/).
include { VG_GIRAFFE }     from '../modules/local/dnaseq/vg_giraffe'
include { VG_PACK }        from '../modules/local/variation/vg_pack'
include { VG_CALL }        from '../modules/local/variation/vg_call'
include { BCFTOOLS_INDEX } from '../modules/local/variation/bcftools_index'
include { BCFTOOLS_MERGE } from '../modules/local/variation/bcftools_merge'
include { BCFTOOLS_QUERY } from '../modules/local/variation/bcftools_query'

workflow DNASEQ {
    take:
    gbz         // path: <outname>.gbz
    dist        // path: <outname>.dist (rebuilt by VG_INDEX_DIST_UPDATE; giraffe -d)
    min         // path: <outname>.shortread.withzip.min
    zipcodes    // path: <outname>.shortread.zipcodes
    snarls      // path: <outname>.snarls
    samples     // tuple [sample, assay, fq1, fq2, idx]

    main:
    ch_dna = samples
        .filter { it[1] == 'dnaseq' }
        .map { [it[0], it[2], it[3], it[4]] }
        .ifEmpty { error 'No dnaseq samples found in the samplesheet' }

    ch_gam  = VG_GIRAFFE(ch_dna, gbz, min, zipcodes, dist)
    ch_pack = VG_PACK(ch_gam.gam, gbz, 'dna')
    // DNA calls against the GBZ: -z tells vg call to use the GBZ haplotypes; -s is
    // the sample id (no override for DNA).
    ch_vcf  = VG_CALL(ch_pack, gbz, snarls, 'dna', true, '')
    ch_csi  = BCFTOOLS_INDEX(ch_vcf, 'dna')

    // Merge order = samplesheet order; the carried row index makes the order
    // deterministic regardless of task completion order. The sorted list is reshaped
    // into a single tuple of lists, which is the input shape BCFTOOLS_MERGE declares.
    ch_merge_in = ch_csi.toSortedList { a, b -> a[3] <=> b[3] }.map { items -> [items.collect { it[0] }, items.collect { it[1] }, items.collect { it[2] }, items.collect { it[3] }] }
    ch_merged   = BCFTOOLS_MERGE(ch_merge_in, 'dna')

    ch_tsv = BCFTOOLS_QUERY(ch_merged.combined_vcf, 'dna')

    emit:
    combined_vcf = ch_merged.combined_vcf
    tsv          = ch_tsv.tsv
}
