// DNA-seq (pool-seq) workflow: per sample giraffe -> pack -> call, then
// bcftools index -> merge -> query.
include { VG_GIRAFFE }        from '../modules/local/dnaseq/vg_giraffe'
include { VG_PACK_DNA }       from '../modules/local/dnaseq/vg_pack_dna'
include { VG_CALL_DNA }       from '../modules/local/dnaseq/vg_call_dna'
include { BCFTOOLS_INDEX }    from '../modules/local/dnaseq/bcftools_index'
include { BCFTOOLS_MERGE_DNA } from '../modules/local/dnaseq/bcftools_merge_dna'
include { BCFTOOLS_QUERY_DNA } from '../modules/local/dnaseq/bcftools_query_dna'

workflow DNASEQ {
    take:
    gbz         // path: <outname>.gbz
    dist        // path: <outname>.dist (rebuilt by VG_INDEX_DIST2; giraffe -d)
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
    ch_pack = VG_PACK_DNA(ch_gam, gbz)
    ch_vcf  = VG_CALL_DNA(ch_pack, gbz, snarls)
    ch_csi  = BCFTOOLS_INDEX(ch_vcf)

    // Merge order = samplesheet order; the carried row index makes the order
    // deterministic regardless of task completion order. The sorted list is
    // reshaped into a single tuple of lists, which is the input shape
    // BCFTOOLS_MERGE_DNA declares.
    ch_merged = ch_csi.toSortedList { a, b -> a[3] <=> b[3] }.map { items -> [items.collect { it[0] }, items.collect { it[1] }, items.collect { it[2] }, items.collect { it[3] }] } | BCFTOOLS_MERGE_DNA

    ch_tsv = BCFTOOLS_QUERY_DNA(ch_merged.combined_vcf)

    emit:
    combined_vcf = ch_merged.combined_vcf
    tsv          = ch_tsv.tsv
}
