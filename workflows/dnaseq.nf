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

    // Index and merge only when the assay has more than one sample: a single sample has
    // nothing to merge, so BCFTOOLS_QUERY reads its VCF directly (bcftools query needs
    // no index unless a region is requested). The samplesheet row index keeps the merge
    // order deterministic regardless of task completion order.
    ch_sorted = ch_vcf.toSortedList { a, b -> a[3] <=> b[3] }   // one emission: every sample, samplesheet order
    ch_multi  = ch_sorted.filter { it.size() > 1 }
    ch_single = ch_sorted.filter { it.size() == 1 }

    // flatMap (not flatten, which recurses into the tuples) turns the collected list
    // back into one item per sample; with a single sample the channel stays empty.
    ch_csi = BCFTOOLS_INDEX(ch_multi.flatMap { items -> items }, 'dna')
    // The sorted (sample, vcf, csi, idx) tuples are reshaped into the single tuple of
    // lists BCFTOOLS_MERGE declares; a collecting operator on an empty channel emits an
    // empty list, so the single-sample case is filtered out before the reshape.
    ch_indexed  = ch_csi.toSortedList { a, b -> a[3] <=> b[3] }.filter { !it.isEmpty() }
    ch_merge_in = ch_indexed.map { items -> [items.collect { it[0] }, items.collect { it[1] }, items.collect { it[2] }, items.collect { it[3] }] }
    ch_merged   = BCFTOOLS_MERGE(ch_merge_in, 'dna')

    // Merged VCF when there was something to merge, the sample VCF when there was not.
    ch_query_in = ch_merged.combined_vcf.mix(ch_single.map { items -> items[0][1] })
    ch_tsv      = BCFTOOLS_QUERY(ch_query_in, 'dna')

    emit:
    combined_vcf = ch_merged.combined_vcf
    tsv          = ch_tsv.tsv
}
