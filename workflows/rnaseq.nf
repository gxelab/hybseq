// RNA-seq workflow: renamed GTF -> spliced graph + indexes, then per sample
// mpmap -> pack -> call, then bcftools index -> merge -> query. Everything from
// pack onwards is shared with the DNA-seq workflow (modules/local/variation/).
include { RENAME_GTF }        from '../modules/local/reference/rename_gtf'
include { VG_RNA }            from '../modules/local/reference/vg_rna'
include { VG_INDEX_XG }       from '../modules/local/reference/vg_index_xg'
include { VG_PRUNE }          from '../modules/local/reference/vg_prune'
include { VG_INDEX_GCSA }     from '../modules/local/reference/vg_index_gcsa'
include { VG_SNARLS }         from '../modules/local/reference/vg_snarls'
include { VG_INDEX_DIST_RNA } from '../modules/local/reference/vg_index_dist_rna'
include { VG_MPMAP }          from '../modules/local/rnaseq/vg_mpmap'
include { VG_PACK }           from '../modules/local/variation/vg_pack'
include { VG_CALL }           from '../modules/local/variation/vg_call'
include { BCFTOOLS_INDEX }    from '../modules/local/variation/bcftools_index'
include { BCFTOOLS_MERGE }    from '../modules/local/variation/bcftools_merge'
include { BCFTOOLS_QUERY }    from '../modules/local/variation/bcftools_query'

workflow RNASEQ {
    take:
    gbz         // path: <outname>.gbz (shared with the DNA-seq workflow)
    gtf         // path: input GTF
    samples     // tuple [sample, assay, fq1, fq2, idx]

    main:
    ch_rna = samples
        .filter { it[1] == 'rnaseq' }
        .map { [it[0], it[2], it[3], it[4]] }
        .ifEmpty { error 'No rnaseq samples found in the samplesheet' }

    ch_gtf = RENAME_GTF(gtf)

    ch_pg      = VG_RNA(gbz, ch_gtf.renamed_gtf)
    ch_xg      = VG_INDEX_XG(ch_pg.pg)
    ch_gcsa    = VG_INDEX_GCSA(VG_PRUNE(ch_pg.pg))
    ch_snarls  = VG_SNARLS(ch_pg.pg)
    ch_dist    = VG_INDEX_DIST_RNA(ch_xg.xg)

    ch_gam  = VG_MPMAP(ch_rna, ch_xg.xg, ch_gcsa.gcsa, ch_gcsa.gcsa_lcp, ch_dist.dist)
    ch_pack = VG_PACK(ch_gam.gam, ch_xg.xg, 'rna')
    // The spliced xg is not a GBZ, so vg call gets no -z (that flag only applies to
    // GBZ input); --rna_call_sample overrides vg call -s, otherwise the sample id.
    ch_vcf  = VG_CALL(ch_pack, ch_xg.xg, ch_snarls.snarls, 'rna', false, params.rna_call_sample ?: '')

    // Index and merge only when the assay has more than one sample, exactly as in the
    // DNA-seq workflow: a single sample has nothing to merge, so BCFTOOLS_QUERY reads
    // its VCF directly.
    ch_sorted = ch_vcf.toSortedList { a, b -> a[3] <=> b[3] }   // one emission: every sample, samplesheet order
    ch_multi  = ch_sorted.filter { it.size() > 1 }
    ch_single = ch_sorted.filter { it.size() == 1 }

    // flatMap (not flatten, which recurses into the tuples) turns the collected list
    // back into one item per sample; with a single sample the channel stays empty.
    ch_csi = BCFTOOLS_INDEX(ch_multi.flatMap { items -> items }, 'rna')
    // Sorted (sample, vcf, csi, idx) tuples -> the single tuple of lists BCFTOOLS_MERGE
    // declares; an empty collection (single sample) is dropped before the reshape.
    ch_indexed  = ch_csi.toSortedList { a, b -> a[3] <=> b[3] }.filter { !it.isEmpty() }
    ch_merge_in = ch_indexed.map { items -> [items.collect { it[0] }, items.collect { it[1] }, items.collect { it[2] }, items.collect { it[3] }] }
    ch_merged   = BCFTOOLS_MERGE(ch_merge_in, 'rna')

    // Merged VCF when there was something to merge, the sample VCF when there was not.
    ch_query_in = ch_merged.combined_vcf.mix(ch_single.map { items -> items[0][1] })
    ch_tsv      = BCFTOOLS_QUERY(ch_query_in, 'rna')

    emit:
    combined_vcf = ch_merged.combined_vcf
    tsv          = ch_tsv.tsv
}
