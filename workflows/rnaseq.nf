// RNA-seq workflow: renamed GTF -> spliced graph + indexes, then per sample
// mpmap -> pack -> call -> per-sample bcftools query.
include { RENAME_GTF }          from '../modules/local/rnaseq/rename_gtf'
include { VG_RNA }              from '../modules/local/rnaseq/vg_rna'
include { VG_INDEX_XG }         from '../modules/local/rnaseq/vg_index_xg'
include { VG_PRUNE }            from '../modules/local/rnaseq/vg_prune'
include { VG_INDEX_GCSA }       from '../modules/local/rnaseq/vg_index_gcsa'
include { VG_SNARLS }           from '../modules/local/rnaseq/vg_snarls'
include { VG_INDEX_DIST_RNA }   from '../modules/local/rnaseq/vg_index_dist_rna'
include { VG_MPMAP }            from '../modules/local/rnaseq/vg_mpmap'
include { VG_PACK_RNA }         from '../modules/local/rnaseq/vg_pack_rna'
include { VG_CALL_RNA }         from '../modules/local/rnaseq/vg_call_rna'
include { BCFTOOLS_QUERY_RNA }  from '../modules/local/rnaseq/bcftools_query_rna'

workflow RNASEQ {
    take:
    gbz         // path: <outname>.gbz (shared with the DNA-seq workflow)
    gtf         // path: input GTF
    samples     // tuple [sample, assay, fq1, fq2, idx]

    main:
    ch_rna = samples
        .filter { it[1] == 'rnaseq' }
        .map { [it[0], it[2], it[4]] }
        .ifEmpty { error 'No rnaseq samples found in the samplesheet' }

    ch_gtf = RENAME_GTF(gtf)

    ch_pg      = VG_RNA(gbz, ch_gtf.renamed_gtf)
    ch_xg      = VG_INDEX_XG(ch_pg.pg)
    ch_gcsa    = VG_INDEX_GCSA(VG_PRUNE(ch_pg.pg))
    ch_snarls  = VG_SNARLS(ch_pg.pg)
    ch_dist    = VG_INDEX_DIST_RNA(ch_xg.xg)

    ch_gam  = VG_MPMAP(ch_rna, ch_xg.xg, ch_gcsa.gcsa, ch_gcsa.gcsa_lcp, ch_dist.dist)
    ch_pack = VG_PACK_RNA(ch_gam, ch_xg.xg)
    ch_vcf  = VG_CALL_RNA(ch_pack, ch_xg.xg, ch_snarls.snarls)
    ch_tsv  = BCFTOOLS_QUERY_RNA(ch_vcf)

    emit:
    tsv = ch_tsv.tsv
}
