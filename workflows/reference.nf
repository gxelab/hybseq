// Reference graph construction and indexing, shared by the DNA-seq and RNA-seq workflows.
include { CACTUS_PANGENOME }     from '../modules/local/reference/cactus_pangenome'
include { VG_INDEX_DIST_UPDATE } from '../modules/local/reference/vg_index_dist_update'

workflow REFERENCE {
    take:
    assemblies               // path: Minigraph-Cactus seqfile

    main:
    ch_cactus      = CACTUS_PANGENOME(assemblies)
    ch_dist_update = VG_INDEX_DIST_UPDATE(ch_cactus.gbz)

    emit:
    gbz      = ch_cactus.gbz
    dist     = ch_dist_update.dist
    min      = ch_cactus.min
    zipcodes = ch_cactus.zipcodes
    snarls   = ch_cactus.snarls
}
