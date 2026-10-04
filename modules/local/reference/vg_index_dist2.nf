// Rebuild the distance index from the GBZ and promote it to the primary
// <outname>.dist, replacing the cactus-pangenome index that CACTUS_PANGENOME
// archived as <outname>.dist.bak. Mirrors the notebook (S05 L26-28):
//   vg index -j ref/<outname>.dist2 ref/<outname>.gbz
//   mv ref/<outname>.dist ref/<outname>.dist.bak
//   mv ref/<outname>.dist2 ref/<outname>.dist
// <outname>.dist2 is an intermediate and is never published.
process VG_INDEX_DIST2 {
    tag "$params.outname"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path gbz

    output:
    path "ref/${params.outname}.dist", emit: dist

    script:
    """
    mkdir -p ref
    vg index -j ref/${params.outname}.dist2 ${gbz}
    mv ref/${params.outname}.dist2 ref/${params.outname}.dist
    """

    stub:
    """
    mkdir -p ref
    touch ref/${params.outname}.dist2
    mv ref/${params.outname}.dist2 ref/${params.outname}.dist
    """
}
