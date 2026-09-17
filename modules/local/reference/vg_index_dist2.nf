// Rebuild the distance index from the GBZ and promote it to the primary
// <outname>.dist, replacing the cactus-pangenome index that CACTUS_PANGENOME
// archived as <outname>.dist.bak. Mirrors the notebook (S05 L26-28):
//   vg index -j gfa/<outname>.dist2 gfa/<outname>.gbz
//   mv gfa/<outname>.dist gfa/<outname>.dist.bak
//   mv gfa/<outname>.dist2 gfa/<outname>.dist
// <outname>.dist2 is an intermediate and is never published.
process VG_INDEX_DIST2 {
    tag "$params.outname"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path gbz

    output:
    path "gfa/${params.outname}.dist", emit: dist

    script:
    """
    mkdir -p gfa
    vg index -j gfa/${params.outname}.dist2 ${gbz}
    mv gfa/${params.outname}.dist2 gfa/${params.outname}.dist
    """

    stub:
    """
    mkdir -p gfa
    touch gfa/${params.outname}.dist2
    mv gfa/${params.outname}.dist2 gfa/${params.outname}.dist
    """
}
