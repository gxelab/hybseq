// Rebuild a distance index from the GBZ as <outname>.dist2.
// NOTE: mapping uses the primary <outname>.dist index produced by cactus-pangenome
// itself; <outname>.dist2 is published but consumed by nothing downstream.
process VG_INDEX_DIST2 {
    tag "$params.outname"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path gbz

    output:
    path "gfa/${params.outname}.dist2", emit: dist2

    script:
    """
    mkdir -p gfa
    vg index -j gfa/${params.outname}.dist2 ${gbz}
    """

    stub:
    """
    mkdir -p gfa
    touch gfa/${params.outname}.dist2
    """
}
