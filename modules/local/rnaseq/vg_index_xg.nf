// Build the xg index of the spliced graph.
process VG_INDEX_XG {
    tag "$params.outname"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path pg

    output:
    path "gfa/${params.outname}_spliced.xg", emit: xg

    script:
    """
    vg index -x gfa/${params.outname}_spliced.xg ${pg}
    """

    stub:
    """
    mkdir -p gfa
    touch gfa/${params.outname}_spliced.xg
    """
}
