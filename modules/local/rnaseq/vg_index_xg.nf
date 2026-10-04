// Build the xg index of the spliced graph.
process VG_INDEX_XG {
    tag "$params.outname"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path pg

    output:
    path "ref/${params.outname}_spliced.xg", emit: xg

    script:
    """
    mkdir -p ref
    vg index -x ref/${params.outname}_spliced.xg ${pg}
    """

    stub:
    """
    mkdir -p ref
    touch ref/${params.outname}_spliced.xg
    """
}
