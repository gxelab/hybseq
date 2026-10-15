// Build the distance index of the spliced xg graph for vg mpmap.
process VG_INDEX_DIST_RNA {
    tag "$params.outname"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path xg

    output:
    path "ref/${params.outname}_spliced.dist", emit: dist

    script:
    """
    mkdir -p ref
    vg index -j ref/${params.outname}_spliced.dist ${xg}
    """

    stub:
    """
    mkdir -p ref
    touch ref/${params.outname}_spliced.dist
    """
}
