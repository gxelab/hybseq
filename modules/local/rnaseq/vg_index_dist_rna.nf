// Build the distance index of the spliced xg graph for vg mpmap.
process VG_INDEX_DIST_RNA {
    tag "$params.outname"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path xg

    output:
    path "gfa/${params.outname}_spliced.dist", emit: dist

    script:
    """
    vg index -j gfa/${params.outname}_spliced.dist ${xg}
    """

    stub:
    """
    mkdir -p gfa
    touch gfa/${params.outname}_spliced.dist
    """
}
