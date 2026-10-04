// Compute snarls for the spliced graph (not required for mpmap when the dist index
// is provided, but required for vg call).
process VG_SNARLS {
    tag "$params.outname"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path pg

    output:
    path "ref/${params.outname}_spliced.snarls", emit: snarls

    script:
    """
    mkdir -p ref
    vg snarls ${pg} > ref/${params.outname}_spliced.snarls
    """

    stub:
    """
    mkdir -p ref
    touch ref/${params.outname}_spliced.snarls
    """
}
