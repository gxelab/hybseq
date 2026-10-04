// Build the GCSA2 index from the pruned spliced graph (`vg index -t -g -b`).
// vg index writes both the .gcsa and the adjacent .gcsa.lcp; both are emitted
// because vg mpmap needs them side by side.
process VG_INDEX_GCSA {
    tag "$params.outname"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path pruned_pg

    output:
    path "ref/${params.outname}_spliced.gcsa",     emit: gcsa
    path "ref/${params.outname}_spliced.gcsa.lcp", emit: gcsa_lcp

    script:
    """
    mkdir -p ref ${params.gcsa_tmpdir}
    vg index -t ${task.cpus} -g ref/${params.outname}_spliced.gcsa -b ${params.gcsa_tmpdir} ${pruned_pg}
    """

    stub:
    """
    mkdir -p ref
    touch ref/${params.outname}_spliced.gcsa ref/${params.outname}_spliced.gcsa.lcp
    """
}
