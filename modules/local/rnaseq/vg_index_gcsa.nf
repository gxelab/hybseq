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
    path "gfa/${params.outname}_spliced.gcsa",     emit: gcsa
    path "gfa/${params.outname}_spliced.gcsa.lcp", emit: gcsa_lcp

    script:
    """
    mkdir -p gfa ${params.gcsa_tmpdir}
    vg index -t ${task.cpus} -g gfa/${params.outname}_spliced.gcsa -b ${params.gcsa_tmpdir} ${pruned_pg}
    """

    stub:
    """
    mkdir -p gfa
    touch gfa/${params.outname}_spliced.gcsa gfa/${params.outname}_spliced.gcsa.lcp
    """
}
