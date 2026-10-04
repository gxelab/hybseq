// Prune the spliced graph before GCSA indexing.
// EPHEMERAL OUTPUT: the pruned graph is only an intermediate for the GCSA build.
// The output has no publishDir, so it stays in the work dir and is never published.
process VG_PRUNE {
    tag "$params.outname"

    container params.vg_container

    input:
    path pg

    output:
    path "${params.outname}_spliced.pruned.pg", emit: pruned_pg

    script:
    """
    vg prune ${pg} > ${params.outname}_spliced.pruned.pg
    """

    stub:
    """
    touch ${params.outname}_spliced.pruned.pg
    """
}
