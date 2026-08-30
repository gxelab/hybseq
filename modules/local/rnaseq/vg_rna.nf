// Build the spliced pangenome graph (PackedGraph) from the GBZ and the renamed
// GTF with vg rna.
process VG_RNA {
    tag "$params.outname"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path gbz
    path gtf

    output:
    path "gfa/${params.outname}_spliced.pg", emit: pg

    script:
    """
    mkdir -p gfa
    vg rna -p --threads ${task.cpus} --transcripts ${gtf} --use-hap-ref --gbz-format ${gbz} \\
        > gfa/${params.outname}_spliced.pg
    """

    stub:
    """
    mkdir -p gfa
    touch gfa/${params.outname}_spliced.pg
    """
}
