// Rewrite NC_/NW_ contig prefixes in the GTF to <ref_name>#0# so they match the
// haplotype-path names inside the Minigraph-Cactus GBZ. Uses the host sed
// (no container directive).
process RENAME_GTF {
    tag "${params.ref_name}.gtf"

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path gtf

    output:
    path "gfa/${params.ref_name}.gtf", emit: renamed_gtf

    script:
    """
    mkdir -p gfa
    rename_gtf_for_vg.sh ${gtf} ${params.ref_name} gfa/${params.ref_name}.gtf
    """

    stub:
    """
    mkdir -p gfa
    cp ${gtf} gfa/${params.ref_name}.gtf
    """
}
