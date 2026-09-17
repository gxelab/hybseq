// Prefix every GTF contig with <ref_name>#0# so the names match the haplotype
// paths inside the Minigraph-Cactus GBZ. The rule is contig-naming-scheme
// agnostic (NC_/NW_, chr*, scaffold*, ...); comment and blank lines are left
// alone. Uses the host sed (no container directive).
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
