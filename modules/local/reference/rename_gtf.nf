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
    path "ref/${params.ref_name}.gtf", emit: renamed_gtf

    script:
    """
    mkdir -p ref
    rename_gtf_for_vg.sh ${gtf} ${params.ref_name} ref/${params.ref_name}.gtf
    """

    stub:
    """
    mkdir -p ref
    cp ${gtf} ref/${params.ref_name}.gtf
    """
}
