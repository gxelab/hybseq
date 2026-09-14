// Compute read support with vg pack (min MAPQ via -Q, default 5).
process VG_PACK_DNA {
    tag "$sample"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(gam), path(gam_log), val(idx)
    path gbz

    output:
    tuple val(sample), path("gam/${sample}.pack"), val(idx), emit: pack

    script:
    """
    mkdir -p gam
    vg pack -x ${gbz} -g ${gam} -o gam/${sample}.pack -t ${task.cpus} -Q ${params.min_mapq}
    """

    stub:
    """
    mkdir -p gam
    touch gam/${sample}.pack
    """
}
