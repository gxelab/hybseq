// Compute read support on the spliced xg graph with vg pack (min MAPQ via -Q, default 5).
process VG_PACK_RNA {
    tag "$sample"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(gam), val(idx)
    path xg

    output:
    tuple val(sample), path("rna/${sample}.pack"), val(idx), emit: pack

    script:
    """
    mkdir -p rna
    vg pack -x ${xg} -g ${gam} -o rna/${sample}.pack -t ${task.cpus} -Q ${params.min_mapq}
    """

    stub:
    """
    mkdir -p rna
    touch rna/${sample}.pack
    """
}
