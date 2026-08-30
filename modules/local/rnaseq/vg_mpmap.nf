// Map single-end RNA-seq reads to the spliced graph with vg mpmap (-n RNA -l short).
process VG_MPMAP {
    tag "$sample"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(fq), val(idx)
    path xg
    path gcsa
    path gcsa_lcp
    path dist

    output:
    tuple val(sample), path("rna/${sample}.gam"), val(idx), emit: gam

    script:
    """
    mkdir -p rna
    vg mpmap -x ${xg} -g ${gcsa} -d ${dist} -n RNA -l short -F GAM -t ${task.cpus} -f ${fq} \\
        > rna/${sample}.gam
    """

    stub:
    """
    mkdir -p rna
    touch rna/${sample}.gam
    """
}
