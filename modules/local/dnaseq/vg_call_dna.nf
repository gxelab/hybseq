// Call variants per sample with vg call (-z -a), compressed with bgzip.
process VG_CALL_DNA {
    tag "$sample"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(pack), val(idx)
    path gbz
    path snarls

    output:
    tuple val(sample), path("gam/${sample}.vcf.gz"), val(idx), emit: vcf

    script:
    """
    mkdir -p gam
    vg call ${gbz} -r ${snarls} -k ${pack} -s ${sample} -z -a -t ${task.cpus} \\
        | bgzip -c > gam/${sample}.vcf.gz
    """

    stub:
    """
    mkdir -p gam
    touch gam/${sample}.vcf.gz
    """
}
