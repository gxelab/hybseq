// Call variants per RNA-seq sample on the spliced graph.
// The -s sample name defaults to the sample id; set --rna_call_sample to override.
process VG_CALL_RNA {
    tag "$sample"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(pack), val(idx)
    path xg
    path snarls

    output:
    tuple val(sample), path("rna/${sample}.vcf.gz"), val(idx), emit: vcf

    script:
    def call_sample = params.rna_call_sample ?: sample
    """
    vg call ${xg} -r ${snarls} -k ${pack} -s ${call_sample} -z -a -t ${task.cpus} \\
        | bgzip -c > rna/${sample}.vcf.gz
    """

    stub:
    """
    mkdir -p rna
    touch rna/${sample}.vcf.gz
    """
}
