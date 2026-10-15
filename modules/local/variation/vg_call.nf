// Call variants per sample with vg call, compressed with bgzip.
// Shared by both assays; the two graph-dependent differences are explicit inputs:
//   graph_is_gbz      adds `-z`, which is only valid for a GBZ input (DNA-seq);
//                     the spliced xg of the RNA-seq workflow must not get it.
//   call_sample_name  `vg call -s` ('' -> the sample id); the RNA-seq workflow
//                     passes params.rna_call_sample.
process VG_CALL {
    tag "$sample"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(pack), val(idx)
    path graph             // <outname>.gbz (DNA-seq) or <outname>_spliced.xg (RNA-seq)
    path snarls
    val assay_dir          // 'dna' | 'rna'
    val graph_is_gbz       // true -> the graph is a GBZ (DNA-seq): pass -z
    val call_sample_name   // '' -> use the sample id

    output:
    tuple val(sample), path("${assay_dir}/${sample}.vcf.gz"), val(idx), emit: vcf

    script:
    def call_sample = call_sample_name ?: sample
    def call_args   = graph_is_gbz ? '-z -a' : '-a'
    """
    mkdir -p ${assay_dir}
    vg call ${graph} -r ${snarls} -k ${pack} -s ${call_sample} ${call_args} -t ${task.cpus} \\
        | bgzip -c > ${assay_dir}/${sample}.vcf.gz
    """

    stub:
    """
    mkdir -p ${assay_dir}
    touch ${assay_dir}/${sample}.vcf.gz
    """
}
