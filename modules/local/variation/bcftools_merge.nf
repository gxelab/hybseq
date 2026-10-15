// Merge the per-sample VCFs of one assay.
// Samples are merged in samplesheet order (the carried row index, sorted upstream),
// which also fixes the sample column order of the merged VCF; the workflow reshapes
// the sorted (sample, vcf, csi, idx) tuples into the tuple of lists declared below.
process BCFTOOLS_MERGE {
    tag "${assay_dir}/combined"

    container params.bcftools_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(vcf), path(csi), val(idx)    // collected lists, sorted by idx
    val assay_dir                                        // 'dna' | 'rna'

    output:
    path "${assay_dir}/combined.vcf.gz", emit: combined_vcf

    script:
    def vcfs = vcf.join(' ')
    """
    mkdir -p ${assay_dir}
    bcftools merge ${vcfs} -O z -o ${assay_dir}/combined.vcf.gz
    """

    stub:
    """
    mkdir -p ${assay_dir}
    touch ${assay_dir}/combined.vcf.gz
    """
}
