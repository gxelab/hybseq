// Merge the per-sample VCFs.
// Samples are merged in samplesheet order (carried row index, sorted upstream),
// which also fixes the sample column order of the merged VCF.
process BCFTOOLS_MERGE_DNA {
    tag 'combined'

    container params.bcftools_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(vcf), path(csi), val(idx)    // collected list, sorted by idx

    output:
    path "gam/combined.vcf.gz", emit: combined_vcf

    script:
    def vcfs = vcf.join(' ')
    """
    mkdir -p gam
    bcftools merge ${vcfs} -O z -o gam/combined.vcf.gz
    """

    stub:
    """
    mkdir -p gam
    touch gam/combined.vcf.gz
    """
}
