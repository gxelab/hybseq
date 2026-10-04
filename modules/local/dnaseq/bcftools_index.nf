// Index each sample VCF before merging (`for i in *.vcf.gz; do bcftools index -f $i; done`).
//
// The input vcf is re-emitted only to carry it into the merge; the publish pattern
// `*.csi` means only the index is published here (the VCF itself is already
// published by VG_CALL_DNA under dna/).
process BCFTOOLS_INDEX {
    tag "$sample"

    container params.bcftools_container

    publishDir [path: "${params.outdir}/dna", pattern: '*.csi', mode: params.publish_dir_mode]

    input:
    tuple val(sample), path(vcf), val(idx)

    output:
    tuple val(sample), path(vcf), path("${vcf}.csi"), val(idx), emit: vcf_csi

    script:
    """
    bcftools index -f ${vcf}
    """

    stub:
    """
    touch ${vcf}.csi
    """
}
