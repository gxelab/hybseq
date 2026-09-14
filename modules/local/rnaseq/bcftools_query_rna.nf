// Tabulate each RNA-seq VCF with bcftools query (one row per record,
// CHROM POS ID REF ALT QUAL FILTER + per-sample GT DP AD{0} AD{1} GQ).
// No merge is performed for RNA-seq; per-sample TSVs are kept as such.
process BCFTOOLS_QUERY_RNA {
    tag "$sample"

    container params.bcftools_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(vcf), val(idx)

    output:
    tuple val(sample), path("rna/${sample}.vcf.tsv.gz"), val(idx), emit: tsv

    script:
    """
    mkdir -p rna
    bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER[\\t%GT\\t%DP\\t%AD{0}\\t%AD{1}\\t%GQ]\\n' \\
        ${vcf} | gzip -c > rna/${sample}.vcf.tsv.gz
    """

    stub:
    """
    mkdir -p rna
    touch rna/${sample}.vcf.tsv.gz
    """
}
