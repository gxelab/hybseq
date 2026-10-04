// Tabulate the merged VCF with bcftools query (one row per record,
// CHROM POS ID REF ALT QUAL FILTER + per-sample GT DP AD{0} AD{1} GQ).
process BCFTOOLS_QUERY_DNA {
    tag 'combined'

    container params.bcftools_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path combined_vcf

    output:
    path "dna/combined.vcf.tsv.gz", emit: tsv

    script:
    """
    mkdir -p dna
    bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER[\\t%GT\\t%DP\\t%AD{0}\\t%AD{1}\\t%GQ]\\n' \\
        ${combined_vcf} | gzip -c > dna/combined.vcf.tsv.gz
    """

    stub:
    """
    mkdir -p dna
    touch dna/combined.vcf.tsv.gz
    """
}
