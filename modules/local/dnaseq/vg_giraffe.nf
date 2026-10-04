// Map reads to the pangenome graph with vg giraffe.
// Single-end (fastq_2 empty) and paired-end (two -f inputs) are both supported;
// the layout is inferred from the samplesheet. fq2 is declared `val` so empty
// strings are accepted (it is still validated via file(..., checkIfExists: true)
// when the samplesheet is parsed).
// The -d index is the rebuilt <outname>.dist that VG_INDEX_DIST2 promotes from
// <outname>.dist2; the cactus-produced index is archived as <outname>.dist.bak.
process VG_GIRAFFE {
    tag "$sample"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(fq1), val(fq2), val(idx)
    path gbz
    path min
    path zipcodes
    path dist

    output:
    tuple val(sample), path("dna/${sample}.gam"), path("dna/${sample}.gam.log"), val(idx), emit: gam

    script:
    def fq_args = fq2 ? "-f ${fq1} -f ${fq2}" : "-f ${fq1}"
    """
    mkdir -p dna
    vg giraffe \\
        -Z ${gbz} \\
        -m ${min} \\
        -z ${zipcodes} \\
        -d ${dist} \\
        ${fq_args} \\
        -t ${task.cpus} \\
        > dna/${sample}.gam 2> dna/${sample}.gam.log
    """

    stub:
    """
    mkdir -p dna
    touch dna/${sample}.gam dna/${sample}.gam.log
    """
}
