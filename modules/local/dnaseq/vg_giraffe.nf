// Map reads to the pangenome graph with vg giraffe.
// Single-end (fastq_2 empty) and paired-end (two -f inputs) are both supported;
// the layout is inferred from the samplesheet. fq2 is declared `val` so empty
// strings are accepted (it is still validated via file(..., checkIfExists: true)
// when the samplesheet is parsed).
// The -d index is the cactus-produced <outname>.dist (not <outname>.dist2, see
// VG_INDEX_DIST2).
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
    tuple val(sample), path("gam/${sample}.gam"), path("gam/${sample}.gam.log"), val(idx), emit: gam

    script:
    def fq_args = fq2 ? "-f ${fq1} -f ${fq2}" : "-f ${fq1}"
    """
    mkdir -p gam
    vg giraffe \\
        -Z ${gbz} \\
        -m ${min} \\
        -z ${zipcodes} \\
        -d ${dist} \\
        ${fq_args} \\
        -t ${task.cpus} \\
        > gam/${sample}.gam 2> gam/${sample}.gam.log
    """

    stub:
    """
    mkdir -p gam
    touch gam/${sample}.gam gam/${sample}.gam.log
    """
}
