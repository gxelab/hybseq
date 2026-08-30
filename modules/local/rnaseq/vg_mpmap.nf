// Map reads to the spliced graph with vg mpmap (-n RNA -l short).
// Single-end (fastq_2 empty) and paired-end (two -f inputs, paired by read name)
// are both supported; the layout is inferred from the samplesheet. fq2 is
// declared `val` so empty strings are accepted (it is still validated via
// file(..., checkIfExists: true) when the samplesheet is parsed).
process VG_MPMAP {
    tag "$sample"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(fq1), val(fq2), val(idx)
    path xg
    path gcsa
    path gcsa_lcp
    path dist

    output:
    tuple val(sample), path("rna/${sample}.gam"), val(idx), emit: gam

    script:
    def fq_args = fq2 ? "-f ${fq1} -f ${fq2}" : "-f ${fq1}"
    """
    mkdir -p rna
    vg mpmap -x ${xg} -g ${gcsa} -d ${dist} -n RNA -l short -F GAM -t ${task.cpus} ${fq_args} \\
        > rna/${sample}.gam
    """

    stub:
    """
    mkdir -p rna
    touch rna/${sample}.gam
    """
}
