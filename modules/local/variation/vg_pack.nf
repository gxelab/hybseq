// Compute read support with vg pack (min MAPQ via -Q, default 5).
// Shared by both assays: DNA-seq packs against the GBZ, RNA-seq against the
// spliced xg; assay_dir ('dna'|'rna') is the published output subdirectory.
process VG_PACK {
    tag "$sample"

    container params.vg_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(gam), val(idx)
    path graph          // <outname>.gbz (DNA-seq) or <outname>_spliced.xg (RNA-seq)
    val assay_dir       // 'dna' | 'rna'

    output:
    tuple val(sample), path("${assay_dir}/${sample}.pack"), val(idx), emit: pack

    script:
    """
    mkdir -p ${assay_dir}
    vg pack -x ${graph} -g ${gam} -o ${assay_dir}/${sample}.pack -t ${task.cpus} -Q ${params.min_mapq}
    """

    stub:
    """
    mkdir -p ${assay_dir}
    touch ${assay_dir}/${sample}.pack
    """
}
