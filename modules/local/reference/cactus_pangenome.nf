// Build the pangenome graph with Minigraph-Cactus (cactus 3.1.4).
// Runs as a managed foreground task; --logFile is written inside --outDir so the
// published layout under results/ref/ is uniform.
// The cactus-produced distance index is archived as <outname>.dist.bak; the
// primary <outname>.dist is replaced by the index VG_INDEX_DIST2 rebuilds.
process CACTUS_PANGENOME {
    tag "$params.outname"

    container params.cactus_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path assemblies

    output:
    path "ref/${params.outname}.gbz",                   emit: gbz
    path "ref/${params.outname}.dist.bak",              emit: dist_bak
    path "ref/${params.outname}.shortread.withzip.min", emit: min
    path "ref/${params.outname}.shortread.zipcodes",    emit: zipcodes
    path "ref/${params.outname}.snarls",                emit: snarls
    path "ref/${params.outname}.log",                   emit: log
    path "ref/${params.outname}.gfa", emit: gfa, optional: true    // --gfa clip filter full
    path "ref/${params.outname}.vcf", emit: vcf, optional: true    // --vcf

    script:
    def toilWorkDir = "${params.outdir}/toil_work"
    """
    mkdir -p ${toilWorkDir}
    mkdir -p ref
    cactus-pangenome js ${assemblies} \\
        --outDir ref \\
        --outName ${params.outname} \\
        --reference ${params.ref_name} \\
        --giraffe clip filter \\
        --gbz clip filter full \\
        --gfa clip filter full \\
        --vcf \\
        --permissiveContigFilter \\
        --haplo \\
        --chrom-vg clip filter \\
        --chrom-og full \\
        --viz \\
        --consCores ${params.cactus_cons_cores} \\
        --indexCores ${params.cactus_index_cores} \\
        --mgCores ${params.cactus_mg_cores} \\
        --mapCores ${params.cactus_map_cores} \\
        --batchSystem single_machine \\
        --logFile ref/${params.outname}.log

    gunzip -c ref/${params.outname}.gfa.gz > ref/${params.outname}.gfa
    gunzip -c ref/${params.outname}.vcf.gz > ref/${params.outname}.vcf
    mv ref/${params.outname}.dist ref/${params.outname}.dist.bak

    rm -rf ${toilWorkDir}
    """

    stub:
    """
    mkdir -p ref
    touch ref/${params.outname}.gbz \\
          ref/${params.outname}.dist.bak \\
          ref/${params.outname}.shortread.withzip.min \\
          ref/${params.outname}.shortread.zipcodes \\
          ref/${params.outname}.snarls \\
          ref/${params.outname}.log \\
          ref/${params.outname}.gfa \\
          ref/${params.outname}.vcf
    """
}
