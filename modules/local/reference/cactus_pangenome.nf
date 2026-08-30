// Build the pangenome graph with Minigraph-Cactus (cactus 3.1.4).
// Runs as a managed foreground task; --logFile is written inside --outDir so the
// published layout under results/gfa/ is uniform.
process CACTUS_PANGENOME {
    tag "$params.outname"

    container params.cactus_container

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path assemblies

    output:
    path "gfa/${params.outname}.gbz",                   emit: gbz
    path "gfa/${params.outname}.dist",                  emit: dist
    path "gfa/${params.outname}.shortread.withzip.min", emit: min
    path "gfa/${params.outname}.shortread.zipcodes",    emit: zipcodes
    path "gfa/${params.outname}.snarls",                emit: snarls
    path "gfa/${params.outname}.log",                   emit: log
    path "gfa/${params.outname}.gfa", emit: gfa, optional: true    // --gfa clip filter full
    path "gfa/${params.outname}.vcf", emit: vcf, optional: true    // --vcf

    script:
    """
    mkdir -p gfa
    cactus-pangenome js ${assemblies} \\
        --outDir gfa \\
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
        --logFile gfa/${params.outname}.log
    """

    stub:
    """
    mkdir -p gfa
    touch gfa/${params.outname}.gbz \\
          gfa/${params.outname}.dist \\
          gfa/${params.outname}.shortread.withzip.min \\
          gfa/${params.outname}.shortread.zipcodes \\
          gfa/${params.outname}.snarls \\
          gfa/${params.outname}.log \\
          gfa/${params.outname}.gfa \\
          gfa/${params.outname}.vcf
    """
}
