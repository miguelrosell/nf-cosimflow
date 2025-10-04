nextflow.enable.dsl=2

process compute_cosine_process {
    tag { expr_file.baseName }
    publishDir "${params.outdir ?: './results'}", mode: 'copy'

    input:
    path expr_file

    output:
    path "cosine_matrix.tsv"                                         , emit: out_cosine_matrix
    path "cosine_heatmap.png"                                        , emit: out_heatmap

    // Use a conda environment or container (comment/uncomment as needed)
    //conda 'envs/r_env.yml'
    //container 'docker://rocker/tidyverse:latest'

    script:
    """
    # Convert xlsx to csv if needed and call R script
    Rscript ${projectDir}/bin/compute_cosine.R \\
        --input ${expr_file} \\
        --out_prefix cosine \\
        --method ${params.method ?: 'cosine'} \\
        --min_gene_mean ${params.min_gene_mean ?: 0.0} \\
        --sample_cols "${params.sample_cols ?: 'auto'}"
    """
}

workflow.onComplete { println "compute_cosine module finished" }
