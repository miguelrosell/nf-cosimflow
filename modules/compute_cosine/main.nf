nextflow.enable.dsl=2

process compute_cosine_process {
    tag { expr_file.baseName }

    // Save final results automatically in /results/
    publishDir "${params.outdir ?: './results'}", mode: 'copy'

    // Use conda environment for reproducibility
    conda "${projectDir}/envs/r_env.yml"

    input:
    path expr_file

    output:
    path "cosine_matrix.csv", emit: out_cosine_matrix
    path "cosine_heatmap.png", emit: out_heatmap

    script:
    """
    Rscript ${projectDir}/bin/compute_cosine.R \\
        --input ${expr_file} \\
        --out_prefix cosine \\
        --method ${params.method ?: 'cosine'} \\
        --min_gene_mean ${params.min_gene_mean ?: 0.0} \\
        --sample_cols "${params.sample_cols ?: 'auto'}"
    """
}

workflow compute_cosine {
    take:
    expr_file_ch

    main:
    compute_cosine_process(expr_file_ch)

    emit:
    out_cosine_matrix = compute_cosine_process.out.out_cosine_matrix
    out_heatmap       = compute_cosine_process.out.out_heatmap
}
