process COMPUTE_COSINE {
    tag "$meta.id"
    label 'process_medium'

    // Conda environment definition.
    // Ideally, for nf-core, we should use a container that includes 'lsa' and 'pheatmap'.
    // Using base-r and installing packages or a custom container is the standard way.
    conda "r-base=4.2 r-optparse r-lsa r-pheatmap r-readr r-readxl r-dplyr"
    
    // Docker/Singularity container definition
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(expression_matrix)

    output:
    tuple val(meta), path("*_matrix.csv") , emit: matrix
    tuple val(meta), path("*_heatmap.png"), emit: heatmap
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # Nextflow automatically adds bin/ to the PATH, so we can call the script directly
    compute_cosine.R \\
        --input $expression_matrix \\
        --out_prefix ${prefix}_cosine \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | sed -n 1p | sed 's/R version //g' | sed 's/ (.*//g')
    END_VERSIONS
    """
}
