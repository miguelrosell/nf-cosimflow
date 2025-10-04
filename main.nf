#workflow entrypoint, DSL2, This pipeline receives the expression file, routes it to the processing module, and retrieves the computed results.
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {
    params.expr = params.expr ?: ''
    if (!params.expr) {
        log.error "Please provide --expr <expression_file.csv|xlsx>"
        System.exit(1)
    }

    Channel
        .fromPath(params.expr)
        .set { expr_file_ch }

    // call module
    compute_cosine = modules.compute_cosine.call(expr_file_ch)

    // collect outputs (for demo)
    compute_cosine.out_cosine_matrix.subscribe { file -> println "Cosine matrix: ${file}" }
    compute_cosine.out_heatmap.subscribe { file -> println "Heatmap: ${file}" }
}
