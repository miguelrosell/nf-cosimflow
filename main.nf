#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================
// nf-cosimflow - pipeline principal (Nextflow DSL2)
// ============================================

include { compute_cosine } from './modules/compute_cosine/main.nf'

workflow {

    // Validate input param (params.expr)
    params.expr = params.expr ?: ''
    if (!params.expr) {
        log.error "❌ Please provide the expression file with --expr <expression_file.csv|xlsx>"
        System.exit(1)
    }

    // Create channel from input file
    expr_file_ch = Channel.fromPath(params.expr)

    // Call the module and capture returned object
    def result = compute_cosine(expr_file_ch)

    // Subscribe to the emitted output channels
    result.out_cosine_matrix.subscribe { file -> println "✅ Cosine matrix generated: ${file}" }
    result.out_heatmap.subscribe { file -> println "✅ Heatmap generated: ${file}" }
}
