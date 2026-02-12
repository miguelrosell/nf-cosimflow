#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import the subworkflow
include { RUN_COSIMFLOW } from './subworkflows/local/run_cosimflow.nf'

workflow {

    // 1. Input validation
    if (params.input) { 
        // 2. Create the channel with the structure [ meta, file ]
        // CRITICAL: nf-core modules always expect a "meta map" as the first element of the tuple.
        ch_input = Channel.fromPath(params.input)
                          .map { file -> [ [id: file.baseName], file ] } 
    } else {
        log.error "❌ Error: Please specify an input file using --input <file.csv>"
        exit 1
    }

    // 3. Run the Subworkflow
    RUN_COSIMFLOW ( ch_input )
}
