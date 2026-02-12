// Import the module
include { COMPUTE_COSINE } from '../../modules/local/compute_cosine.nf'

workflow RUN_COSIMFLOW {
    
    take:
    ch_expression // Expected channel structure: [ [id:'sample_id'], file_path ]

    main:
    
    // Run the cosine similarity module
    COMPUTE_COSINE ( ch_expression )

    emit:
    matrix   = COMPUTE_COSINE.out.matrix
    heatmap  = COMPUTE_COSINE.out.heatmap
    versions = COMPUTE_COSINE.out.versions
}
