// 1. Parameters & Channels Definition
params.input     = null
params.metadata  = null
params.model     = null
params.batch_key = 'Sample_Name'
params.outdir    = null

// 2. The Process Block
process SCANPY_ANALYSIS {
    
    // Route outputs to S3
    publishDir "${params.outdir}", mode: 'copy'

    input: 
    path matrix
    path metadata
    path model
    
    output:
    path '*.h5ad'
    path '*.png'
    path 'marker_genes.csv'
    
    script:
    """
    python ${projectDir}/process_scanpy.py \
        -i ${matrix} \
        -mod ${model} \
        -m ${metadata} \
        -b ${params.batch_key} \
        -t ${task.cpus}
    """
}

// 3. The Workflow Orchestration Block
workflow {
    
    // Initialize channels from S3 paths
    ch_input    = Channel.fromPath(params.input)
    ch_metadata = Channel.fromPath(params.metadata)
    ch_model    = Channel.fromPath(params.model)

    // Execute the process with the staged channels
    SCANPY_ANALYSIS(ch_input, ch_metadata, ch_model)
}