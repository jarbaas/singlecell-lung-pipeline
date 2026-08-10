// 1. Parameters Definition
params.input     = null
params.metadata  = null
params.model     = null
params.batch_key = 'Sample_Name'
params.outdir    = null

// 2. The Process Block
process SCANPY_ANALYSIS {
    
    publishDir "${params.outdir}", mode: 'copy'

    input: 
    path matrix
    path metadata
    path model
    path script_file 
    
    output:
    path '*.h5ad'
    path 'figures/*.png'
    path 'marker_genes.csv'
    
    script:
    """
    python ${script_file} \
        -i ${matrix} \
        -mod ${model} \
        -m ${metadata} \
        -b ${params.batch_key} \
        -t ${task.cpus}
    """
}

// 3. The Workflow Orchestration Block
workflow {
    
    // Initialize channels
    ch_input    = Channel.fromPath(params.input)
    ch_metadata = Channel.fromPath(params.metadata)
    ch_model    = Channel.fromPath(params.model)
    ch_script   = Channel.fromPath("${projectDir}/process_scanpy.py")

    // Execute the process
    SCANPY_ANALYSIS(ch_input, ch_metadata, ch_model, ch_script)
}