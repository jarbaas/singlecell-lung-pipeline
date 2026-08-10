params.input           = null
params.metadata        = null
params.model           = null
params.batch_key       = 'Sample_Name'
params.outdir          = null
params.max_mt          = 15.0
params.container_image = null

// Catch missing containers before launching jobs


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
        --max_mt ${params.max_mt} \
        -t ${task.cpus}
    """
}


workflow {
    
    if (!params.container_image) {
        error "--container_image is not defined. You must provide a valid Docker/ECR URI."
    }

    ch_input    = Channel.fromPath(params.input)
    ch_metadata = Channel.fromPath(params.metadata)
    ch_model    = Channel.fromPath(params.model)
    ch_script   = Channel.fromPath("${projectDir}/process_scanpy.py")

    
    SCANPY_ANALYSIS(ch_input, ch_metadata, ch_model, ch_script)
}