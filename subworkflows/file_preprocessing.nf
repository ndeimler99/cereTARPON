/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CHECK_AND_CONVERT_TO_BAM } from "../bin/process.nf"
include { COMBINE_BAM } from "../bin/process.nf"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/ 
 
 workflow PREPROCESS_FILES {
    
    main:
        Pinguscript.ping_start(nextflow, workflow, params)
        print("Pre-processing Files")
 
        // if params.multiplexed
            // if files are already demultiplexed
            // if files need to be demultiplexed
        // else
        if (file(params.input).isDirectory()){
            // combine files
            // input_ch = Channel.fromPath ( "${params.input}/*" ).map{ it -> [it.baseName, it]}
            // bam_ch = CHECK_AND_CONVERT_TO_BAM(input_ch)
            // bam_ch = COMBINE_BAM(input_ch)
        }
        else {
            input_ch = Channel.fromPath (params.input).map{ it -> [it.baseName, it]}
            bam_ch = CHECK_AND_CONVERT_TO_BAM(input_ch)
        }
     
    emit:
        input = bam_ch
 }


 