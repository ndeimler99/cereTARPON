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
            if (params.demultiplexed){
                input_ch = Channel.fromPath ( "${params.input}/*" ).map{ it -> [it.baseName, it]}
                bam_ch = CHECK_AND_CONVERT_TO_BAM(input_ch)
            }
            else {
                // combine_files
                input_ch = Channel.fromPath ( "${params.input}/*" ).map{ it -> [it.baseName, it]}
                bam_ch = CHECK_AND_CONVERT_TO_BAM(input_ch)
                bam_ch = bam_ch.map { key,value -> value}
                bam_ch = COMBINE_BAM(bam_ch.collect().map { it -> ["${params.sample_name}", it]})
            }
         
        }
        else {
            input_ch = Channel.fromPath (params.input).map{ it -> [it.baseName, it]}
            bam_ch = CHECK_AND_CONVERT_TO_BAM(input_ch)
        }

    emit:
        input = bam_ch
 }


 