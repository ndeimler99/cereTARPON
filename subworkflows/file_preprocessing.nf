/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CHECK_AND_CONVERT_TO_FASTQ } from "../bin/process.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/ 
 
 workflow PREPROCESS_FILES {
    
    main:
        Pinguscript.ping_start(nextflow, workflow, params)
        print("Pre-processing Files")
 
        
        // if input is directory and --demux is true (nanopore demux has already been run)
        if (file(params.input).isDirectory() && params.demuxed){
            // do not combine files
            // map all files to their file_prefix
            input_ch = Channel.fromPath ( "${params.input}/*" ).map{ it -> [it.baseName, it]}
            
            // if file suffix is .bam convert to fastq fil
            fastq_ch = CHECK_AND_CONVERT_TO_FASTQ(input_ch)
        }
        // if --nanopore_barcodes is true
        // else if (params.nanopore_barcodes){
        //     if (file(params.input).isDirectory()){
        //         //combine all files // why am I combining all files can dorado emux operate on directory of files?
        //     }
        //     // demux using nanopore dorado demux
        //         // output as bam file
        //     // map all files to their file_prefix
        //     // isolate telomere reads and return pre-isolated channel, isolate_channel, removed channel
        // }
        // if input is single bam file and barcodes file is provided and capture_probe_sequence is none
        // else if (!params.capture_probe_sequence && params.sample_file == ""){
        //     if (file(params.input).isDirectory()){
        //         //merge files and convert if necesary
        //     }
        //     else {
        //         // convert if necesary
        //     }
        //     // map file to sample_name
        //     // isolate telomere reads and return pre-isolated channel, isolate_channel, removed channel
        // }
        // if input is single bam file and barcodes file is provided and capture_probe_sequence is provided
        // else if (params.capture_probe_sequence != false && params.sample_file != ""){    
        //     if (file(params.input).isDirectory()){
        //         //merge files and convert if necesary
        //     }
        //     else {
        //         // convert if necesary
        //     }
        //     // isolate telomere reads
        //     // demultiplex by adaptor sequence and barcode file
        //     // return pre-isolated channel, isolate_channel, removed channel
        // }
        // else if (!params.capture_probe_sequence && params.sample_file != ""){
        //     if (file(params.input).isDirectory()){
        //         //merge files and convert if necesary
        //     }
        //     else {
        //         //convert if necesary
        //     }
        //     // isolate_telomere_reads
        //     // demultiplex by sample_file barcodes
        //     // end of telomere is now end of read due to demultiplexing logic
        // }
        // else if (params.capture_probe_sequence != false && params.sample_file == ""){
        //     if (file(params.input).isDirectory()){
        //         //merge files and convert if necesaryy
        //     }
        //     else {
        //         //convert if necesaryy
        //     }
        //     // one sample in sequencng run
        //     // end of telomere is adaptor sequence
        // }

    emit:
        input = fastq_ch
 }


 