/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BARCODE_HAMMING_CHECK } from "../bin/process.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow validate_parameters {
    
    main:
        Pinguscript.ping_start(nextflow, workflow, params)
        parameters_passed = true
        print("Checking parameters")

        // Checking to see if report.html already exists in outdir to prevent accidental overwrite of data
        try {
            file("${params.outdir}/report.html", checkIfExists:true)
            if (!params.overwrite_outdir) {
                parameters_passed = false
                println("Out Directory Already Exists, Please Provide New Out Directory Name or Allow Overwriting of Pre-existing directory")
            }
        }
        catch (Exception e) {
            
        }

        // checking to see if input files/directory exist
        try {
            file(params.input, checkIfExists:true)
        }
        catch (Exception e) {
            parameters_passed = false
            println("Error - Input File or Directory Does not Exist")
        }

        try {
            file(params.cere_genome, checkIfExists:true)
            file(params.fsa_idx, checkIfExists: true)
        }
        catch(Exception e){
            parameters_passed = false
            println("Error - Cerevisiae Genome Fasta or Index File Does not Exist")
        }

        // if demux is specified sample file should not be provided and nanopore barcodes should be false
        // if (params.multiplexed && params.sample_file == ''){
        //     parameters_passed = false
        //     println ("Sample File Must Be Specified")
        // }

        // If a sample file is specified, it must be a valid file
        // try {
        //     if (params.sample_file != ""){
        //         file(params.sample_file, checkIfExists:true)
        //     }
        // }
        // catch (Exception e) {
        //     parameters_passed = false
        //     println("Error - Sample File not Found")
        // }

        // //check to ensure barcodes hamming distance is greater than the number of allowable errors in the barcode
        // if (params.sample_file != ''){
        //     try {
        //         barcode_check = BARCODE_HAMMING_CHECK(file(params.sample_file))
        //     }
        //     catch (Exception e){
        //         parameters_passed = false
        //         println "Supplied Barcode Sequences are too Similar for Demultiplexing with ${params.barcode_errors} Errors Allowed. Please reduce error amount."
        //     }
        // }

    emit:
        passed = parameters_passed
}