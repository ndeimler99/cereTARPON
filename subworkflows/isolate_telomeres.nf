/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { alignment } from "../bin/process.nf"
include { isolate } from "../bin/process.nf"
include { telo_filter } from "../bin/process.nf"
include { identify_start } from "../bin/process.nf"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/ 

workflow ISOLATE_TELO_SEQS {

    take:
        input_read_fastq

    main:
        Pinguscript.ping_start(nextflow, workflow, params)

        aln_results = alignment(input_read_fastq)
        
        putative_telomeres = isolate(aln_results.alignment)

        filtered_telomeres = telo_filter(putative_telomeres.telomeric)

        telo_stats = identify_start(filtered_telomeres.telomeric)

        // group all results for html report?

    emit:
        grouped_results

}