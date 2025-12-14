/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ALIGNMENT } from "../bin/process.nf"
include { ISOLATE_TELO_READS } from "../bin/process.nf"
// include { FILTER_TELOMERES } from "../bin/process.nf"
// include { IDENTIFY_TELO_COORDS } from "../bin/process.nf"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/ 

workflow ISOLATE_TELO_SEQS {

    take:
        input_bam

    main:
        Pinguscript.ping_start(nextflow, workflow, params)

        aln_results = ALIGNMENT(input_bam)
        
        putative_telomeres = ISOLATE_TELO_READS(aln_results.alignment)

        // filtered_telomeres = FILTER_TELOMERES(putative_telomeres.telomeric)

        // telo_stats = IDENTIFY_TELO_COORDS(filtered_telomeres.telomeric)

        // group all results for html report?

    // emit:
    //     grouped_results

}