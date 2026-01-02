/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ALIGNMENT } from "../bin/process.nf"
include { ISOLATE_TELO_READS } from "../bin/process.nf"
include { SUBTELO_FILTER } from "../bin/process.nf"
include { IDENTIFY_TELO_COORDS } from "../bin/process.nf"
include { PRETELO_FILTER } from "../bin/process.nf"
include { CALCULATE_TELO_LENGTH } from "../bin/process.nf"
include { GET_ALIGNMENT_STATS } from "../bin/process.nf"
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

        aln_stats = GET_ALIGNMENT_STATS(aln_results.alignment)
        
        putative_telomeres = ISOLATE_TELO_READS(aln_results.input, aln_results.alignment)

        filtered_telomeres = SUBTELO_FILTER(putative_telomeres.putative_reads)

        identified_telomeres = IDENTIFY_TELO_COORDS(filtered_telomeres.filtered)

        final_telomeres = PRETELO_FILTER(identified_telomeres.telomeric)

        telomeric_stats = CALCULATE_TELO_LENGTH(final_telomeres.telomeres)
        // group all results for html report?

    emit:
        aln_stats = aln_stats.alignment_stats
        telo_stats = telomeric_stats.telo_stats

}