import groovy.json.JsonOutput

process COMBINE_BAM {
    label 'cereTARPON'
    tag "$file_type Combining BAM Files"

    input: 
        tuple val(file_type), path(input_files)

    output:
        tuple val(params.run_name), path("${file_type}.bam"), emit:combined

    script:
    """
    samtools merge -o ${file_type}.bam ${input_files} 
    """
}

process BARCODE_HAMMING_CHECK {

    label 'cereTARPON'
    tag "$params.run_name Checking Barcode Hamming Distance"

    input:
        path(sample_file)

    output:
        path("passed.txt"), optional:true

    script:
    """
    check_hamming_distance.py --sample_file ${sample_file} --barcode_errors ${params.barcode_errors}
    """
}



process CHECK_AND_CONVERT_TO_BAM {
    label 'cereTARPON'
    tag "$params.run_name Converting FASTQ to BAM"

    input:
         tuple val(sample), path(input_fh)

    output:
        tuple val(sample), path("${sample}.bam")

    script:
    if (input_fh.extension == "bam")
        """
        echo "Already BAM"
        """
    else
        """
        picard -Xmx100G FastqToSam FASTQ=${input_fh} OUTPUT=${sample}.bam SAMPLE_NAME=${sample}
        """
}


process ALIGNMENT {
    label 'cereTARPON'
    tag 'Alignment'
    cpus Math.min(params.threads as int, Runtime.runtime.availableProcessors())

    input:
        tuple val(sample_name), path(input_bam)

    output:
        tuple val(sample_name), path(input_bam), emit: input
        tuple val(sample_name), path("*.aligned.bam"), emit: alignment

    script:
    """
    samtools fastq -@ ${task.cpus} ${input_bam} > ${sample_name}.fastq
    minimap2 -ax map-ont -t ${task.cpus} ${params.cere_genome}  ${sample_name}.fastq > ${sample_name}.aligned.sam
    samtools view -b ${sample_name}.aligned.sam > ${sample_name}.aligned.bam
    """
}

process ISOLATE_TELO_READS {

    label 'cereTARPON'
    tag 'Isolating Reads'

    input:
        tuple val(sample_name), path(input_bam)
        tuple val(sample_name), path(alignment_file)
    
    output:
//        tuple val(sample_name), path("putative_read_ids.txt"), emit: putative_reads
        tuple val(sample_name), path("*putative_reads.bam"), emit: putative_reads
        path("isolation.results.txt")

    //publishDir "${params.outdir}/${sample_name}/", mode: 'copy', overwrite:true, pattern:"*.coordinates_identified.bam"
    publishDir "${params.outdir}/${sample_name}/", mode: 'copy', overwrite: true, pattern:"isolation.results.txt"

    script:
    """
    isolate_read_ids.py --alignment_file ${alignment_file} \
                                        --karyotype_file ${params.fsa_idx} \
                                        --subtelo_stretch ${params.subtelomeric_ref_stretch} \
                                        --read_ids ${sample_name}.read_ids.txt 
                                        
    samtools view -@ ${task.cpus} -N ${sample_name}.read_ids.txt ${input_bam} > ${sample_name}.filtered_reads.bam

    isolate_telomeric_reads.py --input_file ${sample_name}.filtered_reads.bam \
                                        --minimum_read_length ${params.min_read_length} \
                                        --min_repeat_ratio ${params.min_repeat_ratio} \
                                        --telomeric ${sample_name}.putative_reads.bam

    cp .command.out isolation.results.txt
    """
}

process SUBTELO_FILTER {
    label 'cereTARPON'
    tag 'Filtering Telomeres'

    input:
        tuple val(sample_name), path(input_file)

    output:
        tuple val(sample_name), path("*.subtelo_filtered.bam"), emit: filtered
        path("subtelomere_filter.results.txt")
    
    publishDir "${params.outdir}/${sample_name}/", mode: 'copy', overwrite: true, pattern:"subtelomere_filter.results.txt"

    script:
    """
    filter_by_subtelo.py --telomere_file ${input_file} \
                        --filtered_telomeres ${sample_name}.subtelo_filtered.bam \
                        --removed_telomeres ${sample_name}.removed_subtelo.bam \
                        --min_subtelo_length ${params.minimum_subtelomere_length} \
                        --min_subtelo_ratio ${params.minimum_subtelomere_ratio}

    cp .command.out subtelomere_filter.results.txt
    """
}

process IDENTIFY_TELO_COORDS {

    label 'cereTARPON'
    tag 'Identifying Telomere End and Start'

    input:
        tuple val(sample_name), path(input_file)

    output:
        tuple val(sample_name), path("*coordinates_identified.bam"), emit: telomeric
        path("telomeric_coordinates.results.txt")
    
    publishDir "${params.outdir}/${sample_name}/", mode: 'copy', overwrite: true, pattern:"telomeric_coordinates.results.txt"


    script:
    """
    identify_telomeric_coordinates.py --input_file ${input_file} \
                                        --output_file ${sample_name}.coordinates_identified.bam \
                                        --end_sequence ${params.ligation_adaptor_sequence} \
                                        --sliding_window ${params.sliding_window_size} \
                                        --interval ${params.sliding_window_interval} \
                                        --composition ${params.telomeric_composition}

    cp .command.out telomeric_coordinates.results.txt
    """
}


process PRETELO_FILTER {

    label 'cereTARPON'
    tag 'Filtering by Pre-telomeric Percentage'

    input:
        tuple val(sample_name), path(input_file)

    output:
        tuple val(sample_name), path("*telomeric.bam"), emit: telomeres
        path("pretelomeric_filtering.results.txt")

    publishDir "${params.outdir}/${sample_name}/", mode: 'copy', overwrite: true, pattern:"pretelomeric_filtering.results.txt"
    publishDir "${params.outdir}/${sample_name}/", mode: 'copy', overwrite: true, pattern:"*telomeric.bam"

    script:
    """
    filter_by_pretelomere_composition.py --telomere_file ${input_file} \
                                            --filtered_telomeres ${sample_name}.telomeric.bam \
                                            --distance ${params.pre_telomeric_distance} \
                                            --maximum_pretelo_composition ${params.maximum_pretelo_composition}

    cp .command.out pretelomeric_filtering.results.txt
    """

}

process CALCULATE_TELO_LENGTH {

    label 'cereTARPON'
    tag 'Filtering by Pre-telomeric Percentage'

    input:
        tuple val(sample_name), path(input_file)

    output:
        tuple val(sample_name), path(input_file), emit: telomeres
        tuple val(sample_name), path("*.telo_stats.txt"), emit: telo_stats
        path("*.pdf")

    publishDir "${params.outdir}/${sample_name}/", mode: 'copy', overwrite: true, pattern:"*telo_stats.txt"
    publishDir "${params.outdir}/${sample_name}/FIGURES/", mode: 'copy', overwrite: true, pattern:"*pdf"



    script:
    """
    get_stats.py --input ${input_file} --out_file ${sample_name}.telo_stats.txt
    teloPlots.R ${sample_name}.telo_stats.txt
    """
    //teloPlots.R ${sample_name}.telo_stats.txt
    

}




