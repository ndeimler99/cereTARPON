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
        tuple val(sample_name), path("*.aligned.bam"), emit: alignment

    script:
    """
    samtools fastq ${input_bam} > ${sample_name}.fastq
    minimap2 -ax map-ont -t ${task.cpus} ${sample_name}.fastq ${params.cere_genome} | samtools sort -@ ${task.cpus} | samtools view -b > ${sample_name}.aligned.bam
    """
}

process ISOLATE_TELO_READS {

    label 'cereTARPON'
    tag 'Isolating Reads'

    input:
        tuple val(sample_name), path(input_file, stageAs: "input.bam")
    
    output:
//        tuple val(sample_name), path("putative_read_ids.txt"), emit: putative_reads
        tuple val(sample_name), path("*putative_reads.bam"), emit: putative_reads
        tuple val(sample_name), path("*non_telomeric.bam"), emit: non_telomeric
        tuple val(sample_name), path("*chimeric.bam"), emit: chimeric_reads
        tuple val(sample_name), path("*reverse_complemented.bam"), emit: reverse_complemented
    //publishDir "${params.outdir}/${sample_name}/", mode: 'copy', overwrite:true, pattern:"*.coordinates_identified.bam"

    script:
    """
    isolate_putative_telomeric_reads.py --alignment_file ${input_file} \
                                        --karyotype_file ${params.fsa_idx} \
                                        --minimum_read_length ${params.min_read_length} \
                                        --min_repeat_ratio ${params.min_repeat_ratio} \
                                        --telomeric ${sample_name}.putative_reads.bam \
                                        --non_telomeric ${sample_name}.non_telomeric.bam \
                                        --chimeric_file ${sample_name}.chimeric.bam
    """
}




