import groovy.json.JsonOutput

process ISOLATE_PUTATIVE_TELOMERES {

    label 'pombeTARPON'
    tag 'Isolating Reads'

    input:
        tuple val(sample_name), path(input_file, stageAs: "input.bam")
    
    output:
//        tuple val(sample_name), path("putative_read_ids.txt"), emit: putative_reads
        tuple val(sample_name), path("*putative_reads.bam"), emit: putative_reads
        tuple val(sample_name), path("*non_telomeric.bam"), emit: non_telomeric
        tuple val(sample_name), path("*chimeric.bam"), emit: chimeric_reads
        tuple val(sample_name), path("*reverse_complemented.bam"), emit: reverse_complemented

    script:
    """
    isolate_putative_telomeric_reads.py --input_file ${input_file} \
                                        --repeat ${params.repeat}  \
                                        --repeat_count ${params.repeat_count} \
                                        --fragment_read_length ${params.read_fragment_length} \
                                        --putative_file ${sample_name}.putative_reads.bam \
                                        --non_telomeric ${sample_name}.non_telomeric.bam \
                                        --min_repeat_threshold ${params.minimum_strand_percentage} \
                                        --chimeric_file ${sample_name}.chimeric.bam \
                                        --reverse_complement_file ${sample_name}.reverse_complemented.bam
    """
    // isolate_putative_telomeric_reads.py --input_file ${input_file} \
    //                                         --repeat ${params.repeat}  \
    //                                         --repeat_count ${params.repeat_count} \
    //                                         --fragment_read_length ${params.read_fragment_length} \
    //                                         --putative_file ${sample_name}.putative_reads.bam \
}

process BLAST_FOR_ISOLATION {
    label 'pombeTARPON'
    tag 'Blasting Against TAS and rDNA Sequences'
    stageInMode 'copy'

    input:
        tuple val(sample_name), path(input_file, stageAs: "input.bam")
        path(ref_file)
    
    output:
        tuple val(sample_name), path("*.blast_results.txt"), emit: blast_results

    script:
    """
    samtools fasta ${input_file} > ${sample_name}.fasta
    blastn -query ${sample_name}.fasta -subject ${ref_file} -outfmt "6 qseqid sseqid pident length slen mismatch gapopen qstart qend sstart send evalue bitscore qlen" > ${sample_name}.blast_results.txt
    """
}

process IDENTIFY_TELO_CORDINATES {
    
    label 'pombeTARPON'
    tag 'Isolating Reads'

    input:
        tuple val(sample_name), path(input_file, stageAs: "input.bam")
    
    output:
        tuple val(sample_name), path("*coordinates_identified.bam"), emit: coordinates_identified
        tuple val(sample_name), path("*no_coordinates.bam"), emit: no_coordinates_identified
        
    publishDir "${params.outdir}/${sample_name}/", mode: 'copy', overwrite:true, pattern:"*.coordinates_identified.bam"

    script:
    """
    identify_telo_coord.py --input_file ${input_file} \
                                        --repeat ${params.repeat}  \
                                        --telo_end_repeat_count ${params.telo_end_repeat_count} \
                                        --telo_end_repeat_errors ${params.telo_end_error_count} \
                                        --coords_file ${sample_name}.coordinates_identified.bam \
                                        --no_coords_file ${sample_name}.no_coordinates.bam \
                                        --telo_start_canonical ${params.canonical_start} \
                                        --telo_start_canonical_errors ${params.canonical_start_errors} \
                                        --telo_start_repeat_count ${params.telo_start_repeat_count}\
                                        --telo_start_repeat_errors ${params.telo_start_error_count}

    """
}

process ALIGNMENT {

    label 'pombeTARPON'
    tag "${sample_name} Aligning to Reference"

    input: 
        tuple val(sample_name), path(input_fh)
        path(POMBE_GENOME_FILE)
    output:
        tuple val(sample_name), path("*.stats.txt"), emit:alignment
        path("*.pdf")

    publishDir "${params.outdir}/${sample_name}/", mode: 'copy', overwrite:true, pattern:"*.pdf"

    script: 
    """
    samtools fastq ${input_fh} > ${sample_name}.fastq
    gzip ${sample_name}.fastq
    minimap2 -ax map-ont ${POMBE_GENOME_FILE} ${sample_name}.fastq.gz > ${sample_name}.alignment.sam 2> ${sample_name}.alignment.err
    samtools view -h -b ${sample_name}.alignment.sam > ${sample_name}.alignment.bam
    alignment_stats.py --alignment_file ${sample_name}.alignment.bam --stats_file ${sample_name}.stats.txt --sample ${sample_name}
    alignment_plots.R ${sample_name}.stats.txt
    """
}

process GET_STATS_FROM_BAM {
    label 'pombeTARPON'
    tag 'Extracting Stats from BAM File'

    input:
        tuple val(sample_name), path(input_file)

    output:
        tuple val(sample_name), path(input_file), emit: telo_bam
        tuple val(sample_name), path("*.stats.txt"), emit: telo_stats

    script:
    """
    telo_stats_from_bam.py --input_file ${input_file} --stats_file ${sample_name}.stats.txt
    """
}

process GENERATE_R_PLOTS {
    label 'pombeTARPON'
    tag 'Generating Plots'

    input:
        tuple val(sample_name), val(files)
        //blast is files[0], clusters is files[1], fasta is files[2], stast is files[3]

    output:
        path("*extended_stats.txt"), emit: stats
        path("*.pdf")

    publishDir "${params.outdir}/${sample_name}/", mode: 'copy', overwrite:true, pattern:"*.extended_stats.txt"
    publishDir "${params.outdir}/${sample_name}/FIGURES/", mode: 'copy', overwrite:true, pattern:"*.pdf"

    script:
    """
    append_to_stats.py --stats_file ${files[3]} --cluster_file ${files[1]} \
                        --fasta_file ${files[2]} \
                        --new_stats_file ${sample_name}.extended_stats.txt \
                        --minimum_cluster_size ${params.minimum_cluster_size} \
                        --repeat ${params.repeat}
    plots.R ${sample_name}.extended_stats.txt individual
    """
}

process GENERATE_HTML_REPORT {
    label 'pombeTARPON'
    tag 'Generating HTML Report'

    input:
        path(stats_files)

    output:
        path("*.html")
        path("*.pdf")

    publishDir "${params.outdir}/", mode: 'copy', overwrite:true, pattern:"*.html"
    publishDir "${params.outdir}/FIGURES/", mode: 'copy', overwrite:true, pattern:"*.pdf"

    script:
    """
    appendStats.py --stats_file ${stats_files} --outFile combined.stats.txt
    plots.R combined.stats.txt comparison
    """
}

process CONVERT_BAM_TO_FASTA_AND_TRIM {

    label 'pombeTARPON'
    tag "$file_type Converting BAM to FASTA"

    input:
        tuple val(sample_name), path(input_file)
        tuple val(sample_name), path(stats_file)

    output:
        tuple val(sample_name), path("${sample_name}.fasta"), emit: telo_fasta
        tuple val(sample_name), path("${sample_name}.trimmed.fasta"), emit: trimmed_telo_fasta

    script:
    """
    samtools fasta -@ 4 $input_file > intermediate.fasta
    trimmed_telo_reads.py --input_file intermediate.fasta --stats_file ${stats_file} --no_telo ${sample_name}.fasta \
                        --no_telo_subtelo ${sample_name}.trimmed.fasta --read_length ${params.cluster_telo_length}
    """
}

process BLAST {
    label 'pombeTARPON'
    tag 'Blasting Against TAS and rDNA Sequences'
    stageInMode 'copy'

    input:
        tuple val(sample_name), path(fasta_file)
        path(ref_file)
    
    output:
        tuple val(sample_name), path("*.blast_results.txt"), emit: blast_results

    script:
    """
    blastn -query ${fasta_file} -subject ${ref_file} -outfmt "6 qseqid sseqid pident length slen mismatch gapopen qstart qend sstart send evalue bitscore qlen" > ${sample_name}.blast_results.txt
    """
}

process MESHCLUST {
    label 'meshclust'
    tag 'Clustering using Meshclust'
    stageInMode 'copy'

    input:
        tuple val(sample_name), path(fasta_file)
    
    output:
       tuple val(sample_name), path("*.clusters.txt"), emit: clusters

    script:
    """
    meshclust -d ${fasta_file} -a y -e y -t ${params.cluster_similarity_threshold} -o ${sample_name}.clusters.txt
    """
}

process PLOT_CLUSTERS {
    label 'pycairo'
    tag 'Plotting Clusters'

    input:
        tuple val(sample_name), val(files)
        //blast is files[0], clusters is files[1], fasta is files[2], stast is files[3]

    output:
        path("*.png")

    publishDir "${params.outdir}/${sample_name}/FIGURES/", mode: 'copy', overwrite:true, pattern:"*.png"

    script:
    """
    plotClusters.py --stats_file ${files[3]} --blast_file ${files[0]} --cluster_file ${files[1]} \
        --plot_file_name ${sample_name}.clustered_plotted.png --x_axis_length ${params.cluster_plot_width} \
        --TAS_perc ${params.TAS1_TAS3_percent_identity} --TAS2_perc ${params.TAS2_percent_identity} \
        --TAS_fraction ${params.TAS_length_fraction} --image_width_px ${params.image_width_px} \
        --image_height_px ${params.image_height_px} --minimum_cluster_size ${params.minimum_cluster_size} \
        --color_dict ${params.color_dict}
    """
}

process CONSENSUS_SEQ {
    label 'pombeTARPON'
    tag 'Creating Consensus Sequence and Blasting'

    input:
        tuple val(sample_name), val(files)
        path(ref_file)
        //blast is files[0], clusters is files[1], fasta is files[2], stast is files[3]

    output:
        path("*.consensus_seqs.fa")
        path("*.pdf")

    publishDir "${params.outdir}/${sample_name}/", mode: 'copy', overwrite:true, pattern:"*.consensus_seqs.fa"
    publishDir "${params.outdir}/${sample_name}/CONSENSUS_FIGURES/", mode: 'move', overwrite:true, pattern:"*.pdf"

    script:
    """
    extract_cluster_fa.py --input_file ${files[2]} --cluster_file ${files[1]} --prefix ${sample_name} \
                            --minimum_cluster_size ${params.minimum_cluster_size}
    for cluster in \$(ls -1 *cluster*.fa); 
        do mafft \$cluster > \$cluster.aln;
    done

    create_consensus_seq.py --consensus_files *.aln --output_fasta ${sample_name}.consensus_seqs.fa \
                                --output_stats ${sample_name}.consensus_stats.txt \
                                --seq_percentage ${params.alignment_percentage}
 
    blastn -query ${sample_name}.consensus_seqs.fa -subject ${ref_file} -outfmt "6 qseqid sseqid pident length slen mismatch gapopen qstart qend sstart send evalue bitscore qlen" > ${sample_name}.consensus_blast.txt
    
    consensus_plots.R ${sample_name}.consensus_stats.txt
    

    """
}

process getParams {

    label "pombeTARPON"
    tag "Getting Parameters"

    output:
        path "params.json", emit:params

    script:
    json_str = JsonOutput.toJson(params)
    json_indented = JsonOutput.prettyPrint(json_str)

    """
    echo '${json_indented}' > "params.json"
    """
}

process getVersions {

    label "pombeTARPON"
    tag "Getting Versions"

    output:
        path "versions.txt", emit: versions

    script:
    """
    python --version | sed 's/ /,/' >> versions.txt
    python -c "import regex; print(f'regex,{regex.__version__}')" >> versions.txt
    python -c "import pandas; print(f'pandas,{pandas.__version__}')" >> versions.txt
    seqkit version | sed 's/ /,/' >> versions.txt
    """
}

process getManifest {
    
    label 'pombeTARPON'
    tag "Collecting Manifest Data"

    output:
        path "manifest.json", emit:manifest

    script:
    json_str = JsonOutput.toJson(workflow.manifest)
    json_indented = JsonOutput.prettyPrint(json_str)
    """
    echo '${json_indented}' > "manifest.json"
    """
}


process COMBINE_FASTQ {
    
    label 'pombeTARPON'
    tag "$file_type Concatenating FASTQ Files"

    input:
        tuple val(file_type), path(input_files)

    output:
        tuple val(params.run_name), path("${file_type}.fastq"), emit:combined

    script:
    """
    cat ${input_files} > ${file_type}.fastq
    """
}

process COMBINE_FASTQ_GZ {
    
    label 'pombeTARPON'
    tag "$file_type Concatenating FASTQ Files"

    input:
        tuple val(file_type), path(input_files)

    output:
        tuple val(params.run_name), path("${file_type}.fastq.gz"), emit:combined

    script:
    """
    cat ${input_files} > ${file_type}.fastq.gz
    """
}

process FASTQ_2_FASTQGZ {
    label 'pombeTARPON'
    tag "$file_type Converting FASTQ to FASTQ.GZ"
    stageInMode "copy"

    input:
        tuple val(file_type), path(input_file)
    
    output:
        tuple val(params.run_name), path("${file_type}.fastq.gz"), emit: combined

    script:
    """
    gzip $input_file 
    mv ${input_file}.gz ${file_type}.fastq.gz
    """
}



process COMBINED_FASTQ_GZ {
    
    label 'pombeTARPON'
    tag "$file_type Concatenating FASTQ Files"

    input:
        tuple val(file_type), path(input_files)

    output:
        tuple val(params.run_name), path("${file_type}.fastq.gz"), emit:combined

    script:
    """
    cat ${input_files} > ${file_type}.fastq.gz
    """
}

process COMBINE_FASTQ_AND_ZIP {
    label 'pombeTARPON'
    tag "$file_type Concatenating and Zipping FASTQ Files"

    input:
        tuple val(file_type), path(input_files)

    output:
        tuple val(params.run_name), path("${file_type}.fastq.gz"), emit:combined

    script:
    """
    cat ${input_files} > ${file_type}.fastq
    gzip ${file_type}.fastq
    """
}

process COMBINE_BAM {
    label 'pombeTARPON'
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

    label 'pombeTARPON'
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

process FASTQ_TO_BAM {
    label 'pombeTARPON'
    tag "$params.run_name Converting FASTQ to BAM"

    input:
         tuple val(run_name), path(fastq_file)

    output:
        tuple val(run_name), path("input.bam")

    script:
    """
    picard -Xmx100G FastqToSam FASTQ=${fastq_file} OUTPUT=input.bam SAMPLE_NAME=input
    """
}


process CHECK_AND_CONVERT_TO_BAM {
    label 'pombeTARPON'
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