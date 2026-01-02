# cereTARPON

cereTARPON is a telomere analysis and research pipeline designed for the analysis of Saccharomyces cerevisiae telomeres that originate from whole genome sequencing datasets. This pipeline may also be suitable for the analysis of data using an enrichment methodology; however, that is untested. If this pipeline is used within your research, please cite this GitHub page. This pipeline is still under development and should be used with caution. If interested in more details,the users should contact Nathaniel Deimler at nathanieldeimler.research@gmail.com to learn more.

The NextFlow Pipeline consists of 3 subworkflows. The first validates all input parameters, the second converts files into the appropriate type, and the third analyses the sequencing data to calculate telomere length and in the future chromosome arm specific telomere statistics. This pipeline should only be used for the analysis of S. cerevisiae sequencing data. For organisms with invariant telomeric repeats please see github.com/ndeimler99/TARPON or github.com/ndeimler99/pombeTARPON in the case of S. pombe. 

Prior to the execution of the pipeline, pod5 files should be basecalled using dorado SUP basecalling with the --no-trim parameter specified. After basecalling, dorado can be used to demultiplex reads prior to the execution of cereTARPON. Currently, this pipeline needs to be run with the following parameters --input.

## Validate Parameters Subworkflow

This subworkflow first ensures that data does not exist within the specified --outdir based on the presence of an html report. The input file (--input), reference genome (--cere_genome), and reference index (-fsa_idx) must exist. As the pipeline continues to evolve and more options are added such as demultiplexing and capture probe sequences are used for enrichment, further parameters will be added here to ensure pipeline use.

## Pre-process Files Subworkflow

Currently --input must be a single file. If it is not a ubam file, it will be converted from fastq format to ubam by picard. In the future, as this pipeline expands, files will additionally be demultiplexed at this step

## Telomere Analysis Subworkflow

Sequencing data is first aligned to the cerevisiae reference genome distributed within this pipeline and publicly availabled on SGD. From this alignment, any reads aligning to the reference genome within the first or last 20000 bases (--subtelomeric_ref_stretch) of any chromosome other than the mitochondrial genome are isolated as putative telomeric sequences. If the read is less than 250bp in length (--min_read_length) it is discarded. The frequency of C strand telomeric repeats is then calculated using the first 250bp of the read and the frequency of G strand telomeric repeats is calculated in the terminal 250 bp of the telomere. If the frequency of C strand repeats / c strand + g strand repeats is less than 0.3 (--min_repeat_ratio), the read orients from the G strand. If the ratio is greater than 0.7 (1-0.3), the read is reverse complemented into the G orietnation to simplify further analysis, and if the ratio is between 0.3 and 0.7 the read is discarded as chimeric. The first 250 bp of the read (--minimum_subtelomere_length) must be composed of less than 0.3 (--minimum_subtelomere_ratio) telomeric repeats to ensure the entire read is not a telomere and the telomere start site could be appropriately identified.

The end of the putative telomere is then identified by searching for the first four nucleotides (--ligation_adaptor_sequence) of the nanopore adaptor sequence allowing for one error as this sequence differs greatly from expected telomeric repeats within the terminal 100bp of the telomeric read. The start of the telomere is then identified through a sliding window approach identifying where the frequency of 'T(G){1,3}' repeats exceeds 50% (--telomeric_composition). This window is 20 nucleotides wide (--sliding_window_size) with a 2 nucleotide jump (--sliding_window_interval). Once the telomere start window is identified, the start of the telomere is defined as the first telomeric repeat found. To ensure the start of the telomere was accurately identified the 200 bp (--pre_telomeric_distance) must be composed of less than 0.7 (--maximum_pretelo_composition). Telomere length and statistics are then returned in a .ubam and .txt file.

## Output

A variety of files is generated as output including basic plots of telomere length. All telomeric reads included in the final analys is outputted as a bam file with a corresponding .txt file that lists all descriptive statistics for thes reads and can easily be used to generate custom plots. Additionally, in the future a default HTML report will also be included in the output.




This pipeline utilizes ezcharts to aid in the generation of bokeh plots for the final HTML report. Please see the following disclaimer as well as the License provided in the LICENSE folder for use of this software. This product includes software developed by Oxford Nanopore Technologies Plc.
