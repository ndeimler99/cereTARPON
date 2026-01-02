#!/usr/bin/env python3 

import argparse
import pysam

def main(args):

    aln_file = pysam.AlignmentFile(args.alignment_file, "r", check_sq=False)

    count_dict = {}
    nucl_dict = {}

    for aln in aln_file:
        if aln.is_secondary or aln.is_unmapped:
            continue
        
        if aln.reference_name not in count_dict:
            count_dict[aln.reference_name] = 0
            nucl_dict[aln.reference_name] = 0

        count_dict[aln.reference_name] += 1
        nucl_dict[aln.reference_name] += len(aln.query_sequence)

    
    read_count = sum(count_dict.values())
    nucl_count = sum(nucl_dict.values())


    with open(args.out_file, "w") as out:
        out.write("Chrom\treads\tnucleotides\n")
        for chrom in count_dict:
            out.write("{}\t{}\t{}\n".format(chrom, count_dict[chrom]/read_count * 100, nucl_dict[chrom]/nucl_count * 100))

                
def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment_file", required=True)
    parser.add_argument("--out_file", required=True)

    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
