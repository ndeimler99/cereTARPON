#!/usr/bin/env python3 

import argparse
import pysam
import numpy as np

def main(args):

    aln_file = pysam.AlignmentFile(args.input_file, "r", check_sq=False)


    with open(args.out_file, "w") as out:
        out.write("read_id\tstrand\ttelo_start\ttelo_end\ttelo_length\tread_quality\ttelo_quality\tread_length\n")
        for aln in aln_file.fetch(until_eof=True):
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(aln.query_name, aln.get_tag("XS"), aln.get_tag("TS"), aln.get_tag("TE"),
                                                    aln.get_tag("TE")-aln.get_tag("TS"), np.mean(aln.query_qualities),
                                                    np.mean(aln.query_qualities[aln.get_tag("TS"):aln.get_tag("TE")]), 
                                                    len(aln.query_sequence)))

                
def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--out_file", required=True)

    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
