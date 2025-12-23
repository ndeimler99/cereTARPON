#!/usr/bin/env python3 

import argparse
import pysam
import regex

def main(args):
    
    args.distance = int(args.distance)
    args.maximum_pretelo_composition = float(args.maximum_pretelo_composition)

    aln_file = pysam.AlignmentFile(args.telomere_file, "r", check_sq=False)
    out = pysam.AlignmentFile(args.filtered_telomeres, "w", template=aln_file)
    
    filtered_count = 0
    total_count = 0
    
    for aln in aln_file.fetch(until_eof=True):
        total_count += 1
        telo_start = aln.get_tag("TS")
        seq = aln.query_sequence[(telo_start-args.distance):telo_start]
        matches = list(regex.finditer(r'T(G){1,3}',seq))
        span = 0
        for match in matches:
            span += (match.span()[1] - match.span()[0])
        comp = span/args.distance

        if comp <= args.maximum_pretelo_composition:
            out.write(aln)
        else:
            filtered_count += 1
    
    print("Number of Telomeric Sequences: {}".format(total_count))
    print("Number of Filtered Telomeres: {}".format(filtered_count))


def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--telomere_file", required=True)
    parser.add_argument("--filtered_telomeres", required=True)
    parser.add_argument("--distance", required=True)
    parser.add_argument("--maximum_pretelo_composition", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
