#!/usr/bin/env python3 

import argparse
import pysam
import regex

def main(args):

    args.min_subtelo_length = int(args.min_subtelo_length)
    args.min_subtelo_ratio = float(args.min_subtelo_ratio)

    aln_file = pysam.AlignmentFile(args.telomere_file, "rb", check_sq=False)
    filtered_file = pysam.AlignmentFile(args.filtered_telomeres, "w", template=aln_file)
    removed_file = pysam.AlignmentFile(args.removed_telomeres, "w", template=aln_file)

    reads = 0
    filtered_reads = 0
    for aln in aln_file.fetch(until_eof=True):
        reads += 1
        matches = list(regex.finditer(r'T(G){1,3}T(G){1,3}T(G){1,3}T(G){1,3}T(G){1,3}', aln.query_sequence[0:args.min_subtelo_length], overlapped=True))
        intervals = [(m.start(), m.end()) for m in matches]

        # Merge overlapping intervals
        intervals.sort()
        merged = []
        for start, end in intervals:
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)

        # Compute total covered length
        total_length = sum(end - start for start, end in merged)
        if total_length / args.min_subtelo_length < args.min_subtelo_ratio:
            filtered_file.write(aln)
            filtered_reads += 1
        else:
            removed_file.write(aln)
    
    print("Total Number of Reads: {}".format(reads))
    print("Number of Retained Telomeric Sequences: {}, {}".format(filtered_reads, filtered_reads/reads * 100))
   
def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--telomere_file", required=True)
    parser.add_argument("--filtered_telomeres", required=True)
    parser.add_argument("--removed_telomeres", required = True)
    parser.add_argument("--min_subtelo_length", required=True)
    parser.add_argument("--min_subtelo_ratio", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
