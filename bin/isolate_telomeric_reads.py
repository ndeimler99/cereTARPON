#!/usr/bin/env python3 

import argparse
import pysam
import regex

def rev_comp(seq):
    rev_dict = {"G":"C", "C":"G", "T":"A", "A":"T"}
    return ''.join([rev_dict[nucl] for nucl in seq[::-1]])

def calculate_N50(lengths):
    lengths.sort()
    lsum = sum(lengths)
    rsum = 0
    N50=None
    for l in lengths:
        rsum += l
        if rsum > lsum/2:
            N50 = l
            return N50

def main(args):

    args.minimum_read_length = int(args.minimum_read_length)
    args.min_repeat_ratio = float(args.min_repeat_ratio)

    aln_file = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)
    telomeric = pysam.AlignmentFile(args.telomeric, "w", template=aln_file)

    telo_lengths = []
    too_short = 0
    g_strand = 0
    c_strand = 0
    chimeras = 0

    for aln in aln_file.fetch(until_eof=True):
        if len(aln.query_sequence) < args.minimum_read_length:
            too_short += 1        
            continue
        
        c_matches = len(list(regex.finditer(r'(C){1,3}A(C){1,3}A(C){1,3}A(C){1,3}A(C){1,3}A', aln.query_sequence[0:250], overlapped=True)))
        g_matches = len(list(regex.finditer(r'T(G){1,3}T(G){1,3}T(G){1,3}T(G){1,3}T(G){1,3}', aln.query_sequence[-250:], overlapped=True)))
        if c_matches + g_matches != 0:
            if c_matches / (c_matches + g_matches) > args.min_repeat_ratio and c_matches / (c_matches + g_matches) < 1-args.min_repeat_ratio:
                chimeras += 1
            elif c_matches / (c_matches + g_matches) >= 1-args.min_repeat_ratio:
                c_strand += 1
                q = aln.query_qualities
                aln.query_sequence = rev_comp(aln.query_sequence)
                aln.query_qualities = q[::-1]
                aln.set_tag("XS", "C")
                telomeric.write(aln)
                telo_lengths.append(len(aln.query_sequence))
            elif c_matches / (c_matches + g_matches) <= args.min_repeat_ratio:
                g_strand += 1
                aln.set_tag("XS", "G")
                telomeric.write(aln)
                telo_lengths.append(len(aln.query_sequence))
        else:
            continue

    print("Number of Telomeric Reads: {}".format(g_strand+c_strand))
    print("\tNumber of C Strand Telomeres: {}".format(c_strand))
    print("\tNumber of G Strand Telomeres: {}".format(g_strand))    
    print("Telomeric Read GB: {}".format(sum(telo_lengths)))
    print("Telomeric N50: {}".format(calculate_N50(telo_lengths)))

                
def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--minimum_read_length", required = True)
    parser.add_argument("--min_repeat_ratio", required=True)
    parser.add_argument("--telomeric", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
