#!/usr/bin/env python3 

import argparse
import pysam

def get_karyotype_dict(karyotype):
    chrom_lengths = {}
    with open(karyotype, "r") as karyo:
        for line in karyo:
            line = line.strip().split()
            chrom_lengths[line[0]] = int(line[1])
    return chrom_lengths

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

    args.subtelo_stretch = int(args.subtelo_stretch)

    karyotype = get_karyotype_dict(args.karyotype_file)
    aln_file = pysam.AlignmentFile(args.alignment_file, "rb", check_sq=False)

    telo_reads = []
    all_reads = []
    all_read_lengths = []

    for aln in aln_file:
        
        if aln.is_unmapped or aln.reference_name == "ref|NC_001224|" or aln.reference_name == "chrmt":
            continue

        if aln.query_sequence is not None:
            all_read_lengths.append(len(aln.query_sequence))
        
        all_reads.append(aln.query_name)

        if aln.reference_start < args.subtelo_stretch or aln.reference_end > karyotype[aln.reference_name] - args.subtelo_stretch:
            telo_reads.append(aln.query_name)
            
    with open(args.read_ids, "w") as out_fh:
        for read in set(telo_reads):
            out_fh.write(read + "\n")

    print("Number of Reads: {}".format(len(set(all_reads))))
    print("GB Mapped to Reference (not mt): {}".format(sum(all_read_lengths)))
    print("Sequencing N50: {}".format(calculate_N50(all_read_lengths)))
    print("Number of Reads in Terminal {}bp of Reference: {}".format(args.subtelo_stretch, len(set(telo_reads))))

                
def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment_file", required=True)
    parser.add_argument("--karyotype_file", required=True)
    parser.add_argument("--subtelo_stretch", required=True)
    parser.add_argument("--read_ids", required = True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
