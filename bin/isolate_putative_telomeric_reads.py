#!/usr/bin/env python3 

import argparse

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

    args.minimum_read_length = int(args.minimum_read_length)
    args.min_repeat_ratio = int(args.min_repeat_ratio)

    karyotype = get_karyotype_dict(args.karyotype_file)
    aln_file = pysam.AlignmentFile(args.alignment_file, "rb", check_sq=False)
    chimeric_file = pysam.AlignmentFile(args.chimeric_file, "w", template=aln_file)
    non_telomeric = pysam.AlignmentFile(args.non_telomeric, "w", template=aln_file)
    telomeric = pysam.AlignmentFile(args.telomeric, "w", template=aln_file)

    telo_reads = []
    lengths = []
    telo_lengths = []
    g_strand = 0
    c_strand = 0
    chimeras = 0
    bp = 0
    total_reads = 0

    for aln in aln_file:
        if aln.query_sequence is not None:
            bp += len(aln.query_sequence)
            total_reads += 1      

        if aln.is_unmapped or aln.reference_name == "ref|NC_001224|":
            non_telomeric.write(aln)
            continue
        
        if aln.reference_start < 5000 or aln.reference_end > karyo_dict[aln.reference_name] - 5000:
            if aln.query_name not in telo_reads:
                telo_reads.append(aln.query_name)
                if len(aln.query_sequence) < args.minimum_read_length:
                    non_telomeric.write(aln)
                    continue
                
                c_matches = len(list(regex.finditer(r'(C){1,3}A(C){1,3}A(C){1,3}A(C){1,3}A(C){1,3}A', aln.query_sequence[0:250], overlapped=True)))
                g_matches = len(list(regex.finditer(r'T(G){1,3}T(G){1,3}T(G){1,3}T(G){1,3}T(G){1,3}', aln.query_sequence[-250:], overlapped=True)))
                if c_matches + g_matches != 0:
                    if c_matches / (c_matches + g_matches) > args.min_repeat_ratio and c_matches / (c_matches + g_matches) < 1-args.min_repeat_ratio:
                        chimeric_file.write(aln)
                        chimeras += 1
                    elif c_matches / (c_matches + g_matches) >= 1-args.min_repeat_ratio:
                        c_strand += 1
                        lengths.append(len(aln.query_sequence))
                        q = aln.query_qualities
                        aln.query_sequence = rev_comp(aln.query_sequence)
                        aln.query_qualities = q[::-1]
                        aln.set_tag("XS", "C")
                        telomeric.write(aln)
                        telo_lengths.append(len(aln.query_sequence))
                    elif c_matches / (c_matches + g_matches) <= args.min_repeat_ratio:
                        g_strand += 1
                        lengths.append(len(aln.query_sequence))
                        aln.set_tag("XS", "G")
                        telomeric.write(aln)
                        telo_lengths.append(len(aln.query_sequence))
                else:
                    non_telomeric.write(aln)

    sequencing_N50 = calculate_N50(lengths)
    telo_N50 = calculate_N50(telo_lengths)

    print("Number of Reads: {}".format(total_reads))
    print("GB of Sequencing: {}".format(bp/1000000000))
    print("Number of Chimeric Reads: {}, {}".format(chimeras, chimeras/total_reads * 100))
    print("Number of Telomeric Reads: {}".format(c_strand+g_strand))
        print("Number of C Strand Telomeres: {}".format(c_strand))
        print("Number of G Strand Telomeres: {}".format(g_strand))
    print("Sequencing N50: {}".format(N50))
    print("Telomeric N50: {}".format(telo_N50))
                
def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment_file", required=True)
    parser.add_argument("--karyotype_file", required=True)
    parser.add_argument("--minimum_read_length", required = True)
    parser.add_argument("--non_telomeric", required=True)
    parser.add_argument("--chimeric_file", required=True)
    parser.add_argument("--telomeric", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
