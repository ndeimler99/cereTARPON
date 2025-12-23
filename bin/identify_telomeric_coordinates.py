#!/usr/bin/env python3 

import argparse
import pysam
import regex

def identify_telo_end(sequence, end_sequence):
    
    seq = sequence[-100:]
    matches = list(regex.finditer(r'%s{e<=1}' % end_sequence, seq))
    if len(matches) > 0:
        if matches[0].span()[0] > 50:
            return len(sequence) - 100 + matches[0].span()[0]
        else:
            return None
    else:
        return None

def identify_telo_start(sequence, sliding_window, interval, composition):
    
    telo_start = None
    for i in range(0, len(sequence)-sliding_window, interval):
        seq = sequence[i:i+sliding_window]
        matches = list(regex.finditer(r'T(G){1,3}',seq))
        span = 0
        for match in matches:
            span += (match.span()[1] - match.span()[0])
        comp = span/len(seq)

        if comp >= composition:
            if telo_start is None:
                telo_start = matches[0].span()[0] + i
        else:
            if telo_start is not None:
                telo_start = None

    return telo_start


def main(args):
    
    args.sliding_window = int(args.sliding_window)
    args.interval = int(args.interval)
    args.composition = float(args.composition)

    aln_file = pysam.AlignmentFile(args.input_file, "r", check_sq=False)
    out_file = pysam.AlignmentFile(args.output_file, "w", template=aln_file)

    end_identified = 0
    reads = 0
    g_end = 0
    c_end = 0
    start_identified = 0
    c_start = 0
    g_start = 0

    for aln in aln_file.fetch(until_eof=True):

        reads += 1
        telo_end = identify_telo_end(aln.query_sequence, args.end_sequence)
        
        if telo_end is not None:
            if aln.get_tag("XS") == "C":
                c_end += 1
            else:
                g_end += 1
            end_identified += 1

            telo_start = identify_telo_start(aln.query_sequence[0:telo_end], args.sliding_window, args.interval, args.composition)
            if telo_start is not None:
                if aln.get_tag("XS") == "C":
                    c_start += 1
                else:
                    g_start += 1
                start_identified += 1
                
                aln.set_tag("TE", telo_end)
                aln.set_tag("TS", telo_start)
                out_file.write(aln)
        else:
            pass
        
            
    print("Total Read Count: {}".format(reads))
    print("Adaptor Identified: {}, {}%".format(end_identified, end_identified/reads * 100))
    print("\tG-Strand: {}".format(g_end))
    print("\tC-Strand: {}".format(c_end))
    print("Start Identified: {}, {}%".format(start_identified, start_identified/reads * 100))
    print("\tG-Strand: {}".format(g_start))
    print("\tC-Strand: {}".format(c_start))

                
def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--output_file", required=True)
    parser.add_argument("--end_sequence", required=True)
    parser.add_argument("--sliding_window", required=True)
    parser.add_argument("--interval", required=True)
    parser.add_argument("--composition", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
