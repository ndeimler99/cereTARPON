#!/usr/bin/env python3 

import argparse
import cairo
import numpy as np

def meshclust_clusters(mesh_file):
    cluster_dict = {}
    with open(mesh_file, 'r') as fh:
        for line in fh:
            line = line.strip().split()
            if line == []:
                continue
            if line[0] in cluster_dict:
                cluster_dict[line[0]].append(line[1].strip(">"))
            else:
                cluster_dict[line[0]] = [line[1].strip(">")]
    return cluster_dict

def stats_results(stats_file):
    stats_dict = {}
    linecount = 0
    with open(stats_file, "r") as stats_fh:
        for line in stats_fh:
            if linecount == 0:
                linecount += 1
                continue
            line = line.strip().split()
            stats_dict[line[0]] = {'telo_length':int(line[6]), 'read_length':int(line[4])}
    return stats_dict


def blast_results(blast_file, max_len, TAS_perc, TAS2_perc, TAS_fraction):
    # max_len is x axis length

    # 6 0-qseqid 1-sseqid 2-pident 3-length 4-slen 5-mismatch 6-gapopen 7-qstart 8-qend 9-sstart 10-send evalue bitscore qlen

    blast_dict = {}
    with open(blast_file, 'r') as blast_fh:
        for line in blast_fh:
            line = line.strip().split()
            if line[0] not in blast_dict:
                blast_dict[line[0]] = []
            
            if int(line[13]) - int(line[8]) > max_len:
                #line[8] = max_len
                continue
            if int(line[13]) - int(line[7]) > max_len: 
                line[7] = int(line[13]) - max_len
            

            if "TAS" in line[1] or "TERRA" in line[1]:
                if line[1] == 'TAS2' and float(line[2]) > TAS2_perc:
                    blast_dict[line[0]].append((line[1],int(line[13]) - int(line[8]), int(line[13]) - int(line[7])))
                else:
                    if float(line[2]) > TAS_perc and abs(int(line[10]) - int(line[9])) >= int(line[4]) * TAS_fraction:
                        blast_dict[line[0]].append((line[1], int(line[13]) - int(line[8]), int(line[13]) - int(line[7])))
                        #print(line[1])
                # make sure it is full length rRNA
            else:
                blast_dict[line[0]].append((line[1], int(line[13]) - int(line[8]), int(line[13]) - int(line[7])))  

    return blast_dict

def get_color_dict(color_dict_file):
    color_dict = {}
    with open(color_dict_file, "r") as fh:
        for line in fh:
            line = line.strip().split()
            color_dict[line[0]] = [float(i) for i in line[1].split(",")]
    return color_dict

def flatten(xss):
    return [x for xs in xss for x in xs]

def main(args):

    args.x_axis_length = int(args.x_axis_length)
    args.TAS_perc = float(args.TAS_perc)
    args.TAS2_perc = float(args.TAS2_perc)
    args.TAS_fraction = float(args.TAS_fraction)
    args.image_width_px  = int(args.image_width_px)
    args.image_height_px = int(args.image_height_px)
    args.minimum_cluster_size = int(args.minimum_cluster_size)
    args.color_dict = get_color_dict(args.color_dict)    


    stats_dict = stats_results(args.stats_file)
    if args.x_axis_length == 0:
        for read in stats_dict:
            if stats_dict[read]["read_length"] > args.x_axis_length:
                args.x_axis_length = stats_dict[read]["read_length"]
    else:
        for read in stats_dict:
            if stats_dict[read]["read_length"] > args.x_axis_length:
                stats_dict[read]["read_length"] = args.x_axis_length

    blast_dict = blast_results(args.blast_file, args.x_axis_length, args.TAS_perc, args.TAS2_perc, args.TAS_fraction)
    cluster_dict = meshclust_clusters(args.cluster_file)

    
    surface = cairo.ImageSurface(cairo.FORMAT_RGB24,
                                 args.image_width_px,
                                 args.image_height_px)

    ctx = cairo.Context(surface)
    ctx.rectangle(0, 0, args.image_width_px, args.image_height_px)
    ctx.set_source_rgb(1,1,1)
    ctx.fill()

    mean_sd_dict = {}
    for cluster in cluster_dict:
        telo_lengths = []
        for read in cluster_dict[cluster]:
            telo_lengths.append(stats_dict[read]["telo_length"])
        mean_sd_dict[cluster] = {"mean":np.mean(telo_lengths), "sd":np.std(telo_lengths), "q2": np.quantile(telo_lengths, 0.5), "q3":np.quantile(telo_lengths, 0.75)}
        mean_sd_dict[cluster]["iqr"] = (mean_sd_dict[cluster]["q3"] - mean_sd_dict[cluster]["q2"]) * 3

    telo_lengths_to_plot = []
    for cluster in cluster_dict:
        if len(cluster_dict[cluster]) < args.minimum_cluster_size:
            continue
        for read in cluster_dict[cluster]:
#            if stats_dict[read]["telo_length"] > mean_sd_dict[cluster]["mean"] + mean_sd_dict[cluster]["iqr"]:
#                continue
            telo_lengths_to_plot.append(stats_dict[read]["telo_length"])

    #max_telo_length = max([stats_dict[read]["telo_length"] for read in stats_dict if read in flatten(list(cluster_dict.values())) and stats_dict[read]["telo_length"]])
    max_telo_length = 6500
    telo_plot_width = int(args.image_width_px) * 0.4
    x_offset_left=int(args.image_width_px * 0.05) + telo_plot_width + 50
    x_offset_right=int(args.image_width_px * 0.05)
    y_offset_top=int(args.image_width_px * 0.05)
    y_offset_bottom=int(args.image_width_px * 0.05)

    # ctx.set_source_rgb(0,0,0)
    # ctx.rectangle(y_offset_top, args.image_height_px - y_offset_bottom, args.image_width_px-x_offset_right*2, 1)
    # ctx.fill()

    #read_count = sum([len(cluster_dict[cluster]) for cluster in cluster_dict if len(cluster_dict[cluster]) >= args.minimum_cluster_size])
    read_count = len(telo_lengths_to_plot)
    cluster_gap_size = int((len(cluster_dict) - 1) * args.image_height_px * 0.02)
    sequence_height = (args.image_height_px-y_offset_bottom-y_offset_top-cluster_gap_size) / read_count
    nucl_width = (args.image_width_px - x_offset_left - x_offset_right) / args.x_axis_length
    telo_nucl_width = (telo_plot_width / max_telo_length)

    ctx.set_source_rgb(0,0,0)
    ctx.rectangle(x_offset_left, args.image_height_px-y_offset_bottom, args.image_width_px-x_offset_left - x_offset_right, 3)
    ctx.fill()

    ctx.set_font_size(20)
    ctx.select_font_face("Courier",
                          cairo.FONT_SLANT_NORMAL,
                          cairo.FONT_WEIGHT_NORMAL)

     # draw sequences
    seq_offset=0
    for i,cluster in enumerate(cluster_dict):
        if len(cluster_dict[cluster]) < args.minimum_cluster_size:
            continue
        #print(seq_offset)
        #a += 0.2
        #ctx.set_source_rgb(a, 0 , 0)
        telo_lengths = []
        j = 0
        for k,read in enumerate(cluster_dict[cluster]):

            if read not in stats_dict:
                continue

#            if stats_dict[read]["telo_length"] > mean_sd_dict[cluster]["mean"] + mean_sd_dict[cluster]["iqr"]:
#                continue
            
            j = j + 1
            ctx.set_source_rgb(0.4,0.4,0.4)
            ctx.rectangle(x_offset_left, args.image_height_px - y_offset_bottom - ((j + 1) * sequence_height) - seq_offset, stats_dict[read]["read_length"] * nucl_width, sequence_height/2)
            ctx.fill()

            if read in blast_dict:
                for item in blast_dict[read]:
                    ctx.set_source_rgb(args.color_dict[item[0]][0],args.color_dict[item[0]][1],args.color_dict[item[0]][2])
                    ctx.rectangle(x_offset_left + item[1]*nucl_width, args.image_height_px - y_offset_bottom - ((j + 1)* sequence_height) - seq_offset, item[2]*nucl_width-item[1] * nucl_width, sequence_height)
                    ctx.fill()
                
            if stats_dict[read]["telo_length"] is not None:
                ctx.set_source_rgb(0.1,0.1,0.1)
                ctx.rectangle(x_offset_left - (stats_dict[read]["telo_length"] * telo_nucl_width), args.image_height_px - y_offset_bottom - ((j + 1)* sequence_height) - seq_offset, stats_dict[read]["telo_length"] * telo_nucl_width, sequence_height)
                ctx.fill()
                telo_lengths.append(stats_dict[read]["telo_length"])

        # ctx.set_source_rgb(1,0,0)
        # ctx.rectangle(x_offset_left - sum(telo_lengths) * telo_nucl_width/len(telo_lengths), args.image_height_px - y_offset_bottom - ((j + 1)* sequence_height) - seq_offset, 5, j*sequence_height)
        # ctx.rectangle(x_offset_left - sum(telo_lengths)* telo_nucl_width/len(telo_lengths) - np.std(telo_lengths) * telo_nucl_width, args.image_height_px - y_offset_bottom - ((j + 1)* sequence_height) - seq_offset + j*sequence_height/2, np.std(telo_lengths)*2 * telo_nucl_width, 2)
        # ctx.rectangle(x_offset_left - sum(telo_lengths)* telo_nucl_width/len(telo_lengths) - np.std(telo_lengths) * telo_nucl_width, args.image_height_px - y_offset_bottom - ((j + 1)* sequence_height) - seq_offset + j*sequence_height*0.25, 3, j*sequence_height*0.5)
        # ctx.rectangle(x_offset_left - sum(telo_lengths)* telo_nucl_width/len(telo_lengths) + np.std(telo_lengths) * telo_nucl_width, args.image_height_px - y_offset_bottom - ((j + 1)* sequence_height) - seq_offset + j*sequence_height*0.25, 3, j*sequence_height*0.5)   
        # ctx.fill()

        # ctx.set_source_rgb(0,1,0)
        # ctx.rectangle(x_offset_left - np.quantile(telo_lengths, 0.5) * telo_nucl_width, args.image_height_px - y_offset_bottom - ((j + 1)* sequence_height) - seq_offset, 5, j*sequence_height)
        # ctx.rectangle(x_offset_left - np.quantile(telo_lengths, 0.25) * telo_nucl_width, args.image_height_px - y_offset_bottom - ((j + 1)* sequence_height) - seq_offset + j*sequence_height*0.25, 5, j*sequence_height*0.5)
        # ctx.rectangle(x_offset_left - np.quantile(telo_lengths, 0.75) * telo_nucl_width, args.image_height_px - y_offset_bottom - ((j + 1)* sequence_height) - seq_offset + j*sequence_height*0.25, 5, j*sequence_height*0.5)
        # ctx.rectangle(x_offset_left - np.quantile(telo_lengths, 0.75) * telo_nucl_width, args.image_height_px - y_offset_bottom - ((j + 1)* sequence_height) - seq_offset + j*sequence_height / 2, np.quantile(telo_lengths, 0.75)*telo_nucl_width - np.quantile(telo_lengths, 0.25) * telo_nucl_width, 2)
        # ctx.fill()


        ctx.set_source_rgb(0,0,0)
        ctx.move_to(10, args.image_height_px - y_offset_bottom - ((j/2 + 1)* sequence_height) - seq_offset + sequence_height * 0.5)
        #ctx.text_path('Clust-{}'.format(cluster))
        ctx.fill()
        seq_offset += j*sequence_height
        seq_offset += int(args.image_height_px * 0.02)
        
    ctx.set_source_rgb(0,0,0)
    for i in range(0, args.x_axis_length+100, 1000):
        if i % 10000 == 0:
            a,b,width,height,c,d = ctx.text_extents('{}'.format(i))
            ctx.move_to(x_offset_left + i*nucl_width - width/2 - 2.5, args.image_height_px - y_offset_bottom + 25)
            ctx.text_path('{}'.format(i))
            ctx.fill()
        if i % 2000 == 0:
            ctx.rectangle(x_offset_left + i*nucl_width - 2.5, args.image_height_px-y_offset_bottom, 5, 5)
            ctx.fill()
    ctx.set_source_rgb(0,0,0)
    ctx.rectangle(x_offset_right, args.image_height_px-y_offset_bottom, x_offset_left-x_offset_right, 3)
    ctx.fill()
    if max_telo_length < 1000:
        val = 100
    elif max_telo_length < 2000:
        val = 500
    elif max_telo_length < 5000:
        val = 1000
    else:
        val = 2000
    for i in range(val, max_telo_length, val):
        a,b,width,height,c,d = ctx.text_extents('{}'.format(i))
        ctx.rectangle(x_offset_left - i*telo_nucl_width - 2.5, args.image_height_px-y_offset_bottom, 5, 5)
        ctx.fill()
        ctx.move_to(x_offset_left - i*telo_nucl_width - width/2 - 2.5, args.image_height_px-y_offset_bottom + 25)
        ctx.text_path('{}'.format(i))
        ctx.fill()
        
    for i,seq in enumerate(args.color_dict):
        ctx.set_source_rgb(args.color_dict[seq][0], args.color_dict[seq][1], args.color_dict[seq][2])
        ctx.rectangle(args.image_width_px/2 - 150*(len(args.color_dict)/2) + i*150, 25, 25, 25)
        ctx.fill()
        ctx.set_source_rgb(0,0,0)
        ctx.move_to(args.image_width_px/2 - 150*(len(args.color_dict)/2) + i*150 + 30, 50-25/4)
        ctx.text_path('{}'.format(seq))
        ctx.fill()

    surface.write_to_png(args.plot_file_name)
    
   

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats_file", required=True)
    parser.add_argument("--blast_file", required=True)
    parser.add_argument("--cluster_file", required=True)
    parser.add_argument("--plot_file_name", required=True)
    parser.add_argument("--x_axis_length", required=True)
    parser.add_argument("--TAS_perc", required=True)
    parser.add_argument("--TAS2_perc", required=True)
    parser.add_argument("--TAS_fraction", required=True)
    parser.add_argument("--image_width_px", required=True)
    parser.add_argument("--image_height_px", required=True)
    parser.add_argument("--minimum_cluster_size", required=True)
    parser.add_argument("--color_dict", required=True)
    #   default  {'TAS1':(0.84,0.10,0.37), 'TAS2':(0.1,0.53,0.89), 'TAS3':(1,0.76, 0.03), '28S_rRNA':(0,0.30,0.25), '18S_rRNA':(0.98,0.07,0.92)}")
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
