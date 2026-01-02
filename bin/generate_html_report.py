#!/usr/bin/env python3

import sys
import os
import pandas as pd
#import plotly.express as px
#import plotly.graph_objects as go
import math
#from jinja2 import Template
import argparse
import json

#from bokeh.models import HoverTool
from dominate import tags as html_tags
from dominate.tags import em, p, style
import ezcharts as ezc
from ezcharts.layout.base import IStyles, Snippet
from ezcharts.components.ezchart import EZChart
from ezcharts.components.fastcat import SeqSummary
#from ezcharts.components.reports.labs import LabsReport
from extra.report import LabsReport
#from report import BasicReport
from ezcharts.layout.snippets import DataTable
from ezcharts.layout.snippets import Grid
from ezcharts.layout.snippets import Tabs
import pandas as pd
from ezcharts.components.reports.comp import ComponentReport
from bokeh.models import HoverTool, ColumnDataSource, Segment
import numpy as np
from bokeh.plotting import figure
import bokeh.palettes
from ezcharts.layout.util import css, render_template
from typing import Optional
import extra.report_utils as report_utils

THEME = 'epi2melabs'

def pass_fail_sample(row):

    if row["Number_of_Reads"] >= args.minimum_read_count:
        return "PASS"
    else:
        return "FAIL"
    
def get_nextflow_attributes(attribute_file):

    attribute_file = open(attribute_file)
    attribute = attribute_file.read()
    attribute = json.loads(attribute)
    return attribute

def main(args):

    args.minimum_read_count = int(args.minimum_read_count)

    params = get_nextflow_attributes(args.params)
    manifest = get_nextflow_attributes(args.manifest)

    versions = pd.read_table(args.versions, sep=",", header=None)
    versions.set_axis(["Software", "Version"], axis=1)


    alignment_results = {}
    for i in args.alignment_files:

        label = i.split(".")[0].split("/")[-1]
        alignment_results[label] = {}

        linecount = 0
        with open(i, "r") as fh:
            for line in fh:
                if linecount == 0:
                    linecount += 1
                    continue
                line = line.strip().split()
                alignment_results[label][line[0]] = [float(line[1]), float(line[2])]

    index_dict = {}
    with open(args.fsa_idx, "r") as fh:
        for line in fh:
            line = line.strip().split("\t")
            index_dict[line[0]] = int(line[1])
    
    total_length = sum(index_dict.values())
    for chrom in index_dict:
        index_dict[chrom] = index_dict[chrom]/total_length * 100

    sample_dict = {}
    for i in args.stats_files:
       sample_dict[i.split(".")[0]] = i

    report = LabsReport(
        f"Report for: {params['run_name']}", args.workflow_name,
        args.params, args.versions, manifest["version"], args.manifest)
    
    with report.add_section("Alignment Results", "Alignment Results"):
        tabs = Tabs()
      
        for sample in alignment_results:

            with tabs.add_tab(sample):
                new_tabs = Tabs()
                with new_tabs.add_tab("Read Count"):
                    read_count = report_utils.plot_genome_distribution(alignment_results[sample], 0, index_dict, sample)
                    EZChart(read_count)
                with new_tabs.add_tab("Nucleotide Count"):
                    nucl_count = report_utils.plot_genome_distribution(alignment_results[sample], 1, index_dict, sample)
                    EZChart(nucl_count)

    with report.add_section("Sample Analysis", "Sample Analysis"):
        tabs = Tabs()

        combined_dict = {}
        for sample in sample_dict:
            with tabs.add_tab(sample):
                df = pd.read_table(sample_dict[sample], sep="\t")
                combined_dict[sample] = list(df["telo_length"])
                telo_length_histo = report_utils.telo_length_histogram(df, "telo_length")
                telo_length_boxplot = report_utils.single_sample_boxplot(df, "telo_length", sample)
                summary_stats = pd.DataFrame({"Statistic":["Number of Telomeres", "Mean Telomere Length", "Median Telomere Length", 
                                                "Mean Read Length", "Median Read Length",
                                                "Mean Read Quality", "Median Read Quality"],
                                        "Value":[len(df["telo_length"]), np.mean(df["telo_length"]), np.median(df["telo_length"]),
                                                np.mean(df["read_length"]), np.median(df["read_length"]),
                                                np.mean(df["read_quality"]), np.median(df["read_quality"])]})
                with Grid(columns=3):
                        EZChart(telo_length_histo, THEME)
                        EZChart(telo_length_boxplot, THEME)
                        DataTable.from_pandas(summary_stats, use_index=False)

    with report.add_section("Sample Comparison", "Sample Comparison"):
        tabs = Tabs()

        with tabs.add_tab("Telomere Length Comparison Boxplot"):
            telo_plot = report_utils.multi_sample_boxplot(combined_dict)
            EZChart(telo_plot,THEME)

        with tabs.add_tab("Telomere Stats Comparison"):
            summary_stats = pd.DataFrame({'Sample':[i for i in combined_dict],
                                            'Number of Telomeres':[len(combined_dict[i]) for i in combined_dict],
                                            'Mean Telomere Length':[np.mean(combined_dict[i]) for i in combined_dict],
                                            'Standard Deviation Telo Length':[np.std(combined_dict[i]) for i in combined_dict],
                                            'Median Telomere Length':[np.median(combined_dict[i]) for i in combined_dict],
                                        })
            DataTable.from_pandas(summary_stats, use_index=False)
            
    report.write(args.report)

def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", required=True)
    parser.add_argument("--workflow_name", required=True, help="The name of the workflow.") # works
    parser.add_argument("--params", required=True) #works
    parser.add_argument("--versions", required=True) #works  
    parser.add_argument("--manifest", required=True) #works
    parser.add_argument("--minimum_read_count", required=True)
    parser.add_argument("--fsa_idx", required=True)
    parser.add_argument("--alignment_files", nargs='+', required=True) #works but only for simplex, not multiplex tested yet
    parser.add_argument("--stats_files", nargs='+', required=True) #works but only for simplex, not multiplex tested yet
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)

