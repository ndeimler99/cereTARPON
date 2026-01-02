#!/usr/bin/env python3
from collections import Counter

from bokeh.models import ColumnDataSource, FactorRange, Whisker, Range1d, HoverTool, Segment
#from bokeh.models import VSpan
import bokeh.transform as bkt
import numpy as np
import pandas as pd
from seaborn.categorical import _BarPlotter

from ezcharts.plots import BokehPlot
from seaborn._statistics import Histogram
from scipy.stats import gaussian_kde
from itertools import cycle
from seaborn.relational import _ScatterPlotter
import numbers
import math
from bokeh.models import ColumnDataSource
from bokeh.palettes import Category10
from bokeh.plotting import figure

def plot_genome_distribution(source_dictionary, index, reference, sample_name):
    chroms = list(source_dictionary.keys())
    x_positions = list(range(len(chroms)))  # numeric positions for x-axis

    # Sample and reference values
    sample_values = [source_dictionary[c][index] for c in chroms]
    reference_values = [reference[c] for c in chroms]

    source = ColumnDataSource(data=dict(
        x=x_positions,
        chrom=chroms,
        value=sample_values,
        reference=reference_values
    ))

    plt = BokehPlot(
        x_range=( -0.5, len(chroms)-0.5 ),  # numeric x-range
        height=400,
        width=700,
        y_axis_label="Percentage"
    )

    p = plt._fig

    # Sample bars
    bars = p.vbar(
        x='x',
        top='value',
        width=0.8,
        source=source,
        color='steelblue',
        legend_label=sample_name
    )

    # Horizontal reference lines
    for xi, y_ref in zip(x_positions, reference_values):
        p.line([xi - 0.4, xi + 0.4], [y_ref, y_ref], line_color="red", line_width=2)
    p.tools = [t for t in p.tools if not isinstance(t, HoverTool)]

    # Hover tool
    hover = HoverTool(
        renderers=[bars],
        tooltips=[
            ("Chromosome", "@chrom"),
            (f"{sample_name} Percentage", "@value{0.00}%"),
            ("Reference", "@reference{0.00}%")
        ]
    )
    p.add_tools(hover)

    # X-axis labels
    p.xaxis.ticker = x_positions
    p.xaxis.major_label_overrides = {i: chrom for i, chrom in enumerate(chroms)}

    p.y_range.start = 0
    p.legend.location = "top_right"
    p.xgrid.grid_line_color = None

    return plt
    

def telo_length_histogram(df, column):

    values = list(df[column])
    hist, edges = np.histogram(values, bins=[i*10 for i in range(0, int(max(df[column]/10)) + 1)])

    source = ColumnDataSource({
        'top': hist,
        'left': edges[:-1],
        'right': edges[1:]
    })

    plt = BokehPlot(
        height=400,
        width=700,
        title="Histogram Example",
        x_axis_label="Telomere Length",
        y_axis_label="Number of Sequences"
    )

    p = plt._fig

    r = p.quad(
        bottom=0,
        top='top',
        left='left',
        right='right',
        source=source,
        fill_color="skyblue",
        line_color="black"
    )
    p.tools = [t for t in p.tools if not isinstance(t, HoverTool)]

    hover = HoverTool(
        renderers=[r],
        tooltips=[
            ("Bin", "@left – @right"),
            ("Count", "@top")
        ]
    )

    p.add_tools(hover)

    return plt


def single_sample_boxplot(df, column, sample):

    values = list(df[column])
    q1 = np.percentile(values, 25)
    q2 = np.percentile(values, 50)
    q3 = np.percentile(values, 75)
    iqr = q3 - q1
    upper = min(q3 + 1.5 * iqr, max(values))
    lower = max(q1 - 1.5 * iqr, min(values))

    source = ColumnDataSource({
        "samples": [sample],
        "q1": [q1],
        "q2": [q2],
        "q3": [q3],
        "upper": [upper],
        "lower": [lower]
    })

    plt = BokehPlot(
        x_range=[sample],
        height=400,
        width=700,
        x_axis_label=sample,
        y_axis_label="Telomere Length"
    )

    p = plt._fig
    p.tools = [t for t in p.tools if not isinstance(t, HoverTool)]

    box = p.vbar(
        x="samples",
        width=0.7,
        bottom="q1",
        top="q3",
        fill_color="skyblue",
        line_color="black",
        source=source
    )

    # whiskers
    p.segment("samples", "upper", "samples", "q3", line_color="black", source=source)
    p.segment("samples", "lower", "samples", "q1", line_color="black", source=source)

    # median line
    p.segment(
        "samples", "q2", "samples", "q2",
        line_color="black",
        line_width=3,
        source=source
    )

    hover = HoverTool(
        renderers=[box],
        tooltips=[
            ("Sample", "@samples"),
            ("Q1", "@q1"),
            ("Median", "@q2"),
            ("Q3", "@q3"),
            ("Lower whisker", "@lower"),
            ("Upper whisker", "@upper")
        ]
    )
    p.add_tools(hover)

    return plt



def multi_sample_boxplot(source_dict):
    """
    Create a boxplot for multiple samples.
    
    df: pandas DataFrame
    column: numeric column to plot
    sample_col: column in df with sample names
    """
    
    samples = source_dict.keys()
    
    # Compute boxplot statistics per sample
    data = {
        "samples": [],
        "q1": [],
        "q2": [],
        "q3": [],
        "upper": [],
        "lower": []
    }
    
    for s in samples:
        vals = source_dict[s]
        if len(vals) == 0:
            continue
        q1, q2, q3 = np.percentile(vals, [25, 50, 75])
        iqr = q3 - q1
        upper = min(q3 + 1.5 * iqr, max(vals))
        lower = max(q1 - 1.5 * iqr, max(vals))
        
        data["samples"].append(s)
        data["q1"].append(q1)
        data["q2"].append(q2)
        data["q3"].append(q3)
        data["upper"].append(upper)
        data["lower"].append(lower)
    
    source = ColumnDataSource(data)
    
    # Create BokehPlot
    plt = BokehPlot(
        x_range=list(samples),
        height=400,
        width=700,
        x_axis_label="Sample",
        y_axis_label="Telomere Length"
    )
    p = plt._fig
    p.tools = [t for t in p.tools if not isinstance(t, HoverTool)]

    # Draw boxes
    box = p.vbar(
        x="samples",
        width=0.7,
        bottom="q1",
        top="q3",
        fill_color="skyblue",
        line_color="black",
        source=source
    )
    
    # Whiskers
    p.segment("samples", "upper", "samples", "q3", line_color="black", source=source)
    p.segment("samples", "lower", "samples", "q1", line_color="black", source=source)
    
    # Median line
    p.segment("samples", "q2", "samples", "q2", line_color="black", line_width=3, source=source)
    
    # Hover tool
    hover = HoverTool(
        renderers=[box],
        tooltips=[
            ("Sample", "@samples"),
            ("Q1", "@q1"),
            ("Median", "@q2"),
            ("Q3", "@q3"),
            ("Lower whisker", "@lower"),
            ("Upper whisker", "@upper")
        ]
    )
    p.add_tools(hover)
    
    return plt