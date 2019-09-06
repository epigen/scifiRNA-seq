#! /bin/env python

import os
import sys
from argparse import ArgumentParser
import time

import pandas as pd
import matplotlib
import seaborn as sns


from scifi_utils import (
    load_metrics,
    get_exact_matches,
    plot_efficiency,
    plot_metrics_lineplot,
    plot_species_mixing,
    gather_stats_per_well,
    plot_well_stats,
    cells_per_droplet_stats,
    write_gene_expression_matrix,
    set_args
)


def parse_args():
    parser = ArgumentParser()
    parser.add_argument(
        dest="metric_files",
        nargs="+",
        help="Input metric files from the output of sci-RNA.summarizer.py. "
        "Several can be given.",
    )
    parser.add_argument(
        dest="output_prefix",
        help="Absolute path prefix for output files. "
        "Example: 'results/SCIXX_experiment.'.",
    )
    parser.add_argument(
        "--plotting-attributes",
        dest="plotting_attributes",
        help="Additional attributes to plot metrics grouped by. A comma-delimited list."
        "Several can be given. Example: 'donor_id sex plate'.",
    )
    parser.add_argument(
        "--plate-column",
        dest="plate_column",
        default="plate",
        help="Name of column in input metric files containing r1 plate information. "
        "Defaults to 'plate'.",
    )
    parser.add_argument(
        "--well-column",
        dest="well_column",
        default="plate_well",
        help="Name of column in input metric files containing r1 well plate information. "
        "Defaults to 'plate_well'.",
    )
    parser.add_argument(
        "--droplet-column",
        dest="droplet_column",
        default="r2",
        help="Name of column in input metric files containing r2 (droplet) information. "
        "Defaults to 'r2'.",
    )
    parser.add_argument(
        "--nrows",
        dest="nrows",
        help="An optional number of rows of input to process. " "Useful for testing.",
    )
    parser.add_argument(
        "--species-mixture",
        dest="species_mixture",
        action="store_true",
        help="Whether the experiment is a species mixture. " "Default is false.",
    )
    parser.add_argument(
        "--expected-cell-number",
        dest="expected_cell_number",
        default=200000,
        type=int,
        help="Expected experiment yield number.",
    )
    parser.add_argument(
        "--save-intermediate",
        dest="save_intermediate",
        action="store_true",
        help="Whether to save intermediate pickle files.",
    )
    default = os.path.join("metadata", "737K-cratac-v1.reverse_complement.csv")
    parser.add_argument(
        "--r2-barcodes", dest="r2_barcodes", default=default,
        help="Whilelist file with r2 barcodes."
             f"Defaults to '{default}.")
    choices = ["original", "reverse_complement"]
    parser.add_argument(
        "--barcode-orientation", dest="barcode_orientation", choices=choices, default=choices[0],
        help="Which orientation the r2 barcodes should be read as."
             f"One of '{', '.join(choices)}'."
    )

    # # Example:
    # args = parser.parse_args([
    #     "/scratch/lab_bock/shared/projects/sci-rna/data/PD190_humanmouse/PD190_humanmouse.metrics.csv.gz",
    #     "results/PD190_humanmouse.",
    #     "--plotting-attributes", "plate_well"])
    args = parser.parse_args()

    if args.plotting_attributes is None:
        args.plotting_attributes = list()
    else:
        args.plotting_attributes = args.plotting_attributes.split(",")
    return args


def main():
    global args
    args = set_args(parse_args())

    # read text files
    try:
        r2_barcodes = pd.read_csv(args.r2_barcodes)
    except (OSError, IOError, FileNotFoundError):
        r2_barcodes = None
    metrics = load_metrics(args.metric_files, nrows=args.nrows)

    # # Plot
    # # # Loglog line plot of rank vs abundance
    t = args.expected_cell_number * 5
    plot_metrics_lineplot(metrics, tail=t)
    for attribute in args.plotting_attributes:
        plot_metrics_lineplot(metrics, tail=t, by_group=attribute)

    # # # Efficiency plot (UMI and gene yield in relation to reads per cell)
    t = int(args.expected_cell_number * 5)
    plot_efficiency(metrics, tail=t)
    plot_efficiency(metrics, tail=t, colour_by="unique_fraction")

    if args.species_mixture:
        plot_species_mixing(metrics)

    # Well and droplet inspection
    if args.well_column in metrics.columns:
        well_metrics = gather_stats_per_well(metrics, seq_content=True)
        plot_well_stats(well_metrics, tail=None, suffix="")

        cells_per_droplet = metrics.groupby(args.droplet_column)[
            args.well_column
        ].nunique()
        cells_per_droplet_stats(cells_per_droplet)

    # Barcode inspection
    if r2_barcodes is not None:
        metrics_filtered = get_exact_matches(
            metrics,
            barcodes=['r2'],
            whitelists=[r2_barcodes[args.barcode_orientation]],
            expected_cell_number=args.expected_cell_number)

        t = int(args.expected_cell_number * 5)
        plot_efficiency(metrics_filtered, suffix="filtered.", tail=t)
        plot_efficiency(metrics_filtered, suffix="filtered.", tail=t, colour_by="unique_fraction")

    # Write expression matrix in sparse format as H5AD
    if len(args.metric_files) == 1:
        print(f"# {time.asctime()} - Reading sparse expression matrix.")
        expr = pd.read_csv(args.metric_files[0].replace("metrics", "expression"))
        write_gene_expression_matrix(expr)

    print(f"# {time.asctime()} - Finished!")


if __name__ == "__main__":
    matplotlib.use("Agg")
    sns.set(context="paper", style="ticks", palette="colorblind", color_codes=True)
    matplotlib.rcParams["svg.fonttype"] = "none"
    # Don't use LaTeX for rendering
    matplotlib.rcParams["text.usetex"] = False
    sys.exit(main())
