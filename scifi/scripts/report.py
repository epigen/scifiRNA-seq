#! /bin/env python

import os
import sys
from argparse import ArgumentParser
import time

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


from src.scifi_utils import (
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
        help="Additional attributes to plot metrics grouped by."
        " A comma-delimited list."
        " Several can be given. Example: 'donor_id,sex,plate'.",
    )
    parser.add_argument(
        "--plate-column",
        dest="plate_column",
        default="plate",
        help="Name of column in input metric files"
        "containing r1 plate information. "
        "Defaults to 'plate'.",
    )
    parser.add_argument(
        "--well-column",
        dest="well_column",
        default="plate_well",
        help="Name of column in input metric files"
        "containing r1 well plate information. "
        "Defaults to 'plate_well'.",
    )
    parser.add_argument(
        "--droplet-column",
        dest="droplet_column",
        default="r2",
        help="Name of column in input metric files"
        "containing r2 (droplet) information. "
        "Defaults to 'r2'.",
    )
    parser.add_argument(
        "--nrows",
        dest="nrows",
        help="An optional number of rows of input to process. "
        "Useful for testing.",
    )
    parser.add_argument(
        "--species-mixture",
        dest="species_mixture",
        action="store_true",
        help="Whether the experiment is a species mixture. "
        "Default is false.",
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
    parser.add_argument(
        "--only-matching-barcodes",
        dest="only_matching_barcodes",
        action="store_true",
        help="Whether to use only barcodes matching the reference."
        "Requires a valid value for the '--r2-barcodes' option."
        "Default is false.",
    )
    default = os.path.join("metadata", "737K-cratac-v1.reverse_complement.csv")
    parser.add_argument(
        "--r2-barcodes", dest="r2_barcodes", default=default,
        help="Whilelist file with r2 barcodes."
             f"Defaults to '{default}.")
    choices = ["original", "reverse_complement"]
    parser.add_argument(
        "--barcode-orientation", dest="barcode_orientation",
        choices=choices, default=choices[0],
        help="Which orientation the r2 barcodes should be read as."
             f"One of '{', '.join(choices)}'."
    )

    # # Example:
    # args = parser.parse_args([
    #     "/scratch/lab_bock/shared/projects/sci-rna/data/PD190_humanmouse/PD190_humanmouse.metrics.csv.gz",
    #     "results/PD190_humanmouse.",
    #     "--plotting-attributes", "plate_well"])
    # args = parser.parse_args([
    #     "/scratch/lab_bock/shared/projects/sci-rna/data/PD190_sixlines/PD190_sixlines.metrics.csv.gz",
    #     "results/PD190_sixlines.",
    #     "--plotting-attributes", "plate_well"])
    # args = parser.parse_args([
    #     "/scratch/lab_bock/shared/projects/sci-rna/data/PD193_fivelines_383k/PD193_fivelines_383k.metrics.csv.gz",
    #     "results/PD193_fivelines_383k.",
    #     "--plotting-attributes", "plate_well,cell_line"])
    args = parser.parse_args()

    if args.plotting_attributes is None:
        args.plotting_attributes = list()
    else:
        args.plotting_attributes = args.plotting_attributes.split(",")
    return args


def main():
    global args
    args = set_args(parse_args())

    metrics = load_metrics(args.metric_files, nrows=args.nrows)

    # Filter r2 barcodes if required
    if args.only_matching_barcodes:
        # read text files
        try:
            r2_barcodes = pd.read_csv(args.r2_barcodes)
        except (OSError, IOError, FileNotFoundError):
            msg = (
                "Option --only-matching-barcodes given but value of "
                f"--r2-barcodes could not be read: '{args.r2_barcodes}'")
            raise FileNotFoundError(msg)
        metrics = get_exact_matches(
            metrics,
            barcodes=['r2'],
            whitelists=[r2_barcodes[args.barcode_orientation]],
            expected_cell_number=args.expected_cell_number)

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
    plot_efficiency(
        metrics,
        keys=['unique_fraction'],
        log_scale=False, tail=t, suffix="reads_vs_unique")

    if args.species_mixture:
        plot_species_mixing(metrics)

    # Well and droplet inspection
    if args.well_column in metrics.columns:
        well_metrics = gather_stats_per_well(metrics, seq_content=True)
        plot_well_stats(well_metrics, tail=None, suffix="")

        cells_per_droplet = (
            metrics
            # .tail(int(args.expected_cell_number / 2))
            .query("umi > 100")
            .groupby(args.droplet_column)
            [args.well_column].nunique()
            .sort_values())
        cells_per_droplet.name = "cells_per_droplet"
        cells_per_droplet_stats(cells_per_droplet)

        # Yield as a function of number of cells per droplet
        m = (
            metrics
            # .tail(int(args.expected_cell_number / 2))
            .query("umi > 100")
            .set_index("r2")
            .join(cells_per_droplet)
            .query("cells_per_droplet < 25"))
        m.to_csv(args.output_prefix + 'cells_per_droplet_stats.csv')

        attrs = [
            ('read', True),
            ('umi', True),
            ('gene', True),
            ('unique_fraction', False)]
        fig, axis = plt.subplots(len(attrs), 1, figsize=(9, 3 * len(attrs)))
        for i, (attr, log) in enumerate(attrs):
            sns.violinplot(m['cells_per_droplet'], m[attr], ax=axis[i])
            if log:
                axis[i].set_yscale("log")
        fig.savefig(
            args.output_prefix + f"cells_per_droplet.packaging_vs_yield.violinplot.svg",
            dpi=300, bbox_inches="tight")

        fig, axis = plt.subplots(len(attrs), 1, figsize=(9, 3 * len(attrs)))
        for i, (attr, log) in enumerate(attrs):
            sns.boxplot(m['cells_per_droplet'], m[attr], fliersize=0, ax=axis[i])
            if log:
                axis[i].set_yscale("log")
        fig.savefig(
            args.output_prefix + f"cells_per_droplet.packaging_vs_yield.boxplot.svg",
            dpi=300, bbox_inches="tight")

        if args.species_mixture:
            # Species mixing plot based on round2 only
            r2_metrics = metrics.groupby(args.droplet_column).agg(
                {'human': np.sum, 'mouse': np.sum, 'total': np.sum,
                 'human_norm': np.sum, 'total_norm': np.sum,
                 'sp_ratio': np.mean, 'sp_ratio_norm': np.mean})
            plot_species_mixing(r2_metrics, suffix="only_r2")

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
