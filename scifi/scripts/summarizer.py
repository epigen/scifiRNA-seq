#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script parses mapped, tagged BAM files, producing barcode-wise
statistics and their respective expression measurements.
"""

import os
import pickle
import sys
import time

from argparse import ArgumentParser
from glob import glob

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.sparse import csr_matrix
from scipy.special import factorial
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable

import pysam
import anndata as an
from joblib import Parallel, delayed
import multiprocessing


def parse_args(cli=None):
    parser = ArgumentParser()
    parser.add_argument(
        dest="input_files",
        nargs="+",
        help="Input BAM files to summarize. "
        "Can be several files or a regex that evaluates to several files.",
    )
    parser.add_argument(
        dest="output_prefix",
        help="Absolute path prefix for output files. "
        "Example: 'results/SCIXX_experiment.'.",
    )
    parser.add_argument(
        "--sample-name",
        dest="sample_name",
        help="Simple name to label output files. "
        "Defaults to same as 'output_prefix' option.",
    )
    parser.add_argument(
        "--cell-barcodes",
        dest="cell_barcodes",
        default=["r1", "r2"],
        nargs="+",
        help="Which barcodes to take into consideration to group molecules by. "
        "Defaults to 'r1 r2'.",
    )
    root = os.path.expanduser("~/projects/sci-rna")
    parser.add_argument(
        "--r1-annot",
        dest="r1_annotation_file",
        help="CSV file with annotations of r1 barcodes. ",
        required=True,
    )
    parser.add_argument(
        "--r1-attributes",
        dest="r1_attributes",
        type=str,
        help="Which r1 attributes to annotate cells with. A comma-separated list.",
    )
    default = os.path.join(
        root, "metadata", "737K-cratac-v1.reverse_complement.csv"
    )
    parser.add_argument(
        "--r2-barcodes",
        dest="r2_barcodes",
        default=default,
        help="Whilelist file with r2 barcodes." f"Defaults to '{default}.",
    )
    choices = ["original", "reverse_complement"]
    parser.add_argument(
        "--barcode-orientation",
        dest="barcode_orientation",
        choices=choices,
        default=choices[0],
        help="Which orientation the r2 barcodes should be read as."
        f" Defaults to '{choices[0]}'.",
    )
    default = ["read", "r2", "umi", "gene", "pos"]
    parser.add_argument(
        "--input-header",
        dest="input_header",
        default=default,
        help=f"Columns of input_files. Defautlts to '{''.join(default)}'.",
    )
    parser.add_argument(
        "--no-output-header",
        dest="no_output_header",
        action="store_true",
        help="Whether not to output a header. Defaults to not to.",
    )
    parser.add_argument(
        "--nrows",
        dest="nrows",
        type=int,
        default=int(1e10),
        help="Number of input rows to process. Useful for testing.",
    )
    parser.add_argument(
        "--only-summary",
        dest="only_summary",
        action="store_true",
        help="Only output summary, no plots.",
    )
    parser.add_argument(
        "--species-mixture",
        dest="species_mixture",
        action="store_true",
        help="Whether experiment is a species mixing.",
    )
    parser.add_argument(
        "--expected-cell-number",
        dest="expected_cell_number",
        default=200000,
        type=int,
        help="Number of expected cells. Only used if plotting.",
    )
    parser.add_argument(
        "--no-save-intermediate",
        dest="no_save_intermediate",
        action="store_true",
        help="Whether to skip saving of intermediate results.",
    )
    parser.add_argument(
        "--min-umi-output",
        dest="min_umi_output",
        default=1,
        type=int,
        help="Minimum UMIs per barcode to output.",
    )
    parser.add_argument(
        "--save-gene-expression",
        dest="save_gene_expression",
        action="store_true",
        help="Whether to create a gene expression matrix.",
    )
    parser.add_argument(
        "--correct-r1-barcodes",
        dest="correct_r1_barcodes",
        action="store_true",
        help="Whether set all r1 barcodes to the respective sequence from annotation.",
    )
    parser.add_argument(
        "--correct-r2-barcodes",
        dest="correct_r2_barcodes",
        action="store_true",
        help="Whether to use a mapping of existing barcodes to correct r2 barcodes.",
    )
    parser.add_argument(
        "--correct-r2-barcode-file",
        dest="correct_r2_barcode_file",
        help="File containing mapping between existing barcodes and correct barcodes.",
    )

    # # Example run:
    # cli = [
    #     "--sample-name", "SCI023_Tcell_D77_A01",
    #     "--r1-attributes", "plate", "plate_well", "donor_id", "sex",
    #     "--cell-barcodes", "r2",
    #     "--only-summary",
    #     "--no-save-intermediate",
    #     "--min-umi-output", "20",
    #     "--no-output-header",
    #     "/scratch/lab_bock/shared/projects/sci-rna/data/SCI024_Tcell_s/SCI023_Tcell_D77_A01/SCI023_Tcell_D77_A01.*.STAR.Aligned.out.bam.featureCounts.bam",
    #     "/scratch/lab_bock/shared/projects/sci-rna/data/SCI024_Tcell_s/SCI023_Tcell_D77_A01/SCI023_Tcell_D77_A01"]

    # cli = [
    #     "--sample-name", "PD193_humanmouse_765k_A01",
    #     "--r1-annot", "metadata/sciRNA-seq.PD193_humanmouse_765k.oligos_2019-09-20.csv",
    #     "--r1-attributes", "plate_well",
    #     "--cell-barcodes", "r2",
    #     "--species-mixture",
    #     "--only-summary",
    #     "--no-save-intermediate",
    #     "--min-umi-output", "3",
    #     "--expected-cell-number", "250326",
    #     "--no-output-header",
    #     "--save-gene-expression",
    #     "--correct-r1-barcodes",
    #     "--correct-r2-barcodes",
    #     "--correct-r2-barcode-file", "/scratch/lab_bock/shared/projects/sci-rna/data/PD193_humanmouse_765k/PD193_humanmouse_765k.fixed_barcodes.mapping.tsv",
    #     "/scratch/lab_bock/shared/projects/sci-rna/data/PD193_humanmouse_765k/PD193_humanmouse_765k_A01/PD193_humanmouse_765k_A01.ALL.STAR.Aligned.out.bam.featureCounts.bam",
    #     "/scratch/lab_bock/shared/projects/sci-rna/data/PD193_humanmouse_765k/PD193_humanmouse_765k_A01/PD193_humanmouse_765k_A01"]

    # cli = [
    #     "--sample-name", "PD190_humanmouse",
    #     "--r1-annot", "metadata/sciRNA-seq.PD190_sixlines.oligos_2019-09-05.csv",
    #     "--r1-attributes", "plate_well",
    #     "--cell-barcodes", "r2",
    #     "--only-summary",
    #     "--no-save-intermediate",
    #     "--min-umi-output", "20",
    #     "--no-output-header",
    #     "/scratch/lab_bock/shared/projects/sci-rna/data/PD190_humanmouse/PD190_humanmouse_*/*.STAR.Aligned.out.bam.featureCounts.bam",
    #     "/scratch/lab_bock/shared/projects/sci-rna/data/PD190_humanmouse/PD190_humanmouse"]

    # cli = [
    #     "--r1-annot", "/home/arendeiro/sci-rna/metadata/sciRNA-seq.PD193_fivelines_383k.oligos_2019-09-20.csv",
    #     "--r1-attributes", "plate_well,cell_line",
    #     "--cell-barcodes", "r2",
    #     "--only-summary",
    #     "--no-save-intermediate",
    #     "--min-umi-output", "3",
    #     "--expected-cell-number", "250326",
    #     "--no-output-header",
    #     "--save-gene-expression",
    #     "--correct-r1-barcodes",
    #     "--sample-name", "PD193_fivelines_383k_JurkatCas9TCRlib_A01",
    #     "/home/arendeiro/sci-rna/data/PD193_fivelines_383k/PD193_fivelines_383k_JurkatCas9TCRlib_A01/PD193_fivelines_383k_JurkatCas9TCRlib_A01.*.STAR.Aligned.out.bam.featureCounts.bam",
    #     "/home/arendeiro/sci-rna/data/PD193_fivelines_383k/PD193_fivelines_383k_JurkatCas9TCRlib_A01/PD193_fivelines_383k_JurkatCas9TCRlib_A01"]

    # cli = [
    #     "--r1-annot", "/home/arendeiro/sci-rna/metadata/sciRNA-seq.PD195-1_Tcells_765k-sample1_P7.oligos_2019-10-15.csv",
    #     "--r1-attributes", "donor_id,donor_sex,activation",
    #     "--cell-barcodes", "r2",
    #     "--only-summary",
    #     "--no-save-intermediate",
    #     "--min-umi-output", "3",
    #     "--expected-cell-number", "250326",
    #     "--no-output-header",
    #     "--save-gene-expression",
    #     "--correct-r1-barcodes",
    #     "--correct-r2-barcodes",
    #     "--correct-r2-barcode-file", "/scratch/lab_bock/shared/projects/sci-rna/data/PD195-1_Tcells_765k-sample1_P7/PD195-1_Tcells_765k-sample1_P7.fixed_barcodes.mapping.tsv",
    #     "--sample-name", "PD195-1_Tcells_765k-sample1_P7_A01",
    #     "/home/arendeiro/sci-rna/data/PD195-1_Tcells_765k-sample1_P7/PD195-1_Tcells_765k-sample1_P7_A01/PD195-1_Tcells_765k-sample1_P7_A01.*.STAR.Aligned.out.bam.featureCounts.bam",
    #     "/home/arendeiro/sci-rna/data/PD195-1_Tcells_765k-sample1_P7/PD195-1_Tcells_765k-sample1_P7_A01/PD195-1_Tcells_765k-sample1_P7_A01"]

    # # For whole experiment
    # cli = [
    #     "--sample-name", "PD193_humanmouse_765k",
    #     "--r1-annot", "metadata/sciRNA-seq.PD193_humanmouse_765k.oligos_2019-09-20.csv",
    #     "--r1-attributes", "plate_well",
    #     "--cell-barcodes", "r2",
    #     "--only-summary",
    #     "--no-save-intermediate",
    #     "--min-umi-output", "3",
    #     "--expected-cell-number", "250326",
    #     "--no-output-header",
    #     "--save-gene-expression",
    #     "--correct-r1-barcodes",
    #     "--correct-r2-barcodes",
    #     "--correct-r2-barcode-file", "/scratch/lab_bock/shared/projects/sci-rna/data/PD193_humanmouse_765k/PD193_humanmouse_765k.fixed_barcodes.mapping.tsv",
    #     "/scratch/lab_bock/shared/projects/sci-rna/data/PD193_humanmouse_765k/PD193_humanmouse_765k_*/PD193_humanmouse_765k_*.ALL.STAR.Aligned.out.bam.featureCounts.bam",
    #     "/scratch/lab_bock/shared/projects/sci-rna/data/PD193_humanmouse_765k/PD193_humanmouse_765k/PD193_humanmouse_765k"]

    args = parser.parse_args(cli)
    if args.sample_name is None:
        args.sample_name = args.output_prefix
    if args.r1_attributes is None:
        args.r1_attributes = []
    else:
        args.r1_attributes = args.r1_attributes.split(",")
    if not isinstance(args.r1_attributes, list):
        raise ValueError(
            "Incorrect r1_attributes. Got: '{}'".format(args.r1_attributes)
        )
    args.save_intermediate = not args.no_save_intermediate
    args.output_header = not args.no_output_header
    if not args.output_prefix.endswith("."):
        args.output_prefix += "."

    # expand glob if needed
    if len(args.input_files) == 1:
        if "*" in args.input_files[0]:
            args.input_files = glob(args.input_files[0])

    if not all([os.path.exists(f) for f in args.input_files]):
        raise ValueError("Not all input files exist!")

    return args


def main(cli=None):
    global args
    global attrs

    args = parse_args(cli)
    print(f"# {time.asctime()} - CLI arguments:")
    print(args)

    # barcode annotations
    annotation = pd.read_csv(args.r1_annotation_file)
    r2_barcodes = pd.read_csv(args.r2_barcodes)
    attrs = (
        annotation.query(f"sample_name == '{args.sample_name}'")
        .set_index("combinatorial_barcode")[args.r1_attributes]
        .squeeze(axis=0)
    )

    # read text files
    df = parse_data(args.input_files, nrows=args.nrows)
    # print(f"# {time.asctime()} - Saving all reads to pickle.")
    # to_pickle(df, "df", array=False)

    # correct r1 barcodes
    if args.correct_r1_barcodes:
        print(
            f"# {time.asctime()} - Updating r1 barcodes to the sequence of the well."
        )
        df["r1"] = attrs.name
        # df['r1'] = df['r1'].value_counts().idxmax()

    # correct r2 barcodes
    args.output_suffix = ""
    if args.correct_r2_barcodes:
        print(
            f"# {time.asctime()} - Updating barcodes not matching reference with corrected ones."
        )
        mapping = pd.read_csv(
            args.correct_r2_barcode_file, header=None, sep="\t"
        )
        mapping = mapping.dropna().set_index(0).squeeze().to_dict()

        # update
        unmatch = ~df["r2"].isin(r2_barcodes[args.barcode_orientation])
        df.loc[unmatch, "r2"] = [
            mapping[x] if x in mapping else x for x in df.loc[unmatch, "r2"]
        ]
        args.output_suffix = "_corrected"

    nr = df.isnull().sum().sum()
    if nr > 0:
        print(f"# {time.asctime()} - Found {nr} null reads. Dropping those.")
        df = df.dropna()
    if "int" not in str(df.dtypes["pos"]):
        df.loc[:, "pos"] = df.loc[:, "pos"].astype(int)

    # Save bulk profile
    print(f"# {time.asctime()} - Saving bulk expression profile.")
    (
        df["gene"]
        .value_counts()
        .sort_values(ascending=False)
        .to_frame(name=args.sample_name)
        .to_csv(
            os.path.join(args.output_prefix + "bulk_expression_profile.csv")
        )
    )

    # Gather metrics per cell
    r1_annotation = None
    if "r1" in args.cell_barcodes:
        r1_annotation = annotation.set_index("combinatorial_barcode")[
            args.r1_attributes
        ]
        r1_annotation.index.name = "r1"
    metrics = gather_stats_per_cell(
        df,
        r1_annotation=r1_annotation,
        species_mixture=args.species_mixture,
        save_intermediate=args.save_intermediate,
        suffix=args.output_suffix,
    )

    if "r1" not in args.cell_barcodes:
        # Add required attributes
        metrics = metrics.assign(**dict(zip(attrs.index, attrs.values)))

    # Save
    int_cols = ["read", "unique_umis", "umi", "gene"]
    metrics.loc[:, int_cols] = metrics.loc[:, int_cols].astype(int)
    if not args.output_header:
        print(
            f"# {time.asctime()} - Not outputing header in file, but header is the following:"
        )
        print(f"# {time.asctime()} - {metrics.columns}")
    (
        metrics.query(f"umi > {args.min_umi_output}").to_csv(
            args.output_prefix + "metrics" + args.output_suffix + ".csv.gz",
            float_format="%.3f",
            header=args.output_header,
        )
    )

    if "r1" in args.cell_barcodes:
        # # now the same only for exactly matching barcodes
        metrics_filtered = get_exact_matches(
            metrics,
            r1_whitelist=annotation["combinatorial_barcode"],
            r2_whitelist=r2_barcodes[args.barcode_orientation],
            save_intermediate=args.save_intermediate,
            plot=not args.only_summary,
        )
        metrics_filtered = metrics_filtered.assign(
            **dict(zip(attrs.index, attrs.values))
        )
        metrics_filtered.loc[:, int_cols] = metrics_filtered.loc[
            :, int_cols
        ].astype(int)
        if not args.output_header:
            print(
                f"# {time.asctime()} - Not outputing header in file, but header is the following:"
            )
            print(f"# {time.asctime()} - {metrics_filtered.columns}")
        (
            metrics_filtered.query(f"umi > {args.min_umi_output}").to_csv(
                args.output_prefix
                + "metrics"
                + args.output_suffix
                + "_filtered.csv.gz",
                float_format="%.3f",
                header=args.output_header,
            )
        )

    if args.only_summary:
        return

    # # Plot
    # # # Loglog line plot of rank vs abundance
    plot_metrics_lineplot(metrics)
    if "r1" in args.cell_barcodes:
        plot_metrics_lineplot(metrics_filtered, suffix="exact_match")

    # # # Efficiency plot (UMI and gene yield in relation to reads per cell)
    t = int(args.expected_cell_number * 5)
    plot_efficiency(metrics, tail=t)
    plot_efficiency(metrics, tail=t, colour_by="unique_fraction")
    if "r1" in args.cell_barcodes:
        plot_efficiency(metrics_filtered, tail=t, suffix="exact_match")
        plot_efficiency(
            metrics_filtered,
            tail=t,
            suffix="exact_match",
            colour_by="unique_fraction",
        )

    if args.species_mixture:
        plot_species_mixing(metrics)

    # Well and droplet inspection
    if len(args.cell_barcodes) > 1:
        well_metrics = gather_stats_per_well(metrics, seq_content=True)
        plot_well_stats(well_metrics, tail=None, suffix="")
        if "r1" in args.cell_barcodes:
            well_metrics_filtered = gather_stats_per_well(
                metrics_filtered,
                seq_content=True,
                save_intermediate=args.save_intermediate,
            )
            plot_well_stats(
                well_metrics_filtered, tail=None, suffix="exact_match"
            )


def parse_data(files, nrows=1e10):

    print(f"# {time.asctime()} - Parsing files.")
    pieces = list()
    for file in files:
        print(f"# {time.asctime()} - Parsing file {file}.")

        bam = pysam.AlignmentFile(file)
        for i, read in enumerate(bam):
            if i >= nrows:
                break
            # Filter
            if (
                read.is_qcfail
                or read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
            ):
                continue
            if "Unassigned" in read.get_tag("XS"):
                continue

            gene = read.get_tag("XT", with_value_type=True)
            if gene[1] == "i":
                try:
                    gene = [
                        x[1]
                        for x in read.get_tags(with_value_type=True)
                        if x[0] == "XT" and x[2] == "Z"
                    ]
                    # print(f"# {time.asctime()} - Found read with double XT tag"
                    #       f" (gene assignemnt): {read.qname}. Fixed it.")
                except IndexError:
                    print(
                        f"# {time.asctime()} - Found read with non character"
                        f" gene assignemnt: {read.qname}. Skipping it."
                    )
                    continue

            piece = [
                # save only the read flowcell position
                read.qname.split("#")[0],
                # ":".join(read.qname.split("#")[0].split(":")[1:]),
                read.get_tag("BC")[0:13],
                read.get_tag("r2"),
                read.get_tag("RX"),
                gene[0],
                read.pos,
            ]
            pieces.append(piece)
        print(f"# {time.asctime()} - Done with file {file}. {i} lines.")

    print(f"# {time.asctime()} - Concatenating parts.")
    return pd.DataFrame(
        pieces, columns=["read", "r1", "r2", "umi", "gene", "pos"]
    )


def gather_stats_per_cell(
    df,
    doublet_threshold=0.85,
    save_intermediate=True,
    species_mixture=True,
    norm_species=True,
    r1_annotation=None,
    suffix="",
):
    print(f"# {time.asctime()} - Gathering metrics per cell.")

    # performance metrics
    # # number of unique reads per cell
    reads_per_cell = (
        df.groupby(args.cell_barcodes, sort=False)["read"]
        .nunique()
        .sort_values()
    )
    if save_intermediate:
        to_pickle(reads_per_cell, "reads_per_cell" + suffix)

    # # mapping rate per cell
    # TODO: add mapping rate per cell (needs additional file)

    # # duplication per cell
    reads_per_umi = df.groupby(
        args.cell_barcodes + ["gene", "pos"], sort=False
    )["umi"].size()
    reads_per_umi = reads_per_umi.reset_index(level=["gene", "pos"], drop=True)
    if save_intermediate:
        to_pickle(reads_per_umi, "reads_per_umi" + suffix)

    unique = reads_per_umi == 1
    unique_per_cell = (
        unique.groupby(level=args.cell_barcodes).sum().rename("unique_umis")
    )
    if save_intermediate:
        to_pickle(unique_per_cell, "unique_per_cell" + suffix)

    # # UMI count per cell
    umi_counts = df.groupby(args.cell_barcodes + ["gene", "pos"], sort=False)[
        "umi"
    ].nunique()
    if args.save_gene_expression:
        print(f"# {time.asctime()} - Writing gene expression.")
        (
            umi_counts.reset_index(level="pos", drop=True)
            .reset_index()
            # Add required attributes
            .assign(**dict(zip(attrs.index, attrs.values)))
            .to_csv(
                args.output_prefix + "expression" + suffix + ".csv.gz",
                index=False,
                header=args.output_header,
            )
        )

    if save_intermediate:
        to_pickle(umi_counts, "umi_counts" + suffix)
    umis_per_cell = (
        umi_counts.groupby(args.cell_barcodes, sort=False).sum().sort_values()
    )
    if save_intermediate:
        to_pickle(umis_per_cell, "umis_per_cell" + suffix)

    # # Genes per cell
    genes_per_cell = (
        umi_counts.reset_index(level="gene")
        .groupby(args.cell_barcodes, sort=False)["gene"]
        .nunique()
        .sort_values()
    )
    if save_intermediate:
        to_pickle(genes_per_cell, "genes_per_cell" + suffix)

    # Species mixing
    if species_mixture:
        print(
            f"# {time.asctime()} - Gathering species-specific metrics per cell."
        )
        # # per UMI
        umi_counts2 = umi_counts.reset_index()
        umi_counts2 = umi_counts2.assign(
            species=(
                umi_counts2["gene"]
                .str.startswith("ENSG")
                .replace(True, "human")
                .replace(False, "mouse")
            )
        )
        species_counts = umi_counts2.groupby(args.cell_barcodes + ["species"])[
            "umi"
        ].sum()

        spc = species_counts.reset_index().pivot_table(
            index=args.cell_barcodes,
            columns="species",
            values="umi",
            fill_value=0,
        )
        spc += 1
        spc = spc.assign(total=spc.sum(1), max=spc.max(1))
        spc = spc.assign(
            ratio=spc["max"] / spc["total"],
            sp_ratio=spc["human"] / spc["total"],
        ).sort_values("total")
        spc = spc.assign(
            doublet=(spc["ratio"] < doublet_threshold)
            .astype(int)
            .replace(0, -1)
        )
        if save_intermediate:
            to_pickle(spc, "spc" + suffix)

        # # per read
        read_counts = (
            df.groupby(args.cell_barcodes + ["gene"], sort=False)["read"]
            .count()
            .reset_index()
        )
        read_counts = read_counts.assign(
            species=(
                read_counts["gene"]
                .str.startswith("ENSG")
                .replace(True, "read_human")
                .replace(False, "read_mouse")
            )
        )
        species_counts_read = read_counts.groupby(
            args.cell_barcodes + ["species"]
        )["read"].sum()

        spc_read = species_counts_read.reset_index().pivot_table(
            index=args.cell_barcodes,
            columns="species",
            values="read",
            fill_value=0,
        )
        spc_read += 1
        spc_read = spc_read.assign(
            read_total=spc_read.sum(1), read_max=spc_read.max(1)
        )
        spc_read = spc_read.assign(
            read_ratio=spc_read["read_max"] / spc_read["read_total"],
            read_sp_ratio=spc_read["read_human"] / spc_read["read_total"],
        ).sort_values("read_total")
        spc_read = spc_read.assign(
            read_doublet=(spc_read["read_ratio"] < doublet_threshold)
            .astype(int)
            .replace(0, -1)
        )
        if save_intermediate:
            to_pickle(spc_read, "spc_read" + suffix)

    print(f"# {time.asctime()} - Joining all metrics.")
    # TODO: speed up by assigning to column using equallly sorted indexes
    metrics = (
        reads_per_cell.to_frame()
        .join(unique_per_cell)
        .join(umis_per_cell)
        .join(genes_per_cell)
    )

    if species_mixture:
        metrics = metrics.join(spc).join(spc_read)

    metrics = metrics.sort_values("umi")

    # Calculate unique fraction
    metrics.loc[:, "unique_fraction"] = metrics["unique_umis"] / metrics["umi"]

    if r1_annotation is not None:
        print(f"# {time.asctime()} - Adding well annotation to metrics.")
        r1_annotation.index.name = "r1"
        metrics = metrics.join(r1_annotation)
    if save_intermediate:
        to_pickle(metrics, "metrics" + suffix)

    if (not species_mixture) or (not norm_species):
        return metrics

    # Assess species bias
    # # per UMI
    r = dict()
    for f in [0.1, 0.2, 0.25, 0.5, 0.75] + list(range(1, 1000, 10)) + [1.5]:
        t = metrics.tail(int(args.expected_cell_number * f))
        r[args.expected_cell_number * f] = t["mouse"].mean() / t["human"].mean()
    r = pd.Series(r).sort_index()

    # Add normalized metrics to stats
    metrics.loc[:, "human_norm"] = (
        metrics["human"] * r[int(args.expected_cell_number)]
    )
    metrics.loc[:, "total_norm"] = metrics[["mouse", "human_norm"]].sum(1)
    metrics.loc[:, "max_norm"] = metrics[["mouse", "human_norm"]].max(1)
    metrics.loc[:, "ratio_norm"] = metrics["max_norm"] / metrics["total_norm"]
    metrics.loc[:, "sp_ratio_norm"] = (
        metrics["human_norm"] / metrics["total_norm"]
    )
    metrics.loc[:, "doublet_norm"] = (
        (metrics.loc[:, "ratio_norm"] < doublet_threshold)
        .astype(int)
        .replace(0, -1)
    )

    # Assess species bias
    # # per read
    r = dict()
    for f in [0.1, 0.2, 0.25, 0.5, 0.75] + list(range(1, 1000, 10)) + [1.5]:
        t = metrics.tail(int(args.expected_cell_number * f))
        r[args.expected_cell_number * f] = (
            t["read_mouse"].mean() / t["read_human"].mean()
        )
    r = pd.Series(r).sort_index()

    # Add normalized metrics to stats
    metrics.loc[:, "read_human_norm"] = (
        metrics["read_human"] * r[int(args.expected_cell_number)]
    )
    metrics.loc[:, "read_total_norm"] = metrics[
        ["read_mouse", "read_human_norm"]
    ].sum(1)
    metrics.loc[:, "read_max_norm"] = metrics[
        ["read_mouse", "read_human_norm"]
    ].max(1)
    metrics.loc[:, "read_ratio_norm"] = (
        metrics["read_max_norm"] / metrics["read_total_norm"]
    )
    metrics.loc[:, "read_sp_ratio_norm"] = (
        metrics["read_human_norm"] / metrics["read_total_norm"]
    )
    metrics.loc[:, "read_doublet_norm"] = (
        (metrics.loc[:, "read_ratio_norm"] < doublet_threshold)
        .astype(int)
        .replace(0, -1)
    )

    if save_intermediate:
        to_pickle(metrics, "metrics" + suffix)

    # # # plot ammount of species-specific bias
    # try:
    #     n = 1
    #     fig, axis = plt.subplots(n, 1, figsize=(3, n * 3), tight_layout=True)
    #     axis.plot(r.index, r, "o-")
    #     axis.axvline(args.expected_cell_number, linestyle="--", color="grey")
    #     axis.set_xlabel("Top N barcodes taken into account")
    #     axis.set_ylabel("Ratio of mean UMIs per cell\nbetween mouse and human")
    #     axis.set_xscale("log")
    #     fig.savefig(
    #         args.output_prefix + f"species_bias.lineplot.svg".replace("..", "."),
    #         dpi=300,
    #         bbox_inches="tight",
    #     )
    # except:
    #     print("Couldn't plot species bias plot")

    # # Plot illustration of normalization procedure
    # metrics2 = metrics.tail(int(2 * args.expected_cell_number))

    # for label in ["", ".log"]:
    #     fig, axis = plt.subplots(2, 2, figsize=(2 * 3, 2 * 3), tight_layout=True)
    #     kwargs = {"s": 0.1, "alpha": 0.2, "rasterized": True, "cmap": get_custom_cmap()}
    #     axis[0, 0].scatter(
    #         metrics2["mouse"], metrics2["human"], c=metrics2["sp_ratio"], **kwargs
    #     )
    #     axis[0, 1].scatter(
    #         metrics2["mouse"], metrics2["human"], c=metrics2["sp_ratio_norm"], **kwargs
    #     )
    #     v = metrics2[["mouse", "human"]].quantile(0.999).max()
    #     v += v * 0.1
    #     for ax in axis[0, :]:
    #         if label == "":
    #             ax.set_xlim((-(v / 30.0), v))
    #             ax.set_ylim((-(v / 30.0), v))
    #         else:
    #             ax.loglog()
    #         ax.plot((0, v), (0, v), linestyle="--", color="grey")
    #         ax.plot(
    #             (0, v),
    #             (0, v * (1 / r[int(args.expected_cell_number)])),
    #             linestyle="--",
    #             color="orange",
    #         )
    #         ax.set_ylabel("Human (UMIs)")

    #     axis[1, 0].scatter(
    #         metrics2["mouse"], metrics2["human_norm"], c=metrics2["sp_ratio"], **kwargs
    #     )
    #     axis[1, 1].scatter(
    #         metrics2["mouse"],
    #         metrics2["human_norm"],
    #         c=metrics2["sp_ratio_norm"],
    #         **kwargs,
    #     )
    #     v = metrics2[["mouse", "human_norm"]].quantile(0.999).max()
    #     v += v * 0.1
    #     for ax in axis[1, :]:
    #         if label == "":
    #             ax.set_xlim((-(v / 30.0), v))
    #             ax.set_ylim((-(v / 30.0), v))
    #         else:
    #             ax.loglog()
    #         ax.plot((0, v), (0, v), linestyle="--", color="grey")
    #         ax.plot((0, v), (0, v), linestyle="--", color="orange")
    #         ax.set_ylabel("Human (norm UMIs)")

    #     for ax in axis.flatten():
    #         ax.set_xlabel("Mouse (UMIs)")
    #     axis[0, 0].set_title("Coloured by ratio")
    #     axis[0, 1].set_title("Coloured by norm ratio")
    #     fig.savefig(
    #         args.output_prefix
    #         + f"species_bias.original_vs_corrected.scatterplot{label}.svg".replace(
    #             "..", "."
    #         ),
    #         dpi=300,
    #         bbox_inches="tight",
    #     )

    return metrics


def write_gene_expression_matrix(output_file=None, top_cells=None):

    print(f"# {time.asctime()} - Producing sparse expression matrix.")

    if top_cells is None:
        top_cells = args.expected_cell_number
    if output_file is None:
        output_file = args.output_prefix + "h5ad"

    umi_counts = from_pickle("umi_counts")
    umis_per_cell = from_pickle("umis_per_cell")
    c = umi_counts.reset_index(level="pos", drop=True).reset_index(level="gene")

    print(f"# {time.asctime()} - Factorizing cells and genes.")
    accept = umis_per_cell.tail(top_cells).index
    c2 = c.loc[c.index.isin(accept)].reset_index()
    c2 = c2.assign(cell=c2[["r1", "r2"]].sum(1))
    cell, unique_cell = c2["cell"].factorize()
    gene, unique_gene = c2["gene"].factorize()

    print(f"# {time.asctime()} - Creating AnnData object.")
    a = an.AnnData(
        csr_matrix(
            (c2["umi"], (cell, gene)),
            shape=(len(unique_cell), len(unique_gene)),
        ),
        dtype=np.int,
    )
    a.obs.index = unique_cell
    a.var.index = unique_gene

    print(f"# {time.asctime()} - Writing h5ad object to disk.")
    a.write(output_file)


def write_full_gene_expression_matrix(output_file=None, top_cells=None):

    print(f"# {time.asctime()} - Producing sparse expression matrix.")

    if output_file is None:
        output_file = args.output_prefix + "h5ad"

    c2 = pd.read_csv("data/SCI024_Tcell_s/SCI024_Tcell_s.expression.csv.gz")
    cell, unique_cell = c2["r2"].factorize()
    gene, unique_gene = c2["gene"].factorize()

    print(f"# {time.asctime()} - Creating AnnData object.")
    a = an.AnnData(
        csr_matrix(
            (c2["umi"], (cell, gene)),
            shape=(len(unique_cell), len(unique_gene)),
        ),
        dtype=np.int,
    )
    a.obs.index = unique_cell
    c3 = c2.set_index("r2")[
        ["sample_name", "donor_id", "well"]
    ].drop_duplicates()
    a.obs = a.obs.join(c3)
    a.var.index = unique_gene

    print(f"# {time.asctime()} - Writing h5ad object to disk.")
    a.write("data/SCI024_Tcell_s/SCI024_Tcell_s.expression.h5ad")


def gather_stats_per_cell_as_droplet(
    df,
    cell_barcodes=["r2"],
    doublet_threshold=0.85,
    save_intermediate=True,
    norm_species=True,
    r1_annotation=None,
    suffix="_r2_only",
):
    print(f"# {time.asctime()} - Gathering metrics per cell.")

    # performance metrics
    # # UMI count per cell
    umi_counts = df.groupby(cell_barcodes + ["gene", "pos"], sort=False)[
        "umi"
    ].nunique()
    if save_intermediate:
        to_pickle(umi_counts, "umi_counts" + suffix)

    umis_per_cell = (
        umi_counts.groupby(cell_barcodes, sort=False).sum().sort_values()
    )
    if save_intermediate:
        to_pickle(umis_per_cell, "umis_per_cell" + suffix)

    # species mixing
    print(f"# {time.asctime()} - Gathering species-specific metrics per cell.")
    umi_counts2 = umi_counts.reset_index()
    umi_counts2 = umi_counts2.assign(
        species=(
            umi_counts2["gene"]
            .str.startswith("ENSG")
            .replace(True, "human")
            .replace(False, "mouse")
        )
    )
    species_counts = umi_counts2.groupby(cell_barcodes + ["species"])[
        "umi"
    ].sum()

    spc = species_counts.reset_index().pivot_table(
        index=cell_barcodes, columns="species", values="umi", fill_value=0
    )
    spc += 1
    spc = spc.assign(total=spc.sum(1), max=spc.max(1))
    spc = spc.assign(
        ratio=spc["max"] / spc["total"], sp_ratio=spc["human"] / spc["total"]
    ).sort_values("total")
    spc = spc.assign(
        doublet=(spc["ratio"] < doublet_threshold).astype(int).replace(0, -1)
    )
    if save_intermediate:
        to_pickle(spc, "spc" + suffix)

    print(f"# {time.asctime()} - Joining all metrics.")
    # TODO: speed up by assigning to column using equallly sorted indexes
    metrics_r2 = umis_per_cell.to_frame().join(spc).sort_values("umi")

    if save_intermediate:
        to_pickle(metrics_r2, "metrics_r2" + suffix)

    if not norm_species:
        return metrics_r2

    # Assess species bias
    r = dict()
    for f in [0.1, 0.2, 0.25, 0.5, 0.75] + list(range(1, 1000, 10)) + [1.5]:
        t = metrics_r2.tail(int(args.expected_cell_number * f))
        r[args.expected_cell_number * f] = t["mouse"].mean() / t["human"].mean()
    r = pd.Series(r).sort_index()

    # Add normalized metrics to stats
    metrics_r2.loc[:, "human_norm"] = (
        metrics_r2["human"] * r[int(args.expected_cell_number)]
    )
    metrics_r2.loc[:, "total_norm"] = metrics_r2[["mouse", "human_norm"]].sum(1)
    metrics_r2.loc[:, "max_norm"] = metrics_r2[["mouse", "human_norm"]].max(1)
    metrics_r2.loc[:, "ratio_norm"] = (
        metrics_r2["max_norm"] / metrics_r2["total_norm"]
    )
    metrics_r2.loc[:, "sp_ratio_norm"] = (
        metrics_r2["human_norm"] / metrics_r2["total_norm"]
    )
    metrics_r2.loc[:, "doublet_norm"] = (
        (metrics_r2.loc[:, "ratio_norm"] < doublet_threshold)
        .astype(int)
        .replace(0, -1)
    )

    if save_intermediate:
        to_pickle(metrics_r2, "metrics_r2" + suffix)

    # Plot ammount of species-specific bias
    n = 1
    fig, axis = plt.subplots(n, 1, figsize=(3, n * 3), tight_layout=True)
    axis.plot(r.index, r, "o-")
    axis.axvline(args.expected_cell_number, linestyle="--", color="grey")
    axis.set_xlabel("Top N barcodes taken into account")
    axis.set_ylabel("Ratio of mean UMIs per cell\nbetween mouse and human")
    axis.set_xscale("log")
    fig.savefig(
        args.output_prefix
        + f"species_bias.lineplot.{suffix}.svg".replace("..", "."),
        dpi=300,
        bbox_inches="tight",
    )

    # Plot illustration of normalization procedure
    metrics2 = metrics_r2.tail(int(10 * args.expected_cell_number))

    for label in ["", ".log"]:
        fig, axis = plt.subplots(
            1, 2, figsize=(2 * 3, 1 * 3), tight_layout=True, squeeze=False
        )
        kwargs = {
            "s": 0.1,
            "alpha": 0.2,
            "rasterized": True,
            "cmap": get_custom_cmap(),
        }
        axis[0, 0].scatter(
            metrics2["mouse"],
            metrics2["human"],
            c=metrics2["sp_ratio"],
            **kwargs,
        )
        # axis[0, 1].scatter(metrics2['mouse'], metrics2['human'], c=metrics2['sp_ratio_norm'], **kwargs)
        v = metrics2[["mouse", "human"]].quantile(0.999).max()
        v += v * 0.1
        for ax in axis[0, :]:
            if label == "":
                ax.set_xlim((-(v / 30.0), v))
                ax.set_ylim((-(v / 30.0), v))
            else:
                ax.loglog()
            ax.plot((0, v), (0, v), linestyle="--", color="grey")
            ax.plot(
                (0, v),
                (0, v * (1 / r[int(args.expected_cell_number)])),
                linestyle="--",
                color="orange",
            )
            ax.set_ylabel("Human (UMIs)")

        for ax in axis.flatten():
            ax.set_xlabel("Mouse (UMIs)")
        axis[0, 0].set_title("Coloured by ratio")
        axis[0, 1].set_title("Coloured by norm ratio")
        fig.savefig(
            args.output_prefix
            + f"species_bias.original_vs_corrected.scatterplot{label}.{suffix}.svg".replace(
                "..", "."
            ),
            dpi=300,
            bbox_inches="tight",
        )

    return metrics_r2


def plot_metrics_lineplot(
    metrics,
    keys=["read", "umi", "gene"],
    tail=None,
    suffix="",
    by_group=None,
    always_legend=False,
):
    def min_max(x):
        return (x - x.min()) / (x.max() - x.min())

    print(f"# {time.asctime()} - Plotting metrics per cell.")
    if tail:
        metrics = metrics.tail(tail)
    n = len(keys)

    if by_group is not None:
        groups = sorted(metrics[by_group].dropna().unique())
        suffix += f"by_{by_group}"
    else:
        by_group = "dummy"
        groups = ["dummy"]
        metrics.loc[:, by_group] = groups[0]

    colors = sns.color_palette("husl", n_colors=len(groups))
    styles = ["solid", "dashed", "dashdot", "dotted"]
    num_styles = len(styles)

    fig, axis = plt.subplots(
        n, 2, figsize=(2 * 6, n * 3)
    )  # , tight_layout=True)
    for i, metric in enumerate(keys):
        for j, group in enumerate(groups):
            d = metrics.loc[metrics[by_group] == group, metric].sort_values()
            d_norm = min_max(d)
            rank = d.rank(ascending=False, method="average", pct=False)
            rank_norm = min_max(rank)
            kwargs = {
                "color": colors[j],
                "linestyle": styles[j % num_styles],
                "rasterized": True,
                "label": group if group != "dummy" else None,
            }
            axis[i, 0].plot(rank, d, **kwargs)
            axis[i, 1].plot(rank_norm, d_norm, **kwargs)
        axis[i, 0].axvline(
            args.expected_cell_number, linestyle="--", color="grey"
        )
        for l, ax in enumerate(axis[i, :]):
            ax.loglog()
            ax.set_xlabel("Barcodes")
            ax.set_ylabel(metric.capitalize() + "s")
            # add legend only to one panel if requested
            if (by_group != "dummy") and (
                always_legend or ((i == 0) and (l == 1))
            ):
                ax.legend(
                    bbox_to_anchor=(1.05, 1),
                    loc="upper left",
                    borderaxespad=0.0,
                    ncol=1,
                )
                # ax.legend_.set_in_layout(False)
    fig.savefig(
        args.output_prefix
        + f"metrics_per_cell.lineplot.{suffix}.svg".replace("..", "."),
        dpi=300,
        bbox_inches="tight",
    )


def plot_metrics_distplot(
    metrics, keys=["read", "umi", "gene"], tail=None, suffix="", by_group=None
):
    print(f"# {time.asctime()} - Plotting metrics per cell.")
    if tail:
        metrics = metrics.tail(tail)
    n = len(keys)

    if by_group is not None:
        groups = metrics[by_group].dropna().unique()
        suffix += f"by_{by_group}"
    else:
        by_group = "dummy"
        groups = ["dummy"]
        metrics.loc[:, by_group] = groups[0]

    colors = sns.color_palette("husl", n_colors=len(groups))

    fig, axis = plt.subplots(
        1,
        n,
        figsize=(n * 3, 1 * 3),
        tight_layout=True,
        sharex="col",
        squeeze=True,
    )
    for i, metric in enumerate(keys):
        axis[i].set_yscale("log")
        for j, group in enumerate(groups):
            d = metrics.loc[metrics[by_group] == group, :]
            sns.distplot(
                np.log10(d[metric]),
                kde=False,
                ax=axis[i],
                label=group if group != "dummy" else None,
            )
        axis[i].set_xlabel(metric.capitalize() + "s (log)")
        axis[i].set_ylabel("Barcodes")
        if by_group != "dummy":
            axis[i].legend()
    fig.savefig(
        args.output_prefix
        + f"metrics_per_cell.distplot.{suffix}.svg".replace("..", "."),
        dpi=300,
        bbox_inches="tight",
    )


def plot_efficiency(
    metrics,
    keys=["umi", "gene"],
    tail=25000,
    suffix="",
    by_group=None,
    colour_by=None,
):
    """
    `metrics` must contain a column "read".
    """
    print(f"# {time.asctime()} - Plotting efficiency metrics per cell.")
    if tail:
        metrics = metrics.tail(tail)

    if colour_by is None:
        colour_by = "efficiency"
        if colour_by not in metrics.columns:
            metrics.loc[:, colour_by] = metrics["umi"] / metrics["read"]

    if by_group is not None:
        groups = metrics[by_group].dropna().unique()
        n = len(groups)
        suffix += f"by_{by_group}"
    else:
        by_group = "dummy"
        groups = ["dummy"]
        metrics.loc[:, by_group] = groups[0]
        n = 1

    fig, axis = plt.subplots(
        n,
        2,
        figsize=(2 * 4, n * 4),
        tight_layout=True,
        squeeze=False,
        sharex="col",
        sharey="col",
    )
    fig.suptitle(f"Experiment performance: {args.sample_name}", ha="center")
    for i, group in enumerate(groups):
        d = metrics.loc[metrics[by_group] == group, :]
        if group != "dummy":
            axis[i, 0].set_title(f"{by_group.capitalize()}: {group}")
        for j, metric in enumerate(keys):
            # # fit curve
            try:
                (m, b), pcov = scipy.optimize.curve_fit(
                    lin_func, d["read"], d[metric]
                )
                # lm = LinearRegression().fit(d['read'].values.reshape(-1, 1), d[metric])
                # assert np.allclose(m, lm.coef_)
                # assert np.allclose(b, lm.intercept_)
                axis[i, j].text(
                    0.5e5,
                    args.expected_cell_number,
                    s="y = {:.3f}x + {:.1f}".format(m, b),
                    ha="center",
                )
                # Plot regression line
                x = np.linspace(d["read"].min(), d["read"].max(), num=1000)
                axis[i, j].plot(x, lin_func(x, m, b), color="orange")
            except TypeError:
                print(f"Couldn't find a fit for {group}, {metric}.")
                pass

            # Plot cells
            col = axis[i, j].scatter(
                d["read"],
                d[metric],
                c=d[colour_by],
                rasterized=True,
                alpha=0.2,
                s=2,
            )
            add_colorbar_to_axis(col, label=colour_by)
            axis[i, j].loglog()
            axis[i, j].set_xlabel("Reads per cell")
            axis[i, j].set_ylabel(f"Useful {metric}s per cell")

            # Plot 10k annotation
            y = lin_func(1e4, m, b)
            axis[i, j].text(
                1e4,
                y,
                s=f"{metric.capitalize()}s recovered\nwith 10.000\nreads:\n{y:.2f}",
                ha="left",
            )
            # Plot X == Y
            xmax = d["read"].max()
            xmax += xmax * 0.1
            ymax = d[metric].max()
            ymax += ymax * 0.1
            x = np.linspace(0, ymax, num=2)
            y = lin_func(x, 1, 0)
            axis[i, j].plot(x, y, linestyle="--", color="black", linewidth=0.5)

            # Plot lines at relevant metrics
            for h in [100, 250, 500, 1000, args.expected_cell_number]:
                axis[i, j].axhline(
                    h, linestyle="--", color="black", linewidth=0.5
                )
            for v in [10000, 100000]:
                axis[i, j].axvline(
                    v, linestyle="--", color="black", linewidth=0.5
                )
    fig.savefig(
        args.output_prefix
        + f"performance_per_cell.scatter.coloured_by_{colour_by}.{suffix}.svg".replace(
            "..", "."
        ),
        dpi=300,
        bbox_inches="tight",
    )


def plot_species_mixing(
    metrics,
    norm=False,
    tail=None,
    suffix="",
    cmaps=None,
    zoom_in_area=3000,
    axislims=10000,
):
    print(f"# {time.asctime()} - Plotting species mixtures.")
    if tail is None:
        tail = 2 * args.expected_cell_number
    if tail is not False:
        metrics = metrics.tail(tail)

    if norm:
        suffix += "_norm"

    if cmaps is None:
        cmaps = [
            get_custom_cmap()
        ]  # , plt.get_cmap("coolwarm"), plt.get_cmap("Spectral_r")]
    for attr, label, kwargs in [
        (
            "sp_ratio_norm" if norm else "sp_ratio",
            "coloured_by_ratio",
            {"vmin": 0, "vmax": 1},
        ),
        # ('doublet_norm' if norm else 'doublet', 'coloured_by_doublet', {"vmin": -1, "vmax": 1}),
    ]:
        for cmap in cmaps:
            n_panels = 4
            fig, axis = plt.subplots(
                1, n_panels, figsize=(n_panels * 4, 4), tight_layout=True
            )
            for ax in axis[:3]:
                col = ax.scatter(
                    metrics["mouse"],
                    metrics["human_norm" if norm else "human"],
                    c=metrics[attr],
                    cmap=cmap,
                    s=2,
                    alpha=0.25,
                    rasterized=True,
                    **kwargs,
                )
                ax.set_xlabel("Mouse (UMIs)")
                ax.set_ylabel("Human (UMIs)")
                add_colorbar_to_axis(col, label="Species fraction")
            v = metrics["total_norm" if norm else "total"].max()
            v = min(axislims, v)
            axis[1].set_xlim((-(v / 30.0), v))
            axis[1].set_ylim((-(v / 30.0), v))
            axis[2].set_xlim((-(zoom_in_area / 30.0), zoom_in_area))
            axis[2].set_ylim((-(zoom_in_area / 30.0), zoom_in_area))
            axis[3].scatter(
                metrics["mouse"],
                metrics["human" if norm else "human"],
                c=metrics[attr],
                cmap=cmap,
                s=1,
                alpha=0.1,
                rasterized=True,
                **kwargs,
            )
            axis[3].loglog()
            axis[3].set_xlabel("Mouse (UMIs)")
            axis[3].set_ylabel("Human (UMIs)")
            fig.savefig(
                args.output_prefix
                + f"species_mix.{suffix}.{label}.{cmap.name}.svg".replace(
                    "..", "."
                ),
                dpi=300,
                bbox_inches="tight",
            )


def gather_stats_per_well(metrics, seq_content=False, save_intermediate=True):
    print(
        f"# {time.asctime()} - Gathering metrics per combinatorial indexing well."
    )
    # # # see how many of the round1 barcodes (well) are captured
    umis_per_well = metrics.groupby("r1")["umi"].sum().sort_values()
    umis_per_well.name = "umis"
    if save_intermediate:
        to_pickle(umis_per_well, "umis_per_well")

    # # # see how many droplets per well
    droplets_per_well = (
        metrics.reset_index().groupby("r1")["r2"].nunique().sort_values()
    )
    droplets_per_well.name = "droplets"
    if save_intermediate:
        to_pickle(droplets_per_well, "droplets_per_well")

    well_metrics = umis_per_well.to_frame().join(droplets_per_well)
    if save_intermediate:
        to_pickle(well_metrics, "well_metrics")

    if seq_content:
        s = well_metrics.index.to_series().str
        well_metrics = well_metrics.assign(
            at_content=s.count("AT|TA") / 11,
            gc_content=s.count("GC|CG") / 11,
            a_content=s.count("A|T") / 11,
            c_content=s.count("C|G") / 11,
        )
        if save_intermediate:
            to_pickle(well_metrics, "well_metrics")

    return well_metrics


def get_exact_matches(
    metrics, r1_whitelist, r2_whitelist, save_intermediate=True, plot=True
):
    print(f"# {time.asctime()} - Filtering for exact barcode matches.")

    # Note to self, if metrics has r1_annotations, one could simply do isnull() to get barcodes matching r1

    r1_match = metrics.index.get_level_values("r1").isin(r1_whitelist)
    r2_match = metrics.index.get_level_values("r2").isin(r2_whitelist)

    print(r1_match.sum() / r1_match.shape[0])
    print(r2_match.sum() / r2_match.shape[0])
    print(metrics.loc[r1_match, "read"].sum() / metrics["read"].sum())
    print(metrics.loc[r2_match, "read"].sum() / metrics["read"].sum())
    print(metrics.loc[r1_match, "umi"].sum() / metrics["umi"].sum())
    print(metrics.loc[r2_match, "umi"].sum() / metrics["umi"].sum())

    r1_match = metrics.index.get_level_values("r1").isin(r1_whitelist)
    r2_match = metrics.index.get_level_values("r2").isin(r2_whitelist)

    tr1_match = r1_match[-200000:]
    tr2_match = r2_match[-200000:]
    tmetrics = metrics.tail(200000)
    print(tr1_match.sum() / tr1_match.shape[0])
    print(tr2_match.sum() / tr2_match.shape[0])
    print(tmetrics.loc[tr1_match, "read"].sum() / tmetrics["read"].sum())
    print(tmetrics.loc[tr2_match, "read"].sum() / tmetrics["read"].sum())
    print(tmetrics.loc[tr1_match, "umi"].sum() / tmetrics["umi"].sum())
    print(tmetrics.loc[tr2_match, "umi"].sum() / tmetrics["umi"].sum())

    print(tmetrics["read"].sum() / metrics["read"].sum())
    print(tmetrics["umi"].sum() / metrics["umi"].sum())

    if save_intermediate:
        to_pickle(r1_match, "r1_match", array=False)
        to_pickle(r2_match, "r2_match", array=False)

    if plot:
        plot_barcode_match_fraction(r1_match, r2_match)
        plot_umi_match_fraction(metrics["umi"], r1_match, r2_match)
    metrics_filtered = metrics.loc[r1_match & r2_match]

    if save_intermediate:
        to_pickle(metrics_filtered, "metrics_filtered")
    return metrics_filtered

    # from Levenshtein import distance
    # from joblib import Parallel, delayed

    # query = metrics.index.get_level_values("r1")[~r1_match][:100]
    # reference = r1_whitelist
    # distances = Parallel(n_jobs=-1)(
    #     delayed(distance)(q, r)
    #     for q in query
    #     for r in reference)

    metrics = metrics.assign(r1_match=r1_match)
    m = metrics.tail(2000000).groupby("r1_match").mean()
    mm = m.reset_index().melt(id_vars=["r1_match"])
    grid = sns.catplot(
        data=mm,
        x="r1_match",
        col="variable",
        y="value",
        kind="bar",
        sharey=False,
        col_wrap=4,
        height=2,
        aspect=1,
    )
    grid.savefig(
        args.output_prefix + "metrics_dependent_on_r1_match.svg",
        bbox_inches="tight",
        dpi=300,
    )


def get_exact_matches_droplet(
    metrics, r2_whitelist, save_intermediate=True, plot=True, suffix="_r2_only"
):
    print(f"# {time.asctime()} - Filtering for exact barcode matches.")

    # Note to self, if metrics has r1_annotations, one could simply do isnull() to get barcodes matching r1
    r2_match = metrics.index.get_level_values("r2").isin(r2_whitelist)

    if save_intermediate:
        to_pickle(r2_match, "r2_match" + suffix, array=False)

    metrics_filtered = metrics.loc[r2_match]

    if save_intermediate:
        to_pickle(metrics_filtered, "metrics_filtered" + suffix)
    return metrics_filtered


def get_stats_per_droplet(
    metrics, doublet_threshold=0.85, save_intermediate=True
):
    metrics_droplet = metrics.groupby("r2")[
        "read", "umi", "gene", "human", "mouse", "total", "max"
    ].sum()

    metrics_droplet = metrics_droplet.assign(
        ratio=metrics_droplet["max"] / metrics_droplet["total"],
        sp_ratio=metrics_droplet["human"] / metrics_droplet["total"],
    ).sort_values("total")
    metrics_droplet = metrics_droplet.assign(
        doublet=(metrics_droplet["ratio"] < doublet_threshold)
        .astype(int)
        .replace(0, -1)
    )

    # Assess species bias
    r = dict()
    for f in [0.1, 0.2, 0.25, 0.5, 0.75] + list(range(1, 1000, 10)) + [1.5]:
        t = metrics_droplet.tail(int(args.expected_cell_number * f))
        r[args.expected_cell_number * f] = t["mouse"].mean() / t["human"].mean()
    r = pd.Series(r).sort_index()

    # Add normalized metrics_droplet to stats
    metrics_droplet.loc[:, "human_norm"] = (
        metrics_droplet["human"] * r[int(args.expected_cell_number)]
    )
    metrics_droplet.loc[:, "total_norm"] = metrics_droplet[
        ["mouse", "human_norm"]
    ].sum(1)
    metrics_droplet.loc[:, "max_norm"] = metrics_droplet[
        ["mouse", "human_norm"]
    ].max(1)
    metrics_droplet.loc[:, "ratio_norm"] = (
        metrics_droplet["max_norm"] / metrics_droplet["total_norm"]
    )
    metrics_droplet.loc[:, "sp_ratio_norm"] = (
        metrics_droplet["human_norm"] / metrics_droplet["total_norm"]
    )
    metrics_droplet.loc[:, "doublet_norm"] = (
        (metrics_droplet.loc[:, "ratio_norm"] < doublet_threshold)
        .astype(int)
        .replace(0, -1)
    )

    if save_intermediate:
        to_pickle(metrics_droplet, "metrics_droplet")

    return metrics_droplet


def plot_well_stats(
    well_metrics, keys=["droplets", "umis"], tail=None, suffix=""
):
    print(
        f"# {time.asctime()} - Plotting metrics per combinatorial indexing well."
    )
    if tail is None:
        tail = 384 * 2
    well_metrics = well_metrics.tail(tail)

    n = len(keys)
    v = well_metrics[keys].min().min()
    fig, axis = plt.subplots(n, 1, figsize=(2 * 3, n * 3), tight_layout=True)
    for i, metric in enumerate(keys):
        d = well_metrics[metric].sort_values()
        rank = d.rank(ascending=False, method="average")
        axis[i].plot(rank, d, rasterized=True)
        axis[i].set_ylim((v - (v * 0.5), axis[i].get_ylim()[1] * 1.5))
        axis[i].set_yscale("log")
        axis[i].set_xlabel("Well")
        axis[i].set_ylabel(metric.capitalize())
    fig.savefig(
        args.output_prefix
        + f"metrics_per_well.lineplot.{suffix}.svg".replace("..", "."),
        dpi=300,
        bbox_inches="tight",
    )

    # # # # observe whether barcode sequence content relates with abundance
    if well_metrics.columns.str.endswith("_content").sum() < 1:
        return

    cols = well_metrics.columns[well_metrics.columns.str.endswith("_content")]
    fig, axis = plt.subplots(
        n, len(cols), figsize=(len(cols) * 3, n * 3), tight_layout=True
    )
    for j, col in enumerate(cols):
        for i, key in enumerate(keys):
            axis[i, j].scatter(
                well_metrics[col],
                well_metrics[key],
                rasterized=True,
                s=2,
                alpha=0.5,
            )
            axis[i, j].set_xlim((-0.1, 1.1))
            axis[i, j].set_xlabel(col.replace("_", " ").upper())
            axis[i, j].set_ylabel(key.capitalize() + " per well")
    fig.savefig(
        args.output_prefix
        + f"metrics_per_well.sequence_content.{suffix}.svg".replace("..", "."),
        dpi=300,
        bbox_inches="tight",
    )


def plot_barcode_match_fraction(r1, r2, suffix=""):
    print(f"# {time.asctime()} - Plotting fraction of correct barcodes.")

    d = pd.DataFrame(
        [[r1.sum(), r2.sum(), (r1 & r2).sum()], [r1.shape[0]] * 3],
        index=["correct", "total"],
        columns=["r1", "r2", "both"],
    )
    d.loc["fraction of exact match", :] = d.loc["correct"] / d.loc["total"]

    grid = sns.catplot(
        data=d.reset_index().melt(id_vars="index"),
        x="variable",
        y="value",
        col="index",
        sharey=False,
        kind="bar",
        height=3,
        aspect=1,
    )
    grid.fig.savefig(
        args.output_prefix
        + f"barcode_matching.barplot.{suffix}.svg".replace("..", "."),
        dpi=300,
        bbox_inches="tight",
    )


def plot_umi_match_fraction(umi, r1_match, r2_match, suffix=""):
    print(f"# {time.asctime()} - Plotting fraction of correct barcodes.")

    d = pd.DataFrame(
        [
            [umi.loc[r1_match].sum(), umi.sum()],
            [umi.loc[r2_match].sum(), umi.sum()],
        ],
        index=["correct", "total"],
        columns=["r1", "r2"],
    )
    d.loc["fraction of exact match", :] = d.loc["correct"] / d.loc["total"]

    grid = sns.catplot(
        data=d.reset_index().melt(id_vars="index"),
        x="variable",
        y="value",
        col="index",
        sharey=False,
        kind="bar",
        height=3,
        aspect=1,
    )
    grid.fig.savefig(
        args.output_prefix
        + f"umi_matching.barplot.{suffix}.svg".replace("..", "."),
        dpi=300,
        bbox_inches="tight",
    )


def plot_comparison_to_10x(c, d, suffix=""):
    # # Compare side by side in logster
    for attr, label, kwargs in [
        (
            "sp_ratio",
            "coloured_by_ratio",
            {"cmap": get_custom_cmap(), "vmin": 0, "vmax": 1},
        ),
        # ('doublet', 'coloured_by_doublet', {"cmap": "tab20b", "vmin": -1, "vmax": 1}),
        # ('total', 'coloured_by_total', {"vmin": 100, "vmax": 10000}),
    ]:
        fig, axis = plt.subplots(
            2, 3, figsize=(3 * 4, 2 * 4), tight_layout=True
        )
        for ax in axis[:, 0]:
            ax.set_title("10X")
            ax.scatter(
                d["mouse"],
                d["human"],
                c=d[attr],
                s=1,
                alpha=0.1,
                rasterized=True,
                **kwargs,
            )
        for ax in axis[:, 1]:
            ax.set_title("scifi-RNA-seq")
            ax.scatter(
                c["mouse"],
                c["human"],
                c=c[attr],
                s=1,
                alpha=0.1,
                rasterized=True,
                **kwargs,
            )
        for ax in axis[:, 2]:
            ax.set_title("scifi-RNA-seq coloured as 10X only")
            c.loc[:, "10X_attr"] = d.loc[
                c.index.get_level_values(1), attr
            ].values
            ax.scatter(
                c["mouse"],
                c["human"],
                c=c["10X_attr"],
                s=1,
                alpha=0.05,
                rasterized=True,
                **kwargs,
            )
        for ax in axis.flatten():
            ax.set_xlabel("Mouse (UMIs)")
            ax.set_ylabel("Human (UMIs)")

        for ax in axis[0, :]:
            ax.set_xlim((-(7500 / 30.0), 7500))
            ax.set_ylim((-(7500 / 30.0), 7500))
        for ax in axis[1, :]:
            ax.loglog()
        fig.savefig(
            args.output_prefix
            + f"species_mix.comparison_to_droplet.{label}.{suffix}.svg".replace(
                "..", "."
            ),
            dpi=300,
            bbox_inches="tight",
        )

    fig, axis = plt.subplots(3, 3, figsize=(3 * 3, 3 * 3), tight_layout=True)
    attr, label, kwargs = (
        "doublet",
        "coloured_by_doublet",
        {"cmap": "tab20b", "vmin": 0, "vmax": 1},
    )
    for i, t in enumerate(pd.np.arange(0.6, 0.8, 0.1)):
        # get only duplets from 10X
        p_10x = d.loc[d.loc[:, "ratio"] < t, :]
        # get scifi cells which were in those duplet droplets
        p_spc = c.loc[c.index.get_level_values(1).isin(p_10x.index)]
        axis[i, 0].set_title(f"10X\n(duplets only at {t:.2f} purity)")
        axis[i, 0].scatter(
            p_10x["mouse"],
            p_10x["human"],
            c=p_10x[attr] + 1,
            s=1,
            alpha=0.1,
            rasterized=True,
            **kwargs,
        )
        axis[i, 1].set_title(
            f"scifi-RNA-seq\n(only droplets which are duplets in 10X at {t:.2f} purity)"
        )
        axis[i, 1].scatter(
            p_spc["mouse"],
            p_spc["human"],
            c=p_spc[attr] + 1,
            s=1,
            alpha=0.1,
            rasterized=True,
            **kwargs,
        )
        axis[i, 2].set_title(
            f"scifi-RNA-seq\n(only droplets which are duplets in 10X at {t:.2f} purity)"
        )
        axis[i, 2].scatter(
            p_10x["mouse"],
            p_10x["human"],
            color=plt.get_cmap("Pastel1")(0),
            s=1,
            alpha=0.1,
            rasterized=True,
            label="10X",
        )
        axis[i, 2].scatter(
            p_spc["mouse"],
            p_spc["human"],
            color=plt.get_cmap("Pastel1")(1),
            s=1,
            alpha=0.1,
            rasterized=True,
            label="scifi-RNA-seq",
        )
        axis[i, 2].legend()

    for ax in axis.flatten():
        ax.set_xlim((-(3000 / 30.0), 3000))
        ax.set_ylim((-(3000 / 30.0), 3000))

    fig.savefig(
        args.output_prefix
        + f"species_mix.comparison_to_only_10X.only_doublets.{label}.{suffix}.svg".replace(
            "..", "."
        ),
        dpi=300,
        bbox_inches="tight",
    )


def sequence_content_analysis(metrics, whitelist, barcode="r1", head=5000000):
    from Bio import motifs
    from Bio.Seq import Seq

    r1 = metrics.head(head).index.get_level_values(barcode).to_series()

    in_ = r1.isin(whitelist.tolist())
    r1_right = r1[in_]
    r1_wrong = r1[~in_]

    motif_right = motifs.create([Seq(x) for x in r1_right if "N" not in x])
    motif_wrong = motifs.create([Seq(x) for x in r1_wrong if "N" not in x])

    r = pd.DataFrame(motif_right.pwm)
    w = pd.DataFrame(motif_wrong.pwm)
    c = np.log2(w / r)

    fig, axis = plt.subplots(1, 3, figsize=(3 * 3, 3))
    kwargs = {"square": True}
    sns.heatmap(r.T, ax=axis[0], **kwargs)
    axis[0].set_title("Correct")
    sns.heatmap(w.T, ax=axis[1], **kwargs)
    axis[1].set_title("Wrong")
    kwargs = {"square": True, "cmap": "RdBu_r", "center": 0, "robust": True}
    sns.heatmap(c.T, ax=axis[2], **kwargs)
    axis[2].set_title("Difference")
    fig.savefig(
        args.output_prefix + f"barcode_{barcode}.sequence_content.svg",
        bbox_inches="tight",
        dpi=300,
    )


def to_pickle(obj, name, array=True, only_array=False):
    print(f"# {time.asctime()} - Saving {name} to pickle.")
    if array:
        pickle.dump(
            obj.values,
            open(args.output_prefix + f"{name}.values.pickle", "wb"),
            protocol=-1,
        )
    if only_array:
        return
    pickle.dump(
        obj, open(args.output_prefix + f"{name}.pickle", "wb"), protocol=-1
    )


def from_pickle(key, array=False):
    print(f"# {time.asctime()} - Loading {key} from pickle.")
    if array:
        return pickle.load(
            open(args.output_prefix + f"{key}.values.pickle", "rb")
        )
    return pickle.load(open(args.output_prefix + f"{key}.pickle", "rb"))

    # df = pickle.load(open(args.output_prefix + "all.pickle", 'rb'))


def cells_per_droplet_stats(cells_per_droplet):
    # Observe statistical properties
    def poisson(k, lamb):

        return np.exp(-lamb) * (lamb ** k / factorial(k))

    def log_linear_poisson(k, lamb):

        return np.log(np.exp(-lamb) * (lamb ** k / factorial(k)))

    cells_per_droplet_counts = cells_per_droplet.value_counts().sort_index()

    lamb, cov_matrix = scipy.optimize.curve_fit(
        log_linear_poisson,
        cells_per_droplet_counts.index.astype(int),
        cells_per_droplet_counts,
    )

    print(lamb)

    fig, axis = plt.subplots(1, 2, figsize=(3 * 2, 3), tight_layout=True)
    bins = cells_per_droplet_counts.shape[0]
    for ax in axis:
        sns.distplot(
            cells_per_droplet, bins=bins, kde=False, ax=ax, label="Real"
        )
    x = np.arange(0, cells_per_droplet_counts.index.max())
    y_hat = scipy.stats.poisson(lamb).pmf(x)
    y_hat *= cells_per_droplet.shape[0]
    for ax in axis:
        ax.plot(
            x + 0.5, y_hat
        )  # the 0.5 is just to center on the middle of the histogram bins

    cpd = scipy.stats.poisson(lamb).rvs(cells_per_droplet.shape[0])
    for ax in axis:
        sns.distplot(
            cpd, bins=bins, kde=False, ax=ax, norm_hist=False, label="Poisson"
        )
    for ax in axis:
        ax.set_xlabel("Cells")
        ax.set_ylabel("Droplets")
        ax.legend()
    axis[1].set_yscale("log")
    axis[1].set_ylim(bottom=0.1)
    fig.savefig(
        args.output_prefix + f"cells_per_droplet.poissonian_properties.svg",
        dpi=300,
        bbox_inches="tight",
    )


def add_colorbar_to_axis(
    collection, label=None, position="right", size="5%", pad=0.05
):

    divider = make_axes_locatable(collection.axes)
    cax = divider.append_axes(position, size=size, pad=pad)
    cbar = plt.colorbar(mappable=collection, cax=cax, label=label, alpha=1)
    cbar.solids.set_edgecolor("face")
    cbar.solids.set_rasterized(True)


def lin_func(x, m, b):
    return m * x + b


def get_custom_cmap(vmin=-1.0, vmax=1.0, inner=0.4):
    norm = matplotlib.colors.Normalize(-1, 1)
    colors = [
        [norm(vmin), "cornflowerblue"],
        [norm(-inner), "crimson"],
        [norm(inner), "crimson"],
        [norm(vmax), "burlywood"],
    ]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
    cmap.name = "custom"
    return cmap


def parallel_groupby_apply(df, levels, function):

    g = df.groupby(levels)
    res = Parallel(n_jobs=multiprocessing.cpu_count())(
        delayed(function)(group) for name, group in g
    )
    return pd.concat(res)


if __name__ == "__main__":
    sns.set(
        context="paper", style="ticks", palette="colorblind", color_codes=True
    )
    matplotlib.rcParams["svg.fonttype"] = "none"
    # Don't use LaTeX for rendering
    matplotlib.rcParams["text.usetex"] = False
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        sys.exit(1)
