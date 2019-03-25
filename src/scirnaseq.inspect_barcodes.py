#!/usr/bin/env python

"""
sciRNA-seq barcode inspecting script.
"""


import sys
import os
import time
import re
from argparse import ArgumentParser
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
import scipy


pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc("text", usetex=False)

sns.set_style("ticks")


def sorted_nicely(l):
    """ Sort a given iterable in the way that humans expect."""

    def convert(text):
        return int(text) if text.isdigit() else text

    def alphanum_key(key):
        return [convert(c) for c in re.split("([0-9]+)", key)]

    return sorted(l, key=alphanum_key)


def get_cli_arguments():
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="python inspect_barcodes.py",
        description="\n".join(
            [
                "sciRNA-seq script from Bock lab. "
                "See https://github.com/epigen/sciRNA-seq "
                "for specific documentation."
            ]
        ),
    )
    parser.add_argument(
        dest="run_name",
        help="Name of run of experiment. Must match produced input files e.g. 'sciRNA_016'",
        type=str,
    )
    parser.add_argument(
        "-a",
        "--annotation",
        dest="annotation",
        help="Path to CSV file with barcode annotation. "
        "e.g. metadata/sciRNA-seq.oligos_2018-05-22v2.csv",
        type=str,
    )
    default = ["round1", "round2"]
    parser.add_argument(
        "-b",
        "--barcodes",
        dest="cell_barcodes",
        help="Cell barcodes used. A comma-delimited (no spaces) list. "
        " Default is '{}'.".format("','".join(default)),
        type=str,
    )
    default = "results"
    parser.add_argument(
        "-o",
        "--output-dir",
        dest="output_dir",
        help="Output directory. Default is '{}'.".format(default),
        default=default,
        type=str,
    )
    default = 0
    parser.add_argument(
        "--min-mismatches",
        dest="min_mismatches",
        help="Minimum mismatches input files were allowed to have. "
        "Default {}".format(default),
        default=default,
        type=int,
    )
    default = 1
    parser.add_argument(
        "--max-mismatches",
        dest="max_mismatches",
        help="Maximum mismatches input files were allowed to have. "
        "Default {}".format(default),
        default=default,
        type=int,
    )
    default = "metadata/737K-cratac-v1.reverse_complement.csv"
    parser.add_argument(
        "--10x-whitelist",
        dest="bc_10x",
        help="Path to text file with 10X barcode whitelist. ",
        default=default,
        type=str,
    )
    parser.add_argument(
        "--only-whitelist",
        dest="only_whitelist",
        action="store_true",
        help="Whether to look only at barcodes present in the 10X barcode whitelist. "
        " Will speed processing dramatically.",
    )
    parser.add_argument(
        "--10x-reverse-complement",
        dest="bc_10x_rc",
        help="Whether 10X barcodes should be reverse complemented.",
        action="store_true"
    )
    parser.add_argument(
        "--estimated-cells",
        dest="estimated_cells",
        help="Estimated number of cells in experiment. "
        "Only used for visualization.",
        type=int,
        default=10000,
    )

    # args = parser.parse_args("-a metadata/sciRNA-seq.oligos_2018-11-14.csv -b round1,round2 --max-mismatches 0 -".split(" "))
    args = parser.parse_args(
        "-a metadata/sciRNA-seq.SCI017.oligos_2019-02-11.csv --estimated-cells 4000 --only-whitelist -b round1,round2 --max-mismatches 0 sci-RNA-seq_SCI019".split(
            " "
        )
    )
    print("# " + time.asctime() + " - Parsing command line arguments:")
    # args = parser.parse_args()
    args.cell_barcodes = args.cell_barcodes.split(",")
    print(args)
    return args


def main():
    print("# " + time.asctime() + " - Start.")

    args = get_cli_arguments()

    annotation = pd.read_csv(args.annotation)
    bc = pd.read_csv(args.bc_10x)
    if "reverse_complement" not in bc.columns:
        from Bio.Seq import Seq

        bc.loc[:, "reverse_complement"] = [
            str(Seq(x).reverse_complement()) for x in bc["original"]
        ]

    # # get barcodes with well information
    barcodes_with_well = list()
    for barcode in annotation["barcode_type"].unique():
        if (
            not annotation.loc[annotation["barcode_type"] == barcode, "plate_well"]
            .isnull()
            .all()
        ):
            barcodes_with_well.append(barcode)

    for mismatches in range(args.min_mismatches, args.max_mismatches):
        print(
            "# "
            + time.asctime()
            + " - Processing files with {} mismatches.".format(mismatches)
        )

        output_prefix = "{}.{}mis.".format(args.run_name, mismatches)
        print(output_prefix)

        # get files to concatenate into HDF
        hdf_file = os.path.join(
            "barcodes", args.run_name + ".mis{}.hdf".format(mismatches)
        )
        umi_count_file = os.path.join(
            "barcodes", output_prefix + "barcode_gene_umi_count.clean.csv"
        )
        umi_dup_file = os.path.join(
            "barcodes", output_prefix + "barcode_umi_dups.count.csv"
        )
        cell_dup_file = os.path.join(
            "barcodes", output_prefix + "barcode_umi_dups.per_cell.csv"
        )
        mapping_cell_rate_file = os.path.join(
            "barcodes", output_prefix + "mapping_rate.per_cell.csv"
        )
        mapping_well_rate_file = os.path.join(
            "barcodes", output_prefix + "mapping_rate.per_well.csv"
        )
        cell_metrics_file = os.path.join(
            "barcodes", output_prefix + "all_metrics.per_cell.csv"
        )

        # start
        print("# " + time.asctime() + " - Concatenating barcode files into HDF file.")
        pieces = glob(
            os.path.join(
                "barcodes",
                "{}.part_*.barcodes.*_*.mis_{}.csv.gz".format(
                    args.run_name, mismatches
                ),
            )
        )
        pieces = sorted_nicely(pieces)

        print(
            "## "
            + time.asctime()
            + " - barcode files to read: '{}'.".format(", ".join(pieces))
        )
        print(pieces[0])
        c = pd.read_csv(pieces[0], compression="gzip")
        c.to_hdf(hdf_file, "cells", mode="w", format="table")
        del c  # allow garbage to be collected

        for piece in pieces[1:]:
            print(piece)
            c = pd.read_csv(piece, compression="gzip")
            c.to_hdf(hdf_file, "cells", append=True)
            del c  # allow garbage to be collected

        # read up HDF
        print("# " + time.asctime() + " - Reading up HDF file.")
        cells = pd.read_hdf(hdf_file, "cells")

        # read transcriptome
        print("# " + time.asctime() + " - Concatenating transcriptome files.")
        transcriptome = pd.DataFrame()
        pieces = sorted_nicely(
            glob(
                os.path.join(
                    "star",
                    "{}.STAR.htseq-count_gene.read_gene.part_*.csv".format(
                        args.run_name
                    ),
                )
            )
        )
        for piece in pieces:
            print(piece)
            t = pd.read_csv(piece, header=None, names=["read", "gene"])
            transcriptome = transcriptome.append(t, ignore_index=True)

        msg = " - total reads in experiment: '{}'".format(transcriptome.shape[0])
        print("## " + time.asctime() + msg)

        # join barcodes and transcriptome
        print("# " + time.asctime() + " - Joining barcodes and transcriptome.")
        df = cells.set_index("read").join(transcriptome.set_index("read"))

        msg = " - total with correct combinatorial barcodes: '{}'".format(df.shape[0])
        print("## " + time.asctime() + msg)

        # Filter on 10X whitelist if requested
        if args.only_whitelist:
            if args.bc_10x_rc:
                ll = "reverse_complement"
            else:
                ll = "original"
            i = df["round2"].isin(bc[ll].tolist())
            ic = i.sum()
            msg = " - reads with correct 10X barcodes: '{}'".format(ic)
            print("## " + time.asctime() + msg)
            msg = " - fraction of correct 10X barcodes: '{}'".format(
                ic / float(df.shape[0])
            )
            print("## " + time.asctime() + msg)
            df = df.loc[i.isin([True]), :]

        # investigate mapping rates
        print("# " + time.asctime() + " - Investigating mapping rates per cell.")

        # # per cell
        reads_per_cell = df.groupby(args.cell_barcodes)[
            "gene"
        ].count()  # gene here is just a dummy
        df.loc[:, "unmapped"] = df["gene"].str.startswith("__").astype(int)
        # unmapped = df.groupby(args.cell_barcodes)['gene'].apply(lambda x: x.str.startswith("__").sum())
        unmapped = df.groupby(args.cell_barcodes)["unmapped"].sum()
        mapping_rate = 1 - (unmapped / reads_per_cell)

        cell_rates = pd.DataFrame(
            [reads_per_cell, unmapped, mapping_rate],
            index=["all_reads", "unmapped", "mapping_rate"],
        ).T
        cell_rates.to_csv(mapping_cell_rate_file)

        # Plot read distribution
        print("# " + time.asctime() + " - Investigating read distributions.")
        fig, axis = plt.subplots(2, 2, figsize=(2 * 4, 2 * 4))
        reads_per_cell = reads_per_cell.sort_values()
        fig.suptitle(args.run_name)
        r = reads_per_cell.rank(ascending=False)
        axis[0, 0].plot(r, reads_per_cell)
        if args.estimated_cells is not None:
            axis[0, 1].plot(
                r.tail(args.estimated_cells), reads_per_cell.tail(args.estimated_cells)
            )
            axis[1, 1].plot(
                r.tail(args.estimated_cells), reads_per_cell.tail(args.estimated_cells)
            )
        axis[1, 0].plot(r, reads_per_cell)
        axis[1, 0].set_xscale("log")
        axis[1, 0].set_yscale("log")
        axis[1, 1].set_xscale("log")
        axis[1, 1].set_yscale("log")
        for ax in axis.flatten():
            ax.axvline(args.estimated_cells, linestyle="--", color="black", alpha=0.5)
            ax.set_xlabel("Cells")
            ax.set_ylabel("Reads")
        sns.despine(fig)
        fig.savefig(
            os.path.join(
                args.output_dir, output_prefix + "reads_per_cell.lineplot.svg"
            ),
            bbox_inches="tight",
        )

        # # per well
        print("# " + time.asctime() + " - Investigating mapping rates per well.")
        # # # first add well information
        for barcode in barcodes_with_well:
            a = annotation.loc[
                annotation["barcode_type"] == barcode,
                ["barcode_sequence", "plate_well"],
            ].rename(
                columns={"barcode_sequence": barcode, "plate_well": barcode + "_well"}
            )
            df = pd.merge(df, a, on=barcode, how="left")

        all_reads = df.groupby(["round1_well"])["gene"].count()
        unmapped = df.groupby(["round1_well"])["gene"].apply(
            lambda x: x.str.startswith("__").sum()
        )
        mapping_rate = 1 - (unmapped / all_reads)

        well_rates = pd.DataFrame(
            [all_reads, unmapped, mapping_rate],
            index=["all_reads", "unmapped", "mapping_rate"],
        ).T
        well_rates.to_csv(mapping_well_rate_file)

        # remove unmapped reads
        msg = " - Removing unmapped reads and collapsing to unique UMIs per cell per gene."
        print("# " + time.asctime() + msg)
        df2 = df.loc[df["unmapped"] != 1, :]

        print("# " + time.asctime() + " - Investigating duplicates.")
        duplicates_per_molecule = df2.groupby(args.cell_barcodes)["umi"].value_counts()
        duplicates_per_molecule.columns = args.cell_barcodes + ["umi", "count"]
        duplicates_per_molecule.to_csv(umi_dup_file)

        fig, axis = plt.subplots(1, 2, figsize=(2 * 3, 3), tight_layout=True)
        sns.distplot(duplicates_per_molecule, ax=axis[0], kde=False)
        axis[0].set_xlabel("Reads per UMI")
        axis[0].set_ylabel("UMIs (log)")
        axis[0].set_yscale("log")
        sns.distplot(np.log10(duplicates_per_molecule), ax=axis[1], kde=False)
        axis[1].set_xlabel("Reads per UMI (log)")
        axis[1].set_ylabel("UMIs (log)")
        # axis[1].set_xscale("log")
        axis[1].set_yscale("log")
        reads_per_umi_plot = os.path.join(
            "results", output_prefix + "barcode_umi_dups.per_cell.distplot.svg"
        )
        fig.savefig(reads_per_umi_plot, dpi=300, bbox_inches="tight")

        unique_molecules = (duplicates_per_molecule == 1).astype(int)
        duplicates_per_cell = (
            unique_molecules.groupby(args.cell_barcodes).sum().to_frame(name="unique")
        )
        duplicates_per_cell.loc[:, "all"] = duplicates_per_molecule.groupby(
            args.cell_barcodes
        ).sum()

        duplicates_per_cell["unique_rate"] = (
            duplicates_per_cell["unique"] / duplicates_per_cell["all"]
        )
        duplicates_per_cell.to_csv(cell_dup_file)

        fig, axis = plt.subplots(1, 1, figsize=(1 * 3, 3), tight_layout=True)
        sns.distplot(duplicates_per_cell["unique_rate"], ax=axis, kde=False)
        axis.set_xlabel("Unique rate per Cell")
        axis.set_ylabel("Cells")
        axis.set_yscale("log")
        cell_unique_rate_plot = os.path.join(
            "results", output_prefix + "cell_umi_dups.per_cell.distplot.svg"
        )
        fig.savefig(cell_unique_rate_plot, dpi=300, bbox_inches="tight")

        # get unique UMIs (UMI counts per cell per gene)
        df3 = df2.groupby(args.cell_barcodes + ["gene"])["umi"].nunique()

        print("# " + time.asctime() + " - Writing gene expression matrix to file.")
        df3.to_csv(umi_count_file, index=True, header=True)

        # duplicates_per_molecule = pd.read_csv(umi_dup_file, index_col=[0, 1])
        # duplicates_per_cell = pd.read_csv(cell_dup_file, index_col=[0, 1])
        # cell_rates = pd.read_csv(mapping_cell_rate_file, index_col=[0, 1])
        # df3 = pd.read_csv(umi_count_file)

        umis_per_cell = df3.reset_index().groupby(args.cell_barcodes)["umi"].sum()
        genes_per_cell = df3.reset_index().groupby(args.cell_barcodes)["gene"].nunique()

        # Plot UMI distribution
        print("# " + time.asctime() + " - Investigating UMI distributions.")
        fig, axis = plt.subplots(2, 2, figsize=(2 * 4, 2 * 4))
        umis_per_cell = umis_per_cell.sort_values()
        fig.suptitle(args.run_name)
        r = umis_per_cell.rank(ascending=False)
        axis[0, 0].plot(r, umis_per_cell)
        axis[0, 1].plot(
            r.tail(args.estimated_cells), umis_per_cell.tail(args.estimated_cells)
        )
        axis[1, 0].plot(r, umis_per_cell)
        axis[1, 0].set_xscale("log")
        axis[1, 0].set_yscale("log")
        axis[1, 1].plot(
            r.tail(args.estimated_cells), umis_per_cell.tail(args.estimated_cells)
        )
        axis[1, 1].set_xscale("log")
        axis[1, 1].set_yscale("log")
        for ax in axis.flatten():
            ax.axvline(args.estimated_cells, linestyle="--", color="black", alpha=0.5)
            ax.set_xlabel("Cells")
            ax.set_ylabel("UMIs")
        sns.despine(fig)
        fig.savefig(
            os.path.join(args.output_dir, output_prefix + "umis_per_cell.lineplot.svg"),
            bbox_inches="tight",
        )

        # plot again all metrics, including per-cell duplicate rate
        d = pd.concat(
            [umis_per_cell, genes_per_cell, duplicates_per_cell.loc[:, "unique_rate"]],
            ignore_index=True,
            axis=1,
        )
        p = pd.concat(
            [cell_rates.loc[d.index, :].dropna(), d], ignore_index=True, axis=1
        )
        p.columns = cell_rates.columns.tolist() + ["umi", "gene", "unique_rate"]
        p.to_csv(cell_metrics_file)

        fig, axis = plt.subplots(2, 3, figsize=(3 * 3, 3 * 2), tight_layout=True)

        sns.distplot(p["mapping_rate"], ax=axis[0, 0], kde=False)
        axis[0, 0].set_xlabel("Mapping rate")
        axis[0, 0].set_ylabel("Barcodes")
        axis[0, 0].set_yscale("log")

        axis[1, 0].scatter(
            p["all_reads"].head(int(8 * 50e4)),
            p["gene"].head(int(8 * 50e4)),
            alpha=0.2,
            s=1,
            rasterized=True,
        )
        axis[1, 0].set_xlabel("Reads per barcode")
        axis[1, 0].set_ylabel("Genes per barcode")
        axis[1, 0].set_xscale("log")
        axis[1, 0].set_yscale("log")

        axis[0, 1].scatter(
            p["all_reads"].head(int(8 * 50e4)),
            p["mapping_rate"].head(int(8 * 50e4)),
            alpha=0.2,
            s=1,
            rasterized=True,
        )
        axis[0, 1].set_xlabel("Reads per barcode")
        axis[0, 1].set_ylabel("Mapping rate")
        axis[0, 1].set_xscale("log")

        axis[0, 2].scatter(
            p["all_reads"].head(int(8 * 50e4)),
            p["unique_rate"].head(int(8 * 50e4)),
            alpha=0.2,
            s=1,
            rasterized=True,
        )
        axis[0, 2].set_xscale("log")
        axis[0, 2].set_xlabel("Reads per barcode")
        axis[0, 2].set_ylabel("Unique rate")

        sns.distplot(np.log10(p["all_reads"]), ax=axis[1, 1], kde=False)
        axis[1, 1].set_xlabel("Reads per barcode (log10)")
        axis[1, 1].set_yscale("log")
        axis[1, 1].set_ylabel("Barcodes")

        axis[1, 2].scatter(
            p["mapping_rate"].head(int(8 * 50e4)),
            p["unique_rate"].head(int(8 * 50e4)),
            alpha=0.2,
            s=1,
            rasterized=True,
        )
        axis[1, 2].set_xlabel("Mapping rate")
        axis[1, 2].set_ylabel("Unique rate")

        reads_vs_mapping_rate_plot = os.path.join(
            "results",
            output_prefix
            + "reads_vs_mapping_rate_vs_duplication.per_cell.scatter.res.svg",
        )
        fig.savefig(reads_vs_mapping_rate_plot, dpi=300, bbox_inches="tight")

        # Pairwise plots

        # # add type of material
        bc_annot = annotation.set_index("barcode_sequence").loc[:, "material"]
        bc_annot.index.name = "round2"
        p_annot = p.join(bc_annot, on="round1")

        p2 = p_annot.loc[p_annot["umi"].nlargest(args.estimated_cells).index]
        g = sns.pairplot(
            p2,
            hue="material",
            plot_kws={
                "linewidth": 0,
                "edgecolors": "none",
                "rasterized": True,
                "alpha": 0.2,
            },
        )
        cell_unique_rate_plot = os.path.join(
            "results", output_prefix + "cell_umi_dups.per_cell.top_cells.pairplot.svg"
        )
        g.fig.savefig(cell_unique_rate_plot, dpi=300, bbox_inches="tight")

        p2 = p_annot.loc[p_annot["umi"].nlargest(args.estimated_cells).index]
        p2.loc[:, ~p2.columns.str.contains("_rate|material")] = np.log10(
            p2.loc[:, ~p2.columns.str.contains("_rate|material")] + 1
        )
        g = sns.pairplot(
            p2,
            hue="material",
            plot_kws={
                "linewidth": 0,
                "edgecolors": "none",
                "rasterized": True,
                "alpha": 0.2,
            },
        )
        cell_unique_rate_plot = os.path.join(
            "results",
            output_prefix + "cell_umi_dups.per_cell.top_cells.log.pairplot.svg",
        )
        g.fig.savefig(cell_unique_rate_plot, dpi=300, bbox_inches="tight")

        p2 = p_annot.loc[p_annot["umi"].nsmallest(args.estimated_cells).index]
        g = sns.pairplot(
            p2,
            hue="material",
            plot_kws={
                "linewidth": 0,
                "edgecolors": "none",
                "rasterized": True,
                "alpha": 0.2,
            },
        )
        cell_unique_rate_plot = os.path.join(
            "results",
            output_prefix + "cell_umi_dups.per_cell.bottom_cells.pairplot.svg",
        )
        g.fig.savefig(cell_unique_rate_plot, dpi=300, bbox_inches="tight")

        p2 = p_annot.loc[p_annot["umi"].nsmallest(args.estimated_cells).index]
        p2.loc[:, ~p2.columns.str.contains("_rate|material")] = np.log10(
            p2.loc[:, ~p2.columns.str.contains("_rate|material")] + 1
        )
        g = sns.pairplot(
            p2,
            hue="material",
            plot_kws={
                "linewidth": 0,
                "edgecolors": "none",
                "rasterized": True,
                "alpha": 0.2,
            },
        )
        cell_unique_rate_plot = os.path.join(
            "results",
            output_prefix + "cell_umi_dups.per_cell.bottom_cells.log.pairplot.svg",
        )
        g.fig.savefig(cell_unique_rate_plot, dpi=300, bbox_inches="tight")

        # Species mixing experiment
        print("# " + time.asctime() + " - Species mixing experiment.")

        # # add species ratio and assigned species
        if isinstance(df3, pd.Series):
            df3 = df3.reset_index()
        df3.loc[:, "human"] = df3.loc[:, "gene"].str.startswith("ENSG")

        sp_ratio = df3.pivot_table(
            index=args.cell_barcodes,
            columns="human",
            values="umi",
            fill_value=0,
            aggfunc=sum,
        )
        sp_ratio.columns = ["mouse", "human"]
        sp_ratio.loc[:, "total_umis"] = sp_ratio.sum(axis=1)
        # # normalize to mean per species
        sp_ratio = sp_ratio.assign(
            mouse_norm=sp_ratio["mouse"] / sp_ratio["mouse"].mean()
        )
        sp_ratio = sp_ratio.assign(
            human_norm=sp_ratio["human"] / sp_ratio["human"].mean()
        )

        sp_ratio.loc[:, "ratio"] = sp_ratio.loc[:, "human"] / sp_ratio.loc[:, "mouse"]
        sp_ratio.loc[:, "ratio"] = sp_ratio.loc[:, "ratio"]
        m = sp_ratio.loc[:, "ratio"].replace(np.inf, 0).replace(-np.inf, 0).abs().max()
        sp_ratio.loc[:, "ratio"] = (
            sp_ratio.loc[:, "ratio"].replace(np.inf, m).replace(-np.inf, 0)
        )
        sp_ratio.loc[:, "log_ratio"] = np.log2(sp_ratio.loc[:, "ratio"]).replace(
            -np.inf, -np.log2(m)
        )
        sp_ratio.loc[:, "species"] = np.nan
        sp_ratio.loc[sp_ratio["log_ratio"] > 1.5, "species"] = "human"
        sp_ratio.loc[sp_ratio["log_ratio"] < -1.5, "species"] = "mouse"
        sp_ratio.loc[:, "species"] = sp_ratio["species"].fillna("doublet")

        msg = " - Doublet rate: {}".format(
            (sp_ratio["species"] == "doublet").sum() / float(sp_ratio.shape[0])
        )
        print("# " + time.asctime() + msg)

        q = sp_ratio["species"].value_counts()
        msg = " - Human to mouse ratio: {}".format(q["human"] / float(q["mouse"]))
        print("# " + time.asctime() + msg)

        frac = dict()
        for x in np.linspace(0, 8, 20):
            sp_ratio.loc[:, "species"] = np.nan
            sp_ratio.loc[sp_ratio["log_ratio"] > x, "species"] = "human"
            sp_ratio.loc[sp_ratio["log_ratio"] < -x, "species"] = "mouse"
            sp_ratio.loc[:, "species"] = sp_ratio.loc[:, "species"].fillna("doublet")

            q = sp_ratio.species.value_counts()
            frac[x] = q["doublet"] / float(q.sum())

        doublet_fraction = pd.Series(frac, name="doublet_fraction")

        fig, axis = plt.subplots(1, 2, figsize=(2 * 3, 3))
        for ax in axis:
            ax.plot(doublet_fraction.index, doublet_fraction)
            ax.set_xlabel("Purity required (log fold)")
        axis[0].set_ylim((0, 1))
        axis[0].set_ylabel("Fraction of doublets")
        axis[1].set_ylabel("Fraction of doublets (log)")
        axis[1].set_yscale("log")
        sns.despine(fig)
        fig.savefig(
            os.path.join(
                args.output_dir,
                output_prefix + "species_mix.summary.dependent_purity.svg",
            ),
            bbox_inches="tight",
            dpi=300,
        )

        # # Summary stats
        fig, axis = plt.subplots(1, 2, figsize=(2 * 3, 1 * 3))
        fig.suptitle(args.run_name)
        s = sp_ratio["species"].value_counts()
        sns.barplot(s.index, s, ax=axis[0])
        axis[0].set_xlabel("All barcodes")
        axis[0].set_ylabel("Barcodes")

        s = sp_ratio.loc[sp_ratio["total_umis"] > 100, "species"].value_counts()
        sns.barplot(s.index, s, ax=axis[1])
        axis[1].set_xlabel("Cells with at least 100 UMI")
        axis[1].set_ylabel("Barcodes")
        sns.despine(fig)
        fig.savefig(
            os.path.join(
                args.output_dir, output_prefix + "species_mix.summary.barplot.svg"
            ),
            bbox_inches="tight",
            dpi=300,
        )

        # # Logster plots
        fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 4))
        fig.suptitle(args.run_name)
        for ax in axis:
            ax.axhline(0, linestyle="--", color="black", alpha=0.1, zorder=100)
            ax.axvline(0, linestyle="--", color="black", alpha=0.1, zorder=100)
        m = max(sp_ratio.loc[:, "mouse"].max(), sp_ratio.loc[:, "human"].max())
        axis[0].plot((0, m), (0, m), linestyle="--", alpha=0.5, color="black")
        axis[0].scatter(
            sp_ratio.loc[:, "mouse"],
            sp_ratio.loc[:, "human"],
            c=sp_ratio.loc[:, "log_ratio"],
            cmap="coolwarm",
            s=2,
            alpha=0.8,
            rasterized=True,
        )
        axis[0].set_xlim((0, m))
        axis[0].set_ylim((0, m))
        a = np.log2(1 + sp_ratio.loc[:, "mouse"])
        b = np.log2(1 + sp_ratio.loc[:, "human"])
        m = max(a.max(), b.max())
        axis[1].plot((0, m), (0, m), linestyle="--", alpha=0.5, color="black")
        axis[1].scatter(
            a,
            b,
            c=sp_ratio.loc[:, "log_ratio"],
            cmap="coolwarm",
            s=2,
            alpha=0.8,
            rasterized=True,
        )
        axis[1].set_xlim((0, m))
        axis[1].set_ylim((0, m))
        axis[0].set_xlabel("Mouse (UMIs)")
        axis[0].set_ylabel("Human (UMIs)")
        axis[1].set_xlabel("Mouse (log2 UMIs)")
        axis[1].set_ylabel("Human (log2 UMIs)")
        sns.despine(fig)
        fig.savefig(
            os.path.join(args.output_dir, output_prefix + "species_mix.plot.svg"),
            bbox_inches="tight",
            dpi=300,
        )

        # Efficiency plot
        def f(x, m, b):
            return m * x + b

        p2 = p.loc[p["umi"].nlargest(100000).index]

        fig, axis = plt.subplots(2, 2, figsize=(2 * 3, 2 * 3), tight_layout=True)
        fig.suptitle("Experiment efficiency:\n" + args.run_name, ha="center")

        for j, var in enumerate(["umi", "gene"]):
            for i, l in [(0, ""), (1, " (log)")]:
                (m, b), pcov = scipy.optimize.curve_fit(f, p2["all_reads"], p2[var])
                # lm = LinearRegression().fit(p2['all_reads'].values.reshape(-1, 1), p2[var])
                # assert np.allclose(m, lm.coef_)
                # assert np.allclose(b, lm.intercept_)
                axis[i, j].text(
                    0.5e5,
                    args.estimated_cells,
                    s="y = {:.6f}x + {:.2f}".format(m, b),
                    ha="center",
                )

                axis[i, j].set_xlabel("Total reads sequenced per cell" + l)
                axis[i, j].set_ylabel(
                    "Useful {}s per cell".format(var.capitalize()) + l
                )
                axis[i, j].scatter(
                    p2["all_reads"],
                    p2[var],
                    alpha=0.2,
                    s=2,
                    edgecolors="none",
                    rasterized=True,
                )
                x = np.linspace(p2["all_reads"].min(), p2["all_reads"].max(), num=1000)
                axis[i, j].plot(x, f(x, m, b), color="orange")
                y = f(1e4, m, b)
                axis[i, j].text(
                    1e4,
                    y,
                    s="{}s recovered\nwith 10.000\nreads sequenced:\n{:.2f}".format(
                        var.capitalize(), y
                    ),
                    ha="left",
                )
                # X == Y
                xmax = p2["all_reads"].max()
                xmax += xmax * 0.1
                ymax = p2[var].max()
                ymax += ymax * 0.1
                x = np.linspace(0, ymax, num=2)
                y = f(x, 1, 0)
                axis[i, j].plot(x, y, linestyle="--", color="black", linewidth=0.5)

                # Axlines
                for h in [100, 250, 500, 1000, args.estimated_cells]:
                    axis[i, j].axhline(h, linestyle="--", color="black", linewidth=0.5)
                for v in [10000, 100000]:
                    axis[i, j].axvline(v, linestyle="--", color="black", linewidth=0.5)

                axis[i, j].ticklabel_format(useOffset=False, style="plain")
                if i == 1:
                    axis[i, j].loglog()
        performance_plot = os.path.join(
            "results", output_prefix + "performance_plot.only_10x_correct.svg"
        )
        fig.savefig(performance_plot, dpi=300, bbox_inches="tight")

        # Compare benefit of 1 vs 2 rounds
        r2_only = (
            df3.groupby([args.cell_barcodes[1]] + ["gene", "human"])["umi"]
            .sum()
            .reset_index()
        )

        r2_sp_ratio = r2_only.pivot_table(
            index=args.cell_barcodes[1],
            columns="human",
            values="umi",
            fill_value=0,
            aggfunc=sum,
        )
        r2_sp_ratio.columns = ["mouse", "human"]
        r2_sp_ratio.loc[:, "total_umis"] = r2_sp_ratio.sum(axis=1)
        # # normalize to mean per species
        r2_sp_ratio = r2_sp_ratio.assign(
            mouse_norm=r2_sp_ratio["mouse"] / r2_sp_ratio["mouse"].mean()
        )
        r2_sp_ratio = r2_sp_ratio.assign(
            human_norm=r2_sp_ratio["human"] / r2_sp_ratio["human"].mean()
        )

        r2_sp_ratio.loc[:, "ratio"] = (
            r2_sp_ratio.loc[:, "human"] / r2_sp_ratio.loc[:, "mouse"]
        )
        r2_sp_ratio.loc[:, "ratio"] = r2_sp_ratio.loc[:, "ratio"]
        m = (
            r2_sp_ratio.loc[:, "ratio"]
            .replace(np.inf, 0)
            .replace(-np.inf, 0)
            .abs()
            .max()
        )
        r2_sp_ratio.loc[:, "ratio"] = (
            r2_sp_ratio.loc[:, "ratio"].replace(np.inf, m).replace(-np.inf, 0)
        )
        r2_sp_ratio.loc[:, "log_ratio"] = np.log2(r2_sp_ratio.loc[:, "ratio"]).replace(
            -np.inf, -np.log2(m)
        )
        r2_sp_ratio.loc[:, "species"] = np.nan
        r2_sp_ratio.loc[r2_sp_ratio["log_ratio"] > 3, "species"] = "human"
        r2_sp_ratio.loc[r2_sp_ratio["log_ratio"] < -3, "species"] = "mouse"
        r2_sp_ratio.loc[:, "species"] = r2_sp_ratio.loc[:, "species"].fillna("doublet")

        msg = " - Doublet rate: {}".format(
            (r2_sp_ratio.species == "doublet").sum() / float(r2_sp_ratio.shape[0])
        )
        print("# " + time.asctime() + msg)

        q = r2_sp_ratio["species"].value_counts()
        msg = " - Human to mouse ratio: {}".format(q["human"] / float(q["mouse"]))
        print("# " + time.asctime() + msg)

        # # Summary stats
        fig, axis = plt.subplots(1, 2, figsize=(2 * 3, 1 * 3))
        fig.suptitle(args.run_name)
        s = r2_sp_ratio["species"].value_counts()
        sns.barplot(s.index, s, ax=axis[0])
        axis[0].set_xlabel("All barcodes")
        axis[0].set_ylabel("Barcodes")

        s = r2_sp_ratio.loc[r2_sp_ratio["total_umis"] > 100, "species"].value_counts()
        sns.barplot(s.index, s, ax=axis[1])
        axis[1].set_xlabel("Cells with at least 100 UMI")
        axis[1].set_ylabel("Barcodes")
        sns.despine(fig)
        fig.savefig(
            os.path.join(
                args.output_dir,
                output_prefix + "only_10X_barcodes.species_mix.summary.barplot.svg",
            ),
            bbox_inches="tight",
            dpi=300,
        )

        # # Logster plots
        colors = sns.color_palette("colorblind")
        fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 4))
        fig.suptitle(args.run_name)
        for ax in axis:
            ax.axhline(0, linestyle="--", color="black", alpha=0.1, zorder=100)
            ax.axvline(0, linestyle="--", color="black", alpha=0.1, zorder=100)
        m = max(
            sp_ratio.loc[:, "mouse"].max(),
            r2_sp_ratio.loc[:, "mouse"].max(),
            sp_ratio.loc[:, "human"].max(),
            r2_sp_ratio.loc[:, "human"].max(),
        )
        axis[0].plot((0, m), (0, m), linestyle="--", alpha=0.5, color="black")
        axis[0].scatter(
            sp_ratio.loc[:, "mouse"],
            sp_ratio.loc[:, "human"],
            s=3,
            alpha=0.3,
            rasterized=True,
            color=colors[0],
            label="round1 + 10X",
        )
        axis[0].scatter(
            r2_sp_ratio.loc[:, "mouse"],
            r2_sp_ratio.loc[:, "human"],
            s=3,
            alpha=0.3,
            rasterized=True,
            color=colors[1],
            label="10X only",
        )
        axis[0].legend()
        axis[0].set_xlim((-(m * 0.05), m))
        axis[0].set_ylim((-(m * 0.05), m))
        axis[0].set_xlabel("Mouse (UMIs)")
        axis[0].set_ylabel("Human (UMIs)")

        a = np.log2(1 + sp_ratio.loc[:, "mouse"])
        b = np.log2(1 + sp_ratio.loc[:, "human"])
        r2a = np.log2(1 + r2_sp_ratio.loc[:, "mouse"])
        r2b = np.log2(1 + r2_sp_ratio.loc[:, "human"])
        m = max(r2a.max(), a.max(), r2b.max(), b.max())
        axis[1].plot((0, m), (0, m), linestyle="--", alpha=0.5, color="black")
        axis[1].scatter(
            a, b, s=3, alpha=0.3, rasterized=True, color=colors[0], label="round1 + 10X"
        )
        axis[1].scatter(
            r2a, r2b, s=3, alpha=0.3, rasterized=True, color=colors[1], label="10X only"
        )
        axis[1].legend()
        axis[1].set_xlim((-(m * 0.05), m))
        axis[1].set_ylim((-(m * 0.05), m))
        axis[1].set_xlabel("Mouse (log2 UMIs)")
        axis[1].set_ylabel("Human (log2 UMIs)")
        sns.despine(fig)
        fig.savefig(
            os.path.join(
                args.output_dir,
                output_prefix + "only_10X_barcodes.species_mix.plot.svg",
            ),
            bbox_inches="tight",
            dpi=300,
        )

        # # # Plot barcodes in 10X which are duplets now demultiplexed by round1
        frac = dict()
        for x in np.linspace(0, 8, 20):
            r2_sp_ratio.loc[:, "species"] = np.nan
            r2_sp_ratio.loc[r2_sp_ratio["log_ratio"] > x, "species"] = "human"
            r2_sp_ratio.loc[r2_sp_ratio["log_ratio"] < -x, "species"] = "mouse"
            r2_sp_ratio.loc[:, "species"] = r2_sp_ratio.loc[:, "species"].fillna(
                "doublet"
            )

            q = r2_sp_ratio.species.value_counts()
            frac[x] = q["doublet"] / float(q.sum())

        doublet_fraction = pd.Series(frac, name="doublet_fraction")

        fig, axis = plt.subplots(1, 2, figsize=(2 * 3, 3))
        for ax in axis:
            ax.plot(doublet_fraction.index, doublet_fraction)
            ax.set_xlabel("Purity required (log fold)")
        axis[0].set_ylim((0, 1))
        axis[0].set_ylabel("Fraction of doublets")
        axis[1].set_ylabel("Fraction of doublets (log)")
        axis[1].set_yscale("log")
        sns.despine(fig)
        fig.savefig(
            os.path.join(
                args.output_dir,
                output_prefix
                + ".only_10X_barcodes.species_mix.summary.dependent_purity.svg",
            ),
            bbox_inches="tight",
            dpi=300,
        )

        fig, axis = plt.subplots(5, 2, figsize=(2 * 4, 5 * 4), tight_layout=True)
        fig.suptitle(args.run_name)
        for i, fold in enumerate([0.1, 0.25, 0.5, 1, 2]):
            r2_tmp = r2_sp_ratio.copy()
            r2_tmp = r2_tmp.loc[r2_tmp["total_umis"] >= 50]
            r2_tmp.loc[:, "species"] = np.nan
            r2_tmp.loc[r2_tmp["log_ratio"] > fold, "species"] = "human"
            r2_tmp.loc[r2_tmp["log_ratio"] < -fold, "species"] = "mouse"
            r2_tmp.loc[:, "species"] = r2_tmp.loc[:, "species"].fillna("doublet")

            r1_2_tmp = sp_ratio.copy()
            r1_2_tmp = r1_2_tmp.loc[r1_2_tmp["total_umis"] >= 50]
            r1_2_tmp.loc[:, "species"] = np.nan
            r1_2_tmp.loc[r1_2_tmp["log_ratio"] > fold, "species"] = "human"
            r1_2_tmp.loc[r1_2_tmp["log_ratio"] < -fold, "species"] = "mouse"
            r1_2_tmp.loc[:, "species"] = r1_2_tmp.loc[:, "species"].fillna("doublet")

            r2_doublets = r2_tmp.loc[r2_tmp["species"] == "doublet"].index
            r1_demux = r1_2_tmp.loc[
                r1_2_tmp.index.get_level_values("round2").isin(r2_doublets.tolist())
            ].index

            demux_doublet_fraction = r1_2_tmp.loc[r1_demux, "species"].value_counts()
            demux_doublet_fraction = demux_doublet_fraction["doublet"] / float(
                demux_doublet_fraction.sum()
            )

            a = np.log2(1 + r1_2_tmp.loc[:, "mouse"])
            b = np.log2(1 + r1_2_tmp.loc[:, "human"])
            r2a = np.log2(1 + r2_tmp.loc[:, "mouse"])
            r2b = np.log2(1 + r2_tmp.loc[:, "human"])

            colors = sns.color_palette("colorblind")

            m = max(r2a.max(), a.max(), r2b.max(), b.max())
            for ax in axis[i, :]:
                ax.axhline(0, linestyle="--", color="black", alpha=0.1, zorder=100)
                ax.axvline(0, linestyle="--", color="black", alpha=0.1, zorder=100)
                ax.set_xlim((-(m * 0.05), m))
                ax.set_ylim((-(m * 0.05), m))
                ax.set_xlabel("Mouse (log2 UMIs)")
                ax.set_ylabel("Human (log2 UMIs)")
                ax.plot((0, m), (0, m), linestyle="--", alpha=0.5, color="black")
            axis[i, 0].text(5, 1, s="Purity for doublet selection = {}".format(fold))
            axis[i, 1].text(
                6, 1, s="Doublet fraction = {:2f}".format(demux_doublet_fraction)
            )

            axis[i, 0].set_title("Doublets with 10X")
            axis[i, 0].scatter(
                r2a.loc[r2_doublets],
                r2b.loc[r2_doublets],
                s=5,
                alpha=0.5,
                rasterized=True,
                color=colors[1],
                label="10X only",
            )

            axis[i, 1].set_title("Demultiplexed with scifi-RNA-seq")
            axis[i, 1].plot((0, m), (0, m), linestyle="--", alpha=0.5, color="black")
            axis[i, 1].scatter(
                a.loc[r1_demux],
                b.loc[r1_demux],
                s=5,
                alpha=0.5,
                rasterized=True,
                color=colors[0],
                label="round1 + 10X",
            )
        sns.despine(fig)
        fig.savefig(
            os.path.join(
                args.output_dir,
                output_prefix + "species_mix.plot.doublet_difference.svg",
            ),
            bbox_inches="tight",
            dpi=300,
        )

        #
        # Estimate number of cells per droplet
        df3 = pd.read_csv(umi_count_file)
        df3 = df3.assign(name=df3[args.cell_barcodes].sum(1))
        cells_per_droplet = df3.groupby(args.cell_barcodes[1])[
            args.cell_barcodes[0]
        ].nunique()

        from scipy.optimize import curve_fit

        def poisson(k, lamb):
            from scipy.special import factorial

            return np.exp(-lamb) * (lamb ** k / factorial(k))

        def neg_binom(k, r, p):
            from scipy.special import gammaln
            from scipy.special import factorial

            infinitesimal = np.finfo(np.float).eps

            # MLE estimate based on the formula on Wikipedia:
            # http://en.wikipedia.org/wiki/Negative_binomial_distribution#Maximum_likelihood_estimation
            result = (
                np.sum(gammaln(r + k))
                - np.sum(np.log(factorial(k)))
                - len(k) * gammaln(r)
                + len(k) * r * np.log(p)
                + np.sum(k * np.log(1 - (p if p < 1 else 1 - infinitesimal)))
            )

            return -result

        cells_per_droplet.to_csv(os.path.join(args.output_dir, output_prefix + "cells_per_droplet.csv"), header=False)
        cells_per_droplet = pd.read_csv(os.path.join(args.output_dir, output_prefix + "cells_per_droplet.csv"), index_col=0, header=None, squeeze=True)

        fig, axis = plt.subplots(2, 3, figsize=(3 * 3, 2 * 3), tight_layout=True)
        for ax in axis[0, :]:
            sns.distplot(cells_per_droplet, kde=False, ax=ax, norm_hist=False)

        cells_per_droplet_counts = cells_per_droplet.value_counts()
        print("Fitting Poisson model")
        lamb, cov_matrix = curve_fit(
            poisson, cells_per_droplet_counts.index, cells_per_droplet_counts
        )
        print(lamb)
        x = np.arange(
            0,
            cells_per_droplet_counts.index.max()
        )
        y_hat = scipy.stats.poisson.pmf(x, lamb)
        y_hat *= cells_per_droplet.shape[0]
        for ax in axis[0, :]:
            ax.plot(x + 0.5, y_hat)  # the 0.5 is just to center on the middle of the histogram bins
        axis[0, 1].set_yscale("log")
        axis[0, 2].set_xscale("log")
        axis[0, 2].set_yscale("log")
        v = cells_per_droplet_counts.max()
        v += v * 1
        axis[0, 1].set_ylim((1, v))
        axis[0, 2].set_ylim((1, v))

        # print("Fitting Negative binomial model")
        # params, cov_matrix = curve_fit(
        #     neg_binom, cells_per_droplet_counts.index.values, cells_per_droplet_counts.values
        # )
        # print(params)
        # params = NegativeBinomial().fit(cells_per_droplet_counts.values)

        # now only cells with >= 100 UMIs
        umis_per_cell = df3.groupby("name")["umi"].sum()
        real_cells = umis_per_cell[umis_per_cell >= 100].index
        stringent_cells_per_droplet = (
            df3.loc[df3["name"].isin(real_cells)]
            .groupby(args.cell_barcodes[1])[args.cell_barcodes[0]]
            .nunique()
        )

        for ax in axis[1, :]:
            sns.distplot(stringent_cells_per_droplet, kde=False, ax=ax, norm_hist=False)

        stringent_cells_per_droplet_counts = stringent_cells_per_droplet.value_counts()
        lamb, cov_matrix = curve_fit(
            poisson, stringent_cells_per_droplet_counts.index, stringent_cells_per_droplet_counts
        )
        print(lamb)

        x = np.arange(
            0,
            stringent_cells_per_droplet_counts.index.max()
        )
        y_hat = scipy.stats.poisson.pmf(x, lamb)
        y_hat *= stringent_cells_per_droplet.shape[0]
        for ax in axis[1, :]:
            ax.plot(x + 0.5, y_hat)  # the 0.5 is just to center on the middle of the histogram bins
        axis[1, 1].set_yscale("log")
        axis[1, 2].set_xscale("log")
        axis[1, 2].set_yscale("log")
        v = stringent_cells_per_droplet_counts.max()
        v += v * 1
        axis[1, 1].set_ylim((1, v))
        axis[1, 2].set_ylim((1, v))

        for ax in axis.flatten():
            ax.set_xlabel("Cells per droplet")
            ax.set_ylabel("Droplets")
        fig.savefig(
            os.path.join(args.output_dir, output_prefix + "poissonian_properties.svg"),
            bbox_inches="tight",
            dpi=300,
        )


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
