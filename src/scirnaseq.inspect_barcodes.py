#!/usr/bin/env python

"""
sciRNA-seq barcode inspecting script.
"""

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')


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


pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


def get_cli_arguments():
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="python inspect_barcodes.py",
        description="\n".join([
            "sciRNA-seq script from Bock lab. " +
            "See https://github.com/epigen/sciRNA-seq " +
            "for specific documentation."]))
    parser.add_argument(
        dest="run_name",
        help="Name of run of experiment. Must match produced input files e.g. BSF_0477_HJ7J7BGX2",
        type=str)
    parser.add_argument(
        "-a", "--annotation",
        dest="annotation",
        help="Path to CSV file with barcode annotation. e.g. metadata/sciRNA-seq.oligos_2018-05-22v2.csv",
        type=str)
    default = ["round1", "round2", "round3a", "round3b"]
    parser.add_argument(
        "-b", "--barcodes",
        dest="cell_barcodes",
        help="Cell barcodes used. A comma-delimited (no spaces) list. Default is '{}'.".format("','".join(default)),
        type=str)
    default = "results"
    parser.add_argument(
        "-o", "--output-dir",
        dest="output_dir",
        help="Output directory. Default is '{}'.".format(default),
        default=default,
        type=str)
    default = 0
    parser.add_argument(
        "--min-mismatches",
        dest="min_mismatches",
        help="Minimum mismatches input files were allowed to have. Default {}".format(default),
        default=default,
        type=int)
    default = 1
    parser.add_argument(
        "--max-mismatches",
        dest="max_mismatches",
        help="Maximum mismatches input files were allowed to have. Default {}".format(default),
        default=default,
        type=int)
    # args = parser.parse_args("-a metadata/sciRNA-seq.oligos_2018-11-14.csv -b round1,round2 --max-mismatches 0 -".split(" "))
    print("# " + time.asctime() + " - Parsing command line arguments:")
    args = parser.parse_args()
    args.cell_barcodes = args.cell_barcodes.split(",")
    print(args)
    return args


def main():
    print("# " + time.asctime() + " - Start.")

    args = get_cli_arguments()

    annotation = pd.read_csv(args.annotation)

    for mismatches in range(args.min_mismatches, args.max_mismatches):
        print("# " + time.asctime() + " - Processing files with {} mismatches.".format(mismatches))

        output_prefix = "{}.{}mis.".format(args.run_name, mismatches)
        print(output_prefix)

        # get files to concatenate into HDF
        hdf_file = os.path.join("barcodes", args.run_name + ".mis{}.hdf".format(mismatches))
        umi_count_file = os.path.join("barcodes", output_prefix + "barcode_gene_umi_count.clean.csv")
        umi_dup_file = os.path.join("barcodes", output_prefix + "barcode_umi_dups.count.csv")
        if not os.path.exists(hdf_file):
            print("# " + time.asctime() + " - Concatenating barcode files into HDF file.")
            pieces = glob(os.path.join("barcodes", "{}.part_*.barcodes.*_*.mis_{}.csv.gz".format(args.run_name, mismatches)))
            pieces = sorted_nicely(pieces)

            print("## " + time.asctime() + " - barcode files to read: '{}'.".format(", ".join(pieces)))
            print(pieces[0])
            c = pd.read_csv(
                pieces[0],
                compression="gzip")
            c.to_hdf(hdf_file, "cells", mode='w', format='table')
            del c  # allow garbage to be collected

            for piece in pieces[1:]:
                print(piece)
                c = pd.read_csv(
                    piece,
                    compression="gzip")
                c.to_hdf(hdf_file, "cells", append=True)
                del c  # allow garbage to be collected

            # read up HDF
            print("# " + time.asctime() + " - Reading up HDF file.")
            cells = pd.read_hdf(hdf_file, "cells")

            # read transcriptome
            print("# " + time.asctime() + " - Concatenating transcriptome files.")
            transcriptome = pd.DataFrame()
            pieces = sorted_nicely(
                # sci-RNA-seq_SCI011_gate_more_3_BSF_0513_Jurkat_3T3.STAR.htseq-count.read_gene.part_1.csv
                glob(os.path.join("star", "{}.STAR.htseq-count.read_gene.part_*.csv".format(args.run_name))))
            # glob(os.path.join("star", "{}_STAR_*_part.Aligned.htseq-count.read_gene.csv".format(args.run_name))))
            for piece in pieces:
                print(piece)
                t = pd.read_csv(piece, header=None, names=['read', 'gene'])
                transcriptome = transcriptome.append(t, ignore_index=True)

            # join barcodes and transcriptome
            print("# " + time.asctime() + " - Joining barcodes and transcriptome.")
            df = cells.set_index("read").join(transcriptome.set_index("read"))

            # investigate mapping rates
            print("# " + time.asctime() + " - Investigating mapping rates per well.")
            # # per well
            # add well information
            for barcode in args.cell_barcodes:
                a = annotation.loc[annotation['barcode_type'] == barcode, ['barcode_sequence', 'plate_well']].rename(columns={"barcode_sequence": barcode, "plate_well": barcode + "_well"})
                df = pd.merge(df, a, on=barcode, how='left')
            all_reads = df.groupby(['round1_well'])['gene'].count()
            unmapped = df.groupby(['round1_well'])['gene'].apply(lambda x: x.str.startswith("__").sum())
            mapping_rate = 1 - (unmapped / all_reads)

            rates = pd.DataFrame([all_reads, unmapped, mapping_rate], index=['all_reads', 'unmapped', 'mapping_rate']).T
            mapping_rate_file = os.path.join("barcodes", output_prefix + "mapping_rate.per_well.csv")
            rates.to_csv(mapping_rate_file)

            # remove unmapped reads
            print("# " + time.asctime() + " - Removing unmapped reads and collapsing to unique UMIs per cell per gene.")
            df2 = df.loc[~df['gene'].str.startswith("__"), :]

            print("# " + time.asctime() + " - Investigating duplicates.")
            duplicates_per_molecule = df2.groupby(args.cell_barcodes)['umi'].value_counts()

            duplicates_per_molecule.to_csv(umi_dup_file)

            # get unique UMIs (UMI counts per cell per gene)
            df3 = df2.groupby(args.cell_barcodes + ['gene'])['umi'].nunique()

            print("# " + time.asctime() + " - Writing gene expression matrix to file.")
            df3.to_csv(
                umi_count_file,
                index=True, header=True)
        else:
            print("# " + time.asctime() + " - HDF file exists. Reading up '{}'".format(umi_count_file))
            df3 = pd.read_csv(umi_count_file)
        df3 = df3.reset_index()

        # investigate duplicates
        # duplicates_per_molecule[duplicates_per_molecule == 1].sum()  # number of unique UMIs
        # duplicates_per_molecule[duplicates_per_molecule > 1].shape[0]  # number of duplicate UMIs
        # duplicates_per_molecule[duplicates_per_molecule > 1].sum()  # total duplicate reads

        # # duplicates per well
        print("# " + time.asctime() + " - Investigating duplicates per well.")
        duplicates_per_molecule = pd.read_csv(umi_dup_file, header=None)
        duplicates_per_molecule.columns = ['round1', 'round2', 'umi', 'count']

        # add well information
        for barcode in args.cell_barcodes:
            a = annotation.loc[annotation['barcode_type'] == barcode, ['barcode_sequence', 'plate_well']].rename(columns={"barcode_sequence": barcode, "plate_well": barcode + "_well"})
            duplicates_per_molecule = pd.merge(duplicates_per_molecule, a, on=barcode, how='left')

        all_umis = duplicates_per_molecule.groupby("round1_well")['count'].count()
        dup_umis = duplicates_per_molecule.groupby("round1_well")['count'].apply(lambda x: (x == 1).sum())
        unique_rate = 1 - (dup_umis / all_umis)

        rates = pd.DataFrame([all_umis, dup_umis, unique_rate], index=['all_reads', 'dup_umis', 'unique_rate']).T
        dup_rate_file = os.path.join("barcodes", output_prefix + "barcode_umi_dups.per_well.csv")
        rates.to_csv(dup_rate_file)

        print("## " + time.asctime() + " - mean unique UMI rate: '{}'".format(rates.loc[rates['all_reads'] > 200, 'unique_rate'].mean()))

        # bins = [int(np.round(x, 0)) for x in [0, 1] + list(np.logspace(1, 4.5, 20))]
        # dups = pd.cut(duplicates_per_molecule[3], bins).value_counts()

        # # Plot duplicates
        # fig, axis = plt.subplots(2, 1, figsize=(2 * 4, 1 * 4))
        # sns.barplot(dups.index, dups, ax=axis[0], palette="summer")
        # sns.barplot(dups.index, dups, ax=axis[1], palette="summer")
        # axis[1].set_yscale('log')
        # axis[0].set_xticklabels(axis[0].get_xticklabels(), visible=False)
        # axis[1].set_xlabel('Duplicates per UMI (interval)')
        # for ax in axis:
        #     ax.set_ylabel('Number of UMIs')
        #     ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        # sns.despine(fig)
        # fig.savefig(
        #     os.path.join(args.output_dir, output_prefix + "umis_duplication.barplot.svg"),
        #     bbox_inches="tight")

        # UMIs per cell (regardless of gene)
        umis_per_cell = df3.groupby(args.cell_barcodes)['umi'].sum().sort_values(ascending=False)

        # Plot barcode distribution
        print("# " + time.asctime() + " - Investigating barcode distributions.")
        fig, axis = plt.subplots(2, 2, figsize=(2 * 4, 2 * 4))
        fig.suptitle(args.run_name)
        axis[0, 0].plot(umis_per_cell.rank(ascending=False), umis_per_cell, linestyle='-')
        axis[0, 1].plot(umis_per_cell.head(1000).rank(ascending=False), umis_per_cell.head(1000), linestyle='-')
        axis[1, 0].plot(umis_per_cell.rank(ascending=False), umis_per_cell, linestyle='-')
        axis[1, 0].set_yscale('log')
        axis[1, 1].plot(umis_per_cell.head(1000).rank(ascending=False), umis_per_cell.head(1000), linestyle='-')
        axis[1, 1].set_yscale('log')
        for ax in axis.flatten():
            ax.axvline(1000, linestyle="--", color="black", alpha=0.5)
            ax.set_xlabel("Cells")
            ax.set_ylabel("UMIs")
        sns.despine(fig)
        fig.savefig(
            os.path.join(args.output_dir, output_prefix + "umis_per_cell.lineplot.svg"),
            bbox_inches="tight")

        # Species mixing experiment
        print("# " + time.asctime() + " - Species mixing experiment.")

        df3['species'] = np.nan
        df3.loc[df3['gene'].str.startswith("ENSG"), "species"] = 'human'
        df3.loc[df3['gene'].str.startswith("ENSMUS"), "species"] = 'mouse'

        # Species mixing experiment
        species_count = pd.pivot_table(df3, index=args.cell_barcodes, columns="species", values="umi", aggfunc=sum, fill_value=0)

        s = pd.DataFrame()
        for min_umi in [1, 2, 5, 10, 20, 30, 40, 50, 100, 500, 1000, 2500, 5000, 10000]:
            p = species_count.loc[species_count.sum(axis=1) >= min_umi, :]
            p = (p.T / p.sum(1)).T
            total = float(p.shape[0])
            s.loc['pure human', min_umi] = (p['human'] >= 0.8).sum() / total # pure human cells
            s.loc['pure mouse', min_umi] = (p['mouse'] >= 0.8).sum() / total # pure mouse cells
            s.loc['mixture', min_umi] = (((p['mouse'] < 0.8) & (p['human'] < 0.8))).sum() / total # mixtures
        s *= 100

        fig, axis = plt.subplots(1, 1, figsize=(1 * 4, 4))
        fig.suptitle(args.run_name)
        sns.heatmap(s, ax=axis, cbar_kws={"label": "Percent of total"}, vmin=0, square=True)
        axis.set_xlabel("Minimum UMIs per cell")
        axis.set_yticklabels(axis.get_yticklabels(), rotation=0, ha="right", va="center")
        fig.savefig(os.path.join(args.output_dir, output_prefix + "species_mix.summary.heatmap.svg"), bbox_inches="tight", dpi=300)

        fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 4))
        fig.suptitle(args.run_name)
        sns.barplot(s[1].index, s[1], ax=axis[0])
        sns.barplot(s[5].index, s[100], ax=axis[1])
        axis[0].set_xlabel("Cells with at least 1 UMI")
        axis[1].set_xlabel("Cells with at least 100 UMI")
        axis[0].set_ylabel("Percent of total")
        axis[1].set_ylabel("Percent of total")
        sns.despine(fig)
        fig.savefig(os.path.join(args.output_dir, output_prefix + "species_mix.summary.barplot.svg"), bbox_inches="tight", dpi=300)

        fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 4))
        fig.suptitle(args.run_name)
        for ax in axis:
            ax.axhline(0, linestyle="--", color="black", alpha=0.1, zorder=100)
            ax.axvline(0, linestyle="--", color="black", alpha=0.1, zorder=100)
        m = max(species_count.loc[:, "mouse"].max(), species_count.loc[:, "human"].max())
        axis[0].plot((0, m), (0, m), linestyle="--", alpha=0.5, color="black")
        axis[0].scatter(
            species_count.loc[:, "mouse"],
            species_count.loc[:, "human"], s=2, alpha=0.2, rasterized=True)
        axis[0].set_xlim((0, m))
        axis[0].set_ylim((0, m))
        a = np.log2(1 + species_count.loc[:, "mouse"])
        b = np.log2(1 + species_count.loc[:, "human"])
        m = max(a.max(), b.max())
        axis[1].plot((0, m), (0, m), linestyle="--", alpha=0.5, color="black")
        axis[1].scatter(a, b, s=2, alpha=0.2, rasterized=True)
        axis[1].set_xlim((0, m))
        axis[1].set_ylim((0, m))
        axis[0].set_xlabel("Mouse (UMIs)")
        axis[0].set_ylabel("Human (UMIs)")
        axis[1].set_xlabel("Mouse (log2 UMIs)")
        axis[1].set_ylabel("Human (log2 UMIs)")
        sns.despine(fig)
        fig.savefig(os.path.join(args.output_dir, output_prefix + "species_mix.plot.svg"), bbox_inches="tight", dpi=300)

        grid = sns.jointplot(
            species_count.loc[:, "mouse"],
            species_count.loc[:, "human"], rasterized=True,
            joint_kws={"alpha": 0.2, "s": 2})
        grid.ax_joint.set_xlabel("Mouse (UMIs)")
        grid.ax_joint.set_ylabel("Human (UMIs)")
        grid.fig.suptitle(args.run_name)
        sns.despine(grid.fig)
        grid.fig.savefig(os.path.join(args.output_dir, output_prefix + "species_mix.jointplot.svg"), bbox_inches="tight", dpi=300)
        grid = sns.jointplot(
            np.log2(1 + species_count.loc[:, "mouse"]),
            np.log2(1 + species_count.loc[:, "human"]), rasterized=True,
            joint_kws={"alpha": 0.2, "s": 2})
        grid.ax_joint.set_xlabel("Mouse (log2 UMIs)")
        grid.ax_joint.set_ylabel("Human (log2 UMIs)")
        sns.despine(grid.fig)
        grid.fig.savefig(os.path.join(args.output_dir, output_prefix + "species_mix.jointplot.log2.svg"), bbox_inches="tight", dpi=300)


def sorted_nicely(l):
    """ Sort a given iterable in the way that humans expect."""
    def convert(text): return int(text) if text.isdigit() else text

    def alphanum_key(key): return [convert(c) for c in re.split('([0-9]+)', key)]

    return sorted(l, key=alphanum_key)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
