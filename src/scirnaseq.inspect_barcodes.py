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
import scipy


pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


def sorted_nicely(l):
    """ Sort a given iterable in the way that humans expect."""
    def convert(text): return int(text) if text.isdigit() else text

    def alphanum_key(key): return [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


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
        help="Name of run of experiment. Must match produced input files e.g. 'sciRNA_016'",
        type=str)
    parser.add_argument(
        "-a", "--annotation",
        dest="annotation",
        help="Path to CSV file with barcode annotation. " +
             "e.g. metadata/sciRNA-seq.oligos_2018-05-22v2.csv",
        type=str)
    default = ["round1", "round2"]
    parser.add_argument(
        "-b", "--barcodes",
        dest="cell_barcodes",
        help="Cell barcodes used. A comma-delimited (no spaces) list. " +
             " Default is '{}'.".format("','".join(default)),
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
        help="Minimum mismatches input files were allowed to have. " +
             "Default {}".format(default),
        default=default,
        type=int)
    default = 1
    parser.add_argument(
        "--max-mismatches",
        dest="max_mismatches",
        help="Maximum mismatches input files were allowed to have. "
             "Default {}".format(default),
        default=default,
        type=int)
    # args = parser.parse_args("-a metadata/sciRNA-seq.oligos_2018-11-14.csv -b round1,round2 --max-mismatches 0 -".split(" "))
    args = parser.parse_args("-a metadata/sciRNA-seq.SCI017.oligos_2019-02-11.csv -b round1,round2 --max-mismatches 0 sci-RNA-seq_SCI017_SSIV_6K".split(" "))
    print("# " + time.asctime() + " - Parsing command line arguments:")
    # args = parser.parse_args()
    args.cell_barcodes = args.cell_barcodes.split(",")
    print(args)
    return args


def main():
    print("# " + time.asctime() + " - Start.")

    args = get_cli_arguments()

    annotation = pd.read_csv(args.annotation)
    bc_10x = "/home/arendeiro/workspace/cellranger-atac-1.0.0/cellranger-atac-cs/1.0.0/lib/python/barcodes/737K-cratac-v1.txt"

    # # get barcodes with well information
    barcodes_with_well = list()
    for barcode in annotation['barcode_type'].unique():
        if not annotation.loc[annotation['barcode_type'] == barcode, 'plate_well'].isnull().all():
            barcodes_with_well.append(barcode)

    for mismatches in range(args.min_mismatches, args.max_mismatches):
        print("# " + time.asctime() + " - Processing files with {} mismatches.".format(mismatches))

        output_prefix = "{}.{}mis.".format(args.run_name, mismatches)
        print(output_prefix)

        # get files to concatenate into HDF
        hdf_file = os.path.join("barcodes", args.run_name + ".mis{}.hdf".format(mismatches))
        umi_count_file = os.path.join("barcodes", output_prefix + "barcode_gene_umi_count.clean.csv")
        umi_dup_file = os.path.join("barcodes", output_prefix + "barcode_umi_dups.count.csv")
        cell_dup_file = os.path.join("barcodes", output_prefix + "barcode_umi_dups.per_cell.csv")
        mapping_cell_rate_file = os.path.join("barcodes", output_prefix + "mapping_rate.per_cell.csv")
        mapping_well_rate_file = os.path.join("barcodes", output_prefix + "mapping_rate.per_well.csv")
        cell_metrics_file = os.path.join("barcodes", output_prefix + "all_metrics.per_cell.csv")
        if not os.path.exists(hdf_file):
            print("# " + time.asctime() + " - Concatenating barcode files into HDF file.")
            pieces = glob(os.path.join("barcodes", "{}.part_*.barcodes.*_*.mis_{}.csv.gz".format(
                args.run_name, mismatches)))
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
                glob(os.path.join("star", "{}.STAR.htseq-count_gene.read_gene.part_*.csv".format(args.run_name))))
            for piece in pieces:
                print(piece)
                t = pd.read_csv(piece, header=None, names=['read', 'gene'])
                transcriptome = transcriptome.append(t, ignore_index=True)

            # join barcodes and transcriptome
            print("# " + time.asctime() + " - Joining barcodes and transcriptome.")
            df = cells.set_index("read").join(transcriptome.set_index("read"))

            # investigate mapping rates
            print("# " + time.asctime() + " - Investigating mapping rates per cell.")

            # # per cell
            all_reads = df.groupby(args.cell_barcodes)['gene'].count()
            df.loc[:, 'unmapped'] = df['gene'].str.startswith("__").astype(int)
            # unmapped = df.groupby(args.cell_barcodes)['gene'].apply(lambda x: x.str.startswith("__").sum())
            unmapped = df.groupby(args.cell_barcodes)['unmapped'].sum()
            mapping_rate = 1 - (unmapped / all_reads)

            cell_rates = pd.DataFrame([all_reads, unmapped, mapping_rate], index=['all_reads', 'unmapped', 'mapping_rate']).T
            cell_rates.to_csv(mapping_cell_rate_file)

            # # per well
            print("# " + time.asctime() + " - Investigating mapping rates per well.")
            # # # first add well information
            for barcode in barcodes_with_well:
                a = (
                    annotation.loc[
                        annotation['barcode_type'] == barcode, ['barcode_sequence', 'plate_well']]
                    .rename(columns={"barcode_sequence": barcode, "plate_well": barcode + "_well"}))
                df = pd.merge(df, a, on=barcode, how='left')

            all_reads = df.groupby(['round1_well'])['gene'].count()
            unmapped = df.groupby(['round1_well'])['gene'].apply(lambda x: x.str.startswith("__").sum())
            mapping_rate = 1 - (unmapped / all_reads)

            well_rates = pd.DataFrame([all_reads, unmapped, mapping_rate], index=['all_reads', 'unmapped', 'mapping_rate']).T
            well_rates.to_csv(mapping_well_rate_file)

            # remove unmapped reads
            msg = " - Removing unmapped reads and collapsing to unique UMIs per cell per gene."
            print("# " + time.asctime() + msg)
            df2 = df.loc[df['unmapped'] != 1, :]

            print("# " + time.asctime() + " - Investigating duplicates.")
            duplicates_per_molecule = df2.groupby(args.cell_barcodes)['umi'].value_counts()
            duplicates_per_molecule.columns = args.cell_barcodes + ['umi', 'count']
            duplicates_per_molecule.to_csv(umi_dup_file)

            fig, axis = plt.subplots(1, 2, figsize=(2 * 3, 3), tight_layout=True)
            sns.distplot(duplicates_per_molecule['count'], ax=axis[0], kde=False)
            axis[0].set_xlabel("Reads per UMI")
            axis[0].set_ylabel("UMIs (log)")
            axis[0].set_yscale("log")
            sns.distplot(np.log10(duplicates_per_molecule['count']), ax=axis[1], kde=False)
            axis[1].set_xlabel("Reads per UMI (log)")
            axis[1].set_ylabel("UMIs (log)")
            # axis[1].set_xscale("log")
            axis[1].set_yscale("log")
            reads_per_umi_plot = os.path.join(
                "results", output_prefix + "barcode_umi_dups.per_cell.distplot.svg")
            fig.savefig(reads_per_umi_plot, dpi=300, bbox_inches="tight")

            duplicates_per_molecule.loc[:, 'unique'] = (duplicates_per_molecule['count'] == 1).astype(int)
            duplicates_per_cell = duplicates_per_molecule.groupby(args.cell_barcodes)['unique'].sum().to_frame(name="unique")
            duplicates_per_cell.loc[:, "all"] = duplicates_per_molecule.groupby(args.cell_barcodes)['count'].sum()

            duplicates_per_cell['unique_rate'] = duplicates_per_cell['unique'] / duplicates_per_cell['all']
            duplicates_per_cell.to_csv(cell_dup_file)

            fig, axis = plt.subplots(1, 1, figsize=(1 * 3, 3), tight_layout=True)
            sns.distplot(duplicates_per_cell['unique_rate'], ax=axis, kde=False)
            axis.set_xlabel("Unique rate per Cell")
            axis.set_ylabel("Cells")
            axis.set_yscale("log")
            cell_unique_rate_plot = os.path.join(
                "results", output_prefix + "cell_umi_dups.per_cell.distplot.svg")
            fig.savefig(cell_unique_rate_plot, dpi=300, bbox_inches="tight")

            # get unique UMIs (UMI counts per cell per gene)
            df3 = df2.groupby(args.cell_barcodes + ['gene'])['umi'].nunique()

            print("# " + time.asctime() + " - Writing gene expression matrix to file.")
            df3.to_csv(
                umi_count_file,
                index=True, header=True)

            # duplicates_per_molecule = pd.read_csv(umi_dup_file, index_col=[0, 1])
            # duplicates_per_cell = pd.read_csv(cell_dup_file, index_col=[0, 1])
            # cell_rates = pd.read_csv(mapping_cell_rate_file, index_col=[0, 1])
            # df3 = pd.read_csv(umi_count_file)

            umis_per_cell = df3.reset_index().groupby(args.cell_barcodes)['umi'].sum()
            genes_per_cell = df3.reset_index().groupby(args.cell_barcodes)['gene'].nunique()

            # plot again all metrics, including per-cell duplicate rate
            d = pd.concat([umis_per_cell, genes_per_cell, duplicates_per_cell.loc[:, 'unique_rate']], ignore_index=True, axis=1)
            p = pd.concat([cell_rates.loc[d.index, :].dropna(), d], ignore_index=True, axis=1)
            p.columns = cell_rates.columns.tolist() + ['umi', 'gene', 'unique_rate']
            p.to_csv(cell_metrics_file)

            fig, axis = plt.subplots(2, 3, figsize=(3 * 3, 3 * 2), tight_layout=True)

            sns.distplot(p['mapping_rate'], ax=axis[0, 0], kde=False)
            axis[0, 0].set_xlabel("Mapping rate")
            axis[0, 0].set_ylabel("Barcodes")
            axis[0, 0].set_yscale("log")

            axis[1, 0].scatter(
                p['all_reads'].head(int(8 * 50e4)), p['gene'].head(int(8 * 50e4)),
                alpha=0.2, s=1, rasterized=True)
            axis[1, 0].set_xlabel("Reads per barcode")
            axis[1, 0].set_ylabel("Genes per barcode")
            axis[1, 0].set_xscale("log")
            axis[1, 0].set_yscale("log")

            axis[0, 1].scatter(
                p['all_reads'].head(int(8 * 50e4)), p['mapping_rate'].head(int(8 * 50e4)),
                alpha=0.2, s=1, rasterized=True)
            axis[0, 1].set_xlabel("Reads per barcode")
            axis[0, 1].set_ylabel("Mapping rate")
            axis[0, 1].set_xscale("log")

            axis[0, 2].scatter(
                p['all_reads'].head(int(8 * 50e4)), p['unique_rate'].head(int(8 * 50e4)),
                alpha=0.2, s=1, rasterized=True)
            axis[0, 2].set_xscale("log")
            axis[0, 2].set_xlabel("Reads per barcode")
            axis[0, 2].set_ylabel("Unique rate")

            sns.distplot(np.log10(p['all_reads']), ax=axis[1, 1], kde=False)
            axis[1, 1].set_xlabel("Reads per barcode (log10)")
            axis[1, 1].set_yscale("log")
            axis[1, 1].set_ylabel("Barcodes")

            axis[1, 2].scatter(
                p['mapping_rate'].head(int(8 * 50e4)), p['unique_rate'].head(int(8 * 50e4)),
                alpha=0.2, s=1, rasterized=True)
            axis[1, 2].set_xlabel("Mapping rate")
            axis[1, 2].set_ylabel("Unique rate")

            reads_vs_mapping_rate_plot = os.path.join(
                "results", output_prefix + "reads_vs_mapping_rate_vs_duplication.per_cell.scatter.res.svg")
            fig.savefig(reads_vs_mapping_rate_plot, dpi=300, bbox_inches="tight")

            # Pairwise plots

            # # add type of material
            bc_annot = annotation.set_index("barcode_sequence").loc[:, "material"]
            bc_annot.index.name = "round2"
            p_annot = p.join(bc_annot, on="round1")

            # add species ratio and assigned species
            df3.loc[:, "human"] = df3.loc[:, "gene"].startswith("ENSG")

            sp_count = df3.groupby(args.cell_barcodes + ['human'])['umi'].sum()
            sp_ratio = sp_count.reset_index().groupby(args.cell_barcodes).apply(
                lambda x: x[x['human'] == True]['umi'] / x['umi'].sum())
            # sp_ratio = pd.read_csv("TMP.sp_ratio.csv")
            sp_ratio.loc[:, "species"] = np.nan
            sp_ratio.loc[sp_ratio > 0.75, "species"] = "human"
            sp_ratio.loc[sp_ratio > 0.25, "species"] = "mouse"
            sp_ratio.loc[:, "species"] = sp_ratio.loc[:, "species"].fillna("doublet")

            p2 = p_annot.loc[p_annot['umi'].nlargest(5000).index]
            g = sns.pairplot(p2, hue="material", plot_kws={"linewidth": 0, "edgecolors": "none", "rasterized": True, "alpha": 0.2})
            cell_unique_rate_plot = os.path.join(
                "results", output_prefix + "cell_umi_dups.per_cell.top_cells.pairplot.svg")
            g.fig.savefig(cell_unique_rate_plot, dpi=300, bbox_inches="tight")

            p2 = p_annot.loc[p_annot['umi'].nlargest(5000).index]
            p2.loc[:, ~p2.columns.str.contains("_rate|material")] = np.log10(p2.loc[:, ~p2.columns.str.contains("_rate|material")] + 1)
            g = sns.pairplot(p2, hue="material", plot_kws={"linewidth": 0, "edgecolors": "none", "rasterized": True, "alpha": 0.2})
            cell_unique_rate_plot = os.path.join(
                "results", output_prefix + "cell_umi_dups.per_cell.top_cells.log.pairplot.svg")
            g.fig.savefig(cell_unique_rate_plot, dpi=300, bbox_inches="tight")

            p2 = p_annot.loc[p_annot['umi'].nsmallest(5000).index]
            g = sns.pairplot(p2, hue="material", plot_kws={"linewidth": 0, "edgecolors": "none", "rasterized": True, "alpha": 0.2})
            cell_unique_rate_plot = os.path.join(
                "results", output_prefix + "cell_umi_dups.per_cell.bottom_cells.pairplot.svg")
            g.fig.savefig(cell_unique_rate_plot, dpi=300, bbox_inches="tight")

            p2 = p_annot.loc[p_annot['umi'].nsmallest(5000).index]
            p2.loc[:, ~p2.columns.str.contains("_rate|material")] = np.log10(p2.loc[:, ~p2.columns.str.contains("_rate|material")] + 1)
            g = sns.pairplot(p2, hue="material", plot_kws={"linewidth": 0, "edgecolors": "none", "rasterized": True, "alpha": 0.2})
            cell_unique_rate_plot = os.path.join(
                "results", output_prefix + "cell_umi_dups.per_cell.bottom_cells.log.pairplot.svg")
            g.fig.savefig(cell_unique_rate_plot, dpi=300, bbox_inches="tight")

            # Efficiency plot
            # # add 10X barcodes
            from Bio.Seq import Seq
            bc = pd.read_csv(bc_10x, header=None, squeeze=True)

            p = pd.read_csv(cell_metrics_file, index_col=[0, 1])
            p.loc[:, "round2_rc"] = [str(Seq(x).reverse_complement()) for x in p.index.get_level_values("round2")]
            p.loc[:, "round2_in_10X"] = p["round2_rc"].isin(bc.tolist())
            print("Fraction of correct round2 barcodes: {}".format(p["round2_in_10X"].sum() / float(p.shape[0])))

            # # count transcriptome reads per cell
            # from sklearn.linear_model import LinearRegression
            sns.set_style("ticks")

            def f(x, m, b):
                return m * x + b

            p2 = p.loc[p['umi'].nlargest(100000).index]

            fig, axis = plt.subplots(2, 2, figsize=(2 * 3, 2 * 3), tight_layout=True)
            fig.suptitle("Experiment efficiency:\n" + args.run_name, ha="center")

            for j, var in enumerate(["umi", "gene"]):
                for i, l in [(0, ""), (1, " (log)")]:
                    (m, b), pcov = scipy.optimize.curve_fit(f, p2['all_reads'], p2[var])
                    # lm = LinearRegression().fit(p2['all_reads'].values.reshape(-1, 1), p2[var])
                    # assert np.allclose(m, lm.coef_)
                    # assert np.allclose(b, lm.intercept_)
                    axis[i, j].text(0.5e5, 5000, s="y = {:.6f}x + {:.2f}".format(m, b), ha="center")

                    axis[i, j].set_xlabel("Total reads sequenced per cell" + l)
                    axis[i, j].set_ylabel("Useful {}s per cell".format(var.capitalize()) + l)
                    axis[i, j].scatter(p2['all_reads'], p2[var], alpha=0.2, s=2, edgecolors="none", rasterized=True)
                    x = np.linspace(p2['all_reads'].min(), p2['all_reads'].max(), num=1000)
                    axis[i, j].plot(x, f(x, m, b), color="orange")
                    y = f(1e4, m, b)
                    axis[i, j].text(1e4, y, s="{}s recovered\nwith 10.000\nreads sequenced:\n{:.2f}"
                                              .format(var.capitalize(), y), ha="left")
                    # X == Y
                    xmax = p2['all_reads'].max()
                    xmax += xmax * 0.1
                    ymax = p2[var].max()
                    ymax += ymax * 0.1
                    x = np.linspace(0, ymax, num=2)
                    y = f(x, 1, 0)
                    axis[i, j].plot(x, y, linestyle="--", color="black", linewidth=0.5)

                    # Axlines
                    for h in [100, 250, 500, 1000, 5000]:
                        axis[i, j].axhline(h, linestyle="--", color="black", linewidth=0.5)
                    for v in [10000, 100000]:
                        axis[i, j].axvline(v, linestyle="--", color="black", linewidth=0.5)

                    axis[i, j].ticklabel_format(useOffset=False, style="plain")
                    if i == 1:
                        axis[i, j].loglog()
            performance_plot = os.path.join(
                "results", output_prefix + "performance_plot.only_10x_correct.svg")
            fig.savefig(performance_plot, dpi=300, bbox_inches="tight")

        else:
            print("# " + time.asctime() + " - HDF file exists. Reading up '{}'".format(umi_count_file))
            df3 = pd.read_csv(umi_count_file)
        df3 = df3.reset_index()

        # investigate duplicates
        # duplicates_per_molecule[duplicates_per_molecule == 1].sum()  # number of unique UMIs
        # duplicates_per_molecule[duplicates_per_molecule > 1].shape[0]  # number of duplicate UMIs
        # duplicates_per_molecule[duplicates_per_molecule > 1].sum()  # total duplicate reads
        print("# " + time.asctime() + " - Investigating duplicates per well.")
        duplicates_per_molecule = pd.read_csv(umi_dup_file, header=None)
        duplicates_per_molecule.columns = ['round1', 'round2', 'umi', 'count']

        # # add well information
        for barcode in barcodes_with_well:
            a = (annotation.loc[
                annotation['barcode_type'] == barcode,
                ['barcode_sequence', 'plate_well']]
                .rename(columns={"barcode_sequence": barcode, "plate_well": barcode + "_well"}))
            duplicates_per_molecule = pd.merge(duplicates_per_molecule, a, on=barcode, how='left')

        # # # per cell
        # all_umis = duplicates_per_molecule.groupby("round1_well")['count'].count()
        # dup_umis = duplicates_per_molecule.groupby("round1_well")['count'].apply(lambda x: (x == 1).sum())
        # unique_rate = 1 - (dup_umis / all_umis)

        # cell_rates = pd.DataFrame([all_umis, dup_umis, unique_rate], index=['all_reads', 'dup_umis', 'unique_rate']).T
        # dup_cell_rate_file = os.path.join("barcodes", output_prefix + "barcode_umi_dups.per_cell.csv")
        # cell_rates.to_csv(dup_cell_rate_file)

        # # per well
        all_umis = duplicates_per_molecule.groupby("round1_well")['count'].count()
        dup_umis = duplicates_per_molecule.groupby("round1_well")['count'].apply(lambda x: (x == 1).sum())
        unique_rate = 1 - (dup_umis / all_umis)

        well_rates = pd.DataFrame([all_umis, dup_umis, unique_rate], index=['all_reads', 'dup_umis', 'unique_rate']).T
        dup_well_rate_file = os.path.join("barcodes", output_prefix + "barcode_umi_dups.per_well.csv")
        well_rates.to_csv(dup_well_rate_file)

        msg = " - mean unique UMI rate/well: '{}'".format(well_rates.loc[well_rates['all_reads'] > 200, 'unique_rate'].mean())
        print("## " + time.asctime() + msg)

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
        # fig.suptitle(args.run_name)
        # axis[0, 0].plot(umis_per_cell2.rank(ascending=False), umis_per_cell2, linestyle='-')
        axis[0, 0].plot(umis_per_cell.rank(ascending=False), umis_per_cell, linestyle='-', color="orange")
        # axis[0, 1].plot(umis_per_cell2.head(5000).rank(ascending=False), umis_per_cell2.head(5000), linestyle='-')
        axis[0, 1].plot(umis_per_cell.head(5000).rank(ascending=False), umis_per_cell.head(5000), linestyle='-', color="orange")
        # axis[1, 0].plot(umis_per_cell2.rank(ascending=False), umis_per_cell2, linestyle='-')
        axis[1, 0].plot(umis_per_cell.rank(ascending=False), umis_per_cell, linestyle='-', color="orange")
        axis[1, 0].set_yscale('log')
        # axis[1, 1].plot(umis_per_cell2.head(5000).rank(ascending=False), umis_per_cell2.head(5000), linestyle='-')
        axis[1, 1].plot(umis_per_cell.head(5000).rank(ascending=False), umis_per_cell.head(5000), linestyle='-', color="orange")
        axis[1, 1].set_yscale('log')
        for ax in axis.flatten():
            ax.axvline(5000, linestyle="--", color="black", alpha=0.5)
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
        species_count = pd.pivot_table(
            df3, index=args.cell_barcodes, columns="species", values="umi", aggfunc=sum, fill_value=0)

        species_count = species_count.join(p.set_index(args.cell_barcodes).loc[:, "round2_in_10X"])

        s = pd.DataFrame()
        for min_umi in [1, 2, 5, 10, 20, 30, 40, 50, 100, 500, 1000, 2500, 5000, 10000]:
            p = species_count.loc[species_count.sum(axis=1) >= min_umi, :]
            p = (p.T / p.sum(1)).T
            total = float(p.shape[0])
            s.loc['pure human', min_umi] = (p['human'] >= 0.8).sum() / total  # pure human cells
            s.loc['pure mouse', min_umi] = (p['mouse'] >= 0.8).sum() / total  # pure mouse cells
            s.loc['mixture', min_umi] = (((p['mouse'] < 0.8) & (p['human'] < 0.8))).sum() / total  # mixtures
        s *= 100

        fig, axis = plt.subplots(1, 1, figsize=(1 * 4, 4))
        fig.suptitle(args.run_name)
        sns.heatmap(s, ax=axis, cbar_kws={"label": "Percent of total"}, vmin=0, square=True)
        axis.set_xlabel("Minimum UMIs per cell")
        axis.set_yticklabels(axis.get_yticklabels(), rotation=0, ha="right", va="center")
        fig.savefig(
            os.path.join(
                args.output_dir, output_prefix + "species_mix.summary.heatmap.svg"), bbox_inches="tight", dpi=300)

        fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 4))
        fig.suptitle(args.run_name)
        sns.barplot(s[1].index, s[1], ax=axis[0])
        sns.barplot(s[5].index, s[100], ax=axis[1])
        axis[0].set_xlabel("Cells with at least 1 UMI")
        axis[1].set_xlabel("Cells with at least 100 UMI")
        axis[0].set_ylabel("Percent of total")
        axis[1].set_ylabel("Percent of total")
        sns.despine(fig)
        fig.savefig(
            os.path.join(
                args.output_dir, output_prefix + "species_mix.summary.barplot.svg"), bbox_inches="tight", dpi=300)

        fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 4))
        fig.suptitle(args.run_name)
        for ax in axis:
            ax.axhline(0, linestyle="--", color="black", alpha=0.1, zorder=100)
            ax.axvline(0, linestyle="--", color="black", alpha=0.1, zorder=100)
        m = max(species_count.loc[:, "mouse"].max(), species_count.loc[:, "human"].max())
        axis[0].plot((0, m), (0, m), linestyle="--", alpha=0.5, color="black")

        for bool_, color in zip([True], ["orange"]):
            axis[0].scatter(
                species_count.loc[species_count['round2_in_10X'] == bool_, "mouse"],
                species_count.loc[species_count['round2_in_10X'] == bool_, "human"],
                s=2, alpha=0.2, rasterized=True, label=bool_)
        axis[0].legend()
        axis[0].set_xlim((0, m))
        axis[0].set_ylim((0, m))
        a = np.log2(1 + species_count.loc[:, "mouse"])
        b = np.log2(1 + species_count.loc[:, "human"])
        m = max(a.max(), b.max())
        axis[1].plot((0, m), (0, m), linestyle="--", alpha=0.5, color="black")
        for bool_, color in zip([True], ["orange"]):
            a = np.log2(1 + species_count.loc[species_count['round2_in_10X'] == bool_, "mouse"])
            b = np.log2(1 + species_count.loc[species_count['round2_in_10X'] == bool_, "human"])
            axis[1].scatter(a, b, s=2, alpha=0.2, rasterized=True, label=bool_)
        axis[1].legend()
        axis[1].set_xlim((0, m))
        axis[1].set_ylim((0, m))
        axis[0].set_xlabel("Mouse (UMIs)")
        axis[0].set_ylabel("Human (UMIs)")
        axis[1].set_xlabel("Mouse (log2 UMIs)")
        axis[1].set_ylabel("Human (log2 UMIs)")
        sns.despine(fig)
        fig.savefig(
            os.path.join(
                args.output_dir, output_prefix + "species_mix.plot.svg"), bbox_inches="tight", dpi=300)

        grid = sns.jointplot(
            species_count.loc[:, "mouse"],
            species_count.loc[:, "human"], rasterized=True,
            joint_kws={"alpha": 0.2, "s": 2})
        grid.ax_joint.set_xlabel("Mouse (UMIs)")
        grid.ax_joint.set_ylabel("Human (UMIs)")
        grid.fig.suptitle(args.run_name)
        sns.despine(grid.fig)
        grid.fig.savefig(
            os.path.join(
                args.output_dir, output_prefix + "species_mix.jointplot.svg"), bbox_inches="tight", dpi=300)
        grid = sns.jointplot(
            np.log2(1 + species_count.loc[:, "mouse"]),
            np.log2(1 + species_count.loc[:, "human"]), rasterized=True,
            joint_kws={"alpha": 0.2, "s": 2})
        grid.ax_joint.set_xlabel("Mouse (log2 UMIs)")
        grid.ax_joint.set_ylabel("Human (log2 UMIs)")
        sns.despine(grid.fig)
        grid.fig.savefig(
            os.path.join(
                args.output_dir, output_prefix + "species_mix.jointplot.log2.svg"), bbox_inches="tight", dpi=300)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
