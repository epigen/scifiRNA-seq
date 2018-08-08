#!/usr/bin/env python

"""
sciRNA-seq barcode inspecting script.
"""


import os
import sys
import glob
import time
from argparse import ArgumentParser
import numpy as np
import pandas as pd
import pysam
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob


pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


def sorted_nicely(l): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)] 
    return sorted(l, key = alphanum_key)


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
default = 1
parser.add_argument(
    "--max-mismatches",
    dest="mismatches",
    help="Maximum mismatches input files were allowed to have. Default {}".format(default),
    default=default,
    type=int)
args = parser.parse_args()
print("# " + time.asctime() + " - Start.")
print(args)

annotation = pd.read_csv(args.annotation)
output_prefix = "{}.{}mis.".format(args.run_name, args.mismatches)

# get files to concatenate into HDF
print("# " + time.asctime() + " - Concatenating barcode files into HDF file.")
pieces = sorted_nicely(
    glob(os.path.join("barcodes", "{}_*.barcodes.*_*.mis_{}.csv.gz".format(args.run_name, args.mismatches))))
hdf_file = os.path.join("barcodes", args.run_name + ".mis{}.hdf".format(args.mismatches))
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
    glob(os.path.join("star", "{}.STAR.htseq-count.read_gene.part_*.csv".format(args.run_name))))
    # glob(os.path.join("star", "{}_STAR_*_part.Aligned.htseq-count.read_gene.csv".format(args.run_name))))
for piece in pieces:
    print(piece)
    t = pd.read_csv(piece, header=None, names=['read', 'gene'])
    transcriptome = transcriptome.append(t, ignore_index=True)


# join barcodes and transcriptome
print("# " + time.asctime() + " - Joining barcodes and transcriptome.")
df = cells.set_index("read").join(transcriptome.set_index("read"))

# remove unmapped reads
print("# " + time.asctime() + " - Removing unmapped reads and collapsing to unique UMIs per cell per gene.")
df2 = df.loc[~df['gene'].str.startswith("__"), :]

# get unique UMIs (UMI counts per cell per gene)
df3 = df2.groupby([args.cell_barcodes] + ['gene'])['umi'].nunique()
print("# " + time.asctime() + " - Writing gene expression matrix to file.")
df3.to_csv(
    os.path.join("barcodes", output_prefix + "barcode_gene_umi_count.clean.csv"),
    index=True, header=True)


# UMIs per cell (regardless of gene)
print("# " + time.asctime() + " - Investigating duplicates.")
umis_per_cell = df3.groupby(level=[args.cell_barcodes]).sum().sort_values(ascending=False)


# investigate duplicates
duplicates_per_molecule = df2.groupby([args.cell_barcodes])['umi'].value_counts()

duplicates_per_molecule[duplicates_per_molecule == 1].sum()  # number of unique UMIs
duplicates_per_molecule[duplicates_per_molecule > 1].shape[0]  # number of duplicate UMIs
duplicates_per_molecule[duplicates_per_molecule > 1].sum()  # total duplicate reads

bins = [int(np.round(x, 0)) for x in [0, 1] + list(np.logspace(1, 4.5, 20))]
dups = pd.cut(duplicates_per_molecule, bins).value_counts()

# Plot duplicates
fig, axis = plt.subplots(2, 1, figsize=(2 * 4, 1 * 4))
sns.barplot(dups.index, dups, ax=axis[0], palette="summer")
sns.barplot(dups.index, dups, ax=axis[1], palette="summer")
axis[1].set_yscale('log')
axis[0].set_xticklabels(axis[0].get_xticklabels(), visible=False)
axis[1].set_xlabel('Duplicates per UMI (interval)')
for ax in axis:
    ax.set_ylabel('Number of UMIs')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
sns.despine(fig)
fig.savefig(
    os.path.join(output_dir, output_prefix + "umis_duplication.barplot.svg"),
    bbox_inches="tight")


# Plot barcode distribution
print("# " + time.asctime() + " - Investigating barcode distributions.")
fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 1 * 4))
axis[0].plot(umis_per_cell.rank(ascending=False), umis_per_cell, linestyle='-')
axis[1].plot(umis_per_cell.head(200).rank(ascending=False), umis_per_cell.head(200), linestyle='-')
for ax in axis:
    ax.set_xlabel("Cells")
    ax.set_ylabel("UMIs")
sns.despine(fig)
fig.savefig(
    os.path.join(output_dir, output_prefix + "umis_per_cell.lineplot.svg"),
    bbox_inches="tight")


# Species mixing experiment
print("# " + time.asctime() + " - Species mixing experiment.")
df4 = df3.reset_index()
df4['species'] = df4['gene'].str.startswith("ENSG").replace(True, "human").replace(False, "mouse")


# Species mixing experiment
species_count = pd.pivot_table(df4, index=args.cell_barcodes, columns="species", values="umi", aggfunc=sum, fill_value=0)

s = pd.DataFrame()
for min_umi in list(range(1, 11)) + [20, 30, 40, 50, 100]:
    p = species_count.loc[species_count.sum(axis=1) >= min_umi, :]
    total = float(p.shape[0])
    s.loc['pure human', min_umi] = ((p['human'] >= min_umi) & (p['mouse'] == 0)).sum() / total # pure human cells
    s.loc['pure mouse', min_umi] = ((p['mouse'] >= min_umi) & (p['human'] == 0)).sum() / total # pure mouse cells
    s.loc['mixture', min_umi] = (((p['mouse'] >= min_umi) & (p['human'] >= min_umi))).sum() / total # mixtures

fig, axis = plt.subplots(1, 1, figsize=(1 * 4, 4))
sns.heatmap(s, ax=axis, cbar_kws={"label": "Fraction of total"}, vmin=0, square=True)
axis.set_xlabel("Minimum UMIs per cell")
axis.set_yticklabels(axis.get_yticklabels(), rotation=0, ha="right", va="center")
fig.savefig(os.path.join(output_dir, output_prefix + "species_mix.summary.heatmap.svg"), bbox_inches="tight", dpi=300)

fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 4))
sns.barplot(s[1].index, s[1], ax=axis[0])
sns.barplot(s[5].index, s[5], ax=axis[1])
axis[0].set_xlabel("Cells with at least 1 UMI")
axis[1].set_xlabel("Cells with at least 5 UMI")
axis[0].set_ylabel("Fraction of total")
axis[1].set_ylabel("Fraction of total")
sns.despine(fig)
fig.savefig(os.path.join(output_dir, output_prefix + "species_mix.summary.barplot.svg"), bbox_inches="tight", dpi=300)


fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 4))
for ax in axis:
    ax.axhline(0, linestyle="--", color="black", alpha=0.1, zorder=100)
    ax.axvline(0, linestyle="--", color="black", alpha=0.1, zorder=100)
axis[0].scatter(
    species_count.loc[:, "mouse"],
    species_count.loc[:, "human"], s=2, alpha=0.2, rasterized=True)
axis[1].scatter(
    np.log2(1 + species_count.loc[:, "mouse"]),
    np.log2(1 + species_count.loc[:, "human"]), s=2, alpha=0.2, rasterized=True)
axis[0].set_xlabel("Mouse (UMIs)")
axis[0].set_ylabel("Human (UMIs)")
axis[1].set_xlabel("Mouse (log2 UMIs)")
axis[1].set_ylabel("Human (log2 UMIs)")
sns.despine(fig)
fig.savefig(os.path.join(output_dir, output_prefix + "species_mix.plot.svg"), bbox_inches="tight", dpi=300)


grid = sns.jointplot(
    species_count.loc[:, "mouse"],
    species_count.loc[:, "human"], rasterized=True,
    joint_kws={"alpha": 0.2, "s": 2})
sns.despine(grid.fig)
grid.fig.savefig(os.path.join(output_dir, output_prefix + "species_mix.jointplot.svg"), bbox_inches="tight", dpi=300)
grid = sns.jointplot(
    np.log2(1 + species_count.loc[:, "mouse"]),
    np.log2(1 + species_count.loc[:, "human"]), rasterized=True,
    joint_kws={"alpha": 0.2, "s": 2})
sns.despine(grid.fig)
grid.fig.savefig(os.path.join(output_dir, output_prefix + "species_mix.jointplot.log2.svg"), bbox_inches="tight", dpi=300)
