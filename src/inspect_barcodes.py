#!/usr/bin/env python

"""
sciRNA-seq barcode inspecting script.
"""

import os
import sys
import glob
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
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

run = "BSF_0477_HJ7J7BGX2"
mismatches = 0
pieces = 4
step = 5000000
annotation_file = "metadata/sciRNA-seq.oligos_2018-05-22v2.csv"
annotation = pd.read_csv(annotation_file)
cell_barcodes = ["round1", "round2", "round3a", "round3b"]


# get files to concatenate into HDF
pieces = sorted_nicely(glob("barcodes/BSF_0477_HJ7J7BGX2_*.barcodes.*_*.mis_{}.csv.gz".format(mismatches)))
hdf_file = os.path.join("barcodes", run + ".mis{}.hdf".format(mismatches))
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
cells = pd.read_hdf(hdf_file, "cells")


# read transcriptome
transcriptome = pd.DataFrame()
for piece in range(2, 5):
    print(piece)
    t = pd.read_csv(os.path.join("star", "BSF_0477_HJ7J7BGX2_STAR.htseq-count.read_gene.part_{}.csv".format(piece)), header=None, names=['read', 'gene'])
    transcriptome = transcriptome.append(t, ignore_index=True)


# join barcodes and transcriptome
df = cells.set_index("read").join(transcriptome.set_index("read"))

# remove unmapped reads
df2 = df.loc[~df['gene'].str.startswith("__"), :]

# get unique UMIs (UMI counts per cell per gene)
df3 = df2.groupby(['round1', 'round2', 'round3a', 'round3b', 'gene'])['umi'].nunique()
df3.to_csv(os.path.join("barcodes", run + ".mis{}.barcode_gene_umi_count.clean.csv".format(mismatches)), index=True, header=True)


# UMIs per cell (regardless of gene)
umis_per_cell = df3.groupby(level=['round1', 'round2', 'round3a', 'round3b']).sum().sort_values(ascending=False)


# investigate duplicates
duplicates_per_molecule = df2.groupby(['round1', 'round2', 'round3a', 'round3b'])['umi'].value_counts()

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
    os.path.join("results", run + ".{}mis.umis_duplication.barplot.svg".format(mismatches)),
    bbox_inches="tight")


# Plot barcode distribution
fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 1 * 4))
axis[0].plot(umis_per_cell.rank(ascending=False), umis_per_cell, linestyle='-')
axis[1].plot(umis_per_cell.head(200).rank(ascending=False), umis_per_cell.head(200), linestyle='-')
for ax in axis:
    ax.set_xlabel("Cells")
    ax.set_ylabel("UMIs")
sns.despine(fig)
fig.savefig(
    os.path.join("results", run + ".{}mis.umis_per_cell.lineplot.svg".format(mismatches)),
    bbox_inches="tight")


# Species mixing experiment
df4 = df3.reset_index()
df4['species'] = df4['gene'].str.startswith("ENSG").replace(True, "human").replace(False, "mouse")


# Species mixing experiment
species_count = pd.pivot_table(df4, index=cell_barcodes, columns="species", values="umi", aggfunc=sum, fill_value=0)

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
fig.savefig(os.path.join("results", run + ".{}mis.species_mix.summary.heatmap.svg".format(mismatches)), bbox_inches="tight", dpi=300)

fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 4))
sns.barplot(s[1].index, s[1], ax=axis[0])
sns.barplot(s[5].index, s[5], ax=axis[1])
axis[0].set_xlabel("Cells with at least 1 UMI")
axis[1].set_xlabel("Cells with at least 5 UMI")
axis[0].set_ylabel("Fraction of total")
axis[1].set_ylabel("Fraction of total")
sns.despine(fig)
fig.savefig(os.path.join("results", run + ".{}mis.species_mix.summary.barplot.svg".format(mismatches)), bbox_inches="tight", dpi=300)


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
fig.savefig(os.path.join("results", run + ".{}mis.species_mix.plot.svg".format(mismatches)), bbox_inches="tight", dpi=300)


grid = sns.jointplot(
    species_count.loc[:, "mouse"],
    species_count.loc[:, "human"], rasterized=True,
    joint_kws={"alpha": 0.2, "s": 2})
sns.despine(grid.fig)
grid.fig.savefig(os.path.join("results", run + ".{}mis.species_mix.jointplot.svg".format(mismatches)), bbox_inches="tight", dpi=300)
grid = sns.jointplot(
    np.log2(1 + species_count.loc[:, "mouse"]),
    np.log2(1 + species_count.loc[:, "human"]), rasterized=True,
    joint_kws={"alpha": 0.2, "s": 2})
sns.despine(grid.fig)
grid.fig.savefig(os.path.join("results", run + ".{}mis.species_mix.jointplot.log2.svg".format(mismatches)), bbox_inches="tight", dpi=300)
